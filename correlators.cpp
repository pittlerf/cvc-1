#define MAIN_PROGRAM
#include "cvc_global.h"
#undef MAIN_PROGRAM

#include "DependencyGraph.hpp"
#include "DependencyResolving.hpp"
#include "meta_types.hpp"
#include "algorithms.hpp"
#include "yaml_parsers.hpp"
#include "Core.hpp"
#include "Logger.hpp"
#include "debug_printf.hpp"
#include "init_g_gauge_field.hpp"
#include "cvc_utils.h"
#include "mpi_init.h"
#include "enums.hpp"
#include "ParallelMT19937_64.hpp"
#include "h5utils.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <yaml-cpp/yaml.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <exception>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range.hpp>

using namespace cvc;

std::map< std::string, std::vector< std::vector<int> > > momentum_lists;
std::map< std::string, stoch_prop_meta_t > ts_prop_meta; 

int main(int argc, char ** argv){
  int exitstatus;
  int io_proc = -1;
  double * gauge_field_with_phase = NULL;

  Core core(argc,argv,"correlators");
  if( !core.is_initialised() ){
		debug_printf(g_proc_id,0,"[correlators] Core initialisation failure!\n");
		return(CVC_EXIT_CORE_INIT_FAILURE);
  }

  // handy stopwatch
  Stopwatch sw(g_cart_grid);
  // console logger on proc_id 0
  Logger logger(0,verbosity::basic_progress,std::cout);
  Logger all_logger(0,0,std::cout);

  mpi_init_xchange_contraction(2);
  mpi_init_xchange_eo_spinor ();

  CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      init_g_gauge_field(),
      "[correlators] Error initialising gauge field!",
      true,
      CVC_EXIT_GAUGE_INIT_FAILURE);

  CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      gauge_field_eq_gauge_field_ti_phase ( &gauge_field_with_phase, g_gauge_field, co_phase_up ),
      "[correlators] Error from gauge_field_eq_gauge_field_ti_phase!",
      true,
      CVC_EXIT_UTIL_FUNCTION_FAILURE);

  ParallelMT19937_64 rng( (unsigned long long)(g_seed^Nconf) );
  std::shared_ptr<std::vector<double>> ranspinor = std::make_shared< std::vector<double> >( _GSI(VOLUME) );

  io_proc = get_io_proc ();
  if( io_proc < 0 ) {
    PRINT_STATUS(io_proc, "[correlators] Error, io proc must be ge 0!");
    EXIT(14);
  }
  fprintf(stdout, "# [correlators] proc%.4d has io proc id %d\n", g_cart_id, io_proc );

  all_logger << "# [correlators] Number of source locations: " << g_source_location_number << std::endl;

  for( int isource_location = 0; isource_location < g_source_location_number; isource_location++ ) {

    /***************************************************************************
     * local source timeslice and source process ids
     ***************************************************************************/

    int local_src_ts = -1;
    int src_proc_id  = -1;
    int src_ts       = ( g_source_coords_list[isource_location][0] +  T_global ) %  T_global;

    exitstatus = get_timeslice_source_info ( src_ts, &local_src_ts, &src_proc_id );
    if( exitstatus != 0 ) {
      fprintf(stderr, "[correlators] Error from get_timeslice_source_info status was %d %s %d\n", exitstatus, __FILE__, __LINE__);
      EXIT(123);
    }
   
    // generate a new z2 volume source for each time slice (when the inversions are done,
    // this will be reduced down to a time slice source) 
    rng.gen_z2(ranspinor->data(), 24);

    // the correlation functions will be written to this h5 file
    char corr_h5_filename[100];
    snprintf(corr_h5_filename, 100, "corr.%04d.t%d.h5", Nconf, src_ts);

    // struct for output purposes
    OutputCollection odefs;
    odefs.corr_h5_filename = corr_h5_filename;
    
    // struct for data storage purposes (also temporary)
    DataCollection data;

    // struct for holding meta-information about propagators and correlators
    // as well input data (such as the random vector and source time slice)
    MetaCollection metas;
    metas.gauge_field_with_phases = gauge_field_with_phase;
    metas.ranspinor = ranspinor;
    metas.src_ts = src_ts;

    // now parse the definitions file and build the dependency graphs for basic propagators
    // and the correlation functions
    logger << "# [correlators] Parsing " << core.get_cmd_options()["definitions_yaml"].as<std::string>() <<
      std::endl;
    YAML::Node input_node = YAML::LoadFile( core.get_cmd_options()["definitions_yaml"].as<std::string>() ); 
    yaml::enter_node(input_node, 0, odefs, metas, data); 
    logger << std::endl;

    debug_printf(
        0,verbosity::memory_info,
        "# [correlators] [MEMORY_INFO] Memory required for %d basic propagator%s: %.2f GB\n",
        metas.props_meta.size(),
        metas.props_meta.size() > 1 ? "s" : "",
        (double)metas.props_meta.size()*sizeof(double)*_GSI(g_nproc*VOLUME)*1.0e-9);

    // we first count all the nodes requiered for our momentum projections
    std::vector<ComponentGraph>
      independent_mom_projections(connected_components_subgraphs(metas.phases_graph));
    size_t mom_proj_size = 0;
    for( size_t i_component = 0; i_component < independent_mom_projections.size(); ++i_component ){
      for( auto v : boost::make_iterator_range(boost::vertices(independent_mom_projections[i_component]))) {
        mom_proj_size++;
      }
    } 
    debug_printf(
        0, 0,
        "# [correlators] [JOB_SIZE_INFO] Number of nodes for momentum projection phases: %lu in %lu connected subgraphs\n",
        mom_proj_size,
        independent_mom_projections.size());

    // now we count the basic propagators and sources
    std::vector<ComponentGraph> 
      independent_srcs_and_props(connected_components_subgraphs(metas.props_graph));
    size_t srcs_and_props_size = 0;
    for( size_t i_component = 0; i_component < independent_srcs_and_props.size(); ++i_component){
      for(auto v : boost::make_iterator_range(boost::vertices(independent_srcs_and_props[i_component]))){
        srcs_and_props_size++;
      }
    }
    debug_printf(
        0, 0,
        "# [correlators] [JOB_SIZE_INFO] Number of nodes for sources and propagators: %lu in %lu connected subgraphs\n",
        srcs_and_props_size,
        independent_srcs_and_props.size());

    // and finally the number of nodes in the correlator graph   
    std::vector<ComponentGraph>
      independent_obs(connected_components_subgraphs(metas.corrs_graph));
    size_t obs_size = 0;
    for( size_t i_component = 0; i_component < independent_obs.size(); ++i_component){
      for( auto v : boost::make_iterator_range(boost::vertices(independent_obs[i_component])) ){
        obs_size++;
      }
    }
    debug_printf(
        0, 0,
        "# [correlators] [JOB_SIZE_INFO] Number of nodes for correlation functions: %lu in %lu connected subgraphs\n",
        obs_size,
        independent_obs.size());


    // we begin by generating all momentum projection phases that we require
    for( size_t i_component = 0; i_component < independent_mom_projections.size(); ++i_component ){
      logger << std::endl << "# [correlators] Generating momentum projection phases " <<
        i_component+1 << " of " << independent_mom_projections.size() << std::endl;
      for( auto v : boost::make_iterator_range(boost::vertices(independent_mom_projections[i_component]))) {
        if( !metas.phases_graph[v].resolved ){
          descend_and_resolve<DepGraph>(v, metas.phases_graph);
        }
      }
    }

    // now we generate the basic propagators which have their own dependency graph
    // which is organised by "id" (flavour) to make the inversions efficient
    for( size_t i_component = 0; i_component < independent_srcs_and_props.size(); ++i_component){
      logger << std::endl << "# [correlators] Working on source / propagator set " << 
        i_component+1 << " of " << independent_srcs_and_props.size() << std::endl;
      for(auto v : boost::make_iterator_range(boost::vertices(independent_srcs_and_props[i_component]))){
        if( !metas.props_graph[v].resolved ){
          descend_and_resolve<DepGraph>(v, metas.props_graph);
        }
      }
    }

    // now we perform all contractions, which also includes generating sequential propagators
    // if three point functions have been defined and performing covariant shifts
    // if any covariantly displaced operators are considered
    for( size_t i_component = 0; i_component < independent_obs.size(); ++i_component){
      logger << std::endl << "# [correlators] Working on object set " << i_component+1 << 
        " of " << independent_obs.size() << std::endl;
      for( auto v : boost::make_iterator_range(boost::vertices(independent_obs[i_component])) ){
        if( !metas.corrs_graph[v].resolved )
          descend_and_resolve<DepGraph>(v, metas.corrs_graph);
      }
      sw.reset();

      debug_printf(0, 0, "# [correlators] [MEMORY INFO]: Before clearing, props_data.size(): %d seq_props_data.size(): %d,"
         " cov_displ_props_data.size(): %d, corrs_data.size() %d\n",
         data.props_data.size(), data.seq_props_data.size(), 
         data.cov_displ_props_data.size(), odefs.corrs_data.size() );
      
      // the sequential propagators that were generated for this connected set are
      // no longer required
      data.seq_props_data.clear();

      // the covariantly shifted propagators that were generated for this connected set are
      // no longer required
      data.cov_displ_props_data.clear();

      // write the correlators in this group to file and clear the temporary storage
      h5::write_correlators(corr_h5_filename, odefs.corrs_data);
      sw.elapsed_print("Correlator I/O");
#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
      // now we can clear the temporary correlator storage
      odefs.corrs_data.clear();

      debug_printf(0, 0, "# [correlators] [MEMORY INFO]: After clearing, props_data.size(): %d seq_props_data.size(): %d,"
         " cov_displ_props_data.size(): %d, corrs_data.size() %d\n",
         data.props_data.size(), data.seq_props_data.size(), 
         data.cov_displ_props_data.size(), odefs.corrs_data.size() );

#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
    }
  } // loop over source locations

	free(gauge_field_with_phase);
  mpi_fini_xchange_contraction();
  mpi_fini_xchange_eo_spinor ();

  return 0;
}
