#define MAIN_PROGRAM
#include "global.h"
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

  std::vector<int> src_time_slices{12};

  for(int src_ts : src_time_slices){
    rng.gen_z2(ranspinor->data(), 24);

    char corr_h5_filename[100];
    snprintf(corr_h5_filename, 100, "corr.%04d.t%d.h5", Nconf, src_ts);
    OutputCollection odefs;
    odefs.corr_h5_filename = corr_h5_filename; 
    
    DataCollection data;
    MetaCollection metas;
    metas.ranspinor = ranspinor;
    metas.src_ts = src_ts;

    logger << "# [correlators] Parsing " << core.get_cmd_options()["definitions_yaml"].as<std::string>() <<
      std::endl;
    YAML::Node input_node = YAML::LoadFile( core.get_cmd_options()["definitions_yaml"].as<std::string>() ); 
    yaml::enter_node(input_node, 0, odefs, metas, data); 
    logger << std::endl;

    debug_printf(0,verbosity::memory_info,
                 "# [correlators] [MEMORY_INFO] Memory required for %d basic propagator%s: %.2f GB\n",
                 metas.props_meta.size(),
                 metas.props_meta.size() > 1 ? "s" : "",
                 (double)metas.props_meta.size()*sizeof(double)*_GSI(g_nproc*VOLUME)*1.0e-9);

    boost::property_map<DepGraph, std::string VertexProperties::*>::type name_map = 
      boost::get(&VertexProperties::name, metas.corrs_graph);
    
    std::vector<ComponentGraph> 
      independent_srcs_and_props(connected_components_subgraphs(metas.props_graph));
    for( size_t i_component = 0; i_component < independent_srcs_and_props.size(); ++i_component){
      logger << std::endl << "# [correlators] Working on source / propagator set " << 
        i_component+1 << " of " << independent_srcs_and_props.size() << std::endl;
      for(auto v : boost::make_iterator_range(boost::vertices(independent_srcs_and_props[i_component]))){
        if( !metas.props_graph[v].resolved ){
          descend_and_resolve<DepGraph>(v, metas.props_graph);
        }
      }
    }

    std::vector<ComponentGraph>
      independent_obs(connected_components_subgraphs(metas.corrs_graph));
    for( size_t i_component = 0; i_component < independent_obs.size(); ++i_component){
      logger << std::endl << "# [correlators] Working on object set " << i_component+1 << 
        " of " << independent_obs.size() << std::endl;
      for( auto v : boost::make_iterator_range(boost::vertices(independent_obs[i_component])) ){
        if( !metas.corrs_graph[v].resolved )
          descend_and_resolve<DepGraph>(v, metas.corrs_graph);
      }
      sw.reset();
      // the sequential propagators that were generated for this connected set are
      // no longer required
      data.seq_props_data.clear();
      // the covariantly shifted propagators that were generated for this connected set are
      // no longer required
      data.deriv_props_data.clear();
      // write the correlators in this group to file and clear the temporary storage
      h5::write_correlators(corr_h5_filename, odefs.corrs_data);
      sw.elapsed_print("Correlator I/O");
#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
      odefs.corrs_data.clear();
    }
  }

	free(gauge_field_with_phase);
  mpi_fini_xchange_contraction();
  mpi_fini_xchange_eo_spinor ();

  return 0;
}
