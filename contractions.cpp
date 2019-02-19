#define MAIN_PROGRAM
#include "global.h"
#undef MAIN_PROGRAM

#include "DependencyGraph.hpp"
#include "DependencyFulfilling.hpp"
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

  Core core(argc,argv,"test_parse_yaml");
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
  std::shared_ptr<std::vector<double>> stochastic_source = std::make_shared< std::vector<double> >(_GSI(VOLUME));

  io_proc = get_io_proc ();
  if( io_proc < 0 ) {
    PRINT_STATUS(io_proc, "[correlators] Error, io proc must be ge 0!");
    EXIT(14);
  }
  fprintf(stdout, "# [correlators] proc%.4d has io proc id %d\n", g_cart_id, io_proc );

  std::vector<int> src_time_slices{12};

  for(int src_ts : src_time_slices){
    rng.gen_z2(ranspinor->data(), 24);

    DataCollection data;
    MetaCollection metas;
    metas.ranspinor = ranspinor;
    metas.stochastic_source = stochastic_source;
    metas.src_ts = src_ts;
    YAML::Node input_node = YAML::LoadFile("definitions.yaml");
    yaml::enter_node(input_node, 0, metas, data); 
    logger << std::endl;

    debug_printf(0,verbosity::memory_info,
                 "Memory required for %d basic propagators, %.2f GB\n",
                 metas.props_meta.size(), 
                 (double)metas.props_meta.size()*sizeof(double)*_GSI(g_nproc*VOLUME)*1.0e-9);

    boost::property_map<DepGraph, std::string VertexProperties::*>::type name_map = 
      boost::get(&VertexProperties::name, metas.corrs_graph);
    
    std::vector<ComponentGraph> 
      independent_srcs_and_props(connected_components_subgraphs(metas.props_graph));
    for( size_t i_component = 0; i_component < independent_srcs_and_props.size(); ++i_component){
      logger << std::endl << "Working on source / propagator set " << i_component << std::endl;
      for(auto v : boost::make_iterator_range(boost::vertices(independent_srcs_and_props[i_component]))){
        if( !metas.props_graph[v].fulfilled ){
          descend_and_fulfill<DepGraph>(v, metas.props_graph);
        }
      }
    }

    for( auto v : boost::make_iterator_range(boost::vertices(metas.corrs_graph)) ){
      if( !metas.corrs_graph[v].fulfilled )
        descend_and_fulfill<DepGraph>(v, metas.corrs_graph);
    }
  }

	free(gauge_field_with_phase);
  mpi_fini_xchange_contraction();
  mpi_fini_xchange_eo_spinor ();

  return 0;
}
