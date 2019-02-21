#include "global.h"
#include "types.h"
#include "meta_types.hpp"
#include "Logger.hpp"
#include "constants.hpp"
#include "DependencyGraph.hpp"
#include "SourceCreators.hpp"
#include "yaml_time_slice_propagator.hpp"
#include "yaml_utils.hpp"

#include <yaml-cpp/yaml.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>

namespace cvc {
namespace yaml {

void construct_time_slice_propagator(const YAML::Node &node, 
                                     const int src_ts,
                                     mom_lists_t & mom_lists,
                                     std::map< std::string, ts_stoch_src_meta_t > & srcs_meta,
                                     std::map< std::string, stoch_prop_meta_t > & props_meta,
                                     DepGraph & props_graph,
                                     std::map< std::string, std::vector<double> > & props_data,
                                     const std::vector<double> & ranspinor)
{
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  ::cvc::Logger logger(0, verbosity::input_relay, std::cout);

  validate_nodetype(node, YAML::NodeType::Map, "TimeSlicePropagator");
  
  if( !(node["id"]) || !(node["solver_id"]) || !(node["solver_driver"]) ||
      !(node["g_src"]) || !(node["P_src"]) ){
    throw( std::invalid_argument("for TimeSlicePropagator, the properties 'id', 'solver_id',"
                                 " 'solver_driver', 'g_src' and 'P_src' must be defined!\n") );
  }
  for( const std::string name : {"id", "solver_id", "solver_driver", "P_src"} ){
    validate_nodetype(node[name], YAML::NodeType::Scalar, name);
  }
  validate_nodetype(node["g_src"], YAML::NodeType::Sequence, "g_src");

  {
    for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
      logger << "\n  " << it->first << ": " << it->second;
    }
  }
  if( !mom_lists.count( node["P_src"].as<std::string>() ) ){
    char msg[200];
    snprintf(msg, 200,
             "The momentum list '%s' does not seem to exist!\n",
             node["P_src"].as<std::string>().c_str() );
    throw( std::invalid_argument(msg) );
  }
  const std::string momlist_key = node["P_src"].as<std::string>();
  for( auto & mom : mom_lists[ momlist_key ] ){
    for(size_t i = 0; i < node["g_src"].size(); ++i){
      int g_src = node["g_src"][i].as<int>();

      ts_stoch_src_meta_t src_meta(mom, g_src, src_ts);
      srcs_meta[src_meta.key()] = src_meta;

      // vertex for the source. Because of the way that we deal with the 
      // vertex names, these vertices are unique and multiple insertions of the same
      // vertex will leave the graph unmodified
      // Vertex src_vertex = boost::add_vertex(src_meta.key(),  props_graph);
      // props_graph[src_vertex].resolve.reset( 
      //     new TimeSliceSourceResolve(src_ts, g_src, mom, src_meta.key(), ranspinor, src) ); 

      ::cvc::stoch_prop_meta_t prop_meta(mom, 
                                         node["g_src"][i].as<int>(),
                                         src_ts,
                                         node["id"].as<std::string>(),
                                         node["solver_driver"].as<std::string>(),
                                         node["solver_id"].as<int>());
      props_meta[prop_meta.key()] = prop_meta;

      // multiple insertions of the same propagator leave the graph unmodified
      Vertex prop_vertex = boost::add_vertex(prop_meta.key(), props_graph);
      props_graph[prop_vertex].resolve.reset( 
          new PropResolve(
            prop_meta.key(),
            node["solver_id"].as<int>(),
            props_data,
            new CreateGammaTimeSliceSource(
              src_ts,
              node["g_src"][i].as<int>(),
              mom,
              src_meta.key(),
              ranspinor)
          ) );
      
      // all propagators of a given flavour (id) will be collected in one group
      // that way, the MG setup, if any, can be used optimally
      Vertex id_vertex = boost::add_vertex(node["id"].as<std::string>(), props_graph);
      ::cvc::add_unique_edge(prop_vertex, id_vertex, props_graph);

      {
        logger << "\nAdded stoch_prop_meta_t: " << prop_meta.key();
      }
    }
  }
  { logger << std::endl; }
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
}

} // namespace(yaml)
} // namespace(cvc)

