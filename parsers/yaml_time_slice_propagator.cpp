#include "global.h"
#include "types.h"
#include "meta_types.hpp"
#include "Logger.hpp"
#include "constants.hpp"
#include "DependencyGraph.hpp"
#include "SourceCreators.hpp"
#include "parsers/yaml_time_slice_propagator.hpp"
#include "parsers/yaml_utils.hpp"

#include <yaml-cpp/yaml.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>

namespace cvc {
namespace yaml {

void time_slice_propagator(const YAML::Node &node, 
                           const int src_ts,
                           mom_lists_t & mom_lists,
                           std::map< std::string, ts_stoch_src_meta_t > & srcs_meta,
                           std::map< std::string, stoch_prop_meta_t > & props_meta,
                           DepGraph & props_graph,
                           std::map< std::string, std::vector<double> > & props_data,
                           DepGraph & phases_graph,
                           std::map< std::string, std::vector<::cvc::complex> > & phases_data,
                           const std::vector<double> & ranspinor)
{
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  ::cvc::Logger logger(0, verbosity::input_relay, std::cout);
  ::cvc::Logger all_logger(0, 0, std::cout);

  validate_nodetype(node, YAML::NodeType::Map, "TimeSlicePropagator");

  static const std::vector<std::string> required_nodes {
    "id", "solver_id", "solver_driver", "g", "P" };
  static const std::vector<std::string> scalar_nodes {
    "id", "solver_id", "solver_driver", "P" };
  static const std::vector<std::string> sequence_nodes {
    "g" };
  
  check_missing_nodes(node, required_nodes, 
      "cvc::yaml::time_slice_propagator", "TimeSlicePropagator");
  
  for( const std::string name : {"id", "solver_id", "solver_driver", "P"} ){
    validate_nodetype(node[name], YAML::NodeType::Scalar, name);
  }
  validate_nodetype(node["g"], YAML::NodeType::Sequence, "g");

  {
    for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
      logger << "\n  " << it->first << ": " << it->second;
    }
  }
  if( !mom_lists.count( node["P"].as<std::string>() ) ){
    char msg[200];
    snprintf(msg, 200,
             "The momentum list '%s' does not seem to exist!\n",
             node["P"].as<std::string>().c_str() );
    throw( std::invalid_argument(msg) );
  }
  const std::string momlist_key = node["P"].as<std::string>();
  for( auto & mom : mom_lists[ momlist_key ] ){
    char phase_string[100];
    snprintf(phase_string, 100, "px%dpy%dpz%d", mom.x, mom.y, mom.z);
    const std::string phase_key(phase_string);
    Vertex phase_vertex = boost::add_vertex(phase_key, phases_graph);
    phases_graph[phase_vertex].resolve.reset( new MomentumPhaseResolve(
          phase_key,
          phases_data,
          mom) );

    for(size_t i = 0; i < node["g"].size(); ++i){
      int g = node["g"][i].as<int>();

      ts_stoch_src_meta_t src_meta(mom, g, src_ts);
      srcs_meta[src_meta.key()] = src_meta;

      // vertex for the source. Because of the way that we deal with the 
      // vertex names, these vertices are unique and multiple insertions of the same
      // vertex will leave the graph unmodified
      // Vertex src_vertex = boost::add_vertex(src_meta.key(),  props_graph);
      // props_graph[src_vertex].resolve.reset( 
      //     new TimeSliceSourceResolve(src_ts, g, mom, src_meta.key(), ranspinor, src) ); 

      ::cvc::stoch_prop_meta_t prop_meta(mom, 
                                         node["g"][i].as<int>(),
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
              node["g"][i].as<int>(),
              mom,
              src_meta.key(),
              phases_data,
              phase_key,
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

