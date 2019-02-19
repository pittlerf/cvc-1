#include "global.h"
#include "types.h"
#include "meta_types.hpp"
#include "Logger.hpp"
#include "constants.hpp"

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
                                     const bool verbose,
                                     mom_lists_t & mom_lists,
                                     std::map< std::string, cvc::stoch_prop_meta_t > & props_meta) {
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  cvc::Logger logger(0, verbosity::input_relay, std::cout);

  if( node.Type() != YAML::NodeType::Map ){
    throw( std::invalid_argument("in construct_time_slice_propagator, 'node' must be of type YAML::NodeType::Map\n") );
  }
  if( !(node["id"]) || !(node["solver_id"]) || !(node["solver_driver"]) ||
      !(node["g_src"]) || !(node["P_src"]) ){
    throw( std::invalid_argument("for TimeSlicePropagator, the properties 'id', 'solver_id',"
                                 " 'solver_driver', 'g_src' and 'P_src' must be defined!\n") );
  }

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
    if( node["g_src"].Type() == YAML::NodeType::Scalar ){
      
      int g_src;
      if( node["g_src"].as<std::string>().find("spin_dilution") != std::string::npos ){
        g_src = -1; 
      } else {
        g_src = node["g_src"].as<int>();
      }

      cvc::stoch_prop_meta_t stoch_prop(mom, g_src, node["id"].as<std::string>());
      props_meta[stoch_prop.key()] = stoch_prop;
      {
        logger << "\nAdded stoch_prop_meta_t: " << stoch_prop.key();
      }
    } else {
      for(size_t i = 0; i < node["g_src"].size(); ++i){
        cvc::stoch_prop_meta_t stoch_prop(mom, node["g_src"][i].as<int>(), 
                                          node["id"].as<std::string>());
        props_meta[stoch_prop.key()] = stoch_prop;
        {
          logger << "\nAdded stoch_prop_meta_t: " << stoch_prop.key();
        }
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

