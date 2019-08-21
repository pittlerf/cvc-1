#include "cvc_global.h"
#include "Logger.hpp"
#include "constants.hpp"
#include "algorithms.hpp"
#include "types.h"
#include "exceptions.hpp"

#include <yaml-cpp/yaml.h>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace cvc {
namespace yaml {

std::vector< mom_t > parse_momentum_list(const YAML::Node & node)
{
  if( node.Type() != YAML::NodeType::Sequence ){
    throw( ::cvc::invalid_argument("'node' must be of type YAML::NodeType::Sequence",
                                   "cvc::yaml::parse_momentum_list") );
  }
  std::vector< mom_t > momenta;
  for(size_t i = 0; i < node.size(); ++i){
    mom_t mom{ node[i][0].as<int>(), node[i][1].as<int>(), node[i][2].as<int>() };
    momenta.push_back(mom);
  }
  momentum_compare_t comp;
  std::sort(momenta.begin(), momenta.end(), comp);
  return momenta;
}

std::vector< mom_t > psq_to_momentum_list(const YAML::Node & node)
{
  if( node.Type() != YAML::NodeType::Scalar ){
    throw( ::cvc::invalid_argument("'node' must be of type YAML::NodeType::Scalar",
                                   "cvc::yaml::psq_to_momentum_list") );
  }
  const int psqmax = node.as<int>();
  const int pmax = static_cast<int>(sqrt( node.as<double>() )); 
  std::vector< mom_t > momenta;
  for( int px = -pmax; px <= pmax; ++px ){
    for( int py = -pmax; py <= pmax; ++py ){
      for( int pz = -pmax; pz <= pmax; ++pz ){
        if( (px*px + py*py + pz*pz) <= psqmax ){
          mom_t mom{ px, py, pz };
          momenta.push_back(mom);
        }
      }
    }
  }
  momentum_compare_t comp;
  std::sort(momenta.begin(), momenta.end(), comp);
  return momenta; 
}

void momentum_list(const YAML::Node & node,
                   mom_lists_t & mom_lists )
{
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  cvc::Logger logger(0, verbosity::input_relay, std::cout);

  if( node.Type() != YAML::NodeType::Map ){
    throw( ::cvc::invalid_argument("'node' must be of type YAML::NodeType::Map",
                                   "cvc::yaml::momentum_list") );
  }
  if( !node["id"] || !(node["Psqmax"] || node["Plist"]) ){
    throw( ::cvc::invalid_argument("For 'MomentumList', the 'id' property and one of 'Psqmax' or 'Plist' must be defined!",
                                   "cvc::yaml::momentum_list") );
  }

  std::string id;
  std::vector<mom_t> momentum_list;
  for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
    { logger << "\n  " << it->first << ": " << it->second; }

    if( it->first.as<std::string>() == "id" ){
      id = it->second.as<std::string>();
    } else if( it->first.as<std::string>() == "Plist" ){
      momentum_list = parse_momentum_list( it->second );
    } else if( it->first.as<std::string>() == "Psqmax" ){
      momentum_list = psq_to_momentum_list( it->second );
    } else {
      char msg[200];
      snprintf(msg, 200,
               "%s is not a valid property for a MomentumList!\n",
               it->first.as<std::string>().c_str());
      throw( ::cvc::invalid_argument(msg, "cvc::yaml::momentum_list") );
    }
  }
 
  { 
    logger << std::endl << "   ";
    for( const auto & mom : momentum_list ){
      std::vector<int> mvec{mom.x, mom.y, mom.z};
      logger << "[";
      for(size_t i = 0; i < mvec.size(); ++i){
        logger << mvec[i];
        if( i < mvec.size()-1 ){
          logger << ",";
        } else {
          logger << "] ";
        }
      }
    }
    logger << std::endl;
  }
  
  mom_lists[id] = momentum_list;

#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
}

} //namespace(yaml)
} //namespace(cvc)
