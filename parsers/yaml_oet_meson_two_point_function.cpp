#include "DependencyGraph.hpp"

#include <yaml-cpp/yaml.h>
#include <exception>
#include <stdexcept>

#include "meta_types.hpp"
#include "yaml_utils.hpp"

namespace cvc {
namespace yaml {

void construct_oet_meson_two_point_function(const YAML::Node &node, 
                                            const bool verbose,
                                            const mom_lists_t & mom_lists,
                                            const std::map< std::string, cvc::stoch_prop_meta_t > & props_meta,
                                            DepGraph & g)
{
  validate_nodetype(node, 
                    std::vector<YAML::NodeType::value>{YAML::NodeType::Map},
                    "OetMesonTwoPointFunction" );
  
  if( !node["id"] || !node["fwd_flav"] || !node["bwd_flav"] || !node["g_src"] || !node["g_snk"] || !node["P_src"] ||
      !node["P_snk"] || !node["momentum_conservation"] || !node["dirac_join"] || !node["flav_join"] ){
    throw( std::invalid_argument("for 'OetMesonTwoPointFunction', the properties 'id', 'fwd_flav', "
                                 "'bwd_flav', 'g_src', 'g_snk', 'P_src', 'P_snk', 'momentum_conservation', "
                                 "'dirac_join' and 'flav_join' must be defined!\n") );
  }
  
  for( const std::string name : {"id", "fwd_flav", "bwd_flav", "momentum_conservation", 
                                 "dirac_join", "flav_join", "P_src", "P_snk"} ){
    validate_nodetype(node[name],
                      std::vector<YAML::NodeType::value>{YAML::NodeType::Scalar}, 
                      name);
    if( name == "momentum_conservation" ){
      validate_bool(node[name].as<std::string>(), name);
    } else if ( name == "dirac_join" || name == "flav_join" ){
      validate_join(node[name].as<std::string>(), name);
    }
  }
  for( const auto & name : {"g_src", "g_snk"} ){
    validate_nodetype(node[name], 
                      std::vector<YAML::NodeType::value>{YAML::NodeType::Scalar, 
                                                         YAML::NodeType::Sequence}, 
                      name);
  }
  for( const auto & name : {"P_src", "P_snk"} ){
    validate_mom_lists_key(mom_lists, node[name].as<std::string>(), name,
                           node["id"].as<std::string>());
  }

  if(verbose){
    for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
      std::cout << "\n  " << it->first << ": " << it->second;
    }
  }
}

} //namespace(yaml)
} //namespace(cvc)

