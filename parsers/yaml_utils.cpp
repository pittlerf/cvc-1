#include "global.h"
#include "types.h"
#include "meta_types.hpp"
#include "parsers/yaml_utils.hpp"
#include "cvc_complex.h"
#include "Logger.hpp"
#include "exceptions.hpp"

#include <yaml-cpp/yaml.h>

#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace cvc {
namespace yaml {

void validate_nodetype(const YAML::Node & node,
                       const YAML::NodeType::value & type,
                       const std::string & object_name){
  validate_nodetype(node, std::vector<YAML::NodeType::value>{type}, object_name);
}

void validate_nodetype(const YAML::Node & node,
                       const std::vector<YAML::NodeType::value> & types,
                       const std::string & object_name)
{
  std::string typestring = "";
  for( const auto & type : types ){
    switch(type){
      case YAML::NodeType::Scalar:
        if( typestring.length() == 0 )
          typestring = "Scalar";
        else
          typestring += " or Scalar";
        break;
      case YAML::NodeType::Sequence:
        if( typestring.length() == 0 )
          typestring = "Sequence";
        else
          typestring += " or Sequence";
        break;
      case YAML::NodeType::Map:
        if( typestring.length() == 0 )
          typestring = "Map";
        else
          typestring += " or Map";
        break;
      default:
        char msg[200];
        snprintf(msg, 200,
                 "in 'validate_nodetype', checking for type %d is not foreseen, this is a bug!", type);
        throw( ::cvc::invalid_argument(msg, "validate_nodetype") );
        break;
    }
  }

  if( std::find(types.begin(), types.end(), node.Type()) == types.end() ){
    char msg[200];
    snprintf(msg, 200,
             "for '%s', the node must be of '%s' type!",
             object_name.c_str(), typestring.c_str() );
    throw( ::cvc::invalid_argument(msg, "validate_nodetype") );
  }
}

void validate_bool(const std::string & str, const std::string & name)
{
  if( !( str == "true" || str == "false" ) ){
    char msg[200];
    snprintf(msg, 200, "'%s' must be either 'true' or 'false'! Please check your YAML definitions file!", name.c_str());
    throw( ::cvc::invalid_argument(msg, "validate_bool") );
  }
}

void validate_join(const std::string & str, const std::string & name)
{
  if( !(str == "inner" || str == "outer" ) ){
    char msg[200];
    snprintf(msg, 200, "'%s' must be either 'inner' or 'outer'! Please check your YAML definitions file!", name.c_str());
    throw( ::cvc::invalid_argument(msg, "validate_join") );
  }
}

void validate_mom_lists_key(const mom_lists_t & mom_lists,
                            const std::string & mom_lists_key,
                            const std::string & mom_property_name,
                            const std::string & object_name)
{
  if( !mom_lists.count(mom_lists_key) ){
    char msg[200];
    snprintf(msg, 200,
             "For the definition of momentum '%s' of '%s', the momentum list '%s' does not seem "
             "to have been defined previously! Please check your YAML definitions file!",
             mom_property_name.c_str(), object_name.c_str(), mom_lists_key.c_str());
    throw( ::cvc::invalid_argument(msg, "validate_mom_lists_key") );
  }
}

void validate_prop_key(const std::map<std::string, stoch_prop_meta_t> & props_meta,
                       const std::string & key,
                       const std::string & quarkline_name,
                       const std::string & object_name)
{
  if( !props_meta.count(key) ){
    char msg[200];
    snprintf(msg, 200,
             "For the definition of '%s', quark line '%s', propagator '%s' does not seem "
             "to have been defined previously! Please check your YAML definitions file!",
             object_name.c_str(), quarkline_name.c_str(), key.c_str());
    throw( ::cvc::invalid_argument(msg, "validate_prop_key") );
  }
}

void check_missing_nodes(
    const YAML::Node & node,
    const std::vector<std::string> & required_nodes,
    const std::string & function_name,
    const std::string & object_name )
{
  validate_nodetype( node, YAML::NodeType::Map, object_name );

  ::cvc::Logger all_logger(0,0,std::cout);

  bool nodes_are_missing = false;
  std::vector<std::string> missing_nodes;
  for( auto const & name : required_nodes ){
    if( !node[name] ) {
      nodes_are_missing = true;
      missing_nodes.push_back(name);
    }
  }
  if( nodes_are_missing ){
    all_logger << "[" << function_name << "] The required properties ";
    for( auto const & name : missing_nodes ){
      all_logger << name << " ";
    }
    all_logger << " are missing!" << std::endl;
      
    std::stringstream msg;
    msg << "for '" << object_name << "', the properties ";
    for( auto const & name : required_nodes ){
      msg << "'" << name << "' ";
    }
    msg << " must be defined!";
    throw( ::cvc::invalid_argument(msg.str(), "check_missing_nodes") );
  }
}


} // namespace(yaml)
} // namespace(cvc)
