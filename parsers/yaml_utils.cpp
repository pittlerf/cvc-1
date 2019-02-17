#include <yaml-cpp/yaml.h>

#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "types.h"

namespace cvc {
namespace yaml {

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
                 "in 'validate_nodetype', checking for type %d is not foreseen\n", type);
        throw( std::invalid_argument(msg) );
        break;
    }
  }

  if( std::find(types.begin(), types.end(), node.Type()) == types.end() ){
    char msg[200];
    snprintf(msg, 200,
             "for '%s', the node must be of '%s' type!\n",
             object_name.c_str(), typestring.c_str() );
    throw( std::invalid_argument(msg) );
  }
}

void validate_bool(const std::string & str, const std::string & name)
{
  if( !( str == "true" || str == "false" ) ){
    char msg[200];
    snprintf(msg, 200, "'%s' must be either 'true' or 'false'!\n", name.c_str());
    throw( std::invalid_argument(msg) );
  }
}

void validate_join(const std::string & str, const std::string & name)
{
  if( !(str == "inner" || str == "outer" ) ){
    char msg[200];
    snprintf(msg, 200, "'%s' must be either 'inner' or 'outer'!\n", name.c_str());
    throw( std::invalid_argument(msg) );
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
             "In the definition of momentum '%s' of '%s', the momentum list '%s' does not seem"
             "to exist!\n",
             mom_property_name.c_str(), object_name.c_str(), mom_lists_key.c_str());
    throw( std::invalid_argument(msg) );
  }
}

} // namespace(yaml)
} // namespace(cvc)
