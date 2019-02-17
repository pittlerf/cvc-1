#pragma once

#include "types.h"

#include <yaml-cpp/yaml.h>
#include <vector>
#include <string>

namespace cvc {
namespace yaml {

void validate_nodetype(const YAML::Node & node,
                       const std::vector<YAML::NodeType::value> & types,
                       const std::string & object_name);

void check_mom_lists_key(const mom_lists_t & mom_lists,
                         const std::string & mom_lists_key,
                         const std::string & mom_property_name,
                         const std::string & object_name);

} //namespace(yaml)
} //namespace(cvc)

