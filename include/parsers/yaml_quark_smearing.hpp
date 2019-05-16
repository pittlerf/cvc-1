#pragma once

#include <yaml-cpp/yaml.h>
#include <map>
#include <string>

#include "meta_types.hpp"

namespace cvc {
namespace yaml {

void quark_smearing(YAML::Node const & node, 
                    std::map< std::string, quark_smearing_meta_t > & quark_smearing_meta );

} //namespace(yaml)
} //namespace(cvc)
