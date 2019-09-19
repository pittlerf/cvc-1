#pragma once

#include <yaml-cpp/yaml.h>
#include <map>
#include <string>

#include "meta_types.hpp"

namespace cvc {
namespace yaml {

void gauge_smearing(YAML::Node const & node, 
                    std::map< std::string, gauge_smearing_meta_t > & gauge_smearing_meta );

} //namespace(yaml)
} //namespace(cvc)
