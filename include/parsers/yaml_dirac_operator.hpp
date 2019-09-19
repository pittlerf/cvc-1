#pragma once

#include <yaml-cpp/yaml.h>
#include <map>
#include <string>

#include "meta_types.hpp"

namespace cvc {
namespace yaml {

void dirac_operator(YAML::Node const & node, 
                    std::map< std::string, dirac_op_meta_t > & dirac_op_metas );

} //namespace(yaml)
} //namespace(cvc)
