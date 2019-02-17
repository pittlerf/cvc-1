#pragma once

#include <vector>
#include <map>
#include <yaml-cpp/yaml.h>

#include "types.h"

namespace cvc {
namespace yaml {

std::vector< mom_t > parse_momentum_list(const YAML::Node & node);

std::vector< mom_t > psq_to_momentum_list(const YAML::Node & node);

void construct_momentum_list(const YAML::Node & node, const bool verbose, mom_lists_t & mom_lists );

} //namespace(yaml)
} //namespace(cvc)
