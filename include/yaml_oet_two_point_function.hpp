#pragma once

#include "DependencyGraph.hpp"

#include <map>
#include <string>
#include <vector>

namespace cvc {
namespace yaml {

void construct_oet_two_point_function(const YAML::Node &node, 
                                      const bool verbose,
                                      const mom_lists_t & mom_lists,
                                      const std::map< std::string, cvc::stoch_prop_meta_t > & props_meta,
                                      DepGraph & g);

} //namespace(yaml)
} //namespace(cvc)

