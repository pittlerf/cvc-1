#pragma once

#include "types.h"
#include "DependencyGraph.hpp"

#include <map>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace cvc {
namespace yaml {

/**
* @brief Parse a 'TimeSlicePropagator' object in the object definitions file
*
* @param node Node of the 'TimeSlicePropagator'
* @param verbose Trigger verbose output
* @param mom_lists Reference to map of vector of momentum triplets (input)
* @param props_meta Reference to map of meta info for the propagators (output)
*/
void construct_time_slice_propagator(const YAML::Node &node, 
                                     const bool verbose,
                                     const int src_ts,
                                     mom_lists_t & mom_lists,
                                     std::map< std::string, cvc::ts_stoch_src_meta_t > & srcs_meta,
                                     std::map< std::string, cvc::stoch_prop_meta_t > & props_meta,
                                     DepGraph & props_graph);

} //namespace(yaml)
} //namespace(cvc)

