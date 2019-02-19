#pragma once

#include "types.h"
#include "DependencyGraph.hpp"

#include <map>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace cvc {
namespace yaml {

void construct_time_slice_propagator(const YAML::Node &node, 
                                     const bool verbose,
                                     const int src_ts,
                                     mom_lists_t & mom_lists,
                                     std::map< std::string, ts_stoch_src_meta_t > & srcs_meta,
                                     std::map< std::string, stoch_prop_meta_t > & props_meta,
                                     DepGraph & props_graph,
                                     const std::vector<double> & ranspinor,
                                     std::vector<double> & src); 

} //namespace(yaml)
} //namespace(cvc)

