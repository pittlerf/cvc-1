#pragma once

#include "DependencyGraph.hpp"

#include <map>
#include <string>
#include <vector>

namespace cvc {
namespace yaml {

void construct_oet_meson_two_point_function(const YAML::Node &node, 
                                            mom_lists_t & mom_lists,
                                            const std::string & output_filename,
                                            std::map< std::string, ::cvc::stoch_prop_meta_t > & props_meta,
                                            std::map< std::string, std::vector<double> > & props_data,
                                            DepGraph & g);

} //namespace(yaml)
} //namespace(cvc)

