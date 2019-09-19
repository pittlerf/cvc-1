#pragma once

#include "DependencyGraph.hpp"
#include "meta_types.hpp"

#include <map>
#include <string>
#include <vector>

namespace cvc {
namespace yaml {

void oet_meson_two_point_function(const YAML::Node &node, 
                                  mom_lists_t & mom_lists,
                                  const int src_ts,
                                  std::map< std::string, ::cvc::stoch_prop_meta_t > & props_meta,
                                  DepGraph & props_graph,
                                  std::map< std::string, std::vector<double> > & props_data,
                                  std::map< std::string, ::cvc::H5Correlator > & corrs_data, 
                                  DepGraph & corrs_graph,
                                  std::map< std::string, std::vector<::cvc::complex> > & phases_data,
                                  DepGraph & phases_graph,
                                  std::map< std::string, ::cvc::dirac_op_meta_t > & dirac_ops_meta,
                                  const std::vector<double> & ranspinor);

} //namespace(yaml)
} //namespace(cvc)

