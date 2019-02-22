#pragma once

#include "DependencyGraph.hpp"
#include "meta_types.hpp"
#include "types.h"

#include <yaml-cpp/yaml.h>

namespace cvc {
namespace yaml {

void construct_oet_meson_three_point_function(
    const YAML::Node &node,
    mom_lists_t & mom_lists,
    const int src_ts,
    std::map< std::string, ::cvc::stoch_prop_meta_t > & props_meta,
    std::map< std::string, std::vector<double> > & props_data,
    std::map< std::string, ::cvc::H5Correlator > & corrs_data,
    std::map< std::string, std::vector<double> > & seq_props_data,
    DepGraph & g);

} // namespace(yaml)
} // naemspace(cvc)

