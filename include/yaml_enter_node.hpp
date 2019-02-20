#pragma once

#include "meta_types.hpp"

#include <yaml-cpp/yaml.h>
#include <string>

namespace cvc {
namespace yaml {

  /**
   * @brief Driver for parsing a YAML node hierarchy for CVC object definitions
   *
   * @param node Starting node.
   * @param depth Current depth.
   * @param metas Collection of maps of meta types.
   */
  void enter_node(const YAML::Node &node,
                  const unsigned int depth,
                  const OutputDefinitions & odefs, 
                  MetaCollection & metas,
                  DataCollection & data);

} //namespace(yaml)
} //namespace(cvc)
