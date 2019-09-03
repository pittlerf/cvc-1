#include "cvc_global.h"
#include "types.h"
#include "meta_types.hpp"
#include "Logger.hpp"
#include "constants.hpp"
#include "DependencyGraph.hpp"
#include "SourceCreators.hpp"
#include "parsers/yaml_time_slice_propagator.hpp"
#include "parsers/yaml_utils.hpp"
#include "exceptions.hpp"

#include <yaml-cpp/yaml.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>

namespace cvc {
namespace yaml {

void dirac_operator(const YAML::Node &node, 
                    std::map< std::string, dirac_op_meta_t > & dirac_op_metas)
{
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  ::cvc::Logger logger(0, verbosity::input_relay, std::cout);
  ::cvc::Logger all_logger(0, 0, std::cout);

  validate_nodetype(node, YAML::NodeType::Map, "DiracOperator");

  static const std::vector<std::string> required_nodes {
    "id", "solver_id", "solver_driver" };
  static const std::vector<std::string> scalar_nodes {
    "id", "solver_id", "solver_driver" };
  
  check_missing_nodes(node, required_nodes, 
      "cvc::yaml::dirac_operator", "DiracOperator");
  
  for( auto const name : required_nodes ){
    validate_nodetype(node[name], YAML::NodeType::Scalar, name);
  }

  {
    for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
      logger << "\n  " << it->first << ": " << it->second;
    }
  }
 
  if( dirac_op_metas.count( node["id"].as<std::string>() ) == 0 ){
    dirac_op_metas[ node["id"].as<std::string>() ] =
      dirac_op_meta_t(node["solver_driver"].as<std::string>(),
                      node["id"].as<std::string>(),
                      node["solver_id"].as<int>());
  }

#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
}

} // namespace(yaml)
} // namespace(cvc)

