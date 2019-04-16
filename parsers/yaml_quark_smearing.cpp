#include "global.h"
#include "Logger.hpp"
#include "constants.hpp"
#include "algorithms.hpp"
#include "types.h"
#include "meta_types.hpp"
#include "parsers/yaml_utils.hpp"

#include <yaml-cpp/yaml.h>
#include <vector>
#include <map>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <string>

namespace cvc {
namespace yaml {

  void quark_smearing(YAML::Node const & node,
                      std::map< std::string, quark_smearing_meta_t > & quark_smearing_metas)
  {
    ::cvc::Logger logger(0, verbosity::input_relay, std::cout);
    {
      for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
        logger << "\n" << it->first << ": " << it->second;
      }
      logger << std::endl;
    }

    validate_nodetype(node, YAML::NodeType::Map, "QuarkSmearing" );

    // to begin with, we require just two node types
    static const std::vector<std::string> required_nodes{
      "id", "type" };
    // the existence of which we check
    check_missing_nodes(node, required_nodes,
        "quark_smearing", "QuarkSmearing");
    // and for which we also require that they are scalar
    std::vector<std::string> scalar_nodes{"id", "type"};
    for( auto const & name : scalar_nodes ){
      validate_nodetype(node[name], YAML::NodeType::Scalar, name);
    }
    // now we determine the smearing type and adjust the requirements as appropriate
    QuarkSmearingType_t type;
    // for which we first need to transform to lowercase
    std::string type_string = node["type"].as<std::string>();
    std::transform(type_string.begin(), type_string.end(), type_string.begin(), ::tolower);
    if( type_string == "jacobi" ){
      type = SMEAR_JACOBI;
      scalar_nodes.push_back("kappa");
      scalar_nodes.push_back("n_iter");
    } else if ( type_string == "mom_jacobi" ){
      type = SMEAR_MOM_JACOBI;
      scalar_nodes.push_back("mom_scale_factor");
    } else {
      char msg[200];
      snprintf(msg, 200,
          "[yaml::quark_smearing] '%s' is not a valid quark field smearing type!\n",
          node["type"].as<std::string>().c_str() ); 
      EXIT_WITH_MSG(CVC_EXIT_INVALID_INPUT, msg);
    }
    // and make sure that all required node types exist and have the right type
    check_missing_nodes(node, scalar_nodes,
        "quark_smearing", "QuarkSmearing");
    for( auto const & name : scalar_nodes ){
      validate_nodetype(node[name], YAML::NodeType::Scalar, name);
    }

    // now we're ready to build the meta information for this smearing instance
    std::string id = node["id"].as<std::string>();
    quark_smearing_metas[id] = quark_smearing_meta_t();

    quark_smearing_metas[id].type = type;
    if( type == SMEAR_JACOBI || type == SMEAR_MOM_JACOBI ){
      quark_smearing_metas[id].n_iter = node["n_iter"].as<unsigned int>();
      quark_smearing_metas[id].kappa = node["kappa"].as<double>();
    }
    if( type == SMEAR_MOM_JACOBI ){
      quark_smearing_metas[id].mom_scale_factor = node["mom_scale_factor"].as<double>();
    }
  }

} // namespace(yaml)
} // namespace(cvc)
