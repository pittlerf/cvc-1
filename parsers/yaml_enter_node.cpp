#include "global.h"
#include "Logger.hpp"
#include "meta_types.hpp"
#include "yaml_parsers.hpp"
#include "constants.hpp"

#include <yaml-cpp/yaml.h>
#include <iostream>

namespace cvc {
namespace yaml {

void enter_node(const YAML::Node &node, 
                const unsigned int depth,
                MetaCollection & metas,
                const bool verbose){
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  cvc::Logger logger(0, verbosity::input_relay, std::cout);

  YAML::NodeType::value type = node.Type();
  std::string indent( 2*(size_t)depth, ' ');
  switch(type){
    case YAML::NodeType::Scalar:
      { logger << node; }
      break;
    case YAML::NodeType::Sequence:
      for(unsigned int i = 0; i < node.size(); ++i){
        if(depth <= 2){
          logger << "[ ";
        }
        const YAML::Node & subnode = node[i];
        enter_node(subnode, depth+1, metas, verbose);
        if( depth <= 2 ){
          logger << " ]";
        } else if( i < node.size()-1 ) {
          logger << ",";
        }
      }
      break;
    case YAML::NodeType::Map:
      for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
        {
          if(depth != 0){ logger << std::endl << indent; }
          logger << it->first << ": ";
        }
        
        if( it->first.as<std::string>() == "MomentumList" ){
          construct_momentum_list(it->second, verbose, metas.mom_lists );
        } else if ( it->first.as<std::string>() == "TimeSlicePropagator" ){
          construct_time_slice_propagator(it->second, verbose, metas.src_ts, metas.mom_lists, 
                                          metas.srcs_meta, metas.props_meta, metas.props_graph,
                                          *(metas.ranspinor), *(metas.stochastic_source) );
        } else if ( it->first.as<std::string>() == "OetMesonTwoPointFunction" ){
          construct_oet_meson_two_point_function(it->second, verbose, metas.mom_lists, 
                                                 metas.props_meta, metas.corrs_graph); 
        } else {
          char msg[200];
          snprintf(msg, 200,
                   "%s is not a valid Object name\n",
                   it->first.as<std::string>().c_str());
          throw( std::invalid_argument(msg) );
        }
      }
    case YAML::NodeType::Null:
      { logger << std::endl; }
      break;
    default:
      break;
  }
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
}

} // namespace(yaml)
} // namespace(cvc)
