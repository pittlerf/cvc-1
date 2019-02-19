#include "global.h"

#include "DependencyGraph.hpp"
#include "meta_types.hpp"
#include "yaml_utils.hpp"
#include "constants.hpp"
#include "Logger.hpp"

#include <yaml-cpp/yaml.h>
#include <exception>
#include <stdexcept>

namespace cvc {
namespace yaml {

void construct_oet_meson_two_point_function(const YAML::Node &node, 
                                            const bool verbose,
                                            mom_lists_t & mom_lists,
                                            const std::map< std::string, cvc::stoch_prop_meta_t > & props_meta,
                                            DepGraph & g)
{
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  cvc::Logger logger(0, verbosity::input_relay, std::cout);

  { 
    for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
      logger << "\n  " << it->first << ": " << it->second;
    }
    logger << std::endl;
  }

  validate_nodetype(node, 
                    std::vector<YAML::NodeType::value>{YAML::NodeType::Map},
                    "OetMesonTwoPointFunction" );
  
  if( !node["id"] || !node["fwd_flav"] || !node["bwd_flav"] || !node["gi"] || !node["gf"] || !node["gb"] ||
      !node["Pi"] || !node["Pf"] || !node["momentum_conservation"] || !node["dirac_join"] || !node["flav_join"] ){
    throw( std::invalid_argument("for 'OetMesonTwoPointFunction', the properties 'id', 'fwd_flav', "
                                 "'bwd_flav', 'gi', 'gf', 'gb', 'Pi', 'Pf', 'momentum_conservation', "
                                 "'dirac_join' and 'flav_join' must be defined!\n") );
  }
  
  for( const std::string name : {"id", "fwd_flav", "bwd_flav", "momentum_conservation", 
                                 "dirac_join", "flav_join", "Pi", "Pf"} ){
    validate_nodetype(node[name],
                      std::vector<YAML::NodeType::value>{YAML::NodeType::Scalar}, 
                      name);
    if( name == "momentum_conservation" ){
      validate_bool(node[name].as<std::string>(), name);
    } else if ( name == "dirac_join" || name == "flav_join" ){
      validate_join(node[name].as<std::string>(), name);
    }
  }
  for( const auto & name : {"gi", "gf", "gb"} ){
    validate_nodetype(node[name], 
                      std::vector<YAML::NodeType::value>{YAML::NodeType::Sequence}, 
                      name);
  }
  for( const auto & name : {"Pi", "Pf"} ){
    validate_mom_lists_key(mom_lists, node[name].as<std::string>(), name,
                           node["id"].as<std::string>());
  }


  for( const auto & pi : mom_lists[ node["Pi"].as<std::string>() ] ){
    for( const auto & pf : mom_lists[ node["Pf"].as<std::string>() ] ){
      mom_t mom_xchange{ -(pi.x+pf.x), -(pi.y+pf.y), -(pi.z+pf.z) };
      if( node["momentum_conservation"].as<std::string>() == "true" ){
        if( mom_xchange.x != 0 || mom_xchange.y != 0 || mom_xchange.z != 0 ){
          continue;
        }
      }

      {
        std::vector<int> pivec{ pi.x, pi.y, pi.z };
        std::vector<int> pfvec{ pf.x, pf.y, pf.z };
        logger << "Momentum: (" << pivec[0];
        for( size_t i_pi = 1; i_pi < pivec.size(); ++i_pi ){
          if( i_pi < 3 ) logger << ",";
          logger << pivec[i_pi];
        }
        logger << " ; " << pfvec[0];
        for( size_t i_pf = 1; i_pf < pfvec.size(); ++i_pf ){
          if( i_pf < 3 ) logger << ",";
          logger << pfvec[i_pf];
        }
        logger << ")" << std::endl;
      }


      YAML::Node gb = node["gb"];
      YAML::Node gi = node["gi"];
      for( size_t i_gi = 0; i_gi < node["gi"].size(); ++i_gi ){
        YAML::Node gf;
        if( node["dirac_join"].as<std::string>() == "inner" ){
          gf = gi;
        } else {
          gf = node["gf"];
        }
        for( size_t i_gf = 0; i_gf < gf.size(); ++i_gf ){
          for( size_t i_gb = 0; i_gb < gb.size(); ++i_gb ){
            {
              logger << "Dirac: (" << gf[i_gf].as<int>() << 
                "," << gi[i_gi].as<int>() << ")  ";
              logger << "BwdDirac: " << gb[i_gb].as<int>() << std::endl;
            }

            std::string fwd_prop_key(stoch_prop_meta_t::key(pi,
                                                            gi[i_gi].as<int>(),
                                                            node["fwd_flav"].as<std::string>()));
            std::string bwd_prop_key(stoch_prop_meta_t::key(zero_mom,
                                                            gb[i_gb].as<int>(),
                                                            node["bwd_flav"].as<std::string>()));

            validate_prop_key(props_meta, fwd_prop_key, "fwd_flav", node["id"].as<std::string>());
            validate_prop_key(props_meta, bwd_prop_key, "bwd_flav", node["id"].as<std::string>());

            char corrkey[500];
            snprintf(corrkey, 500,
                     "%s+-g-%s-g/gf%d/pfx%dpfy%dpfz%d/gi%d/pix%dpiy%dpz%d",
                     node["bwd_flav"].as<std::string>().c_str(),
                     node["fwd_flav"].as<std::string>().c_str(),
                     gf[i_gf].as<int>(),
                     pf.x, pf.y, pf.z,
                     gi[i_gi].as<int>(),
                     pi.x, pi.y, pi.z);

            Vertex corrvertex = add_vertex(corrkey, g);
            g[corrvertex].fulfill.reset( new 
                CorrFulfill(fwd_prop_key, bwd_prop_key, pf, gf[i_gf].as<int>()) ); 

          } // gb
        } // gf
      } // gi
    } // pf
  } // pi
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
}

} //namespace(yaml)
} //namespace(cvc)

