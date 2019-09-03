#include "cvc_global.h"
#include "cvc_complex.h"

#include "DependencyResolving.hpp"
#include "DependencyGraph.hpp"
#include "meta_types.hpp"
#include "parsers/yaml_utils.hpp"
#include "constants.hpp"
#include "Logger.hpp"
#include "types.h"

#include <yaml-cpp/yaml.h>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>

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
                                  const std::vector<double> & ranspinor)
{
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  ::cvc::Logger logger(0, verbosity::input_relay, std::cout);

  { 
    for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
      logger << "\n  " << it->first << ": " << it->second;
    }
    logger << std::endl;
  }

  validate_nodetype(node, YAML::NodeType::Map, "OetMesonTwoPointFunction" );

  static const std::vector<std::string> required_nodes{
    "id", "fwd_flav", "bwd_flav", "gi", "gf", "gb", "Pi", "Pf",
    "momentum_conservation", "dirac_join", "smearing_join", "smearing" };
  static const std::vector<std::string> scalar_nodes{
    "id", "fwd_flav", "bwd_flav", "Pi", "Pf",
    "momentum_conservation", "dirac_join", "smearing_join" };
  static const std::vector<std::string> sequence_nodes{
    "gi", "gf", "gb", "smearing" };

  check_missing_nodes(node, required_nodes, 
      "cvc::yaml::oet_meson_two_point_function", "OetMesonTwoPointFunction");
 
  for( auto const & name : scalar_nodes ){ 
    validate_nodetype(node[name], YAML::NodeType::Scalar, name);
    if( name == "momentum_conservation" ){
      validate_bool(node[name].as<std::string>(), name);
    } else if ( name == "dirac_join" || name == "smearing_join" ){
      validate_join(node[name].as<std::string>(), name);
    }
  }
  for( const auto & name : sequence_nodes ){
    validate_nodetype(node[name], YAML::NodeType::Sequence, name);
  }
  for( const auto & name : {"Pi", "Pf"} ){
    validate_mom_lists_key(mom_lists, node[name].as<std::string>(), name,
                           node["id"].as<std::string>());
  }
  for( const auto & name : {"fwd_flav", "bwd_flav"} ){
    validate_dirac_op_key(dirac_ops_meta, node[name].as<std::string>(), name,
                          node["id"].as<std::string>() );
  }

  const std::string fwd_flav = node["fwd_flav"].as<std::string>();
  const std::string bwd_flav = node["bwd_flav"].as<std::string>();

  // loop over the momentum and gamma combiations to construct all instances of this
  // correlation function
  for( const auto & pi : mom_lists[ node["Pi"].as<std::string>() ] ){
    char phase_string[100];
    snprintf(phase_string, 100, "px%dpy%dpz%d", pi.x, pi.y, pi.z);
    const std::string pi_key(phase_string);
    Vertex phase_vertex = boost::add_vertex(std::string(phase_string), phases_graph);
    phases_graph[phase_vertex].resolve.reset( new MomentumPhaseResolve(
          std::string(phase_string),
          phases_data,
          pi) );

    for( const auto & pf : mom_lists[ node["Pf"].as<std::string>() ] ){
      // this is a two point function, so there's no momentum exchange, but we
      // call the "relative" momentum "mom_exchange"
      mom_t mom_xchange{ -(pi.x+pf.x), -(pi.y+pf.y), -(pi.z+pf.z) };
      if( node["momentum_conservation"].as<std::string>() == "true" ){
        if( mom_xchange.x != 0 || mom_xchange.y != 0 || mom_xchange.z != 0 ){
          continue;
        }
      }

      {
        std::vector<int> pivec{ pi.x, pi.y, pi.z };
        std::vector<int> pfvec{ pf.x, pf.y, pf.z };
        logger << "# [yaml::oet_meson_two_point_function] Momentum: (" << pivec[0];
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
      
      snprintf(phase_string, 100, "px%dpy%dpz%d", pf.x, pf.y, pf.z);
      Vertex phase_vertex = boost::add_vertex(std::string(phase_string), phases_graph);
      phases_graph[phase_vertex].resolve.reset( new MomentumPhaseResolve(
            std::string(phase_string),
            phases_data,
            pf) );

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
              logger << "# [yaml::oet_meson_two_point_function] Dirac: (" << gf[i_gf].as<int>() << 
                "," << gi[i_gi].as<int>() << ")  ";
              logger << "BwdDirac: " << gb[i_gb].as<int>() << std::endl;
            }

            stoch_prop_meta_t fwd_prop_meta(pi,
                                            gi[i_gi].as<int>(),
                                            src_ts,
                                            node["fwd_flav"].as<std::string>());
            const std::string fwd_prop_key = fwd_prop_meta.key();
            if( props_meta.count(fwd_prop_key) == 0 ){
              ts_stoch_src_meta_t fwd_prop_src_meta(pi, gi[i_gi].as<int>(), src_ts);
              props_meta[fwd_prop_key] = fwd_prop_meta;
              Vertex fwd_prop_vertex = boost::add_vertex(fwd_prop_key, props_graph);
              props_graph[fwd_prop_vertex].resolve.reset(
                  new PropResolve(
                    fwd_prop_key,
                    dirac_ops_meta[fwd_flav].solver_id,
                    props_data,
                    new CreateGammaTimeSliceSource(
                      src_ts,
                      gi[i_gi].as<int>(),
                      pi,
                      fwd_prop_src_meta.key(),
                      phases_data,
                      pi_key,
                      ranspinor)
                    )
                  );

              // all simple propagators of a given flavour (id) are collected in one
              // group, that way, the MG setup, if any, can be used optimally
              Vertex fwd_flav_vertex = boost::add_vertex(fwd_flav, props_graph);
              ::cvc::add_unique_edge(fwd_prop_vertex, fwd_flav_vertex, props_graph);
            }

            stoch_prop_meta_t bwd_prop_meta(zero_mom,
                                            gb[i_gb].as<int>(),
                                            src_ts,
                                            node["bwd_flav"].as<std::string>());
            const std::string bwd_prop_key = bwd_prop_meta.key();
            if( props_meta.count(bwd_prop_key) == 0 ){
              ts_stoch_src_meta_t bwd_prop_src_meta(zero_mom, gb[i_gb].as<int>(), src_ts);
              props_meta[bwd_prop_key] = bwd_prop_meta;
              Vertex bwd_prop_vertex = boost::add_vertex(bwd_prop_key, props_graph);
              props_graph[bwd_prop_vertex].resolve.reset(
                  new PropResolve(
                    bwd_prop_key,
                    dirac_ops_meta[bwd_flav].solver_id,
                    props_data,
                    new CreateGammaTimeSliceSource(
                      src_ts,
                      gb[i_gb].as<int>(),
                      zero_mom,
                      bwd_prop_src_meta.key(),
                      phases_data,
                      "px0pz0py0",
                      ranspinor)
                    )
                  );

              Vertex bwd_flav_vertex = boost::add_vertex(bwd_flav, props_graph);
              ::cvc::add_unique_edge(bwd_prop_vertex, bwd_flav_vertex, props_graph);
            }

            char corrtype[100];
            snprintf(corrtype, 100, "%s+-g-%s-g",
                     bwd_flav.c_str(),
                     fwd_flav.c_str());

            char subpath[100];
            std::list<std::string> path_list;
            path_list.push_back(corrtype);
            snprintf(subpath, 100, "t%d", src_ts);
            path_list.push_back(subpath);
            snprintf(subpath, 100, "gf%d", gf[i_gf].as<int>());
            path_list.push_back(subpath);
            snprintf(subpath, 100, "pfx%dpfy%dpfz%d", pf.x, pf.y, pf.z);
            path_list.push_back(subpath);
            snprintf(subpath, 100, "gi%d", gi[i_gi].as<int>());
            path_list.push_back(subpath);
            snprintf(subpath, 100, "pix%dpiy%dpiz%d", pi.x, pi.y, pi.z);
            path_list.push_back(subpath);

            Vertex corrvertex = boost::add_vertex(h5::path_list_to_key(path_list), corrs_graph);
            // for the two point function, both propagators are stored in the same
            // map
            // when the contraction is performed, the bwd_prop is daggered and
            // a gamma5 is explicitly added to account for gamma5 hermiticity
            // -> contract_twopoint_gamma5_gamma_snk_only_snk_momentum
            // -> we need the full gamma structure at the sink
            corrs_graph[corrvertex].resolve.reset( new 
                CorrResolve(fwd_prop_key,
                            bwd_prop_key,
                            pf, 
                            gf[i_gf].as<int>(),
                            path_list,
                            props_data,
                            props_data,
                            corrs_data,
                            phases_data,
                            ::cvc::complex{1.0, 0.0} ) );
            

            // adding a vertex for the correlator type allows correlators to be
            // processed in bunches for efficient I/O
            Vertex corrtypevertex = boost::add_vertex(std::string(corrtype), corrs_graph);
            ::cvc::add_unique_edge(corrvertex, corrtypevertex, corrs_graph);

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

