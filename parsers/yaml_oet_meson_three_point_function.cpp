#include "global.h"
#include "cvc_complex.h"

#include "DependencyResolving.hpp"
#include "DependencyGraph.hpp"
#include "meta_types.hpp"
#include "yaml_utils.hpp"
#include "constants.hpp"
#include "Logger.hpp"
#include "types.h"
#include "deriv_tools.hpp"

#include <yaml-cpp/yaml.h>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>

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
    std::map< std::string, std::vector<double> > & deriv_props_data,
    DepGraph & g)
{
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  ::cvc::Logger logger(0, verbosity::input_relay, std::cout);
  ::cvc::Logger all_logger(0, 0, std::cout);

  { 
    for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
      logger << "\n  " << it->first << ": " << it->second;
    }
    logger << std::endl;
  }

  validate_nodetype(node, YAML::NodeType::Map, "OetMesonThreePointFunction" );
  
  static const std::vector<std::string> required_nodes{
    "id", "dt", "fwd_flav", "bwd_flav", "seq_flav", "gi", "gf", "gb", "gc", "Pi", "Pf",
    "momentum_exchange", "dirac_join", "Dc", "seq_solver_id", "seq_solver_driver"};

  static const std::vector<std::string> scalar_nodes{
    "id", "fwd_flav", "bwd_flav", "seq_flav", "dirac_join", "momentum_exchange",
    "Pi", "Pf", "seq_solver_id", "seq_solver_driver" };

  static const std::vector<std::string> sequence_nodes{
    "dt", "gi", "gf", "gb", "gc", "Dc" };

  check_missing_nodes(node, required_nodes,
      "construct_oet_meson_three_point_function", "OetMesonThreePointFunction");
 
  for( const std::string name : scalar_nodes ){
    validate_nodetype(node[name],
                      std::vector<YAML::NodeType::value>{YAML::NodeType::Scalar}, 
                      name);
    if( name == "momentum_exchange" ){
      validate_bool(node[name].as<std::string>(), name);
    } else if ( name == "dirac_join" ){
      validate_join(node[name].as<std::string>(), name);
    }
  }
  for( const auto & name : sequence_nodes ){
    validate_nodetype(node[name], 
                      std::vector<YAML::NodeType::value>{YAML::NodeType::Sequence}, 
                      name);
  }
  for( const auto & name : {"Pi", "Pf"} ){
    validate_mom_lists_key(mom_lists, node[name].as<std::string>(), name,
                           node["id"].as<std::string>());
  }

  // loop over the source-sink-separation as well as the momentum combinations at this level
  // further below over all gamma and derivative combinations
  // to construct all instances of this correlation function
  for( size_t i_dt = 0; i_dt < node["dt"].size(); ++i_dt ){
    const int dt = node["dt"][i_dt].as<int>();
    {
      logger << "# [construct_oet_meson_three_point_function] Source-sink-separation: " <<
        dt << std::endl;
    }
    const int seq_src_ts = (src_ts + dt + T_global) % T_global;

    for( const auto & pi : mom_lists[ node["Pi"].as<std::string>() ] ){
      for( const auto & pf : mom_lists[ node["Pf"].as<std::string>() ] ){
        mom_t mom_xchange{ -(pi.x+pf.x), -(pi.y+pf.y), -(pi.z+pf.z) };
        if( node["momentum_exchange"].as<std::string>() == "false" ){
          if( mom_xchange.x != 0 || mom_xchange.y != 0 || mom_xchange.z != 0 ){
            continue;
          }
        }

        {
          std::vector<int> pivec{ pi.x, pi.y, pi.z };
          std::vector<int> pfvec{ pf.x, pf.y, pf.z };
          std::vector<int> pcvec{ mom_xchange.x, mom_xchange.y, mom_xchange.z };
          logger << "# [construct_oet_meson_three_point_function] Momentum: (" << pivec[0];
          for( size_t i_pi = 1; i_pi < pivec.size(); ++i_pi ){
            if( i_pi < 3 ) logger << ",";
            logger << pivec[i_pi];
          }
          logger << " ; " << pcvec[0];
          for( size_t i_pc = 1; i_pc < pcvec.size(); ++i_pc ){
            if( i_pc < 3 ) logger << ",";
            logger << pcvec[i_pc];
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
        YAML::Node gc = node["gc"];
        YAML::Node Dc = node["Dc"];
        for( size_t i_gi = 0; i_gi < gi.size(); ++i_gi ){
          YAML::Node gf;
          if( node["dirac_join"].as<std::string>() == "inner" ){
            gf = gi;
          } else {
            gf = node["gf"];
          }
          for( size_t i_gf = 0; i_gf < gf.size(); ++i_gf ){
            for( size_t i_gb = 0; i_gb < gb.size(); ++i_gb ){
              for( size_t i_gc = 0; i_gc < gc.size(); ++i_gc ){
                for( size_t i_Dc = 0; i_Dc < Dc.size(); ++i_Dc ){
                  {
                    logger << "# [construct_oet_meson_three_point_function] Dirac: (" << 
                      gf[i_gf].as<int>() << "," <<
                      gc[i_gc].as<int>() <<
                      "," << gi[i_gi].as<int>() << ")  " <<
                      "BwdDirac: " << gb[i_gb].as<int>() << "  " <<
                      "Derivative: " << Dc[i_Dc].as<int>() <<
                      std::endl;
                  }

                  std::vector< std::vector<deriv_t> > deriv_chains;
                  deriv_chains = create_derivatives(0, Dc[i_Dc].as<int>(), deriv_chains);
                  
                  if( deriv_chains.size() == 0 ){
                    std::string fwd_prop_key(stoch_prop_meta_t::key(pi,
                                                                    gi[i_gi].as<int>(),
                                                                    src_ts,
                                                                    node["fwd_flav"].as<std::string>()));

                    std::string bwd_prop_key(stoch_prop_meta_t::key(zero_mom,
                                                                    gb[i_gb].as<int>(),
                                                                    src_ts,
                                                                    node["bwd_flav"].as<std::string>()));

                    validate_prop_key(props_meta, fwd_prop_key, "fwd_flav", node["id"].as<std::string>());
                    validate_prop_key(props_meta, bwd_prop_key, "bwd_flav", node["id"].as<std::string>());

                    ::cvc::seq_stoch_prop_meta_t seq_prop(mom_xchange,
                                                          gf[i_gf].as<int>(),
                                                          (src_ts + dt + T_global) % T_global,
                                                          src_ts,
                                                          node["seq_flav"].as<std::string>(),
                                                          pi,
                                                          gb[i_gb].as<int>(),
                                                          node["bwd_flav"].as<std::string>()
                                                          );

                    const std::string seq_prop_key = seq_prop.key();

                    char seq_src_string[200];
                    snprintf(seq_src_string, 200, "g%d_px%dpy%dpz%d::ts%d::%s",
                        gf[i_gf].as<int>(), pf.x, pf.y, pf.z,
                        seq_src_ts,
                        bwd_prop_key.c_str());
                    const std::string seq_src_key(seq_src_string);

                    Vertex seq_prop_vertex = boost::add_vertex(seq_prop_key, g);
                    g[seq_prop_vertex].resolve.reset( new PropResolve(
                          seq_prop_key,
                          node["seq_solver_id"].as<int>(),
                          seq_props_data,
                          new CreateSequentialGammaTimeSliceSource(
                            src_ts,
                            seq_src_ts,
                            gf[i_gf].as<int>(),
                            pf,
                            seq_src_key,
                            bwd_prop_key,
                            props_data) ) );

                    char corrtype[100];
                    snprintf(corrtype, 100, "s%s%s+-g-%s-g",
                             node["bwd_flav"].as<std::string>().c_str(),
                             node["seq_flav"].as<std::string>().c_str(),
                             node["fwd_flav"].as<std::string>().c_str());

                    char subpath[100];
                    std::list<std::string> path_list;
                    path_list.push_back(corrtype);
                    snprintf(subpath, 100, "t%d", src_ts);
                    path_list.push_back(subpath);
                    snprintf(subpath, 100, "dt%d", dt);
                    path_list.push_back(subpath);
                    snprintf(subpath, 100, "gf%d", gf[i_gf].as<int>());
                    path_list.push_back(subpath);
                    snprintf(subpath, 100, "pfx%dpfy%dpfz%d", pf.x, pf.y, pf.z);
                    path_list.push_back(subpath);
                    snprintf(subpath, 100, "gc%d", gc[i_gc].as<int>());
                    path_list.push_back(subpath); 
                    snprintf(subpath, 100, "gi%d", gi[i_gi].as<int>());
                    path_list.push_back(subpath);
                    snprintf(subpath, 100, "pix%dpiy%dpiz%d", pi.x, pi.y, pi.z);
                    path_list.push_back(subpath);

                    Vertex corrvertex = boost::add_vertex(h5::path_list_to_key(path_list), g);
                    ::cvc::add_unique_edge(corrvertex, seq_prop_vertex, g);

                    // for the three point function, the forward and daggered propagator
                    // are in separate maps
                    g[corrvertex].resolve.reset( new 
                        CorrResolve(fwd_prop_key,
                                    seq_prop_key,
                                    mom_xchange, 
                                    gc[i_gc].as<int>(),
                                    path_list,
                                    props_data,
                                    seq_props_data,
                                    corrs_data,
                                    ::cvc::complex{1.0, 0.0} ) );
                  
                  } else { // if( deriv_chains.size() > 0 )
                  
                    for( auto const & deriv_chain : deriv_chains ){
                      std::string fwd_prop_key(stoch_prop_meta_t::key(pi,
                                                                      gi[i_gi].as<int>(),
                                                                      src_ts,
                                                                      node["fwd_flav"].as<std::string>()));

                      std::string bwd_prop_key(stoch_prop_meta_t::key(zero_mom,
                                                                      gb[i_gb].as<int>(),
                                                                      src_ts,
                                                                      node["bwd_flav"].as<std::string>()));

                      validate_prop_key(props_meta, fwd_prop_key, "fwd_flav", node["id"].as<std::string>());
                      validate_prop_key(props_meta, bwd_prop_key, "bwd_flav", node["id"].as<std::string>());

                      ::cvc::seq_stoch_prop_meta_t seq_prop(mom_xchange,
                                                            gf[i_gf].as<int>(),
                                                            (src_ts + dt + T_global) % T_global,
                                                            src_ts,
                                                            node["seq_flav"].as<std::string>(),
                                                            pi,
                                                            gb[i_gb].as<int>(),
                                                            node["bwd_flav"].as<std::string>()
                                                            );

                      const std::string seq_prop_key = seq_prop.key();

                      char seq_src_string[200];
                      snprintf(seq_src_string, 200, "g%d_px%dpy%dpz%d::ts%d::%s",
                          gf[i_gf].as<int>(), pf.x, pf.y, pf.z,
                          seq_src_ts,
                          bwd_prop_key.c_str());
                      const std::string seq_src_key(seq_src_string);

                      Vertex seq_prop_vertex = boost::add_vertex(seq_prop_key, g);
                      g[seq_prop_vertex].resolve.reset( new PropResolve(
                            seq_prop_key,
                            node["seq_solver_id"].as<int>(),
                            seq_props_data,
                            new CreateSequentialGammaTimeSliceSource(
                              src_ts,
                              seq_src_ts,
                              gf[i_gf].as<int>(),
                              pf,
                              seq_src_key,
                              bwd_prop_key,
                              props_data) ) );

                      char corrtype[100];
                      snprintf(corrtype, 100, "s%s%s+-g-%s-g",
                               node["bwd_flav"].as<std::string>().c_str(),
                               node["seq_flav"].as<std::string>().c_str(),
                               node["fwd_flav"].as<std::string>().c_str());

                      char subpath[100];
                      std::list<std::string> path_list;
                      path_list.push_back(corrtype);
                      snprintf(subpath, 100, "t%d", src_ts);
                      path_list.push_back(subpath);
                      snprintf(subpath, 100, "dt%d", dt);
                      path_list.push_back(subpath);
                      snprintf(subpath, 100, "gf%d", gf[i_gf].as<int>());
                      path_list.push_back(subpath);
                      snprintf(subpath, 100, "pfx%dpfy%dpfz%d", pf.x, pf.y, pf.z);
                      path_list.push_back(subpath);
                      snprintf(subpath, 100, "gc%d", gc[i_gc].as<int>());
                      path_list.push_back(subpath);
                      for( auto const & deriv : deriv_chain ){
                        snprintf(subpath, 100, "Ddim%d:dir%d", deriv.dim, deriv.dir);
                        path_list.push_back(subpath);
                      }
                      snprintf(subpath, 100, "gi%d", gi[i_gi].as<int>());
                      path_list.push_back(subpath);
                      snprintf(subpath, 100, "pix%dpiy%dpiz%d", pi.x, pi.y, pi.z);
                      path_list.push_back(subpath);

                      Vertex corrvertex = boost::add_vertex(h5::path_list_to_key(path_list), g);
                      ::cvc::add_unique_edge(corrvertex, seq_prop_vertex, g);

                      // add the various derivative vertices
                      std::vector<Vertex> deriv_vertices;
                      std::string deriv_prop_key = fwd_prop_key;
                      for( auto const & deriv : deriv_chain ){
                        char deriv_prop_string[200];
                        snprintf(deriv_prop_string, 200,
                            "Ddim%d:dir%d::%s", deriv.dim, deriv.dir, deriv_prop_key.c_str());

                        std::string new_deriv_prop_key(deriv_prop_string);
                        Vertex deriv_vertex = boost::add_vertex(new_deriv_prop_key, g);
                        g[deriv_vertex].resolve.reset( new CovDevResolve(
                              deriv_prop_key, deriv.dir, deriv.dim) );
                        deriv_vertices.push_back( deriv_vertex );
                        deriv_prop_key = new_deriv_prop_key;
                      }
                      // now connect the correlator and the derivatives 
                      for(std::vector<Vertex>::reverse_iterator rit = deriv_vertices.rbegin();
                          rit != deriv_vertices.rend(); ++rit  ){
                        // only the very highest derivative is directly connected to the correlator
                        if( rit == deriv_vertices.rbegin() ){
                          ::cvc::add_unique_edge(corrvertex, *rit, g);
                        }
                        // while all derivatives are connected to each other up until the lowest derivative
                        if( rit+1 != deriv_vertices.rend() ){
                          ::cvc::add_unique_edge(*rit, *(rit+1), g);
                        }
                        // note that if other correlators get connected to some of the lower derivatives
                        // from here, the graph will be correctly connected and the processing order will
                        // be correct
                      }

                      // for the three point function, the forward and daggered propagator
                      // are in separate maps
                      g[corrvertex].resolve.reset( new 
                          CorrResolve(deriv_prop_key,
                                      seq_prop_key,
                                      mom_xchange, 
                                      gc[i_gc].as<int>(),
                                      path_list,
                                      deriv_props_data,
                                      seq_props_data,
                                      corrs_data,
                                      ::cvc::complex{1.0, 0.0} ) );
                    } // end of for( deriv_chain in deriv_chains )
                  } // end of if( deriv_chains.size() > 0 ) 
                  

        //          // adding a vertex for the correlator type allows correlators to be
        //          // processed in bunches for efficient I/O
        //          Vertex corrtypevertex = boost::add_vertex(std::string(corrtype), g);
        //          ::cvc::add_unique_edge(corrvertex, corrtypevertex, g);
                } // Dc
              } // gc
            } // gb
          } // gf
        } // gi
      } // pf
    } // pi
  } // i_dt
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
}

} //namespace(yaml)
} //namespace(cvc)

