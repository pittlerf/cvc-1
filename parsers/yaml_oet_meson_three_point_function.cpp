#include "cvc_global.h"
#include "cvc_complex.h"

#include "DependencyResolving.hpp"
#include "DependencyGraph.hpp"
#include "meta_types.hpp"
#include "parsers/yaml_utils.hpp"
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

void oet_meson_three_point_function(
    const YAML::Node &node,
    mom_lists_t & mom_lists,
    const int src_ts,
    std::map< std::string, ::cvc::stoch_prop_meta_t > & props_meta,
    DepGraph & props_graph,
    std::map< std::string, std::vector<double> > & props_data,
    std::map< std::string, ::cvc::H5Correlator > & corrs_data,
    DepGraph & corrs_graph,
    std::map< std::string, std::vector<::cvc::complex> > & phases_data,
    DepGraph & phases_graph,
    std::map< std::string, std::vector<double> > & seq_props_data,
    std::map< std::string, std::vector<double> > & cov_displ_props_data,
    double * const gauge_field_with_phases,
    const std::vector<double> & ranspinor)
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
    "momentum_exchange", "dirac_join", "Dc", "dagger_sequential",
    "fwd_solver_id", "bwd_solver_id", "seq_solver_id",
    "solver_driver"};

  static const std::vector<std::string> scalar_nodes{
    "id", "fwd_flav", "bwd_flav", "seq_flav", "dirac_join", "momentum_exchange",
    "Pi", "Pf", "dagger_sequential", "solver_driver",
    "seq_solver_id", "fwd_solver_id", "bwd_solver_id" };

  static const std::vector<std::string> sequence_nodes{
    "dt", "gi", "gf", "gb", "gc", "Dc" };

  check_missing_nodes(node, required_nodes,
      "cvc::yaml::oet_meson_three_point_function", "OetMesonThreePointFunction");
 
  for( const std::string name : scalar_nodes ){
    validate_nodetype(node[name],
                      std::vector<YAML::NodeType::value>{YAML::NodeType::Scalar}, 
                      name);
    if( name == "momentum_exchange" || name == "dagger_sequential" ){
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

  const bool dagger_seq = node["dagger_sequential"].as<bool>();

  // loop over the source-sink-separation as well as the momentum combinations at this level
  // further below over all gamma and cov_displative combinations
  // to construct all instances of this correlation function
  for( size_t i_dt = 0; i_dt < node["dt"].size(); ++i_dt ){
    const int dt = node["dt"][i_dt].as<int>();
    {
      logger << "# [yaml::oet_meson_three_point_function] Source-sink-separation: " <<
        dt << std::endl;
    }
    const int seq_src_ts = (src_ts + dt + T_global) % T_global;

    for( const auto & pi : mom_lists[ node["Pi"].as<std::string>() ] ){
      mom_t psrc{ pi.x, pi.y, pi.z };
      // when the forward propagator rather than the sequentual propagator is
      // daggered in the contraction, we need to supply the source momentum
      // with a minus sign
      if( !dagger_seq ){
        psrc = mom_t{ -pi.x, -pi.y, -pi.z };
      }

      for( const auto & pf : mom_lists[ node["Pf"].as<std::string>() ] ){
        const mom_t mom_xchange{ -(pi.x+pf.x), -(pi.y+pf.y), -(pi.z+pf.z) };
        // we have to be careful, if the sequential propagator is daggered as a whole,
        // the momentum that is actually injected at the sink must be supplied with a minus sign
        mom_t pseq{ pf.x, pf.y, pf.z };
        if( dagger_seq ){
          pseq = mom_t{ -pf.x, -pf.y, -pf.z };
        }

        if( node["momentum_exchange"].as<std::string>() == "false" ){
          if( mom_xchange.x != 0 || mom_xchange.y != 0 || mom_xchange.z != 0 ){
            continue;
          }
        }

        char phase_string[100];
        snprintf(phase_string, 100, "px%dpy%dpz%d", psrc.x, psrc.y, psrc.z);
        const std::string psrc_key(phase_string);
        Vertex psrc_vertex = boost::add_vertex(psrc_key, phases_graph);
        phases_graph[psrc_vertex].resolve.reset( new MomentumPhaseResolve(
              psrc_key,
              phases_data,
              psrc) );

        snprintf(phase_string, 100, "px%dpy%dpz%d", mom_xchange.x, mom_xchange.y, mom_xchange.z);
        const std::string mom_xchange_key(phase_string);
        Vertex mom_xchange_vertex = boost::add_vertex(mom_xchange_key, phases_graph);
        phases_graph[mom_xchange_vertex].resolve.reset( new MomentumPhaseResolve(
              mom_xchange_key,
              phases_data,
              mom_xchange) );
        
        snprintf(phase_string, 100, "px%dpy%dpz%d", pseq.x, pseq.y, pseq.z);
        const std::string pseq_key(phase_string);
        Vertex pseq_vertex = boost::add_vertex(pseq_key, phases_graph);
        phases_graph[pseq_vertex].resolve.reset( new MomentumPhaseResolve(
              pseq_key,
              phases_data,
              pseq) );

        {
          std::vector<int> pivec{ pi.x, pi.y, pi.z };
          std::vector<int> pfvec{ pf.x, pf.y, pf.z };
          std::vector<int> pcvec{ mom_xchange.x, mom_xchange.y, mom_xchange.z };
          logger << "# [yaml::oet_meson_three_point_function] Momentum: (" << pivec[0];
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
              // take care of the forward and backward propagators
              stoch_prop_meta_t fwd_prop_meta(psrc,
                                              gi[i_gi].as<int>(),
                                              src_ts,
                                              node["fwd_flav"].as<std::string>(),
                                              node["solver_driver"].as<std::string>(),
                                              node["fwd_solver_id"].as<int>() );
              const std::string fwd_prop_key = fwd_prop_meta.key();
              
              if( props_meta.count(fwd_prop_key) == 0 ){
                ts_stoch_src_meta_t fwd_prop_src_meta(psrc, gi[i_gi].as<int>(), src_ts);
                props_meta[fwd_prop_key] = fwd_prop_meta;
                Vertex fwd_prop_vertex = boost::add_vertex(fwd_prop_key, props_graph);
                props_graph[fwd_prop_vertex].resolve.reset(
                    new PropResolve(
                      fwd_prop_key,
                      node["fwd_solver_id"].as<int>(),
                      props_data,
                      new CreateGammaTimeSliceSource(
                        src_ts,
                        gi[i_gi].as<int>(),
                        psrc,
                        fwd_prop_src_meta.key(),
                        phases_data,
                        psrc_key,
                        ranspinor)
                      )
                    );

                // all simple propagators of a given flavour (id) will be collected in one group
                // that way, the MG setup, if any, can be used optimally
                Vertex fwd_flav_vertex = boost::add_vertex(node["fwd_flav"].as<std::string>(), props_graph);
                ::cvc::add_unique_edge(fwd_prop_vertex, fwd_flav_vertex, props_graph); 
              }


              // the backward propagator for the sequential propagator is
              // always a zero momentum one
              stoch_prop_meta_t bwd_prop_meta(zero_mom,
                                              gb[i_gb].as<int>(),
                                              src_ts,
                                              node["bwd_flav"].as<std::string>(),
                                              node["solver_driver"].as<std::string>(),
                                              node["bwd_solver_id"].as<int>() );
              const std::string bwd_prop_key = bwd_prop_meta.key();

              if( props_meta.count(bwd_prop_key) == 0 ){
                ts_stoch_src_meta_t bwd_prop_src_meta(zero_mom, gb[i_gb].as<int>(), src_ts);
                props_meta[bwd_prop_key] = bwd_prop_meta;
                Vertex bwd_prop_vertex = boost::add_vertex(bwd_prop_key, props_graph);
                props_graph[bwd_prop_vertex].resolve.reset(
                    new PropResolve(
                      bwd_prop_key,
                      node["bwd_solver_id"].as<int>(),
                      props_data,
                      new CreateGammaTimeSliceSource(
                        src_ts,
                        gb[i_gb].as<int>(),
                        zero_mom,
                        bwd_prop_src_meta.key(),
                        phases_data,
                        "px0py0pz0",
                        ranspinor)
                      )
                    );

                Vertex bwd_flav_vertex = boost::add_vertex(node["bwd_flav"].as<std::string>(), props_graph);
                ::cvc::add_unique_edge(bwd_prop_vertex, bwd_flav_vertex, props_graph);
              }

              for( size_t i_gc = 0; i_gc < gc.size(); ++i_gc ){
                // we have to be careful here too: in the case dagger_seq is true,
                // the Dirac structure at the sink is implicitly daggered
                // -> there may be additional signs or phase changes
                // which are not taken into account!
                // what is taken into account, however, is the sign change
                // on the sequential (final) momentum if dagger_seq is true
                // -> see definition of pseq above
                ::cvc::seq_stoch_prop_meta_t seq_prop(pseq,
                                                      gf[i_gf].as<int>(),
                                                      seq_src_ts,
                                                      src_ts,
                                                      node["seq_flav"].as<std::string>(),
                                                      zero_mom, // this is the momentum carried by the backward
                                                                // propagator
                                                      gb[i_gb].as<int>(),
                                                      node["bwd_flav"].as<std::string>()
                                                      );

                const std::string seq_prop_key = seq_prop.key();

                char seq_src_string[200];
                snprintf(seq_src_string, 200, "g%d_px%dpy%dpz%d::ts%d::%s",
                    gf[i_gf].as<int>(), pseq.x, pseq.y, pseq.z,
                    seq_src_ts,
                    bwd_prop_key.c_str());
                const std::string seq_src_key(seq_src_string);

                Vertex seq_prop_vertex = boost::add_vertex(seq_prop_key, corrs_graph);
                corrs_graph[seq_prop_vertex].resolve.reset( new PropResolve(
                      seq_prop_key,
                      node["seq_solver_id"].as<int>(),
                      seq_props_data,
                      new CreateSequentialGammaTimeSliceSource(
                        src_ts,
                        seq_src_ts,
                        gf[i_gf].as<int>(),
                        pseq,
                        seq_src_key,
                        bwd_prop_key,
                        props_data,
                        pseq_key,
                        phases_data) ) );

                char corrtype[100];
                if( dagger_seq ){
                  snprintf(corrtype, 100, "s%s%s+-g-%s-g",
                           node["bwd_flav"].as<std::string>().c_str(),
                           node["seq_flav"].as<std::string>().c_str(),
                           node["fwd_flav"].as<std::string>().c_str());
                } else {
                  snprintf(corrtype, 100, "s%s%s-g-%s+-g",
                           node["bwd_flav"].as<std::string>().c_str(),
                           node["seq_flav"].as<std::string>().c_str(),
                           node["fwd_flav"].as<std::string>().c_str());
                }
                for( size_t i_Dc = 0; i_Dc < Dc.size(); ++i_Dc ){
                  {
                    logger << "# [yaml::oet_meson_three_point_function] Dirac: (" << 
                      gf[i_gf].as<int>() << "," <<
                      gc[i_gc].as<int>() <<
                      "," << gi[i_gi].as<int>() << ")  " <<
                      "BwdDirac: " << gb[i_gb].as<int>() << "  " <<
                      "Covariant Displacements: " << Dc[i_Dc].as<int>() <<
                      std::endl;
                  }

                  std::vector< std::vector<deriv_t> > cov_displ_chains;
                  cov_displ_chains = create_derivative_chains(0, Dc[i_Dc].as<int>(), cov_displ_chains);
                 
                  if( cov_displ_chains.size() == 0 ){
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

                    Vertex corrvertex = boost::add_vertex(h5::path_list_to_key(path_list), corrs_graph);
                    ::cvc::add_unique_edge(corrvertex, seq_prop_vertex, corrs_graph);

                    // for the three point function, the forward and daggered propagator
                    // are in separate maps
                    // note that when the sequential and forward propagators are contracted
                    // the gamma5 from employing gamma5-hermiticity is explicitly included
                    // in the contraction routine contract_twopoint_gamma5_gamma_snk_only_snk_momentum
                    if( dagger_seq ){
                      corrs_graph[corrvertex].resolve.reset( new 
                          CorrResolve(fwd_prop_key,
                                      seq_prop_key,
                                      mom_xchange, 
                                      gc[i_gc].as<int>(),
                                      path_list,
                                      props_data,
                                      seq_props_data,
                                      corrs_data,
                                      phases_data,
                                      ::cvc::complex{1.0, 0.0} ) );
                    } else {
                      corrs_graph[corrvertex].resolve.reset( new 
                          CorrResolve(seq_prop_key,
                                      fwd_prop_key,
                                      mom_xchange, 
                                      gc[i_gc].as<int>(),
                                      path_list,
                                      seq_props_data,
                                      props_data,
                                      corrs_data,
                                      phases_data,
                                      ::cvc::complex{1.0, 0.0} ) );
                    }
                  } else { // --> cov_displ_chains.size() > 0
                  
                    for( auto const & cov_displ_chain : cov_displ_chains ){
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
                      // in the creation of the key for multiple covariant displacements,
                      // we reverse the iteration order over the elements in the derivative
                      // chain as the first element of the chain is applied first,
                      // followed by the second and then by the third
                      // Thus, in order to obtain operator notation (with the first application
                      // appearing *right-most* in the key), we need to push_back the elements
                      // from the last to the first.
                      for( auto cov_displ_i = cov_displ_chain.rbegin(); cov_displ_i != cov_displ_chain.rend(); ++cov_displ_i ){
                        snprintf(subpath, 100, "Ddim%d_dir%d", cov_displ_i->dim, cov_displ_i->dir);
                        path_list.push_back(subpath);
                      }
                      snprintf(subpath, 100, "gi%d", gi[i_gi].as<int>());
                      path_list.push_back(subpath);
                      snprintf(subpath, 100, "pix%dpiy%dpiz%d", pi.x, pi.y, pi.z);
                      path_list.push_back(subpath);

                      Vertex corrvertex = boost::add_vertex(h5::path_list_to_key(path_list), corrs_graph);
                      ::cvc::add_unique_edge(corrvertex, seq_prop_vertex, corrs_graph);

                      // add the various covariantly displaced vertices
                      std::vector<Vertex> cov_displ_vertices;

                      // generate a key for the covariantly displaced quark propagator
                      // in operator ordering (see snprintf below, displacements are added
                      // to the key from the left)
                      std::string cov_displ_prop_key = fwd_prop_key;
                      for( size_t i_cov_displ = 0; i_cov_displ < cov_displ_chain.size(); ++i_cov_displ ){
                        const deriv_t & cov_displ = cov_displ_chain[i_cov_displ];

                        char cov_displ_prop_string[300];
                        int nchar = snprintf(cov_displ_prop_string, 300,
                                             "Ddim%d_dir%d::%s", cov_displ.dim, cov_displ.dir,
                                             cov_displ_prop_key.c_str());
                        if( nchar >= 300 ){
                          char msg[200];
                          snprintf(msg, 200,
                                   "[yaml::oet_meson_three_point_function] Exceeded maximum number of characters (299)"
                                   "in the writing of 'cov_displ_prop_string'. You should increase the number"
                                   "of characters that can be written here! Sorry about the poor implementation!\n");
                          EXIT_WITH_MSG(CVC_EXIT_SNPRINTF_OVERFLOW, msg);
                        }

                        std::string new_cov_displ_prop_key(cov_displ_prop_string);
                        Vertex cov_displ_vertex = boost::add_vertex(new_cov_displ_prop_key, corrs_graph);
                        corrs_graph[cov_displ_vertex].resolve.reset( new CovDisplResolve(
                              cov_displ_prop_key, // src key, in the first iteration, this
                                                  // is sequal to fwd_prop_key
                              new_cov_displ_prop_key, // displaced prop key
                              props_data,
                              cov_displ_props_data,
                              i_cov_displ == 0,
                              cov_displ.dir,
                              cov_displ.dim,
                              gauge_field_with_phases) );
                        if( i_cov_displ == 0 ){
                          corrs_graph[cov_displ_vertex].independent = true;
                        }
                        cov_displ_vertices.push_back( cov_displ_vertex );

                        // replace the key with the next one in line
                        cov_displ_prop_key = new_cov_displ_prop_key;
                      }

                      // now connect the correlator and the covariantly displaced propagators 
                      for(std::vector<Vertex>::reverse_iterator rit = cov_displ_vertices.rbegin();
                          rit != cov_displ_vertices.rend(); ++rit  ){
                        // all displacement are connected to each other up until the lowest one
                        // which is not connected because the basic propagators are assumed
                        // to exist
                        if( rit+1 != cov_displ_vertices.rend() ){
                          debug_printf(0,verbosity::graph_connections,"Connecting %s to %s\n",
                              corrs_graph[(*rit)].name.c_str(), corrs_graph[*(rit+1)].name.c_str() );
                          ::cvc::add_unique_edge(*rit, *(rit+1), corrs_graph);
                        }
                        // only the very highest displacement is directly connected to the correlator
                        // note that it's important to connect the correlator AFTER connecting the
                        // displacement vertices, such that that the correlator ends up
                        // a level higher in the hierarchy
                        if( rit == cov_displ_vertices.rbegin() ){
                          debug_printf(0,verbosity::graph_connections,"Connecting %s to %s\n",
                              corrs_graph[corrvertex].name.c_str(), corrs_graph[(*rit)].name.c_str() );
                          ::cvc::add_unique_edge(corrvertex, *rit, corrs_graph);
                        }
                        // note that if other correlators get connected to some of the lower displacements
                        // from here, the graph will be correctly connected and the processing order will
                        // be correct
                      }
 
                      // now just complete the task by creating the actual correlator vertex
                      if( dagger_seq ){
                        corrs_graph[corrvertex].resolve.reset( new 
                            CorrResolve(cov_displ_prop_key,
                                        seq_prop_key,
                                        mom_xchange, 
                                        gc[i_gc].as<int>(),
                                        path_list,
                                        cov_displ_props_data,
                                        seq_props_data,
                                        corrs_data,
                                        phases_data,
                                        ::cvc::complex{1.0, 0.0} ) );
                      } else {
                        // FIXME: when the covariantly displaced propagator
                        // is daggered, one has to of course think carefully about what this means
                        // for the displacements!
                        corrs_graph[corrvertex].resolve.reset( new 
                            CorrResolve(seq_prop_key,
                                        cov_displ_prop_key,
                                        mom_xchange, 
                                        gc[i_gc].as<int>(),
                                        path_list,
                                        seq_props_data,
                                        cov_displ_props_data,
                                        corrs_data,
                                        phases_data,
                                        ::cvc::complex{1.0, 0.0} ) );
                      }
                    } // end of for( cov_displ_chain in cov_displ_chains )
                  } // end of if( cov_displ_chains.size() > 0 ) 
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

