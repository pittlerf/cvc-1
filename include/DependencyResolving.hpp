#pragma once

#include "cvc_global.h"
#include "debug_printf.hpp"
#include "enums.hpp"
#include "types.h"
#include "constants.hpp"
#include "prepare_source.h"
#include "Logger.hpp"
#include "Stopwatch.hpp"
#include "dummy_solver.h"
#include "h5utils.hpp"
#include "cvc_complex.h"
#include "propagator_io.h"
#include "SourceCreators.hpp"
#include "Q_phi.h"
#include "make_phase_field.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <highfive/H5File.hpp>

#include <string>
#include <iostream>
#include <cstdio>
#include <algorithm>

namespace cvc {

struct ResolveDependency {
    virtual void operator()() const = 0;
};

struct NullResolve : public ResolveDependency {
  NullResolve() {}
  void operator()() const 
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }
};

struct MomentumPhaseResolve : public ResolveDependency {
  const std::string mom_key;
  std::map< std::string, std::vector<::cvc::complex> > & phases_data;
  const mom_t p;

  MomentumPhaseResolve(
      const std::string & mom_key_in,
      std::map< std::string, std::vector<::cvc::complex> > & phases_data_in,
      const mom_t & p_in) :
    mom_key(mom_key_in),
    phases_data(phases_data_in),
    p(p_in) {}
  
  void operator()() const 
  {
    debug_printf(0,verbosity::resolve,
        "[MomentumPhaseResolve] Making phase field for px:%d py:%d pz:%d\n",
        p.x, p.y, p.z);

    if( !phases_data.count(mom_key) ){
      size_t const VOL3 = LX*LY*LZ;
      phases_data.emplace( std::make_pair(mom_key,
                                          std::vector<::cvc::complex>(VOL3) ) );
      make_phase_field(phases_data[mom_key], p);
    }
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }
}; // MomentumPhaseResolve

struct TimeSliceSourceResolve : public ResolveDependency {
  const int src_ts;
  const int gamma;
  const mom_t p;
  const std::string src_key;
  const std::vector<double> & ranspinor;
  std::vector<double> & src;

  TimeSliceSourceResolve(const int src_ts_in,
                         const int gamma_in, 
                         const mom_t & p_in, 
                         const std::string & src_key_in,
                         const std::vector<double> & ranspinor_in,
                         std::vector<double> & src_in) :
    src_ts(src_ts_in),
    gamma(gamma_in),
    p(p_in),
    src_key(src_key_in),
    ranspinor(ranspinor_in),
    src(src_in) {}

  void operator ()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
    int exitstatus;
    debug_printf(0,verbosity::resolve,
                 "# [TimeSliceSourceResolve] Creating source %s\n", src_key.c_str());
    CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      prepare_gamma_timeslice_oet(
        src.data(),
        ranspinor.data(),
        gamma,
        src_ts,
        p),
      "[TimeSliceSourceResolve] Error from prepare_gamma_timeslice_oet!",
      true,
      CVC_EXIT_MALLOC_FAILURE);
    

    if( g_write_source ){
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
      std::string filename = src_key;
      std::replace( filename.begin(), filename.end(), '/', '_');
      filename += ".lime";
      write_propagator(src.data(), filename.c_str(), 0, 64);
    } 

    if(g_verbose >= verbosity::detailed_progress){
      const bool have_source = ( (src_ts / T) == g_proc_coords[0] );
      const unsigned int vol3 = LX * LY * LZ;
      const unsigned int src_index = _GSI((src_ts % T)*vol3);
      
      // no process will ever have proc_id == -1, so we can limit output to the appropriate
      // subset as follows:
      Logger logger(have_source ? g_proc_id : -1, verbosity::detailed_progress, std::cout);
      logger << "# [TimeSliceSourceResolve] process " << g_proc_id <<
        " elements 0-24" << std::endl;

      for(size_t i = src_index; i < src_index+24; ++i ){
        logger << src[i] << " ";
      }
      logger << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }
};

struct SeqSourceResolve : public ResolveDependency {
  int src_ts;
  mom_t pf;
  std::string src_prop_key;

  SeqSourceResolve(const int _src_ts, const mom_t & _pf, const std::string& _src_prop_key) :
    src_ts(_src_ts), pf(_pf), src_prop_key(_src_prop_key) {}

  void operator()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
    debug_printf(0,verbosity::resolve,
                 "# [SeqSourceResolve] Creating source on ts %d of %s\n", src_ts, src_prop_key.c_str());
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }
};

struct PropResolve : public ResolveDependency {
  const int solver_id;
  const std::string prop_key;
  std::map< std::string, std::vector<double> > & props_data;

  // source creation functor for this propagator
  std::shared_ptr<CreateSource> create_source;

  PropResolve(const std::string & prop_key_in,
              const int solver_id_in,
              std::map< std::string, std::vector<double> > & props_data_in,
              CreateSource * create_source_in) : 
    prop_key(prop_key_in),
    solver_id(solver_id_in), 
    props_data(props_data_in),
    create_source(create_source_in) {}

  void operator()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif

    debug_printf(0,verbosity::resolve,
                 "# [PropResolve] Inverting to generate propagator %s\n", prop_key.c_str());
    
    Stopwatch sw(g_cart_grid);
    std::vector<double> workspinorA(_GSI(VOLUMEPLUSRAND), 0.0);
    std::vector<double> workspinorB(_GSI(VOLUMEPLUSRAND), 0.0);

    // call the source creation functor
    (*create_source)(workspinorA);

    int exitstatus = 0;
    CHECK_EXITSTATUS_NEGATIVE(
      exitstatus,
      _TMLQCD_INVERT(workspinorB.data(), workspinorA.data(), solver_id),
      "[PropResolve] Error from TMLQCD_INVERT",
      true,
      CVC_EXIT_UTIL_FUNCTION_FAILURE);

    if(g_verbose >= verbosity::resolve) sw.elapsed_print("TMLQCD_INVERT"); 

    if( !props_data.count(prop_key) ){
      props_data.emplace( std::make_pair(prop_key,
                                         std::vector<double>(_GSI(VOLUME) ) ) );
    }
    memcpy( props_data[ prop_key ].data(), workspinorB.data(), _GSI(VOLUME)*sizeof(double) );

#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif

    if( g_write_propagator ){
      std::string filename = prop_key + ".lime";
      write_propagator(workspinorB.data(), filename.c_str(), 0, 64);
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
    }

  }
};

struct CovDisplResolve : public ResolveDependency {
  const std::string cov_displ_prop_key;
  const std::string source_prop_key;
  std::map< std::string, std::vector<double> > & props_data;
  std::map< std::string, std::vector<double> > & cov_displ_props_data;
  const bool src_in_props_data;
  const int dir;
  const int dim;
  double * const gauge_field;

  CovDisplResolve(const std::string& source_prop_key_in, 
                  const std::string& cov_displ_prop_key_in,
                  std::map< std::string, std::vector<double> > & props_data_in,
                  std::map< std::string, std::vector<double> > & cov_displ_props_data_in, 
                  const bool src_in_props_data_in,
                  const int dir_in,
                  const int dim_in,
                  double * const gauge_field_in) :
    source_prop_key(source_prop_key_in),
    cov_displ_prop_key(cov_displ_prop_key_in),
    props_data(props_data_in),
    cov_displ_props_data(cov_displ_props_data_in),
    src_in_props_data(src_in_props_data_in),
    dir(dir_in),
    dim(dim_in),
    gauge_field(gauge_field_in) {}

  void operator()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif

    Stopwatch sw(g_cart_grid);
    debug_printf(0,verbosity::resolve,
        "# [CovDisplResolve] Applying covariant displacement in dir %c, dim %c on %s\n", 
        shift_dir_names[dir],
        latDim_names[dim], 
        source_prop_key.c_str());

    if( !cov_displ_props_data.count( cov_displ_prop_key ) ){
      cov_displ_props_data.emplace( std::make_pair(cov_displ_prop_key,
                                                   std::vector<double>(_GSI(VOLUME) ) ) );
    }
    // either the source is a basic propagator or has already been displaced
    double * const src = src_in_props_data ? props_data[ source_prop_key ].data() :
                                             cov_displ_props_data[ source_prop_key ].data();
                                             
    spinor_field_eq_cov_displ_spinor_field(cov_displ_props_data[ cov_displ_prop_key ].data(),
                                           src,
                                           dim,
                                           dir,
                                           gauge_field);

    if( g_verbose >= verbosity::resolve ) sw.elapsed_print("covariant displacement");

#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }

};

struct CorrResolve : public ResolveDependency {
  const ::cvc::complex normalisation;
  const std::string propkey;
  const std::string dagpropkey;
  const std::list<std::string> path_list;
  std::map< std::string, std::vector<double> > & props_data;
  std::map< std::string, std::vector<double> > & dag_props_data;
  std::map< std::string, H5Correlator > & corrs_data;
  std::map< std::string, std::vector<::cvc::complex> > & phases_data; 
  const mom_t p;
  const int gamma;

  CorrResolve(const std::string & propkey_in,
              const std::string & dagpropkey_in, 
              const mom_t& p_in,
              const int gamma_in, 
              const std::list<std::string> & path_list_in,
              std::map< std::string, std::vector<double> > & props_data_in,
              std::map< std::string, std::vector<double> > & dag_props_data_in,
              std::map< std::string, H5Correlator > & corrs_data_in,
              std::map< std::string, std::vector<::cvc::complex> > & phases_data_in, 
              const ::cvc::complex & normalisation_in) :
    propkey(propkey_in),
    dagpropkey(dagpropkey_in), 
    p(p_in),
    gamma(gamma_in),
    path_list(path_list_in),
    props_data(props_data_in),
    dag_props_data(dag_props_data_in),
    corrs_data(corrs_data_in),
    phases_data(phases_data_in),
    normalisation(normalisation_in) {}

  void operator()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
    const std::string key = h5::path_list_to_key(path_list);

    debug_printf(0, verbosity::resolve, 
        "# [CorrResolve] Contracting [(%s)+ g5 | g%d px%dpy%dpz%d | %s] to form %s\n",
        dagpropkey.c_str(), 
        gamma, p.x, p.y, p.z,
        propkey.c_str(), key.c_str() );

    if( !corrs_data.count(key) ){
      corrs_data.emplace( std::make_pair(key,
                                         H5Correlator(path_list, 2*T) ) );
    }

    char mom_key[30];
    snprintf(mom_key, 30, "px%dpy%dpz%d", p.x, p.y, p.z);
    const int mom_vec[3]={p.x, p.y, p.z}; 

    Stopwatch sw(g_cart_grid);
    contract_twopoint_gamma5_gamma_only_ext_momentum(
        corrs_data[key].storage.data(), 
        gamma,
        phases_data[std::string(mom_key)].data(),
        dag_props_data[ dagpropkey ].data(),
        props_data[ propkey ].data() );
    scale_cplx( corrs_data[key].storage.data(), T, normalisation );
    if(g_verbose >= verbosity::resolve){
      sw.elapsed_print_and_reset("contract_twopoint_gamma5_gamma_only_ext_momentum local current and normalisation");
    }

#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }
};

/**
 * @brief recursively resolve dependencies starting at a particular vertex
 *
 * @tparam Graph
 * @param v
 * @param g
 */
template <typename Graph>
static inline void descend_and_resolve(typename boost::graph_traits<Graph>::vertex_descriptor v,
                                       Graph & g)
{
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  debug_printf(0,verbosity::resolve,"# [descend_and_resolve] Entered %s\n",g[v].name.c_str());
  
  // if we hit a vertex which can be resolved immediately, let's do so
  // this will break one class of infinite recursions
  if( g[v].independent && !g[v].resolved ){
    debug_printf(0,verbosity::resolve,"# [descend_and_resolve] Calling resolve of %s\n",g[v].name.c_str());
    (*(g[v].resolve))();
    g[v].resolved = true;
  }
 
  // otherwise we descend further, but never up the level hierarchy
  typename boost::graph_traits<Graph>::out_edge_iterator e, e_end;
  for( boost::tie(e, e_end) = boost::out_edges(v, g); e != e_end; ++e)
    if( g[boost::target(*e, g)].resolved == false && 
        g[boost::target(*e, g)].level < g[v].level ){
      debug_printf(0,verbosity::resolve,"# [descend_and_resolve] Descending into %s\n",g[boost::target(*e, g)].name.c_str());
      descend_and_resolve( boost::target(*e, g), g );
    }

  debug_printf(0,verbosity::resolve,"# [descend_and_resolve] Came up the hierarchy, ready to resolve depedency!\n");

  // in any case, when we come back here, we are ready to resolve
  if( g[v].resolved == false ){
    debug_printf(0,verbosity::resolve,"# [descend_and_resolve] Calling resolve of %s\n",g[v].name.c_str());
    (*(g[v].resolve))();
    g[v].resolved = true;
  }
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
}

} // namespace(cvc)
