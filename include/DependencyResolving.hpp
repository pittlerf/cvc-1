#pragma once

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

#include <string>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <highfive/H5File.hpp>

#include <cstdio>

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
    src_ts(src_ts_in), gamma(gamma_in), p(p_in), src_key(src_key_in),
    ranspinor(ranspinor_in), src(src_in) {}

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
      "[contractions] Error from prepare_gamma_timeslice_oet!",
      true,
      CVC_EXIT_MALLOC_FAILURE);
    
    if(g_verbose >= verbosity::detailed_progress){
      const bool have_source = ( (src_ts / T) == g_proc_coords[0] );
      const unsigned int vol3 = LX * LY * LZ;
      const unsigned int src_index = _GSI((src_ts % T)*vol3);
      
      Logger logger(have_source ? g_proc_id : -1, verbosity::detailed_progress, std::cout);
      logger << "TimeSliceSourceResolve: process " << g_proc_id <<
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
  const std::vector<double> & src;


  PropResolve(const std::string & prop_key_in,
              const int solver_id_in,
              const std::vector<double> & src_in,
              std::map< std::string, std::vector<double> > & props_data_in) : 
    prop_key(prop_key_in), solver_id(solver_id_in), src(src_in), 
    props_data(props_data_in) {}

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

    memcpy(workspinorA.data(), src.data(), _GSI(VOLUME)*sizeof(double));

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
  }
};

struct CovDevResolve : public ResolveDependency {
  std::string spinor_key;
  int dir;
  int dim;

  CovDevResolve(const std::string& _spinor_key, const int _dir, const int _dim) :
    spinor_key(_spinor_key), dir(_dir), dim(_dim) {}

  void operator()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
    debug_printf(0,verbosity::resolve,
        "# [CovDevResolve] Applying CovDev in dir %c, dim %c on %s\n", 
        shift_dir_names[dir],
        latDim_names[dim], 
        spinor_key.c_str());
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
  const std::string output_filename;
  std::map< std::string, std::vector<double> > & props_data;
  const mom_t p;
  const int gamma;

  CorrResolve(const std::string & propkey_in,
              const std::string & dagpropkey_in, 
              const mom_t& p_in,
              const int gamma_in, 
              const std::list<std::string> & path_list_in,
              const std::string & output_filename_in, 
              std::map< std::string, std::vector<double> > & props_data_in, 
              const ::cvc::complex & normalisation_in) :
    propkey(propkey_in),
    dagpropkey(dagpropkey_in), 
    p(p_in),
    gamma(gamma_in),
    path_list(path_list_in),
    output_filename(output_filename_in),
    props_data(props_data_in),
    normalisation(normalisation_in) {}

  void operator()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
    debug_printf(0, verbosity::resolve, 
        "# [CorrFullfill] Contracting %s\n", h5::path_list_to_key(path_list).c_str());

    std::vector<int> mom = {p.x, p.y, p.z};
    std::vector<double> corr(2*T);
    Stopwatch sw(g_cart_grid);
    contract_twopoint_gamma5_gamma_snk_only_snk_momentum(
        corr.data(), gamma, mom.data(), props_data[ dagpropkey ].data(),
        props_data[ propkey ].data());
    scale_cplx( corr.data(), T, normalisation );
    if(g_verbose >= verbosity::resolve){
      sw.elapsed_print_and_reset("contract_twopoint_gamma5_gamma_snk_only_snk_momentum local current and normalisation");
    }
    
    sw.reset();
    h5::write_t_dataset(output_filename, path_list, corr);
    if( g_verbose >= verbosity::resolve ){
      sw.elapsed_print("write_t_data in CorrFullFill");
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

  debug_printf(0,verbosity::resolve,"# [descend_and_resolve] Came up the hierarchy, ready to resolve depdency!\n");

  // in any case, when we come back here, we are ready to resolve
  //if( g[v].resolved == false ){
    debug_printf(0,verbosity::resolve,"# [descend_and_resolve] Calling resolve of %s\n",g[v].name.c_str());
    (*(g[v].resolve))();
    g[v].resolved = true;
  //}
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
}

} // namespace(cvc)
