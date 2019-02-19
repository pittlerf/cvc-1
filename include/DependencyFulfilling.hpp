#pragma once

#include "debug_printf.hpp"
#include "enums.hpp"
#include "types.h"
#include "constants.hpp"
#include "prepare_source.h"
#include "Logger.hpp"

#include <string>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <cstdio>

namespace cvc {

struct FulfillDependency {
    virtual void operator()() const = 0;
};

struct NullFulfill : public FulfillDependency {
  NullFulfill() {}
  void operator()() const 
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }
};

struct TimeSliceSourceFulfill : public FulfillDependency {
  int src_ts;
  int gamma;
  mom_t p;
  std::string src_key;
  const std::vector<double> & ranspinor;
  std::vector<double> & src;

  TimeSliceSourceFulfill(const int src_ts_in,
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
    debug_printf(0,verbosity::fulfill,
                 "TimeSliceSourceFulfill: Creating source %s\n", src_key.c_str());
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
      logger << "TimeSliceSourceFulfill: process " << g_proc_id <<
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

struct SeqSourceFulfill : public FulfillDependency {
  int src_ts;
  mom_t pf;
  std::string src_prop_key;

  SeqSourceFulfill(const int _src_ts, const mom_t & _pf, const std::string& _src_prop_key) :
    src_ts(_src_ts), pf(_pf), src_prop_key(_src_prop_key) {}

  void operator()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
    debug_printf(0,verbosity::fulfill,
                 "SeqSourceFulfill: Creating source on ts %d of %s\n", src_ts, src_prop_key.c_str());
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }
};

struct PropFulfill : public FulfillDependency {
  std::string src_key;
  std::string flav;

  PropFulfill(const std::string& _src_key, const std::string& _flav) : 
    src_key(_src_key), flav(_flav) {}

  void operator()() const
  {
    debug_printf(0,verbosity::fulfill,
                 "PropFulfill: Inverting %s on %s\n", flav.c_str(), src_key.c_str());
  }
};

struct CovDevFulfill : public FulfillDependency {
  std::string spinor_key;
  int dir;
  int dim;

  CovDevFulfill(const std::string& _spinor_key, const int _dir, const int _dim) :
    spinor_key(_spinor_key), dir(_dir), dim(_dim) {}

  void operator()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
    debug_printf(0,verbosity::fulfill,
        "CovDevFulfill: Applying CovDev in dir %c, dim %c on %s\n", 
        shift_dir_names[dir],
        latDim_names[dim], 
        spinor_key.c_str());
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }

};

struct CorrFulfill : public FulfillDependency {
  std::string propkey;
  std::string dagpropkey;
  mom_t p;
  int gamma;

  CorrFulfill(const std::string& _propkey, const std::string& _dagpropkey, const mom_t& _p, const int _gamma) :
    propkey(_propkey), dagpropkey(_dagpropkey), p(_p), gamma(_gamma) {}

  void operator()() const
  {
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
    debug_printf(0,verbosity::fulfill, 
        "CorrFullfill: Contracting %s+-g%d/px%dpy%dpz%d-%s\n",
        dagpropkey.c_str(), gamma, p.x, p.y, p.z, propkey.c_str());
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }
};

/**
 * @brief recursively fulfill dependencies starting at a particular vertex
 *
 * @tparam Graph
 * @param v
 * @param g
 */
template <typename Graph>
static inline void descend_and_fulfill(typename boost::graph_traits<Graph>::vertex_descriptor v,
                                       Graph & g)
{
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
  debug_printf(0,verbosity::fulfill,"Entered %s\n",g[v].name.c_str());
  
  // if we hit a vertex which can be fulfilled immediately, let's do so
  // this will break one class of infinite recursions
  if( g[v].independent && !g[v].fulfilled ){
    debug_printf(0,verbosity::fulfill,"Calling fulfill of %s\n",g[v].name.c_str());
    (*(g[v].fulfill))();
    g[v].fulfilled = true;
  }
 
  // otherwise we descend further, but never up the level hierarchy
  typename boost::graph_traits<Graph>::out_edge_iterator e, e_end;
  for( boost::tie(e, e_end) = boost::out_edges(v, g); e != e_end; ++e)
    if( g[boost::target(*e, g)].fulfilled == false && 
        g[boost::target(*e, g)].level < g[v].level ){
      debug_printf(0,verbosity::fulfill,"Descending into %s\n",g[boost::target(*e, g)].name.c_str());
      descend_and_fulfill( boost::target(*e, g), g );
    }

  debug_printf(0,verbosity::fulfill,"Came up the hierarchy, ready to fulfill!\n");

  // in any case, when we come back here, we are ready to fulfill
  //if( g[v].fulfilled == false ){
    debug_printf(0,verbosity::fulfill,"Calling fulfill of %s\n",g[v].name.c_str());
    (*(g[v].fulfill))();
    g[v].fulfilled = true;
  //}
#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif
}

} // namespace(cvc)
