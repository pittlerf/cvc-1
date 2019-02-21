#pragma once

#include "global.h"
#include "Logger.hpp"
#include "debug_printf.hpp"
#include "propagator_io.h"
#include "Stopwatch.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <vector>

namespace cvc {

  struct CreateSource {
    virtual void operator()(std::vector<double> & src) const = 0;
  };
  
  struct CreateNullSource : public CreateSource {
    void operator()(std::vector<double> & src) const {}
  };

  struct CreateGammaTimeSliceSource : public CreateSource {
    const int src_ts;
    const int gamma;
    const mom_t p;
    const std::string src_key;
    const std::vector<double> & ranspinor;
  
  	CreateGammaTimeSliceSource(const int src_ts_in,
  												     const int gamma_in,
                               const mom_t & p_in,
                               const std::string & src_key_in,
                               const std::vector<double> & ranspinor_in) :
      src_ts(src_ts_in),
      gamma(gamma_in),
      p(p_in), 
      src_key(src_key_in),
      ranspinor(ranspinor_in) {}
  
    void operator()(std::vector<double> & src) const 
    {
#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
      Stopwatch sw(g_cart_grid);
      int exitstatus;
      debug_printf(0,verbosity::resolve,
                   "# [CreateGammaTimeSliceSource] Creating source %s\n", src_key.c_str());
      CHECK_EXITSTATUS_NONZERO(
        exitstatus,
        prepare_gamma_timeslice_oet(
          src.data(),
          ranspinor.data(),
          gamma,
          src_ts,
          p),
        "[CreateGammaTimeSliceSource] Error from prepare_gamma_timeslice_oet!",
        true,
        CVC_EXIT_MALLOC_FAILURE);

      if( g_verbose >= verbosity::detailed_progress ) sw.elapsed_print("prepare_gamma_timeslice_oet");
  
      if( g_write_source ){
#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
        std::string filename = src_key;
        std::replace( filename.begin(), filename.end(), '/', '_');
        filename += ".lime";
        write_propagator(src.data(), filename.c_str(), 0, 64);
      } 
  
#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
    }
  };

} // namespace(cvc)
