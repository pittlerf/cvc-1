#pragma once

#include "cvc_global.h"
#include "Logger.hpp"
#include "debug_printf.hpp"
#include "propagator_io.h"
#include "Stopwatch.hpp"
#include "prepare_source.h"

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
  
  struct CreateSequentialGammaTimeSliceSource : public CreateSource {
    const int src_ts;
    const int seq_src_ts;
    const int seq_gamma;
    const mom_t seq_p;
    const std::string seq_src_key;
    const std::string bwd_prop_key;
    std::map< std::string, std::vector<double> > & props_data;
    const std::string phase_key;
    std::map< std::string, std::vector<::cvc::complex> > & phases_data;
  
  	CreateSequentialGammaTimeSliceSource(const int src_ts_in,
                                         const int seq_src_ts_in,
  												               const int seq_gamma_in,
                                         const mom_t & seq_p_in,
                                         const std::string & seq_src_key_in,
                                         const std::string & bwd_prop_key_in,
                                         std::map<std::string, std::vector<double> > & props_data_in,
                                         const std::string & phase_key_in,
                                         std::map<std::string, std::vector<::cvc::complex> > & phases_data_in) :
      seq_src_ts(seq_src_ts_in),
      src_ts(src_ts_in),
      seq_gamma(seq_gamma_in),
      seq_p(seq_p_in), 
      seq_src_key(seq_src_key_in),
      bwd_prop_key(bwd_prop_key_in),
      props_data(props_data_in),
      phases_data(phases_data_in),
      phase_key(phase_key_in) {}
  
    void operator()(std::vector<double> & src) const 
    {
#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
      Stopwatch sw(g_cart_grid);
      int exitstatus;
      debug_printf(0,verbosity::resolve,
                   "# [CreateSequentialGammaTimeSliceSource] Creating source %s\n", seq_src_key.c_str());
      CHECK_EXITSTATUS_NONZERO(
        exitstatus,
        prepare_gamma_timeslice_oet(
          src.data(),
          props_data[bwd_prop_key].data(),
          seq_gamma,
          seq_src_ts,
          seq_p),
        "[CreateSequentialGammaTimeSliceSource] Error from prepare_gamma_timeslice_oet!",
        true,
        CVC_EXIT_MALLOC_FAILURE);

      if( g_verbose >= verbosity::detailed_progress ) sw.elapsed_print("prepare_gamma_timeslice_oet");
  
      if( g_write_source ){
#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
        std::string filename = seq_src_key;
        std::replace( filename.begin(), filename.end(), '/', '_');
        filename += ".lime";
        write_propagator(src.data(), filename.c_str(), 0, 64);
      } 
  
#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
    }
  };

  struct CreateGammaTimeSliceSource : public CreateSource {
    const int src_ts;
    const int gamma;
    const mom_t p;
    const std::string src_key;
    const std::vector<double> & ranspinor;
    const std::string phase_key;
    std::map< std::string, std::vector<::cvc::complex> > & phases_data;
  
  	CreateGammaTimeSliceSource(const int src_ts_in,
  												     const int gamma_in,
                               const mom_t & p_in,
                               const std::string & src_key_in,
                               std::map< std::string, std::vector<::cvc::complex> > & phases_data_in,
                               const std::string & phase_key_in,
                               const std::vector<double> & ranspinor_in) :
      src_ts(src_ts_in),
      gamma(gamma_in),
      p(p_in), 
      src_key(src_key_in),
      phase_key(phase_key_in),
      phases_data(phases_data_in),
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
