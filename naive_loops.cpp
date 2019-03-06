#define MAIN_PROGRAM
#include "global.h"
#undef MAIN_PROGRAM

#include "DependencyGraph.hpp"
#include "DependencyResolving.hpp"
#include "meta_types.hpp"
#include "algorithms.hpp"
#include "yaml_parsers.hpp"
#include "Core.hpp"
#include "Logger.hpp"
#include "debug_printf.hpp"
#include "init_g_gauge_field.hpp"
#include "cvc_utils.h"
#include "mpi_init.h"
#include "enums.hpp"
#include "ParallelMT19937_64.hpp"
#include "h5utils.hpp"
#include "make_phase_field.hpp"
#include "loop_tools.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <yaml-cpp/yaml.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range.hpp>

using namespace cvc;

std::map< std::string, std::vector< std::vector<int> > > momentum_lists;
std::map< std::string, stoch_prop_meta_t > ts_prop_meta; 

int main(int argc, char ** argv){
  int exitstatus;
  int io_proc = -1;
  double * gauge_field_with_phase = NULL;

  Core core(argc,argv,"naive_loops");
  if( !core.is_initialised() ){
		debug_printf(g_proc_id,0,"[naive_loops] Core initialisation failure!\n");
		return(CVC_EXIT_CORE_INIT_FAILURE);
  }

  // handy stopwatch
  Stopwatch sw(g_cart_grid);
  // console logger on proc_id 0
  Logger logger(0,verbosity::basic_progress,std::cout);
  Logger all_logger(0,0,std::cout);

  mpi_init_xchange_contraction(2);
  mpi_init_xchange_eo_spinor ();

  CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      init_g_gauge_field(),
      "[naive_loops] Error initialising gauge field!",
      true,
      CVC_EXIT_GAUGE_INIT_FAILURE);

  CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      gauge_field_eq_gauge_field_ti_phase ( &gauge_field_with_phase, g_gauge_field, co_phase_up ),
      "[naive_loops] Error from gauge_field_eq_gauge_field_ti_phase!",
      true,
      CVC_EXIT_UTIL_FUNCTION_FAILURE);

  ParallelMT19937_64 rng( (unsigned long long)(g_seed^Nconf) );
  std::vector<double> ranspinor( _GSI(VOLUME) );

  io_proc = get_io_proc ();
  if( io_proc < 0 ) {
    PRINT_STATUS(io_proc, "[naive_loops] Error, io proc must be ge 0!");
    EXIT(14);
  }
  fprintf(stdout, "# [naive_loops] proc%.4d has io proc id %d\n", g_cart_id, io_proc );

  all_logger << "# [naive_loops] Number of source locations: " << g_source_location_number << std::endl;

  // the loops will be written to this file
  char loops_h5_filename_cstring[100];
  snprintf(loops_h5_filename_cstring, 100, "loops.%04d.h5", Nconf);
  const std::string loops_h5_filename(loops_h5_filename_cstring);

  std::vector<double> propagator(_GSI(VOLUMEPLUSRAND));
  std::vector<double> source(_GSI(VOLUMEPLUSRAND));

  // pre-generate the phase fields that we need
  std::vector< std::vector<::cvc::complex> > phase_fields(g_source_momentum_number);
  for(unsigned int i_mom = 0; i_mom < g_source_momentum_number; ++i_mom){
    phase_fields[i_mom].resize(VOLUME);
    mom_t mom{g_source_momentum_list[i_mom][0],
              g_source_momentum_list[i_mom][1],
              g_source_momentum_list[i_mom][2]};
    make_phase_field(phase_fields[i_mom], mom);
  }


  char istoch_cstring[20];
  char gamma_cstring[10];
  char p_cstring[100];
  for( int i_sample = 0; i_sample < g_nsample; i_sample++ ) {

    // build h5 path lists for the output
    std::vector< std::list<std::string> > path_lists(16*g_source_momentum_number);
    snprintf(istoch_cstring, 20, "istoch_%04d", i_sample);
    for(unsigned int i_gamma = 0; i_gamma < 16; ++i_gamma){
      snprintf(gamma_cstring, 10, "g%d", i_gamma);
      for(unsigned int i_mom = 0; i_mom < g_source_momentum_number; ++i_mom){
        snprintf(p_cstring, 100, "px%dpy%dpz%d", 
                                 g_source_momentum_list[i_mom][0],
                                 g_source_momentum_list[i_mom][1],
                                 g_source_momentum_list[i_mom][2]);

        path_lists[g_source_momentum_number*i_gamma + i_mom].push_back(istoch_cstring);
        path_lists[g_source_momentum_number*i_gamma + i_mom].push_back(gamma_cstring);
        path_lists[g_source_momentum_number*i_gamma + i_mom].push_back(p_cstring);
      }
    }

    // allocate and zero out memory for gamma and momentum projected loops
    std::vector<double> loops(2*T*16*g_source_momentum_number, 0.0);

    sw.reset();
    // generate a random source with the properties 
    //   z[i]*conj(z[i]) = 1.0 \forall i
    //   \lim_{n -> \infty} 1\n \sum z_n[i] = 0
    //   \lim_{n -> \infty} 1\n \sum z_n[i] conj(z_n[j]) = \delta_{i,j}
    if( g_noise_type == GAUSSIAN_NOISE ){
      rng.gen_gauss(ranspinor.data(), 24);
    } else if ( g_noise_type == Z2_NOISE ){
      rng.gen_z2(ranspinor.data(), 24);
    } else {
      throw( std::invalid_argument("The only supported noise types are 'Gaussian' and 'Z2'!\n") );
    }
    sw.elapsed_print_and_reset("volume source generation");

    // copy the random spinor to the source vector to make sure that it's preserved
    double * const src_data = source.data();
    double * const ran_data = ranspinor.data();
 
    memcpy(src_data, ran_data, sizeof(double)*_GSI(VOLUME));
  
    int exitstatus = 0;
    CHECK_EXITSTATUS_NEGATIVE(
        exitstatus,
        _TMLQCD_INVERT(propagator.data(), source.data(), 0),
        "[naive_loops] Error from TMQLCD_INVERT",
        true,
        CVC_EXIT_UTIL_FUNCTION_FAILURE);

    sw.elapsed_print_and_reset("TMLQCD_INVERT");

    contract_twopoint_all_gammas_ext_momentum_list(loops.data(), phase_fields, ranspinor.data(), propagator.data());
    sw.elapsed_print_and_reset("contract_twopoint_all_gammas_ext_momentum_list");

    h5::write_loops(loops_h5_filename, loops, path_lists, 16, g_source_momentum_number);
    sw.elapsed_print("loop_io");

#ifdef HAVE_MPI
      MPI_Barrier(g_cart_grid);
#endif
  } // loop over samples

	free(gauge_field_with_phase);
  mpi_fini_xchange_contraction();
  mpi_fini_xchange_eo_spinor ();

  return 0;
}
