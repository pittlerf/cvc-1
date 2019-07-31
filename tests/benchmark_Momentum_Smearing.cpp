
#define MAIN_PROGRAM
#include "global.h"
#undef MAIN_PROGRAM

#include "cvc_utils.h"
#include "Q_phi.h"
#include "Stopwatch.hpp"
#include "Core.hpp"
#include "ParallelMT19937_64.hpp"
#include "enums.hpp"
#include "init_g_gauge_field.hpp"
#include "debug_printf.hpp"
#include "smearing_techniques.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <vector>
#include <cstring>

#include <omp.h>

using namespace cvc;

int main(int argc, char** argv)
{
  Core core(argc,argv);
  if( !(core.is_initialised()) ){
    std::cout << "Core initialisation failed!\n";
    return(CVC_EXIT_CORE_INIT_FAILURE);
  }

  int exitstatus = 0; 

#ifdef HAVE_MPI
  Stopwatch sw(g_cart_grid);
#else
  Stopwatch sw(0);
#endif
  
#pragma omp parallel
  {
    if( omp_get_thread_num() == 0 ){
      char msg[400];
      snprintf(msg, 400, "### [proc %3d]: omp_num_threads = %d\n", g_proc_id, omp_get_num_threads());
      debug_printf(g_proc_id, 0, msg); 
    }
  }

  // initialize and *load* the gauge field
  CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      init_g_gauge_field(),
      "[benchmark_Momentum_Smearing] Error initialising gauge field!",
      true, 
      CVC_EXIT_GAUGE_INIT_FAILURE);

  double * gauge_field_smeared;
  CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      alloc_gauge_field(&gauge_field_smeared, VOLUMEPLUSRAND),
      "[benchmark_Momentum_Smearing] Error allocating memory for smeared gauge field!",
      true,
      CVC_EXIT_GAUGE_INIT_FAILURE);

  memcpy( gauge_field_smeared, g_gauge_field, 72*VOLUME*sizeof(double) ); 

  sw.reset();
  CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      APE_Smearing( gauge_field_smeared, 2.5, 15 ),
      "[benchmark_Momentum_Smearing] Error from APE_Smearing!",
      true,
      CVC_EXIT_UTIL_FUNCTION_FAILURE);
  sw.elapsed_print("APE smearing");

  const unsigned long long global_volume = g_nproc*VOLUME;

  ParallelMT19937_64 rng(897872134ULL);

  // let's use one particular momentum configuration, doesn't matter which one
  const mom_t mom = { 1, 2, -3 };

  const unsigned int outer_iterations = 32;
  const unsigned int inner_iterations = 90;
  
  std::vector<std::vector<double>> spinor1(outer_iterations);
  for( auto & elem : spinor1 ) elem.resize(24*(VOLUME+RAND));

  // fill all spinor1 spinors with random numbers
  for( auto & elem : spinor1 ) rng.gen_z2(elem.data(),24);
  
  sw.reset();

  for( int outer_iter = 0; outer_iter < outer_iterations; ++outer_iter) {
    Momentum_Smearing(gauge_field_smeared,
                      mom,
                      0.8 /* momentum factor for mom smearing phase */,
                      spinor1[outer_iter].data(),
                      inner_iterations,
                      0.25);
  }

  char msg[400];
  snprintf(msg, 400, "%d x %d iterations of momentum smearing:",
           outer_iterations, inner_iterations);
  duration bench_dur = sw.elapsed_print(msg);

  debug_printf(
    0,0,
    "At a cost of 432 flops (per site) to prepare the gauge field and 990 flops (per site)\n"
    "per Jacobi iteration, this corresponds to %lf GFlop/s\n",
    ( ((double)432.0)*g_nproc*VOLUME*1.0e-9 + /* flops to prepare the mom_phase gauge field  */ 
    ( (double)990.0)*g_nproc*VOLUME*outer_iterations*inner_iterations*1.0e-9) /* flops for the Jacobi iterations */
      / bench_dur.mean
  );

  free(gauge_field_smeared);

  return 0;
}
