
#define MAIN_PROGRAM
#include "global.h"
#undef MAIN_PROGRAM

#include "ParallelMT19937_64.hpp"
#include "enums.hpp"
#include "index_tools.hpp"
#include "Stopwatch.hpp"
#include "Core.hpp"
#include "debug_printf.hpp"
#include "loop_tools.h"
#include "propagator_io.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <unistd.h>


constexpr unsigned int startsite = 0;
constexpr unsigned int no_testsites = 1000;
constexpr unsigned int no_testvals = 50000;

using namespace cvc;

int main(int argc, char** argv)
{
  Core core(argc,argv);
  if( !(core.is_initialised()) ){
    std::cout << "Core initialisation failed!\n";
    return(CVC_EXIT_CORE_INIT_FAILURE);
  } 

#ifdef HAVE_MPI
  Stopwatch sw(g_cart_grid);
#else
  Stopwatch sw(0);
#endif

  debug_printf(0, 0, "\n\n############ TESTING ParallelMT19937_64 ###############\n\n");

  sw.reset();
  ParallelMT19937_64 rangen(982932ULL);
  sw.elapsed_print("ParallelMT19937_64 initialisation");

  std::vector< double > ranspinor(24*VOLUME);
  // we we also convert the seeds into double and write these out
  std::vector< double > local_seeds(24*VOLUME);

  for(unsigned long long i = 0; i < VOLUME; ++i){
    local_seeds[24*i] = static_cast<double>(local_to_global_site_index(i));
  }

  //// generate some z2 cross z2 random numbers
  sw.reset();
  rangen.gen_z2(ranspinor.data(),
                24);
  sw.elapsed_print("gen_z2");

#ifdef HAVE_MPI
  MPI_Barrier(g_cart_grid);
#endif

  char filename[100];
  snprintf(filename, 100,
           "local_seeds.npt%d_npx%d_npy%d_npz%d.lime",
           g_nproc_t, g_nproc_x, g_nproc_y, g_nproc_z);
  
  int exitstatus = 0;
  CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      write_propagator(local_seeds.data(), filename, 0, 64),
      "error in write_propagator",
      true,
      1);
  
  snprintf(filename, 100,
           "ranspinor.npt%d_npx%d_npy%d_npz%d.lime",
           g_nproc_t, g_nproc_x, g_nproc_y, g_nproc_z);
  
  CHECK_EXITSTATUS_NONZERO(
      exitstatus,
      write_propagator(ranspinor.data(), filename, 0, 64),
      "error in write_propagator",
      true,
      1);
  
  return 0;
}
