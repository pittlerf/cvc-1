
#define MAIN_PROGRAM
#include "cvc_global.h"
#undef MAIN_PROGRAM

#include "index_tools.hpp"
#include "Stopwatch.hpp"
#include "Core.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <vector>
#include <random>
#include <cstdio>

using namespace cvc;

int main(int argc, char** argv)
{
  Core core(argc,argv);
  if( !(core.is_initialised()) ){
    std::cout << "Core initialisation failed!\n";
    return(CVC_EXIT_CORE_INIT_FAILURE);
  } 

  std::default_random_engine generator;
  generator.seed(12345);
  std::uniform_int_distribution<unsigned long long> distribution(0,g_nproc*VOLUME);

  unsigned long long gi;
  for(int i = 0; i < 100; ++i){
    gi = distribution(generator);
    for( int proc_id = 0; proc_id < g_nproc; ++proc_id){
      if( g_proc_id == proc_id ){
        std::array<int,4> proc_coords;
        global_site_get_proc_coords(gi, proc_coords);
        std::string local_string( is_local(gi) ? "true" : "false" );
        std::cout << "Global index " << gi << " is local:" << local_string <<
          " on rank " << g_proc_id << " and should be on rank (t,x,y,z) " << 
         proc_coords[DIM_T] << " " << proc_coords[DIM_X] << " " <<
         proc_coords[DIM_Y] << " " << proc_coords[DIM_Z] << std::endl;
        fflush(stdout);
      }
      MPI_Barrier(g_cart_grid);
    }
    MPI_Barrier(g_cart_grid);
    if( g_proc_id == 0 ){
      std::cout << std::endl;
      fflush(stdout);
    }
  }

  return 0;
}
