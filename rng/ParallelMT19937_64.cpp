#include "index_tools.hpp"
#include "cvc_global.h"
#include "ParallelMT19937_64.hpp"
#include "debug_printf.hpp"
#include "SequenceOfUnique.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "loop_tools.h"

#include <unistd.h>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <exception>
#include <stdexcept>

namespace cvc {
  // 311ULL is just an offset, we could use any other
  // It's the seed which decides where in the sequence we start
  ParallelMT19937_64::ParallelMT19937_64(const unsigned long long seed) :
    seed_gen(seed, 311ULL)
  { 
    const unsigned long long local_volume = VOLUME;
    const unsigned long long global_volume = VOLUME*g_nproc;
    local_seeds.resize( VOLUME );
    local_rngs.resize( VOLUME );
    
    unsigned long long ran64;
    for(unsigned long long gi = 0; gi < global_volume; ++gi){
      // use the sequence generator to produce our seeds for the global lattice
      ran64 = seed_gen.next();
      if( is_local(gi) ){
        local_seeds[ global_to_local_site_index(gi) ] = ran64;
      }
    }
    #pragma omp parallel for
    for(unsigned int li = 0; li < VOLUME; ++li){
      local_rngs[li].init(local_seeds[li]);
    }
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }

  void
  ParallelMT19937_64::gen_z2(double * const buffer, const unsigned int n_per_site)
  {
    const double one_ov_sqrt_2 = 1.0/sqrt(2.0);
    PARALLEL_AND_FOR(li, 0, VOLUME)
    {
      for(unsigned int i_per_site = 0; i_per_site < n_per_site; ++i_per_site){
        buffer[n_per_site*li + i_per_site] = local_rngs[li].gen_real() >= 0.5 ?
          one_ov_sqrt_2 : -one_ov_sqrt_2;
      }
    }
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }

  void
  ParallelMT19937_64::gen_gauss(double * const buffer, const unsigned int n_per_site)
  {
    if( n_per_site % 2 != 0 ){
      throw( std::invalid_argument("[ParallelMT19937_64::gen_gauss] n_per_site must be even!") );
    }
    const double TWO_MPI = 2.0 * M_PI;
    const unsigned int half_n_per_site = n_per_site / 2;

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
      double x1, x2, norm;
      FOR_IN_PARALLEL(li, 0, VOLUME)
      {
        for(unsigned int i_per_site = 0; i_per_site < half_n_per_site; ++i_per_site){
          x1 = sqrt( -log(1.0 - local_rngs[li].gen_real()) );
          x2 = TWO_MPI * local_rngs[li].gen_real();

          // real and imaginary components are in consecutive elements
          buffer[n_per_site*li + 2*i_per_site  ] = x1 * sin(x2);
          buffer[n_per_site*li + 2*i_per_site+1] = x1 * cos(x2);
 
          // our aim is to generate RNs such that v[i]*conj(v[i]) = 1.0 -> need to normalise
          norm = 1.0 / sqrt( buffer[n_per_site*li + 2*i_per_site  ]*buffer[n_per_site*li + 2*i_per_site  ] +
                             buffer[n_per_site*li + 2*i_per_site+1]*buffer[n_per_site*li + 2*i_per_site+1] );

          buffer[n_per_site*li + 2*i_per_site  ] *= norm;
          buffer[n_per_site*li + 2*i_per_site+1] *= norm;
        }
      }
    } // omp parallel closing brace
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }

  void
  ParallelMT19937_64::gen_test(double * const buffer, const unsigned int n_per_site)
  {
    PARALLEL_AND_FOR(li, 0, VOLUME)
    {
      for(unsigned int i_per_site = 0; i_per_site < n_per_site; ++i_per_site){
        buffer[n_per_site*li + i_per_site] = (double)li;
      }
    }
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }

  void
  ParallelMT19937_64::gen_real(double * const buffer, const unsigned int n_per_site) 
  {
    PARALLEL_AND_FOR(li, 0, VOLUME)
    {
      for(unsigned int i_per_site = 0; i_per_site < n_per_site; ++i_per_site){
        buffer[n_per_site*li + i_per_site] = local_rngs[li].gen_real();
      }
    }
#ifdef HAVE_MPI
    MPI_Barrier(g_cart_grid);
#endif
  }

  double
  ParallelMT19937_64::gen_real_at_site(const size_t local_site_idx){
    return(local_rngs[local_site_idx].gen_real());
  }

  unsigned long long
  ParallelMT19937_64::get_local_seed(const size_t local_site_idx){
    return(local_seeds[local_site_idx]);
  }

} // namespace(cvc)
