#include "global.h"
#include "loop_tools.h"
#include "cvc_complex.h"
#include "cvc_linalg.h"

#include <vector>
#include <cstring>
#include <iostream>

namespace cvc {

void contract_twopoint_all_gammas_ext_momentum_list(
    double * const contr, 
    const std::vector< std::vector<::cvc::complex> > & phase_fields,
    double const * const chi, double const * const psi) 
{
  size_t const VOL3 = LX*LY*LZ;
  unsigned int const n_momenta = phase_fields.size();
  size_t const n_elem = 2*T*16*n_momenta;

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // thread-local temporary for accumulation of momentum projected correlator
    std::vector<double> contr_tmp(n_elem, 0.0);

    // pointer into data of contr_tmp to make sure that accumulation
    // into "contr" below can be done atomically, although I'm not sure
    // it's even neessary
    double * contr_tmp_ptr = contr_tmp.data();

    size_t iix, ix3d, o_idx;
    unsigned int x0;
    double spinor1[24];

    ::cvc::complex w1, w2;

    // perform sink gamma multiplication, contraction and momentum projection
    FOR_IN_PARALLEL(ix, 0, VOLUME)
    {
      iix = _GSI(ix);
      ix3d = ix % VOL3;
      x0 = g_lexic2coords[ix][0];

      for(unsigned int i_mom = 0; i_mom < n_momenta; ++i_mom){
        ::cvc::complex const * const phases = phase_fields[i_mom].data();
        for(unsigned int i_gamma = 0; i_gamma < 16; ++i_gamma){
          // t is the slowest running index here
          o_idx = 2*(i_gamma + 16*(i_mom + n_momenta*x0));
          
          _fv_eq_gamma_ti_fv(spinor1, i_gamma, psi+iix);
          _co_eq_fv_dag_ti_fv(&w1, chi+iix, spinor1);
          w2.re = w1.re * phases[ix3d].re - w1.im * phases[ix3d].im;
          w2.im = w1.re * phases[ix3d].im + w1.im * phases[ix3d].re;
         
          // no locking or atomicity required here, contr_tmp is thread-local
          contr_tmp_ptr[o_idx    ] += w2.re;
          contr_tmp_ptr[o_idx + 1] += w2.im;
        }
      }
    }

    // accumulate results from all threads 
    // atomic write should be faster than having a lock
    for(size_t idx = 0; idx < n_elem; ++idx){
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
      contr[idx] += contr_tmp_ptr[idx];
    }
  } // end of parallel section if OpenMP in use

#ifdef HAVE_MPI
#if (defined PARALLELTX) || (defined PARALLELTXY) || (defined PARALLELTXYZ)
  std::vector<double> mpi_buffer(n_elem);
  memcpy(mpi_buffer.data(), contr, n_elem*sizeof(double) );
  if( MPI_Reduce(mpi_buffer.data(), contr, n_elem, MPI_DOUBLE, MPI_SUM, 0, g_ts_comm ) != MPI_SUCCESS ) {
    fprintf( stderr, "[contract_twopoint_all_gammas_ext_momentum_list] Error from MPI_Reduce %s %d\n", __FILE__, __LINE__);
  }
#endif 
#endif

}

} // namespace(cvc)
