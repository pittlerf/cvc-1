#include "global.h"
#include "loop_tools.h"
#include "cvc_complex.h"
#include "cvc_linalg.h"

#include <vector>
#include <cstring>

namespace cvc {

void contract_twopoint_gamma5_gamma_only_ext_momentum(
    double * const contr, const int gamma_id, ::cvc::complex const * const phases,
    double const * chi, double const * psi) 
{
  size_t const VOL3 = LX*LY*LZ;
  // for safety, we can set 'contr' to zero
  memset((void*)contr, 0, sizeof(double)*2*T);

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // thread-local temporary for accumulation of momentum projected correlator
    std::vector<double> contr_tmp(2*T, 0.0);

    // pointer into data of contr_tmp to make sure that accumulation
    // into "contr" below can be done atomically, although I'm not sure
    // it's even neessary
    double * contr_tmp_ptr = contr_tmp.data();

    size_t iix, ix3d;
    unsigned int x0;
    double spinor1[24], spinor2[24];

    ::cvc::complex w1, w2;

    // perform sink gamma multiplication, contraction and momentum projection
    FOR_IN_PARALLEL(ix, 0, VOLUME)
    {
      iix = _GSI(ix);
      ix3d = ix % VOL3;
      x0 = g_lexic2coords[ix][0];

      _fv_eq_gamma_ti_fv(spinor1, gamma_id, psi+iix);
      _fv_eq_gamma_ti_fv(spinor2, 5, spinor1);
      _co_eq_fv_dag_ti_fv(&w1, chi+iix, spinor2);
      w2.re = w1.re * phases[ix3d].re - w1.im * phases[ix3d].im;
      w2.im = w1.re * phases[ix3d].im + w1.im * phases[ix3d].re;
      
      // no locking or atomicity required here, contr_tmp is thread-local
      contr_tmp[2*x0  ] += w2.re;
      contr_tmp[2*x0+1] += w2.im;
    }

    // accumulate results from all threads 
    // atomic write should be faster than having a lock
    for(int t = 0; t < 2*T; t++){
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
      contr[t] += contr_tmp_ptr[t];
    }
  } // end of parallel section if OpenMP in use

#ifdef HAVE_MPI
#if (defined PARALLELTX) || (defined PARALLELTXY) || (defined PARALLELTXYZ)
  std::vector<double> mpi_buffer(2*T);
  memcpy(mpi_buffer.data(), contr, 2*T*sizeof(double) );
  if( MPI_Reduce(mpi_buffer.data(), contr, 2*T, MPI_DOUBLE, MPI_SUM, 0, g_ts_comm ) != MPI_SUCCESS ) {
    fprintf( stderr, "[contract_twopoint_gamma5_gamma_only_ext_momentum] Error from MPI_Reduce %s %d\n", __FILE__, __LINE__);
  }
#endif 
#endif

}

} // namespace(cvc)
