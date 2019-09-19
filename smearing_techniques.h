/**********************************
 * smearing_techniques.h  
 *
 * originally:
 * Author: Marc Wagner
 * Date: September 2007
 *
 * February 2010
 * taken over to cvc_parallel, parallelized
 *
 * Fri Dec  9 09:24:13 CET 2016
 * now copied from cvc_parallel_5d
 *
 **********************************/
#ifndef _SMEARING_TECHNIQUES_H
#define _SMEARING_TECHNIQUES_H

#include "cvc_linalg.h"
#include "types.h"

namespace cvc {

int APE_Smearing(double * const smeared_gauge_field, double const APE_smearing_alpha, int const APE_smearing_niter);
int Jacobi_Smearing(double * const smeared_gauge_field, double * const psi, int const N, double const kappa);

/**
 * @brief Wrapper for Jacobi_Smearing to implement momentum smearing
 *
 * @param smeared_gauge_field smeared gauge field which will have momentum phases
 *                            applied and passed on to Jacobi_Smearing (note:
 *                            original will not be modified!)
 * @param momentum momentum vector in lattice units (will be multiplied by 2*pi/L_{x,y,z})
 * @param mom_scale_factor scale factor to turn momentum phase into smearing phase. This
 *          depends on the problem at hand. For mesons a value of ~ 0.8 was shown to be
 *          ideal while for Baryons a value around 0.45 gave the best results.
 * @param psi spinor field to be smeared
 * @param N number of Jacobi iterations
 * @param kappa strength of Jacobi smearing, normalisation is 1/(1+6*kappa)
 *
 * @return 
 */
int Momentum_Smearing(double const * const smeared_gauge_field,
    mom_t const & momentum, double const mom_scale_factor,
    double * const psi, int const N, double const kappa);

}
#endif
