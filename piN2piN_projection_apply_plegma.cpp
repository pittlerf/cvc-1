/***************************************************
 * piN2piN_projection
 * 
 * PURPOSE:
 *   originally copied from piN2piN_correlators.cpp
 * TODO:
 * DONE:
 *
 ****************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef HAVE_MPI
#  include <mpi.h>
#endif
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <getopt.h>

#ifdef HAVE_LHPC_AFF
#include "lhpc-aff.h"
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef __cplusplus
}
#endif

#define MAIN_PROGRAM
#include <hdf5.h>


#include "cvc_complex.h"
#include "ilinalg.h"
#include "icontract.h"
#include "global.h"
#include "cvc_geometry.h"
#include "cvc_utils.h"
#include "cvc_timer.h"
#include "mpi_init.h"
#include "io.h"
#include "read_input_parser.h"
#include "contractions_io.h"
#include "matrix_init.h"
#include "table_init_z.h"
#include "table_init_2pt.h"
#include "table_init_d.h"
#include "table_init_i.h"
#include "contract_diagrams.h"
#include "aff_key_conversion.h"
#include "zm4x4.h"
#include "gamma.h"
#include "twopoint_function_utils.h"
#include "rotations.h"
#include "group_projection.h"
#include "little_group_projector_set.h"

#define MAX_UDLI_NUM 100000
#define EPS 1e-14
#define ANNIHILATION 1 
#define CREATION 0
using namespace cvc;


/* ************************************************************************
 *
 * routine that multiplies with the external gamma5 structure at the sink
 * used for gammas like C, Cg4, cg1g4g5 etc.
 *
 * ************************************************************************/

static inline void mult_with_gamma5_matrix_adj_source ( double ** buffer_write ) {

    double **buffer_temporary=init_2level_dtable(16,2);

    buffer_temporary[0][0]=-1*buffer_write[2][0]; buffer_temporary[0][1]=-1*buffer_write[2][1];
    buffer_temporary[1][0]=-1*buffer_write[3][0]; buffer_temporary[1][1]=-1*buffer_write[3][1];
    buffer_temporary[2][0]=-1*buffer_write[0][0]; buffer_temporary[2][1]=-1*buffer_write[0][1];
    buffer_temporary[3][0]=-1*buffer_write[1][0]; buffer_temporary[3][1]=-1*buffer_write[1][1];

    buffer_temporary[4][0]=-1*buffer_write[6][0]; buffer_temporary[4][1]=-1*buffer_write[6][1];
    buffer_temporary[5][0]=-1*buffer_write[7][0]; buffer_temporary[5][1]=-1*buffer_write[7][1];
    buffer_temporary[6][0]=-1*buffer_write[4][0]; buffer_temporary[6][1]=-1*buffer_write[4][1];
    buffer_temporary[7][0]=-1*buffer_write[5][0]; buffer_temporary[7][1]=-1*buffer_write[5][1];

    buffer_temporary[8][0] =-1*buffer_write[10][0]; buffer_temporary[8][1] =-1*buffer_write[10][1];
    buffer_temporary[9][0] =-1*buffer_write[11][0]; buffer_temporary[9][1] =-1*buffer_write[11][1];
    buffer_temporary[10][0]=-1*buffer_write[ 8][0]; buffer_temporary[10][1]=-1*buffer_write[ 8][1];
    buffer_temporary[11][0]=-1*buffer_write[ 9][0]; buffer_temporary[11][1]=-1*buffer_write[ 9][1];

    buffer_temporary[12][0]=-1*buffer_write[14][0]; buffer_temporary[12][1]=-1*buffer_write[14][1];
    buffer_temporary[13][0]=-1*buffer_write[15][0]; buffer_temporary[13][1]=-1*buffer_write[15][1];
    buffer_temporary[14][0]=-1*buffer_write[12][0]; buffer_temporary[14][1]=-1*buffer_write[12][1];
    buffer_temporary[15][0]=-1*buffer_write[13][0]; buffer_temporary[15][1]=-1*buffer_write[13][1];

    for (int i=0; i<16; ++i){
      buffer_write[i][0]=buffer_temporary[i][0];
    }

    fini_2level_dtable(&buffer_temporary);
}



static inline void mult_with_gamma5_matrix_sink ( double ** buffer_write ) {

    double **buffer_temporary=init_2level_dtable(16,2);

    buffer_temporary[0][0]=buffer_write[8][0]; buffer_temporary[0][1]=buffer_write[8][1];
    buffer_temporary[1][0]=buffer_write[9][0]; buffer_temporary[1][1]=buffer_write[9][1];
    buffer_temporary[2][0]=buffer_write[10][0]; buffer_temporary[2][1]=buffer_write[10][1];
    buffer_temporary[3][0]=buffer_write[11][0]; buffer_temporary[3][1]=buffer_write[11][1];

    buffer_temporary[4][0]=buffer_write[12][0]; buffer_temporary[4][1]=buffer_write[12][1];
    buffer_temporary[5][0]=buffer_write[13][0]; buffer_temporary[5][1]=buffer_write[13][1];
    buffer_temporary[6][0]=buffer_write[14][0]; buffer_temporary[6][1]=buffer_write[14][1];
    buffer_temporary[7][0]=buffer_write[15][0]; buffer_temporary[7][1]=buffer_write[15][1];

    buffer_temporary[8][0]=buffer_write[ 0][0]; buffer_temporary[8][1]=buffer_write[0][1];
    buffer_temporary[9][0]=buffer_write[ 1][0]; buffer_temporary[9][1]=buffer_write[1][1];
    buffer_temporary[10][0]=buffer_write[2][0]; buffer_temporary[10][1]=buffer_write[2][1];
    buffer_temporary[11][0]=buffer_write[3][0]; buffer_temporary[11][1]=buffer_write[3][1];

    buffer_temporary[12][0]=buffer_write[4][0]; buffer_temporary[12][1]=buffer_write[4][1];
    buffer_temporary[13][0]=buffer_write[5][0]; buffer_temporary[13][1]=buffer_write[5][1];
    buffer_temporary[14][0]=buffer_write[6][0]; buffer_temporary[14][1]=buffer_write[6][1];
    buffer_temporary[15][0]=buffer_write[7][0]; buffer_temporary[15][1]=buffer_write[7][1];

    for (int i=0; i<16; ++i){
      buffer_write[i][0]=buffer_temporary[i][0];
    }

    fini_2level_dtable(&buffer_temporary);
}





/* ************************************************************************
 *
 * routine that applies parity matrix  to the current
 * spin matrix 
 *
 * ************************************************************************/

static inline void mult_parity_matrix ( double ** buffer_target , double ** buffer_source ) {

    buffer_target[0][0]=   buffer_source[0][0];   buffer_target[0][1]=   buffer_source[0][1];
    buffer_target[1][0]=   buffer_source[1][0];   buffer_target[1][1]=   buffer_source[1][1];
    buffer_target[2][0]=-1*buffer_source[2][0];   buffer_target[2][1]=-1*buffer_source[2][1];
    buffer_target[3][0]=-1*buffer_source[3][0];   buffer_target[3][1]=-1*buffer_source[3][1];

    buffer_target[4][0]=   buffer_source[4][0];   buffer_target[4][1]=   buffer_source[4][1];
    buffer_target[5][0]=   buffer_source[5][0];   buffer_target[5][1]=   buffer_source[5][1];
    buffer_target[6][0]=-1*buffer_source[6][0];   buffer_target[6][1]=-1*buffer_source[6][1];
    buffer_target[7][0]=-1*buffer_source[7][0];   buffer_target[7][1]=-1*buffer_source[7][1];

    buffer_target[8][0] =-1*buffer_source[8][0];   buffer_target[8][1]= -1*buffer_source[8][1];
    buffer_target[9][0] =-1*buffer_source[9][0];   buffer_target[9][1]= -1*buffer_source[9][1];
    buffer_target[10][0]=   buffer_source[10][0];  buffer_target[10][1]=   buffer_source[10][1];
    buffer_target[11][0]=   buffer_source[11][0];  buffer_target[11][1]=   buffer_source[11][1];

    buffer_target[12][0]=-1*buffer_source[12][0];   buffer_target[12][1]=-1*buffer_source[12][1];
    buffer_target[13][0]=-1*buffer_source[13][0];   buffer_target[13][1]=-1*buffer_source[13][1];
    buffer_target[14][0]=   buffer_source[14][0];   buffer_target[14][1]=   buffer_source[14][1];
    buffer_target[15][0]=   buffer_source[15][0];   buffer_target[15][1]=   buffer_source[15][1];

}


/* ************************************************************************
 *
 * routine that applies time reversal matrix  to the current
 * spin matrix 
 *
 * ************************************************************************/

static inline void mult_time_reversal_matrix ( double ** buffer_target , double ** buffer_source ) {

    buffer_target[0][0]=   buffer_source[10][0];  buffer_target[0][1]=   buffer_source[10][1];
    buffer_target[1][0]=   buffer_source[11][0];  buffer_target[1][1]=   buffer_source[11][1];
    buffer_target[2][0]=-1*buffer_source[8][0];   buffer_target[2][1]=-1*buffer_source[8][1];
    buffer_target[3][0]=-1*buffer_source[9][0];   buffer_target[3][1]=-1*buffer_source[9][1];

    buffer_target[4][0]=   buffer_source[14][0];   buffer_target[4][1]=   buffer_source[14][1];
    buffer_target[5][0]=   buffer_source[15][0];   buffer_target[5][1]=   buffer_source[15][1];
    buffer_target[6][0]=-1*buffer_source[12][0];   buffer_target[6][1]=-1*buffer_source[12][1];
    buffer_target[7][0]=-1*buffer_source[13][0];   buffer_target[7][1]=-1*buffer_source[13][1];

    buffer_target[8][0] =-1*buffer_source[2][0];   buffer_target[8][1]= -1*buffer_source[2][1];
    buffer_target[9][0] =-1*buffer_source[3][0];   buffer_target[9][1]= -1*buffer_source[3][1];
    buffer_target[10][0]=   buffer_source[0][0];  buffer_target[10][1]=    buffer_source[0][1];
    buffer_target[11][0]=   buffer_source[1][0];  buffer_target[11][1]=    buffer_source[1][1];

    buffer_target[12][0]=-1*buffer_source[12][0];   buffer_target[12][1]=-1*buffer_source[6][1];
    buffer_target[13][0]=-1*buffer_source[13][0];   buffer_target[13][1]=-1*buffer_source[7][1];
    buffer_target[14][0]=   buffer_source[14][0];   buffer_target[14][1]=   buffer_source[4][1];
    buffer_target[15][0]=   buffer_source[15][0];   buffer_target[15][1]=   buffer_source[5][1];

}

/*************************************************************************
 *
 * routine that applies CP  matrix  to the current
 * spin matrix 
 *
 * ************************************************************************/

static inline void mult_cp_matrix ( double ** buffer_target , double ** buffer_source ) {

    buffer_target[0][0]=   buffer_source[15][0];   buffer_target[0][1]=   buffer_source[15][1];
    buffer_target[1][0]=-1*buffer_source[11][0];   buffer_target[1][1]=-1*buffer_source[11][1];
    buffer_target[2][0]=-1*buffer_source[ 7][0];   buffer_target[2][1]=-1*buffer_source[7][1];
    buffer_target[3][0]=   buffer_source[ 3][0];   buffer_target[3][1]=   buffer_source[3][1];

    buffer_target[4][0]=-1*buffer_source[14][0];   buffer_target[4][1]=-1*buffer_source[14][1];
    buffer_target[5][0]=   buffer_source[10][0];   buffer_target[5][1]=   buffer_source[10][1];
    buffer_target[6][0]=   buffer_source[ 6][0];   buffer_target[6][1]=   buffer_source[ 6][1];
    buffer_target[7][0]=-1*buffer_source[ 2][0];   buffer_target[7][1]=-1*buffer_source[ 2][1];

    buffer_target[8][0] =-1*buffer_source[13][0];   buffer_target[8][1]=-1*buffer_source[13][1];
    buffer_target[9][0] =   buffer_source[ 9][0];   buffer_target[9][1]=   buffer_source[ 9][1];
    buffer_target[10][0]=   buffer_source[ 5][0];  buffer_target[10][1]=   buffer_source[ 5][1];
    buffer_target[11][0]=-1*buffer_source[ 1][0];  buffer_target[11][1]=-1*buffer_source[ 1][1];

    buffer_target[12][0]=   buffer_source[12][0];   buffer_target[12][1]=   buffer_source[12][1];
    buffer_target[13][0]=-1*buffer_source[ 8][0];   buffer_target[13][1]=-1*buffer_source[ 8][1];
    buffer_target[14][0]=-1*buffer_source[ 4][0];   buffer_target[14][1]=-1*buffer_source[ 4][1];
    buffer_target[15][0]=   buffer_source[ 0][0];   buffer_target[15][1]=   buffer_source[ 0][1];

}


/*************************************************************************
 *
 * routine that applies CT  matrix  to the current
 * spin matrix 
 *
 * ************************************************************************/

static inline void mult_ct_matrix ( double ** buffer_target , double ** buffer_source ) {

    buffer_target[0][0]=   buffer_source[ 5][0];   buffer_target[0][1]=   buffer_source[ 5][1];
    buffer_target[1][0]=-1*buffer_source[ 1][0];   buffer_target[1][1]=-1*buffer_source[ 1][1];
    buffer_target[2][0]=-1*buffer_source[13][0];   buffer_target[2][1]=-1*buffer_source[13][1];
    buffer_target[3][0]=   buffer_source[ 9][0];   buffer_target[3][1]=   buffer_source[ 9][1];

    buffer_target[4][0]=-1*buffer_source[ 4][0];   buffer_target[4][1]=-1*buffer_source[ 4][1];
    buffer_target[5][0]=   buffer_source[ 0][0];   buffer_target[5][1]=   buffer_source[ 0][1];
    buffer_target[6][0]=   buffer_source[12][0];   buffer_target[6][1]=   buffer_source[12][1];
    buffer_target[7][0]=-1*buffer_source[ 8][0];   buffer_target[7][1]=-1*buffer_source[ 8][1];

    buffer_target[8][0] =-1*buffer_source[ 7][0];   buffer_target[8][1]=-1*buffer_source[ 7][1];
    buffer_target[9][0] =   buffer_source[ 3][0];   buffer_target[9][1]=   buffer_source[ 3][1];
    buffer_target[10][0]=   buffer_source[15][0];  buffer_target[10][1]=   buffer_source[15][1];
    buffer_target[11][0]=-1*buffer_source[11][0];  buffer_target[11][1]=-1*buffer_source[11][1];

    buffer_target[12][0]=   buffer_source[ 6][0];   buffer_target[12][1]=   buffer_source[ 6][1];
    buffer_target[13][0]=-1*buffer_source[ 2][0];   buffer_target[13][1]=-1*buffer_source[ 2][1];
    buffer_target[14][0]=-1*buffer_source[14][0];   buffer_target[14][1]=-1*buffer_source[14][1];
    buffer_target[15][0]=   buffer_source[10][0];   buffer_target[15][1]=   buffer_source[10][1];

}

static inline void mult_pt_matrix ( double ** buffer_target , double ** buffer_source ) {

    buffer_target[0][0]=   buffer_source[10][0];   buffer_target[0][1]=   buffer_source[10][1];
    buffer_target[1][0]=   buffer_source[11][0];   buffer_target[1][1]=   buffer_source[11][1];
    buffer_target[2][0]=   buffer_source[ 8][0];   buffer_target[2][1]=   buffer_source[ 8][1];
    buffer_target[3][0]=   buffer_source[ 9][0];   buffer_target[3][1]=   buffer_source[ 9][1];

    buffer_target[4][0]=   buffer_source[14][0];   buffer_target[4][1]=   buffer_source[14][1];
    buffer_target[5][0]=   buffer_source[15][0];   buffer_target[5][1]=   buffer_source[15][1];
    buffer_target[6][0]=   buffer_source[12][0];   buffer_target[6][1]=   buffer_source[12][1];
    buffer_target[7][0]=   buffer_source[13][0];   buffer_target[7][1]=   buffer_source[13][1];

    buffer_target[8][0] =   buffer_source[ 2][0];  buffer_target[8][1]=   buffer_source[ 2][1];
    buffer_target[9][0] =   buffer_source[ 3][0];  buffer_target[9][1]=   buffer_source[ 3][1];
    buffer_target[10][0]=   buffer_source[ 0][0];  buffer_target[10][1]=  buffer_source[ 0][1];
    buffer_target[11][0]=   buffer_source[ 1][0];  buffer_target[11][1]=  buffer_source[ 1][1];

    buffer_target[12][0]=   buffer_source[ 6][0];   buffer_target[12][1]=   buffer_source[ 6][1];
    buffer_target[13][0]=   buffer_source[ 7][0];   buffer_target[13][1]=   buffer_source[ 7][1];
    buffer_target[14][0]=   buffer_source[ 4][0];   buffer_target[14][1]=   buffer_source[ 4][1];
    buffer_target[15][0]=   buffer_source[ 5][0];   buffer_target[15][1]=   buffer_source[ 5][1];

}


/*************************************************************************
 *
 * routine that applies CPT  matrix  to the current
 * spin matrix 
 *
 * ************************************************************************/

static inline void mult_cpt_matrix ( double ** buffer_target , double ** buffer_source ) {

    buffer_target[0][0]=   buffer_source[ 5][0];   buffer_target[0][1]=   buffer_source[ 5][1];
    buffer_target[1][0]=-1*buffer_source[ 1][0];   buffer_target[1][1]=-1*buffer_source[ 1][1];
    buffer_target[2][0]=   buffer_source[13][0];   buffer_target[2][1]=   buffer_source[13][1];
    buffer_target[3][0]=-1*buffer_source[ 9][0];   buffer_target[3][1]=-1*buffer_source[ 9][1];

    buffer_target[4][0]=-1*buffer_source[ 4][0];   buffer_target[4][1]=-1*buffer_source[ 4][1];
    buffer_target[5][0]=   buffer_source[ 0][0];   buffer_target[5][1]=   buffer_source[ 0][1];
    buffer_target[6][0]=-1*buffer_source[12][0];   buffer_target[6][1]=-1*buffer_source[12][1];
    buffer_target[7][0]=   buffer_source[ 8][0];   buffer_target[7][1]=   buffer_source[ 8][1];

    buffer_target[8][0] =   buffer_source[ 7][0];  buffer_target[8][1]=    buffer_source[ 7][1];
    buffer_target[9][0] =-1*buffer_source[ 3][0];  buffer_target[9][1]= -1*buffer_source[ 3][1];
    buffer_target[10][0]=   buffer_source[15][0];  buffer_target[10][1]=   buffer_source[15][1];
    buffer_target[11][0]=-1*buffer_source[11][0];  buffer_target[11][1]=-1*buffer_source[11][1];

    buffer_target[12][0]=-1*buffer_source[ 6][0];   buffer_target[12][1]=-1*buffer_source[ 6][1];
    buffer_target[13][0]=   buffer_source[ 2][0];   buffer_target[13][1]=   buffer_source[ 2][1];
    buffer_target[14][0]=-1*buffer_source[14][0];   buffer_target[14][1]=-1*buffer_source[14][1];
    buffer_target[15][0]=   buffer_source[10][0];   buffer_target[15][1]=   buffer_source[10][1];

}


static inline void apply_signs_discrete_symmetry ( int *sign, char * source_gamma ) {
    if (strcmp(source_gamma,"cg1")){
      sign[0]=+1;
      sign[1]=+1;
      sign[2]=+1;
      sign[3]=+1;
      sign[4]=+1;
      sign[5]=+1;
      sign[6]=-1;
      sign[7]=-1;
      sign[8]=-1;
    }
    else if (strcmp(source_gamma,"cg2")){
      sign[0]=+1;
      sign[1]=-1;
      sign[2]=+1;
      sign[3]=-1;
      sign[4]=+1;
      sign[5]=-1;
      sign[6]=-1;
      sign[7]=+1;
      sign[8]=+1;
    }
    else if (strcmp(source_gamma,"cg3")){
      sign[0]=+1;
      sign[1]=+1;
      sign[2]=+1;
      sign[3]=+1;
      sign[4]=+1;
      sign[5]=+1;
      sign[6]=-1;
      sign[7]=-1;
      sign[8]=-1;
    }
    else if (strcmp(source_gamma,"cg1g4")){
      sign[0]=+1;
      sign[1]=-1;
      sign[2]=+1;
      sign[3]=+1;
      sign[4]=-1;
      sign[5]=-1;
      sign[6]=+1;
      sign[7]=+1;
      sign[8]=-1;
    }
    else if (strcmp(source_gamma,"cg2g4")){
      sign[0]=+1;
      sign[1]=+1;
      sign[2]=+1;
      sign[3]=-1;
      sign[4]=-1;
      sign[5]=+1;
      sign[6]=+1;
      sign[7]=-1;
      sign[8]=+1;
    }
    else if (strcmp(source_gamma,"cg3g4")){
      sign[0]=+1;
      sign[1]=-1;
      sign[2]=+1;
      sign[3]=+1;
      sign[4]=-1;
      sign[5]=-1;
      sign[6]=+1;
      sign[7]=+1;
      sign[8]=-1;
    }
    else if (strcmp(source_gamma,"cg1g4g5")){
      sign[0]=-1;
      sign[1]=-1;
      sign[2]=-1;
      sign[3]=+1;
      sign[4]=+1;
      sign[5]=+1;
      sign[6]=+1;
      sign[7]=-1;
      sign[8]=+1;
    }
    else if (strcmp(source_gamma,"cg2g4g5")){
      sign[0]=-1;
      sign[1]=+1;
      sign[2]=-1;
      sign[3]=+1;
      sign[4]=+1;
      sign[5]=-1;
      sign[6]=+1;
      sign[7]=+1;
      sign[8]=-1;
    }
    else if (strcmp(source_gamma,"cg3g4g5")){
      sign[0]=-1;
      sign[1]=-1;
      sign[2]=-1;
      sign[3]=+1;
      sign[4]=+1;
      sign[5]=+1;
      sign[6]=+1;
      sign[7]=-1;
      sign[8]=+1;
    }
    else if (strcmp(source_gamma,"Cg5")){
      sign[0]=+1;
      sign[1]=+1;
      sign[2]=+1;
      sign[3]=+1;
      sign[4]=+1;
      sign[5]=+1;
      sign[6]=+1;
      sign[7]=-1;
      sign[8]=-1;
    }
    else if (strcmp(source_gamma,"Cg5g4")){
      sign[0]=+1;
      sign[1]=+1;
      sign[2]=+1;
      sign[3]=+1;
      sign[4]=-1;
      sign[5]=-1;
      sign[6]=-1;
      sign[7]=-1;
      sign[8]=-1;

    }
    else if (strcmp(source_gamma,"Cg4")){
      sign[0]=+1;
      sign[1]=+1;
      sign[2]=-1;
      sign[3]=+1;
      sign[4]=+1;
      sign[5]=-1;
      sign[6]=-1;
      sign[7]=+1;
      sign[8]=-1;
    }
    else if (strcmp(source_gamma,"C")){
      sign[0]=+1;
      sign[1]=+1;
      sign[2]=-1;
      sign[3]=-1;
      sign[4]=-1;
      sign[5]=-1;
      sign[6]=+1;
      sign[7]=-1;
      sign[8]=+1;
    }
    else if (strcmp(source_gamma,"g5")){
      sign[0]=+1;
      sign[1]=+1;
      sign[2]=-1;
      sign[3]=-1;
      sign[4]=-1;
      sign[5]=-1;
      sign[6]=+1;
      sign[7]=+1;
      sign[8]=-1;
    }
    else {
      fprintf(stderr, "# [apply_signs_discrete_symmetry] Not recognized gamma in discrete symmetry signs\n");
      exit(1);
    }

}


/* ************************************************************************
 *
 * routine that applies charge conjugation matrix and transposition to the current
 * spin matrix 
 *
 * ************************************************************************/

static inline void mult_with_charge_conjugation ( double ** buffer_target , double ** buffer_source ) {

    buffer_target[0][0]=   buffer_source[15][0];  buffer_target[0][1]=   buffer_source[15][1];
    buffer_target[1][0]=-1*buffer_source[11][0];  buffer_target[1][1]=-1*buffer_source[11][1];
    buffer_target[2][0]=   buffer_source[7][0];   buffer_target[2][1]=   buffer_source[7][1];
    buffer_target[3][0]=-1*buffer_source[3][0];   buffer_target[3][1]=-1*buffer_source[3][1];

    buffer_target[4][0]=-1*buffer_source[14][0];  buffer_target[4][1]=-1*buffer_source[14][1];
    buffer_target[5][0]=   buffer_source[10][0];  buffer_target[5][1]=   buffer_source[10][1];
    buffer_target[6][0]=-1*buffer_source[6][0];   buffer_target[6][1]=   buffer_source[6][1];
    buffer_target[7][0]=   buffer_source[2][0];   buffer_target[7][1]=-1*buffer_source[2][1];

    buffer_target[8][0] =   buffer_source[13][0];  buffer_target[8][1]=    buffer_source[13][1];
    buffer_target[9][0] =-1*buffer_source[9][0];   buffer_target[9][1]= -1*buffer_source[9][1];
    buffer_target[10][0]=   buffer_source[5][0];   buffer_target[10][1]=   buffer_source[5][1];
    buffer_target[11][0]=-1*buffer_source[1][0];   buffer_target[11][1]=-1*buffer_source[1][1];

    buffer_target[12][0]=-1*buffer_source[12][0];  buffer_target[12][1]=-1*buffer_source[12][1];
    buffer_target[13][0]=   buffer_source[8][0];   buffer_target[13][1]=   buffer_source[8][1];
    buffer_target[14][0]=-1*buffer_source[4][0];   buffer_target[14][1]=   buffer_source[4][1];
    buffer_target[15][0]=   buffer_source[0][0];   buffer_target[15][1]=-1*buffer_source[0][1];

}



static inline void mult_with_c( double ***buffer_accum,  double *** buffer_write, twopoint_function_type * tp, char *gamma_string_source, char *gamma_string_sink ){
   int *sign_table_source=init_1level_itable(9);
   int *sign_table_sink=init_1level_itable(9);

   int *sign_table_gamma5=init_1level_itable(9);

   apply_signs_discrete_symmetry ( sign_table_source, gamma_string_source );
   apply_signs_discrete_symmetry ( sign_table_sink,   gamma_string_sink  );

   apply_signs_discrete_symmetry ( sign_table_gamma5,   "g5"  );

   double ***buffer_charge_conjugated= init_3level_dtable(tp->T,tp->d*tp->d,2);

   mult_with_charge_conjugation(buffer_charge_conjugated[0],buffer_write[0]);

   for (int spin_source=0; spin_source<tp->d; ++spin_source) {
     for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {
       /* The minus sign comes from the last column in Marcus's notes */
       buffer_charge_conjugated[0][spin_sink*4+spin_source][0]*=sign_table_source[1]*sign_table_sink[1]*sign_table_source[3]*sign_table_sink[3]*-1;
       buffer_charge_conjugated[0][spin_sink*4+spin_source][1]*=sign_table_source[1]*sign_table_sink[1]*sign_table_source[3]*sign_table_sink[3]*-1;
       if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
          buffer_charge_conjugated[0][spin_sink*4+spin_source][0]*=sign_table_gamma5[1]*sign_table_gamma5[3];
          buffer_charge_conjugated[0][spin_sink*4+spin_source][1]*=sign_table_gamma5[1]*sign_table_gamma5[3];

       }
       if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
          buffer_charge_conjugated[0][spin_sink*4+spin_source][0]*=sign_table_gamma5[1]*sign_table_gamma5[3];
          buffer_charge_conjugated[0][spin_sink*4+spin_source][1]*=sign_table_gamma5[1]*sign_table_gamma5[3];
       }

       buffer_accum[0][spin_sink*4+spin_source][0]=buffer_charge_conjugated[0][spin_sink*4+spin_source][0];
       buffer_accum[0][spin_sink*4+spin_source][1]=buffer_charge_conjugated[0][spin_sink*4+spin_source][1];
     }
   }
   for (int time_extent = 1; time_extent < tp->T ; ++ time_extent ){
     mult_with_charge_conjugation(buffer_charge_conjugated[time_extent],buffer_write[tp->T-time_extent]);
     for (int spin_source=0; spin_source<tp->d; ++spin_source) {
       for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {

         buffer_charge_conjugated[time_extent][spin_sink*4+spin_source][0]*=sign_table_source[1]*sign_table_sink[1]*sign_table_source[3]*sign_table_sink[3]*-1;
         buffer_charge_conjugated[time_extent][spin_sink*4+spin_source][1]*=sign_table_source[1]*sign_table_sink[1]*sign_table_source[3]*sign_table_sink[3]*-1;

         if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
             buffer_charge_conjugated[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[1]*sign_table_gamma5[3];
             buffer_charge_conjugated[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[1]*sign_table_gamma5[3];

         }
         if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
             buffer_charge_conjugated[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[1]*sign_table_gamma5[3];
             buffer_charge_conjugated[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[1]*sign_table_gamma5[3];
         }
         buffer_accum[time_extent][spin_sink*4+spin_source][0]=buffer_charge_conjugated[time_extent][spin_sink*4+spin_source][0];
         buffer_accum[time_extent][spin_sink*4+spin_source][1]=buffer_charge_conjugated[time_extent][spin_sink*4+spin_source][1];

       }
     }
   }

   fini_3level_dtable(&buffer_charge_conjugated);
   free(sign_table_source);
   free(sign_table_sink);
   free(sign_table_gamma5);


}


static inline void mult_with_ct( double ***buffer_accum,  double *** buffer_write, twopoint_function_type * tp, char *gamma_string_source, char *gamma_string_sink ){


    /* CT */
    /* tagname suffix  is CT, abbreviation from charge conjugation + time reversal */
   int *sign_table_source=init_1level_itable(9);
   int *sign_table_sink=init_1level_itable(9);

   int *sign_table_gamma5=init_1level_itable(9);

   apply_signs_discrete_symmetry ( sign_table_source, gamma_string_source );
   apply_signs_discrete_symmetry ( sign_table_sink,   gamma_string_sink  );

   apply_signs_discrete_symmetry ( sign_table_gamma5,   "g5"  );

   double ***buffer_CT= init_3level_dtable(tp->T,tp->d*tp->d,2);

   for (int time_extent = 0; time_extent < tp->T ; ++time_extent ){

     mult_ct_matrix(buffer_CT[time_extent],buffer_write[time_extent]);
     for (int spin_source=0; spin_source<tp->d; ++spin_source) {
       for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {

         buffer_CT[time_extent][spin_sink*4+spin_source][0]*=sign_table_source[6]*sign_table_sink[6]*sign_table_source[3]*sign_table_sink[3];
         buffer_CT[time_extent][spin_sink*4+spin_source][1]*=sign_table_source[6]*sign_table_sink[6]*sign_table_source[3]*sign_table_sink[3];

         if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
             buffer_CT[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[6]*sign_table_gamma5[3];
             buffer_CT[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[6]*sign_table_gamma5[3];

         }
         if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
             buffer_CT[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[6]*sign_table_gamma5[3];
             buffer_CT[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[6]*sign_table_gamma5[3];
         }

         buffer_accum[time_extent][spin_sink*4+spin_source][0]=buffer_CT[time_extent][spin_sink*4+spin_source][0];
         buffer_accum[time_extent][spin_sink*4+spin_source][1]=buffer_CT[time_extent][spin_sink*4+spin_source][1];


       }
     }
   }

   fini_3level_dtable(&buffer_CT);
   free(sign_table_source);
   free(sign_table_sink);
   free(sign_table_gamma5);

}

static inline void mult_with_p( double ***buffer_accum,  double *** buffer_write, twopoint_function_type * tp, char *gamma_string_source, char *gamma_string_sink ){

   /* parity */
   /* tagname suffix  is P, abbreviation from parity */

   int *sign_table_source=init_1level_itable(9);
   int *sign_table_sink=init_1level_itable(9);

   int *sign_table_gamma5=init_1level_itable(9);

   apply_signs_discrete_symmetry ( sign_table_source, gamma_string_source );
   apply_signs_discrete_symmetry ( sign_table_sink,   gamma_string_sink  );

   apply_signs_discrete_symmetry ( sign_table_gamma5,   "g5"  );

   double ***buffer_P= init_3level_dtable(tp->T,tp->d*tp->d,2);

   for (int time_extent = 0; time_extent < tp->T ; ++ time_extent ){

     mult_parity_matrix(buffer_P[time_extent],buffer_write[time_extent]);

     for (int spin_source=0; spin_source<tp->d; ++spin_source) {
       for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {

         buffer_P[time_extent][spin_sink*4+spin_source][0]*=sign_table_source[2]*sign_table_sink[2];
         buffer_P[time_extent][spin_sink*4+spin_source][1]*=sign_table_source[2]*sign_table_sink[2];

         if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
            buffer_P[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[2];
            buffer_P[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[2];

         }
         if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
            buffer_P[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[2];
            buffer_P[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[2];
         }

         buffer_accum[time_extent][spin_sink*4+spin_source][0]=buffer_P[time_extent][spin_sink*4+spin_source][0];
         buffer_accum[time_extent][spin_sink*4+spin_source][1]=buffer_P[time_extent][spin_sink*4+spin_source][1];

       }
     }
   }

   fini_3level_dtable(&buffer_P);
   free(sign_table_source);
   free(sign_table_sink);
   free(sign_table_gamma5);

}

static inline void mult_with_pt( double ***buffer_accum,  double *** buffer_write, twopoint_function_type * tp, char *gamma_string_source, char *gamma_string_sink ){


   /* PT */
   /* tagname suffix  is PT, abbreviation from  parity plus time-reversal */
   int *sign_table_source=init_1level_itable(9);
   int *sign_table_sink=init_1level_itable(9);

   int *sign_table_gamma5=init_1level_itable(9);

   apply_signs_discrete_symmetry ( sign_table_source, gamma_string_source );
   apply_signs_discrete_symmetry ( sign_table_sink,   gamma_string_sink  );

   apply_signs_discrete_symmetry ( sign_table_gamma5,   "g5"  );

   double ***buffer_PT= init_3level_dtable(tp->T,tp->d*tp->d,2);

   mult_pt_matrix(buffer_PT[0],buffer_write[0]);

   for (int spin_source=0; spin_source<tp->d; ++spin_source) {
     for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {

       /* The minus sign comes from Marcus's notes last column */
       buffer_PT[0][spin_sink*4+spin_source][0]*=sign_table_source[7]*sign_table_sink[7]*-1;
       buffer_PT[0][spin_sink*4+spin_source][1]*=sign_table_source[7]*sign_table_sink[7]*-1;
       if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
         buffer_PT[0][spin_sink*4+spin_source][0]*=sign_table_gamma5[7];
         buffer_PT[0][spin_sink*4+spin_source][1]*=sign_table_gamma5[7];

       }
       if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
          buffer_PT[0][spin_sink*4+spin_source][0]*=sign_table_gamma5[7];
          buffer_PT[0][spin_sink*4+spin_source][1]*=sign_table_gamma5[7];
       }
       buffer_accum[0][spin_sink*4+spin_source][0]=buffer_PT[0][spin_sink*4+spin_source][0];
       buffer_accum[0][spin_sink*4+spin_source][1]=buffer_PT[0][spin_sink*4+spin_source][1];

     }
   }
   for (int time_extent = 1; time_extent < tp->T ; ++ time_extent ){

     mult_pt_matrix(buffer_PT[time_extent],buffer_write[tp->T - time_extent]);
     for (int spin_source=0; spin_source<tp->d; ++spin_source) {
       for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {

         buffer_PT[time_extent][spin_sink*4+spin_source][0]*=sign_table_source[7]*sign_table_sink[7]*-1;
         buffer_PT[time_extent][spin_sink*4+spin_source][1]*=sign_table_source[7]*sign_table_sink[7]*-1;

         if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
           buffer_PT[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[7];
           buffer_PT[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[7];

         }
         if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
           buffer_PT[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[7];
           buffer_PT[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[7];
         }

         buffer_accum[time_extent][spin_sink*4+spin_source][0]=buffer_PT[time_extent][spin_sink*4+spin_source][0];
         buffer_accum[time_extent][spin_sink*4+spin_source][1]=buffer_PT[time_extent][spin_sink*4+spin_source][1];


       }
     }
   }

   fini_3level_dtable(&buffer_PT);
   free(sign_table_source);
   free(sign_table_sink);
   free(sign_table_gamma5);

}

static inline void mult_with_cp( double ***buffer_accum,  double *** buffer_write, twopoint_function_type * tp, char *gamma_string_source, char *gamma_string_sink ){


   double ***buffer_CP= init_3level_dtable(tp->T,tp->d*tp->d,2);
   int *sign_table_source=init_1level_itable(9);
   int *sign_table_sink=init_1level_itable(9);

   int *sign_table_gamma5=init_1level_itable(9);

   apply_signs_discrete_symmetry ( sign_table_source, gamma_string_source );
   apply_signs_discrete_symmetry ( sign_table_sink,   gamma_string_sink  );

   apply_signs_discrete_symmetry ( sign_table_gamma5,   "g5"  );


   mult_cp_matrix(buffer_CP[0],buffer_write[0]);

   for (int spin_source=0; spin_source<tp->d; ++spin_source) {
     for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {
       /* minus sign here is the last column from Marcus's note */
       buffer_CP[0][spin_sink*4+spin_source][0]*=sign_table_source[5]*sign_table_sink[5]*sign_table_source[3]*sign_table_sink[3]*-1;
       buffer_CP[0][spin_sink*4+spin_source][1]*=sign_table_source[5]*sign_table_sink[5]*sign_table_source[3]*sign_table_sink[3]*-1;
       if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
         buffer_CP[0][spin_sink*4+spin_source][0]*=sign_table_gamma5[5]*sign_table_gamma5[3];
         buffer_CP[0][spin_sink*4+spin_source][1]*=sign_table_gamma5[5]*sign_table_gamma5[3];

       }
       if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
         buffer_CP[0][spin_sink*4+spin_source][0]*=sign_table_gamma5[5]*sign_table_gamma5[3];
         buffer_CP[0][spin_sink*4+spin_source][1]*=sign_table_gamma5[5]*sign_table_gamma5[3];
       }

       buffer_accum[0][spin_sink*4+spin_source][0]=buffer_CP[0][spin_sink*4+spin_source][0];
       buffer_accum[0][spin_sink*4+spin_source][1]=buffer_CP[0][spin_sink*4+spin_source][1];
     }
   }


   for (int time_extent = 1; time_extent < tp->T ; ++ time_extent ){

     mult_cp_matrix(buffer_CP[time_extent],buffer_write[tp->T - time_extent]);
     for (int spin_source=0; spin_source<tp->d; ++spin_source) {
       for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {

         buffer_CP[time_extent][spin_sink*4+spin_source][0]*=sign_table_source[5]*sign_table_sink[5]*sign_table_source[3]*sign_table_sink[3]*-1;
         buffer_CP[time_extent][spin_sink*4+spin_source][1]*=sign_table_source[5]*sign_table_sink[5]*sign_table_source[3]*sign_table_sink[3]*-1;

         if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
           buffer_CP[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[5]*sign_table_gamma5[3];
           buffer_CP[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[5]*sign_table_gamma5[3];

         }
         if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
           buffer_CP[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[5]*sign_table_gamma5[3];
           buffer_CP[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[5]*sign_table_gamma5[3];
         }

         buffer_accum[time_extent][spin_sink*4+spin_source][0]=buffer_CP[time_extent][spin_sink*4+spin_source][0];
         buffer_accum[time_extent][spin_sink*4+spin_source][1]=buffer_CP[time_extent][spin_sink*4+spin_source][1];

       }
     }
   }

   fini_3level_dtable(&buffer_CP);
   free(sign_table_source);
   free(sign_table_sink);
   free(sign_table_gamma5);


}
static inline void mult_with_cpt( double ***buffer_accum,  double *** buffer_write, twopoint_function_type * tp, char *gamma_string_source, char *gamma_string_sink ){


   /* CPT */
   /* tagname suffix  is CPT, abbreviation from charge conjugation + parity + time reversal */
   int *sign_table_source=init_1level_itable(9);
   int *sign_table_sink=init_1level_itable(9);

   int *sign_table_gamma5=init_1level_itable(9);

   apply_signs_discrete_symmetry ( sign_table_source, gamma_string_source );
   apply_signs_discrete_symmetry ( sign_table_sink,   gamma_string_sink  );

   apply_signs_discrete_symmetry ( sign_table_gamma5,   "g5"  );

   double ***buffer_CPT= init_3level_dtable(tp->T,tp->d*tp->d,2);

   for (int time_extent = 0; time_extent < tp->T ; ++ time_extent ){

     mult_cpt_matrix(buffer_CPT[time_extent],buffer_write[time_extent]);

     for (int spin_source=0; spin_source<tp->d; ++spin_source) {

       for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {

         buffer_CPT[time_extent][spin_sink*4+spin_source][0]*=sign_table_source[8]*sign_table_sink[8]*sign_table_sink[3]*sign_table_source[3]*(1.);
         buffer_CPT[time_extent][spin_sink*4+spin_source][1]*=sign_table_source[8]*sign_table_sink[8]*sign_table_sink[3]*sign_table_source[3]*(1.);

         if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
           buffer_CPT[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[8]*sign_table_gamma5[3];
           buffer_CPT[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[8]*sign_table_gamma5[3];
         }
         if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
           buffer_CPT[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[8]*sign_table_gamma5[3];
           buffer_CPT[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[8]*sign_table_gamma5[3];
         }

         buffer_accum[time_extent][spin_sink*4+spin_source][0]=buffer_CPT[time_extent][spin_sink*4+spin_source][0];
         buffer_accum[time_extent][spin_sink*4+spin_source][1]=buffer_CPT[time_extent][spin_sink*4+spin_source][1];

       }
     }
   }
   fini_3level_dtable(&buffer_CPT);
   free(sign_table_source);
   free(sign_table_sink);
   free(sign_table_gamma5);


}




static inline void mult_with_t( double ***buffer_accum,  double *** buffer_write, twopoint_function_type * tp, char *gamma_string_source, char *gamma_string_sink){
   int *sign_table_source=init_1level_itable(9);
   int *sign_table_sink=init_1level_itable(9);

   int *sign_table_gamma5=init_1level_itable(9);

   apply_signs_discrete_symmetry ( sign_table_source, gamma_string_source );
   apply_signs_discrete_symmetry ( sign_table_sink,   gamma_string_sink  );

   apply_signs_discrete_symmetry ( sign_table_gamma5,   "g5"  );

   /* time reversal */
   /* tagname suffix  is T, abbreviation from time reversal */
   double ***buffer_time_reversal= init_3level_dtable(tp->T,tp->d*tp->d,2);

   mult_time_reversal_matrix(buffer_time_reversal[0],buffer_write[0]);

   for (int spin_source=0; spin_source<tp->d; ++spin_source) {
     for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {

       buffer_time_reversal[0][spin_sink*4+spin_source][0]*=sign_table_source[4]*sign_table_sink[4]*-1;
       buffer_time_reversal[0][spin_sink*4+spin_source][1]*=sign_table_source[4]*sign_table_sink[4]*-1;

       buffer_accum[0][spin_sink*4+spin_source][0]=buffer_time_reversal[0][spin_sink*4+spin_source][0];
       buffer_accum[0][spin_sink*4+spin_source][1]=buffer_time_reversal[0][spin_sink*4+spin_source][1];
     }
   }

   for (int time_extent = 1; time_extent < tp->T ; ++ time_extent ){
     mult_time_reversal_matrix(buffer_time_reversal[time_extent],buffer_write[tp->T-time_extent]);
     for (int spin_source=0; spin_source<tp->d; ++spin_source) {
       for (int spin_sink=0; spin_sink < tp->d; ++spin_sink ) {

         buffer_time_reversal[time_extent][spin_sink*4+spin_source][0]*=sign_table_source[4]*sign_table_sink[4]*-1;
         buffer_time_reversal[time_extent][spin_sink*4+spin_source][1]*=sign_table_source[4]*sign_table_sink[4]*-1;

         if ((strcmp(gamma_string_source,"C")==0) || (strcmp(gamma_string_source,"Cg4")==0) || (strcmp(gamma_string_source,"cg1g4g5")==0) || (strcmp(gamma_string_source,"cg2g4g5")==0) || (strcmp(gamma_string_source,"cg3g4g5")==0) ){
            buffer_time_reversal[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[4];
            buffer_time_reversal[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[4];
         }
         if ((strcmp(gamma_string_sink,"C")==0) || (strcmp(gamma_string_sink,"Cg4")==0) || (strcmp(gamma_string_sink,"cg1g4g5")==0) || (strcmp(gamma_string_sink,"cg2g4g5")==0) || (strcmp(gamma_string_sink,"cg3g4g5")==0) ){
            buffer_time_reversal[time_extent][spin_sink*4+spin_source][0]*=sign_table_gamma5[4];
            buffer_time_reversal[time_extent][spin_sink*4+spin_source][1]*=sign_table_gamma5[4];

         }

         buffer_accum[time_extent][spin_sink*4+spin_source][0]=buffer_time_reversal[time_extent][spin_sink*4+spin_source][0];
         buffer_accum[time_extent][spin_sink*4+spin_source][1]=buffer_time_reversal[time_extent][spin_sink*4+spin_source][1];

       }
     }
   }

   fini_3level_dtable(&buffer_time_reversal);
   free(sign_table_source);
   free(sign_table_sink);
   free(sign_table_gamma5);


}


static inline void sink_and_source_gamma_list( char *filename, char *tagname, int number_of_gammas_source, int number_of_gammas_sink, char ***gamma_string_source, char ***gamma_string_sink, int pion_source, char *diagram_name){

#if 0
    hid_t file_id, dataset_id, attr_id; /* identifiers */
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_id = H5Dopen2 ( file_id, tagname, H5P_DEFAULT );
    attr_id = H5Aopen_idx (dataset_id, 0);
    hid_t atype = H5Aget_type(attr_id);
    hsize_t sz = H5Aget_storage_size(attr_id);
    char* date_source = (char *)malloc(sizeof(char)*(sz+1));
    H5Aread( attr_id,atype, (void*)date_source);
    printf("%s\n", date_source);
    H5Aclose(attr_id );
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    printf("data secr %s %s %s\n",date_source, filename,tagname );
#endif
    *gamma_string_source =(char **)malloc(sizeof(char *)*number_of_gammas_source);
    for (int j=0; j<number_of_gammas_source; j++){
      (*gamma_string_source)[j]=(char *)malloc(sizeof(char)*100);
    }

    *gamma_string_sink =(char **)malloc(sizeof(char *)*number_of_gammas_sink);
    for (int j=0; j<number_of_gammas_sink; j++){
      (*gamma_string_sink)[j]=(char *)malloc(sizeof(char)*100);
    }

    if (strcmp(diagram_name, "T")==0){
      snprintf((*gamma_string_source)[0],100,"%s","Cg5");
      snprintf((*gamma_string_source)[1],100,"%s","C");
      snprintf((*gamma_string_source)[2],100,"%s","Cg5g4");
      snprintf((*gamma_string_source)[3],100,"%s","Cg4");

      snprintf((*gamma_string_sink)[0],100,"%s","cg1");
      snprintf((*gamma_string_sink)[1],100,"%s","cg2");
      snprintf((*gamma_string_sink)[2],100,"%s","cg3");
      snprintf((*gamma_string_sink)[3],100,"%s","cg1g4");
      snprintf((*gamma_string_sink)[4],100,"%s","cg2g4");
      snprintf((*gamma_string_sink)[5],100,"%s","cg3g4");
      snprintf((*gamma_string_sink)[6],100,"%s","cg1g4g5");
      snprintf((*gamma_string_sink)[7],100,"%s","cg2g4g5");
      snprintf((*gamma_string_sink)[8],100,"%s","cg3g4g5");

    }
    else if (strcmp(diagram_name ,"piN")==0){
      snprintf((*gamma_string_source)[0],100,"%s","Cg5");
      snprintf((*gamma_string_source)[1],100,"%s","C");
      snprintf((*gamma_string_source)[2],100,"%s","Cg5g4");
      snprintf((*gamma_string_source)[3],100,"%s","Cg4");
   
      snprintf((*gamma_string_sink)[0],100,"%s","Cg5");
      snprintf((*gamma_string_sink)[1],100,"%s","C");
      snprintf((*gamma_string_sink)[2],100,"%s","Cg5g4");
      snprintf((*gamma_string_sink)[3],100,"%s","Cg4");
    }
    else if (strcmp(diagram_name ,"D")==0){
      snprintf((*gamma_string_source)[0],100,"%s","cg1");
      snprintf((*gamma_string_source)[1],100,"%s","cg2");
      snprintf((*gamma_string_source)[2],100,"%s","cg3");
      snprintf((*gamma_string_source)[3],100,"%s","cg1g4");
      snprintf((*gamma_string_source)[4],100,"%s","cg2g4");
      snprintf((*gamma_string_source)[5],100,"%s","cg3g4");
      snprintf((*gamma_string_source)[6],100,"%s","cg1g4g5");
      snprintf((*gamma_string_source)[7],100,"%s","cg2g4g5");
      snprintf((*gamma_string_source)[8],100,"%s","cg3g4g5");


      snprintf((*gamma_string_sink)[0],100,"%s","cg1");
      snprintf((*gamma_string_sink)[1],100,"%s","cg2");
      snprintf((*gamma_string_sink)[2],100,"%s","cg3");
      snprintf((*gamma_string_sink)[3],100,"%s","cg1g4");
      snprintf((*gamma_string_sink)[4],100,"%s","cg2g4");
      snprintf((*gamma_string_sink)[5],100,"%s","cg3g4");
      snprintf((*gamma_string_sink)[6],100,"%s","cg1g4g5");
      snprintf((*gamma_string_sink)[7],100,"%s","cg2g4g5");
      snprintf((*gamma_string_sink)[8],100,"%s","cg3g4g5");

    }
    else if (strcmp(diagram_name ,"N")==0){
      snprintf((*gamma_string_source)[0],100,"%s","Cg5");
      snprintf((*gamma_string_source)[1],100,"%s","C");
      snprintf((*gamma_string_source)[2],100,"%s","Cg5g4");
      snprintf((*gamma_string_source)[3],100,"%s","Cg4");

      snprintf((*gamma_string_sink)[0],100,"%s","Cg5");
      snprintf((*gamma_string_sink)[1],100,"%s","C");
      snprintf((*gamma_string_sink)[2],100,"%s","Cg5g4");
      snprintf((*gamma_string_sink)[3],100,"%s","Cg4");
    }

#if 0
    char *token;
    printf("Tag %s\n", date_source);
    token=strtok(date_source,",{");
    /* walk through other tokens */
    /* we do ignore the external gammas here because they are always one */
    token = strtok(NULL, ",{");
    token = strtok(NULL, ",{");
    for (int j=0; j<number_of_gammas_source;++j){
      token = strtok(NULL, ",{");
      if (token == NULL){
        fprintf(stderr,"# [sink_and_source_gamma_list] Error in obtaining the source gammas, check out your hdf5 file\n");
        exit(1);
      }
      if (j==(number_of_gammas_source-1)){
        snprintf((*gamma_string_source)[j],strlen(token),"%s",token);
      }
      else{
        snprintf((*gamma_string_source)[j],100,"%s",token);
      }
      fprintf(stdout,"# [sink_and_source_gamma_list] Gamma string list at the source %d %s\n", j, (*gamma_string_source)[j]);
    }
    if (pion_source ==1) {
      /* if the source contains the pion we have to jump over its gamma structure */
      token=strtok(NULL,",{");
    }
    for (int j=0; j<number_of_gammas_sink;++j){
      token=strtok(NULL, ",{");
      if (token == NULL){
        fprintf(stderr,"# [sink_and_source_gamma_list] Error in obtaining the source gammas, check out your hdf5 file\n");
        exit(1);
      }
      if (j==(number_of_gammas_sink-1)){
        char *temp=strtok(token, "}");
        snprintf((*gamma_string_sink)[j],strlen(temp)+1,"%s",temp);
      }
      else{
        snprintf((*gamma_string_sink)[j],100,"%s",token);
      }
      fprintf(stdout,"# [sink_and_source_gamma_list] Gamma string list at the sink %d %s\n", j, (*gamma_string_sink)[j]);
    }
    free(date_source);

#endif
}

int i2index( int pi2x, int pi2y, int pi2z ){
   int ret_index=-1;
   for (int i=0; i<g_seq_source_momentum_number; ++i){
     if ((g_seq_source_momentum_list[i][0]==pi2x) && (g_seq_source_momentum_list[i][1]==pi2y) && (g_seq_source_momentum_list[i][2]==pi2z)){
       ret_index=i;
       break;
     }
   }
   return ret_index;
}

int f1f2index_piN( int ** buffer_mom, int pf1x, int pf1y, int pf1z, int pf2x, int pf2y, int pf2z ){
   int ret_index=-1;
   for (int i=0; i<343; ++i){
     if ((buffer_mom[i][3]==pf1x) && (buffer_mom[i][4]==pf1y) && (buffer_mom[i][5]==pf1z) && (buffer_mom[i][6]==pf2x) && (buffer_mom[i][7]==pf2y) && (buffer_mom[i][8]==pf2z)){
       ret_index=i; 
       break; 
     }
   }
   return ret_index;
}


int f1f2index_T( int ** buffer_mom, int pf1x, int pf1y, int pf1z ){
   int ret_index=-1;
   for (int i=0; i<27; ++i){
     if ((buffer_mom[i][3]==pf1x) && (buffer_mom[i][4]==pf1y) && (buffer_mom[i][5]==pf1z) ){
       ret_index=i;
       break;
     }
   }
   return ret_index;
}


int f1f2index_D( int ** buffer_mom, int pf1x, int pf1y, int pf1z ){
   int ret_index=-1;
   for (int i=0; i<27; ++i){
     if ((buffer_mom[i][0]==pf1x) && (buffer_mom[i][1]==pf1y) && (buffer_mom[i][2]==pf1z) ){
       ret_index=i;
       break;
     }
   }
   return ret_index;
}




int gammaindex( int sink , int source, char *diagram_name ) {

   int sink_index;
   int source_index;
   if (strcmp(diagram_name, "T")==0){
   if ( sink == 14 ){
       sink_index= 0;
   }
   else if (sink == 11) {
       sink_index= 1;
   }
   else if (sink == 8){
       sink_index= 2;
   }
   else if (sink == 2){
       sink_index=3; 
   }
   else if (sink == 9){
       sink_index=0;
   }
   else if (sink == 0){
       sink_index=1;
   }
   else if (sink == 7){
       sink_index=2;
   }
   else if (sink ==13){
       sink_index=3;
   }
   else if (sink ==4){
       sink_index=4;
   }
   else if (sink == 15) {
       sink_index=5;
   }
   else if (sink == 12){
       sink_index=6;
   }
   else if (sink == 5){
       sink_index=7;
   }     
   else if (sink == 10){
       sink_index=8;
   }
   else { 
    fprintf(stderr, "Unregocnized sink gamma %d\n",sink);
    exit(1);
   }
   if ( source == 14 ){ 
       source_index= 0; 
   }
   else if (source == 11) { 
       source_index= 1; 
   } 
   else if (source == 8){ 
       source_index= 2; 
   } 
   else if (source == 2){ 
       source_index=3;  
   } 
   else {  
    fprintf(stderr, "Unregocnized source gamma %d\n",source); 
    exit(1); 
   } 
   return ( sink_index* 4 + source_index);
   }
   else if ((strcmp(diagram_name, "piN")==0) && (strcmp(diagram_name, "N")==0)) {
   if ( source == 14 ){
       source_index= 0;
   }
   else if (source == 11) {
       source_index= 1;
   }
   else if (source == 8){
       source_index= 2;
   }
   else if (source == 2){
       source_index=3;
   }
   else {
    fprintf(stderr, "Unregocnized source gamma %d\n",source);
    exit(1);
   }
   if ( sink == 14 ){
       sink_index= 0;
   }
   else if (sink == 11) {
       sink_index= 1;
   }
   else if (sink == 8){
       sink_index= 2;
   }
   else if (sink == 2){
       sink_index=3;
   }
   else {
    fprintf(stderr, "Unregocnized sink gamma %d\n",sink);
    exit(1);
   }
   return ( sink_index* 4 + source_index);
   }
   else if (strcmp(diagram_name, "D")==0){
   if ( sink == 14 ){
       sink_index= 0;
   }
   else if (sink == 11) {
       sink_index= 1;
   }
   else if (sink == 8){
       sink_index= 2;
   }
   else if (sink == 2){
       sink_index=3;
   }
   else if (sink == 9){
       sink_index=0;
   }
   else if (sink == 0){
       sink_index=1;
   }
   else if (sink == 7){
       sink_index=2;
   }
   else if (sink ==13){
       sink_index=3;
   }
   else if (sink ==4){
       sink_index=4;
   }
   else if (sink == 15) {
       sink_index=5;
   }
   else if (sink == 12){
       sink_index=6;
   }
   else if (sink == 5){
       sink_index=7;
   }
   else if (sink == 10){
       sink_index=8;
   }
   else {
    fprintf(stderr, "Unregocnized sink gamma %s %d\n",diagram_name, sink);
    exit(1);
   }
   if ( source == 14 ){
       source_index= 0;
   }
   else if (source == 11) {
       source_index= 1;
   }
   else if (source == 8){
       source_index= 2;
   }
   else if (source == 2){
       source_index=3;
   }
   else if (source == 9){
       source_index=0;
   }
   else if (source == 0){
       source_index=1;
   }
   else if (source == 7){
       source_index=2;
   }
   else if (source ==13){
       source_index=3;
   }
   else if (source ==4){
       source_index=4;
   }
   else if (source == 15) {
       source_index=5;
   }
   else if (source == 12){
       source_index=6;
   }
   else if (source == 5){
       source_index=7;
   }
   else if (source == 10){
       source_index=8;
   }
   else {
    fprintf(stderr, "Unregocnized source gamma %s %d\n", diagram_name, source);
    exit(1);
   }
   return ( sink_index* 9 + source_index);
   }
   else {
     fprintf(stdout,"Diagram is not implemented\n");
     exit(1);
   }

}

int udli_id_lookup_piNpiN( char udli_list[][500], int udli_count, int *momtable_sink_f2, char *gamma_string_sink_member, int *momtable_sink_f1, int *momtable_source_i2, char *gamma_string_source_member){
    char udliname[500];
    int udliid=-1;
    snprintf(udliname, 500, "%d%d%d%s%d%d%d%d%d%d%s", momtable_sink_f2[0],
                                                      momtable_sink_f2[1],
                                                      momtable_sink_f2[2],
                                                      gamma_string_sink_member,
                                                      momtable_sink_f1[0],
	                                              momtable_sink_f1[1],
                                                      momtable_sink_f1[2],
                                                      momtable_source_i2[0],
                                                      momtable_source_i2[1],
                                                      momtable_source_i2[2],
                                                      gamma_string_source_member);
    for (int i=0; i<udli_count; ++i){
      if (strcmp(udli_list[i], udliname)==0){
        udliid=i;
        break;
      }
    }
    return udliid;
}
void udli_id_store_piNpiN( char udli_list[][500], int udli_count, int *momtable_sink_f2, char *gamma_string_sink_member, int *momtable_sink_f1, int *momtable_source_i2, char *gamma_string_source_member){
    snprintf(udli_list[udli_count], 500, "%d%d%d%s%d%d%d%d%d%d%s", momtable_sink_f2[0],
                                                      momtable_sink_f2[1],
                                                      momtable_sink_f2[2],
                                                      gamma_string_sink_member,
                                                      momtable_sink_f1[0],
                                                      momtable_sink_f1[1],
                                                      momtable_sink_f1[2],
                                                      momtable_source_i2[0],
                                                      momtable_source_i2[1],
                                                      momtable_source_i2[2],
                                                      gamma_string_source_member);
    
}
int udli_id_lookup_BB( char udli_list[][500], int udli_count, int *momtable_sink, char *gamma_string_sink_member, char *gamma_string_source_member){
    char udliname[500];
    int udliid=-1;
    snprintf(udliname, 500, "%d%d%d%s%s", momtable_sink[0],
                                          momtable_sink[1],
                                          momtable_sink[2],
                                          gamma_string_sink_member,
                                          gamma_string_source_member);
    for (int i=0; i<udli_count; ++i){
      if (strcmp(udli_list[i], udliname)==0){
        udliid=i;
        break;
      }
    }
    return udliid;
} 

void udli_id_store_BB( char udli_list[][500], int udli_count, int *momtable_sink, char *gamma_string_sink_member, char *gamma_string_source_member){
    snprintf(udli_list[udli_count], 500, "%d%d%d%s%s", momtable_sink[0],
                                          momtable_sink[1],
                                          momtable_sink[2],
                                          gamma_string_sink_member,
                                          gamma_string_source_member);
}


int udli_id_lookup_piND( char udli_list[][500], int udli_count, int *momtable_sink, char* gamma_string_sink_member,int *momtable_source_1, int *momtable_source_2, char *gamma_string_source_member){
    char udliname[500];
    int udliid=-1;
    snprintf(udliname, 500, "%d%d%d%s%d%d%d%d%d%d%s", momtable_sink[0],
                                          momtable_sink[1],
                                          momtable_sink[2],
                                          gamma_string_sink_member,
                                          momtable_source_1[0],
                                          momtable_source_1[1],
                                          momtable_source_1[2],
                                          momtable_source_2[0],
                                          momtable_source_2[1],
                                          momtable_source_2[2],
                                          gamma_string_source_member);
    for (int i=0; i<udli_count; ++i){
      if (strcmp(udli_list[i], udliname)==0){
        udliid=i;
        break;
      }
    }
    return udliid;
}

void udli_id_store_piND( char udli_list[][500], int udli_count, int *momtable_sink, char* gamma_string_sink_member,int *momtable_source_1, int *momtable_source_2, char *gamma_string_source_member){
    snprintf(udli_list[udli_count], 500, "%d%d%d%s%d%d%d%d%d%d%s", momtable_sink[0],
                                          momtable_sink[1],
                                          momtable_sink[2],
                                          gamma_string_sink_member,
                                          momtable_source_1[0],
                                          momtable_source_1[1],
                                          momtable_source_1[2],
                                          momtable_source_2[0],
                                          momtable_source_2[1],
                                          momtable_source_2[2],
                                          gamma_string_source_member);
}
int udli_id_lookup_DpiN( char udli_list[][500], int udli_count, int *momtable_sink_1, char *gamma_string_sink_member, int *momtable_sink_2, int *momtable_source_1, char *gamma_string_source_member){
    char udliname[500];
    int udliid=-1;
    snprintf(udliname, 500, "%d%d%d%s%d%d%d%d%d%d%s", momtable_sink_1[0],
                                          momtable_sink_1[1],
                                          momtable_sink_1[2],
                                          gamma_string_sink_member,
                                          momtable_sink_2[0],
                                          momtable_sink_2[1],
                                          momtable_sink_2[2],
                                          momtable_source_1[0],
                                          momtable_source_1[1],
                                          momtable_source_1[2],
                                          gamma_string_source_member);
    for (int i=0; i<udli_count; ++i){
      if (strcmp(udli_list[i], udliname)==0){
        udliid=i;
        break;
      }
    }
    return udliid;

}
void udli_id_store_DpiN( char udli_list[][500], int udli_count, int *momtable_sink_1, char *gamma_string_sink_member, int *momtable_sink_2, int *momtable_source_1, char *gamma_string_source_member){
  snprintf(udli_list[udli_count], 500, "%d%d%d%s%d%d%d%d%d%d%s", momtable_sink_1[0],
                                          momtable_sink_1[1],
                                          momtable_sink_1[2],
                                          gamma_string_sink_member,
                                          momtable_sink_2[0],
                                          momtable_sink_2[1],
                                          momtable_sink_2[2],
                                          momtable_source_1[0],
                                          momtable_source_1[1],
                                          momtable_source_1[2],
                                          gamma_string_source_member);
}


char *convert_gamma_to_string ( int gamma ) {

   char *ret_string=(char *)malloc(sizeof(char)*100);
   if ( gamma == 9){
     snprintf(ret_string,100,"cg1,1");
   }
   else if (gamma== 0){
     snprintf(ret_string,100,"cg2,1");
   }
   else if (gamma==7){
     snprintf(ret_string,100,"cg3,1");
   }
   else if (gamma==13){
     snprintf(ret_string,100,"cg1g4,1");
   }
   else if (gamma==4){
     snprintf(ret_string,100,"cg2g4,1");
   }
   else if (gamma==15){
     snprintf(ret_string,100,"cg3g4,1");
   }
   else if ( gamma==12){
     snprintf(ret_string,100,"cg1g4g5,5");
   }
   else if ( gamma==5) {
     snprintf(ret_string,100,"cg2g4g5,5");
   }
   else if (gamma==10) {
     snprintf(ret_string,100,"cg3g4g5,5");
   }
   else if (gamma==14) {
     snprintf(ret_string,100, "Cg5,1");
   }
   else if (gamma==11) {
     snprintf(ret_string,100, "C,5");
   }
   else if (gamma==8) {
     snprintf(ret_string,100, "Cg5g4,1");
   }
   else if (gamma==2) {
     snprintf(ret_string,100, "Cg4,5");
   }
   else {
     fprintf(stderr, "# [convert_gamma_to_string] Non recognized gamma in conversion\n");}
   return ret_string;

}

char *tagname_forgamma_multiplets ( int *gamma_multiplet, int length ){

   char *ret_string=(char *)malloc(sizeof(char)*100);
   if (length==3){
     if ((gamma_multiplet[0]==9) && (gamma_multiplet[1]==0)  && (gamma_multiplet[2]==7)){
       snprintf(ret_string, 100, "cgxyz");
     }
     else if ((gamma_multiplet[0]==13) && (gamma_multiplet[1]==4)  && (gamma_multiplet[2]==15)){
       snprintf(ret_string, 100, "cgxyzg0");
     }
     else if ((gamma_multiplet[0]==12) && (gamma_multiplet[1]==5)  && (gamma_multiplet[2]==10)){
       snprintf(ret_string, 100, "cgxyzg0g5");

     }
     else {
       fprintf(stderr, "# [tagname_forgamma_multiplets] No recognized multiplett %d %d %d\n", gamma_multiplet[0],gamma_multiplet[1],gamma_multiplet[2]);
       exit(1);
     }
   }
   else if (length==1) {
     if (gamma_multiplet[0]==14){
       snprintf(ret_string,100, "cg5");
     }
     else if  (gamma_multiplet[0]==11){
       snprintf(ret_string,100, "c");
     }
     else if  (gamma_multiplet[0]==8){
       snprintf(ret_string,100, "cg5g4");
     }
     else if  (gamma_multiplet[0]==2){
       snprintf(ret_string,100, "cg4");
     }
     else {
       fprintf(stderr, "# [tagname_forgamma_multiplets] No recognized multiplett length 1 the gamma structure is the following %d \n",gamma_multiplet[0]);
       exit(1);
     }
  }
  else {
       fprintf(stderr, "# [tagname_forgamma_multiplets] No recognized length of mupliplett other length \n");
       exit(1);
  }
  return ret_string;

}

/***********************************************************
 * main program
 * This code is supposed to be produce table of projected correlation
 * functions for piN-D system
 ***********************************************************/
int main(int argc, char **argv) {
 
#define _ZCOEFF_EPS 8.e-12


  int c;
  int filename_set = 0;
  int exitstatus;
  int check_reference_rotation = 0;
  char filename[400];
  char tagname[400];

  
  /***********************************************************
   * initialize MPI if used
   ***********************************************************/
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  /***********************************************************
   * evaluate command line arguments
   ***********************************************************/
  while ((c = getopt(argc, argv, "ch?f:")) != -1) {
    switch (c) {
    case 'f':
      strcpy(filename, optarg);
      filename_set=1;
      break;
    case 'c':
      check_reference_rotation = 1;
      fprintf ( stdout, "# [piN2piN_projection_apply_plegma] check_reference_rotation set to %d\n", check_reference_rotation );
      break;
    case 'h':
    case '?':
    default:
      exit(1);
      break;
    }
  }

  /***********************************************************
   * set the default values
   ***********************************************************/
  if(filename_set==0) strcpy(filename, "cvc.input");
  read_input_parser(filename);

#ifdef HAVE_OPENMP
  omp_set_num_threads(g_num_threads);
#else
  fprintf(stdout, "# [piN2piN_projection_apply_plegma] Warning, resetting global thread number to 1\n");
  g_num_threads = 1;
#endif

  /***********************************************************
   * package-own initialization of MPI parameters
   ***********************************************************/
  mpi_init(argc, argv);

  /***********************************************************
   * report git version
   ***********************************************************/
  if ( g_cart_id == 0 ) {
    fprintf(stdout, "# [piN2piN_projection_apply_plegma] git version = %s\n", g_gitversion);
  }

  /***********************************************************
   * set geometry
   ***********************************************************/
  exitstatus = init_geometry();
  if( exitstatus != 0 ) {
    fprintf(stderr, "# [piN2piN_projection_apply_plegma] Error from init_geometry, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
    EXIT(1);
  }
  geometry();

  /***********************************************************
   * set io process
   ***********************************************************/
  int const io_proc = get_io_proc ();

  /***********************************************************
   * TEST: report size of twopoint_function_type
   ***********************************************************/
  if ( io_proc == 0 ) {
    fprintf ( stdout, "# [piN2piN_projection_apply_plegma] sizeof twopoint_function_type = %lu\n", sizeof ( twopoint_function_type ) );
  }

  int const momentum_orbit_pref[4][3] = { {0,0,0}, {0,0,1}, {1,1,0}, {1,1,1} };

  int udli_count = 0;
  int udliid;
  char udli_list[MAX_UDLI_NUM][500];
  char udli_name[500];
  double ***correlation_function[MAX_UDLI_NUM];
//  double ***correlation_function_symmetrized[MAX_UDLI_NUM];


  /******************************************************
   * check source coords list
   ******************************************************/
  for ( int i = 0; i < g_source_location_number; i++ ) {
    g_source_coords_list[i][0] = ( g_source_coords_list[i][0] +  T_global ) %  T_global;
    g_source_coords_list[i][1] = ( g_source_coords_list[i][1] + LX_global ) % LX_global;
    g_source_coords_list[i][2] = ( g_source_coords_list[i][2] + LY_global ) % LY_global;
    g_source_coords_list[i][3] = ( g_source_coords_list[i][3] + LZ_global ) % LZ_global;
  }

  /******************************************************
   * total number of source locations, 
   * base x coherent
   *
   * fill a list of all source coordinates
   ******************************************************/
  int const source_location_number = g_source_location_number * g_coherent_source_number;
  int ** source_coords_list = init_2level_itable ( source_location_number, 4 );
  if( source_coords_list == NULL ) {
    fprintf ( stderr, "# [piN2piN_projection_apply_plegma] Error from init_2level_itable %s %d\n", __FILE__, __LINE__ );
    EXIT( 43 );
  }
  for ( int ib = 0; ib < g_source_location_number; ib++ ) {
    g_source_coords_list[ib][0] = ( g_source_coords_list[ib][0] +  T_global ) %  T_global;
    g_source_coords_list[ib][1] = ( g_source_coords_list[ib][1] + LX_global ) % LX_global;
    g_source_coords_list[ib][2] = ( g_source_coords_list[ib][2] + LY_global ) % LY_global;
    g_source_coords_list[ib][3] = ( g_source_coords_list[ib][3] + LZ_global ) % LZ_global;

    int const t_base = g_source_coords_list[ib][0];

    for ( int ic = 0; ic < g_coherent_source_number; ic++ ) {
      int const ibc = ib * g_coherent_source_number + ic;

      int const t_coherent = ( t_base + ( T_global / g_coherent_source_number ) * ic ) % T_global;
      source_coords_list[ibc][0] = t_coherent;
      source_coords_list[ibc][1] = ( g_source_coords_list[ib][1] + (LX_global/2) * ic ) % LX_global;
      source_coords_list[ibc][2] = ( g_source_coords_list[ib][2] + (LY_global/2) * ic ) % LY_global;
      source_coords_list[ibc][3] = ( g_source_coords_list[ib][3] + (LZ_global/2) * ic ) % LZ_global;
    }
  }


  /***********************************************************
   * initialize gamma matrix algebra and several
   * gamma basis matrices
   ***********************************************************/
  init_gamma_matrix ("plegma");

  /******************************************************
   * set gamma matrices
   *   tmLQCD counting
   ******************************************************/
  gamma_matrix_type gamma[16];
  for ( int i = 0; i < 16; i++ ) {
    gamma_matrix_set ( &(gamma[i]), i, 1. );
  }


  /*****************************************************
   * 
   * Loop over source positions
   *
   *****************************************************/
 
  for ( int k = 0; k < source_location_number; k++ ) {

    int const pi2[3] = {
            g_seq_source_momentum_list[0][0],
            g_seq_source_momentum_list[0][1],
            g_seq_source_momentum_list[0][2] };
#ifdef HAVE_HDF5
    /***********************************************************
     * read data block from h5 file
     ***********************************************************/

    if ( (strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"piN")==0)) {
      snprintf ( filename, 400, "plegma_format_piN_%04d.h5",
                         Nconf);
    }
    else if ((strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {
      snprintf ( filename, 400, "mod_Diagramm%04d_sx%02dsy%02dsz%02dst%03d_T.h5", Nconf,
                                                                                source_coords_list[k][1],
                                                                                source_coords_list[k][2],
                                                                                source_coords_list[k][3],
                                                                                source_coords_list[k][0]);
    }
    else if ((strcmp(g_twopoint_function_list[0].particlename_source,"D")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {
      snprintf ( filename, 400, "mod_Diagramm%04d_sx%02dsy%02dsz%02dst%03d_D.h5", Nconf,
                                                                                source_coords_list[k][1],
                                                                                source_coords_list[k][2],
                                                                                source_coords_list[k][3],
                                                                                source_coords_list[k][0]);
    }
    else if ((strcmp(g_twopoint_function_list[0].particlename_source,"N")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"N")==0)) {
      snprintf ( filename, 400, "mod_Diagramm%04d_sx%02dsy%02dsz%02dst%03d_N.h5", Nconf,
                                                                                source_coords_list[k][1],
                                                                                source_coords_list[k][2],
                                                                                source_coords_list[k][3],
                                                                                source_coords_list[k][0]);
    }

    int ** buffer_mom ;
    if ( (strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"piN")==0)) {

      buffer_mom=init_2level_itable ( 343, 9 );

    }
    else if ( (strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {

      buffer_mom=init_2level_itable ( 27, 6 );

    }
    else if ( (strcmp(g_twopoint_function_list[0].particlename_source,"D")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {

      buffer_mom=init_2level_itable ( 27, 3 );

    }
    else if ( (strcmp(g_twopoint_function_list[0].particlename_source,"N")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"N")==0)) {

      buffer_mom=init_2level_itable ( 27, 3 );

    }

    if ( buffer_mom == NULL ) {
       fprintf(stderr, "# [piN2piN_projection_apply_plegma]  Error from ,init_2level_itable %s %d\n", __FILE__, __LINE__ );
       EXIT(12);
    }

    if (((strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"piN")==0)) ||
        ((strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0))) {
       snprintf(tagname, 400, "sx%.02dsy%.02dsz%.02dst%.02d/pi2=%d_%d_%d/mvec",
                         source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0],
                         pi2[0],pi2[1],pi2[2]);
    }
    else if (((strcmp(g_twopoint_function_list[0].particlename_source,"N")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"N")==0)) ||
             ((strcmp(g_twopoint_function_list[0].particlename_source,"D")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0))) {
       snprintf(tagname, 400, "sx%.02dsy%.02dsz%.02dst%.02d/mvec",
                         source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0]);
    }

    fprintf(stdout, "# [piN2piN_projection_apply_plegma] Try to read in momenta %s %s\n",filename, tagname);
    exitstatus = read_from_h5_file ( (void*)(buffer_mom[0]), filename, tagname, io_proc, 1 );
    if ( exitstatus != 0 ) {
      fprintf(stderr, "[piN2piN_projection_apply_plegma] Error from read_from_h5_file, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
      EXIT(12);
    }
#endif


    char **gamma_string_list_source;
    char **gamma_string_list_sink;

/*
    snprintf ( filename, 400, "plegma_format_piN_%04d.h5",
                         Nconf
                          );
*/

    if ((strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"piN")==0)) {
       sink_and_source_gamma_list( filename, tagname, 4, 4, &gamma_string_list_source, &gamma_string_list_sink, 1,"piN");

    }
    else if ((strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {
       sink_and_source_gamma_list( filename, tagname, 4, 9, &gamma_string_list_source, &gamma_string_list_sink, 1,"T");
    }
    else if ((strcmp(g_twopoint_function_list[0].particlename_source,"D")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {
       sink_and_source_gamma_list( filename, tagname, 9, 9, &gamma_string_list_source, &gamma_string_list_sink, 1,"D");

    } 
    else if ((strcmp(g_twopoint_function_list[0].particlename_source,"N")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"N")==0)) {
       sink_and_source_gamma_list( filename, tagname, 4, 4, &gamma_string_list_source, &gamma_string_list_sink, 1,"N");
    }



    double ******buffer_sum;

    if ((strcmp(g_twopoint_function_list[0].particlename_source,"D")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {
      buffer_sum  = init_6level_dtable(1,  g_twopoint_function_list[0].T, 27, 81, 16, 2  );
    }
    else if ((strcmp(g_twopoint_function_list[0].particlename_source,"N")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"N")==0)) {
      buffer_sum  = init_6level_dtable(1,  g_twopoint_function_list[0].T, 27, 16, 16, 2  );
    }
    else if ((strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {
      buffer_sum  = init_6level_dtable(27,  g_twopoint_function_list[0].T, 27, 36, 16, 2  );
    }
    else if ((strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"piN")==0)) {
      buffer_sum  = init_6level_dtable(27,  g_twopoint_function_list[0].T, 343, 16, 16, 2  );
    }






    for ( int ipi2 = 0; ipi2 < g_seq_source_momentum_number; ipi2++ ) {


      if ((strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"piN")==0)) {
       snprintf ( tagname, 400, "/%04d/sx%.02dsy%.02dsz%.02dst%.02d/pi2=%d_%d_%d/mxb-mxb",
                         Nconf,
                         source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0],
                         g_seq_source_momentum_list[ipi2][0],
                         g_seq_source_momentum_list[ipi2][1],
                         g_seq_source_momentum_list[ipi2][2]);
      }
      else if ((strcmp(g_twopoint_function_list[0].particlename_source,"piN")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {
       snprintf ( tagname, 400, "/sx%.02dsy%.02dsz%.02dst%.02d/pi2=%d_%d_%d/T",
                         source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0],
                         g_seq_source_momentum_list[ipi2][0],
                         g_seq_source_momentum_list[ipi2][1],
                         g_seq_source_momentum_list[ipi2][2]);
      }
      else if ((strcmp(g_twopoint_function_list[0].particlename_source,"D")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"D")==0)) {
       snprintf ( tagname, 400, "/sx%.02dsy%.02dsz%.02dst%.02d/D",
                         source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0]);
      }  
      else if ((strcmp(g_twopoint_function_list[0].particlename_source,"N")==0) && (strcmp(g_twopoint_function_list[0].particlename_sink,"N")==0)) {
       snprintf ( tagname, 400, "/sx%.02dsy%.02dsz%.02dst%.02d/N",
                         source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0]);
      }


      exitstatus = read_from_h5_file ( (void*)(buffer_sum[ipi2][0][0][0][0]), filename, tagname, io_proc );

      if ( exitstatus != 0 ) {
         fprintf(stderr, "# [piN2piN_projection_apply_plegma] Error from read_from_h5_file, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
         EXIT(12);
      }

    }/* loop over pi2 */


    fprintf(stdout, "# [piN2piN_projection_apply_plegma] reading done\n");
   
  /******************************************************
   * loop on 2-point functions
   ******************************************************/

    for ( int i2pt = 0; i2pt < g_twopoint_function_number; i2pt++ ) {

      printf("# [piN2piN_projection_apply] start analyzing twopoint function index %d\n", i2pt);

      for (int ii=0; ii< g_twopoint_function_list[i2pt].nlistmomentumf1; ++ii){

        printf("# [piN2piN_projection_apply] pf1 (%d %d %d)\n", g_twopoint_function_list[i2pt].pf1list[ii][0] , g_twopoint_function_list[i2pt].pf1list[ii][1] , g_twopoint_function_list[i2pt].pf1list[ii][2] );
        printf("# [piN2piN_projection_apply] pf2 (%d %d %d)\n", g_twopoint_function_list[i2pt].pf2list[ii][0] , g_twopoint_function_list[i2pt].pf2list[ii][1] , g_twopoint_function_list[i2pt].pf2list[ii][2] );

      }


      const int numberofnplets_source = g_twopoint_function_list[i2pt].number_of_gammas_i1/g_twopoint_function_list[i2pt].contniuum_spin_particle_source;

      const int spin1dimension_source = g_twopoint_function_list[i2pt].contniuum_spin_particle_source;

      const int numberofnplets_sink = g_twopoint_function_list[i2pt].number_of_gammas_f1/g_twopoint_function_list[i2pt].contniuum_spin_particle_sink;

      const int spin1dimension_sink   = g_twopoint_function_list[i2pt].contniuum_spin_particle_sink;

      const int spin1212dimension = g_twopoint_function_list[i2pt].d;

      hid_t file_id, group_id, dataset_id, dataspace_id;  /* identifiers */
      herr_t      status;


      snprintf ( filename, 400, "Projected_corrfunc_sink_%s_source_%s_group_%s_irrep_%s_%.04d.h5",
                         g_twopoint_function_list[i2pt].particlename_sink,
                         g_twopoint_function_list[i2pt].particlename_source,
                         g_twopoint_function_list[i2pt].group,
                         g_twopoint_function_list[i2pt].irrep,
                         Nconf);

      struct stat fileStat;
      if(stat( filename, &fileStat) < 0 ) {

        /* Open an existing file. */
        fprintf ( stdout, "# [test_hdf5] create new file %s\n",filename );
        file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      } else {
        fprintf ( stdout, "# [test_hdf5] open existing file %s\n", filename );
        file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
      }


      /* Determining the dimension of the irrep */

      int const Ptot[3] = {
       g_twopoint_function_list[i2pt].pf1list[0][0] + g_twopoint_function_list[i2pt].pf2list[0][0] , 
       g_twopoint_function_list[i2pt].pf1list[0][1] + g_twopoint_function_list[i2pt].pf2list[0][1] , 
       g_twopoint_function_list[i2pt].pf1list[0][2] + g_twopoint_function_list[i2pt].pf2list[0][2]  };

      snprintf ( filename, 400, "projection_coefficients_%s_group_%s_irrep_%s.h5",
                            g_twopoint_function_list[i2pt].particlename_source,
                            g_twopoint_function_list[i2pt].group,
                            g_twopoint_function_list[i2pt].irrep );

      fprintf(stdout, "# [piN2piN_projection_apply]  determining dimension of the irrep\n");

      if ( strcmp(g_twopoint_function_list[i2pt].particlename_source, "piN") == 0) {
        snprintf( tagname, 400, "/pfx%dpfy%dpfz%d_pi%dN%d/mu_%d/beta_%d/Nreplicas_Nps_Ndimirrep",  Ptot[0],Ptot[1], Ptot[2], g_twopoint_function_list[i2pt].total_momentum_pion_source[0], g_twopoint_function_list[i2pt].total_momentum_nucleon_source[0], 0, 0 );}
      else { 
        snprintf( tagname, 400, "/pfx%dpfy%dpfz%d/mu_%d/beta_%d/Nreplicas_Nps_Ndimirrep",  Ptot[0],Ptot[1], Ptot[2],  0, 0 );       }

      int *tmp= init_1level_itable( 3 );

      exitstatus = read_from_h5_file ( (void*)(tmp), filename, tagname, io_proc, 1 );
      if (exitstatus != 0){
        fprintf(stderr, "Error in opening file for determining the dimension of irrep\n");
        exit(1);
      }

      const int irrep_dim=tmp[2];

      fini_1level_itable(&tmp);
      fprintf(stdout,"# [piN2piN_projection_apply] Irrep dimension=%d\n", irrep_dim);


      snprintf ( filename, 400, "projection_coefficients_%s_group_%s_irrep_%s.h5",
        g_twopoint_function_list[i2pt].particlename_source,
        g_twopoint_function_list[i2pt].group,
        g_twopoint_function_list[i2pt].irrep );


      int ****Nps_Nreplica_source= init_4level_itable( 1, irrep_dim, g_twopoint_function_list[i2pt].ncombination_total_momentum_source, 3 );

      int ****Nps_Nreplica_sink= init_4level_itable( 1, irrep_dim, g_twopoint_function_list[i2pt].ncombination_total_momentum_sink, 3 );
      int *****momtable_source=(int *****)malloc(sizeof(int****)*1);
      momtable_source[0]=(int ****)malloc(sizeof(int***)*irrep_dim);
      for (int i=0; i<irrep_dim; ++i){
        momtable_source[0][i]=(int***)malloc(sizeof(int**)*g_twopoint_function_list[i2pt].ncombination_total_momentum_source);
      }

      int *****momtable_sink=(int *****)malloc(sizeof(int****)*1);
      momtable_sink[0]=(int ****)malloc(sizeof(int***)*irrep_dim);
      for (int i=0; i<irrep_dim; ++i){
        momtable_sink[0][i]=(int***)malloc(sizeof(int**)*g_twopoint_function_list[i2pt].ncombination_total_momentum_sink);
      }
      double ******projection_coeff_c_ORT=(double ******)malloc(sizeof(double*****)*1);
      projection_coeff_c_ORT[0]=(double *****)malloc(sizeof(double****)*irrep_dim);
      for (int i=0; i<irrep_dim; ++i){
        projection_coeff_c_ORT[0][i]=(double****)malloc(sizeof(double***)*g_twopoint_function_list[i2pt].ncombination_total_momentum_source);
      }

      double ******projection_coeff_a_ORT=(double ******)malloc(sizeof(double*****)*1);
      projection_coeff_a_ORT[0]=(double *****)malloc(sizeof(double****)*irrep_dim);
      for (int i=0; i<irrep_dim; ++i){
        projection_coeff_a_ORT[0][i]=(double****)malloc(sizeof(double***)*g_twopoint_function_list[i2pt].ncombination_total_momentum_sink);
      }

      for ( int ibeta = 0; ibeta < 1; ibeta++ ) {
 
        for ( int imu = 0; imu < irrep_dim; imu++ ) {

          for ( int icombination_source=0; icombination_source < g_twopoint_function_list[i2pt].ncombination_total_momentum_source ; ++icombination_source ) {

            snprintf ( filename, 400, "projection_coefficients_%s_group_%s_irrep_%s.h5",
              g_twopoint_function_list[i2pt].particlename_source,
              g_twopoint_function_list[i2pt].group,
              g_twopoint_function_list[i2pt].irrep );


            if (strcmp( g_twopoint_function_list[i2pt].particlename_source, "piN") == 0){
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d_pi%dN%d/mu_%d/beta_%d/Nreplicas_Nps_Ndimirrep",  Ptot[0],Ptot[1], Ptot[2], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], imu, ibeta );
            }
            else {
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d/mu_%d/beta_%d/Nreplicas_Nps_Ndimirrep",  Ptot[0],Ptot[1], Ptot[2], imu, ibeta );
            }

            exitstatus = read_from_h5_file ( (void*)(Nps_Nreplica_source[ibeta][imu][icombination_source]), filename, tagname, io_proc, 1 );
            if (exitstatus != 0){
              fprintf(stderr, "Error in opening file %s tag %s\n", filename, tagname );
              exit(1);
            }
            printf("Nps_Nreplica_source[ibeta][imu][icombination_source][1]%d\n",  Nps_Nreplica_source[ibeta][imu][icombination_source][1]);
            momtable_source[ibeta][imu][icombination_source] = init_2level_itable ( Nps_Nreplica_source[ibeta][imu][icombination_source][1], 6 );
            if (strcmp ( g_twopoint_function_list[i2pt].particlename_source, "piN")==0){
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d_pi%dN%d/mu_%d/beta_%d/momlist_f1f2",  Ptot[0],Ptot[1], Ptot[2], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], imu, ibeta );
            }
            else {
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d/mu_%d/beta_%d/momlist_f1f2",  Ptot[0],Ptot[1], Ptot[2], imu, ibeta );
            }


            exitstatus = read_from_h5_file ( (void*)(momtable_source[ibeta][imu][icombination_source][0]), filename, tagname, io_proc, 1 );
            if (exitstatus != 0){
               fprintf(stderr, "# [piN2piN_projection_apply] Error in opening %s %s\n",filename, tagname );
               exit(1); 
            }

            fprintf(stdout, "# [piN2piN_projection_apply] momtable source loaded %s %s\n", tagname, filename );

            projection_coeff_c_ORT[ibeta][imu][icombination_source]= init_3level_dtable( Nps_Nreplica_source[ibeta][imu][icombination_source][1]*spin1212dimension*spin1dimension_source, Nps_Nreplica_source[ibeta][imu][icombination_source][0], 2);
            if (strcmp ( g_twopoint_function_list[i2pt].particlename_source, "piN")==0){
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d_pi%dN%d/mu_%d/beta_%d/c_data_ort",  Ptot[0],Ptot[1], Ptot[2], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], imu, ibeta );
            } else {
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d/mu_%d/beta_%d/c_data_ort",  Ptot[0],Ptot[1], Ptot[2], imu, ibeta );
            }
            fprintf(stdout, "# [piN2piN_projection_apply] projection coeff source loaded %s %s\n", filename, tagname);

            exitstatus = read_from_h5_file ( (void*)(projection_coeff_c_ORT[ibeta][imu][icombination_source][0][0]), filename, tagname, io_proc, 0 );
            if (exitstatus != 0){
              fprintf(stderr, "# [piN2piN_projection_apply] Error in opening %s %s\n",filename, tagname );
               exit(1);
            }
          }
          for ( int icombination_sink=0; icombination_sink < g_twopoint_function_list[i2pt].ncombination_total_momentum_sink ; ++icombination_sink ) {

            snprintf ( filename, 400, "projection_coefficients_%s_group_%s_irrep_%s.h5",
              g_twopoint_function_list[i2pt].particlename_sink,
              g_twopoint_function_list[i2pt].group,
              g_twopoint_function_list[i2pt].irrep );


            if (strcmp( g_twopoint_function_list[i2pt].particlename_sink, "piN") == 0){
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d_pi%dN%d/mu_%d/beta_%d/Nreplicas_Nps_Ndimirrep",  Ptot[0],Ptot[1], Ptot[2], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], imu, ibeta );
            }
            else {
             snprintf( tagname, 400, "/pfx%dpfy%dpfz%d/mu_%d/beta_%d/Nreplicas_Nps_Ndimirrep",  Ptot[0],Ptot[1], Ptot[2], imu, ibeta );
            }

            exitstatus = read_from_h5_file ( (void*)(Nps_Nreplica_sink[ibeta][imu][icombination_sink]), filename, tagname, io_proc, 1 );
            if (exitstatus != 0){
              fprintf(stderr, "Error in opening file %s tag %s\n", filename, tagname );
              exit(1);
            }

            momtable_sink[ibeta][imu][icombination_sink] = init_2level_itable ( Nps_Nreplica_sink[ibeta][imu][icombination_sink][1], 6 );
            if (strcmp ( g_twopoint_function_list[i2pt].particlename_sink, "piN")==0){
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d_pi%dN%d/mu_%d/beta_%d/momlist_f1f2",  Ptot[0],Ptot[1], Ptot[2], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], imu, ibeta );
            }
            else {
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d/mu_%d/beta_%d/momlist_f1f2",  Ptot[0],Ptot[1], Ptot[2], imu, ibeta );
            }


            exitstatus = read_from_h5_file ( (void*)(momtable_sink[ibeta][imu][icombination_sink][0]), filename, tagname, io_proc, 1 );
            if (exitstatus != 0){
               fprintf(stderr, "# [piN2piN_projection_apply] Error in opening %s %s\n",filename, tagname );
               exit(1);
            }

            fprintf(stdout, "# [piN2piN_projection_apply] momtable sink loaded %s %s\n", tagname, filename );


            projection_coeff_a_ORT[ibeta][imu][icombination_sink]= init_3level_dtable(Nps_Nreplica_sink[ibeta][imu][icombination_sink][0], Nps_Nreplica_sink[ibeta][imu][icombination_sink][1]*spin1212dimension*spin1dimension_sink,2);

            if (strcmp(g_twopoint_function_list[i2pt].particlename_sink , "piN" )==0) {
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d_pi%dN%d/mu_%d/beta_%d/a_data_ort",  Ptot[0],Ptot[1], Ptot[2], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], imu, ibeta );
            } else {
              snprintf( tagname, 400, "/pfx%dpfy%dpfz%d/mu_%d/beta_%d/a_data_ort",  Ptot[0],Ptot[1], Ptot[2],  imu, ibeta );

            }
            printf("Filename %s Tagname %s\n", filename, tagname );
            exitstatus = read_from_h5_file ( (void*)(projection_coeff_a_ORT[ibeta][imu][icombination_sink][0][0]), filename, tagname, io_proc, 0 );
            if (exitstatus != 0){
              fprintf(stderr, "# [piN2piN_projection_apply] Error in opening %s %s\n",filename, tagname );
              exit(1);
            }
            fprintf(stdout, "# [piN2piN_projection_apply] coeff sink loaded %s %s\n", tagname, filename );

           }
         
         }
    
       }


      /******************************************************
       * Open the table for creation and annihilation interpolating operator coefficients 
       *****************************************************/
       snprintf(tagname,400,"/sx%0.2dsy%0.2dsz%0.2dst%03d", source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0]);
       status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

       status = H5Gget_objinfo (file_id, tagname, 0, NULL);
   
       if (status != 0){

         /* Create a group named "/MyGroup" in the file. */
         group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

         /* Close the group. */
         status = H5Gclose(group_id);

       }

       snprintf(tagname,400,"/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d", source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2]);
       status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

       status = H5Gget_objinfo (file_id, tagname, 0, NULL);

       if (status < 0){

         /* Create a group named "/MyGroup" in the file. */
         group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

         /* Close the group. */
         status = H5Gclose(group_id);

       }

       for ( int ibeta = 0; ibeta < 1; ibeta++ ) {


         snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d", source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta );

         status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

         status = H5Gget_objinfo (file_id, tagname, 0, NULL);
   
         if (status != 0){

           /* Create a group named "/MyGroup" in the file. */
           group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

           /* Close the group. */
           status = H5Gclose(group_id);

         } 

         for ( int imu = 0; imu < irrep_dim; imu++ ) {

           snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d/mu_%d",source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta, imu );

           status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

           status = H5Gget_objinfo (file_id, tagname, 0, NULL);
           if (status != 0){

             /* Create a group named "/MyGroup" in the file. */
             group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

             /* Close the group. */
             status = H5Gclose(group_id);

           }
 
           for ( int icombination_source=0; icombination_source < g_twopoint_function_list[i2pt].ncombination_total_momentum_source ; ++icombination_source ) {
             snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d/mu_%d/source_totmomN%dp%d/", source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta, imu, g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source] );

             status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

             status = H5Gget_objinfo (file_id, tagname, 0, NULL);
             if (status != 0){

               /* Create a group named "/MyGroup" in the file. */
               group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

               /* Close the group. */
               status = H5Gclose(group_id);

             }    
             for ( int icombination_sink=0; icombination_sink < g_twopoint_function_list[i2pt].ncombination_total_momentum_sink ; ++icombination_sink ) {

               snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d/mu_%d/source_totmomN%dp%d/sink_totmomN%dp%d", source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta, imu, g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source],g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink]  );

               status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

               status = H5Gget_objinfo (file_id, tagname, 0, NULL);
               if (status != 0){

                 /* Create a group named "/MyGroup" in the file. */
                 group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                 /* Close the group. */
                 status = H5Gclose(group_id);

               }


               for (int gamma_nplettid_source=0; gamma_nplettid_source < numberofnplets_source; ++gamma_nplettid_source) {

                 int *gamma_table_source=init_1level_itable(spin1dimension_source);
                 for (int iii=0; iii<spin1dimension_source;++iii){
                   gamma_table_source[iii]=g_twopoint_function_list[i2pt].list_of_gammas_i1[gamma_nplettid_source*spin1dimension_source+iii][0];
                 }
                 char *gamma_string_source=tagname_forgamma_multiplets ( gamma_table_source, spin1dimension_source );

                 snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d/mu_%d/source_totmomN%dp%d/sink_totmomN%dp%d/source_%s/",source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta, imu, g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source],g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink], gamma_string_source);

                 status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

                 status = H5Gget_objinfo (file_id, tagname, 0, NULL);
                 if (status != 0){

                   /* Create a group named "/MyGroup" in the file. */
                   group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                   /* Close the group. */
                   status = H5Gclose(group_id);

                 }

                 for (int gamma_nplettid_sink=0; gamma_nplettid_sink < numberofnplets_sink; ++gamma_nplettid_sink) {

                   int *gamma_table_sink=init_1level_itable(spin1dimension_sink);
                   for (int iii=0; iii<spin1dimension_sink;++iii){
                     gamma_table_sink[iii]=g_twopoint_function_list[i2pt].list_of_gammas_f1[gamma_nplettid_sink*spin1dimension_sink+iii][0];
                   }
                   char *gamma_string_sink=tagname_forgamma_multiplets ( gamma_table_sink, spin1dimension_sink );

                   snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d/mu_%d/source_totmomN%dp%d/sink_totmomN%dp%d/source_%s/sink_%s",source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta, imu, g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source],g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink], gamma_string_source, gamma_string_sink);

                   status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

                   status = H5Gget_objinfo (file_id, tagname, 0, NULL);
                   if (status != 0){

                     /* Create a group named "/MyGroup" in the file. */
                     group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                     /* Close the group. */
                     status = H5Gclose(group_id);

                   }
                   int npsso_length=Nps_Nreplica_source[ibeta][imu][icombination_source][0]; 
                   for (int nreplicum_source=0; nreplicum_source < npsso_length ; ++nreplicum_source ) {

                     snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d/mu_%d/source_totmomN%dp%d/sink_totmomN%dp%d/source_%s/sink_%s/Replicasource_%d", source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta, imu, g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source],g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink], gamma_string_source, gamma_string_sink, nreplicum_source);

                     status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

                     status = H5Gget_objinfo (file_id, tagname, 0, NULL);
                     if (status != 0){

                       /* Create a group named "/MyGroup" in the file. */
                       group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                       /* Close the group. */
                       status = H5Gclose(group_id);

                     }

                     int npssi_length=Nps_Nreplica_sink[ibeta][imu][icombination_sink][0];
 
                     for ( int nreplicum_sink=0; nreplicum_sink < npssi_length ; ++nreplicum_sink ) {

                       snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d/mu_%d/source_totmomN%dp%d/sink_totmomN%dp%d/source_%s/sink_%s/Replicasource_%d/Replicasink_%d", source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta, imu, g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source],g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink], gamma_string_source, gamma_string_sink, nreplicum_source, nreplicum_sink);

                       status = H5Eset_auto(NULL, H5P_DEFAULT, NULL);

                       status = H5Gget_objinfo (file_id, tagname, 0, NULL);
                       if (status != 0){
 
                         /* Create a group named "/MyGroup" in the file. */
                         group_id = H5Gcreate2(file_id, tagname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                         /* Close the group. */
                         status = H5Gclose(group_id);

                       }
 
                       double **projected_correlation_function=init_2level_dtable ( g_twopoint_function_list[i2pt].T, 2);

                       double **projected_correlation_function_symmetrized=init_2level_dtable ( g_twopoint_function_list[i2pt].T, 2);

                       /* loop for possible individial momenta leading to the same tot momentum  at the source */

                       for ( int i1_source =0; i1_source < Nps_Nreplica_source[ibeta][imu][icombination_source][1] ; ++i1_source ) {
      
                         /* loop for possible gamma structures for the delta the triplet for example */
                          
                         for ( int i2_source=0; i2_source < spin1dimension_source ; ++i2_source ) {

                           char *gamma_string_source_member=convert_gamma_to_string( gamma_table_source[i2_source]);

                           /* loop for possible individial momenta leading to the same tot momentum  at the sink */

                           for ( int i1_sink =0; i1_sink < Nps_Nreplica_sink[ibeta][imu][icombination_sink][1] ; ++i1_sink ) {

                             /* loop for possible gamma structures for the delta the triplet for example */

                             for ( int i2_sink =0; i2_sink < spin1dimension_sink ; ++i2_sink ) {

                               char *gamma_string_sink_member=convert_gamma_to_string( gamma_table_sink[i2_sink]);

                               char *correlation_function_filename_suffix=(char *)malloc(sizeof(char)*100);
                               char *correlation_function_tagname_suffix =(char *)malloc(sizeof(char)*100);

                               if ( (strcmp(g_twopoint_function_list[i2pt].particlename_sink, "D") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "D") == 0) ){

                                snprintf(correlation_function_filename_suffix, 100, "D");
                                snprintf(correlation_function_tagname_suffix, 100, "D"); 

                               } else if ( (strcmp(g_twopoint_function_list[i2pt].particlename_sink, "D") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "piN") == 0 )){

                                snprintf(correlation_function_filename_suffix, 100, "T");
                                snprintf(correlation_function_tagname_suffix, 100, "T"); 

                               } else if ( (strcmp(g_twopoint_function_list[i2pt].particlename_sink, "piN") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "D") == 0)){

                                snprintf(correlation_function_filename_suffix, 100, "TpiNsink");
                                snprintf(correlation_function_tagname_suffix, 100, "TpiNsink");

                               } else if ( (strcmp(g_twopoint_function_list[i2pt].particlename_sink, "piN") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "piN") == 0)){

                                snprintf(correlation_function_filename_suffix, 100, "piN");
                                snprintf(correlation_function_tagname_suffix, 100, "piN");

                               } else if ( (strcmp(g_twopoint_function_list[i2pt].particlename_sink, "N") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "N") == 0)){

                                snprintf(correlation_function_filename_suffix, 100, "N");
                                snprintf(correlation_function_tagname_suffix, 100, "N");

                               }
                               else { 
                                fprintf(stderr, "# [piN2piN_projection_apply_plegma] only piN -delta system is implemented\n");
                                exit(1);
                               }

 
                               if ((strcmp(g_twopoint_function_list[i2pt].particlename_sink, "piN") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "piN") == 0)){


                                 int pf1x=momtable_sink[ibeta][imu][icombination_sink][i1_sink][0];
                                 int pf1y=momtable_sink[ibeta][imu][icombination_sink][i1_sink][1];
                                 int pf1z=momtable_sink[ibeta][imu][icombination_sink][i1_sink][2];
                                 int pf2x=momtable_sink[ibeta][imu][icombination_sink][i1_sink][3];
                                 int pf2y=momtable_sink[ibeta][imu][icombination_sink][i1_sink][4];
                                 int pf2z=momtable_sink[ibeta][imu][icombination_sink][i1_sink][5];

                                 int pi1x=momtable_source[ibeta][imu][icombination_source][i1_source][0];
                                 int pi1y=momtable_source[ibeta][imu][icombination_source][i1_source][1];
                                 int pi1z=momtable_source[ibeta][imu][icombination_source][i1_source][2];
                                 int pi2x=momtable_source[ibeta][imu][icombination_source][i1_source][3];
                                 int pi2y=momtable_source[ibeta][imu][icombination_source][i1_source][4];
                                 int pi2z=momtable_source[ibeta][imu][icombination_source][i1_source][5];

                                 
                                 int applied_transformations=0;
 
                                 udliid=udli_id_lookup_piNpiN( udli_list, udli_count, &momtable_sink[ibeta][imu][icombination_sink][i1_sink][3], gamma_string_sink_member,&momtable_sink[ibeta][imu][icombination_sink][i1_sink][0],&momtable_source[ibeta][imu][icombination_source][i1_source][3],gamma_string_source_member);
                                 if (udliid==-1){
                                   snprintf ( tagname, 400, "/sx%.02dsy%.02dsz%.02dst%03d/gf25/pf2x%.02dpf2y%.02dpf2z%.02d/gf1%s/pf1x%.02dpf1y%.02dpf1z%.02d/gi25/pi2x%.02dpi2y%.02dpi2z%.02d/gi1%s/pi1x%.02dpi1y%.02dpi1z%.02d/%s", source_coords_list[k][1],
                                                        source_coords_list[k][2],
                                                        source_coords_list[k][3],
                                                        source_coords_list[k][0],
                                                        pf2x,
                                                        pf2y,
                                                        pf2z,
                                                        gamma_string_sink_member,
                                                        pf1x,
                                                        pf1y,
                                                        pf1z,
                                                        pi2x,
                                                        pi2y,
                                                        pi2z,
                                                        gamma_string_source_member,
                                                        pi1x,
                                                        pi1y,
                                                        pi1z,
                                                        correlation_function_tagname_suffix);

                                   const int i_total_momentum=Ptot[0]*Ptot[0]+Ptot[1]*Ptot[1]+Ptot[2]*Ptot[2];
                             
                                   snprintf ( filename, 400, "%s%04d_PX%.02dPY%.02dPZ%.02d_%s.h5",
                                                        filename_prefix,
                                                        Nconf,
                                                        momentum_orbit_pref[i_total_momentum][0],
                                                        momentum_orbit_pref[i_total_momentum][1],
                                                        momentum_orbit_pref[i_total_momentum][2],
                                                        correlation_function_filename_suffix);
                                   const int spin1212_dimension=g_twopoint_function_list[i2pt].d*g_twopoint_function_list[i2pt].d ;
                                   const int time_extent=g_twopoint_function_list[i2pt].T;

                                   correlation_function[udli_count]=init_3level_dtable( time_extent, spin1212_dimension ,2 );
//                                   correlation_function_symmetrized[udli_count]=init_3level_dtable( time_extent, spin1212_dimension ,2 );
                                     
                                   int i2idx=i2index( pi2x, pi2y, pi2z);
                                   int f1f2idx=f1f2index_piN(buffer_mom, pf1x, pf1y, pf1z, pf2x, pf2y, pf2z);
                                   int gammaidx=gammaindex(gamma_table_sink[i2_sink],gamma_table_source[i2_source],"piN");
                                   
                                   
                                   for (int time_c=0; time_c<time_extent; ++time_c){
                                     for (int spin_index=0; spin_index < spin1212_dimension; ++ spin_index){
                                       for (int realimag=0; realimag < 2; ++realimag ){
                                         correlation_function[udli_count][time_c][spin_index][realimag]=buffer_sum[i2idx][time_c][f1f2idx][gammaidx][spin_index][realimag];
                                       }
                                     }
                                     if ((gamma_table_source[i2_source]==11) || (gamma_table_source[i2_source]==2)){
                                       mult_with_gamma5_matrix_adj_source(correlation_function[udli_count][time_c]);
                                     }
                                     if ((gamma_table_sink[i2_sink]==11) || (gamma_table_sink[i2_sink]==2)){
                                       mult_with_gamma5_matrix_sink(correlation_function[udli_count][time_c]);
                                     }
/*                                     for (int spin_index=0; spin_index < spin1212_dimension; ++ spin_index){
                                       for (int realimag=0; realimag < 2; ++realimag ){
                                         correlation_function_symmetrized[udli_count][time_c][spin_index][realimag]=correlation_function[udli_count][time_c][spin_index][realimag];

                                       }
                                     }*/
                                   }

                                   /* applying time reversal */
#if 0
                                   double ***buffer_t_reversal= init_3level_dtable(time_extent, spin1212_dimension ,2);
                                   mult_with_t( buffer_t_reversal,  correlation_function[udli_count], &g_twopoint_function_list[i2pt], gamma_string_source_member, gamma_string_sink_member );
                                   
                                   for (int time_coord=0; time_coord < time_extent; ++time_coord ){
                                     for (int spin_inner=0; spin_inner < 16; ++spin_inner) {
                                       for (int realimag=0; realimag < 2; ++realimag){
                                         correlation_function_symmetrized[udli_count][time_coord][spin_inner][realimag]+=buffer_t_reversal[time_coord][spin_inner][realimag];
                                       }
                                     }
                                   }
                                   fini_3level_dtable( &buffer_t_reversal );

                                   applied_transformations+=2;

                                   /* charge conjugated, the sink and the source gamma is interchanged, the momenta left unchanged*/
                                   /* The corresponding discrete symmetry transformations                           */
                                   /*                                                   (1) charge conjugation (C)  */
                                   /*                                                   (2) charge conjugation + time-reversal (CT) */
     
                                   const int pi1osq=pi1x*pi1x+pi1y*pi1y+pi1z*pi1z;
                                   const int pi2osq=pi2x*pi2x+pi2y*pi1y+pi2z*pi2z;
                                   const int pf1osq=pf1x*pf1x+pi1y*pi1y+pi1z*pi1z;
                                   const int pf2osq=pf2x*pf2x+pi1y*pi1y+pi1z*pi1z;


                                   if ( (pi1osq<=3) && (pi2osq<=3) && (pf1osq<=3) && (pf2osq<=3)) {

                                     double ***buffer_write= init_3level_dtable(time_extent, spin1212_dimension,2);

                                     gammaidx=gammaindex(gamma_table_source[i2_source],gamma_table_sink[i2_sink]);
                                     i2idx=i2index( pf2x, pf2y, pf2z);
                                     f1f2idx=f1f2index( buffer_mom, pi1x, pi1y, pi1z, pi2x, pi2y, pi2z );

                                     for (int time = 0; time < time_extent ; ++ time ){
                                       for (int spin_inner=0; spin_inner < spin1212_dimension; ++spin_inner) {
                                         for (int realimag=0; realimag < 2; ++realimag){
                                           buffer_write[time][spin_inner][realimag]=buffer_sum[i2idx][time][f1f2idx][gammaidx][spin_inner][realimag];
                                         }
                                       }
                                       if ((gamma_table_source[i2_source]==11) || (gamma_table_source[i2_source]==2)){
                                         mult_with_gamma5_matrix_sink(buffer_write[time_extent]);

                                       }
                                       if ((gamma_table_sink[i2_sink]==11) || (gamma_table_source[i2_sink]==2)){
                                         mult_with_gamma5_matrix_adj_source(buffer_write[time_extent]);

                                       }

                                     }
                                     double ***buffer_c= init_3level_dtable(time_extent, spin1212_dimension,2);
                
                                     mult_with_c( buffer_c,  buffer_write, &g_twopoint_function_list[i2pt], gamma_string_source_member, gamma_string_sink_member);

                                     for (int time_coord=0; time_coord < time_extent; ++time_coord ){
                                       for (int spin_inner=0; spin_inner < spin1212_dimension; ++spin_inner){
                                         for (int realimag=0; realimag < 2; ++realimag){
                                           correlation_function_symmetrized[udli_count][time_coord][spin_inner][realimag]+=buffer_c[time_coord][spin_inner][realimag];
                                         }
                                       }
                                     }
                                     fini_3level_dtable(&buffer_c);

                                     double ***buffer_ct= init_3level_dtable(time_extent, spin1212_dimension, 2);
                                     mult_with_ct( buffer_ct,  buffer_write, &g_twopoint_function_list[i2pt], gamma_string_source_member, gamma_string_sink_member );

                                     for (int time_coord=0; time_coord < time_extent; ++time_coord ){
                                       for (int spin_inner=0; spin_inner < spin1212_dimension ; ++spin_inner) {
                                         for (int realimag=0; realimag < 2; ++realimag){
                                           correlation_function_symmetrized[udli_count][time_coord][spin_inner][realimag]+=buffer_ct[time_coord][spin_inner][realimag];
                                         }
                                       }
                                     }
                                     fini_3level_dtable( &buffer_ct);

                                     applied_transformations+=2;

                                     /* parity */

                                     i2idx=i2index( -pi2x, -pi2y, -pi2z);
                                     f1f2idx=f1f2index( buffer_mom, -pf1x, -pf1y, -pf1z, -pf2x, -pf2y, -pf2z );
                                     gammaidx=gammaindex(gamma_table_sink[i2_sink],gamma_table_source[i2_source]);


                                     for (int time = 0; time < time_extent ; ++ time ){
                                       for (int spin_inner=0; spin_inner < spin1212_dimension; ++spin_inner) {
                                         for (int realimag=0; realimag < 2; ++realimag){
                                           buffer_write[time][spin_inner][realimag]=buffer_sum[i2idx][time][f1f2idx][gammaidx][spin_inner][realimag];
                                         }
                                       }

                                       if ((gamma_table_source[i2_source]==11) || (gamma_table_source[i2_source]==2)){
                                         mult_with_gamma5_matrix_adj_source(buffer_write[time]);
                                       }
                                       if ((gamma_table_sink[i2_sink]==11) || (gamma_table_sink[i2_sink]==2)){
                                         mult_with_gamma5_matrix_sink(buffer_write[time]);
                                       }

                                     }
                  

                                     double ***buffer_p= init_3level_dtable(time_extent,spin1212_dimension, 2);
                                     mult_with_p( buffer_p,  buffer_write, &g_twopoint_function_list[i2pt], gamma_string_source_member, gamma_string_sink_member  );
                                     for (int time_coord=0; time_coord < time_extent; ++time_coord ){
                                       for (int spin_inner=0; spin_inner < spin1212_dimension; ++spin_inner) {
                                         for (int realimag=0; realimag < 2; ++realimag){
                                           correlation_function_symmetrized[udli_count][time_coord][spin_inner][realimag]+=buffer_p[time_coord][spin_inner][realimag];
                                         }
                                       }
                                     }

                                     fini_3level_dtable( &buffer_p );

                                     double ***buffer_pt= init_3level_dtable(time_extent, spin1212_dimension, 2);
                                     mult_with_pt( buffer_pt,  buffer_write, &g_twopoint_function_list[i2pt], gamma_string_source_member, gamma_string_sink_member  );
                                     for (int time_coord=0; time_coord < time_extent; ++time_coord ){
                                       for (int spin_inner=0; spin_inner < spin1212_dimension; ++spin_inner) {
                                         for (int realimag=0; realimag < 2; ++realimag){
                                           correlation_function_symmetrized[udli_count][time_coord][spin_inner][realimag]+=buffer_pt[time_coord][spin_inner][realimag];
                                         }
                                       }
                                     }
                                     fini_3level_dtable( &buffer_pt );

                                     /* Charge+parity conjugation */
                                     i2idx=i2index( -pf2x, -pf2y, -pf2z);
                                     f1f2idx=f1f2index( buffer_mom, -pi1x, -pi1y, -pi1z, -pi2x, -pi2y, -pi2z );
                                     gammaidx=gammaindex(gamma_table_source[i2_source],gamma_table_sink[i2_sink]);

                                     for (int time  = 0; time  < time_extent ; ++ time ){
                                       for (int spin_inner=0; spin_inner < spin1212_dimension; ++spin_inner) {
                                         for (int realimag=0; realimag < 2; ++realimag){
                                           buffer_write[time][spin_inner][realimag]=buffer_sum[i2idx][time][f1f2idx][gammaidx][spin_inner][realimag];
                                         }
                                       }
                                       if ((gamma_table_source[i2_source]==11) || (gamma_table_source[i2_source]==2) ){
                                         mult_with_gamma5_matrix_sink(buffer_write[time]);
                                       }

                                       if ((gamma_table_sink[i2_sink]==11) || (gamma_table_sink[i2_sink]==2)){
                                         mult_with_gamma5_matrix_adj_source(buffer_write[time]);
                                       }
 
                                     }

                                     double ***buffer_cp= init_3level_dtable(time_extent,spin1212_dimension,2);
                                     mult_with_cp( buffer_cp,  buffer_write, &g_twopoint_function_list[i2pt],  gamma_string_source_member, gamma_string_sink_member   );
                   
                                     for (int time=0; time  < time_extent; ++time  ){
                                       for (int spin_inner=0; spin_inner < spin1212_dimension; ++spin_inner) {
                                         for (int realimag=0; realimag < 2; ++realimag){
                                           correlation_function_symmetrized[udli_count][time][spin_inner][realimag]+=buffer_cp[time][spin_inner][realimag];
                                         }
                                       }
                                     }
                                     fini_3level_dtable( &buffer_cp );

                                     double ***buffer_cpt= init_3level_dtable(time_extent, spin1212_dimension, 2);
                                     mult_with_cpt( buffer_cpt,  buffer_write, &g_twopoint_function_list[i2pt], gamma_string_source_member, gamma_string_sink_member  );
                  
                                     for (int time=0; time < time_extent; ++time ){
                                       for (int spin_inner=0; spin_inner < spin1212_dimension; ++spin_inner) {
                                         for (int realimag=0; realimag < 2; ++realimag){
                                           correlation_function_symmetrized[udli_count][time][spin_inner][realimag]+=buffer_cpt[time][spin_inner][realimag];
                                         }
                                       }
                                     }

                                     applied_transformations+=2;
                                     fini_3level_dtable( &buffer_cpt );
                                     fini_3level_dtable( &buffer_write );

                                   }


                                   /* Perform normalization */
                                   for (int time_coord=0; time_coord < time_extent; ++time_coord ){
                                     for (int spin_inner=0; spin_inner < spin1212_dimension; ++spin_inner) {
                                       for (int realimag=0; realimag < 2; ++realimag){
                                         correlation_function_symmetrized[udli_count][time_coord][spin_inner][realimag]/=applied_transformations;
                                       }
                                     }
                                   }
#endif
                                   udli_id_store_piNpiN( udli_list, udli_count, &momtable_sink[ibeta][imu][icombination_sink][i1_sink][3], gamma_string_sink_member,&momtable_sink[ibeta][imu][icombination_sink][i1_sink][0],&momtable_source[ibeta][imu][icombination_source][i1_source][3],gamma_string_source_member);
                                   udliid=udli_count;
                                   ++udli_count;
                                   if (udli_count==MAX_UDLI_NUM){
                                     fprintf(stderr,"# [piN2piN_projection_apply] Error in the size of udli %d list\n", udli_count );  
                                     exit(1);
                                   }
                                 }
                                 /*else {
                                   fprintf(stdout,"# [piN2piN_projection_apply] udliid found %d udlicount %d\n", udliid, udli_count);
                                 }*/
                               } else if (((strcmp(g_twopoint_function_list[i2pt].particlename_sink, "N") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "N") == 0)) || (((strcmp(g_twopoint_function_list[i2pt].particlename_sink, "D") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "D") == 0)))){
                                 udliid=udli_id_lookup_BB( udli_list, udli_count, &momtable_sink[ibeta][imu][icombination_sink][i1_sink][0], gamma_string_sink_member,gamma_string_source_member);
                                 if (udliid==-1){
                                   int pf1x=momtable_sink[ibeta][imu][icombination_sink][i1_sink][0];
                                   int pf1y=momtable_sink[ibeta][imu][icombination_sink][i1_sink][1];
                                   int pf1z=momtable_sink[ibeta][imu][icombination_sink][i1_sink][2];
                                   int pi1x=momtable_source[ibeta][imu][icombination_source][i1_source][0];
                                   int pi1y=momtable_source[ibeta][imu][icombination_source][i1_source][1];
                                   int pi1z=momtable_source[ibeta][imu][icombination_source][i1_source][2];

                                   const int i_total_momentum=Ptot[0]*Ptot[0]+Ptot[1]*Ptot[1]+Ptot[2]*Ptot[2];

                                   correlation_function[udli_count]=init_3level_dtable( g_twopoint_function_list[i2pt].T,g_twopoint_function_list[i2pt].d*g_twopoint_function_list[i2pt].d ,2 );
				   const int spin1212_dimension=g_twopoint_function_list[i2pt].d*g_twopoint_function_list[i2pt].d ;
                                   const int time_extent=g_twopoint_function_list[i2pt].T;

                                   
                                   int f1f2idx=f1f2index_D( buffer_mom, pf1x, pf1y, pf1z );
                                   int gammaidx=gammaindex(gamma_table_sink[i2_sink],gamma_table_source[i2_source],"D");

                                   for (int time_c=0; time_c<time_extent; ++time_c){
                                     for (int spin_index=0; spin_index < spin1212_dimension; ++ spin_index){
                                       for (int realimag=0; realimag < 2; ++realimag ){
                                         correlation_function[udli_count][time_c][spin_index][realimag]=buffer_sum[0][time_c][f1f2idx][gammaidx][spin_index][realimag];
                                       }
                                     }
                                     if ((gamma_table_source[i1_source]==12) || (gamma_table_source[i1_source]==5) || (gamma_table_source[i1_source]==10) ){
                                       mult_with_gamma5_matrix_adj_source(correlation_function[udli_count][time_c]);
                                     }
                                     if ((gamma_table_sink[i1_sink]==12) || (gamma_table_sink[i1_sink]==5) || (gamma_table_sink[i1_sink]==10)) {
                                       mult_with_gamma5_matrix_sink(correlation_function[udli_count][time_c]);
                                     }
/*                                     for (int spin_index=0; spin_index < spin1212_dimension; ++ spin_index){
                                       for (int realimag=0; realimag < 2; ++realimag ){
                                         correlation_function_symmetrized[udli_count][time_c][spin_index][realimag]=correlation_function[udli_count][time_c][spin_index][realimag];

                                       }
                                     }*/
                                   }


                                   udli_id_store_BB( udli_list, udli_count, &momtable_sink[ibeta][imu][icombination_sink][i1_sink][0], gamma_string_sink_member,gamma_string_source_member);
                                   udliid=udli_count;
                                   ++udli_count;
                                   if (udli_count==MAX_UDLI_NUM){
                                     fprintf(stderr,"# [piN2piN_projection_apply] Error in the size of udli %d list\n", udli_count );
                                     exit(1);
                                   }

                                 }
 
                               } else  if ( (strcmp(g_twopoint_function_list[i2pt].particlename_sink, "D") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "piN") == 0 )){
                                 udliid=udli_id_lookup_piND( udli_list, udli_count, &momtable_sink[ibeta][imu][icombination_sink][i1_sink][0], gamma_string_sink_member,&momtable_source[ibeta][imu][icombination_source][i1_source][0],&momtable_source[ibeta][imu][icombination_source][i1_source][3],gamma_string_source_member);

                                 if (udliid==-1){
                                   int pf1x=momtable_sink[ibeta][imu][icombination_sink][i1_sink][0];
                                   int pf1y=momtable_sink[ibeta][imu][icombination_sink][i1_sink][1];
                                   int pf1z=momtable_sink[ibeta][imu][icombination_sink][i1_sink][2];
  
                                   int pi1x=momtable_source[ibeta][imu][icombination_source][i1_source][0];
                                   int pi1y=momtable_source[ibeta][imu][icombination_source][i1_source][1];
                                   int pi1z=momtable_source[ibeta][imu][icombination_source][i1_source][2];
                                   int pi2x=momtable_source[ibeta][imu][icombination_source][i1_source][3];
                                   int pi2y=momtable_source[ibeta][imu][icombination_source][i1_source][4];
                                   int pi2z=momtable_source[ibeta][imu][icombination_source][i1_source][5];


                                   const int i_total_momentum=Ptot[0]*Ptot[0]+Ptot[1]*Ptot[1]+Ptot[2]*Ptot[2];


                                   correlation_function[udli_count]=init_3level_dtable( g_twopoint_function_list[i2pt].T,g_twopoint_function_list[i2pt].d*g_twopoint_function_list[i2pt].d ,2 );

                                   const int spin1212_dimension=g_twopoint_function_list[i2pt].d*g_twopoint_function_list[i2pt].d ;
                                   const int time_extent=g_twopoint_function_list[i2pt].T;

                                   int i2idx=i2index( pi2x, pi2y, pi2z);
                                   int f1f2idx=f1f2index_T( buffer_mom, pf1x, pf1y, pf1z );
                                   int gammaidx=gammaindex(gamma_table_sink[i2_sink],gamma_table_source[i2_source],"T");

                                   for (int time_c=0; time_c<time_extent; ++time_c){
                                     for (int spin_index=0; spin_index < spin1212_dimension; ++ spin_index){
                                       for (int realimag=0; realimag < 2; ++realimag ){
                                         correlation_function[udli_count][time_c][spin_index][realimag]=buffer_sum[i2idx][time_c][f1f2idx][gammaidx][spin_index][realimag];
                                       }
                                     }
                                     if ((gamma_table_source[i1_source]==11) || (gamma_table_source[i1_source]==2)){
                                       mult_with_gamma5_matrix_adj_source(correlation_function[udli_count][time_c]);
                                     }
                                     if ((gamma_table_sink[i1_sink]==12) || (gamma_table_sink[i1_sink]==5) || (gamma_table_sink[i1_sink]==10)) {
                                       mult_with_gamma5_matrix_sink(correlation_function[udli_count][time_c]);
                                     }
/*                                     for (int spin_index=0; spin_index < spin1212_dimension; ++ spin_index){
                                       for (int realimag=0; realimag < 2; ++realimag ){
                                         correlation_function_symmetrized[udli_count][time_c][spin_index][realimag]=correlation_function[udli_count][time_c][spin_index][realimag];

                                       }
                                     }*/
                                   }

                                   udli_id_store_piND( udli_list, udli_count, &momtable_sink[ibeta][imu][icombination_sink][i1_sink][0], gamma_string_sink_member,&momtable_source[ibeta][imu][icombination_source][i1_source][0],&momtable_source[ibeta][imu][icombination_source][i1_source][3],gamma_string_source_member);
                                   udliid=udli_count;
                                   ++udli_count;
                                   if (udli_count==MAX_UDLI_NUM){
                                     fprintf(stderr,"# [piN2piN_projection_apply] Error in the size of udli %d list\n", udli_count );
                                     exit(1);
                                   }

                                 }

                               } else if ( (strcmp(g_twopoint_function_list[i2pt].particlename_sink, "piN") == 0) && (strcmp(g_twopoint_function_list[i2pt].particlename_source, "D") == 0) ){

                                 int pf1x=momtable_sink[ibeta][imu][icombination_sink][i1_sink][0];
                                 int pf1y=momtable_sink[ibeta][imu][icombination_sink][i1_sink][1];
                                 int pf1z=momtable_sink[ibeta][imu][icombination_sink][i1_sink][2];
                                 int pf2x=momtable_sink[ibeta][imu][icombination_sink][i1_sink][3];
                                 int pf2y=momtable_sink[ibeta][imu][icombination_sink][i1_sink][4];
                                 int pf2z=momtable_sink[ibeta][imu][icombination_sink][i1_sink][5];

                                 int pi1x=momtable_source[ibeta][imu][icombination_source][i1_source][0];
                                 int pi1y=momtable_source[ibeta][imu][icombination_source][i1_source][1];
                                 int pi1z=momtable_source[ibeta][imu][icombination_source][i1_source][2];



                                 udliid=udli_id_lookup_DpiN( udli_list, udli_count, &momtable_sink[ibeta][imu][icombination_sink][i1_sink][0], gamma_string_sink_member,&momtable_sink[ibeta][imu][icombination_sink][i1_sink][3],&momtable_source[ibeta][imu][icombination_source][i1_source][0],gamma_string_source_member);
                                 if (udliid==-1){
                                    snprintf ( tagname, 400, "/sx%.02dsy%.02dsz%.02dst%03d/gf25/pf2x%.02dpf2y%.02dpf2z%.02d/gf1%s/pf1x%.02dpf1y%.02dpf1z%.02d/gi1%s/pi1x%.02dpi1y%.02dpi1z%.02d/%s", source_coords_list[k][1],
                                                        source_coords_list[k][2],
                                                        source_coords_list[k][3],
                                                        source_coords_list[k][0],
                                                        pf2x,
                                                        pf2y,
                                                        pf2z,
                                                        gamma_string_sink_member,
                                                        pf1x,
                                                        pf1y,
                                                        pf1z,
                                                        gamma_string_source_member,
                                                        pi1x,
                                                        pi1y,
                                                        pi1z,
                                                        correlation_function_tagname_suffix);

                                   const int i_total_momentum=Ptot[0]*Ptot[0]+Ptot[1]*Ptot[1]+Ptot[2]*Ptot[2];

                                   snprintf ( filename, 400, "%s%04d_PX%.02dPY%.02dPZ%.02d_%s.h5",
                                                        filename_prefix,
                                                        Nconf,
                                                        momentum_orbit_pref[i_total_momentum][0],
                                                        momentum_orbit_pref[i_total_momentum][1],
                                                        momentum_orbit_pref[i_total_momentum][2],
                                                        correlation_function_filename_suffix);

                                   correlation_function[udli_count]=init_3level_dtable( g_twopoint_function_list[i2pt].T,g_twopoint_function_list[i2pt].d*g_twopoint_function_list[i2pt].d ,2 );

                                   udli_id_store_DpiN( udli_list, udli_count, &momtable_sink[ibeta][imu][icombination_sink][i1_sink][0], gamma_string_sink_member,&momtable_sink[ibeta][imu][icombination_sink][i1_sink][3],&momtable_source[ibeta][imu][icombination_source][i1_source][0],gamma_string_source_member);


                                   exitstatus = read_from_h5_file ( (void*)(correlation_function[udli_count][0][0]), filename, tagname, io_proc, 0 );
                                   if (exitstatus!=0){
                                     fprintf(stderr,"# [piN2piN_projection_apply] Error in opening file %s %s\n", filename, tagname );
                                     exit(1);
                                   }
                                   udliid=udli_count;
                                   ++udli_count;
                                   if (udli_count==MAX_UDLI_NUM){
                                     fprintf(stderr,"# [piN2piN_projection_apply] Error in the size of udli %d list\n", udli_count );
                                     exit(1);
                                   }

                                 }
                                 
                               } else {
                                  fprintf(stderr, "Projector for particle at sink %s and at source %s is not implemented yet",g_twopoint_function_list[i2pt].particlename_sink,g_twopoint_function_list[i2pt].particlename_source);
                                  exit(1);
                               }

                               //printf("udliid %d limit: spin12 %d t %d , sink %d source %d \n", udliid,spin1212dimension,g_twopoint_function_list[i2pt].T,Nps_Nreplica_sink[ibeta][imu][icombination_source][1],Nps_Nreplica_source[ibeta][imu][icombination_sink][1]);

                               for ( int i3_source =0; i3_source < spin1212dimension; ++i3_source ) {

                                 double proj_source_real=projection_coeff_c_ORT[ibeta][imu][icombination_source][i1_source*spin1dimension_source*spin1212dimension+i2_source*spin1212dimension+i3_source][nreplicum_source][0];
                                 double proj_source_imag=projection_coeff_c_ORT[ibeta][imu][icombination_source][i1_source*spin1dimension_source*spin1212dimension+i2_source*spin1212dimension+i3_source][nreplicum_source][1];

                                 if ( (abs(proj_source_real)<1e-8 ) && ( abs(proj_source_imag)<1e-8 )){
                                   continue;
                                 }

                                 for ( int i3_sink =0; i3_sink < spin1212dimension ; ++i3_sink ) {

                                   double proj_sink_real=projection_coeff_a_ORT[ibeta][imu][icombination_sink][nreplicum_sink][i1_sink*spin1dimension_sink*spin1212dimension+i2_sink*spin1212dimension+i3_sink][0];
                                   double proj_sink_imag=projection_coeff_a_ORT[ibeta][imu][icombination_sink][nreplicum_sink][i1_sink*spin1dimension_sink*spin1212dimension+i2_sink*spin1212dimension+i3_sink][1];

                                   if ( (abs(proj_sink_real)<1e-8 ) && ( abs(proj_sink_imag)<1e-8 )){
                                     continue;
                                   }


                                   for (int t=0; t<g_twopoint_function_list[i2pt].T; ++t){


                                     double tmp_real_part=+proj_sink_real*proj_source_real*correlation_function[udliid][t][i3_source+4*i3_sink][0]
                                                          -proj_sink_real*proj_source_imag*correlation_function[udliid][t][i3_source+4*i3_sink][1]
                                                          -proj_sink_imag*proj_source_imag*correlation_function[udliid][t][i3_source+4*i3_sink][0]
                                                          -proj_sink_imag*proj_source_real*correlation_function[udliid][t][i3_source+4*i3_sink][1];
                                     projected_correlation_function[t][0]+=tmp_real_part;

                                     double tmp_imag_part=-proj_sink_imag*proj_source_imag*correlation_function[udliid][t][i3_source+4*i3_sink][1]
	                                                  +proj_sink_real*proj_source_real*correlation_function[udliid][t][i3_source+4*i3_sink][1]
							  +proj_sink_real*proj_source_imag*correlation_function[udliid][t][i3_source+4*i3_sink][0]
							  +proj_sink_imag*proj_source_real*correlation_function[udliid][t][i3_source+4*i3_sink][0];
                                     projected_correlation_function[t][1]+=tmp_imag_part;
#if 0
                                     tmp_real_part=+proj_sink_real*proj_source_real*correlation_function_symmetrized[udliid][t][i3_source+4*i3_sink][0]
                                                          -proj_sink_real*proj_source_imag*correlation_function_symmetrized[udliid][t][i3_source+4*i3_sink][1]
                                                          -proj_sink_imag*proj_source_imag*correlation_function_symmetrized[udliid][t][i3_source+4*i3_sink][0]
                                                          -proj_sink_imag*proj_source_real*correlation_function_symmetrized[udliid][t][i3_source+4*i3_sink][1];
                                     projected_correlation_function_symmetrized[t][0]+=tmp_real_part;

                                     tmp_imag_part=-proj_sink_imag*proj_source_imag*correlation_function_symmetrized[udliid][t][i3_source+4*i3_sink][1]
                                                          +proj_sink_real*proj_source_real*correlation_function_symmetrized[udliid][t][i3_source+4*i3_sink][1]
                                                          +proj_sink_real*proj_source_imag*correlation_function_symmetrized[udliid][t][i3_source+4*i3_sink][0]
                                                          +proj_sink_imag*proj_source_real*correlation_function_symmetrized[udliid][t][i3_source+4*i3_sink][0];
                                     projected_correlation_function_symmetrized[t][1]+=tmp_imag_part;
#endif
/*
                                     if (t == 0) {
                                       fprintf(stdout,"# [piN2piN_projection_apply] (beta=%d) (mu=%d) (repl source %d) (repl sink %d) (mom_source_id %d) (gamma_source_id %d) (mom_sink_id %d) (gamma_sink_id %d) (inner spin source %d) (inner spin sink %d) %e %e\n", ibeta, imu, nreplicum_source, nreplicum_sink, i1_source, i2_source, i1_sink, i2_sink, i3_source, i3_sink, tmp_real_part, tmp_imag_part ); 
                                     }*/

                                   }/*time*/

                                 }/*i3 sink */

                               }/*i3_source */

                               free(correlation_function_filename_suffix);
     
                               free(correlation_function_tagname_suffix);
  
                               free(gamma_string_sink_member);

                             }/*i2_sink */

                           }/*i1_sink */

                           free(gamma_string_source_member);

                         }/*i2_source*/

                       }/*i1_source*/


                       snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d/mu_%d/source_totmomN%dp%d/sink_totmomN%dp%d/source_%s/sink_%s/Replicasource_%d/Replicasink_%d/data", source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta, imu, g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source],g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink], gamma_string_source, gamma_string_sink, nreplicum_source, nreplicum_sink);

                       hsize_t dims[2];
                       dims[0]=g_twopoint_function_list[i2pt].T;
                       dims[1]=2;
                       dataspace_id = H5Screate_simple(2, dims, NULL);


                       /* Create a dataset in group "MyGroup". */
                       dataset_id = H5Dcreate2(file_id, tagname, H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                       /* Write the first dataset. */
                       status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(projected_correlation_function[0][0]));

                       /* Close the data space for the first dataset. */
                       status = H5Sclose(dataspace_id);

                       /* Close the first dataset. */
                       status = H5Dclose(dataset_id);

                       fini_2level_dtable(&projected_correlation_function);
#if 0
                       snprintf ( tagname, 400, "/sx%0.2dsy%0.2dsz%0.2dst%03d/px%dpy%dpz%d/beta_%d/mu_%d/source_totmomN%dp%d/sink_totmomN%dp%d/source_%s/sink_%s/Replicasource_%d/Replicasink_%d/data_symmetrized",source_coords_list[k][1],
                         source_coords_list[k][2],
                         source_coords_list[k][3],
                         source_coords_list[k][0], Ptot[0],Ptot[1],Ptot[2], ibeta, imu, g_twopoint_function_list[i2pt].total_momentum_nucleon_source[icombination_source], g_twopoint_function_list[i2pt].total_momentum_pion_source[icombination_source],g_twopoint_function_list[i2pt].total_momentum_nucleon_sink[icombination_sink], g_twopoint_function_list[i2pt].total_momentum_pion_sink[icombination_sink], gamma_string_source, gamma_string_sink, nreplicum_source, nreplicum_sink);

                       dims[0]=g_twopoint_function_list[i2pt].T;
                       dims[1]=2;
                       dataspace_id = H5Screate_simple(2, dims, NULL);


                       /* Create a dataset in group "MyGroup". */
                       dataset_id = H5Dcreate2(file_id, tagname, H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                       /* Write the first dataset. */
                       status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(projected_correlation_function_symmetrized[0][0]));

                       /* Close the data space for the first dataset. */
                       status = H5Sclose(dataspace_id);

                       /* Close the first dataset. */
                       status = H5Dclose(dataset_id);

                       fini_2level_dtable(&projected_correlation_function_symmetrized);
#endif

                     } /*nreplicum sink */


                   } /*nreplicum source */

                   fini_1level_itable(&gamma_table_sink);
                   free(gamma_string_sink);


                 } /*ngamma multiplett sink*/

                 fini_1level_itable(&gamma_table_source);
                 free(gamma_string_source);

               } /*ngamma muplitplett source*/

             } /* momentum combination sink */
 
           } /* momentum combination source */

         } /*imu */

       } /*ibeta*/ 

       for (int iudli=0; iudli<udli_count; ++iudli){
         fini_3level_dtable(&correlation_function[udli_count]);
//         fini_3level_dtable(&correlation_function_symmetrized[udli_count]);
       }
       udli_count=0;



       fini_4level_itable(&Nps_Nreplica_source);
       fini_4level_itable(&Nps_Nreplica_sink);
       for (int ibeta=0; ibeta<1; ++ibeta){
         for (int imu=0; imu<irrep_dim; ++imu){
           for (int icombinations_source=0; icombinations_source < g_twopoint_function_list[i2pt].ncombination_total_momentum_source;icombinations_source++){
             fini_3level_dtable(&projection_coeff_c_ORT[ibeta][imu][icombinations_source]);
             fini_2level_itable(&momtable_source[ibeta][imu][icombinations_source]);
           }
           free(projection_coeff_c_ORT[ibeta][imu]);
           free(momtable_source[ibeta][imu]);
         }
         free(projection_coeff_c_ORT[ibeta]);
         free(momtable_source[ibeta]);
       }
       free(projection_coeff_c_ORT);
       free(momtable_source);
       for (int ibeta=0; ibeta<1; ++ibeta){
         for (int imu=0; imu<irrep_dim; ++imu){
           for (int icombinations_sink=0; icombinations_sink < g_twopoint_function_list[i2pt].ncombination_total_momentum_sink;icombinations_sink++){
             fini_2level_itable(&momtable_sink[ibeta][imu][icombinations_sink]);
             fini_3level_dtable(&projection_coeff_a_ORT[ibeta][imu][icombinations_sink]);
           }
           free(projection_coeff_a_ORT[ibeta][imu]);
           free(momtable_sink[ibeta][imu]);
         }
         free(projection_coeff_a_ORT[ibeta]);
         free(momtable_sink[ibeta]);
       }
       free(projection_coeff_a_ORT);
       free(momtable_sink);



       /* Close the file. */
       status = H5Fclose(file_id);

     } /* end of loop on 2-point functions*/

 
     fini_2level_itable(&buffer_mom);

     fini_6level_dtable(&buffer_sum);


   }  /* end of loop on source positions */



   /******************************************************/
   /******************************************************/
#if 0
#endif  /* of if 0 */

   /******************************************************
    * finalize
    *
    * free the allocated memory, finalize
   ******************************************************/
  free_geometry();

#ifdef HAVE_MPI
   MPI_Finalize();
#endif
   if(g_cart_id == 0) {
    g_the_time = time(NULL);
    fprintf(stdout, "# [piN2piN_projection_apply_plegma] %s# [piN2piN_projection_apply_plegma] end fo run\n", ctime(&g_the_time));
    fflush(stdout);
    fprintf(stderr, "# [piN2piN_projection_apply_plegma] %s# [piN2piN_projection_apply_plegma] end fo run\n", ctime(&g_the_time));
    fflush(stderr);
   }

  return(0);
}