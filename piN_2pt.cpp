/****************************************************
 * piN_2pt.c
 * 
 * Wed Nov 16 15:21:03 CET 2016
 *
 * PURPOSE:
 * - originally copied from proton_2pt

 * TODO:
 * DONE:
 *
 ****************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef HAVE_MPI
#  include <mpi.h>
#endif
#ifdef OPENMP
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

#  ifdef HAVE_TMLQCD_LIBWRAPPER
#    include "tmLQCD.h"
#  endif

#ifdef __cplusplus
}
#endif

#define MAIN_PROGRAM

#include "cvc_complex.h"
#include "ilinalg.h"
#include "icontract.h"
#include "global.h"
#include "cvc_geometry.h"
#include "cvc_utils.h"
#include "mpi_init.h"
#include "io.h"
#include "propagator_io.h"
#include "gauge_io.h"
#include "read_input_parser.h"
/* #include "smearing_techniques.h" */
#include "contractions_io.h"
#include "matrix_init.h"
#include "project.h"

using namespace cvc;

void usage() {
  fprintf(stdout, "Code to perform contractions for proton 2-pt. function\n");
  fprintf(stdout, "Usage:    [options]\n");
  fprintf(stdout, "Options: -f input filename [default cvc.input]\n");
  fprintf(stdout, "         -a write ascii output too [default no ascii output]\n");
  fprintf(stdout, "         -F fermion type, must be set [default -1, no type]\n");
  fprintf(stdout, "         -q/Q/p/P p[i,f][1,2] source and sink momenta [default 0]\n");
  fprintf(stdout, "         -h? this help\n");
#ifdef HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Finalize();
#endif
  exit(0);
}


int main(int argc, char **argv) {
  
  const int n_c=3;
  const int n_s=4;
  const char outfile_prefix[] = "piN_2pt";

  int c, i, k;
  int filename_set = 0;
  int exitstatus;
  int it, ir, is;
  int gsx[4], sx0, sx1, sx2, sx3;
  int write_ascii=0;
  int fermion_type = -1;
  int write_xspace = 0;
  int source_proc_id = 0, source_proc_coords[4];
  char filename[200], contype[1200];
  double ratime, retime;
  double plaq_m, plaq_r;
  double *spinor_work[2];
  unsigned int ix, iix;
  unsigned int VOL3;
  size_t sizeof_spinor_field = 0;
  spinor_propagator_type *connq=NULL;
  double ****connt[10], **connt_p=NULL, **connt_n=NULL;
  double **buffer=NULL;
  int io_proc = -1;

#ifdef HAVE_LHPC_AFF
  struct AffWriter_s *affw = NULL;
  struct AffNode_s *affn = NULL, *affdir=NULL;
  char * aff_status_str;
  double _Complex *aff_buffer = NULL;
  char aff_buffer_path[200];
  /*  uint32_t aff_buffer_size; */
#endif

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  while ((c = getopt(argc, argv, "xah?f:F:")) != -1) {
    switch (c) {
    case 'f':
      strcpy(filename, optarg);
      filename_set=1;
      break;
    case 'F':
      if(strcmp(optarg, "Wilson") == 0) {
        fermion_type = _WILSON_FERMION;
      } else if(strcmp(optarg, "tm") == 0) {
        fermion_type = _TM_FERMION;
      } else {
        fprintf(stderr, "[piN_2pt] Error, unrecognized fermion type\n");
        exit(145);
      }
      fprintf(stdout, "# [piN_2pt] will use fermion type %s ---> no. %d\n", optarg, fermion_type);
      break;
    case 'x':
      write_xspace = 1;
      fprintf(stdout, "# [piN_2pt] will write xspace correlator\n");
      break;
    case 'h':
    case '?':
    default:
      usage();
      break;
    }
  }

  /* set the default values */
  if(filename_set==0) strcpy(filename, "cvc.input");
  fprintf(stdout, "# reading input from file %s\n", filename);
  read_input_parser(filename);

  if(fermion_type == -1 ) {
    fprintf(stderr, "# [piN_2pt] fermion_type must be set\n");
    exit(1);
  }

#ifdef HAVE_TMLQCD_LIBWRAPPER

  fprintf(stdout, "# [piN_2pt] calling tmLQCD wrapper init functions\n");

  /*********************************
   * initialize MPI parameters for cvc
   *********************************/
  exitstatus = tmLQCD_invert_init(argc, argv, 1);
  if(exitstatus != 0) {
    EXIT(14);
  }
  exitstatus = tmLQCD_get_mpi_params(&g_tmLQCD_mpi);
  if(exitstatus != 0) {
    EXIT(15);
  }
  exitstatus = tmLQCD_get_lat_params(&g_tmLQCD_lat);
  if(exitstatus != 0) {
    EXIT(16);
  }
#endif

#ifdef OPENMP
  omp_set_num_threads(g_num_threads);
#else
  fprintf(stdout, "[piN_2pt] Warning, resetting global thread number to 1\n");
  g_num_threads = 1;
#endif

  /* initialize MPI parameters */
  mpi_init(argc, argv);

  /******************************************************
   * initialize geometry
   ******************************************************/

  if(init_geometry() != 0) {
    fprintf(stderr, "[piN_2pt] Error from init_geometry\n");
    EXIT(1);
  }

  geometry();

  VOL3 = LX*LY*LZ;
  sizeof_spinor_field = 24*VOLUME*sizeof(double);

#ifdef HAVE_MPI
  /***********************************************
   * set io process
   ***********************************************/
  if( g_proc_coords[0] == 0 && g_proc_coords[1] == 0 && g_proc_coords[2] == 0 && g_proc_coords[3] == 0) {
    io_proc = 2;
    fprintf(stdout, "# [piN_2pt] proc%.4d is io process\n", g_cart_id);
  } else {
    if( g_proc_coords[1] == 0 && g_proc_coords[2] == 0 && g_proc_coords[3] == 0) {
      io_proc = 1;
      fprintf(stdout, "# [piN_2pt] proc%.4d is send process\n", g_cart_id);
    } else {
      io_proc = 0;
    }
  }
#else
  io_proc = 2;
#endif


#ifndef HAVE_TMLQCD_LIBWRAPPER
  /******************************************************
   * allocate and read the gauge field
   ******************************************************/
  alloc_gauge_field(&g_gauge_field, VOLUMEPLUSRAND);
  switch(g_gauge_file_format) {
    case 0:
      sprintf(filename, "%s.%.4d", gaugefilename_prefix, Nconf);
      if(g_cart_id==0) fprintf(stdout, "reading gauge field from file %s\n", filename);
      exitstatus = read_lime_gauge_field_doubleprec(filename);
      break;
    case 1:
      sprintf(filename, "%s.%.5d", gaugefilename_prefix, Nconf);
      if(g_cart_id==0) fprintf(stdout, "\n# [piN_2pt] reading gauge field from file %s\n", filename);
      exitstatus = read_nersc_gauge_field(g_gauge_field, filename, &plaq_r);
      break;
  }
  if(exitstatus != 0) {
    fprintf(stderr, "[piN_2pt] Error, could not read gauge field\n");
    EXIT(21);
  }
#  ifdef HAVE_MPI
  xchange_gauge();
#  endif
#else
  /******************************************************
   * retrieve gauge field from tmLQCD
   ******************************************************/
  Nconf = g_tmLQCD_lat.nstore;
  if(g_cart_id== 0) fprintf(stdout, "[piN_2pt] Nconf = %d\n", Nconf);

  exitstatus = tmLQCD_read_gauge(Nconf);
  if(exitstatus != 0) {
    EXIT(3);
  }

  exitstatus = tmLQCD_get_gauge_field_pointer(&g_gauge_field);
  if(exitstatus != 0) {
    EXIT(4);
  }
  if(&g_gauge_field == NULL) {
    fprintf(stderr, "[piN_2pt] Error, &g_gauge_field is NULL\n");
    EXIT(5);
  }
#endif

  /* measure the plaquette */
  plaquette(&plaq_m);
  if(g_cart_id==0) {
    if(g_gauge_file_format == 1) fprintf(stdout, "# [piN_2pt] read plaquette value    : %25.16e\n", plaq_r);
    fprintf(stdout, "# [piN_2pt] measured plaquette value: %25.16e\n", plaq_m);
  }


  /***********************************************************
   * allocate memory for the spinor fields
   *
   * of size VOLUME
   * u, (d,) tfii
   * for one timeslice tffi xor pffii xor pfifi
   ***********************************************************/
  g_spinor_field = NULL;
  no_fields = n_s*n_c;
  if(fermion_type == _TM_FERMION) {
    no_fields *= 4;
  } else {
    no_fields *= 3;
  }
  no_fields += 2;

  g_spinor_field = (double**)calloc(no_fields, sizeof(double*));
  for(i=0; i<no_fields-2; i++) alloc_spinor_field(&g_spinor_field[i], VOLUME);
  alloc_spinor_field(&g_spinor_field[no_fields-2], VOLUMEPLUSRAND);
  alloc_spinor_field(&g_spinor_field[no_fields-1], VOLUMEPLUSRAND);
  spinor_work[0] = g_spinor_field[no_fields-2];
  spinor_work[1] = g_spinor_field[no_fields-1];

  /***********************************************************
   * allocate memory for the contractions
   **********************************************************/
  connq = create_sp_field( (size_t)VOL3 );
  if(connq == NULL) {
    fprintf(stderr, "[piN_2pt] Error, could not alloc connq\n");
    EXIT(2);
  }

  /***********************************************************
   * determine source coordinates, find out, if source_location is in this process
   ***********************************************************/
  /* global source coordinates */
  gsx[0] = g_source_location / ( LX_global * LY_global * LZ_global);
  gsx[1] = (g_source_location % ( LX_global * LY_global * LZ_global)) / (LY_global * LZ_global);
  gsx[2] = (g_source_location % ( LY_global * LZ_global)) / LZ_global;
  gsx[3] = (g_source_location % LZ_global);
  /* local source coordinates */
  sx0 = gsx[0] % T;
  sx1 = gsx[1] % LX;
  sx2 = gsx[2] % LY;
  sx3 = gsx[3] % LZ;
  source_proc_id = 0;
#ifdef HAVE_MPI
  source_proc_coords[0] = gsx[0] / T;
  source_proc_coords[1] = gsx[1] / LX;
  source_proc_coords[2] = gsx[2] / LY;
  source_proc_coords[3] = gsx[3] / LZ;

  if(g_cart_id == 0) {
    fprintf(stdout, "# [piN_2pt] global source coordinates: (%3d,%3d,%3d,%3d)\n",  gsx[0], gsx[1], gsx[2], gsx[3]);
    fprintf(stdout, "# [piN_2pt] source proc coordinates: (%3d,%3d,%3d,%3d)\n",  source_proc_coords[0], source_proc_coords[1], source_proc_coords[2], source_proc_coords[3]);
  }

  exitstatus = MPI_Cart_rank(g_cart_grid, source_proc_coords, &source_proc_id);
  if(exitstatus !=  MPI_SUCCESS ) {
    fprintf(stderr, "[piN_2pt] Error from MPI_Cart_rank, status was %d\n", exitstatus);
    EXIT(9);
  }
#endif
  if( source_proc_id == g_cart_id) {
    fprintf(stdout, "# [piN_2pt] process %2d has source location\n", source_proc_id);
  }


  /***********************************************************
   * up-type propagator
   ***********************************************************/
  ratime = _GET_TIME;
  for(is=0;is<n_s*n_c;is++) {

    memset(spinor_work[0], 0, sizeof_spinor_field);
    memset(spinor_work[1], 0, sizeof_spinor_field);
    if(source_proc_id == g_cart_id)  {
      spinor_work[0][_GSI(g_ipt[sx0][sx1][sx2][sx3])+2*is] = 1.;
    }

    exitstatus = tmLQCD_invert(spinor_work[1], spinor_work[0], 0, 0);
    if(exitstatus != 0) {
      fprintf(stderr, "[piN_2pt] Error from tmLQCD_invert, status was %d\n", exitstatus);
      EXIT(12);
    }
    memcpy( g_spinor_field[is], spinor_work[1], sizeof_spinor_field);
  }
  retime = _GET_TIME;
  if(g_cart_id == 0) fprintf(stderr, "# [piN_2pt] time for up propagator = %e seconds\n", retime-ratime);

  /***********************************************************
   * dn-type propagator
   ***********************************************************/
  if(fermion_type == _TM_FERMION) {
    ratime = _GET_TIME;
    for(is=0;is<n_s*n_c;is++) {

      memset(spinor_work[0], 0, sizeof_spinor_field);
      memset(spinor_work[1], 0, sizeof_spinor_field);
      if(source_proc_id == g_cart_id)  {
        spinor_work[0][_GSI(g_ipt[sx0][sx1][sx2][sx3])+2*is] = 1.;
      }

      exitstatus = tmLQCD_invert(spinor_work[1], spinor_work[0], 1, 0);
      if(exitstatus != 0) {
        fprintf(stderr, "[piN_2pt] Error from tmLQCD_invert, status was %d\n", exitstatus);
        EXIT(12);
      }
      memcpy( g_spinor_field[n_s*n_c+is], spinor_work[1], sizeof_spinor_field);
    }
    retime = _GET_TIME;
    if(g_cart_id == 0) fprintf(stderr, "# [piN_2pt] time for dn propagator = %e seconds\n", retime-ratime);
  }

  /***********************************************************
   * Tud sequential propagator
   ***********************************************************/
  ratime = _GET_TIME;
  for(is=0;is<n_s*n_c;is++) {

    memset(spinor_work[0], 0, sizeof_spinor_field);
    memset(spinor_work[1], 0, sizeof_spinor_field);

    exitstatus = init_sequential_source(spinor_work[0], g_spinor_field[is], gsx[0], g_mxb_mxb_2pt_list[0].pi2, 5);
    if(exitstatus != 0) {
      fprintf(stderr, "[piN_2pt] Error from init_sequential_source, status was %d\n", exitstatus);
      EXIT(12);
    }

    exitstatus = tmLQCD_invert(spinor_work[1], spinor_work[0], 1, 0);
    if(exitstatus != 0) {
      fprintf(stderr, "[piN_2pt] Error from tmLQCD_invert, status was %d\n", exitstatus);
      EXIT(12);
    }
    memcpy( g_spinor_field[n_s*n_c*(1+(fermion_type == _TM_FERMION))+is], spinor_work[1], sizeof_spinor_field);
  }
  retime = _GET_TIME;
  if(g_cart_id == 0) fprintf(stderr, "# [piN_2pt] time for T_ud propagator = %e seconds\n", retime-ratime);
  

  /******************************************************
   * connt
   ******************************************************/
  for(i=0; i<10; i++) {
    init_4level_buffer(connt+i, T, g_sink_momentum_number, g_sv_dim, 2*g_sv_dim);
  }


  /******************************************************
   * contractions
   ******************************************************/

  for(dt=0; dt<=dt_source_sink_max; dt++) {
    sink_timeslice = ( gsx[0] + dt ) %T_global;
    if(g_cart_id == 0) fprintf(stdout, "# [] processing sink timeslice %d\n", sink_timeslice);

    /******************************************************
     * B diagrams
     ******************************************************/

    /***********************************************************
     * Tud sequential propagator
     ***********************************************************/
    ratime = _GET_TIME;
    for(is=0;is<n_s*n_c;is++) {

      memset(spinor_work[0], 0, sizeof_spinor_field);
      memset(spinor_work[1], 0, sizeof_spinor_field);

      exitstatus = init_sequential_source(spinor_work[0], g_spinor_field[is+n_s*n_c*(1+(fermion_type == _TM_FERMION))], sink_timeslice, g_mxb_mxb_2pt_list[0].pf2, 5);
      if(exitstatus != 0) {
        fprintf(stderr, "[piN_2pt] Error from init_sequential_source, status was %d\n", exitstatus);
        EXIT(12);
      }

      exitstatus = tmLQCD_invert(spinor_work[1], spinor_work[0], 1, 0);
      if(exitstatus != 0) {
        fprintf(stderr, "[piN_2pt] Error from tmLQCD_invert, status was %d\n", exitstatus);
        EXIT(12);
      }
      memcpy( g_spinor_field[n_s*n_c*(2+(fermion_type == _TM_FERMION))+is], spinor_work[1], sizeof_spinor_field);
    }
    retime = _GET_TIME;
    if(g_cart_id == 0) fprintf(stderr, "# [piN_2pt] time for Pffii propagator = %e seconds\n", retime-ratime);

    ratime = _GET_TIME;
#ifdef HAVE_OPENMP
#pragma omp parallel 
{
#endif
  /* variables */
  fermion_propagator_type fp1, fp2, fp3, uprop, dprop;
  spinor_propagator_type sp1, sp2;
#ifdef HAVE_OPENMP
#pragma omp for
#endif
  for(ix=0; ix<VOLUME; ix++) {

    // assign the propagators
    _assign_fp_point_from_field(uprop, g_spinor_field, ix);
    if(fermion_type==_TM_FERMION) {
      _assign_fp_point_from_field(dprop, g_spinor_field+n_s*n_c, ix);
    } else {
      _fp_eq_fp(dprop, uprop);
    }

    // flavor rotation for twisted mass fermions
    if(fermion_type == _TM_FERMION) {
      _fp_eq_rot_ti_fp(fp1, uprop, +1, fermion_type, fp2);
      _fp_eq_fp_ti_rot(uprop, fp1, +1, fermion_type, fp2);
      _fp_eq_rot_ti_fp(fp1, dprop, -1, fermion_type, fp2);
      _fp_eq_fp_ti_rot(dprop, fp1, -1, fermion_type, fp2);
    }

    // S_u x Cg5
    _fp_eq_fp_ti_Cg5(fp1, uprop, fp3);

    // Cg5 x S_d
    _fp_eq_Cg5_ti_fp(fp2, dprop, fp3);
    
    /******************************************************
     * first contribution
     ******************************************************/

    // reduce
    _fp_eq_zero(fp3);
    _fp_eq_fp_eps_contract13_fp(fp3, fp1, fp2);

    // reduce to spin propagator
    _sp_eq_zero( sp1 );
    _sp_eq_fp_del_contract34_fp(sp1, uprop, fp3);

    /******************************************************
     * second contribution
     ******************************************************/

    // reduce
    _fp_eq_zero(fp3);
    _fp_eq_fp_eps_contract24_fp(fp3, fp1, fp2);

    // reduce to spin propagator
    _sp_eq_zero( sp2 );
    _sp_eq_fp_del_contract23_fp(sp2, fp3, uprop);

    /******************************************************
     * add both contribution, write to field
     * - the second contribution receives additional minus
     *   from rearranging the epsilon tensor indices
     ******************************************************/

    _sp_pl_eq_sp(sp1, sp2);

    _sp_eq_sp( connq[ix], sp1);
  }  // end of loop on VOLUME
#ifdef OPENMP
}  /* end of parallel region */
#endif
  retime = _GET_TIME;
  if(g_cart_id == 0)  fprintf(stdout, "# [piN_2pt] time for contractions = %e seconds\n", retime-ratime);

  /***********************************************
   * free gauge fields and spinor fields
   ***********************************************/
  if(g_gauge_field != NULL) {
    free(g_gauge_field);
    g_gauge_field=(double*)NULL;
  }
  if(g_spinor_field!=NULL) {
    for(i=0; i<no_fields; i++) free(g_spinor_field[i]);
    free(g_spinor_field); g_spinor_field=(double**)NULL;
  }


  /***********************************************
   * finish calculation of connq
   ***********************************************/
  ratime = _GET_TIME;
  if(g_propagator_bc_type == 0) {
    // multiply with phase factor
    fprintf(stdout, "# [piN_2pt] multiplying with boundary phase factor\n");
    iix = 0;
    for(it=0;it<T_global;it++) {
      ir = (it - gsx[0] + T_global) % T_global;
      const complex w1 = { cos( 3. * M_PI*(double)ir / (double)T_global ), sin( 3. * M_PI*(double)ir / (double)T_global ) };
#ifdef HAVE_OPENMP
#pragma omp parallel private(ix,iix) shared(connq)
{
#endif
      spinor_propagator_type sp1;
#ifdef HAVE_OPENMP
#pragma omp for
#endif
      for(ix=0;ix<VOL3;ix++) {
        _sp_eq_sp(sp1, connq[iix] );
        _sp_eq_sp_ti_co( connq[iix], sp1, w1);
        iix++;
      }
#ifdef HAVE_OPENMP
}  /* end of parallel region */
#endif
    }
  } else if (g_propagator_bc_type == 1) {
    // multiply with step function
    fprintf(stdout, "# [piN_2pt] multiplying with boundary step function\n");
    for(ir=0; ir<T; ir++) {
      it = ir + g_proc_coords[0] * T;  // global t-value, 0 <= t < T_global
      if(it < gsx[0]) {
#ifdef HAVE_OPENMP
#pragma omp parallel private(ix,iix) firstprivate(it)
{
#endif
        spinor_propagator_type sp1;
        iix = it * VOL3;
#ifdef HAVE_OPENMP
#pragma omp for
#endif
        for(ix=0;ix<VOL3;ix++) {
          _sp_eq_sp(sp1, connq[iix] );
          _sp_eq_sp_ti_re( connq[iix], sp1, -1.);
          iix++;
        }
#ifdef HAVE_OPENMP
}  /* end of parallel region */
#endif
      }
    }  /* end of if it < gsx[0] */
  }
  retime = _GET_TIME;
  if(g_cart_id == 0)  fprintf(stdout, "# [piN_2pt] time for boundary phase = %e seconds\n", retime-ratime);

  if(write_ascii) {
    /***********************************************
     * each MPI process dump his part in ascii format
     ***********************************************/
    ratime = _GET_TIME;
    sprintf(filename, "%s_x.%.4d.t%.2dx%.2dy%.2dz%.2d.proct%.2dprocx%.2dprocy%.2dprocz%.2d.ascii", outfile_prefix, Nconf, gsx[0], gsx[1], gsx[2], gsx[3],
       g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3]);
    FILE *ofs = fopen(filename, "w");
    if(ofs == NULL) {
      fprintf(stderr, "[piN_2pt] Error opening file %s\n", filename);
      EXIT(56);
    }
    for(ix=0; ix<VOLUME; ix++) {
      sprintf(contype, "# ix = %u", ix);
      printf_sp(connq[ix], contype, ofs);
    }
    fclose(ofs);
    retime = _GET_TIME;
    if(g_cart_id == 0)  fprintf(stdout, "# [piN_2pt] time for writing ascii = %e seconds\n", retime-ratime);
  }

  /***********************************************
   * write to file
   ***********************************************/
  if(write_xspace) {
    ratime = _GET_TIME;
    sprintf(contype, "\n<description> proton 2pt spinor propagator position space\n"\
      "<components>%dx%d</components>\n"\
      "<data_type>%s</data_type>\n"\
      "<precision>%d</precision>\n"\
      "<source_coords_t>%2d</source_coords_t>\n"\
      "<source_coords_x>%2d</source_coords_x>\n"\
      "<source_coords_y>%2d</source_coords_y>\n"\
      "<source_coords_z>%2d</source_coords_z>\n"\
      "<spin_structure>Cg5-Cg5</spin_structure>\n",\
      g_sv_dim, g_sv_dim, "complex", 64, gsx[0], gsx[1], gsx[2], gsx[3]);

    sprintf(filename, "%s_x.%.4d.t%.2dx%.2dy%.2dz%.2d", outfile_prefix, Nconf, gsx[0], gsx[1], gsx[2], gsx[3]);
    write_lime_contraction(connq[0][0], filename, 64, g_sv_dim*g_sv_dim, contype, Nconf, 0);
    retime = _GET_TIME;
    if(g_cart_id == 0) {
      fprintf(stdout, "# [piN_2pt] time for writing xspace = %e seconds\n", retime-ratime);
    }
  }

  /***********************************************
   * momentum projections
   ***********************************************/
  init_4level_buffer(&connt, T, g_sink_momentum_number, g_sv_dim, 2*g_sv_dim);
  for(it=0; it<T; it++) {
    exitstatus = momentum_projection2 (connq[it*VOL3][0], connt[it][0][0], g_sv_dim*g_sv_dim, g_sink_momentum_number, g_sink_momentum_list, &(gsx[1]) );
  }

  init_2level_buffer(&connt_p, g_sink_momentum_number, 2*T);
  init_2level_buffer(&connt_n, g_sink_momentum_number, 2*T);
#ifdef HAVE_OPENMP
#pragma omp parallel
{
#endif
#ifdef HAVE_OPENMP
#pragma omp for
#endif
  for(it=0; it<T; it++) {
    spinor_propagator_type sp1, sp2;
    complex w;
    for(k=0; k<g_sink_momentum_number; k++) {
      _sp_eq_sp(sp1, connt[it][k]);
      _sp_eq_gamma_ti_sp(sp2, 0, sp1);
      _sp_pl_eq_sp(sp1, sp2);
      _co_eq_tr_sp(&w, sp1);
      connt_p[k][2*it  ] = w.re * 0.25;
      connt_p[k][2*it+1] = w.im * 0.25;
      _sp_eq_sp(sp1, connt[it][k]);
      _sp_eq_gamma_ti_sp(sp2, 0, sp1);
      _sp_mi_eq_sp(sp1, sp2);
      _co_eq_tr_sp(&w, sp1);
      connt_n[k][2*(T+it)  ] = w.re * 0.25;
      connt_n[k][2*(T+it)+1] = w.im * 0.25;
    }  /* end of loop on sink momenta */
  }  /* end of loop on T */
#ifdef HAVE_OPENMP
}  /* end of parallel region */
#endif
  fini_4level_buffer(&connt);

  ratime = _GET_TIME;
#ifdef HAVE_MPI
  if(io_proc == 2) {
    init_2level_buffer(&buffer, g_sink_mommentum_number, 2*T_global);
  }
  if(io_proc>0) {
    exitstatus = MPI_Gather(connt_p[0], 2*g_sink_mommentum_number*T, MPI_DOUBLE, buffer[0],  2*g_sink_mommentum_number*T, MPI_DOULBE, 0, g_tr_comm);
    if(exitstatus != MPI_SUCCESS) {
      EXIT(124);
    }
  }
#else
  buffer = connt_p;
#endif
  if(io_proc == 2) {
#ifdef HAVE_LHPC_AFF
    for(k=0; k<g_sink_momentum_number; k++) {
      sprintf(aff_buffer_path, "/%s/P+/px%.2dpy%.2dpz%.2d/t%.2dx%.2dy%.2dz%.2d", outfile_prefix, 
          g_sink_momentum_list[k][0], g_sink_momentum_list[k][1], g_sink_momentum_list[k][2],
          gsx[0], gsx[1], gsx[2], gsx[3]);
      fprintf(stdout, "# [piN_2pt] current aff path = %s\n", aff_buffer_path);
      affdir = aff_writer_mkpath(affw, affn, aff_buffer_path);
      memcpy(aff_buffer, buffer[k], 2*T_global*sizeof(double));
      int status = aff_node_put_complex (affw, affdir, aff_buffer, (uint32_t)T_global);
      if(status != 0) {
        fprintf(stderr, "[piN_2pt] Error from aff_node_put_double, status was %d\n", status);
        EXIT(8);
      }
    }
#endif
  }

#ifdef HAVE_MPI
  if(io_proc>0) {
    exitstatus = MPI_Gather(connt_n[0], 2*g_sink_mommentum_number*T, MPI_DOUBLE, buffer[0],  2*g_sink_mommentum_number*T, MPI_DOULBE, 0, g_tr_comm);
    if(exitstatus != MPI_SUCCESS) {
      EXIT(124);
    }
  }
#else
  buffer = connt_n;
#endif
  if(io_proc == 2) {
#ifdef HAVE_LHPC_AFF
    for(k=0; k<g_sink_momentum_number; k++) {
      sprintf(aff_buffer_path, "/%s/P-/px%.2dpy%.2dpz%.2d/t%.2dx%.2dy%.2dz%.2d", outfile_prefix, 
          g_sink_momentum_list[k][0], g_sink_momentum_list[k][1], g_sink_momentum_list[k][2],
          gsx[0], gsx[1], gsx[2], gsx[3]);
      fprintf(stdout, "# [piN_2pt] current aff path = %s\n", aff_buffer_path);
      affdir = aff_writer_mkpath(affw, affn, aff_buffer_path);
      memcpy(aff_buffer, connt_n[k], 2*T_global*sizeof(double));
      int status = aff_node_put_complex (affw, affdir, aff_buffer, (uint32_t)T_global);
      if(status != 0) {
        fprintf(stderr, "[piN_2pt] Error from aff_node_put_double, status was %d\n", status);
        EXIT(8);
      }
    }  /* end of loop on sink momenta */
#endif
  }  /* end of if io_proc == 2 */
  retime = _GET_TIME;
  if(io_proc == 2) fprintf(stdout, "# [piN_2pt] time for writing = %e seconds\n", retime - ratime);

#ifdef HAVE_MPI
  if(io_proc == 2) {
    fini_2level_buffer(&buffer);
  }
#endif
#ifdef HAVE_LHPC_AFF
  if(io_proc == 2) {
    aff_status_str = (char*)aff_writer_close (affw);
    if( aff_status_str != NULL ) {
      fprintf(stderr, "[piN_2pt] Error from aff_writer_close, status was %s\n", aff_status_str);
      return(11);
    }
    if(aff_buffer != NULL) free(aff_buffer);
  }  /* end of if io_proc == 2 */
#endif  /* of ifdef HAVE_LHPC_AFF */


  fini_2level_buffer(&connt_p);
  fini_2level_buffer(&connt_n);

  /***********************************************
   * free the allocated memory, finalize
   ***********************************************/
  free_geometry();
  free_sp_field(&connq);

#ifdef HAVE_TMLQCD_LIBWRAPPER
    tmLQCD_finalise();
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  g_the_time = time(NULL);
  fprintf(stdout, "# [piN_2pt] %s# [piN_2pt] end fo run\n", ctime(&g_the_time));
  fflush(stdout);
  fprintf(stderr, "# [piN_2pt] %s# [piN_2pt] end fo run\n", ctime(&g_the_time));
  fflush(stderr);

  return(0);
}