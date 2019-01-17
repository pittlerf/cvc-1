/****************************************************
 * loop_analyse
 *
 * PURPOSE:
 * DONE:
 * TODO:
 ****************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef HAVE_MPI
#  include <mpi.h>
#endif
#ifdef HAVE_OPENMP
#  include <omp.h>
#endif
#include <getopt.h>

#ifdef HAVE_LHPC_AFF
#include "lhpc-aff.h"
#endif

#define MAIN_PROGRAM

#include "cvc_complex.h"
#include "cvc_linalg.h"
#include "global.h"
#include "cvc_geometry.h"
#include "cvc_utils.h"
#include "mpi_init.h"
#include "set_default.h"
#include "io.h"
#include "propagator_io.h"
#include "read_input_parser.h"
#include "contractions_io.h"
#include "Q_clover_phi.h"
#include "contract_cvc_tensor.h"
#include "prepare_source.h"
#include "prepare_propagator.h"
#include "project.h"
#include "table_init_z.h"
#include "table_init_d.h"
#include "dummy_solver.h"
#include "Q_phi.h"
#include "clover.h"
#include "contract_loop.h"
#include "ranlxd.h"

#define _OP_ID_UP 0
#define _OP_ID_DN 1
#define _OP_ID_ST 2

using namespace cvc;

void usage() {
  fprintf(stdout, "Code to analyse loop data\n");
  fprintf(stdout, "Usage:    [options]\n");
  fprintf(stdout, "Options:  -f <input filename> : input filename for cvc      [default cpff.input]\n");
  EXIT(0);
}

int main(int argc, char **argv) {
  
  const char outfile_prefix[] = "loop";

  /* const char fbwd_str[2][4] =  { "fwd", "bwd" }; */

  int const conf_traj = 4;  /* meaning of this ? from HMC run ? */

  int c;
  int filename_set = 0;
  int exitstatus;
  int io_proc = -1;
  char filename[400];
  int Qsq = -1;

  struct timeval ta, tb;
  long unsigned int seconds, useconds;

  char output_filename[400];

  char data_tag[400];

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  while ((c = getopt(argc, argv, "h?f:Q:")) != -1) {
    switch (c) {
    case 'f':
      strcpy(filename, optarg);
      filename_set=1;
      break;
    case 'Q':
      Qsq = atoi ( optarg );
      break;
    case 'h':
    case '?':
    default:
      usage();
      break;
    }
  }

  g_the_time = time(NULL);

  /* set the default values */
  if(filename_set==0) snprintf ( filename, 400, "%s.input", outfile_prefix );
  /* fprintf(stdout, "# [loop_analyse] Reading input from file %s\n", filename); */
  read_input_parser(filename);

  /***************************************************************************
   * initialize MPI parameters for cvc
   ***************************************************************************/
  mpi_init(argc, argv);

  /***************************************************************************
   * report git version
   ***************************************************************************/
  if ( g_cart_id == 0 ) {
    fprintf(stdout, "# [loop_analyse] git version = %s\n", g_gitversion);
  }


  /***************************************************************************
   * set number of openmp threads
   ***************************************************************************/
  set_omp_number_threads ();

  /***************************************************************************
   * initialize geometry fields
   ***************************************************************************/
  if ( init_geometry() != 0 ) {
    fprintf(stderr, "[loop_analyse] Error from init_geometry %s %d\n", __FILE__, __LINE__);
    EXIT(4);
  }

  geometry();

  /***************************************************************************
   * some additional xchange operations
   ***************************************************************************/
  mpi_init_xchange_contraction(2);

  /***************************************************************************
   * set io process
   ***************************************************************************/
  io_proc = get_io_proc ();
  if( io_proc < 0 ) {
    fprintf(stderr, "[loop_analyse] Error, io proc must be ge 0 %s %d\n", __FILE__, __LINE__);
    EXIT(14);
  }
  fprintf(stdout, "# [loop_analyse] proc%.4d has io proc id %d\n", g_cart_id, io_proc );
  if (g_filtered_qsq > Qsq){
    if (g_proc_id == 0){
     fprintf(stderr, "# [loop_analyse] We do not have momenta till %f\n",g_filtered_qsq);
     EXIT(1);
    }
  }
  if (g_proc_id == 0){
    fprintf(stdout, "# [loop_analyse] Filtering value of Qsq %f\n", g_filtered_qsq);
  }
  if (g_loop_number == 0){
   if (g_proc_id == 0){
    fprintf(stderr, "# [loop_analyse] Specify the type of loops you are interested in e.g.\n");
    fprintf(stderr, "# [loop_analyse] filterlooptype = Scalar, Loops\n");
    EXIT(1);
   }
  }
  if (g_proc_id == 0)
   fprintf(stdout, "# [loop_analyse] Following loops will be filtered\n");
  for (int i=0; i<g_loop_number; ++i){
    switch( g_loop_type[i] ){
     case 0: if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] Scalar\n");
             break;
     case 1: if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] dOp\n");
             break;
     case 2: if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] Loops\n");
             break;
     case 3: if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] LpsDw\n");
             break;
     case 4: if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] LoopsCV\n");
             break;
     case 5: if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] LpsDwCv\n");
             break;
    }
  }

  /***************************************************************************
   * loop data filename
   ***************************************************************************/
  snprintf ( filename, 400, "%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d.h5", filename_prefix, Nconf, filename_prefix2, g_nsample, Nsave, Qsq );
  if ( io_proc == 2 && g_verbose > 2 ) fprintf ( stdout, "# [loop_analyse] loop filename = %s\n", filename );

  /***************************************************************************
   * count momenta and build momentum list
   ***************************************************************************/
  g_sink_momentum_number = 0;
  for( int x1 = -LX_global/2+1; x1 < LX_global/2; x1++ ) {
  for( int x2 = -LY_global/2+1; x2 < LY_global/2; x2++ ) {
  for( int x3 = -LZ_global/2+1; x3 < LZ_global/2; x3++ ) {
    int const qq = x1*x1 + x2*x2 + x3*x3;
    if ( qq <= Qsq ) {
      /*
      g_sink_momentum_list[g_sink_momentum_number][0] = x1;
      g_sink_momentum_list[g_sink_momentum_number][1] = x2;
      g_sink_momentum_list[g_sink_momentum_number][2] = x3;
      */
      g_sink_momentum_number++;
    }
  }}}
  if ( g_sink_momentum_number <= 0 ) {
    fprintf ( stderr, "[loop_analyse] Error, momentum list is empty %s %d\n", __FILE__, __LINE__ );
    EXIT(1);
  } else {
    if (io_proc == 2 && g_verbose > 1 ) fprintf ( stdout, "# [loop_analyse] number of momenta <= %3d is %3d\n", Qsq, g_sink_momentum_number );
  }

  exitstatus = loop_get_momentum_list_from_h5_file ( g_sink_momentum_list, filename, g_sink_momentum_number, io_proc );
  if ( exitstatus != 0 ) {
    fprintf ( stderr, "[] Error from loop_get_momentum_list_from_h5_file, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
    EXIT(1);

  }

  int filtered_sink_momentum_list[MAX_MOMENTUM_NUMBER][3];
  int filtered_sink_momentum_index[MAX_MOMENTUM_NUMBER];
  int index=0;
  for (int i=0;i<g_sink_momentum_number; ++i)
  {
    int momentum = g_sink_momentum_list[i][0]*g_sink_momentum_list[i][0]+g_sink_momentum_list[i][1]*g_sink_momentum_list[i][1]+g_sink_momentum_list[i][2]*g_sink_momentum_list[i][2];
    if (momentum< g_filtered_qsq){
      filtered_sink_momentum_list[index][0]=g_sink_momentum_list[i][0];
      filtered_sink_momentum_list[index][1]=g_sink_momentum_list[i][1];
      filtered_sink_momentum_list[index][2]=g_sink_momentum_list[i][2];
      filtered_sink_momentum_index[index]=i;
      index++;
    }
  }
  int filtered_sink_momentum_number= index;

  char *attribute_correlator;
  get_attribute_from_h5_file (&attribute_correlator, filename, "Correlator-info", io_proc ) ;

  char *attribute_ensemble_info;
  get_attribute_from_h5_file (&attribute_ensemble_info, filename, "Ensemble-info", io_proc ) ;
  char attribute_nmoms[100];
  snprintf(attribute_nmoms,100, "%d",filtered_sink_momentum_number);
  char attribute_qsq[100];
  snprintf(attribute_qsq,100,"%d",(int)g_filtered_qsq);


  snprintf ( filename, 400, "filtered_%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d.h5", filename_prefix, Nconf, filename_prefix2, g_nsample, Nsave, (int)g_filtered_qsq );
  if ( io_proc == 2 && g_verbose > 2 ) fprintf ( stdout, "# [loop_analyse] loop filename = %s\n", filename );

  exitstatus = loop_write_momentum_list_to_h5_file ( filtered_sink_momentum_list, filename, filtered_sink_momentum_number, io_proc );
  if ( exitstatus != 0 ) {
    fprintf ( stderr, "[] Error from loop_write_momentum_list_to_h5_file, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
    EXIT(1);

  }
  set_attribute_in_h5_file (attribute_correlator, filename, "Correlator-info", io_proc ) ;
  set_attribute_in_h5_file (attribute_ensemble_info, filename, "Ensemble-info", io_proc ) ;
  set_attribute_in_h5_file (attribute_nmoms, filename, "Nmoms", io_proc ) ;
  set_attribute_in_h5_file (attribute_qsq, filename, "Qsq", io_proc ) ;

  free(attribute_correlator);
  free(attribute_ensemble_info);





  if ( g_verbose > 2 && io_proc == 2 ) {
    for ( int imom = 0; imom < filtered_sink_momentum_number; imom++ ) {
      fprintf ( stdout, " %3d  ( %3d, %3d, %3d)\n", imom, filtered_sink_momentum_list[imom][0], filtered_sink_momentum_list[imom][1], filtered_sink_momentum_list[imom][2] );
    }
  } 

  /***************************************************************************
   * allocate memory for contractions
   ***************************************************************************/
  double **** loop = init_4level_dtable ( g_nsample, T, g_sink_momentum_number, 32 );
  if ( loop == NULL ) {
    fprintf(stderr, "[loop_analyse] Error from init_4level_dtable %s %d\n", __FILE__, __LINE__ );;
    EXIT(48);
  }

  /***************************************************************************
   * allocate memory for filtered_contractions
   ***************************************************************************/
  double **** loop_filtered = init_4level_dtable ( g_nsample, T, filtered_sink_momentum_number, 32 );
  if ( loop == NULL ) {
    fprintf(stderr, "[loop_analyse] Error from init_4level_dtable %s %d\n", __FILE__, __LINE__ );;
    EXIT(48);
  }

  /***************************************************************************
   * loop on different loop types
   ***************************************************************************/

  for ( int iloop_type=0; iloop_type < g_loop_number ; ++iloop_type){

  /***************************************************************************
   * loop on stochastic oet samples
   ***************************************************************************/
    for ( int isample = g_sourceid; isample < g_sourceid2; isample += g_sourceid_step )
    {

      int const Nstoch = isample * Nsave + 1;
      char loop_type[100];
      char loop_name[100];

      int inner_loop_length;

      switch( g_loop_type[iloop_type] ){
       case 0: snprintf ( loop_type, 100, "%s", "Scalar" );
               inner_loop_length=1;
               if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] Processing Scalar\n");
               break;
       case 1: snprintf ( loop_type, 100, "%s", "dOp" );
               inner_loop_length=1;
               if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] dOp\n");
               break;
       case 2: snprintf ( loop_type, 100, "%s", "Loops" );
               inner_loop_length=4;
               if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] Loops\n");
               break;
       case 3: snprintf ( loop_type, 100, "%s", "LpsDw" );
               inner_loop_length=4;
               if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] LpsDw\n");
               break;
       case 4: if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] LoopsCV\n");
               EXIT(1); 
               break;
       case 5: snprintf ( loop_type, 100, "%s", "LpsDwCv" );
               inner_loop_length=4;
               if (g_proc_id == 0) fprintf(stdout, "# [loop_analyse] LpsDwCv\n");
               break;
      }

      snprintf ( loop_name, 100, "%s", "loop" );

      for (int loop_direction=0; loop_direction < inner_loop_length ; ++loop_direction){

        char direction_part[100];
        snprintf(direction_part, 100, "dir_%02d", loop_direction);
     
        if (inner_loop_length > 1)
          snprintf ( data_tag, 100, "/conf_%.4d/Nstoch_%.4d/%s/%s/%s", Nconf, Nstoch, loop_type, direction_part, loop_name );
        else
          snprintf ( data_tag, 100, "/conf_%.4d/Nstoch_%.4d/%s/%s", Nconf, Nstoch, loop_type, loop_name );

        if ( io_proc == 2 && g_verbose > 2 ) fprintf( stdout, "# [loop_analyse] data_tag = %s\n", data_tag);

        snprintf ( filename, 400, "%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d.h5", filename_prefix, Nconf, filename_prefix2, g_nsample, Nsave, Qsq );

        exitstatus = loop_read_from_h5_file ( loop[isample], filename, data_tag, g_sink_momentum_number, 16, io_proc );
        if ( exitstatus != 0 ) {
          fprintf ( stderr, "[loop_analyse] Error from loop_read_from_h5_file, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
          EXIT(1);
        }

        if ( io_proc > 0 && g_verbose > 4 ) {
          /*****************************************************************
           * write in ASCII format
           *****************************************************************/
          for ( int iproc = 0; iproc < g_nproc_t; iproc++ ) {
            if ( g_tr_id == iproc ) {
              char output_filename[400];
              sprintf ( output_filename, "Nconf_%.4d.Nstoch_%.4d.%s.%s", Nconf, Nstoch, loop_type, loop_name );
              FILE * ofs = fopen ( output_filename, "w" );
              if ( ofs == NULL ) {
                fprintf ( stderr, "[loop_analyse] Error from fopen %s %d\n", __FILE__, __LINE__ );
                EXIT(1);
              }

              fprintf ( ofs, "# [loop_analyse] %s\n", data_tag );
              for ( int x0 = 0; x0 < T; x0++ ) {
                int const y0 = x0 + g_proc_coords[0] * T;

                for ( int imom = 0; imom < filtered_sink_momentum_number; imom++ ) {

                  for( int ic = 0; ic < 16; ic++ ) {

                    loop_filtered[isample][x0][imom][2*ic  ] = loop[isample][x0][filtered_sink_momentum_index[imom]][2*ic  ];
                    loop_filtered[isample][x0][imom][2*ic+1] = loop[isample][x0][filtered_sink_momentum_index[imom]][2*ic+1];

                    fprintf ( ofs, " %3d %4d   %3d% 3d% 3d   %d %d  %25.16e %25.16e\n", Nstoch, y0, 
                      g_sink_momentum_list[imom][0], g_sink_momentum_list[imom][1], g_sink_momentum_list[imom][2],
                      ic/4, ic%4, loop[isample][x0][imom][2*ic], loop[isample][x0][imom][2*ic+1] );

                  }  /* end of loop on components */
                }  /* end of loop on momenta */
              }  /* end of loop on time slices */

              fclose ( ofs );

            }  /* end of if g_tr_id == iproc */
#ifdef HAVE_MPI
          MPI_Barrier ( g_tr_comm );
#endif
          }  /* end of loop on procs in time direction */
        }  /* end of if io_proc > 0 and verbosity high level enough */
#if 0
#endif  /* of if 0 */

        /*****************************************************************/
        /*****************************************************************/

        if (inner_loop_length > 1)
          snprintf ( data_tag, 100, "/conf_%.4d/Nstoch_%.4d/%s/%s", Nconf, Nstoch, loop_type, direction_part );
        else
          snprintf ( data_tag, 100, "/conf_%.4d/Nstoch_%.4d/%s", Nconf, Nstoch, loop_type);

        if ( io_proc == 2 && g_verbose > 2 ) fprintf( stdout, "# [loop_analyse] data_tag = %s\n", data_tag);

        snprintf ( filename, 400, "filtered_%s.%.4d_%s_Ns%.4d_step%.4d_Qsq%d.h5", filename_prefix, Nconf, filename_prefix2, g_nsample, Nsave, (int)g_filtered_qsq );
        if ( io_proc == 2 && g_verbose > 2 ) fprintf ( stdout, "# [loop_analyse] loop filename = %s\n", filename );

        exitstatus = contract_loop_write_to_h5_file ( loop_filtered[isample], filename, data_tag, filtered_sink_momentum_number, 16, io_proc );
        if ( exitstatus != 0 ) {
          fprintf ( stderr, "[loop_analyse] Error from loop_write_to_h5_file, status was %d %s %d\n", exitstatus, __FILE__, __LINE__ );
          EXIT(1);
        }


        /*****************************************************************
         * set loop type
         *****************************************************************/
        /*
         * sprintf ( loop_type, "%s", "dOp" );
         * same for gen-oet
         */
      }  /* end of loop for the different loop directions */

    }  /* end of loop on oet samples */

  }  /* loop on different types */

  /***************************************************************************
   * decallocate fields
   ***************************************************************************/
  fini_4level_dtable ( &loop );
  fini_4level_dtable ( &loop_filtered ); 

  /***************************************************************************
   * free the allocated memory, finalize
   ***************************************************************************/

  free_geometry();

#ifdef HAVE_MPI
  mpi_fini_xchange_contraction();
  mpi_fini_datatypes();
  MPI_Finalize();
#endif

  if(g_cart_id==0) {
    g_the_time = time(NULL);
    fprintf(stdout, "# [loop_analyse] %s# [loop_analyse] end of run\n", ctime(&g_the_time));
    fprintf(stderr, "# [loop_analyse] %s# [loop_analyse] end of run\n", ctime(&g_the_time));
  }

  return(0);

}
