#ifndef _TABLE_INIT_Z_H
#define _TABLE_INIT_Z_H

/****************************************************
 * table_init_z.h
 *
 * PURPOSE:
 * DONE:
 * TODO:
 ****************************************************/

namespace cvc {

static inline double _Complex * init_1level_ztable ( int const N0 ) {
  return( ( double _Complex *) calloc ( N0 , sizeof( double _Complex ) ) );
}  // end of init_1level_ztable

/************************************************************************************/
/************************************************************************************/

static inline void fini_1level_ztable ( double _Complex **s  ) {
  if ( *s != NULL ) free ( *s );
  // fprintf ( stdout, "# [fini_1level_ztable] active\n");
  *s = NULL;
}  // end of fini_1level_ztable

/************************************************************************************/
/************************************************************************************/

static inline double _Complex ** init_2level_ztable (int const N0, int const N1 ) {
  double _Complex * s__ = NULL;
  s__ = init_1level_ztable ( N0*N1);
  if ( s__ == NULL ) return( NULL );

  double _Complex ** s_ = ( double _Complex **) malloc( N0 * sizeof( double _Complex *) );
  if ( s_ == NULL ) return ( NULL );

  for ( int i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_2level_ztable

/************************************************************************************/
/************************************************************************************/


static inline void fini_2level_ztable ( double _Complex *** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_2level_ztable] active\n");
    fini_1level_ztable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_2level_ztable

/************************************************************************************/
/************************************************************************************/


static inline double _Complex *** init_3level_ztable (int const N0, int const N1, int const N2 ) {
  double _Complex ** s__ = NULL;
  s__ = init_2level_ztable ( N0*N1, N2);
  if ( s__ == NULL ) return( NULL );

  double _Complex *** s_ = ( double _Complex ***) malloc( N0 * sizeof( double _Complex **) );
  if ( s_ == NULL ) return ( NULL );

  for ( int i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_3level_ztable

/************************************************************************************/
/************************************************************************************/


static inline void fini_3level_ztable ( double _Complex **** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_3level_ztable] active\n");
    fini_2level_ztable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_3level_ztable

/************************************************************************************/
/************************************************************************************/


static inline double _Complex **** init_4level_ztable (int const N0, int const N1, int const N2, int const N3 ) {
  double _Complex *** s__ = NULL;
  s__ = init_3level_ztable ( N0*N1, N2, N3);
  if ( s__ == NULL ) return( NULL );

  double _Complex **** s_ = ( double _Complex ****) malloc( N0 * sizeof( double _Complex ***) );
  if ( s_ == NULL ) return ( NULL );

  for ( int i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_4level_ztable

/************************************************************************************/
/************************************************************************************/


static inline void fini_4level_ztable ( double _Complex ***** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_4level_ztable] active\n");
    fini_3level_ztable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_4level_ztable

/************************************************************************************/
/************************************************************************************/


static inline double _Complex ***** init_5level_ztable (int const N0, int const N1, int const N2, int const N3, int const N4 ) {
  double _Complex **** s__ = NULL;
  s__ = init_4level_ztable ( N0*N1, N2, N3, N4);
  if ( s__ == NULL ) return( NULL );

  double _Complex ***** s_ = ( double _Complex *****) malloc( N0 * sizeof( double _Complex ****) );
  if ( s_ == NULL ) return ( NULL );

  for ( int i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_5level_ztable

/************************************************************************************/
/************************************************************************************/


static inline void fini_5level_ztable ( double _Complex ****** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_5level_ztable] active\n");
    fini_4level_ztable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_5level_ztable

/************************************************************************************/
/************************************************************************************/


static inline double _Complex ****** init_6level_ztable (int const N0, int const N1, int const N2, int const N3, int const N4, int const N5 ) {
  double _Complex ***** s__ = NULL;
  s__ = init_5level_ztable ( N0*N1, N2, N3, N4, N5);
  if ( s__ == NULL ) return( NULL );

  double _Complex ****** s_ = ( double _Complex ******) malloc( N0 * sizeof( double _Complex *****) );
  if ( s_ == NULL ) return ( NULL );

  for ( int i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_6level_ztable

/************************************************************************************/
/************************************************************************************/


static inline void fini_6level_ztable ( double _Complex ******* s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_6level_ztable] active\n");
    fini_5level_ztable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_6level_ztable

/************************************************************************************/
/************************************************************************************/


static inline double _Complex ******* init_7level_ztable (int const N0, int const N1, int const N2, int const N3, int const N4, int const N5, int const N6 ) {
  double _Complex ****** s__ = NULL;
  s__ = init_6level_ztable ( N0*N1, N2, N3, N4, N5, N6);
  if ( s__ == NULL ) return( NULL );

  double _Complex ******* s_ = ( double _Complex *******) malloc( N0 * sizeof( double _Complex ******) );
  if ( s_ == NULL ) return ( NULL );

  for ( int i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_7level_ztable

/************************************************************************************/
/************************************************************************************/


static inline void fini_7level_ztable ( double _Complex ******** s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_7level_ztable] active\n");
    fini_6level_ztable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_7level_ztable

/************************************************************************************/
/************************************************************************************/


static inline double _Complex ******** init_8level_ztable (int const N0, int const N1, int const N2, int const N3, int const N4, int const N5, int const N6, int const N7 ) {
  double _Complex ******* s__ = NULL;
  s__ = init_7level_ztable ( N0*N1, N2, N3, N4, N5, N6, N7);
  if ( s__ == NULL ) return( NULL );

  double _Complex ******** s_ = ( double _Complex ********) malloc( N0 * sizeof( double _Complex *******) );
  if ( s_ == NULL ) return ( NULL );

  for ( int i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_8level_ztable

/************************************************************************************/
/************************************************************************************/


static inline void fini_8level_ztable ( double _Complex ********* s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_8level_ztable] active\n");
    fini_7level_ztable ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_8level_ztable

/************************************************************************************/
/************************************************************************************/


}  /* end of namespace cvc */

#endif