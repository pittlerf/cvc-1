#!/bin/bash

TYPE="double _Complex"
TAG="z"

TTAG=$( echo $TAG | tr '[:lower:]' '[:upper:]')

FILE=table_init_${TAG}.h

cat << EOF > $FILE
#ifndef _TABLE_INIT_${TTAG}_H
#define _TABLE_INIT_${TTAG}_H

/****************************************************
 * table_init_${TAG}.h
 *
 * PURPOSE:
 * DONE:
 * TODO:
 ****************************************************/

namespace cvc {
EOF


cat << EOF >> $FILE

static inline $TYPE * init_1level_${TAG}table ( int const N0 ) {
  return( ( $TYPE *) calloc ( N0 , sizeof( $TYPE ) ) );
}  // end of init_1level_${TAG}table

/************************************************************************************/
/************************************************************************************/

static inline void fini_1level_${TAG}table ( ${TYPE} **s  ) {
  if ( *s != NULL ) free ( *s );
  // fprintf ( stdout, "# [fini_1level_${TAG}table] active\n");
  *s = NULL;
}  // end of fini_1level_${TAG}table

/************************************************************************************/
/************************************************************************************/

EOF



for LEVEL in $(seq 2 8 ); do

  PTR2=""
  for ((k=1; k < $LEVEL; k++ )) do
    PTR2="${PTR2}*"
  done
  PTR="${PTR2}*"
  PTR3="${PTR}*"

  printf "static inline %s %s init_%dlevel_%stable (" "$TYPE" "$PTR"  $LEVEL "$TAG"
for ((k=1; k < $LEVEL; k++ )) do
  printf "int const N%d, " $(($k - 1 ))
done
printf "int const N%d ) {\n" $(($LEVEL - 1 ))

printf "  %s %s s__ = NULL;\n" "${TYPE}" "$PTR2"
printf "  s__ = init_%dlevel_%stable ( N0*N1" $(( $LEVEL - 1 )) "$TAG"

for ((k=2; k < $LEVEL; k++ )) do
  printf ", N%d" $k
done
printf ");\n"

cat << EOF
  if ( s__ == NULL ) return( NULL );

  ${TYPE} $PTR s_ = ( ${TYPE} $PTR) malloc( N0 * sizeof( ${TYPE} $PTR2) );
  if ( s_ == NULL ) return ( NULL );

  for ( int i = 0; i < N0; i++ ) s_[i] = s__ + i * N1;
  return( s_ );
}  // end of init_${LEVEL}level_${TAG}table

/************************************************************************************/
/************************************************************************************/


static inline void fini_${LEVEL}level_${TAG}table ( ${TYPE} ${PTR3} s  ) {
  if ( *s != NULL ) {
    // fprintf ( stdout, "# [fini_${LEVEL}level_${TAG}table] active\n");
    fini_$(( $LEVEL - 1))level_${TAG}table ( *s );
    free ( *s );
    *s = NULL;
  }
}  // end of fini_${LEVEL}level_${TAG}table

/************************************************************************************/
/************************************************************************************/


EOF

done >> $FILE

cat << EOF >> $FILE
}  /* end of namespace cvc */

#endif
EOF

exit 0