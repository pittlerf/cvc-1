#pragma once

#include "types.h"

namespace cvc {

static inline void describe_mom(const mom_t & mom, int & pnonzero, int & psq, int & pmax)
{
  psq = mom.x*mom.x + mom.y*mom.y + mom.z*mom.z;
  pnonzero = static_cast<int>( mom.x != 0 ) + static_cast<int>( mom.y != 0 ) + 
             static_cast<int>( mom.z != 0 );
  pmax = 0;
  for( const auto & p : {mom.x, mom.y, mom.z} ){
    pmax = abs(p) > pmax ? abs(p) : pmax;
  }
}

/**
 * @brief Comparison functor for momentum triplets, suitable for std::sort
 *
 * we sort momenta by psq first.
 * For equal psq, we sort by maximum component: (2,2,1) < (0,0,3)
 * and then by the number of non-zeros (not sure if this last condition would ever occur)
 */
typedef struct momentum_compare_t {
  bool operator()( const mom_t & a, const mom_t & b ){
    int asq, anonzero, amax;
    int bsq, bnonzero, bmax;
    describe_mom(a, anonzero, asq, amax);
    describe_mom(b, bnonzero, bsq, bmax);

    if( asq == bsq ){
      if( amax == bmax ){
        return anonzero < bnonzero;
      } else {
        return amax < bmax;
      }
    } else {
      return asq < bsq;
    }
  }
} momentum_compare_t;

} //namespace(cvc)

