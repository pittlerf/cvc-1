#pragma once

#include "types.h"
#include "enums.hpp"

#include <vector>

namespace cvc {

static inline
std::vector< std::vector<deriv_t> >
create_derivatives(
  const int depth,
  const int max_depth,
  const std::vector< std::vector<deriv_t> > curr_derivs)
{
  std::vector< std::vector<deriv_t> > derivs;
  if( depth == 0 && max_depth > 0 ){
    for( int dim = DIM_T; dim < DIM_NDIM; ++dim ){
      for( int dir = DIR_FWD; dir < DIR_NDIR; ++dir ){
        std::vector<deriv_t> deriv_chain{ deriv_t{dim, dir} };
        derivs.push_back( deriv_chain );
      }
    }
  } else {
    for( auto const & curr_deriv_chain : curr_derivs ){
      for( int dim = DIM_T; dim < DIM_NDIM; ++dim ){
        for( int dir = DIR_FWD; dir < DIR_NDIR; ++dir ){
          std::vector<deriv_t> deriv_chain = curr_deriv_chain;
          deriv_chain.push_back( deriv_t{dim, dir} );
          derivs.push_back( deriv_chain );
        }
      }
    }
  }
  
  if( depth+1 < max_depth ){
    return(create_derivatives(depth+1, max_depth, derivs));
  } else {
    return(derivs);
  }
}

} // namespace(cvc)

