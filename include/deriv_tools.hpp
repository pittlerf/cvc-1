#pragma once

#include "types.h"
#include "enums.hpp"

#include <vector>

namespace cvc {

static inline
std::vector< std::vector<deriv_t> >
create_derivative_chains(
  const int depth,
  const int max_depth,
  const std::vector< std::vector<deriv_t> > curr_deriv_chains)
{
  std::vector< std::vector<deriv_t> > deriv_chains;
  if( depth == 0 && max_depth > 0 ){
    for( int dim = DIM_T; dim < DIM_NDIM; ++dim ){
      for( int dir = DIR_FWD; dir < DIR_NDIR; ++dir ){
        std::vector<deriv_t> deriv_chain{ deriv_t{dim, dir} };
        deriv_chains.push_back( deriv_chain );
      }
    }
  } else {
    for( auto const & curr_deriv_chain : curr_deriv_chains ){
      for( int dim = DIM_T; dim < DIM_NDIM; ++dim ){
        for( int dir = DIR_FWD; dir < DIR_NDIR; ++dir ){
          std::vector<deriv_t> deriv_chain = curr_deriv_chain;
          deriv_chain.push_back( deriv_t{dim, dir} );
          deriv_chains.push_back( deriv_chain );
        }
      }
    }
  }
  
  if( depth+1 < max_depth ){
    return(create_derivative_chains(depth+1, max_depth, deriv_chains));
  } else {
    return(deriv_chains);
  }
}

} // namespace(cvc)

