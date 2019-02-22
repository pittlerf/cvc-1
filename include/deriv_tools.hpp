#pragma once

#include "types.h"
#include "enums.hpp"

#include <vector>
#include <iostream>

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
    std::cout << "[create_derivatives] Recursed from depth " << depth-1 << " to depth " << depth << 
      ", max_depth:" << max_depth << std::endl;
    for( auto const & curr_deriv_chain : curr_derivs ){
      std::cout << "[create_derivatives] curr_deriv_chain.size(): " << curr_deriv_chain.size() << std::endl;
      for( int dim = DIM_T; dim < DIM_NDIM; ++dim ){
        for( int dir = DIR_FWD; dir < DIR_NDIR; ++dir ){
          std::vector<deriv_t> deriv_chain = curr_deriv_chain;
          deriv_chain.push_back( deriv_t{dim, dir} );
          std::cout << "[create_derivatives] deriv_chain.size(): " << deriv_chain.size() << std::endl;
          derivs.push_back( deriv_chain );
        }
      }
    }
  }
  std::cout << "[create_derivatives] derivs.size(): " << derivs.size() << std::endl;
  
  if( depth+1 < max_depth ){
    return(create_derivatives(depth+1, max_depth, derivs));
  } else {
    return(derivs);
  }
}

} // namespace(cvc)

