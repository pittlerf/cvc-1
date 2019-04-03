#pragma once

#include "types.h"
#include "enums.hpp"

#include <vector>

namespace cvc {

  /**
   * @brief create a vector of vectors to represent chains of covariant derivatives
   *
   * A covariant derivative is applied in a given direction (forward or backward)
   * and dimension (t,x,y,z <-> 0,1,2,3).
   *
   * A chain of such derivatives can also be applied in order, for example:
   * 
   *    D_0^f D_1^b D_2^f
   *
   * If taken as operator notation with the subscript denoting the space-time dimension
   * and the superscript denoting the direction (f)orward, (b)ackward, this is read
   * right to left.
   *
   * First, D_2^f is applied, then D_1^b and finally D_0^f.
   *
   * This function recursively creates a vector of vectors of pairs, each element 
   * of type 'deriv_t'y  with elements .dim and .dir for the dimension and direction,
   * respectively.
   *
   * To use it, one would create an empty std::vector< std::vector< cvc::deriv_t > >
   * externally and pass it as the third argument to this function. One would
   * specify 'depth = 0' to indicate that this is the first level and
   * 'max_depth' specifies the number of derivative applications.
   *
   * The function will then create all combinations of (0,1,2,3) and (0,1) for
   * the dimensions and directions, respectively, until the maximum depth is
   * reached.
   *
   * For example (see tests/test_create_derivative_chains.cpp)
   *
   * std::vector< std::vector<deriv_t> > derivs;
   * derivs = create_derivative_chains(0, 3, derivs);
   *
   * creates a vector of vectors with 512 elements of 3 elements each:
   *
   * (0,0) (0,0) (0,0) 
   * (0,0) (0,0) (0,1) 
   * (0,0) (0,0) (1,0) 
   * (0,0) (0,0) (1,1) 
   * (0,0) (0,0) (2,0) 
   * (0,0) (0,0) (2,1)
   * [...]
   * (2,0) (3,1) (2,1) 
   * (2,0) (3,1) (3,0) 
   * (2,0) (3,1) (3,1) 
   * (2,1) (0,0) (0,0) 
   * (2,1) (0,0) (0,1) 
   * (2,1) (0,0) (1,0)  
   * [...]
   * (3,1) (3,1) (1,1) 
   * (3,1) (3,1) (2,0) 
   * (3,1) (3,1) (2,1) 
   * (3,1) (3,1) (3,0) 
   * (3,1) (3,1) (3,1) 
   *
   * @param depth Current recursion level.
   * @param max_depth Maximal recusrion level.
   * @param curr_deriv_chains Current state of derivative chains vector.
   *
   * @return vector of vectors with derivative chains of the requested depth 
   */
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

