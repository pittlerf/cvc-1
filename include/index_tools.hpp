#ifndef INDEX_TOOLS_HPP
#define INDEX_TOOLS_HPP

#include "cvc_global.h"
#include "enums.hpp"

#include <array>

namespace cvc {

  /**
   * @brief Convert a local to a global site index
   *
   * @param local_site_index
   *
   * @return global site index
   */
static inline
unsigned long long local_to_global_site_index(const unsigned int local_site_index)
{
  unsigned long long z_global, y_global, x_global, t_global;
  unsigned long long idx = local_site_index;
  z_global = g_proc_coords[DIM_Z]*LZ + (idx % LZ);
  idx /= LZ;
  y_global = g_proc_coords[DIM_Y]*LY + (idx % LY);
  idx /= LY;
  x_global = g_proc_coords[DIM_X]*LX + (idx % LX);
  idx /= LX;
  t_global = g_proc_coords[DIM_T]*T + (idx % T); 

  return z_global +
         LZ_global*y_global +
         LY_global*LZ_global*x_global +
         LX_global*LY_global*LZ_global*z_global;
}

static inline
int global_to_local_site_index(const unsigned long long global_site_index)
{
  unsigned int t,x,y,z;
  unsigned long long idx = global_site_index;
  z = static_cast<int>( (idx % LZ_global) % LZ ); 
  idx /= LZ_global;
  y = static_cast<int>( (idx % LY_global) % LY );
  idx /= LY_global;
  x = static_cast<int>( (idx % LX_global) % LX );
  idx /= LX_global;
  t = static_cast<int>( (idx % T_global) % T );

  return g_ipt[t][x][y][z];
}

static inline 
void global_site_get_proc_coords(unsigned long long global_site_index,
                                 std::array<int,4> & proc_coords)
{
  unsigned long long idx = global_site_index;
  proc_coords[DIM_Z] = static_cast<int>( (idx % LZ_global)/LZ );
  idx /= LZ_global;
  proc_coords[DIM_Y] = static_cast<int>( (idx % LY_global)/LY );
  idx /= LY_global;
  proc_coords[DIM_X] = static_cast<int>( (idx % LX_global)/LX );
  idx /= LX_global;
  proc_coords[DIM_T] = static_cast<int>( (idx % T_global)/T );
}

/**
 * @brief check if global site index is within the local lattice
 *
 * @param global_site_index
 *
 * @return
 */
static inline bool is_local(unsigned long long global_site_index)
{
  std::array<int,4> proc_coords;
  global_site_get_proc_coords(global_site_index, proc_coords);

  return ( (proc_coords[DIM_Z] == g_proc_coords[DIM_Z]) &&
           (proc_coords[DIM_Y] == g_proc_coords[DIM_Y]) &&
           (proc_coords[DIM_X] == g_proc_coords[DIM_X]) &&
           (proc_coords[DIM_T] == g_proc_coords[DIM_T]) 
         );
}

static inline int global_site_get_rank(unsigned long long global_site_index)
{
  std::array<int,4> proc_coords;
  global_site_get_proc_coords(global_site_index, proc_coords);
  int rank;
#ifdef HAVE_MPI
  MPI_Cart_rank(g_cart_grid, proc_coords.data(), &rank);
  return(rank);
#else
  return(0);
#endif
}

} // namespace(cvc)

#endif // INDEX_TOOLS_HPP
