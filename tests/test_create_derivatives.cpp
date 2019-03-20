#include "deriv_tools.hpp"

#include "types.h"

#include <vector>
#include <iostream>

using namespace cvc;

int main(int argc, char ** argv)
{
  std::vector< std::vector<deriv_t> > derivs;
  derivs = create_derivatives(0, 3, derivs);

  std::cout << "[test_create_derivatives] derivs.size(): " << derivs.size() << std::endl;
  for(auto const & deriv_chain : derivs){
    for( auto const & deriv : deriv_chain ){
      std::cout << "(" << deriv.dim << "," << deriv.dir << ") ";
    }
    std::cout << std::endl;
  }
  return 0;
}
