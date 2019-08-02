#include "exceptions.hpp"

#include <stdexcept>

void throw_error()
{
  throw cvc::runtime_error("bleh", "main");
}

int main(void){
  try
  {
    throw_error();
  } 
  
  catch (const std::exception& e)
  {
    std::cout << "Caught exception " << e.what() << std::endl;
  }

  return 0;
}
