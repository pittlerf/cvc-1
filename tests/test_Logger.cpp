#define MAIN_PROGRAM
#include "cvc_global.h"
#undef MAIN_PROGRAM

#include "Logger.hpp"
#include "Core.hpp"

#include <iostream>

int main(int argc, char ** argv)
{
  cvc::Core core(argc,argv);
  cvc::Logger logger(cvc::g_proc_id,0,std::cout);

  logger << std::endl;
  logger << "test " << cvc::g_proc_id << std::endl;
  logger << std::endl;
 
  return 0;
}
