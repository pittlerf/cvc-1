#define MAIN_PROGRAM
#include "cvc_global.h"
#undef MAIN_PROGRAM

#include "h5utils.hpp"
#include "Core.hpp"

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Group.hpp>

#include <list>
#include <string>

int main(int argc, char ** argv)
{
  cvc::Core core(argc,argv);

  std::vector<double> data(cvc::T, cvc::g_proc_coords[0]);

  std::list<std::string> path_list;
  path_list.push_back("group2");
  path_list.push_back("subgroup2");
  path_list.push_back("subsubgroup");
  path_list.push_back("dataset_two");

  cvc::h5::write_t_dataset("test.h5", path_list, data); 

  return 0;
}
