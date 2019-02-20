#include "h5utils.hpp"

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Group.hpp>

#include <list>
#include <string>

int main(int argc, char ** argv)
{
  HighFive::File file("test.h5", 
    HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

  std::vector<double> data(50, 3.1415);
  HighFive::Group group = file.createGroup("/group");
  HighFive::Group subgroup = file.createGroup("/group/subgroup");
  HighFive::DataSet dataset = file.createDataSet<double>("/group/subgroup/dataset_one", HighFive::DataSpace::From(data));
  dataset.write(data);
  file.flush();

  std::list<std::string> path_list;
  path_list.push_back("group2");
  path_list.push_back("subgroup2");
  path_list.push_back("subsubgroup");

  cvc::h5::write_dataset(file, path_list, "dataset_two", data); 

  return 0;
}
