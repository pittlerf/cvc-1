#include <global.h>
#include <enums.hpp>
#include <cvc_utils.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#include <hdf5.h>
#include <cstring>
#include <cstdio>
#include <string>
#include <iostream>
#include <string>
#include <list>

namespace cvc {
namespace h5 {

  static inline std::string path_list_to_key(const std::list<std::string> & path_list)
  {
    std::string key;
    for( auto const & subpath : path_list ){
      key += "/" + subpath;
    }
    return(key);
  }

  /**
   * @brief Simple writing of time-dependent data
   *
   *  Data should be ordered such that t is the slowest dimension and any
   *  reductions over the time slice communicator should already have been
   *  performed as the function below will be a no-op for any rank that
   *  is not at the origin of the XYZ cartesian grid of each time slice
   *  group.
   *
   * @param file
   * @param path_list
   * @param local_data
   */
  static inline void write_t_dataset(const std::string & filename, 
      const std::list<std::string> & path_list, 
      const std::vector<double> & local_data)
  {
    const int io_proc = get_io_proc();
  
    if(io_proc > 0){
#ifdef HAVE_MPI
      std::vector<double> buffer( local_data.size()*(T_global/T) );
      int existatus;
      CHECK_EXITSTATUS_NOT(
          existatus,
          MPI_SUCCESS,
          MPI_Gather(local_data.data(), local_data.size(), MPI_DOUBLE, 
                     buffer.data(), local_data.size(), MPI_DOUBLE, 0, g_tr_comm),
          "[cvc::h5::write_t_dataset] Failure in MPI_Gather\n",
          true,
          CVC_EXIT_MPI_FAILURE);
#else
      std::vector<double> & buffer = data;
#endif
    
      if(io_proc == 2){
        HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create); 
    
        std::string path;
        for( auto const & subpath : path_list ){
          path += "/" + subpath;
          // up until the last element, these are groups
          // the final element is instead the name of the dataset
          // and we will not create a group with this name
          if( !file.exist(path) && subpath != *(--path_list.end()) ){
            file.createGroup(path);
          }
        }
        // we just attempt to create the dataset. If it already exists, things will fail
        // and HighFive will give us a useful exception
        HighFive::DataSet dataset = file.createDataSet<double>(path, HighFive::DataSpace::From(buffer));
        dataset.write(buffer);
        file.flush();
      }
    }
  }


  #define H5UTILS_MAX_KEY_LENGTH 500
  
  /**
   * @brief check if a particular group / dataset exists
   * Traverses the h5 tree and checks if the requested group
   * hierarchy exists.
   *
   * @param loc_id HDF5 File or group id below which the hierarchy
   * should be traversed
   * @param key Group / Dataset key of the form "/grp1/grp2/grp3/dataset1"
   * @param fail_path First link path of searched hierarchy which
   *                  could not be found. Note pass as reference.
   *
   * @return true if path was found, false if it was not. Path of
   *         failure is returned vial fail_path argument.
   */
  static inline bool check_key_exists(hid_t loc_id,
                                      char const * const key,
                                      std::string & fail_path,
                                      bool full_path = true)
  {
    char key_copy[H5UTILS_MAX_KEY_LENGTH];
    if( strlen(key) < H5UTILS_MAX_KEY_LENGTH ){
      strcpy( key_copy, key );
    } else {
      char message[500];
      snprintf(message, 
               500, 
               "[h5_check_key_exists] length of key exceeds %d characters. %s line %d", 
               H5UTILS_MAX_KEY_LENGTH, __FILE__, __LINE__);
      EXIT_WITH_MSG(CVC_EXIT_H5_KEY_LENGTH, message);
    }
  
    htri_t status;
    std::string curr_path("");
  
    if( full_path ){
      curr_path = "/";
      status = H5Lexists(loc_id, curr_path.c_str(), H5P_DEFAULT);
      if( status <= 0 ){
        fail_path = "/";
        return false;
      }
    }
  
    const char grp_sep[] = "/";
    char * curr_grp_substr = strtok( key_copy, grp_sep );
    while( curr_grp_substr != NULL ){
      curr_path += std::string(curr_grp_substr);
      status = H5Lexists(loc_id, curr_path.c_str(), H5P_DEFAULT);
      if( status <= 0 ){
        fail_path = curr_path;
        return false;
      }
      curr_grp_substr = strtok( NULL, grp_sep );
      curr_path += std::string("/");
    }
  
    // traversal was successful
    fail_path = "";
    return true;
  }

} // namespace(h5)
} // namespace(cvc)
