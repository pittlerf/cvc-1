#pragma once

#include "cvc_global.h"
#include "enums.hpp"
#include "cvc_utils.h"
#include "debug_printf.hpp"
#include "constants.hpp"

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

  static inline std::string recursive_path_create(
      HighFive::File & file,
      const std::list<std::string> & path_list,
      const bool is_dataset_path = true)
  {
    std::string path;
    for( auto const & subpath : path_list ){
      path += "/" + subpath;
      // up until the last element, these are groups
      // the final element is instead the name of the dataset
      // and we will not create a group with this name
      auto end_elem = path_list.end();
      if( is_dataset_path ){
        end_elem = (--path_list.end());
      }
      if( !file.exist(path) && subpath != *(end_elem) ){
        debug_printf(0, verbosity::detailed_progress,
            "# [cvc::h5::recursive_path_create] Creating H5 path %s\n", path.c_str());

        file.createGroup(path);
      }
    }
    return path;
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
      int exitstatus;
      CHECK_EXITSTATUS_NOT(
          exitstatus,
          MPI_SUCCESS,
          MPI_Gather(local_data.data(), local_data.size(), MPI_DOUBLE, 
                     buffer.data(), local_data.size(), MPI_DOUBLE, 0, g_tr_comm),
          "[cvc::h5::write_t_dataset] Failure in MPI_Gather\n",
          true,
          CVC_EXIT_MPI_FAILURE);
#else
      const std::vector<double> & buffer = data;
#endif
    
      if(io_proc == 2){
        HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create); 
    
        // create the path hierarchy all the way up until the dataset name
        std::string path = recursive_path_create(file, path_list, true);

        // we check if the dataset path exists and if it does not, we create it
        // and we attempt to write the dataset, if it already exsists it should
        // just be overwritten, as long as the buffer has the same size
        if( !file.exist(path_list_to_key(path_list)) ){
          HighFive::DataSet dataset = file.createDataSet<double>(path, HighFive::DataSpace::From(buffer));
          dataset.write(buffer);
        } else {
          HighFive::DataSet dataset = file.getDataSet(path);
          dataset.write(buffer);
        }
        file.flush();
      }
    }
  }

  static inline void write_correlators(const std::string & filename,
      std::map< std::string, H5Correlator > & corrs)
  {
    const int io_proc = get_io_proc();

    if(io_proc > 0){
      // we use a shared ptr to the file such that we only have to open it once on io_proc==2, 
      // rather than opening it again and again in the loop below
      std::shared_ptr<HighFive::File> file_ptr;
      std::string path;
      std::vector<double> buffer;

      if(io_proc == 2){
        file_ptr.reset( new HighFive::File(filename, HighFive::File::ReadWrite | HighFive::File::Create ) );
      }
      for( auto & elem : corrs ){
        int elem_size = elem.second.storage.size();
        
        if(io_proc == 2){
          buffer.resize( elem_size * (T_global/T) );
          // we create the path to the dataset up to the dataset name itself
          path = recursive_path_create( *(file_ptr), elem.second.path_list, true);
          if( g_verbose >= verbosity::detailed_progress ){
            std::cout << "# [cvc::h5::write_correlators] Dataset path : " << path << std::endl;
          }
        } else {
          // on all other processes, make sure that buffer.data() points to a valid
          // address, even if it is not used
          buffer.resize(1);
        }

        std::vector<double> & local_data = elem.second.storage;
#ifdef HAVE_MPI
        int exitstatus;
        CHECK_EXITSTATUS_NOT(
            exitstatus,
            MPI_SUCCESS,
            MPI_Gather(local_data.data(), elem_size, MPI_DOUBLE, 
                       buffer.data(), elem_size, MPI_DOUBLE, 0, g_tr_comm),
            "[cvc::h5::write_correlators] Failure in MPI_Gather\n",
            true,
            CVC_EXIT_MPI_FAILURE);
#else
        buffer = data;
#endif
        if(io_proc == 2){
          // if the dataset does not exist it, we create it
          // otherwise we open it and attempt to write it
          if( !file_ptr->exist(path) ){
            HighFive::DataSet dataset = file_ptr->createDataSet<double>(path, HighFive::DataSpace::From(buffer));
            dataset.write(buffer);
          } else {
            HighFive::DataSet dataset = file_ptr->getDataSet(path);
            dataset.write(buffer);
          }
        }
      } // loop over corrs
      if( io_proc == 2 ){
        file_ptr->flush();
      }
    } // if(io_proc > 0)
  }

  static inline void write_loops(const std::string & filename,
      const std::vector<double> & local_loops,
      const std::vector< std::list<std::string> > & path_lists,
      const unsigned int n_gamma, const unsigned int n_mom)
  {
    const int io_proc = get_io_proc();

    // only the accumulator processes in T perform any operations here
    if( io_proc > 0 ){
      std::vector<double> buffer;
      std::vector<double> reorder_buffer;
      const size_t n_elem = 2*T_global*path_lists.size();
      
      // on the true origin process, we need to accumulate the loops and reshuffle them
      // so we need some storage
      if(io_proc == 2){
        buffer.resize( n_elem );
        reorder_buffer.resize( n_elem );
      } else {
        // on all other processes, make sure that buffer.data() points to a valid
        // address, even if it is not used
        buffer.resize(1);
      }

      // gather the loops from all accumulator processes on the I/O process
#ifdef HAVE_MPI
      int exitstatus;
      CHECK_EXITSTATUS_NOT(
          exitstatus,
          MPI_SUCCESS,
          MPI_Gather(local_loops.data(), n_elem*T/T_global, MPI_DOUBLE, 
                     buffer.data(), n_elem*T/T_global, MPI_DOUBLE, 0, g_tr_comm),
          "[cvc::h5::write_loops] Failure in MPI_Gather\n",
          true,
          CVC_EXIT_MPI_FAILURE);
#else
      buffer = local_loops;
#endif

      if(io_proc == 2){
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
        {
          // now that we have the data from all T processes, we can reorder
          // appropriately for writeout
          size_t in_idx;
          size_t out_idx; 
#ifdef HAVE_OPENMP
#pragma omp for collapse(3)
#endif
          for(unsigned int t = 0; t < T_global; ++t){
            for(unsigned int i_mom = 0; i_mom < n_mom; ++i_mom){
              for(unsigned int i_gamma = 0; i_gamma < n_gamma; ++i_gamma){
                // in slowest to fastest: t, mom, gamma, complex
                in_idx = 2*(i_gamma + 16*(i_mom + n_mom*t));
                // out slowest to fastest: gamma, mom, T, complex
                out_idx = 2*(t + T_global*(i_mom + n_mom*i_gamma));
                reorder_buffer[out_idx  ] = buffer[in_idx];
                reorder_buffer[out_idx+1] = buffer[in_idx+1];
              }
            }
          }
        } // OpenMP parallel closing brace

        std::string path;
        HighFive::File loop_file(filename, HighFive::File::ReadWrite | HighFive::File::Create );
        // each loop that we are about to write is 2*T_global in size
        for( size_t i_elem = 0; i_elem < path_lists.size(); ++i_elem ){
          path = recursive_path_create(loop_file, path_lists[i_elem]);
          if( g_verbose >= verbosity::detailed_progress ){
            std::cout << "# [cvc::h5::write_loops] Dataset path : " << path << std::endl;
          }
          // currently we make a copy but there must be a better way for doing this
          // using HigFive, couldn't figure it out though...
          std::vector<double> temp(reorder_buffer.begin() + i_elem*2*T_global, 
                                   reorder_buffer.begin() + (i_elem+1)*2*T_global);
        
          // if the data set does not exist, we create it, otherwise we simply attempt to write
          // it!
          if( !loop_file.exist(path) ){
            HighFive::DataSet dataset = loop_file.createDataSet<double>(path, HighFive::DataSpace::From(temp));
            dataset.write(temp);
          } else {
            HighFive::DataSet dataset = loop_file.getDataSet(path);
            dataset.write(temp);
          }
        }
        loop_file.flush();
      } // if(io_proc == 2)
    } // if(io_proc > 0)
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
   * @param full_path Pass 'true' if the key refers to a complete
   *                  path with a leading forward slash.
   *
   * @return true if path was found, false if it was not. Path of
   *         failure is returned vial fail_path argument.
   */
  static inline bool check_key_exists(hid_t loc_id,
                                      char const * const key,
                                      std::string & fail_path,
                                      bool full_path = true)
  {
    if( full_path ){
      char key_copy[H5UTILS_MAX_KEY_LENGTH];
      if( strlen(key) < H5UTILS_MAX_KEY_LENGTH ){
        strcpy( key_copy, key );
      } else {
        char message[500];
        snprintf(message, 
                 500, 
                 "[cvc::h5::check_key_exists] length of key exceeds %d characters. %s line %d", 
                 H5UTILS_MAX_KEY_LENGTH, __FILE__, __LINE__);
        EXIT_WITH_MSG(CVC_EXIT_H5_KEY_LENGTH, message);
      }
  
      htri_t status;
      std::string curr_path("");
  
      curr_path = "/";
      status = H5Lexists(loc_id, curr_path.c_str(), H5P_DEFAULT);
      if( status <= 0 ){
        fail_path = "/";
        return false;
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
    } else {
      if( H5Lexists(loc_id, key, H5P_DEFAULT) <= 0 ){
        return false;
      } else {
        fail_path = "";
        return true;
      }
    } // if(full_path)
  }

} // namespace(h5)
} // namespace(cvc)
