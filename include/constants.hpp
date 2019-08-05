#pragma once

namespace cvc {

  extern const char latDim_names[4];
  extern const char shift_dir_names[2];
  extern const mom_t zero_mom;
  
  namespace verbosity {
    /**
     * @brief verbosity level at which an extreme level of detail is passed to the log
     */
    constexpr int detailed_progress = 5;
    /**
     * @brief verbosity level at which all graph connection actions will be output in detail
     */
    constexpr int graph_connections = 4;
    /**
     * @brief verbosity level at which input parsers will relay what they read
     */
    constexpr int input_relay = 3;
    /**
     * @brief verbosity level at which [descend_and_resolve] reports in detail
     */
    constexpr int resolve = 3;
    /**
     * @brief verbosity level at which detailed memory info is output before and after
     *        contractions are done and written to disk
     */
    constexpr int memory_info = 2;
    constexpr int basic_progress = 1;
  } //verbosity
} //namespace(cvc)
