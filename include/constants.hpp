#pragma once

namespace cvc {

  extern const char latDim_names[4];
  extern const char shift_dir_names[2];
  extern const mom_t zero_mom;
  
  namespace verbosity {
    /**
     * @brief verbosity level at which input parsers will relay what they read
     */
    constexpr int input_relay = 3;
    constexpr int graph_connections = 4;
    constexpr int resolve = 3;
    constexpr int memory_info = 2;
    constexpr int basic_progress = 1;
    constexpr int detailed_progress = 5;
  } //verbosity
} //namespace(cvc)
