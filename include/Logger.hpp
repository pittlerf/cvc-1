#pragma once

#include "global.h"
#include <iostream>
#include <string>
#include <ios>

namespace cvc {

class Logger {
private:
  int output_proc_id;
  int verbosity_threshold;
  std::ostream out;

public: 
  Logger() : output_proc_id(0), verbosity_threshold(0), out(std::cout.rdbuf()) {}

  Logger(const int proc_id_in, const int verbosity_threshold_in, std::ostream & out_in) :
    output_proc_id(proc_id_in), verbosity_threshold(verbosity_threshold_in), out(out_in.rdbuf()) {}

  void set_verbosity_threshold(const int level)
  {
    verbosity_threshold = level;
  }
  
  void set_output_proc_id(const int id)
  {
    output_proc_id = id;
  }
 
  template<typename T>
  std::ostream & operator<<(T val)
  { 
    if( output_proc_id == g_proc_id && verbosity_threshold <= g_verbose ){
      return (out << val);
    } else {
      return (out);
    }
  }

  // support for io maniupulators and std::endl
  std::ostream & operator<<( std::ios_base & (*func)(std::ios_base &) )
  {
    if( output_proc_id == g_proc_id && verbosity_threshold <= g_verbose ){
      return (out << func);
    } else {
      return (out);
    }
  }

  std::ostream & operator<<( std::ostream& (*func)(std::ostream &) )
  {
    if( output_proc_id == g_proc_id && verbosity_threshold <= g_verbose ){
      return (out << func);
    } else {
      return (out);
    }
  }

};

} //namespace(cvc)

