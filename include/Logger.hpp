#pragma once

#include "global.h"
#include <iostream>
#include <string>
#include <ios>

namespace cvc {


  /**
   * @brief A streambuf which does not produce any output
   */
class NullBuffer : public std::basic_streambuf<char>
{
public:
  int overflow(int c) { return c; }
};

class Logger {
private:
  int output_proc_id;
  int verbosity_threshold;
  NullBuffer nulbuf;
  std::ostream nul_stream;
  std::ostream out_stream;

public: 
  Logger() : 
    out_stream(std::cout.rdbuf()),
    nul_stream(&nulbuf),
    verbosity_threshold(0),
    output_proc_id(0) {}

  Logger(const int proc_id_in, const int verbosity_threshold_in) :
    out_stream(std::cout.rdbuf()),
    nul_stream(&nulbuf),
    verbosity_threshold(verbosity_threshold_in),
    output_proc_id(proc_id_in) {}

  Logger(const int proc_id_in, const int verbosity_threshold_in, std::ostream & out_in) :
    out_stream(out_in.rdbuf()),
    nul_stream(&nulbuf),
    verbosity_threshold(verbosity_threshold_in),
    output_proc_id(proc_id_in) {}

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
      return (out_stream << val);
    } else {
      return (nul_stream << val);
    }
  }

  // support for io maniupulators and std::endl
  std::ostream & operator<<( std::ios_base & (*func)(std::ios_base &) )
  {
    if( output_proc_id == g_proc_id && verbosity_threshold <= g_verbose ){
      return (out_stream << func);
    } else {
      return (nul_stream << func);
    }
  }

  std::ostream & operator<<( std::ostream& (*func)(std::ostream &) )
  {
    if( output_proc_id == g_proc_id && verbosity_threshold <= g_verbose ){
      return (out_stream << func);
    } else {
      return (nul_stream << func);
    }
  }

};

} //namespace(cvc)

