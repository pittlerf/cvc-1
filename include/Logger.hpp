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
  const std::ostream out_target;
  std::ostream out;

public: 
  Logger() : 
    out_target(std::cout.rdbuf()),
    out(std::cout.rdbuf())
  {
    set_verbosity_threshold(0);
    set_output_proc_id(0);
  }

  Logger(const int proc_id_in, const int verbosity_threshold_in) :
    out_target(std::cout.rdbuf()),
    out(std::cout.rdbuf())
  {
    set_verbosity_threshold(verbosity_threshold_in);
    set_output_proc_id(proc_id_in);
  }

  Logger(const int proc_id_in, const int verbosity_threshold_in, std::ostream & out_in) :
    out_target(out_in.rdbuf()),
    out(out_in.rdbuf())
  {
    set_verbosity_threshold(verbosity_threshold_in);
    set_output_proc_id(proc_id_in);
  }

  void set_verbosity_threshold(const int level)
  {
    verbosity_threshold = level;
  }
  
  void set_output_proc_id(const int id)
  {
    output_proc_id = id;
    // only the target process will actually direct its output
    // to the target stream while all others will use the
    // nulbuf
    if( output_proc_id == g_proc_id ){
      out.rdbuf(out_target.rdbuf());
    } else {
      out.rdbuf(&nulbuf);
    }
  }
 
  template<typename T>
  std::ostream & operator<<(T val)
  {
    if( verbosity_threshold <= g_verbose ){
      return (out << val);
    }
  }

  // support for io maniupulators and std::endl
  std::ostream & operator<<( std::ios_base & (*func)(std::ios_base &) )
  {
    if( verbosity_threshold <= g_verbose ){
      return (out << func);
    }
  }

  std::ostream & operator<<( std::ostream& (*func)(std::ostream &) )
  {
    if( verbosity_threshold <= g_verbose ){
      return (out << func);
    }
  }

};

} //namespace(cvc)

