#include <stdexcept>
#include <string>
#include <iostream>

namespace cvc {

class runtime_error : public std::runtime_error {
  public:
    runtime_error( const char * what_arg, const char * function_name ) :
      runtime_error( std::string(what_arg), std::string(function_name) ) {}

    runtime_error( const char * what_arg, const std::string& function_name ) :
      runtime_error( std::string(what_arg), function_name ) {}

    runtime_error( const std::string& what_arg, const char * function_name ) :
      runtime_error( what_arg, std::string(function_name) ) {} 

    runtime_error( const std::string& what_arg, const std::string& function_name ) :
      std::runtime_error(what_arg)
    {
      std::cout << "[cvc::runtime_error]: exception in " << 
        function_name << " : " <<
        what_arg << std::endl; 
    }
};

class invalid_argument : public std::invalid_argument {
  public:
    invalid_argument( const char * what_arg, const char * function_name ) :
      invalid_argument( std::string(what_arg), std::string(function_name) ) {}

    invalid_argument( const char * what_arg, const std::string& function_name ) :
      invalid_argument( std::string(what_arg), function_name ) {}

    invalid_argument( const std::string& what_arg, const char * function_name ) :
      invalid_argument( what_arg, std::string(function_name) ) {} 

    invalid_argument( const std::string& what_arg, const std::string& function_name ) :
      std::invalid_argument(what_arg)
    {
      std::cout << "[cvc::invalid_argument]: exception in " <<
        function_name << " : " <<
        what_arg << std::endl;
    }
};

} // namespace(cvc)
