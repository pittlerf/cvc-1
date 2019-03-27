/***********************************************************************
 *  
 * Copyright (C) 2018 Bartosz Kostrzewa
 *
 ***********************************************************************/

#ifndef DEBUG_PRINTF_HPP
#define DEBUG_PRINTF_HPP

#include "global.h"

#include <cstdarg>
#include <cstdio>

/* Function along the lines of printf which produces output on a single
 * or all MPI tasks (unordered) when g_debug_level is at or
 * above the provided threshold 
 * to have output by all MPI tasks, simply pass g_proc_id for proc_id */

namespace cvc {

static inline void debug_printf(int const proc_id,
                                int const verbosity_level_threshold,
                                char const * format,
                                ...)
{
  if( g_proc_id == proc_id && g_verbose >= verbosity_level_threshold ){
    va_list arglist;
    va_start(arglist, format);
    vprintf(format, arglist);
    va_end(arglist);
    fflush(stdout);
  }
}

/**
 * @brief printf-like function which prepends MPI task ID to error message
 *
 * @param format format string
 * @param ... variables for placeholders
 */
static inline void error_printf(char const * format, ...)
{
  size_t fmtlen = strlen(format);
  if( fmtlen > 10000 ){
    fmtlen = 10000;
  }
  char * newformat = (char)malloc(fmtlen+100);
  if( newformat == NULL ){
    EXIT_WITH_MSG(CVC_EXIT_MALLOC_FAILURE, "malloc failure in error_printf!\n");
  }
  snprintf(newformat, fmtlen+100, "ERROR on MPI Task: %d -- %s", g_proc_id, format);
  va_list arglist;
  va_start(arglist, format);
  vprintf(newformat, arglist);
  va_end(arglist);
  fflush(stdout);
}

} // namespace(cvc)

#endif
