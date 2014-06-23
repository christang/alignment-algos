/**
 *  Package HMAP2.1
 *  File: application.h
 *  Desc: Application parameters
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_APPLICATION
#define _HMAP2_APPLICATION

#include <string>

#include "pstore.h"

enum output_format_t {
  oHMAP                = 0,
  oPIR                 = 1,
  oFASTA               = 2
};

class ApplicationParams
{

 public:
  
  ApplicationParams ();
  
  void read (ParamStore* p);
  
  static const output_format_t default_output_format;
  static const int             default_line_length;
  static const int             default_verbosity;
  static const string          default_log_file;
  
  output_format_t output_format;
  int             line_length;
  int             verbosity;
  string          log_file;
  
};

#endif  // _HMAP2_APPLICATION
