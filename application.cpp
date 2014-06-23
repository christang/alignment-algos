/**
 *  Package HMAP2.1
 *  File: application.cpp
 *  Desc: Application parameters
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */


#include "application.h"

const output_format_t ApplicationParams::default_output_format = oFASTA;
const int             ApplicationParams::default_line_length   = 60;
const int             ApplicationParams::default_verbosity     = 0;
const string          ApplicationParams::default_log_file      = "";

ApplicationParams::ApplicationParams ()
  : output_format (default_output_format),
    line_length   (default_line_length),
    verbosity     (default_verbosity),
    log_file      (default_log_file)
{
  // Empty
}
  
void ApplicationParams::read (ParamStore* p)
{
  string s;
  
  s="OUTPUT_FORMAT";
  if (p->find(s)) {
    int v=default_output_format; p->getValue(s) >> v;
    output_format = static_cast<output_format_t> (v);
  }

  s="OUTPUT_LINE_LENGTH";
  if (p->find(s)) p->getValue(s) >> line_length;

  s="VERBOSE";
  if (p->find(s)) p->getValue(s) >> verbosity;

  s="LOG_FILE";
  if (p->find(s)) p->getValue(s) >> log_file;
}
