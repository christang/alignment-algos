/**
 *  Package HMAP2.1
 *  File: rcfile.h
 *  Desc: Resource file parameters
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_RCFILE
#define _HMAP2_RCFILE

#include "pstore.h"

class RCfile :
  public ParamStore {
  
public:
  
  RCfile ();
  RCfile (const string& fname);

  static const std::string default_RC_fname;

private:
  
  void init ();
  void resolve (string& fname);

  string fname;

};

template <class param_t>
RCfile& operator>> (RCfile& arg, param_t& p)
{
  p.read (&arg);
  return arg;
}

#endif  // _HMAP2_RCFILE
