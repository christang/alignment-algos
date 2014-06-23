/**
 *  Package HMAP2.1
 *  File: alib.h
 *  Desc: Alignment parameters
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_ALIB
#define _HMAP2_ALIB

#include <string>

#include "pstore.h"

enum align_t {      // alignment overhang treatment
  global_local         = 0,  // overhangs penalized in template not query
  global               = 1,  // overhangs penalized (local_const = min_score)
  local_global         = 2,  // overhangs penalized in query not template
  local                = 3,  // local alignment (local_const = 0)
  semi_local           = 4   // not penalized   (local_const = min_score)
};

class AliParams 
{
  
public:
  
  AliParams ();
  
  void read (ParamStore* p);
  
  static const align_t         default_align_type;
  static const float           default_gap_init_penalty;
  static const float           default_gap_extn_penalty;

  align_t         align_type;
  float           gap_init_penalty;
  float           gap_extn_penalty;
  string          submatrix_fn;

};

#endif  // _HMAP2_ALIB
