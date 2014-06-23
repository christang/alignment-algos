/**
 *  Package HMAP2.1
 *  File: noalib.h
 *  Desc: Near-optimal alignment parameters
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_NOALIB
#define _HMAP2_NOALIB

#include "pstore.h"

class NOaliParams {

public:

  NOaliParams ();
  NOaliParams (const NOaliParams& p);

  void read (ParamStore* ps);

  static const int   default_number_suboptimal;
  static const float default_delta_ratio;
  static const unsigned int default_k_limit;
  static const unsigned int default_sort_limit;
  static const unsigned int default_user_limit;
  static const float default_max_overlap;
  static const unsigned int default_rounds;

  int   number_suboptimal;
  int   subopt_per_round;
  float delta_ratio;
  unsigned int k_limit;
  unsigned int sort_limit;
  unsigned int user_limit;
  float max_overlap;
  float final_overlap;
  unsigned int rounds;

};

#endif  // _HMAP2_NOALIB
