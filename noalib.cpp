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

#include "noalib.h"

const int   NOaliParams::default_number_suboptimal = 200;
const float NOaliParams::default_delta_ratio       = 0.01f;
const unsigned int NOaliParams::default_k_limit    = 16;
const unsigned int NOaliParams::default_sort_limit = 100;
const unsigned int NOaliParams::default_user_limit = 100000;
const float NOaliParams::default_max_overlap       = 0.30;
const unsigned int NOaliParams::default_rounds     = 4;

NOaliParams::NOaliParams ()
  : number_suboptimal (default_number_suboptimal),
    subopt_per_round  (default_number_suboptimal),
    delta_ratio       (default_delta_ratio),
    k_limit    (default_k_limit),
    sort_limit (default_sort_limit),
    user_limit (default_user_limit),
    max_overlap  (default_max_overlap),
    final_overlap(default_max_overlap),
    rounds     (default_rounds)
{
  // Empty
}

NOaliParams::NOaliParams (const NOaliParams& p)
  : number_suboptimal (p.number_suboptimal),
    subopt_per_round  (p.subopt_per_round),
    delta_ratio       (p.delta_ratio),
    k_limit           (p.k_limit),
    sort_limit        (p.sort_limit),
    user_limit        (p.user_limit),
    max_overlap       (p.max_overlap),
    final_overlap     (p.final_overlap),
    rounds            (p.rounds)
{
  // Empty
}

void NOaliParams::read (ParamStore* p)
{
  string s;
  s="NUM_SUBOPT"; if (p->find(s)) p->getValue(s) >> number_suboptimal;
  s="NUM_ROUND_SUBOPT"; if (p->find(s)) p->getValue(s) >> subopt_per_round;
  s="DELTA_RATIO"; if (p->find(s)) p->getValue(s) >> delta_ratio;
  s="K_LIMIT";     if (p->find(s)) p->getValue(s) >> k_limit;
  s="USER_LIMIT";  if (p->find(s)) p->getValue(s) >> user_limit;
  s="SORT_LIMIT";  if (p->find(s)) p->getValue(s) >> sort_limit;
  s="MAX_OVERLAP"; if (p->find(s)) p->getValue(s) >> max_overlap;
  s="FINAL_OVERLAP"; if (p->find(s)) p->getValue(s) >> final_overlap;
  s="ROUNDS";      if (p->find(s)) p->getValue(s) >> rounds; 
}
