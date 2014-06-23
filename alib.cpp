/**
 *  Package HMAP2.1
 *  File: alib.cpp
 *  Desc: Alignment parameters
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include <iostream>
#include "alib.h"

const align_t         AliParams::default_align_type        = semi_local;
const float           AliParams::default_gap_init_penalty  = 4.73;
const float           AliParams::default_gap_extn_penalty  = 0.34;

AliParams::AliParams ()
  : align_type       (default_align_type),
    gap_init_penalty (default_gap_init_penalty),
    gap_extn_penalty (default_gap_extn_penalty),
    submatrix_fn     ("")
{
  // Empty
}

void AliParams::read (ParamStore* p)
{
  string s;

  s="ALIGN_MODE";
  if (p->find(s)) {
    int v=default_align_type; p->getValue(s) >> v;
    align_type = static_cast<align_t> (v);
  }
  
  s="GAP_INIT_PENALTY";
  if (p->find(s)) p->getValue(s) >> gap_init_penalty;
  
  s="GAP_EXTN_PENALTY";
  if (p->find(s)) p->getValue(s) >> gap_extn_penalty;

  s="SUB_MATRIX";
  if (p->find(s)) p->getValue(s) >> submatrix_fn;
}
