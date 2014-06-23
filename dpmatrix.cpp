/**
 *  Package HMAP2.1
 *  File: dpmatrix.cpp
 *  Desc: Template class for dynamic programming matrix. 
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include "dpmatrix.h"

const int DPCell::null = -1;

DPCell::DPCell() 
  : prev_query_idx    (null),
    prev_template_idx (null),
    query_idx         (null),
    template_idx      (null),
    score             (0.0f)
{
  // Empty
}

void DPCell::setTB (int pq, int pt, float s)
{
  prev_query_idx = pq;
  prev_template_idx = pt;
  score = s;
}
