/**
 *  Package HMAP2.1
 *  File: gstrings.h
 *  Desc: Generate gapped strings for displaying alignments
 *
 *  Created on 11/26/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include "gstrings.h"

const char SequenceGaps::gapChar = '-';

void SequenceGaps::build (const string& seq, 
			  string& value, 
			  char gc)
{
  string result = "";
  assert (template_len==(int)seq.size());
  for (int i=0; i<template_len-1; ++i) {
    result.append(1,seq[i]);
    result.append(anchors[i],gc);
  }
  result.append(1,seq[template_len-1]);
  value = result;
}
