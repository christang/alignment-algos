/**
 *  Package HMAP2.1
 *  File: sequence.cpp
 *  Desc: Base class describing sequences.
 *
 *  Created on 9/11/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include "sequence.h"

const char SequenceElem::Head = '^';
const char SequenceElem::Tail = '$';

SequenceElem::SequenceElem ()
  : index (-1),
    olc   (' ')  {}

SequenceElem::SequenceElem (const SequenceElem& e)
  : index (e.index),
    olc   (e.olc)  {}

SequenceElem::SequenceElem (int i, char o)
  : index (i),
    olc   (o) {}

SequenceElem& SequenceElem::operator= (const SequenceElem& e)
{
  if (&e == this) return *this;
  index = e.index;
  olc   = e.olc;
  return *this;
}
