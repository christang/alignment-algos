/**
 *  Package HMAP2.1
 *  File: enumerator.h
 *  Desc: Base template class for enumerators in this library.
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_ENUMERATOR
#define _HMAP2_ENUMERATOR

template <class S1, class S2, class Etype> class DPMatrix;
template <class S1, class S2, class Etype> class AlignmentSet;

template <class S1, class S2, class Etype>
class Enumerator {
public:
  virtual int estimateSize () const = 0;
  virtual void enumerate (DPMatrix<S1,S2,Etype>& dpm,
			  AlignmentSet<S1,S2,Etype>& as) = 0;
};

#endif  // _HMAP2_ENUMERATOR
