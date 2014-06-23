/**
 *  Package HMAP2.1
 *  File: significance.h
 *  Desc: Derive significance values from defined distribution.
 *
 *  Created on 11/27/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_SIGNIFICANCE
#define _HMAP2_SIGNIFICANCE

template <class Stype>
class Significance {

public:
  
  float significance (float score) const ;

protected:
  
  Significance () {};
  inline const Stype& Derived() const { return static_cast<const Stype&>(*this); } ;
  
};

template <class Stype>
float Significance<Stype>::significance (float score) const
{ return Derived().significance(score); }

#endif  //_HMAP2_SIGNIFICANCE
