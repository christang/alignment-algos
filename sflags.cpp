/**
 *  Package HMAP2.1
 *  File: sflags.cpp
 *  Desc: Handle suboptimal region flags.
 *
 *  Created on 11/26/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include "sflags.h"

SuboptFlags::SuboptFlags (bool f, size_t len)
  : subopt(f,len),last_elem(0) 
{
  Base::seq_name="Flags=suboptimal region"; 
  Base::seq_string="";
  setString();
}

void SuboptFlags::append (const string& s)
{
  for (string::const_iterator it=s.begin(); it!=s.end(); ++it) {
    if (last_elem>=subopt.size()) 
      throw string ("Sequence flags longer than template!");
    if (*it == '0') subopt[last_elem++]=false;
    else subopt[last_elem++]=true;
  }
}

void SuboptFlags::append (const char* cs)
{ 
  string s(cs);
  append(s);
}

size_t SuboptFlags::size () const 
{ return subopt.size(); }

void SuboptFlags::setString ()
{
  seq_string=="";
  for (int i=0; i<(int)subopt.size(); ++i) {
    if (subopt[i]) seq_string.append("1");
    else seq_string.append("0");
  }
}

void SuboptFlags::Set (unsigned int i,bool b)
{ 
  if(i>subopt.size()) throw string("Subopt index out of range");
  
  subopt[i]=b; 
  if (b) seq_string.replace(i,1,"1");
  else seq_string.replace(i,1,"0");
}
