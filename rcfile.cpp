/**
 *  Package HMAP2.1
 *  File: rcfile.cpp
 *  Desc: Resource file parameters
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include <fstream>
#include <iostream>

#include "rcfile.h"

const string RCfile::default_RC_fname = "~/.hmaprc";

RCfile::RCfile ()
  : fname (default_RC_fname)
{
  resolve(fname);
  ifstream fin (fname.c_str());
  if (!fin.good ()) 
    cerr << "No defaults file (" << default_RC_fname
         << ").  Using programmed defaults." << endl;
     //throw default_RC_fname + " defaults not found";
  else init();
}

RCfile::RCfile (const string& fn)
  : fname (fn)
{
  resolve(fname);
  init();
}

void RCfile::init ()
{
  ifstream fin (fname.c_str());
  if (fin.good ()) 
    cerr << "Loading params from " << fname << endl;
  else throw fname + " file not found";

  string key;
  string value;
  while (extract(fin,key,value)) map_store[key] = value;
}

void RCfile::resolve (string& fn)
{
  if (fn.substr(0,1)=="~") {
    char *hdir=getenv("HOME");
    string file;
    if(hdir==NULL) { file="."; }
    else file=hdir;
    file+=fn.substr(1);
    fn=file;
  }
}
