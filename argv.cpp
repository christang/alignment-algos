/**
 *  Package HMAP2.1
 *  File: argv.cpp
 *  Desc: Process command line arguments
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include <stdlib.h>
#include "argv.h"

Argv::Argv (int argc,
	    const char** argv) 
  : dohelp (false)
{
  init(argc,argv);
}

int Argv::count()
{
  return (int)args_store.size();
}

stringstream& Argv::getArg(int c,
			   bool cleanbuff,
			   bool eraseafter)
{
  if (c>=(int)args_store.size()) throw string ("Command line arg missing");
  if (cleanbuff) { argv_buffer.str(""); argv_buffer.clear(); }
  argv_buffer << args_store[c] << flush;
  if (eraseafter) args_store.erase(args_store.begin()+c);
  return argv_buffer;
}

bool Argv::getArg(int c, string& arg)
{
  if (c>=(int)args_store.size()) return false;
  arg = args_store[c];
  return true;
}

bool Argv::getSwitch(const char* sw,
		     bool eraseafter)
{
  for (vector<string>::iterator i=args_store.begin();
       i!=args_store.end();i++) {
    if (*i==sw) {
      if (eraseafter) args_store.erase(i);
      return true;
    }
  }
  return false;
}

stringstream& Argv::getSwitch(const char* sw, 
			      int c,
			      bool cleanbuff,
			      bool eraseafter)
{
  bool s = false;
  vector<string>::iterator i,beg_i;

  if (cleanbuff) { argv_buffer.str(""); argv_buffer.clear(); }
  for (i=args_store.begin();i!=args_store.end();i++) {
    if (!s && *i==sw) {
      beg_i = i;
      s = true;
    } else if (s && c>0) {
      argv_buffer << *i << " ";
      --c;
    } else if (s && c==0) break;
  }
  if (c>0) throw string("Switch arg missing for ")+sw ;
  if (s && eraseafter) args_store.erase(beg_i,i);
  return argv_buffer;
}

bool Argv::help()
{
  return dohelp;
}

void Argv::init(int argc,
		const char** argv)
{
  for (int i=1; i<argc; ++i) {
    string s = argv[i];
    if (s == "-help") dohelp = true;
  }

  for (int i=1; i<argc; ++i) {
    string s = argv[i];
    if (s.find ("--") == 0) {
      if (i+1>=argc) throw string ("Argument missing for ")+s;
      string key(argv[i++]);
      string value(argv[i]);
      key = key.substr(2);
      setValue(key,value);
    }
    else {
      args_store.push_back(s);
    }
  }
}
