/**
 *  Package HMAP2.1
 *  File: argv.h
 *  Desc: Process command line arguments
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_ARGV
#define _HMAP2_ARGV

#include <string>
#include <sstream>
#include <vector>

#include "pstore.h"

using namespace std;

class Argv 
  : public ParamStore
{

public:

  Argv (int argc, const char** argv);

  int count();
  bool getArg(int c, string& arg);
  stringstream& getArg(int c,
		       bool cleanbuff=true,
		       bool eraseafter=false);
  bool getSwitch(const char* sw,
		 bool eraseafter=true);
  stringstream& getSwitch(const char* sw, int c,
			  bool cleanbuff=true,
			  bool eraseafter=true);
  bool help();

private:

  void init(int argc, const char** argv);

  bool dohelp;
  vector<string> args_store;
  stringstream argv_buffer;

};

template <class param_t>
Argv& operator>> (Argv& arg, param_t& p)
{
  p.read (&arg);
  return arg;
}

#endif  // _HMAP2_ARGV
