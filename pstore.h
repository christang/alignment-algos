/**
 *  Package HMAP2.1
 *  File: pstore.cpp
 *  Desc: Parameter store
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_PSTORE
#define _HMAP2_PSTORE

#include <map>
#include <string>
#include <sstream>

using namespace std;

class ParamStore
{
  
public:

  void clear();

  bool find(const string& key);
  stringstream& getValue(const string& key,
			 bool cleanbuff=true,
			 bool eraseafter=false
			 );
  bool setValue(const string& key,
		const string& value
		);

protected:
  
  ParamStore ();
  bool extract(istream& in,
	       string& key,
	       string& value
	       );
  void parseline(string& line,
		 string& key,
		 string& value
		 );

  map <string, string> map_store;
  stringstream pstore_buffer;

};

#endif  // _HMAP2_PSTORE
