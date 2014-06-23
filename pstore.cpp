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

#include <iostream>
#include "pstore.h"

ParamStore::ParamStore ()
  : pstore_buffer ()
{
  // Empty
}

void ParamStore::clear ()
{
  map_store.clear ();
}

bool ParamStore::find (const string& key)
{
  return map_store.find(key) != map_store.end();
}

stringstream& ParamStore::getValue (const string& key,
				    bool cleanbuff,
				    bool eraseafter)
{
  if (cleanbuff) { pstore_buffer.str(""); pstore_buffer.clear(); }
  pstore_buffer << map_store[key] << flush;
  if (eraseafter) map_store.erase(map_store.find(key));
  return pstore_buffer;
}

bool ParamStore::setValue (const string& key,
			   const string& value)
{
  map_store[key] = value;
  return true;
}

bool ParamStore::extract (istream& in,
			  string& key,
			  string& value)
{
  string line("");
  while (line == "" || line[0] == '#') {
    if (in.eof()) return false;
    getline(in,line);
  }
  key = "";
  value = "";
  parseline (line, key, value);
  return true;
}

void ParamStore::parseline (string& line,
			    string& key,
			    string& value)
{
  if (line.empty()) return;

  unsigned int i0 = line.find_first_of(':',0);
  unsigned int i1 = line.find_first_not_of(" \t",i0+1);
  if (i0==string::npos)
    throw string ("Param parse error");
  key = line.substr(0,i0);
  if (i1==string::npos) value = "";
  else value = line.substr(i1);
}
