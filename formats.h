
#ifndef _HMAP2_FORMATS
#define _HMAP2_FORMATS

#include <iostream>
#include <sstream>

#include "alignment.h"
#include "gstrings.h"

struct Formats
{
  struct FastaOut {
    FastaOut (int len=60);
    int line_length;
  };
  struct FastaIn {
    FastaIn (const char* cs="", bool flag=true);
    bool head_tail;
    string find_me;
  };
  struct FastaAlignmentIn {
    FastaAlignmentIn (bool flag=true);
    bool head_tail;
  };
  struct PIROut {
    PIROut (int len=60);
    int line_length;
  };
  struct PIRIn {
    PIRIn (bool flag=true);
    bool head_tail;
  };
  struct HMAPOut {
    HMAPOut (const char* submatrix="", int len=60);
    int line_length;
    string submatrix;
  };
};

#endif  //_HMAP2_FORMATS
