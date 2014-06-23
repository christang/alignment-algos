
#include "formats.h"

Formats::FastaOut::FastaOut (int len) : line_length (len) {}

Formats::FastaIn::FastaIn (const char* cs, bool ht) 
  : head_tail(ht),find_me(cs) {}

Formats::FastaAlignmentIn::FastaAlignmentIn (bool ht)
  : head_tail(ht) {}

Formats::PIROut::PIROut (int len) : line_length (len) {}

Formats::PIRIn::PIRIn (bool ht)
  : head_tail(ht) {}

Formats::HMAPOut::HMAPOut (const char* submtx,int len) 
  : line_length (len), submatrix(submtx) {}

