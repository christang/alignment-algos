
#include "fastaio.h"

FastaWrite::FastaWrite (ostream& o, int len) 
  : output (&o), line_length (len) {}

FastaWrite operator<< (ostream& o, Formats::FastaOut p)
{ return FastaWrite (o,p.line_length); }

void FastaWrite::write (const string& seq)
{
  int size = (int)seq.size();
  for (int i=0; i<size; i+=line_length) {
    *output << seq.substr(i,line_length);
    *output << endl;
  }
}

FastaRead::FastaRead (istream& i, const string& s, bool flag)
  : input (&i), head_tail (flag), find_me (s) {}

FastaRead operator>> (istream& i, Formats::FastaIn p)
{ return FastaRead (i,p.find_me,p.head_tail); }

FastaAlignmentRead::FastaAlignmentRead (istream& i, bool flag)
  : Base (i,"",flag) {}

FastaAlignmentRead operator>> (istream& i, Formats::FastaAlignmentIn p)
{ return FastaAlignmentRead (i,p.head_tail); }
