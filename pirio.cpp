
#include "pirio.h"

PIRWrite::PIRWrite (ostream& o, int len) 
  : output (&o), line_length (len) {}

PIRWrite operator<< (ostream& o, Formats::PIROut p)
{ return PIRWrite (o,p.line_length); }

void PIRWrite::write (const string& seq)
{
  int size = seq.size();

  for (int i=0; i<size; i+=line_length) {
    *output << seq.substr(i,line_length);
    *output << endl;
  }
}

void PIRWrite::fix_ends (string& seq)
{
  if (seq[0]==SequenceElem::Head)
    seq.erase(0,1);
  if (seq[seq.length()-1]==SequenceElem::Tail)
    seq.erase(seq.length()-1);
}

PIRRead::PIRRead (istream& i, bool flag)
  : input (&i), head_tail (flag) { }

PIRRead operator>> (istream& i, Formats::PIRIn p)
{ return PIRRead (i,p.head_tail); }
