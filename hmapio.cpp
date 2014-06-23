
#include "hmapio.h"

HMAPWrite::HMAPWrite (ostream& o, int len, const string& sfn) 
  : output (&o), line_length (len), submatrix_fn (sfn) {}

HMAPWrite operator<< (ostream& o, Formats::HMAPOut p)
{ return HMAPWrite (o,p.line_length,p.submatrix); }

void HMAPWrite::write (const string& seq_templ_sse,
		       const string& seq_templ,
		       const string& seq_marks,
		       const string& seq_query,
		       const string& seq_query_sse)
{
  int size = seq_templ.size();

  for (int i=0; i<size; i+=line_length) {
    *output << endl;
    *output << "       " << seq_templ_sse.substr(i,line_length) << endl;
    *output << "model: " <<     seq_templ.substr(i,line_length) << endl;
    *output << "       " <<     seq_marks.substr(i,line_length) << endl;
    *output << "query: " <<     seq_query.substr(i,line_length) << endl;
    *output << "       " << seq_query_sse.substr(i,line_length) << endl;
  }
}

void HMAPWrite::fix_ends (string& seq)
{
  if (seq[0]==SequenceElem::Head)
    seq.erase(0,1);
  if (seq[seq.length()-1]==SequenceElem::Tail)
    seq.erase(seq.length()-1);
}
