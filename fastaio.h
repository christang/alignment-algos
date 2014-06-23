
#ifndef _HMAP2_FASTAIO
#define _HMAP2_FASTAIO

#include "formats.h"
#include "sequence.h"

//----------------------------------------------------------------------------

class FastaWrite 
{
public:
  FastaWrite (ostream& o, int len);
  
  template <class S1, class S2, class Etype>
  void write (AlignmentSet<S1,S2,Etype>& as);

  template <class elem_t>
  void write (Sequence<elem_t>& seq);

  void write (const string& s);

  template <class S1, class S2>
  void makeAnnotation (AlignedPairList<S1,S2>& ali,string& s);

  ostream* output;
  int line_length;
};

FastaWrite operator<< (ostream& o, Formats::FastaOut p);

template <class S1, class S2, class Etype>
ostream& operator<< (FastaWrite w, AlignmentSet<S1,S2,Etype>& as)
{ w.write (as); return *w.output; };

template <class elem_t>
ostream& operator<< (FastaWrite w, Sequence<elem_t>& seq)
{ w.write (seq); return *w.output; };

template <class elem_t>
void FastaWrite::write (Sequence<elem_t>& seq)
{
  *output << "> ";  
  *output << seq.seq_name;
  *output << endl;

  write (*seq.getString());
}

template <class S1, class S2, class Etype>
void FastaWrite::write (AlignmentSet<S1,S2,Etype>& as)
{
  SequenceGaps gaps(as);

  *output << "> ";
  *output << as.getTemplateSequence()->seq_name;
  *output << endl;

  string gapped_string; 
  gaps.build (*as.getTemplateSequence()->getString(),gapped_string);
  write (gapped_string);
  
  int count = 0;
  for (typename AlignmentSet<S1,S2,Etype>::iterator it=as.begin(); 
       it!=as.end(); ++it) {
    *output << "> ";
    *output << as.getQuerySequence()->seq_name;
    *output << "_" << count++;
    string annot; makeAnnotation (*it,annot);
    if (annot!="") *output << " " << annot;
    *output << endl;

    gaps.build (*as.getQuerySequence()->getString(),*it,gapped_string);
    write (gapped_string);
  }
}

template <class S1, class S2>
void FastaWrite::makeAnnotation (AlignedPairList<S1,S2>& ali,string& s)
{
  stringstream buff("");
  buff << "(sc=";
  buff << ali.score;
  buff << ",ev=";
  buff << ali.significance;
  buff << ",id=";
  buff << ali.identity;
  buff << "%)";
  s.assign(buff.str());
}

//----------------------------------------------------------------------------

class FastaRead {

public:
  FastaRead (istream& i, const string& s="", bool flag=true);
  
  template <class S>
  void readInto (S& sequence);

  template <class S>
  bool search_failed (S& sequence);

  istream* input;
  bool head_tail;
  string find_me;
};

FastaRead operator>> (istream& i, Formats::FastaIn p);

template <class S>
istream& operator>> (FastaRead r, S& s)
{ r.readInto (s); return *r.input; }

template <class S>
void FastaRead::readInto (S& s)
{
  while (input->good() && search_failed(s)) {}

  if (!input->good()) {
    if (find_me=="") throw string ("Error reading fasta file");
    else throw string ("Could not find search string: ") + find_me;
  }

  if (head_tail) s.append ("^");

  char c;
  string buff;
  while (input->good()) {
    c = input->peek();
    if (c=='>') break;
    getline (*input,buff);
    s.append (buff);
  }

  if (head_tail) s.append ("$");
}

template <class S>
bool FastaRead::search_failed(S& s) 
{
  char c = input->peek();
  string buff = "";
  if (c=='>') {
    c = input->get();
    if (input->peek()==' ') input->ignore(32,' ');
    getline (*input,buff);
    if (find_me=="" || buff.find(find_me)!=string::npos) {
      s.seq_name = buff;
      return false;
    }
  } else getline (*input,buff);
  return true;
}

template <>
inline bool FastaRead::search_failed<string> (string& s) 
{
  char c = input->peek();
  string buff = "";
  if (c=='>') {
    c = input->get();
    if (input->peek()==' ') input->ignore(32,' ');
    getline (*input,buff);
    return false;
  } else getline (*input,buff);
  return true;
}

//----------------------------------------------------------------------------

class FastaAlignmentRead : public FastaRead {

  typedef FastaRead Base;

 public:
  FastaAlignmentRead (istream& i, bool flag=true);

  template <class S>
  void readEntries (S& apl);

};

FastaAlignmentRead operator>> (istream& i, Formats::FastaAlignmentIn p);

template <class S>
istream& operator>> (FastaAlignmentRead r, S& apl)
{ r.readEntries (apl); return *r.input; }

template <class S>
void FastaAlignmentRead::readEntries (S& apl)
{
  // S must be an AlignedPairlist

  string templ_tmp_ali, query_tmp_ali;
  
  Base::readInto (templ_tmp_ali);
  Base::readInto (query_tmp_ali);

  apl.readFrom (query_tmp_ali, templ_tmp_ali);

}

#endif  //_HMAP2_FASTAIO
