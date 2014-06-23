
#ifndef _HMAP2_HMAPIO
#define _HMAP2_HMAPIO

#include "formats.h"
#include "sequence.h"
#include "submatrix.h"

//----------------------------------------------------------------------------

class HMAPWrite
{
public:
  HMAPWrite (ostream& o, int len, const string& sfn);
  
  template <class S1, class S2, class Etype>
  void write (AlignmentSet<S1,S2,Etype>& as);
  
  void write (const string& seq_templ_sse,
	      const string& seq_templ,
	      const string& seq_marks,
	      const string& seq_query,
	      const string& seq_query_sse);

  void fix_ends (string& s);

  template <class S1, class S2>
  void makeAnnotation (AlignedPairList<S1,S2>& ali,string& s);

  template <class S1, class S2, class Etype>
  void generateMarks (AlignedPairList<S1,S2>& ali,
		      AlignmentSet<S1,S2,Etype>& as,
		      string& s);

  ostream* output;
  int line_length;
private:
  string submatrix_fn;
};

HMAPWrite operator<< (ostream& o, Formats::HMAPOut p);

template <class S1, class S2, class Etype>
ostream& operator<< (HMAPWrite w, AlignmentSet<S1,S2,Etype>& as)
{ w.write (as); return *w.output; };

template <class S1, class S2, class Etype>
void HMAPWrite::write (AlignmentSet<S1,S2,Etype>& as)
{
  int count= 0;
  string gapped_templ_sse,gapped_templ,gapped_marks;
  string gapped_query,gapped_query_sse; 
  valarray<bool> mask (false,as.size());
  for (typename AlignmentSet<S1,S2,Etype>::iterator it=as.begin(); 
       it!=as.end(); ++it) {
    mask[count]=true;
    SequenceGaps gaps(as,mask);

    *output << ">" << as.getQuerySequence()->seq_name;
    *output << "_" << count;
    string annot; makeAnnotation (*it,annot);
    if (annot!="") *output << " " << annot;
    *output << endl << endl;

    *output << "model: length " << as.getTemplateSequence()->size()-2 << endl;
    *output << "query: length " << as.getQuerySequence()->size()-2 << endl;

    gaps.build (*as.getTemplateSequence()->getSSEString(),
		gapped_templ_sse,' ');
    fix_ends (gapped_templ_sse);

    gaps.build (*as.getTemplateSequence()->getString(),gapped_templ);
    fix_ends (gapped_templ);

    generateMarks (*it, as, gapped_marks);

    gaps.build (gapped_marks,*it,gapped_marks,' ');

    fix_ends (gapped_marks);

    gaps.build (*as.getQuerySequence()->getString(),*it,gapped_query);
    fix_ends (gapped_query);

    gaps.build (*as.getQuerySequence()->getSSEString(),*it,
		gapped_query_sse,' ');
    fix_ends (gapped_query_sse);

    write(gapped_templ_sse,gapped_templ,gapped_marks,
	  gapped_query,gapped_query_sse);

    *output << endl;

    mask[count++]=false;
  }  

}

template <class S1, class S2>
void HMAPWrite::makeAnnotation (AlignedPairList<S1,S2>& ali,string& s)
{
  stringstream buff("");
  buff << "(sc=";
  buff << ali.score;
  buff << ",ev=";
  buff << ali.significance;
  buff << ",id=";
  buff << ali.identity;
  buff << "%)";
  buff << "  UID=" << ali.uid;
  s.assign(buff.str());
}

template <class S1, class S2, class Etype>
void HMAPWrite::generateMarks (AlignedPairList<S1,S2>& ali,
			       AlignmentSet<S1,S2,Etype>& as,
			       string& marks)
{
  BlosumMatrix* bm = 0;
  if (submatrix_fn!="") {
    bm = new BlosumMatrix (submatrix_fn.c_str());
  }

  int qp = -1;
  stringstream buffer("");
  const string* q_seq = as.getQuerySequence()->getString();
  const string* t_seq = as.getTemplateSequence()->getString();

  for (typename AlignedPairList<S1,S2>::iterator it=ali.begin(); 
       it!=ali.end(); ++it) {

    int qi = it->query_idx();
    int ti = it->template_idx();

    char qc = (*q_seq)[qi];
    char tc = (*t_seq)[ti];

    float s = as.getDPMatrix()->getSim(qi,ti);

    //    cerr << "qi: " << qi << "\tqp: " << qp << endl;

    buffer << string (qi-qp-1,' ');
    qp = qi;

    if( ( qc == SequenceElem::Head ) ||
	( qc == SequenceElem::Tail) ) {
      buffer << qc;
    }
    else if( qc == tc ) {
      buffer << '|';
    }
    else if( bm && ( bm->score(qc,tc)>0) ) {
      buffer << ':';
    }
    else if( s > 0 ) {
      buffer << '.';
    }
    else {
      buffer << ' ';
    }

  }

  marks = buffer.str();
}

#endif  //_HMAP2_HMAPIO
