/**
 *  Package HMAP2.1
 *  File: hmapalib_seq.h
 *  Desc: Sequence class for HMAP profiles.
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include <fstream>

#include "hmapalib_seq.h"
#include "hmath.h"

#ifdef WIN32
double vs_erfc(double);
#endif

HMAPElem::HMAPElem ()
  : aa_profile   (20), 
    gap_values   (4), 
    sse_values   (3)

{
  // Empty
}

HMAPElem::HMAPElem (istream& in)
  : aa_profile   (20), 
    gap_values   (4), 
    sse_values   (3)
{
  readHMAP (in);
}

HMAPElem::HMAPElem (const HMAPElem& e) 
  : SequenceElem   (e),
    aa_profile     (e.aa_profile), 
    gap_values     (e.gap_values), 
    motif_value    (e.motif_value), 
    motif_confid   (e.motif_confid),
    sse_values     (e.sse_values), 
    sse_confid     (e.sse_confid),
    surfacc_value  (e.surfacc_value), 
    surfacc_confid (e.surfacc_confid)
{
  // Empty
}


HMAPElem& HMAPElem::operator= (const HMAPElem& e)
{
  if (&e == this) return *this;
  aa_profile = e.aa_profile;
  gap_values = e.gap_values;
  motif_value = e.motif_value;
  motif_confid = e.motif_confid;
  sse_values = e.sse_values; 
  sse_confid = e.sse_confid;
  surfacc_value = e.surfacc_value; 
  surfacc_confid = e.surfacc_confid;
  return *this;
}

void HMAPElem::readHMAP (istream& in)
{
  //  cerr << "readHMAP: top" << endl;

  int  temp_i;
  char temp_c;

  in >> temp_i; // this is not used
  in >> olc;
  for (int i=0; i<20; i++) {
    in >> aa_profile[i];
    aa_profile[i]/=100.0f;
  }
  calcHydropathy();
  in >> temp_c;
  if (temp_c != '-') 
    throw string ("Parse error before '-'");
  for (int i=0; i<4; i++) {
    in >> gap_values[i];
  }
  in >> motif_value;
  in >> motif_confid;
  in >> temp_c;
  if (temp_c != '*') 
    throw string ("Parse error before '*'");
  for (int i=0; i<3; i++) {
    in >> sse_values[i];
  }
  in >> sse_confid;
  in >> surfacc_value;
  in >> surfacc_confid;

  // Assign lods type for Hmap sequence

  unsigned int idxtype = 3;
  if (sse_values[0]>.5f) idxtype = 0;
  if (sse_values[1]>.5f) idxtype = 1;
  if (sse_values[2]>.5f) idxtype = 2;

  unsigned int idxconf = 0;
  if (sse_confid > .33f) idxconf = 1;
  if (sse_confid > .66f) idxconf = 2;

  lods_type = idxtype * 3 + idxconf;

  char buff[10];
  in.getline(buff,10); // skip to next line

  //  cerr << "readHMAP: end" << endl;
}

void HMAPElem::calcHydropathy ()
{
  valarray<float> hpath (20);

  hpath [0] =  0.5f;  // A 
  hpath [1] = -2.2f;  // R 
  hpath [2] = -1.0f;  // N 
  hpath [3] = -1.3f;  // D 
  hpath [4] =  1.0f;  // C 
  hpath [5] = -1.4f;  // Q 
  hpath [6] = -2.1f;  // E 
  hpath [7] =  0.0f;  // G 
  hpath [8] = -0.5f;  // H 
  hpath [9] =  0.9f;  // I 

  hpath [10] =  0.8f;  // L 
  hpath [11] = -3.5f;  // K 
  hpath [12] =  0.6f;  // M 
  hpath [13] =  0.7f;  // F 
  hpath [14] = -0.8f;  // P 
  hpath [15] = -0.3f;  // S 
  hpath [16] = -0.2f;  // T 
  hpath [17] =  0.3f;  // W 
  hpath [18] =  0.1f;  // Y 
  hpath [19] =  0.8f;  // V 

  hydropathy = 0.f;
  for (unsigned int j=0; j<20; ++j)
    hydropathy += aa_profile [j] * hpath [j];
}

HMAPSequence::HMAPSequence ()
  : evd1_field (0.0f),
    evd2_field (0.0f),
    sse_string ("")
{
  // Empty
}

HMAPSequence::HMAPSequence (const char* fn)
  : evd1_field (0.0f),
    evd2_field (0.0f),
    sse_string ("")
{
  ifstream fin(fn);
  readHMAP (fin);
}

HMAPSequence::HMAPSequence (istream& in)
  : evd1_field (0.0f),
    evd2_field (0.0f),
    sse_string ("")

{
  readHMAP (in);
}

const string* HMAPSequence::getSSEString () const
{
  if (sse_string=="") buildSSEString();
  return & sse_string;
}

void HMAPSequence::readHMAP (istream& in)
{
  char buff[256];
 
  if (!in.good()) throw string ("Error reading file");
 
  in.getline(buff,256,':');
  if (strcmp(buff,"PDB")==0) {
    in.getline(buff,256); // ignore rest of line
    in.getline(buff,256,':');
  }
  
  if (strcmp(buff,"ID ")==0) in >> seq_name;
  else throw string ("Parse error before 'ID'");
  in.getline(buff,256); // ignore rest of line

  in.getline(buff,256,':');
  if (strcmp(buff,"DE ")==0) in >> de_field;
  else throw string ("Parse error before 'DE'");
  in.getline(buff,256); // ignore rest of line

  in.getline(buff,256,':');
  if (strcmp(buff,"SR ")==0) in >> sr_field;
  else throw string ("Parse error before 'SR'");
  in.getline(buff,256); // ignore rest of line

  in.getline(buff,256,':');
  if (strcmp(buff,"EVD")==0) in >> evd1_field >> evd2_field;
  else throw string ("Parse error before 'EVD'");
  in.getline(buff,256); // ignore rest of line
  
  in.getline(buff,256,':');
  if (strcmp(buff,"LEN")==0) in >> seq_length;
  else throw string ("Parse error before 'LEN'");
  in.getline(buff,256); // ignore rest of line
  reserve(seq_length+2);

  int index = 0;

  HMAPElem* head = new HMAPElem ();
  head->olc = SequenceElem::Head;
  head->index = index++;
  push_back (head);

  for (unsigned int i=0; i<seq_length; i++) {
    HMAPElem* node = new HMAPElem(in);
    node->index = index++;
    push_back (node);
  }
  
  HMAPElem* tail = new HMAPElem ();
  tail->olc   = SequenceElem::Tail;
  tail->index = index++;
  push_back (tail);
  
  head->gap_values = at(1)->gap_values;
  tail->gap_values = at(seq_length)->gap_values;

  in.getline(buff,256);
  if (strcmp(buff,"//")!=0) 
    throw string ("end of profile '//' not found"); 
}

void HMAPSequence::buildSSEString () const
{
  stringbuf buffer("");
  for (vector<HMAPElem*>::const_iterator it=begin();it!=end();++it) {
    char sse;
    float helix   = (*it)->p_helix();
    float strand  = (*it)->p_strand();
    float coil    = (*it)->p_coil();
    float confid  = (*it)->sse_confid;
    if ((*it)->isHead()) {
      sse = SequenceElem::Head;
    } else if ((*it)->isTail()) {
      sse = SequenceElem::Tail;
    } else if (helix>strand && helix>coil) {
      if (helix < .5 || confid < .5) sse = 'h';
      else sse = 'H';
    } else if (strand>helix && strand>coil) {
      if (strand < .5 || confid < .5) sse = 'e';
      else sse = 'E';
    } else {
      sse = ' ';
    }
    buffer.sputc (sse);
  }
  sse_string = buffer.str();
}

void HMAPSequence::getDefaultFlags (SuboptFlags& subopt_flags) 
{
  subopt_flags.Set(0,true);
  for(unsigned int i=1;i<=seq_length;++i) {
    if (at(i)->p_coil() > 0.3f)
      subopt_flags.Set(i,false);     // set loop regions to false
    else 
      subopt_flags.Set(i,true);
  }
  subopt_flags.Set(seq_length+1,true);
}

LogisticNormal::LogisticNormal (float qp, float qw,
				float tp, float tw,
				float eff)
  : q_peak (qp), q_width (qw),
    t_peak (tp), t_width (tw),
    eff_num(eff)  {}

float LogisticNormal::significance (float score) const
{
  float ev1 = oneSidedSignificance(score, t_peak, t_width);
  float ev2 = oneSidedSignificance(score, q_peak, q_width);
  if (ev1 >= 0 && ev2 >= 0) return sqrt(ev1*ev2);
  else if (ev1 >= 0) return ev1;
  else if (ev2 >= 0) return ev2;
  else return 9999.0f;
}

float LogisticNormal::oneSidedSignificance(float score, float peak, 
					   float width) const
{
  // clt 6.11.2002 
  // There are two transforms used below, one based on a normal pdf and
  // one based on the logistic function.  The current implementation is
  // based on the observation that the p-value follows a normal more
  // closely when z-score < 0 but follows the logistic when z-score > 0.
  // Note: the e-value for z-score > 0 behaves like an upper bound.

  // Also note: this e-value will only be accurate for global-local 
  // alignments (and probably global) but not local

  if (width<=0) return -1.0f;

  double zscore = (score-peak)/width;

  if (zscore<0) {
    // clt 5.30.2002 e-value computed from a normal distribution
#ifdef WIN32
    float pvalue = (float)(vs_erfc(zscore/1.41421356)/2.0);
#else
    float pvalue = (float)(erfc(zscore/1.41421356)/2.0);
#endif
    float evalue = (float)(eff_num) * pvalue;
    return evalue;
  } else {
    // lei's e-value corrected w/beta constant
    float pvalue = float(1.0/(exp(zscore*1.81379936)+1.0));
    float evalue = float(eff_num) * pvalue;
    return evalue;
  }
  
}

ostream& operator<< (ostream& os, HMAPElem& n)

{

int i;
char buffer[256];

sprintf(buffer,"%4i ",n.index);
os << buffer
   << n.olc << " ";

os.setf(ios::fixed);
for(i=0; i<20; i++) os << n.aa_profile[i]*100.0f << " ";
os << endl;

os << "   -   ";
for (i=0;i<4;i++) os << n.gap_values[i] << " ";

os << n.motif_value << " "
   << n.motif_confid << endl;

os << "   *   ";
for(i=0; i<3; i++) os << n.sse_values[i] << " ";
os << n.sse_confid << " "
   << n.surfacc_value << " "
   << n.surfacc_confid << endl;

return os;

}

ostream& operator<<(ostream& os,HMAPSequence& s)

{
unsigned int i;

os << "ID : " << s.seq_name << endl
   << "DE : " << s.de_field << endl
   << "SR : " << s.sr_field << endl
   << "EVD: " << s.evd1_field << " " << s.evd2_field << endl
   << "LEN: " << s.size()-2 << endl;

for(i=1;i<s.size()-1;i++) {
  s[i]->index=i;
  os << *s[i];
  }

os << "//\n";

return os;

}
