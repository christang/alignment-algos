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

#ifndef _HMAP2_HMAPALIB_SEQ
#define _HMAP2_HMAPALIB_SEQ

#include <istream>
#include <string>
#include <valarray>
#include <vector>

#include "sequence.h"
#include "sflags.h"
#include "significance.h"
#include "struct.h"

using namespace HMAP;

class HMAPElem : public SequenceElem
{

public: 

  valarray<float> aa_profile;
  valarray<float> gap_values;
  float motif_value;
  float motif_confid;
  valarray<float> sse_values;
  float sse_confid;
  float surfacc_value;
  float surfacc_confid;
  HM_Data rdata;
  unsigned int lods_type;
  float hydropathy;


  inline float gap_init() const { return gap_values[0]; };
  inline float gap_extn() const { return gap_values[1]; };
  inline  void gap_init(float gi) { gap_values[0]=gi; };
  inline  void gap_extn(float ge) { gap_values[1]=ge; };
  
  inline float p_helix()  const { return sse_values[0]; };
  inline float p_strand() const { return sse_values[1]; };
  inline float p_coil()   const { return sse_values[2]; };
  
  HMAPElem();
  HMAPElem(istream& in);
  HMAPElem(const HMAPElem& e);
  
  void readHMAP (istream& in);
  void calcHydropathy ();

  HMAPElem& operator= (const HMAPElem& e);

 
};

ostream& operator<< (ostream& is, HMAPElem&);


class HMAPSequence : public Sequence<HMAPElem*>
{

public:

  HMAPSequence (const char* fn);
  HMAPSequence (istream& in);

  //HMAPSequence (const HMAPSequence& s);
  //HMAPSequence& operator= (const HMAPSequence& s);
  //~HMAPSequence ();

  string de_field;
  string sr_field;
  float  evd1_field;
  float  evd2_field;

  const string* getSSEString () const;
  friend ostream& operator<< (ostream& os, HMAPSequence& s);

  void getDefaultFlags( SuboptFlags& sof );

protected:
  
  HMAPSequence ();
  void readHMAP (istream& in);
  void buildSSEString () const;
  mutable string sse_string;

};

class LogisticNormal : public Significance<LogisticNormal> 
{

public:
 
  LogisticNormal (float qp, float qw, float tp, float tw, float eff=5000.f);
  
  float significance(float score) const; 
  
private:

  float oneSidedSignificance(float score, float peak, float width) const;

  float q_peak;
  float q_width;

  float t_peak;
  float t_width;

  float eff_num;

};

#endif  // _HMAP2_HMAPALIB_SEQ
