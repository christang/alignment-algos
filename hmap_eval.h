
#include "alib.h"
#include "argv.h"
#include "evaluator.h"
#include "hmapalib_seq.h"
#include "hmath.h"
#include "noalib.h"
#include "rcfile.h"
#include "simmatrix.h"

#ifndef _HMAP2_HMAPALIB
#define _HMAP2_HMAPALIB

class HMAPaliParams :
  public AliParams,
  public NOaliParams  
{

public:
  
  HMAPaliParams ();

  void read (ParamStore* p);

  static const float default_alpha;
  static const float default_beta ;
  static const float default_gamma;
  static const bool  default_normalize_mtx;
  static const float default_zero_shift;

  float alpha;
  float beta;
  float gamma;
  bool  normalize_mtx;
  float zero_shift;

};

class HMAPaliEval :
  public Evaluator<HMAPSequence,HMAPSequence,HMAPaliEval> 
{

public:

  HMAPaliEval (HMAPaliParams& p);

  inline float similarity (const HMAPSequence& q,
			   const HMAPSequence& t,
			   int q_pos, int t_pos) const 
  {
    float ip = dot_product  (q[q_pos]->aa_profile,t[t_pos]->aa_profile);
    float pc = pearson_corr (q[q_pos]->sse_values,t[t_pos]->sse_values);
    
cout << q[q_pos]->sse_values[0] << " " << q[q_pos]->sse_values[1] << " " << q[q_pos]->sse_values[2] << " " <<
     t[t_pos]->sse_values[0] << " " << t[t_pos]->sse_values[1] << " " << t[t_pos]->sse_values[2] << " " << pc << endl; 
    float sim =
      ip * exp (params->alpha * pc *
		q[q_pos]->sse_confid *
		t[t_pos]->sse_confid);
    return sim;
  }
  
  inline float deletion   (const HMAPSequence& q,
			   const HMAPSequence& t,
			   int q_pos1, int q_pos2, 
			   int t_pos1, int t_pos2) const
  {
    int dist = t_pos2 - t_pos1;
    if (dist < 2) return 0;

    float gi, ge;
    gi = min(t[t_pos1]->gap_init(),t[t_pos2]->gap_init());
    ge = min(t[t_pos1]->gap_extn(),t[t_pos2]->gap_extn());

    switch (params->align_type) {
    case global:                   // overhangs penalized
    case global_local:             // overhangs penalized in template
      return gi + ge*(dist-2);
    case local:                    // overhangs not penalized
    case semi_local:               // overhangs not penalized
    case local_global:             // overhangs penalized in query
      if (t[t_pos1]->isHead() || t[t_pos2]->isTail())
	return 0;
      else
	return gi + ge*(dist-2);
    default:
      throw string ("Illegal gap style");
    }
  }
  
  inline float insertion  (const HMAPSequence& q,
			   const HMAPSequence& t,
			   int q_pos1, int q_pos2, 
			   int t_pos1, int t_pos2) const
  {
    int dist = q_pos2 - q_pos1;
    if (dist < 2) return 0;

    float gi, ge;
    gi = min(t[t_pos1]->gap_init(),t[t_pos2]->gap_init());
    ge = min(t[t_pos1]->gap_extn(),t[t_pos2]->gap_extn());

    switch (params->align_type) {
    case global:                   // overhangs penalized
    case local_global:             // overhangs penalized in query
      return gi + ge*(dist-2);
    case local:                    // overhangs not penalized
    case semi_local:               // overhangs not penalized
    case global_local:             // overhangs penalized in template
      if (q[q_pos1]->isHead() || q[q_pos2]->isTail())
	return 0;
      else
	return gi + ge*(dist-2);
    default:
      throw string ("Illegal gap style");
    }
  }

  void pre_calculate (const HMAPSequence& s1, const HMAPSequence& s2) const;
  void post_process (SimilarityMatrix& s) const;

private:
  
  HMAPaliParams* params;

};

#endif  // HMAP2_HMAPALIB
