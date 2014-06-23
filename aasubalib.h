
#include "alib.h"
#include "evaluator.h"
#include "sequence.h"
#include "submatrix.h"

template <class S1, class S2>
class AASubstitutionEval :
  public Evaluator<S1,S2,AASubstitutionEval<S1,S2> >
{
 
public:

  AASubstitutionEval (AliParams& p, SubstitutionMatrix& m)
    : params(&p), sub_matrix(&m) {};

  inline float similarity (const S1& q, const S2& t,
                           int q_pos, int t_pos) const
  {
    if (q[q_pos]->isHead()||q[q_pos]->isTail()||
	t[t_pos]->isHead()||t[t_pos]->isTail()) return 0.f;
    char aa1 = q[q_pos]->olc;
    char aa2 = t[t_pos]->olc;
    return sub_matrix->score (aa1,aa2);
  };

  inline float deletion   (const S1& q, const S2& t,
                           int q_pos1, int q_pos2,
			   int t_pos1, int t_pos2) const
  { 
    int len = t_pos2-t_pos1-1;
    if (len < 1) return 0.f;
    else {
      switch (params->align_type) {
      case global:
      case global_local:
	return params->gap_init_penalty + 
	  params->gap_extn_penalty * (len-1);
      case local:
      case semi_local:
      case local_global:
	if (t[t_pos1]->isHead() || t[t_pos2]->isTail())
	  return 0;
	else
	  return params->gap_init_penalty + 
	    params->gap_extn_penalty * (len-1);
      default:
	throw string ("Illegal gap style");
      }
    }
  }
  
  inline float insertion  (const S1& q, const S2& t,
                           int q_pos1, int q_pos2,
			   int t_pos1, int t_pos2) const
  { 
    int len = q_pos2-q_pos1-1;
    if (len < 1) return 0.f;
    else {
      switch (params->align_type) {
      case global:
      case local_global:
	return params->gap_init_penalty + 
	  params->gap_extn_penalty * (len-1);
      case local:
      case semi_local:
      case global_local:
	if (q[q_pos1]->isHead() || q[q_pos2]->isTail())
	  return 0;
	else
	  return params->gap_init_penalty + 
	    params->gap_extn_penalty * (len-1);
      default:
	throw string ("Illegal gap style");
      }
    }
  }
 
  inline void pre_calculate (const S1& q, const S2& t) const {}
  inline void post_process (SimilarityMatrix& s) const {}

private:
  
  AliParams* params;
  SubstitutionMatrix* sub_matrix;

};
