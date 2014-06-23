
#include "hmap_eval.h"

const float HMAPaliParams::default_alpha         = 0.5f;
const float HMAPaliParams::default_beta          = 1.0f;
const float HMAPaliParams::default_gamma         = 0.1f;
const float HMAPaliParams::default_zero_shift    = 0.12f;
const bool  HMAPaliParams::default_normalize_mtx = true;

HMAPaliParams::HMAPaliParams ()
  : AliParams(),
    NOaliParams (),
    alpha         (default_alpha),
    beta          (default_beta),
    gamma         (default_gamma),
    normalize_mtx (default_normalize_mtx),
    zero_shift    (default_zero_shift)
{
  // Empty
}

void HMAPaliParams::read (ParamStore* p)
{
  string s;

  s="CORE_MATCH_WEIGHT"; if (p->find(s)) p->getValue(s) >> alpha;
  s="CORE_GAP_WEIGHT"; if (p->find(s)) p->getValue(s) >> beta;
  s="MOTIF_MATCH_WEIGHT"; if (p->find(s)) p->getValue(s) >> gamma;
  s="NORMALIZE_SIM_MTX"; if (p->find(s)) p->getValue(s) >> normalize_mtx;
  s="ZERO_SHIFT"; if (p->find(s)) p->getValue(s) >> zero_shift;
  
  this->NOaliParams::read(p);
  this->AliParams::read(p);
}

HMAPaliEval::HMAPaliEval (HMAPaliParams& p) : params(&p) {}

void HMAPaliEval::pre_calculate (const HMAPSequence& s1, const HMAPSequence& s2) const
{
  for (unsigned int i=0; i<s2.size(); ++i) {
    float Pi = exp(params->beta * (1.f - 1.25f * s2[i]->p_coil()));
    s2[i]->gap_init(params->gap_init_penalty * Pi);
    s2[i]->gap_extn(params->gap_extn_penalty * Pi);
  }
}

void HMAPaliEval::post_process (SimilarityMatrix& s) const
{
  norm_elements (s,s,1,s.rows()-1,1,s.cols()-1);
  shift_elements (s,s,1,s.rows()-1,1,s.cols()-1,-params->zero_shift);
}
