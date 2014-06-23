#ifndef _SSSS_SHARED
#define _SSSS_SHARED

struct SSE_Data {  // holds basic info on template SSEs
  int sse_id;
  int ss_type; // 329 = helix, 330 = strand
  int beg_id;
  int end_id;
};

struct Frag_ID {
  int sse_idx;
  int frag_idx;
};

struct Aligned_Pair {
  int templ_res;
  int query_res;
};

struct Frag_Connection {
  Frag_ID prev_frag;
  Frag_ID next_frag;
  int prev_end_res_idx;
  int next_beg_res_idx;
  float connection_score;
};

struct Res_Pair {
  float t;
  float q;
  int rel_pos;
};

#endif  //_SSSS_SHARED
