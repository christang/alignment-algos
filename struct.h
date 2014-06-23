/**
 *  Package HmapII
 *  File: structure.h
 *  Desc: Wrapper for Troll structure classes
 *
 *  Created on 6/29/2003
 *  Author: cltang @ honig lab
 *
 */

#ifndef _WSTRUCTURE
#define _WSTRUCTURE

#include <map>
#include <valarray>
#include <string>

// from Troll
#include "app.h"
#include "structure.h"
#include "pdbfile.h"
#include "troll.h"

using namespace std;
using namespace Troll;

namespace HMAP {

  enum geom_t {
    AA = 0,  // Calpha atom distance
    CN = 1  // C_i to N_i+1 atom distance
  };

  class HM_Data {

  public:
  HM_Data() : prev_sse(NULL),next_sse(NULL),isse(-1),sse_type(TC_Coil),accessibility(0.0) { };
    int deletedHBs (HM_Data& r, int off1=0, int off2=0);
    float cosSSEAxes (HM_Data& r);
    SSE* prev_sse;
    SSE* next_sse;
    int isse;
    unsigned long sse_type;
    double accessibility;
    Vector n,ca,c,cb;

  };

} // Namespace

#endif // _WSTRUCTURE

