#ifndef AA_SEQUENCE
#define AA_SEQUENCE

#include <string>
#include <vector>

#include "sequence.h"

class AASequence : public Sequence<SequenceElem*> 
{
  typedef Sequence<SequenceElem*> Base;
public:
  AASequence ();
  //AASequence (const AASequence& s);
  //AASequence& operator= (const AASequence& s);
  ~AASequence ();
  void append (const string& s);
  void append (const char* s);
  void cleargaps (char c);

  vector<vector<float> > distance;
  void calcPrimaryDistances ();
};


#endif
