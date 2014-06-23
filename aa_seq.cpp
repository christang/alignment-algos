
#include "aa_seq.h"

AASequence::AASequence() : Base() {}

AASequence::~AASequence() 
{
  for (vector<SequenceElem*>::reverse_iterator it=rbegin(); it!=rend(); ++it)
    delete *it;
}

void AASequence::append(const string& s)
{
  unsigned int index = size();
  for (string::const_iterator it=s.begin(); it!=s.end(); ++it) 
    push_back (new SequenceElem(index++,*it));
}

void AASequence::append(const char* cs)
{
  string s(cs);
  append(s);
}

void AASequence::cleargaps(char c)
{
  for (vector<SequenceElem*>::iterator it=vector<SequenceElem*>::begin();
       it!=vector<SequenceElem*>::end(); ) {
    if ((*it)->olc == c) vector<SequenceElem*>::erase(it);
    else ++it;
  }
  seq_string="";
}

void AASequence::calcPrimaryDistances ()
{
  //if (verbose>1) cout << "Distance of residues" << endl;
  seq_length = size()-2;
  distance.resize(seq_length);
  for (unsigned int i=2; i<seq_length+2; ++i) {
    distance[i-2].resize(i-1);
    //HMAPElem* t2 = at(i);
    //if (verbose>1) cout << t2->residue->getSecondary() << t2->olc << i << ":";
    for (unsigned int j=0; j<i-1; ++j) {
      //HMAPElem* t1 = at(j);
      //float rd = t1->residue->distanceTo(*t2->residue);
      distance[i-2][j] = 0; //rd;
      //if (verbose>1) cout << " " << rd;
    }
    //if (verbose>1) cout << endl;
  }
}
