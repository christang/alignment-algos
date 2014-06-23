/**
 *  Package HMAP2.1
 *  File: gn2lib_seq.h
 *  Desc: Derived sequence class for SMAP profiles
 *
 *  Created on 11/10/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_GN2LIB_SEQ
#define _HMAP2_GN2LIB_SEQ

#include "alignment.h"
#include "hmapalib_seq.h"

enum {
  HM_FileNotFound
} ;

class SMAPSequence;
class Gn2Eval;

typedef AlignedPair<HMAPSequence,SMAPSequence> GnAlignedPair;
typedef list<GnAlignedPair> LGnAlignedPair;
typedef AlignedPairList<HMAPSequence,SMAPSequence> GnAlignment;
typedef AlignmentSet<HMAPSequence,SMAPSequence,Gn2Eval> GnAlignmentSet;
typedef vector<GnAlignment> VGnAlignment;

class SMAPSequence : public HMAPSequence
{

public:

  SMAPSequence (const char* fn, int v=0, bool f=false);
  SMAPSequence (istream& in, int v=0, bool f=false);
  ~SMAPSequence ();
  
  vector<vector<unsigned long> > brokenhb;
  vector<vector<float> > distance;
  vector<float> weighted_contact_number;
  vector<vector <unsigned long> > intra_hb_table;

  int verbose;
  bool gn2;                // 11.19.07 set true to make gn2 faster

  vector<vector<float> > distance2;         // 11.11.07 not used in gn2
  vector<vector<float> > angle;             // 11.11.07 not used in gn2

  void readSMAP(istream&);
  friend ostream& operator<< (ostream& os, SMAPSequence& s);

  bool get_backbone_HB_contact( int, int ) const;
  Structure *getStructure() { return &structure; }

  void updateCore (const GnAlignmentSet& as, float ratio);

private:

  void calcStructProperties ();
  void calcHBondContactMap ();
  void calcBroken_HBs ();
  void calcCBDistance2 ();
  void calcPrimaryDistances ();
  void calcSecondaryDistances ();
  void calcWeightedContactNumber ();
  void calcAccessibility ();
  void calcSSAngles ();

  Structure structure;
  string pdb_id,pdb_chain;

  bool** HBondContactMap;
  vector<HBond> Backbone_HBonds;
  float** cb_dist2;
};

#endif  // _HMAP2_GN2LIB_SEQ
