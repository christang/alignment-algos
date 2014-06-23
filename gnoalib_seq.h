/**
 *  Package HMAP2.1
 *  File: gnoalib_seq.h
 *  Desc: Derived sequence class for GNOALI profiles
 *
 *  Created on 11/10/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_GNOALIB_SEQ
#define _HMAP2_GNOALIB_SEQ

#include "hmapalib_seq.h"

enum {
  HM_FileNotFound
  } ;

class SMAPSequence : public HMAPSequence
{
  
public:

  SMAPSequence (const char* fn, int v=0);
  SMAPSequence (istream& in, int v=0);
  
  vector<vector<unsigned long> > brokenhb;
  vector<vector<float> > distance;
  vector<vector<float> > distance2;
  vector<vector<float> > angle;

  void readSMAP(istream&);
  friend ostream& operator<< (ostream& os, SMAPSequence& s);

  bool get_backbone_HB_contact( int, int ) const;
  Structure *getStructure() { return &structure; }

  int verbose;

private:

  void calcStructProperties ();
  void calcHBondContactMap ();
  void calcBroken_HBs ();
  void calcPrimaryDistances ();
  void calcSecondaryDistances ();
  void calcSSAngles ();

  Structure structure;
  string pdb_id,pdb_chain;
  bool** HBondContactMap;

  vector<HBond> Backbone_HBonds;

};

#endif  // _HMAP2_GNOALIB_SEQ
