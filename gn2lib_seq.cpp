/**
 *  Package HMAP2.1
 *  File: gn2lib_seq.cpp
 *  Desc: Derived sequence class for GNOALI profiles
 *
 *  Created on 11/10/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include "surfvsurface.h"
#include "gn2lib_seq.h"
#include <sstream>

using namespace std;

SMAPSequence::SMAPSequence (const char* fn, int v, bool f)
  : verbose (v), gn2 (f)
{
  ifstream fin(fn);
  readSMAP(fin);
}


SMAPSequence::SMAPSequence (istream& in, int v, bool f)
  : verbose (v), gn2 (f), HBondContactMap (0), cb_dist2 (0)
{
  readSMAP(in);
}

SMAPSequence::~SMAPSequence ()
{
  unsigned int nr ( structure.residue.size() );

  // delete the inter-residue square distance matrix
  if (cb_dist2!=0) {
    for(unsigned int i=0; i<nr; ++i)
      delete[] cb_dist2[i];
    delete[] cb_dist2;
  }
  
  // delete the hbond contact map matrix
  if (HBondContactMap!=0) {
    for(unsigned int i=0; i<nr+1; ++i)
      delete[] HBondContactMap[i];
    delete[] HBondContactMap;
  }
}

void SMAPSequence::readSMAP (istream& in)

{

  char buff[256];
 
  in.getline(buff,256,':');
  if (strcmp(buff,"PDB")) throw string ("SMAP file before 'PDB'");

  in >> pdb_id >> pdb_chain;
  in.getline(buff,256); // ignore rest of line

  // Create and read pdb.
  try {
  PDBFile pdbfile((char *)pdb_id.c_str());
  string subdef("chain "+pdb_chain);
  Subset chsub((char *)subdef.c_str());
  pdbfile.selection=chsub; // Read only necessary chain.
  pdbfile >> structure;
  } catch (...) { throw string("Can't read PDB file."); }

  if( (int)structure.chain.size()!=1 ) {
    string emsg("Chain ");
    stringstream ss;
    ss << structure.chain.size();
    emsg+=pdb_chain;
    emsg+=" (size = ";
    emsg+=ss.str();
    emsg+=")";
    emsg+=" is ambiguous.  Exiting.\n";
    throw emsg;
  }

  readHMAP (in);

  //cerr << "readSMAP: about to calcStructProperties" << endl;

  calcStructProperties ();
  if(seq_length!=structure.residue.size()) 
    throw string("Length of profile and length of PDB file do not match.");

  //cerr << "readSMAP: end" << endl;
}

void SMAPSequence::calcStructProperties ()
{
  Chain *ch=structure.chain[0];
  unsigned int nr=ch->residue.size();

  Residue **r=&ch->residue[0];
  for (unsigned int i=0;i<nr;++i) r[i]->param=i+1;

  // Reassign lods type for gn2 sequence.

  for (unsigned int i=0;i<nr;++i) {
    at(i+1)->lods_type = 0;
    if (at(i+1)->sse_values[1]>.5f) at(i+1)->lods_type = 1;
    if (at(i+1)->sse_values[2]>.5f) at(i+1)->lods_type = 2;
  }

  // Assign sse data.
  
  SSE **sse=&ch->ssel[0];
  unsigned int nsse=ch->ssel.size();

  for(unsigned int i=0;i<nsse;++i) {
    Residue **res=&sse[i]->res[0];
    for (unsigned int j=0;j<sse[i]->res.size();++j) {
      at(res[j]->param)->rdata.isse = i;
      at(res[j]->param)->rdata.sse_type = sse[i]->type;
      at(res[j]->param)->rdata.prev_sse = sse[i];
      at(res[j]->param)->rdata.next_sse = sse[i];
    }
  }

  for(unsigned int i=0;i<nr;++i) { 
    // For coiled regions:
    if (at(i+1)->rdata.isse==-1) {
      // Find the SSE of the first residue before i that has SSE
      for (unsigned int j=i;j>0;--j) {
	if (at(j)->rdata.isse!=-1) {
	  at(i+1)->rdata.prev_sse=sse[at(j)->rdata.isse];
	  break;
	}
      }
      // Ditto, for the SSE of the first residue after i.
      for (unsigned int j=i+1;j<nr;++j) {
	if (at(j+1)->rdata.isse!=-1) {
	  at(i+1)->rdata.next_sse=sse[at(j+1)->rdata.isse];
	  break;
	}
      }
    }
    //cout << at(i+1)->olc << i+1 << " " << at(i+1)->rdata.prev_sse << " " << at(i+1)->rdata.next_sse << "" << endl;
  }  

  // Assign residue-specific data 

  for(unsigned int i=0;i<nr;++i) {

    //cerr << "calcStructProperties: for loop: i: " << i << ", nr: " << nr << endl;
    //cerr << "current residue type: " << at(i+1)->olc << endl;

    try {
      at(i+1)->rdata.n=r[i]->atom["N"]->pos;
      at(i+1)->rdata.ca=r[i]->atom["CA"]->pos;
      at(i+1)->rdata.c=r[i]->atom["C"]->pos;
    } catch (...) {
      // there must be at least one atom in residue...
      at(i+1)->rdata.n=r[i]->atom.front()->pos;
      at(i+1)->rdata.ca=r[i]->atom.front()->pos;
      at(i+1)->rdata.c=r[i]->atom.front()->pos;

      cerr << "***missing atoms***\nresidue: " << at(i+1)->olc << i+1 ;
      cerr << ", atoms in residue: " << r[i]->atom.size() << endl;
    }

    try {
      at(i+1)->rdata.cb=r[i]->atom["CB"]->pos;
    } catch (...) {
      if( at(i+1)->olc != 'G' )
	cerr << "residue: " << at(i+1)->olc << i+1 << ", CB missing" << endl;
      at(i+1)->rdata.cb=at(i+1)->rdata.ca;
    }

  }

  Backbone_HBonds = structure.GetHBonds (aabackbone,aabackbone);
  calcHBondContactMap ();
  calcBroken_HBs ();
  calcPrimaryDistances ();
  calcWeightedContactNumber ();
  if (!gn2) calcAccessibility ();       // only calculate these for gnoali
  if (!gn2) calcSecondaryDistances ();
  if (!gn2) calcSSAngles ();

  at(0)->rdata.n=at(1)->rdata.n;
  at(0)->rdata.ca=at(1)->rdata.ca;
  at(0)->rdata.cb=at(1)->rdata.cb;
  at(0)->rdata.c=at(1)->rdata.c;
  at(0)->rdata.accessibility=at(1)->rdata.accessibility;

  at(nr+1)->rdata.n=at(nr)->rdata.n;
  at(nr+1)->rdata.ca=at(nr)->rdata.ca;
  at(nr+1)->rdata.cb=at(nr)->rdata.cb;
  at(nr+1)->rdata.c=at(nr)->rdata.c;
  at(nr+1)->rdata.accessibility=at(nr)->rdata.accessibility;

  //cerr << "calcStructProperties: end" << endl;
}

void SMAPSequence::calcAccessibility ()
{
  // Assign residue solvent accessibility

  SurfvSurface stotal(structure.atom,2);

  if (verbose>1) {
    cout << "Residue solvent accessibility." << endl;
  }
  unsigned int nr = structure.residue.size();
  Chain *ch=structure.chain[0];
  Residue **r=&ch->residue[0];

  for (unsigned int i=0; i<nr; ++i) {
    float sa(0.0),ra(0.0);
    ResidueDescription *rd=r[i]->description;
    Atom **ap=&r[i]->atom[0];

    unsigned int na=r[i]->atom.size();
    for(unsigned int j=0; j<na; ++j) {
      sa+=(float)stotal.Surface(ap[j]->pos,ap[j]->radius+
                                Troll::Application::probe_radius);
      AtomDescription *ad=rd->FindAtom(ap[j]->name);
      if(ad->refa<0.0) {
        cerr << "Undefined reference area for " << *ap[j] << endl;
        continue;
      }
      ra+=ad->refa;
    }

    sa/=ra;

    at(i+1)->rdata.accessibility = sa<1.f?sa:1.f;

    if (verbose>1) {
      cout << at(i+1)->olc << " " << at(i+1)->rdata.accessibility << endl;
    }
  }
}

// Setup the matrix of C-beta squared distances, cb_dists**

void SMAPSequence::calcCBDistance2 ()
{
  unsigned int nr( structure.residue.size() );

  cb_dist2 = new float* [ nr ];

  for( unsigned int i=0; i<nr; i++ ) {
    cb_dist2[i] = new float [ nr ];
  }

  for( unsigned int i=0; i<nr; i++ ) {
    cb_dist2[i][i] = 0;
    for( unsigned int j=0; j<i; j++ ) {
      cb_dist2[i][j] = ( at(i+1)->rdata.cb - at(j+1)->rdata.cb ).norm();
      cb_dist2[i][j]*= cb_dist2[i][j];  // square it
      cb_dist2[j][i] = cb_dist2[i][j];  // make it symmetric
    }
  }
}

void SMAPSequence::calcWeightedContactNumber () 
{
  unsigned int ao ( 0 );
  unsigned int nr ( structure.residue.size() );
  float span = ao * 2 + 1;
  
  weighted_contact_number.resize(nr+2,0.f);
  weighted_contact_number[0]=0.f;
  weighted_contact_number[nr+1]=0.f;

  calcCBDistance2();

  for (unsigned int i=0; i<nr; ++i) {
    for (unsigned int j=0; j<nr; ++j) {
      if (cb_dist2[i][j] > 14.5f && cb_dist2[i][j] < 256.f) 
	for (unsigned int z=max(0u,i-ao); z<=min(nr-1,i+ao); ++z)
	  weighted_contact_number[z+1] += ( 0.722f / cb_dist2[i][j] ) / span;
    }
  }
  
  //for (unsigned int i=1; i<=nr; ++i) cerr << "WCN: " << at(i)->olc << i << " " << weighted_contact_number[i] << endl;

}

void SMAPSequence::updateCore (const GnAlignmentSet& as,float ratio)
{

  unsigned int ao ( 0 );
  unsigned int nr( structure.residue.size() );
  float span ( ao * 2 + 1 );
  float len ( (float) as.size() );

  vector<float> model_cn(nr, 0.f);

  // work out model_cn for each alignment

  for (VGnAlignment::const_iterator al_it=as.begin();
       al_it!=as.end();  al_it++) {
    vector<bool> occupancy (nr+2, false);
    for (LGnAlignedPair::const_iterator ap_it=al_it->begin();
	 ap_it!=al_it->end();  ap_it++)
      occupancy [ap_it->template_idx()] = true;
    //for (unsigned int i=0; i<nr+2; ++i) cerr << occupancy [i]; cerr << endl;
    for (unsigned int i=0; i<nr; ++i) 
      for (unsigned int j=0; j<nr; ++j)
	if (cb_dist2[i][j] > 14.5f && cb_dist2[i][j] < 256.f && occupancy[j+1])
	  model_cn [i] += ( 0.722f / cb_dist2[i][j] ) / len;
  }
  
  // update core values from model_cn
 
  for (unsigned int i=1; i<=nr; ++i) cerr << "OLD: " << at(i)->olc << i << " " << weighted_contact_number[i] << "\t" << model_cn [i-1] << endl;

  for (unsigned int i=1; i<=nr; ++i) {
    weighted_contact_number[i] *= ratio;
    for (unsigned int z=max(1u,i-ao); z<=min(nr,i+ao); ++z)
      weighted_contact_number[z] += (1.f - ratio) * model_cn[i-1] / span;
  }

  for (unsigned int i=1; i<=nr; ++i) cerr << "NEW: " << at(i)->olc << i << " " << weighted_contact_number[i] << endl;

}

void SMAPSequence::calcHBondContactMap ()
{
  unsigned int nr ( structure.residue.size() );

  HBondContactMap = new bool* [ nr+1 ];

  for( unsigned int i=0; i<=nr; ++i ) {
    HBondContactMap[i] = new bool [ i+1 ];
  }

  for( unsigned int i=0; i<=nr; ++i ) {
    for( unsigned int j=0; j<=i; ++j ) {
      HBondContactMap[i][j] = false;
    }
  }

  for( unsigned int i=0; i<nr; i++ ) {
    structure.residue[i]->param = i+1;
  }

  for( unsigned int i=0; i<Backbone_HBonds.size(); ++i ) {
    int r1 ( Backbone_HBonds[i].donor.n->residue->param );
    int r2 ( Backbone_HBonds[i].acceptor.n->residue->param );

    //    cerr << "r1: " << r1 << "\t";
    //    cerr << "r2: " << r2 << endl;

    if( r1 >= r2 ) {
      HBondContactMap[r1][r2] = true;
    }
    else {
      HBondContactMap[r2][r1] = true;
    }

  }

}

bool SMAPSequence::get_backbone_HB_contact( int i, int j ) const {

  int nr( (int) structure.residue.size() );

  //  cerr << "get_backbone_HB_contact: nr: " << nr << endl;

  if( ( i >= (nr+1) ) || ( j >= (nr+1) ) ) {
    cerr << "Index out of bounds when checking backbone H-bond contact. Exiting." << endl;
    exit(-1);
  }

  if( i>=j ) {
    return HBondContactMap[i][j];
  }
  else {
    return HBondContactMap[j][i];
  }

}


void SMAPSequence::calcBroken_HBs ()
{
  unsigned int i, j, nr ( structure.residue.size() );
  
  HBond *hbp=&Backbone_HBonds[0];
  unsigned int nhb ( Backbone_HBonds.size() );

  // Each value hb_table[a][b] is the number of
  // bb hbs between residue a and b.
  vector<vector <unsigned long> > hb_table(nr);
  for(i=0;i<nr;++i) hb_table[i].resize(nr,0);
  for(i=0;i<nhb;++i) {
    int r1( hbp[i].donor.n->residue->param-1 );
    int r2( hbp[i].acceptor.n->residue->param-1 );
    if(r1==r2) continue;
    hb_table[r1][r2]=1;
    hb_table[r2][r1]=1;
  }

  vector<unsigned long> hb_table_row_sum(nr,0);
  for(i=0;i<nr;++i) {
    for(j=0;j<nr;++j) {
      hb_table_row_sum[i] += hb_table[i][j];
    }
  }

  // Each value intra_hb_table[a][b] is the number of 
  // bb hbs internal to residues a to b
  intra_hb_table.resize(nr);
  for(i=0;i<nr;++i) intra_hb_table[i].resize(nr,0);
  for(i=1;i<nr;++i) intra_hb_table[i][i-1] = 2*hb_table[i][i-1];
  for(i=2;i<nr;++i) {
    for (int j=i-2;j>=0;--j) {
      intra_hb_table[i][j] = intra_hb_table[i-1][j] + intra_hb_table[i][j+1] 
	- intra_hb_table[i-1][j+1] + 2*hb_table[i][j];
    }
  }

  // Each value brokenbh[a][b] is the number of bb hbs
  // between deleted by removing residues a to b, minus
  // the hydrogen bonds internal to residues a to b.
  brokenhb.resize(nr,vector<unsigned long>(nr,0));
  for(i=0;i<nr;++i) brokenhb[i][i] = hb_table_row_sum[i];
  for(i=1;i<nr;++i) {
    for(int j=i-1;j>=0;--j) {
      brokenhb[i][j] = brokenhb[i-1][j] + brokenhb[i][j+1] 
	- brokenhb[i-1][j+1];
    }
  }

  for(i=1;i<nr;++i) {
    for(int j=i-1;j>=0;--j) {
      brokenhb[i][j] -= intra_hb_table[i][j];
    }
  }

  // Output values if verbose
  if (verbose>2) {
    cout << "Per Residue Intra HBs" << endl;
    for(i=0;i<nr;++i) {
      HMAPElem* t1 = at(i+1);
      cout << t1->olc << i+1 << ":";
      for (j=0;j<i;++j) cout << " " << intra_hb_table[i][j] + brokenhb[i][j];
      cout << endl;
    }
  }

  if (verbose>2) {
    cout << "Intra HBs" << endl;
    for(i=0;i<nr;++i) {
      HMAPElem* t1 = at(i+1);
      cout << t1->olc << i+1 << ":";
      for (j=0;j<i+1;++j) cout << " " << intra_hb_table[i][j];
      cout << endl;
    }
  }

  if (verbose>1) {
    cout << "Broken HBs" << endl;
    for(i=0;i<nr;++i) {
      HMAPElem* t1 = at(i+1);
      cout << t1->olc << i+1 << ":";
      for (j=0;j<i+1;++j) cout << " " << brokenhb[i][j];
      cout << endl;
    }
  }
}


void SMAPSequence::calcPrimaryDistances ()
{
  if (verbose>1) cout << "Distance of residues" << endl;
  distance.resize(seq_length);
  for (unsigned int i=2; i<seq_length+2; ++i) {
    distance[i-2].resize(i-1);
    HMAPElem* t2 = at(i);
    if (verbose>1) cout << t2->rdata.sse_type << t2->olc << i << ":";
    for (unsigned int j=0; j<i-1; ++j) {
      HMAPElem* t1 = at(j);
      //double rd=(t2->rdata.n-t1->rdata.c).norm();
      double rd=(t2->rdata.cb-t1->rdata.cb).norm();
      distance[i-2][j] = rd;
      if (verbose>1) cout << " " << rd;
    }
    if (verbose>1) cout << endl;
  }
} 

void SMAPSequence::calcSecondaryDistances ()
{
  if (verbose>1) cout << "Distance of residues once removed" << endl;
  distance2.resize(seq_length);
  for (unsigned int i=2; i<seq_length+2; ++i) {
    distance2[i-2].resize(i-1);
    unsigned int ii=i<seq_length+1?i+1:i;
    HMAPElem* t2 = at(ii);
    if (verbose>1) cout << at(i)->rdata.sse_type << at(i)->olc << i << ":";
    for (unsigned int j=0; j<i-1; ++j) {
      unsigned int jj=j>0?j-1:j;
      HMAPElem* t1 = at(jj);
      double rd=(t2->rdata.n-t1->rdata.c).norm();
      distance2[i-2][j] = rd;
      if (verbose>1) cout << " " << rd;
    }
    if (verbose>1) cout << endl;
  }
}


void SMAPSequence::calcSSAngles ()
{
  if (verbose>1) cout << "Cosine of angle between SSEs." << endl;
  angle.resize(seq_length);
  for(unsigned int i=2; i<seq_length+2; ++i) {
    angle[i-2].resize(i-1,-1.f);
    HMAPElem* t2 = at(i);
    if (verbose>1) cout << t2->rdata.sse_type << t2->olc << i << ":";
    for (unsigned int j=0; j<i-1; ++j) {
      HMAPElem* t1 = at(j);
      if(t1->rdata.prev_sse && t2->rdata.next_sse) {
	Vector* temp = t1->rdata.prev_sse->Axis();
	Vector a = temp[1]-temp[0];
	temp = t2->rdata.next_sse->Axis();
	Vector b = temp[1]-temp[0];
	double ad;
	if (a.norm() == 0 || b.norm() == 0) ad=1.0;
	else ad=a*b/a.norm()/b.norm();
	angle[i-2][j]=ad;
      }
      if (verbose>1) cout << " " << angle[i-2][j];
    }
    if (verbose>1) cout << endl;
  }
}


ostream& operator<<(ostream& os,SMAPSequence& s)
{
  
  os << "PDB: " << s.pdb_id << " " << s.pdb_chain << endl;
  os << (HMAPSequence&)s;
  
  return os;
  
}
