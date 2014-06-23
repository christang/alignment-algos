/**
 *  Package HMAP2.1
 *  File: gnoalib_seq.cpp
 *  Desc: Derived sequence class for GNOALI profiles
 *
 *  Created on 11/10/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include "surfvsurface.h"
#include "gnoalib_seq.h"
#include <sstream>

using namespace std;

SMAPSequence::SMAPSequence (const char* fn, int v)
  : verbose (v)
{
  //  cerr << "SMAPSequence constructor: top" << endl;
  ifstream fin(fn);
  //  cerr << "SMAPSequence constructor: loaded " << fn << endl;
  //  cerr << "SMAPSequence constructor: about to read SMAP" << endl;
  readSMAP(fin);
  //  cerr << "SMAPSequence constructor: top" << endl;
}


SMAPSequence::SMAPSequence (istream& in, int v)
  : verbose (v)
{
  readSMAP(in);
}

/*
// AK
vector<vector<unsigned long> > SMAPSequence::export_backbone_hb_table() const {
}
*/

void SMAPSequence::readSMAP (istream& in)

{
  //  cerr << "readSMAP: top" << endl;

  char buff[256];
 
  in.getline(buff,256,':');
  if (strcmp(buff,"PDB")) throw string ("SMAP file before 'PDB'");

  in >> pdb_id >> pdb_chain;
  in.getline(buff,256); // ignore rest of line

  //  cerr << "readSMAP: about to read pdb" << endl;

  // Create and read pdb.
  PDBFile pdbfile((char *)pdb_id.c_str());
  string subdef("chain "+pdb_chain);
  Subset chsub((char *)subdef.c_str());
  pdbfile.selection=chsub; // Read only necessary chain.
  pdbfile >> structure;

  //  cerr << "readSMAP: created pdbfile" << endl;

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

  //  cerr << "readSMAP: about to readHMAP" << endl;

  readHMAP (in);

  //  cerr << "readSMAP: about to calcStructProperties" << endl;

  calcStructProperties ();
  if(seq_length!=structure.residue.size()) 
    throw string("Length of profile and length of PDB file do not match.");

  //  cerr << "readSMAP: end" << endl;

}


void SMAPSequence::calcStructProperties ()
{
  //  cerr << "calcStructProperties: top" << endl;

  int i,j;

  Chain *ch=structure.chain[0];
  unsigned int nr=ch->residue.size();

  Residue **r=&ch->residue[0];
  for (i=0;i<(int)nr;i++) r[i]->param=i+1;

  // Assign sse data.
  
  SSE **sse=&ch->ssel[0];
  unsigned int nsse=ch->ssel.size();

  //  cerr << "here! 1" << endl;

  for(i=0;i<(int)nsse;i++) {
    Residue **res=&sse[i]->res[0];
    for (j=0;j<(int)sse[i]->res.size();j++) {
      at(res[j]->param)->rdata.isse = i;
      at(res[j]->param)->rdata.sse_type = sse[i]->type;
      at(res[j]->param)->rdata.prev_sse = sse[i];
      at(res[j]->param)->rdata.next_sse = sse[i];
    }
  }

  //  cerr << "here! 2" << endl;

  for(i=0;i<(int)nr;i++) { 
    // For coiled regions:
    if (at(i+1)->rdata.isse==-1) {
      // Find the SSE of the first residue before i that has SSE
      for (j=i-1;j>=0;--j) {
	if (at(j+1)->rdata.isse!=-1) {
	  at(i+1)->rdata.prev_sse=sse[at(j+1)->rdata.isse];
	  break;
	}
      }
      // Ditto, for the SSE of the first residue after i.
      for (j=i+1;j<(int)nr;++j) {
	if (at(j+1)->rdata.isse!=-1) {
	  at(i+1)->rdata.next_sse=sse[at(j+1)->rdata.isse];
	  break;
	}
      }
    }
    //cout << at(i+1)->olc << i+1 << " " << at(i+1)->rdata.prev_sse << " " << at(i+1)->rdata.next_sse << "" << endl;
  }  

  //  cerr << "here! 3" << endl;

  // Assign residue specific data 

  SurfvSurface stotal(structure.atom,2);
  
  if (verbose>1) {
    cout << "Residue solvent accessibility." << endl;
  }

  for(i=0;i<(int)nr;i++) {

    //    cerr << "calcStructProperties: for loop: i: " << i << ", nr: " << nr << endl;
    //    cerr << "current residue type: " << at(i+1)->olc << endl;

    at(i+1)->rdata.n=r[i]->atom["N"]->pos;
    at(i+1)->rdata.ca=r[i]->atom["CA"]->pos;
    at(i+1)->rdata.c=r[i]->atom["C"]->pos;

    //    cerr << "stored positions in rdata fields" << endl;

    if( at(i+1)->olc == 'G' ) {
      //      cerr << "in if portion" << endl;
      at(i+1)->rdata.cb=r[i]->atom["CA"]->pos;    
    }
    else {
      //      cerr << "in else portion" << endl;
      at(i+1)->rdata.cb=r[i]->atom["CB"]->pos;    
    }

    //    cerr << "here! 4" << endl;

    float sa(0.0),ra(0.0);
    ResidueDescription *rd=r[i]->description;
    Atom **ap=&r[i]->atom[0];

    //    cerr << "here! 4.1" << endl;

    unsigned int na=r[i]->atom.size();
    for(j=0;j<(int)na;j++) {
      sa+=(float)stotal.Surface(ap[j]->pos,ap[j]->radius+
				Troll::Application::probe_radius);
      AtomDescription *ad=rd->FindAtom(ap[j]->name);
      if(ad->refa<0.0) {
	cerr << "Undefined reference area for " << *ap[j] << endl;
	continue;
      }
      ra+=ad->refa;
    }

    //    cerr << "here! 4.2" << endl;

    sa/=ra;

    //    cerr << "here! 4.3" << endl;

    at(i+1)->rdata.accessibility = sa<1.f?sa:1.f;

    //    cerr << "here! 4.4" << endl;

    if (verbose>1) {
      cout << at(i+1)->olc << " " << at(i+1)->rdata.accessibility << endl;
    }

    //    cerr << "here! 4.5" << endl;
  }

  //  cerr << "here! 5" << endl;

  at(0)->rdata.n=r[0]->atom["N"]->pos;
  at(0)->rdata.ca=r[0]->atom["CA"]->pos;
  at(0)->rdata.c=r[0]->atom["C"]->pos;
  at(0)->rdata.accessibility = at(1)->rdata.accessibility;

  at(nr+1)->rdata.n=r[nr-1]->atom["N"]->pos;
  at(nr+1)->rdata.ca=r[nr-1]->atom["CA"]->pos;
  at(nr+1)->rdata.c=r[nr-1]->atom["C"]->pos;
  at(nr+1)->rdata.accessibility = at(nr)->rdata.accessibility;

  //  cerr << "calcStructProperties: loaded n,ca, c positions" << endl;

  Backbone_HBonds = structure.GetHBonds(aabackbone,aabackbone);
  //  cerr << "calcStructProperties: got Backbone_HBonds" << endl;
  calcHBondContactMap();
  //  cerr << "calcStructProperties: did calcHBondContactMap" << endl;
  calcBroken_HBs ();
  calcPrimaryDistances ();
  calcSecondaryDistances ();
  calcSSAngles ();

  //  cerr << "calcStructProperties: end" << endl;
}


//AK
void SMAPSequence::calcHBondContactMap () {

  //  cerr << "calcHBondContactMap: top" << endl;

  int nr( (int)structure.residue.size() );

  //  cerr << "calcHBondContactMap: nr: " << nr << endl;

  HBondContactMap = new bool* [ nr+1 ];

  for( int i=0; i<(nr+1); i++ ) {
    HBondContactMap[i] = new bool [ i+1 ];
  }

  //  cerr << "calcHBondContactMap: reserved space for HBondContactMap**" << endl;

  for( int i=0; i<(nr+1); i++ ) {
    for( int j=0; j<=i; j++ ) {
      HBondContactMap[i][j] = false;
    }
  }

  //  cerr << "calcHBondContactMap: initialized HBondContactMap** to false" << endl;

  for( int i=0; i<nr; i++ ) {
    structure.residue[i]->param = i+1;
  }

  //  cerr << "calcHBondContactMap: reset param values in structure" << endl;

  for( int i=0; i<(int)Backbone_HBonds.size(); i++ ) {
    int r1=Backbone_HBonds[i].donor.n->residue->param;
    int r2=Backbone_HBonds[i].acceptor.n->residue->param;

    //    cerr << "r1: " << r1 << "\t";
    //    cerr << "r2: " << r2 << endl;

    if( r1 >= r2 ) {
      HBondContactMap[r1][r2] = true;
    }
    else {
      HBondContactMap[r2][r1] = true;
    }

  }

  //  cerr << "calcHBondContactMap: end" << endl;
}

bool SMAPSequence::get_backbone_HB_contact( int i, int j ) const {

  int nr((int)structure.residue.size());

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
  int i,j,nr((int)structure.residue.size());
  
  HBond *hbp=&Backbone_HBonds[0];
  unsigned int nhb=Backbone_HBonds.size();

  // Each value hb_table[a][b] is the number of
  // bb hbs between residue a and b.
  vector<vector <unsigned long> > hb_table(nr);
  for(i=0;i<nr;++i) hb_table[i].resize(nr,0);
  for(i=0;i<(int)nhb;++i) {
    int r1=hbp[i].donor.n->residue->param-1;
    int r2=hbp[i].acceptor.n->residue->param-1;
    if(r1==r2) continue;
    hb_table[r1][r2]++;
    hb_table[r2][r1]++;
  }

  vector<unsigned long> hb_table_row_sum(nr,0);
  for(i=0;i<nr;++i) {
    for(j=0;j<nr;++j) {
      hb_table_row_sum[i] += hb_table[i][j];
    }
  }

  // Each value intra_hb_table[a][b] is the number of 
  // bb hbs internal to residues a to b
  vector<vector <unsigned long> > intra_hb_table(nr);
  for(i=0;i<nr;++i) intra_hb_table[i].resize(nr,0);
  for(i=1;i<nr;++i) intra_hb_table[i][i-1] = 2*hb_table[i][i-1];
  for(i=2;i<nr;++i) {
    for (j=i-2;j>=0;--j) {
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
    for(j=i-1;j>=0;--j) {
      brokenhb[i][j] = brokenhb[i-1][j] + brokenhb[i][j+1] 
	- brokenhb[i-1][j+1];
    }
  }

  for(i=1;i<nr;++i) {
    for(j=i-1;j>=0;--j) {
      brokenhb[i][j] -= intra_hb_table[i][j];
    }
  }

  // Output values if verbose
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
      double rd=(t2->rdata.n-t1->rdata.c).norm();
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
