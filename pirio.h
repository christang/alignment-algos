
#ifndef _HMAP2_PIRIO
#define _HMAP2_PIRIO

#include "formats.h"
#include "sequence.h"

//----------------------------------------------------------------------------

class PIRWrite
{
public:
  PIRWrite (ostream& o, int len);

  template <class S1, class S2, class Etype>
  void write (AlignmentSet<S1,S2,Etype>& as);
  
  void write (const string& s);

  void fix_ends (string& s);

  ostream* output;
  int line_length;
};

PIRWrite operator<< (ostream& o, Formats::PIROut p);

template <class S1, class S2, class Etype>
ostream& operator<< (PIRWrite w, AlignmentSet<S1,S2,Etype>& as)
{ w.write (as); return *w.output; };

template <class S1, class S2, class Etype>
void PIRWrite::write (AlignmentSet<S1,S2,Etype>& as)
{
  int count = 0;
  string gapped_string; 
  valarray<bool> mask (false,as.size());
  for (typename AlignmentSet<S1,S2,Etype>::iterator it=as.begin(); 
       it!=as.end(); ++it) {
    mask[count]=true;
    SequenceGaps gaps(as,mask);

    *output << "#start" << endl << endl;
    *output << ">P1;";
    *output << as.getTemplateSequence()->seq_name << endl;

    *output << "structureN:";
    *output << as.getTemplateSequence()->seq_name;
    *output << "::::" << endl;

    gaps.build (*as.getTemplateSequence()->getString(),gapped_string);
    fix_ends (gapped_string);
    write (gapped_string);

    *output << endl;
    *output << ">P1;";
    *output << as.getQuerySequence()->seq_name << endl;

    *output << "sequence:";
    *output << as.getQuerySequence()->seq_name;
    *output << "::::" << endl;

    gaps.build (*as.getQuerySequence()->getString(),*it,gapped_string);
    fix_ends (gapped_string);
    write (gapped_string);
    
    *output << endl;
    *output << "#end" << endl;
    mask[count++]=false;
  }
}

//----------------------------------------------------------------------------

class PIRRead {

public:
  PIRRead (istream& i, bool flag=true);

  template <class S>
  void readInto (S& apl);

  istream* input;
  bool head_tail;

};
  
PIRRead operator>> (istream& i, Formats::PIRIn p);

template <class S>
istream& operator>> (PIRRead r, S& s)
{ r.readInto (s); return *r.input; }

template <class S>
void PIRRead::readInto (S& apl)
{
  // S must be an AlignedPairlist

  string line;
  string templ_tmp_ali, query_tmp_ali;

  // Key up to #start
  
  while( line.find( "#start", 0 ) == string::npos ) {
    getline( *input, line );
    if( input->eof() ) { throw string ("Error (1) parsing PIR"); }
  }
  
  if( input->eof() ) { throw string ("Error (2) parsing PIR"); }

  // line now contains "#start..."

  while( line.find( "structure", 0 ) == string::npos ) {
    getline( *input, line );
  }
  
  getline( *input, line );
  // line now contains the first line of the template sequence
  
  while( true ) { // copy the template sequence lines
    templ_tmp_ali.append( line );
    
    if( line.length() == 0 || line[ line.length() - 1 ] == '*' ) { // check for the last line of the sequence
      break;
    }
    
    getline( *input, line );
  }
  
  // line now contains query sequence header line

  while( line.find( "sequence", 0 ) == string::npos ) {
    getline( *input, line );
  }

  getline( *input, line );
  // line now contains the first line of the query sequence

  while( true ) { // copy the query sequence lines
    query_tmp_ali.append( line );
    
    if( line.length() == 0 ||  line[ line.length() - 1 ] == '*' ) { // check for the last line of the sequence
      break;
    }
    
    getline( *input, line );
  }

  // remove * from end of seqs and add fasta characters
  if( templ_tmp_ali[ templ_tmp_ali.length() - 1 ] == '*' ) {
    templ_tmp_ali.resize( templ_tmp_ali.length() - 1 );
  }

  if( query_tmp_ali[ query_tmp_ali.length() - 1 ] == '*' ) {
    query_tmp_ali.resize( query_tmp_ali.length() - 1 );
  }

  if ( head_tail ) {

    if( *templ_tmp_ali.begin() != '^' ) { 
      templ_tmp_ali.insert( templ_tmp_ali.begin(), '^' ); }
    if( templ_tmp_ali[ templ_tmp_ali.length() - 1 ] != '$' ) {
      templ_tmp_ali.insert( templ_tmp_ali.end(), '$' );
    }
    
    if( *query_tmp_ali.begin() != '^' ) { 
      query_tmp_ali.insert( query_tmp_ali.begin(), '^' ); }
    if( query_tmp_ali[ query_tmp_ali.length() - 1 ] != '$' ) {
      query_tmp_ali.insert( query_tmp_ali.end(), '$' );
    }

  }
  
  apl.readFrom (query_tmp_ali, templ_tmp_ali);

}

#endif  //_HMAP2_PIRIO
