/**
 *  Package HMAP2.1
 *  File: alignment.h
 *  Desc: Storage template classes for aligned sequences. 
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_ALIGNMENT
#define _HMAP2_ALIGNMENT

#include <algorithm>
#include <assert.h>
#include <list>
#include <utility>
#include <vector>
#include "math.h"

#include "dpmatrix.h"
#include "enumerator.h"
#include "sequence.h"
#include "significance.h"
#include "sflags.h"

#include <iostream>

using namespace std;

//----------------------------------------------------------------------------
typedef pair<int,int> _AlignedPair;

template <class S1, class S2>
class AlignedPair : public _AlignedPair
{
  typedef _AlignedPair Base;
public:
  AlignedPair():Base(-1,-1) {};
  AlignedPair(int i,int j):Base(i,j) {};
  AlignedPair(const AlignedPair& a):Base(a) {};
  inline int query_idx() const { return first; }
  inline int template_idx() const { return second; }
};

//----------------------------------------------------------------------------
template <class S1, class S2>
class AlignedPairList : public list<AlignedPair<S1,S2> > 
{
  typedef AlignedPair<S1,S2> Pair;
  typedef list<Pair> Base;
public:
 AlignedPairList() 
   : Base(0),score(0.f),identity(0.f),significance(9999.f),uid(-1) {};

 AlignedPairList( const Pair& a, float s )
   : score(s),identity(0.f),significance(9999.f),uid(-1) { push_back(a); };

 AlignedPairList( const AlignedPairList& a ) 
   : Base(a),score(a.score),identity(a.identity),
    significance(a.significance),SSE_CO(a.SSE_CO),coverage(a.coverage),uid(a.uid) {};

  AlignedPairList( const S1&, const S2& );

  AlignedPairList& operator=( const AlignedPairList& rhs ); 

  void append(int i,int j);
  void prepend(int i,int j);
  void calcIdentity (const string& a, const string& b);
  void readFrom (const string& q, const string& t);

  //AK
  void print_frag( string, string ) const;
  string get_templ_string( string ) const;
  string get_query_string( string ) const;
  void print_pairs() const;
  void remove_first_pair();
  void remove_last_pair();
  void remove_ends();
  void combine ( AlignedPairList& );
  int get_last_query_idx();
  int get_first_query_idx();
  int get_last_template_idx();
  int get_first_template_idx();
  bool frag_follows( AlignedPairList& );
  int get_shift ( AlignedPairList&, string&, SuboptFlags&, int *len=0 );
  float get_simple_shift ( AlignedPairList&, SuboptFlags&, int *len=0 );
  void get_q_all ( AlignedPairList&, SuboptFlags&, int *n_agree,
		   float *q_mod=0, float *q_dev=0, float *q_comb=0);
  float get_area_diff( const AlignedPairList& ) const;
  void compare_segments( bool&, bool&,
			 const float&, const float&, const float&, const float&,
			 const float&, const float&, const float&, const float&,
			 float&, float&,
			 float&, float&, float&, float&, float&, float&, float&, float& ) const;
  float get_trapezoid_area( float, float, float, float, float, float, float, float ) const;
  void fix_zigzag();
  float find_perpendicular_distance( int, int, int, int );

  template <class Stype>
  void calcSignificance (const Significance<Stype>& s);
  inline bool operator< (const AlignedPairList& a) const 
    { return score > a.score; } ;
  
  float score;
  float identity;
  float significance;
  float SSE_CO; //AK
  float coverage; // AK
  int uid;
};

template <class S1, class S2>
void AlignedPairList<S1,S2>::readFrom (const string& query, 
				       const string& templ) 
{
  score = 0.f;
  identity = 0.f;
  significance = 9999.f;
  uid = -1;

  Base::clear();

  int seq1(-1), seq2(-1);
  unsigned int i = 0;

  if (query.size()!=templ.size()) 
    throw ("readFrom error: query and templ not equal length");
  
  float aligned = 0.f;

  while( i < query.size() ) {

    if( query [ i ] != '-' ) { seq1++; }
    if( templ [ i ] != '-' ) { seq2++; }

    if( ( query [ i ] != '-' ) && ( templ [ i ] != '-' ) ) {
      append( seq1, seq2 );
      
      if ( query [ i ] != '^' && query [ i ] != '$' &&
	   templ [ i ] != '^' && templ [ i ] != '$')
	aligned += 1.f;
      if ( query [ i ] == templ [i] ) 
	if ( query [ i ] != '^' && query [ i ] != '$' &&
	     templ [ i ] != '^' && templ [ i ] != '$')
	  identity += 1.f;
    }
    
    i++;
  }

  identity /= aligned;
  identity *= 100.f;
}

//AK
template <class S1, class S2>
AlignedPairList<S1,S2>::AlignedPairList( const S1& s1, const S2& s2 )
{
  score = 0.f;
  identity = 0.f;
  significance = 9999.f;
  uid = -1;

  int seq1(-1), seq2(-1);
  unsigned int i = 0;
  // s1 and s2 are assumed to be of the same length

  while( i < s1.size() ) {

    if( s1.olc( i ) != '-' ) { seq1++; }
    if( s2.olc( i ) != '-' ) { seq2++; }

    if( ( s1.olc( i ) != '-' ) && ( s2.olc( i ) != '-' ) ) {
      append( seq2, seq1 );
    }

    i++;
  }

}

template <class S1, class S2>
AlignedPairList<S1,S2>& AlignedPairList<S1,S2>::operator= (const AlignedPairList<S1,S2>& rhs)
{
  if (this == &rhs) return *this;
  
  Base::operator=(rhs);
  score = rhs.score;
  identity = rhs.identity;
  significance = rhs.significance;
  uid = rhs.uid;  

  return *this;
}

template <class S1, class S2>
void AlignedPairList<S1,S2>::print_frag( string query, string templ ) const
{

  cerr << get_templ_string( templ ) << "\t" << "(" << Base::front().template_idx() << "-" << Base::back().template_idx() << ")" << endl;
  cerr << get_query_string( query ) << "\t" << "(" << Base::front().query_idx() << "-" << Base::back().query_idx() << ")";

}



template <class S1, class S2>
string AlignedPairList<S1,S2>::get_templ_string( string templ_str ) const
{
  string res;

  // print template, with gaps where needed
  typename list<Pair>::const_iterator pair_it = list<Pair>::begin();
  typename list<Pair>::const_iterator pair_end = list<Pair>::end();

  //  cerr << templ_str[pair_it->template_idx()];
  //  res.append( templ_str[pair_it->template_idx()] ); // print the first template residue
  res.append( templ_str.substr( pair_it->template_idx(), 1 ) ); // print the first template residue
                                                    // first pair of a fragment always matched

  const Pair* prev = &*pair_it; ++pair_it;
  
  for( ; pair_it!=pair_end; ++pair_it) {
    if( pair_it->template_idx() == prev->template_idx()+1 ) { // no insertion in template
                                                                  // possible insertion in query
      int x = pair_it->query_idx();
      int y = x - prev->query_idx();
      int gap = y-1;
      for( int k=0; k<gap; k++ ) {
	res.append( "-" );
	//	cerr << "-"; // print the gap character
      }
    }
    else { // insertion in template (or zigzag)

      // print the inserted template residues
      int templ_res_id = prev->template_idx()+1;
      while( templ_res_id < pair_it->template_idx() ) {
	res.append( templ_str.substr( templ_res_id, 1 ) );
	//	cerr << templ_str[templ_res_id];
	templ_res_id++;
      }

      if( pair_it->query_idx() > prev->query_idx()+1 ) { // zigzag case
	// print gaps to match the query residues      

	int query_res_id = prev->query_idx()+1;
	while( query_res_id < pair_it->query_idx() ) {
	  res.append( "-" );
	  //	  cerr << "-";
	  query_res_id++;
	}

      }

    }
    
    res.append( templ_str.substr( pair_it->template_idx(), 1 ) );
    //    cerr << templ_str[pair_it->template_idx()];
    
    prev = &*pair_it;
  }

  //  cerr << "\t" << "(" << front().template_idx() << "-" << back().template_idx() << ")" << endl;

  return res;

}


template <class S1, class S2>
string AlignedPairList<S1,S2>::get_query_string( string query_str ) const
{

  string res;

  // print query, with gaps where needed
  typename list<Pair>::const_iterator pair_it = list<Pair>::begin();
  typename list<Pair>::const_iterator pair_end = list<Pair>::end();

  //  cerr << query_str[pair_it->query_idx()];
  //  res.append( query_str[pair_it->query_idx()] ); // print the first query residue
  res.append( query_str.substr( pair_it->query_idx(), 1 ) ); // print the first query residue
                                                    // first pair of a fragment always matched

  const Pair* prev = &*pair_it; ++pair_it;
  
  for( ; pair_it!=pair_end; ++pair_it) {
    if( pair_it->query_idx() == prev->query_idx()+1 ) { // no insertion in query
                                                                  // possible insertion in query
      int x = pair_it->template_idx();
      int y = x - prev->template_idx();
      int gap = y-1;
      for( int k=0; k<gap; k++ ) {
	res.append( "-" );
	//	cerr << "-"; // print the gap character
      }
    }
    else { // insertion in query (or zigzag)

      // print the inserted query residues
      int query_res_id = prev->query_idx()+1;
      while( query_res_id < pair_it->query_idx() ) {
	res.append( query_str.substr( query_res_id, 1 ) );
	//	cerr << query_str[query_res_id];
	query_res_id++;
      }

      if( pair_it->template_idx() > prev->template_idx()+1 ) { // zigzag case
	// print gaps to match the template residues      

	int templ_res_id = prev->template_idx()+1;
	while( templ_res_id < pair_it->template_idx() ) {
	  res.append( "-" );
	  //	  cerr << "-";
	  templ_res_id++;
	}

      }

    }
    
    res.append( query_str.substr( pair_it->query_idx(), 1 ) );
    //    cerr << query_str[pair_it->query_idx()];
    
    prev = &*pair_it;
  }

  //  cerr << "\t" << "(" << front().query_idx() << "-" << back().query_idx() << ")" << endl;

  return res;


}


template <class S1, class S2>
void AlignedPairList<S1,S2>::get_q_all ( AlignedPairList& native, 
					 SuboptFlags& core, 
					 int *n_agree, float *q_mod, 
					 float *q_dev, float *q_comb )
{ 
  assert (n_agree != 0);
  if (get_last_query_idx() != (int)core.size()-1) 
    throw string ("Core file length does not match alignment");
  
  *n_agree = -2;  // account for head and tail
  typename list<Pair>::const_iterator nat_it = native.begin();
  for ( typename list<Pair>::const_iterator cur_it = list<Pair>::begin(); 
	cur_it!=list<Pair>::end(); ) {
    while (nat_it->query_idx() < cur_it->query_idx() 
	   && nat_it != native.end()) ++nat_it;

    while (cur_it->query_idx() < nat_it->query_idx()
	   && cur_it != list<Pair>::end()) ++cur_it;
    
    if (nat_it == native.end() || cur_it == list<Pair>::end())
      break;

    if (nat_it->query_idx() != cur_it->query_idx())
      continue;

    if (core[cur_it->query_idx()])
      if (nat_it->template_idx() == cur_it->template_idx())
	(*n_agree)++;
    
    nat_it++; cur_it++;
  }

  vector<bool> seen(core.size(),false);

  int d_mod = -2;  // account for head and tail, and count core positions
  for ( typename list<Pair>::const_iterator cur_it = list<Pair>::begin(); 
	cur_it!=list<Pair>::end(); cur_it++ )
    if (core[cur_it->query_idx()]) 
      { ++d_mod; seen[cur_it->query_idx()]=true; }
  
  int d_dev = -2;
  for ( typename list<Pair>::const_iterator nat_it = native.begin(); 
	nat_it!=native.end(); nat_it++ )
    if (core[nat_it->query_idx()]) 
      { ++d_dev; seen[nat_it->query_idx()]=true; }
  
  int d_comb = -2;
  for ( typename vector<bool>::const_iterator it = seen.begin();
	it!=seen.end(); it++ )
    if ( *it ) ++d_comb;

  if (q_mod != 0)     *q_mod  = float(*n_agree) / float(d_mod);

  if (q_dev != 0)     *q_dev  = float(*n_agree) / float(d_dev);

  if (q_comb != 0)    *q_comb = float(*n_agree) / float(d_comb);
}


template <class S1, class S2>
float AlignedPairList<S1,S2>::get_simple_shift ( AlignedPairList<S1,S2>& apl, 
						 SuboptFlags& core, int *len)
{
  if (get_last_query_idx() != (int)core.size()-1) 
    throw string ("Core file length does not match alignment");

  //cerr << "# start" << endl;

  int al(0), ts(0);

  typename list<Pair>::const_iterator idx(apl.begin());

  //cerr << idx->query_idx() << " - " << idx->template_idx() << endl;

  for ( typename list<Pair>::const_iterator it = list<Pair>::begin(); 
	it!=list<Pair>::end(); it++ ) {
    while ( idx->query_idx()<it->query_idx() ) { 
      idx++; 
      if (idx==apl.end()) break;
    }
    if (idx==apl.end()) break;
    if ( idx->query_idx() == it->query_idx() &&
	 core[ it->query_idx() ] ) {
	ts += abs (idx->template_idx()-it->template_idx());
	al ++;
    }
    //cerr << it->query_idx() << " : " << ts << " @ " << al << endl;
  }

  // Return length of alignment to *len

  if (len!=0) *len = al;
  if (!al) throw string ("No residues aligned");
  return float(ts)/float(al);

}


//AK
template <class S1, class S2>
void AlignedPairList<S1,S2>::print_pairs() const
{
  for ( typename list<Pair>::const_iterator it = list<Pair>::begin(); 
	it!=list<Pair>::end(); it++ ) {
    cerr << it->template_idx() << "\t" << it->query_idx() << endl;
  }
  cerr << endl;
}

//AK
template <class S1, class S2>
void AlignedPairList<S1,S2>::remove_first_pair()
{
  erase( list<Pair>::begin() );
}

//AK
template <class S1, class S2>
void AlignedPairList<S1,S2>::remove_last_pair()
{
  erase( --list<Pair>::end() );
}


//AK
template <class S1, class S2>
void AlignedPairList<S1,S2>::remove_ends()
{
  erase( list<Pair>::begin() );
  erase( --list<Pair>::end() );
}

//AK
template <class S1, class S2>
void AlignedPairList<S1,S2>::combine ( AlignedPairList& a )
{
  splice( list<Pair>::end(), a );
  score += a.score;
}

//AK
template <class S1, class S2>
int AlignedPairList<S1,S2>::get_last_query_idx()
{
  typename list<Pair>::iterator it = list<Pair>::end();
  it--;

  return it->query_idx();
}

//AK
template <class S1, class S2>
int AlignedPairList<S1,S2>::get_first_query_idx()
{
  return list<Pair>::begin()->query_idx();
}

//AK
template <class S1, class S2>
int AlignedPairList<S1,S2>::get_last_template_idx()
{
  typename list<Pair>::iterator it = list<Pair>::end();
  it--;

  return it->template_idx();
}

//AK
template <class S1, class S2>
int AlignedPairList<S1,S2>::get_first_template_idx()
{
  return list<Pair>::begin()->template_idx();
}

//AK
template <class S1, class S2>
bool AlignedPairList<S1,S2>::frag_follows( AlignedPairList& a )
{
  // the '+1' is to ensure that the loop between the sse frags ...
  // is at least 1 residue long
  return( get_last_query_idx()+1 < a.get_first_query_idx() );
}

//AK
template <class S1, class S2>
float AlignedPairList<S1,S2>::get_area_diff( const AlignedPairList& a2 ) const
{

  //  cerr << "entering get_area_diff" << endl;

  typename list<Pair>::const_iterator it1 = Base::begin(); 
  const Pair* prev1 = &*it1;
  it1++;

  typename list<Pair>::const_iterator it2 = a2.begin(); 
  const Pair* prev2 = &*it2;
  it2++;

  float xp, yp, xa1, xa2, ya1, ya2, xb1, xb2, yb1, yb2;

  float x_ali_a1, y_ali_a1, x_ali_a2, y_ali_a2;
  float x_ali_b1, y_ali_b1, x_ali_b2, y_ali_b2;

  bool exists, main_is_former;
  bool has_area = false;
  float area_diff = 0;

  const Pair *former, *latter, *former_prev, *latter_prev;

  //  cerr << "last temple res main: " << back().template_idx() << endl;
  //  cerr << "last temple res other: " << a2.back().template_idx() << endl;

  while( ( it1 != Base::end() ) || ( it2 != a2.end() ) ) {

    //    cerr << "it1 templ res: " << it1->template_idx() << endl;
    //    cerr << "it2 templ res: " << it2->template_idx() << endl;

    if( it1->template_idx() <= it2->template_idx() ) {
      //      cerr << "top of if" << endl;
      main_is_former = true;
      former = &*it1;
      former_prev = prev1;
      latter = &*it2;
      latter_prev = prev2;
      prev1 = &*it1;
      it1++;
    }
    else {
      main_is_former = false;
      former = &*it2;
      former_prev = prev2;
      latter = &*it1;
      latter_prev = prev1;
      prev2 = &*it2;
      it2++;
    }

    xa1 = former_prev->template_idx();   ya1 = former_prev->query_idx();
    xa2 = former->template_idx();        ya2 = former->query_idx();
    xb1 = latter_prev->template_idx();   yb1 = latter_prev->query_idx();
    xb2 = latter->template_idx();        yb2 = latter->query_idx();

    compare_segments( exists, has_area,
		      xa1, ya1, xa2, ya2,
		      xb1, yb1, xb2, yb2,
		      xp, yp,
		      x_ali_a1, y_ali_a1, x_ali_a2, y_ali_a2,
		      x_ali_b1, y_ali_b1, x_ali_b2, y_ali_b2 );

    if( has_area ) {

      if( !exists ) { // no intersection inside these segments

	area_diff += abs( get_trapezoid_area( x_ali_a1, y_ali_a1, x_ali_a2, y_ali_a2,
					      x_ali_a1, Base::back().query_idx(),
					      x_ali_a2, Base::back().query_idx() ) -

			  get_trapezoid_area( x_ali_b1, y_ali_b1, x_ali_b2, y_ali_b2,
					      x_ali_b1, Base::back().query_idx(),
					      x_ali_b2, Base::back().query_idx() ) );
      }
      else { // there is an intersection inside these segments

	// calc area before intersection
	area_diff += abs( get_trapezoid_area( x_ali_a1, y_ali_a1, xp, yp,
					      x_ali_a1, Base::back().query_idx(),
					      xp,       Base::back().query_idx() ) -

			  get_trapezoid_area( x_ali_b1, y_ali_b1, xp, yp,
					      x_ali_b1, Base::back().query_idx(),
					      xp,       Base::back().query_idx() ) );

	// calc area after intersection
	area_diff += abs( get_trapezoid_area( xp, yp, x_ali_a2, y_ali_a2,
					      xp,       Base::back().query_idx(),
					      x_ali_a2, Base::back().query_idx() ) -

			  get_trapezoid_area( xp, yp, x_ali_b2, y_ali_b2,
					      xp,       Base::back().query_idx(),
					      x_ali_b2, Base::back().query_idx() ) );

      }
    }

    if( xa2 == xb2 ) { // handle the case where the segments end at the same x-coordinate
      if( main_is_former ) {
	prev2 = &*it2;
	it2++;
      }
      else {
	prev1 = &*it1;
	it1++;
      }
    }

  }
  
  //  cerr << "about to return area" << endl;

  return area_diff;
}

template <class S1, class S2>
void AlignedPairList<S1,S2>::compare_segments( bool& exists, bool& has_area,
					       const float& xa1, const float& ya1,
					       const float& xa2, const float& ya2,
					       const float& xb1, const float& yb1,
					       const float& xb2, const float& yb2,
					       float& xp, float& yp,
					       float& x_ali_a1, float& y_ali_a1,
					       float& x_ali_a2, float& y_ali_a2,
					       float& x_ali_b1, float& y_ali_b1,
					       float& x_ali_b2, float& y_ali_b2 ) const
{

  bool same_p1, same_p2;

  // check for same start points
  same_p1 = (xa1 == xb1) && (ya1 == yb1);
  same_p2 = (xa2 == xb2) && (ya2 == yb2);

  if( same_p1 && same_p2 ) { // identical segments
    exists = true;
    has_area = false;

    x_ali_a1 = xa1;   y_ali_a1 = ya1;
    x_ali_a2 = xa2;   y_ali_a2 = ya2;
    x_ali_b1 = xb1;   y_ali_b1 = yb1;
    x_ali_b2 = xb2;   y_ali_b2 = yb2;

    return;
  }

  float x_inter_min = xa1 > xb1 ? xa1 : xb1;
  float x_inter_max = xa2 < xb2 ? xa2 : xb2;

  float m_a = ( ya2 - ya1 ) / ( xa2 - xa1 );
  float m_b = ( yb2 - yb1 ) / ( xb2 - xb1 );

  float intercept_a = ya1 - ( m_a * xa1 );
  float intercept_b = yb1 - ( m_b * xb1 );

  if( same_p1 && !same_p2 ) { // segments start together, but end separately
    exists = true;
    has_area = m_a != m_b;

    xp = xa1;
    yp = ya1;

    x_ali_a1 = xa1;           y_ali_a1 = ya1;
    x_ali_a2 = x_inter_max;   y_ali_a2 = ( m_a * x_inter_max ) + intercept_a;
    x_ali_b1 = xb1;           y_ali_b1 = yb1;
    x_ali_b2 = x_inter_max;   y_ali_b2 = ( m_b * x_inter_max ) + intercept_b;

    return;
  }

  if( !same_p1 && same_p2 ) { // segments start together, but end separately
    exists = true;
    has_area = m_a != m_b;

    xp = xa2;
    yp = ya2;

    x_ali_a1 = x_inter_min;   y_ali_a1 = ( m_a * x_inter_min ) + intercept_a;
    x_ali_a2 = xa2;           y_ali_a2 = ya2;
    x_ali_b1 = x_inter_min;   y_ali_b1 = ( m_b * x_inter_min ) + intercept_b;
    x_ali_b2 = xb2;           y_ali_b2 = yb2;

    return;
  }

  // segments don't start or end together
  if( m_a == m_b ) {
    if( intercept_a == intercept_b ) { // segments line on same line
      exists = true;
      has_area = false;

      x_ali_a1 = x_inter_min;   y_ali_a1 = ( m_a * x_inter_min ) + intercept_a;
      x_ali_a2 = x_inter_max;   y_ali_a2 = ( m_a * x_inter_max ) + intercept_a;
      x_ali_b1 = x_inter_min;   y_ali_b1 = ( m_b * x_inter_min ) + intercept_b;
      x_ali_b2 = x_inter_max;   y_ali_b2 = ( m_b * x_inter_max ) + intercept_b;

      return;
    }
    else { // intercept_a != intercept_b --> segments lie on parallel lines
      exists = false;
      has_area = true;

      x_ali_a1 = x_inter_min;   y_ali_a1 = ( m_a * x_inter_min ) + intercept_a;
      x_ali_a2 = x_inter_max;   y_ali_a2 = ( m_a * x_inter_max ) + intercept_a;
      x_ali_b1 = x_inter_min;   y_ali_b1 = ( m_b * x_inter_min ) + intercept_b;
      x_ali_b2 = x_inter_max;   y_ali_b2 = ( m_b * x_inter_max ) + intercept_b;

      return;
    }
  }
  else { // slopes not equal, so find intersection point
    xp = ( intercept_a - intercept_b ) / ( m_b - m_a );

    if( x_inter_min <= xp && xp <= x_inter_max ) { // interesection occurs in both segments
      exists = true;
      has_area = true;
      yp = intercept_a + ( m_a * xp );

      x_ali_a1 = x_inter_min;   y_ali_a1 = ( m_a * x_inter_min ) + intercept_a;
      x_ali_a2 = x_inter_max;   y_ali_a2 = ( m_a * x_inter_max ) + intercept_a;
      x_ali_b1 = x_inter_min;   y_ali_b1 = ( m_b * x_inter_min ) + intercept_b;
      x_ali_b2 = x_inter_max;   y_ali_b2 = ( m_b * x_inter_max ) + intercept_b;

      return;
    }
    else { // intersection occurs outside one or both segments
      exists = false;
      has_area = true;

      x_ali_a1 = x_inter_min;   y_ali_a1 = ( m_a * x_inter_min ) + intercept_a;
      x_ali_a2 = x_inter_max;   y_ali_a2 = ( m_a * x_inter_max ) + intercept_a;
      x_ali_b1 = x_inter_min;   y_ali_b1 = ( m_b * x_inter_min ) + intercept_b;
      x_ali_b2 = x_inter_max;   y_ali_b2 = ( m_b * x_inter_max ) + intercept_b;

      return;
    }


  }

}

template <class S1, class S2>
float AlignedPairList<S1,S2>::get_trapezoid_area( float x_ali1, float y_ali1,
						  float x_ali2, float y_ali2, 
						  float x_base1, float y_base1,
						  float x_base2, float y_base2 ) const
{

  return( ( ( ( y_base1 - y_ali1 ) + ( y_base2 - y_ali2 ) ) / 2 ) * ( x_ali2 - x_ali1 ) );

}


template <class S1, class S2>
void AlignedPairList<S1,S2>::fix_zigzag()
{
  typename list<Pair>::iterator it = Base::begin();
  Pair* prev = &*it;
  it++;

  // find each zigzag in the alignment
  for ( ; it!=Base::end(); it++ ) {
    if( ( it->template_idx() - prev->template_idx() > 1 ) &&
	( it->query_idx() - prev->query_idx() > 1 ) ) { // found a zigzag

      int q_beg = prev->query_idx();
      int t_beg = prev->template_idx();
      int q_end = it->query_idx();
      int t_end = it->template_idx();

      //      cerr << "found zigzag at: (" << q_beg << "," << t_beg << "), (" << q_end << "," << t_end << ")" << endl;

      int q_new = q_beg;
      int t_new = t_beg;

      while( ( ( q_end - q_new ) > 1 ) && ( ( t_end - t_new ) > 1 ) ) {

	q_new++; t_new++;

	while( find_perpendicular_distance( q_end-q_beg, t_end-t_beg, (q_new+1)-q_beg, t_new-t_beg ) <
	       find_perpendicular_distance( q_end-q_beg, t_end-t_beg, q_new-q_beg,     t_new-t_beg ) ) {

	  q_new++;
	  //	  cerr << "new point: (" << q_new << "," << t_new << ")" << endl; 
	}

	while( find_perpendicular_distance( q_end-q_beg, t_end-t_beg, q_new-q_beg, (t_new+1)-t_beg ) <
	       find_perpendicular_distance( q_end-q_beg, t_end-t_beg, q_new-q_beg,     t_new-t_beg ) ) {

	  t_new++;
	  //	  cerr << "new point: (" << q_new << "," << t_new << ")" << endl; 
	}

	// at this point, either q_new or t_new has changed (but not both) to the point closest
	// to the segment joining (q_beg, t_beg) and (q_end, t_end)
	insert( it, Pair( q_new, t_new ) );
      }

    }
    prev = &*it;
  }
}


template <class S1, class S2>
float AlignedPairList<S1,S2>::find_perpendicular_distance( int x1p, int y1p, int xp, int yp )
{

  float dist_a_sq = x1p*x1p + y1p*y1p;
  float dist_b_sq = xp*xp + yp*yp;

  float cos_theta_sq = ( ( (x1p*xp) + (y1p*yp) ) * ( (x1p*xp) + (y1p*yp) ) ) / ( dist_a_sq * dist_b_sq );
  float sin_theta_sq = 1 - cos_theta_sq;

  return( sqrt( dist_b_sq * sin_theta_sq ) );
}


template <class S1, class S2>
inline void AlignedPairList<S1,S2>::append(int i,int j)
{ AlignedPair<S1,S2> p(i,j); this->push_back(p); };

template <class S1, class S2>
inline void AlignedPairList<S1,S2>::prepend(int i,int j)
{ AlignedPair<S1,S2> p(i,j); this->push_front(p); };

template <class S1, class S2>
void AlignedPairList<S1,S2>::calcIdentity(const string& query,
					  const string& templ)
{
  int same = - 2;  // adjust for head and tail elements
  int total = (int) min (query.size(),templ.size()) - 2;
  for (typename list<Pair>::iterator it=list<Pair>::begin(); 
       it!=list<Pair>::end(); ++it)
    if (query[it->query_idx()] == templ[it->template_idx()]) ++same;
  identity = float(same)/float(total) * 100.f;
}

template <class S1, class S2>
template <class Stype>
void AlignedPairList<S1,S2>::
calcSignificance(const Significance<Stype>& s) 
{
  significance = s.significance(score);
}

//----------------------------------------------------------------------------
template <class S1, class S2, class Etype>
class AlignmentSet : public vector<AlignedPairList<S1,S2> >
{
  typedef AlignedPairList<S1,S2> Alignment;
  typedef vector<Alignment> Base;
public:
  AlignmentSet (const Alignment& apl):Base(1,apl) {};
  AlignmentSet (const AlignmentSet& as);
  AlignmentSet (DPMatrix<S1,S2,Etype>& dpm, 
		Enumerator<S1,S2,Etype>& en);

  inline const S1*    getQuerySequence () const 
    { return dpmatrix->getQuerySequence(); } ;
  inline const S2* getTemplateSequence () const 
    { return dpmatrix->getTemplateSequence(); } ;
  inline const DPMatrix<S1,S2,Etype>* getDPMatrix () const
    { return dpmatrix; }

  void sortSet (int max);
  
  void assignIdentity ();
  
  template <class Stype>
  void assignSignificance (const Significance<Stype>& s);

private:
  void build ();
  DPMatrix<S1,S2,Etype>* dpmatrix;
  Enumerator<S1,S2,Etype>* enumerator;
};

template <class S1, class S2, class Etype>
AlignmentSet<S1,S2,Etype>::AlignmentSet (const AlignmentSet<S1,S2,Etype>& as)
  : Base (as),
    dpmatrix (&as.dpmatrix),
    enumerator (&as.enumerator) {}

template <class S1, class S2, class Etype>
AlignmentSet<S1,S2,Etype>::AlignmentSet (DPMatrix<S1,S2,Etype>& dpm, 
					 Enumerator<S1,S2,Etype>& en)
  : dpmatrix (&dpm),
    enumerator (&en)
{
  build ();
}

template <class S1, class S2, class Etype>
void AlignmentSet<S1,S2,Etype>::sortSet (int max)
{
  if (max>=(int)vector<Alignment>::size()) {
    sort (vector<Alignment>::begin(),vector<Alignment>::end());
  } else if (max>0) {
    typename vector<Alignment>::iterator it=vector<Alignment>::begin();
    partial_sort (it,it+max,vector<Alignment>::end());
    this->erase (it+max,vector<Alignment>::end());
  }
}

template <class S1, class S2, class Etype>
void AlignmentSet<S1,S2,Etype>::build ()
{
  this->reserve(enumerator->estimateSize());
  enumerator->enumerate(*dpmatrix,*this);
  assignIdentity();
}

template <class S1, class S2, class Etype>
void AlignmentSet<S1,S2,Etype>::assignIdentity ()
{
  for (typename vector<Alignment>::iterator it=vector<Alignment>::begin(); 
       it!=vector<Alignment>::end(); ++it)
    it->calcIdentity(*dpmatrix->getQuerySequence()->getString(),
		     *dpmatrix->getTemplateSequence()->getString());
}

template <class S1, class S2, class Etype>
template <class Stype>
void AlignmentSet<S1,S2,Etype>::
assignSignificance (const Significance<Stype>& s)
{
  for (typename vector<Alignment>::iterator it=vector<Alignment>::begin(); 
       it!=vector<Alignment>::end(); ++it)
    it->calcSignificance (s);
}


//----------------------------------------------------------------------------

struct IdentityComparator {
  template <class S1, class S2>
  inline bool operator() (const AlignedPairList<S1,S2>& a,
			  const AlignedPairList<S1,S2>& b)
  { return a.identity > b.identity; };
};

struct ScoreComparator {
  template <class S1, class S2>
  inline bool operator() (const AlignedPairList<S1,S2>& a,
			  const AlignedPairList<S1,S2>& b)
  { return a.score > b.score; };
};

#endif  // _HMAP2_ALIGNMENT
