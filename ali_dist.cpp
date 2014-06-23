#include "ali_dist.h"
#include <string>
#include <iostream>

using namespace std;

/// PRIVATE FUNCTIONS ///

// convert different alignment formats to vector<Res_Pair>
void Ali_Dist::convert_strings_to_VRP( const string& templ, const string& query, vector<Res_Pair>& vrp )
{

  vrp.clear();
  
  if( templ.length() != query.length() ) {
    cerr << "Sequences are of unequal lengths." << endl;
    exit(-1);
  }
  
  int templ_idx(0), query_idx(0);
  
  for( int i=0; i<(int)templ.length(); i++ ) {
    if( ( templ[i] != '-' ) && ( query[i] != '-' ) ) {

      Res_Pair tmp;
      tmp.t = (float)templ_idx;
      tmp.q = (float)query_idx;
      tmp.rel_pos = -2; // default initial value

      vrp.push_back( tmp );
      templ_idx++;
      query_idx++;
    }
    else if( ( templ[i] != '-' ) && ( query[i] == '-' ) ) {
      templ_idx++;
    }
    else if( ( templ[i] == '-' ) && ( query[i] != '-' ) ) {
      query_idx++;
    }
  }

}


int Ali_Dist::get_sequence_length( const string& s )
{
  // find length of sequence in s
  int length(0);

  for( int i=0; i<(int)s.length(); i++ ) {
    if( ( s[i] != '-' ) && ( s[i] != '^' ) && ( s[i] != '$' ) ) {
      length++;
    }
  }

  return length;
}


int Ali_Dist::get_sequence_length( const vector<Res_Pair>& vrp, char s )
{
  // s is T(emplate) or Q(uery)

  if( s == 'T' ) {
    return( (int)vrp.back().t - 1 );
  }
  else if( s == 'Q' ) {
    return( (int)vrp.back().q - 1 );
  }
  else { // error checking
    cerr << s << " is an invalid value for the variable s.  Exiting." << endl;
    exit(-1);
  }

}


bool Ali_Dist::extract_next_ali( ifstream& file, string& templ, string& query )
{
  // extract the next ali from 'file' and store it in 'templ' and 'query'
  // return true if is one or more sequences left in 'file'
  // return false if file is at eof

  templ.clear(); // erase previous sequence
  query.clear();

  string line;

  getline( file, line );

  while( line.find( "#start", 0 ) == string::npos ) { // while line does NOT contain "#start"
    getline( file, line );
    
    if( file.eof() ) { return false; }
  }

  if( file.eof() ) { return false; }
  
  // line now contains "#start..."
  
  while( line.find( "structure", 0 ) == string::npos ) {
    getline( file, line );
  }

  getline( file, line );
  // line now contains the first line of the template sequence
  
  while( true ) { // copy the template sequence lines
    templ.append( line );
    
    if( line.length() == 0 || line[ line.length() - 1 ] == '*' ) { // check for the last line of the sequence
      break;
    }
      
    getline( file, line );
  }

  // line now contains query sequence header line
  
  while( line.find( "sequence", 0 ) == string::npos ) {
    getline( file, line );
  }
    
  getline( file, line );
  // line now contains the first line of the query sequence
  
  while( true ) { // copy the query sequence lines
    query.append( line );
    
    if( line.length() == 0 || line[ line.length() - 1 ] == '*' ) { // check for the last line of the sequence
      break;
    }
    
    getline( file, line );
  }

  format_string_ends( templ );
  format_string_ends( query );

  return true;
}


void Ali_Dist::format_string_ends( string& s )
{
  // remove * from end of sequence, if it's there
  if( s[ s.length() - 1 ] == '*' ) {
    s.resize( s.length() - 1 );
  }

  // insert initial and final fasta characters (^,$), if not already present
  if( *s.begin() != '^' ) { s.insert( s.begin(), '^' ); }
  if( s[ s.length() - 1 ] != '$' ) {
    s.insert( s.end(), '$' );
  }
}


// area calculation
int Ali_Dist::get_relative_position( float t, float q, const vector<Res_Pair>& ali_pts )
{

  // Purpose: Determine whether the point (t,q) is above(+1), below(-1), or on(0) alignment ali_pts

  // find the segment of ali_pts which bounds (t,q)
  vector<Res_Pair>::const_iterator next = ali_pts.begin();
  vector<Res_Pair>::const_iterator prev = next;
  next++;
  
  while( next != ali_pts.end() && ( (*next).t < t ) ) {
    next++;
    prev++;
  }

  if( prev == ali_pts.end() ) { // error catching
    cerr << "get_rel_pos: Could not find point within range of alignment. Exiting." << endl;
    cerr << "get_rel_pos: First template idx: " << ali_pts.front().t << endl;
    cerr << "get_rel_pos: Last template idx:  " << ali_pts.back().t << endl;
    cerr << "get_rel_pos: Current point template idx: " << t << endl;
    exit(-1);
  }

  // next and prev now bound (t,q).  next might be even with t.
  if( t == (*next).t ) {  // next and t are even

    if( q == (*next).q ) { // next and (t,q) are the same point
      return 0;
    }
    else if( q > (*next).q ) { // next is directly above (t,q)
      return 1;
    }
    else { // next is directly below (t,q)
      return -1;
    }
  }
  else { // t is between prev and next

    // determine the slope (m) and intercept (b) of the line connecting prev and next
    // NOTE: don't need to check for divide by zero since x values are always different in next and prev
    float m = ( (*next).q - (*prev).q ) / ( (*next).t - (*prev).t );
    float b = (*prev).q - ( m * (*prev).t );

    // determine the 'shadow' of (t,q) on this line above or below (t,q)
    float q_shadow = ( m * t ) + b;

    if( q == q_shadow ) { // (t,q) is on the line joining next and prev
      return 0;
    }
    else if( q > q_shadow ) { // (t,q) is above the line joining next and prev
      return 1;
    }
    else { // (t,q) is below the line joining next and prev
      return -1;
    }

  }

}


void Ali_Dist::insert_intersections( vector<Res_Pair>& a1, vector<Res_Pair>& a2 )
{

  vector<Res_Pair>::iterator prev1 = a1.begin();
  vector<Res_Pair>::iterator prev2 = a2.begin();
  vector<Res_Pair>::iterator next1 = prev1; next1++;
  vector<Res_Pair>::iterator next2 = prev2; next2++;

  while( ( next1 != a1.end() ) && ( next2 != a2.end() ) ) {

    // check for an intersection
    if( ( (*prev1).rel_pos * (*next1).rel_pos == -1.f ) ||   // points in a1 flip sign
	( (*prev2).rel_pos * (*next2).rel_pos == -1.f ) ) {  // points in a2 flip sign

      // setup variables
      float x1i = (*prev1).t;  float y1i = (*prev1).q;
      float x2i = (*prev2).t;  float y2i = (*prev2).q;
      float x1f = (*next1).t;  float y1f = (*next1).q;
      float x2f = (*next2).t;  float y2f = (*next2).q;

      float m1 = ( y1f - y1i ) / ( x1f - x1i );
      float m2 = ( y2f - y2i ) / ( x2f - x2i );

      // check if there is an intersection point
      if( m1 == m2 ) { // slopes equal - no intersection

	// move up the 'next' pointer which is further behind
	if( (*next1).t < (*next2).t ) {
	  prev1 = next1; next1++;
	}
	else if( (*next1).t > (*next2).t ) {
	  prev2 = next2; next2++;
	}
	else { // they're even, move up both
	  prev1 = next1; next1++;
	  prev2 = next2; next2++;
	}

	continue;
      }

      // find the intersection point (xp,yp)
      float xp = ( ( y1i - y2i ) - ( m1*x1i - m2*x2i ) ) / ( m2 - m1 );
      float yp = y1i + m1*( xp - x1i );

      // check if intersection point is inside segment bounds
      if( !( ( xp > x1i ) && ( xp < x1f ) && ( xp > x2i ) && ( xp < x2f ) ) ) { // out of bounds

	// move up the 'next' pointer which is further behind
	if( (*next1).t < (*next2).t ) {
	  prev1 = next1; next1++;
	}
	else if( (*next1).t > (*next2).t ) {
	  prev2 = next2; next2++;
	}
	else { // they're even, move up both
	  prev1 = next1; next1++;
	  prev2 = next2; next2++;
	}

	continue;

      }

      // insert the intersection into both a1 and a2
      Res_Pair tmp;
      tmp.t = xp;
      tmp.q = yp;
      tmp.rel_pos = 0;

      next1 = a1.insert( next1, tmp ); // next pointers set to new intersection point
      next2 = a2.insert( next2, tmp );
    }
    else { // no possible intersection here
      // move up the 'next' pointer which is further behind
      if( (*next1).t < (*next2).t ) {
	prev1 = next1; next1++;
      }
      else if( (*next1).t > (*next2).t ) {
	prev2 = next2; next2++;
      }
      else { // they're even, move up both
	prev1 = next1; next1++;
	prev2 = next2; next2++;
      }

    }

  }



}


void Ali_Dist::insert_matching_points( vector<Res_Pair>& a1, vector<Res_Pair>& a2 )
{
  
  vector<Res_Pair>::iterator prev1 = a1.begin();
  vector<Res_Pair>::iterator prev2 = a2.begin();
  vector<Res_Pair>::iterator next1 = prev1; next1++;
  vector<Res_Pair>::iterator next2 = prev2; next2++;
  
  while( ( next1 != a1.end() ) && ( next2 != a2.end() ) ) {

    // find a point that is in a1 but not a2 or vice versa
    if( (*next1).t != (*next2).t ) {

      float m, b, q_shadow;

      if( (*next1).t < (*next2).t ) { // add point to a2

	// determine the slope (m) and intercept (b) of the line connecting prev2 and next2
	// NOTE: don't need to check for divide by zero since x values are always different in next and prev
	m = ( (*next2).q - (*prev2).q ) / ( (*next2).t - (*prev2).t );
	b = (*prev2).q - ( m * (*prev2).t );

	// determine the 'shadow' of (t,q) on this line above or below (t,q)
	q_shadow = ( m * (*next1).t ) + b;

	// insert new point into a2
	Res_Pair tmp;
	tmp.t = (*next1).t;
	tmp.q = q_shadow;
	tmp.rel_pos = -1 * (*next1).rel_pos;

	next2 = a2.insert( next2, tmp ); // prev pointers set to new intersection point
      }
      else { // add point to a1

	m = ( (*next1).q - (*prev1).q ) / ( (*next1).t - (*prev1).t );
	b = (*prev1).q - ( m * (*prev1).t );
	q_shadow = ( m * (*next2).t ) + b;

	// insert new point into a1
	Res_Pair tmp;
	tmp.t = (*next2).t;
	tmp.q = q_shadow;
	tmp.rel_pos = -1 * (*next2).rel_pos;

	next1 = a1.insert( next1, tmp ); // prev pointers set to new intersection point
      }

    }
    else {
      // move up pointers
      prev1 = next1; next1++;
      prev2 = next2; next2++;
    }

  }

}


float Ali_Dist::calculate_area_between_VRPs( vector<Res_Pair>& a1, vector<Res_Pair>& a2 )
{
  // get the area between a1 and a2

  float total_area(0);

  if( a1.size() != a2.size() ) { // error checking
    cerr << "Alignments must be the same size before calculating area. ERROR!" << endl;
    cerr << "a1(main).size(): " << a1.size() << endl;
    cerr << "a1(test).size(): " << a2.size() << endl;
    exit(-1);
  }

  for( int i=1; i<(int)a2.size(); i++ ) {
      
    if( ( a1[i-1].rel_pos == 0 ) && ( a1[i].rel_pos == 0.f ) ) { // no area difference
      continue;
    }
    else {
      // get areas
      float area1 = ( ( a1[i].q + a1[i-1].q ) / 2.0f ) *
	( a1[i].t - a1[i-1].t );
      
      float area2 = ( ( a2[i].q + a2[i-1].q ) / 2.0f ) * 
	( a2[i].t - a2[i-1].t );
      
      if( ( a1[i-1].rel_pos > 0 ) || ( a1[i].rel_pos > 0 ) ) { // a1 is on top
	total_area += ( area1 - area2 );
      }
      else { // a2 must be on top
	total_area += ( area2 - area1 );
      }
      
    }
    
  }
  
  return total_area;
}


// mutual coverage calculation
float Ali_Dist::calc_templ_mutual_coverage()
{
  int common_templ_residues(0);

  int num_pairs_main = main_ali.size() - 2;
  int num_pairs_test = test_ali.size() - 2;
  float avg_num_pairs = (float)( num_pairs_main + num_pairs_test ) / 2.f;

  // initialize iterators to first pair of each alignment (not including ^,^)
  vector<Res_Pair>::iterator main_it = main_ali.begin(); main_it++;
  vector<Res_Pair>::iterator test_it = test_ali.begin(); test_it++;

  while( ( main_it != main_ali.end() ) && ( test_it != test_ali.end() ) ) { // don't count $,$

    if( main_it->t == test_it->t ) { // they are even, move both up one
      common_templ_residues++;
      main_it++;
      test_it++;
    }
    else if( main_it->t < test_it->t ) { // main_it is trailing, move main_it up one
      main_it++;
    }
    else { // main_it->t > test_it->t ... test_it is trailing, move test_it up one
      test_it++;
    }

  }

  return( (float)common_templ_residues / avg_num_pairs );
}


float Ali_Dist::calc_query_mutual_coverage()
{
  int common_query_residues(0);

  int num_pairs_main = main_ali.size() - 2;
  int num_pairs_test = test_ali.size() - 2;
  float avg_num_pairs = (float)( num_pairs_main + num_pairs_test ) / 2.f;

  // initialize iterators to first pair of each alignment (not including ^,^)
  vector<Res_Pair>::iterator main_it = main_ali.begin(); main_it++;
  vector<Res_Pair>::iterator test_it = test_ali.begin(); test_it++;

  while( ( main_it != main_ali.end() ) && ( test_it != test_ali.end() ) ) { // don't count $,$

    if( main_it->q == test_it->q ) { // they are even, move both up one
      common_query_residues++;
      main_it++;
      test_it++;
    }
    else if( main_it->q < test_it->q ) { // main_it is trailing, move main_it up one
      main_it++;
    }
    else { // main_it->t > test_it->t ... test_it is trailing, move test_it up one
      test_it++;
    }

  }

  return( (float)common_query_residues / avg_num_pairs );
}


/// PUBLIC FUNCTIONS ///
// constructor/destructor
Ali_Dist::Ali_Dist()
{

}


Ali_Dist::~Ali_Dist()
{


}


// GENERAL
// load alignments from different source formats
void Ali_Dist::load_main( string fn )
{
 // load a fasta file (fn), set template length, set main_ali

  ifstream main_ali_stream( fn.c_str() );
  string templ_main_ali, query_main_ali;
  string line;

  // extract the template and query strings from the native alignment
  getline( main_ali_stream, line );
  
  while( line[0] != '>' ) { // forward to first line of template sequence
    getline( main_ali_stream, line );
  }
  
  getline( main_ali_stream, line );  // 'line' now contains first line of template sequence
  
  while( line[0] != '>' ) { // forward to first line of query sequence
    templ_main_ali.append( line );
    getline( main_ali_stream, line );
  }
  
  getline( main_ali_stream, line );  // 'line' now contains first line of query sequence
  
  while( !main_ali_stream.eof() ) {
    query_main_ali.append( line );
    getline( main_ali_stream, line );
  }

  format_string_ends( templ_main_ali );
  format_string_ends( query_main_ali );

  templ_length = get_sequence_length( templ_main_ali );
  query_length = get_sequence_length( query_main_ali );

  // convert templ_main_ali and query_main_ali into vector<Res_Pair> main_ali
  convert_strings_to_VRP( templ_main_ali, query_main_ali, main_ali );

  //  cerr << "load_main: end" << endl;
  //  cerr << "main_ali" << endl;
  //  print_vrp( main_ali );

}


void Ali_Dist::load_main( const vector<Res_Pair>& vrp )
{
  main_ali = vrp;
  templ_length = get_sequence_length( vrp, 'T' );
  query_length = get_sequence_length( vrp, 'Q' );
}


void Ali_Dist::load_test( string fn )
{

}

void Ali_Dist::load_test( const vector<Res_Pair>& vrp )
{
  test_ali = vrp;
}

void Ali_Dist::load_PIR_file( string fn )
{
  // load a file with one or many alignments in PIR format 

}

void Ali_Dist::batch_compare_to_main_ali( string fn )
{
  // open fn, extract each PIR file in it one at a time,
  // calculate the distance between each PIR ali and main_ali,
  // store distances in batch_dists
  
  int counter(1);

  ifstream pir_file( fn.c_str() );
  
  string templ_test_ali, query_test_ali;

  while( extract_next_ali( pir_file, templ_test_ali, query_test_ali ) ) {

    // set test_ali to new alignment
    convert_strings_to_VRP( templ_test_ali, query_test_ali, test_ali );

    // get the mutual coverages
    float tmc = calc_templ_mutual_coverage();
    float qmc = calc_query_mutual_coverage();
    float product_mc = tmc * qmc;

    vector<float> tmp;
    // get the area between the two alignments and store the distance
    tmp.push_back( get_area_between_main_and_test() / (float)templ_length );

    // store the mutual coverage values
    tmp.push_back( tmc );
    tmp.push_back( qmc );
    tmp.push_back( product_mc );

    batch_dists.push_back( tmp ); // store in batch_dists

    counter++;
  }

}


float Ali_Dist::get_area_between_main_and_test()
{
  // calculate the area between main_ali and test_ali

  vector<Res_Pair> main_tmp = main_ali; // create a tmp copy of main_ali to work with

  // for each point in main_tmp, determine whether it is above (+1), below(-1), or on(0) tmp_ali
  for( int i=0; i<(int)main_tmp.size(); i++ ) {
    main_tmp[i].rel_pos = get_relative_position( main_tmp[i].t, main_tmp[i].q, test_ali );
  }

  // for each point in tmp_ali_pts, determine whether it is above (+1), below(-1), or on(0) nat_ali
  for( int i=0; i<(int)test_ali.size(); i++ ) {
    test_ali[i].rel_pos = get_relative_position( test_ali[i].t, test_ali[i].q, main_tmp );
  }
  
  // prepare alis for area calculation
  insert_intersections( main_tmp, test_ali );
  insert_matching_points( main_tmp, test_ali );

  //  cerr << "get_area: about to calculuate_area" << endl;

  return( calculate_area_between_VRPs( main_tmp, test_ali ) );
}


float Ali_Dist::get_dist_between_main_and_test()
{
  //  cerr << "get_dist_between: top" << endl;

  return( get_area_between_main_and_test() / (float)templ_length );
}


vector<Res_Pair> Ali_Dist::get_local_native_ali( int t_beg, int t_end )
{
  vector<Res_Pair> res;

  float t_res_beg = (float)t_beg;
  float t_res_end = (float)t_end;

  int idx(0);

  // advance idx so it is pointing to the first pair in main_ali within the region
  while( ( idx < (int)main_ali.size() ) && ( main_ali[idx].t < t_res_beg ) ) { idx++; }

  if( main_ali[idx].t < t_res_end ) { // there are pairs to be found in the region

    // advance idx until it leaves the region, saving the pairs it finds 
    while( ( idx < (int)main_ali.size() ) && ( main_ali[idx].t <= t_res_end ) ) { 
      res.push_back( main_ali[idx] );
      idx++; 
    }

  }

  return res;
}


float Ali_Dist::get_local_qt_shift( int t_beg, int t_end )
{

  vector<Res_Pair> local_nat_ali = get_local_native_ali( t_beg, t_end );

  if( (int)local_nat_ali.size() <= 0 ) { // error checking
    cerr << "No native pairs between template residues " << t_beg << " and  " << t_end << ". Exiting." << endl;
    exit(-1);
  }

  float sum_qt_shift(0);

  for( int i=0; i<(int)local_nat_ali.size(); i++ ) {
    sum_qt_shift += ( local_nat_ali[i].q - local_nat_ali[i].t );
  }

  sum_qt_shift /= (float)local_nat_ali.size(); // turn the sum into an average

  return sum_qt_shift;
}

void Ali_Dist::print_batch_dists( ostream& os )
{
  os << "ali#\tshift\tmin_shift" << endl;

  int min_ali_idx = -1;
  float min_dist = 100000000.f;
  
  // step through all_ali_file, parse out each pair of sequences, setup tmp_ali_pts, and find area
  
  for( int i=0; i<(int)batch_dists.size(); i++ ) {

    if( batch_dists[i][0] < min_dist ) {
      min_dist = batch_dists[i][0];
      min_ali_idx = i;
    }

    // print line
    os << i+1 << "\t" << batch_dists[i][0] << "\t" << min_dist << "\t"
       << batch_dists[i][1] << "\t" << batch_dists[i][2] << "\t" << batch_dists[i][3]
       << endl;
  }

  os << "Rank of closest:  " << min_ali_idx+1 << endl;
  os << "Shift of closest: " << min_dist << endl;

}



// OPERATORS

// ACCESS

// DEBUGGING
void Ali_Dist::print_vrp( const vector<Res_Pair>& vrp )
{
  for( int i=0; i<(int)vrp.size(); i++ ) {
    cerr << "(" << vrp[i].t << "," << vrp[i].q << ")" << endl;
  }
}

void Ali_Dist::pause()
{
  cerr << "Press a key and hit Enter to continue." << endl;
  string pause;
  cin >> pause;
}
