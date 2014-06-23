// written by Andy Kuziemko (Jul 30, 2007)

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "ali_dist.h"

using namespace std;

void usage();

int main ( int argc, const char** argv ) {

  Ali_Dist X;

  X.load_main( (string)argv[2] );

  X.batch_compare_to_main_ali( (string)argv[1] );

  X.print_batch_dists();
  
}

void usage() {

  cerr << "Calculates the area of difference between two alignments " << endl;
  cerr << "with identical template and query sequences." << endl;
  exit(0);
  
}


