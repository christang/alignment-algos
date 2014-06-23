// Report residue characteristics from structure-based seq ali.
// Print contact number, accessibility, and alignment state of
// template and give template reference sequence.

#include <fstream>

#include "aa_seq.h"
#include "fastaio.h"
#include "gn2lib_seq.h"

void usage();

typedef AlignedPairList<AASequence,AASequence> Alignment;

int main ( int argc, char** argv ) {

  if (argc!=4) usage();

  try{

  Alignment ali;
  
  ifstream ali_file ( argv[1] );

  // in FastaAlignment, template is on top

  ali_file >> Formats::FastaAlignmentIn() >> ali;
  ali.remove_ends();

// query seq, resetting and rewinding stream
//   AASequence tmp, query;
//   ali_file. clear ();
//   ali_file. seekg (0, ios::beg);
//   ali_file >> Formats::FastaIn() >> tmp;
//   ali_file >> Formats::FastaIn() >> query;
//   query.cleargaps('-');

  // cerr << Formats::FastaOut() << query;

  // setup gn2 profile
  Troll::Application app;
  parameters.Load( Troll::Application::topology_file );
  SMAPSequence prof ( argv[2] );

  // query profile
  HMAPSequence hmap ( argv[3] );

  // start from first aligned position
  Alignment::const_iterator it=ali.begin();
  int idx = it->query_idx();
  int ali_idx = it->template_idx();

  for (;  it!=ali.end();  ++it) {
    while (idx < it->query_idx()) {
      cout << 2 << "\t(" << it->query_idx()-idx << ")\t" 
	   << "-\t-" << endl;
      idx = it->query_idx();
    }
    while (ali_idx < it->template_idx()) {
      cout << 0 << "\t" << prof.weighted_contact_number[ali_idx] << "\t" 
	   << prof.at(ali_idx)->rdata.accessibility << "\t"
	   << "-" << "\t"
	   << prof.at(ali_idx)->olc << endl;
      ++ali_idx;
    }
    cout << 1 << "\t" << prof.weighted_contact_number[ali_idx] << "\t" 
	 << prof.at(ali_idx)->rdata.accessibility << "\t"
	 << hmap.at(idx)->hydropathy << "\t"
	 << prof.at(ali_idx)->hydropathy << "\t";
    if (prof.at(ali_idx)->p_coil()>prof.at(ali_idx)->p_strand()&&
        prof.at(ali_idx)->p_coil()>prof.at(ali_idx)->p_helix()) 
      cout << "c" << "\t";
    else if (prof.at(ali_idx)->p_strand()>prof.at(ali_idx)->p_coil()&&
             prof.at(ali_idx)->p_strand()>prof.at(ali_idx)->p_helix()) 
      cout << "e" << "\t";
    else if (prof.at(ali_idx)->p_helix()>prof.at(ali_idx)->p_strand()&&
             prof.at(ali_idx)->p_helix()>prof.at(ali_idx)->p_coil())
      cout << "h" << "\t";
    else { cerr << "error" << endl; exit(1); }

    cout << hmap.at(idx)->olc << "\t"
	 << prof.at(ali_idx)->olc << endl;
    ++idx;
    ++ali_idx;
  }

  }
  catch (string s) {
    cerr << s << endl;
  }
  catch (...) {
    cerr << "Error" << endl;
  }
  
}

void usage() {

  cerr << "Prints three columns of numbers" << endl;
  cerr << "Del/Ali/Ins (0/1)  Weighted contact number/(Ins length)  Accessibility  query Hydropathy  templ Hydropathy  query AA  templ AA" << endl;
  cerr << "Usage: cn_acc_analysis <ali> <templ prof> <query prof>" << endl;
  cerr << endl;
  exit(-1);
  
}

