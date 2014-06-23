// Variation on get_area_diffs

#include <fstream>
#include <sstream>

#include "aa_seq.h"
#include "ali_dist.h"
#include "argv.h"
#include "fastaio.h"
#include "gn2lib_seq.h"
#include "pirio.h"
#include "sflags.h"

void usage();

typedef AlignedPairList<AASequence,AASequence> Alignment;

template <class S1, class S2>
int AlignedPairList<S1,S2>::get_shift ( AlignedPairList<S1,S2>& native, 
          string& qstr, SuboptFlags& core,
          int *ali_len)
{
  if (get_last_query_idx() != (int)core.size()-1) 
    throw string ("Core file length does not match alignment");
  
  AlignmentSet<S1,S2,AASubstitutionEval> as (native);
  as.push_back (*this);

  SequenceGaps sg (as,native.get_last_query_idx()+1,
       native.get_last_template_idx()+1);
  
  string tpl_gapped = "";
  string nat_gapped = "";
  string ali_gapped = "";

  // Mask non-core sequence
  // Only care about core parts of the query
  for (unsigned int i=0; i<qstr.length(); ++i)
    if (!core[i]) qstr[i]='.';

  // Mask zigzag regions
  // Pretend that any part of the query in the
  // zigzag part of a native ali is never wrongly aligned
  // in the test ali
  typename list<Pair>::const_iterator pair_it  = native.begin();
  typename list<Pair>::const_iterator pair_end = native.end();
  const Pair* prev = &*pair_it; ++pair_it;
  for (; pair_it!=pair_end; ++pair_it) {
    // check for zigzag region
    if (pair_it->query_idx()-prev->query_idx()>1 &&
  pair_it->template_idx()-prev->template_idx()>1)
      for (int i=prev->query_idx()+1; i<pair_it->query_idx();++i)
  qstr[i]='.';
    prev = &*pair_it;
  }

  string tstr = string (native.get_last_template_idx()+1,'*');

  sg.build (tstr, tpl_gapped, '-');
  sg.build (qstr, native, nat_gapped, '-');
  sg.build (qstr,  *this, ali_gapped, '-');

  //cerr << "templ : " << tpl_gapped << endl;
  //cerr << "Native: " << nat_gapped << endl;
  //cerr << "   Ali: " << ali_gapped << endl;

  if (ali_len!=0) *ali_len=-2;
  int diff(0), shift(0);

  // The Difference-based Shift score
    // This is equivalent to seq_Sum_i ( abs ( nat_R_i - ali_R_i ) )
    // where i is the index of summation, 
    // nat_R_i is the position of residue i in the native ali
    // and ali_R_i is the position of the same residue in the test ali

  for (unsigned int i=0; i<nat_gapped.length(); ++i) {
    // shift
    if (nat_gapped[i]!='-'&&nat_gapped[i]!='.') ++diff;
    if (ali_gapped[i]!='-'&&ali_gapped[i]!='.') --diff;
    shift += abs (diff);

    // alignment length
    if (ali_len!=0 && ali_gapped[i]!='-' && tpl_gapped[i]!='-') 
      (*ali_len)++;

    //cerr << tpl_gapped[i] << nat_gapped[i] << ali_gapped[i] << " : " 
    //   << diff << " : " << shift << " l=" << *ali_len << endl;
  }
  return shift;
}

int main ( int argc, const char** argv ) {

  if (argc<3) usage();

  Argv args (argc,argv);
  bool use_all = args.getSwitch("-all",true);

  // Do area-based shift calculation

  Ali_Dist X;
  X.load_main( (string)argv[2] );
  X.batch_compare_to_main_ali( (string)argv[1] );

  // Do residue-based shift calculation

  ifstream all_ali_file( argv[1] );
  ifstream nat_ali_file( argv[2] );

  Alignment seq_ali, nat_ali, opt_ali;

  nat_ali_file >> Formats::FastaAlignmentIn() >> nat_ali;
  
  int rank = 0;
  float min_area_based=999999999.f;
  int min_res_based=999999999;
  int max_agree=-1;
  float max_q_mod=-1.f;
  float max_q_dev=-1.f;
  float max_q_comb=-1.f;
  float max_from_opt=-1.f;
  float len=(float) nat_ali.get_last_template_idx()-1;

  char buffer[50];
  stringstream partII;
  if (argc>3) printf ("Using core definitions\n");
  else printf ("Using all residues\n");
  printf ("Native alignment length: %d\n",nat_ali.size());
  printf ("Native alignment %%ID: %4.2f\n",nat_ali.identity);
  printf ("\nRunning statistics\n");
  printf ("Rank \t%%ID\t#ali'd\tshift_r\tshift_a\t#agree\tQ_mod\tQ_dev\tQ_comb\n");

  partII << "\nCummulative statistics\n";
  partII << "Rank \t%ID\t#ali'd\tshift_r\tshift_a\t#agree\tQ_mod\tQ_dev\tQ_comb\n";

  int bd_idx = 0;
  
  string q_seq (nat_ali.get_last_query_idx()+1,'*');
 
  SuboptFlags *allr = new SuboptFlags (true, nat_ali.get_last_query_idx()+1);
  SuboptFlags *core = new SuboptFlags (true, nat_ali.get_last_query_idx()+1);
  if (argc>3) {
    Troll::Application app;
    parameters.Load( Troll::Application::topology_file );
    SMAPSequence smap (argv[3],0,true);
    q_seq = *smap.getString();

    if (!use_all) smap.getDefaultFlags( *core );
  }
  
  bool opt = true;

  while (1) {

    try {
      all_ali_file >> Formats::PIRIn() >> seq_ali;
      if (opt) opt_ali = seq_ali;
    } catch (...) {
      break;
    }

    int ali_len;
    //float area = seq_ali.get_area_diff (nat_ali);
    float area_based = X.batch_dists[bd_idx++][0];
    int res_based = seq_ali.get_shift (nat_ali, q_seq, *core, &ali_len);
    float pctid = seq_ali.identity;
    
    int n_agree;
    float q_mod, q_dev, q_comb;
    seq_ali.get_q_all (nat_ali, *allr, &n_agree, &q_mod, &q_dev, &q_comb);

    min_area_based = min ( min_area_based, area_based );
    min_res_based = min ( min_res_based, res_based );
    max_agree = max ( max_agree, n_agree );
    max_q_mod = max ( max_q_mod, q_mod );
    max_q_dev = max ( max_q_dev, q_dev );
    max_q_comb = max ( max_q_comb, q_comb );

    //int tf = 10* n_agree - res_based;

    cout << rank << "\t";
    partII << rank++ << "\t";

    sprintf (buffer,"%4.2f",pctid);
    cout << buffer << "\t";
    partII << buffer << "\t";

    cout << ali_len << "\t";
    partII << ali_len << "\t";

    cout << res_based << "\t";
    partII << min_res_based << "\t";

    sprintf (buffer,"%4.2f",area_based);
    cout << buffer << "\t";
    sprintf (buffer,"%4.2f",min_area_based);
    partII << buffer << "\t";

    cout << n_agree << "\t";
    partII << max_agree << "\t";

    sprintf (buffer,"%4.2f",q_mod * 100.f);
    cout << buffer << "\t";
    sprintf (buffer,"%4.2f",max_q_mod * 100.f);
    partII << buffer << "\t";
    
    sprintf (buffer,"%4.2f",q_dev * 100.f);
    cout << buffer << "\t";
    sprintf (buffer,"%4.2f",max_q_dev * 100.f);
    partII << buffer << "\t";
    
    sprintf (buffer,"%4.2f",q_comb * 100.f);
    cout << buffer << "\t";
    sprintf (buffer,"%4.2f",max_q_comb * 100.f);
    partII << buffer << "\t";

    if (!opt) {
      float from_opt = seq_ali.get_area_diff (opt_ali);
      sprintf (buffer,"%4.2f",from_opt/len);
      cout << buffer;

      max_from_opt = max ( max_from_opt, from_opt );
      sprintf (buffer,"%4.2f",max_from_opt/len);
      partII << buffer;
    }

    cout << "\t[R]" << endl;
    partII << "\t[C]" << endl;

    opt = false;

  }
  
  cout << partII.str();

    //cerr << "native:"  << endl;
    //nat_ali.print_pairs();
    
    //cerr << "alignment:" << endl;
    //seq_ali.print_pairs();
    
    //cerr << "core:" << endl; 
    //cerr << Formats::FastaOut() << *core << endl;
    
}


void usage()
{
  cerr << "Calculates the shift between two alignments " << endl;
  cerr << "with identical template and query sequences." << endl;
  cerr << endl;
  cerr << "get_shifts <seq ali> <nat ali> [core flags]" << endl;
  exit(0);
}
