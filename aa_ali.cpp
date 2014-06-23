
#include <time.h>
#include <fstream>

#include "aa_seq.h"
#include "aasubalib.h"
#include "application.h"
#include "argv.h"
#include "cw.h"
#include "dpmatrix.h"
#include "fastaio.h"
#include "formats.h"
#include "noalib.h"
#include "optimal.h"
#include "pirio.h"
#include "rcfile.h"
#include "sequence.h"
#include "sflags.h"

void usage();

typedef AASubstitutionEval<AASequence,AASequence> AAEval;

int main ( int argc, const char** argv ) {

  try {

    clock_t t0 = clock();

    // Check command line arguments
    if (argc==0) usage();
    Argv args (argc,argv);
    if (args.help()) usage();
    
    // Check command line dashed arguments 
    // (& Remove them from arg list)
    string topfile;
    if (args.getSwitch("-top",false))
      args.getSwitch ("-top", 1) >> topfile;

    bool optflag = args.getSwitch("-opt",true);
    
    // Initialize program parameters
    AliParams ali_params;
    ApplicationParams app_params;
    NOaliParams noa_params;
    RCfile default_rc;
    default_rc >> ali_params >> app_params >> noa_params; 

    // Initialize topfile parameters
    if (!topfile.empty()) {
      RCfile top_rc (topfile);
      top_rc >> ali_params >> app_params;
    }
    
    // Initialize command line overrides
    args >> ali_params >> app_params >> noa_params;
    
    // Finally, check for correct number of regular arguments
    if (args.count()!=1) usage();
    
    // Setup alignment; read in sequences
    AASequence query, templ;
    ifstream seqs (args.getArg(0).str().c_str());
    cerr << "Reading in query profile" << endl;
    seqs >> Formats::FastaIn() >> templ;

    cerr << "Reading in template profile" << endl;
    seqs >> Formats::FastaIn() >> query;

    BlosumMatrix blosum (ali_params.submatrix_fn.c_str());
    AAEval ge (ali_params, blosum);
    
    DPMatrix<AASequence,AASequence,AAEval> dpm (query, templ, 
						ge,fwd, 
						ali_params.align_type);

    cout << dpm << endl;
    
    clock_t t1 = clock();
    
    Optimal<AASequence,AASequence,AAEval> opt(ali_params.align_type);
    AlignmentSet<AASequence,AASequence,AAEval> alignments (dpm, opt);
    
    if (!optflag) {
      SuboptFlags subopt(templ.size(),true);
      ConstrainedNearOptimal<AASequence,AASequence,AAEval> 
	cno (noa_params,subopt);
      cno.enumerate (dpm,alignments);
    }

    alignments.assignIdentity ();

    clock_t t2 = clock();

    switch (app_params.output_format) {
    case oFASTA:
      cout << Formats::FastaOut(app_params.line_length) << alignments;
      break;
    case oPIR:
      cout << Formats::PIROut(app_params.line_length) << alignments;
      break;
    case oHMAP:
      cerr << "Cannot use this format!\n";
      exit (-1);
    }

    double algorithm_time = (t2-t1)/(double)CLOCKS_PER_SEC;
    double cpu_time = (t2-t0)/(double)CLOCKS_PER_SEC;

    cout << "time for alignment was (sec) " << algorithm_time << endl;
    cout << "total cpu time was (sec) " << cpu_time << endl << endl;

  }
  catch (string e) { cerr << e << endl; exit(-1); }
  
}

void usage() {

  cerr << endl << "Usage: ";
  cerr << "aaa fasta_seqs" << endl;
  cerr << endl;
  cerr << "   Create optimal and near-optimal alignments using a substitution scoring matrix" << endl << endl;
  cerr << "   template.flag  specify regions for suboptimal alignment" << endl;
  cerr << "   -opt           just do an optimal alignment (template.flag is ignored)" << endl;
  cerr << "   -top <file>    specify a parameter file" << endl;
  //cerr << "   -help          get list of available parameters & default parameters" << endl; 
  cerr << "      It is possible to change these values by using the general format:" << endl;
  cerr << "      --PARAMETER_NAME value" << endl << endl;
  cerr << "      To specify the substitution matrix (Blosum-format), use:" << endl;
  cerr << "      --SUB_MATRIX filename" << endl;
  cerr << endl;
  exit(0);
  
}


