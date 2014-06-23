
#include <fstream>
#include <time.h>

#include "application.h"
#include "cw.h"
#include "dpmatrix.h"
#include "fastaio.h"
#include "formats.h"
#include "hmap_eval.h"
#include "hmapio.h"
#include "optimal.h"
#include "pirio.h"
#include "sflags.h"
#include "ucw.h"

void usage();

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
    bool ucwflag = args.getSwitch("-ucw",true);
    
    // Initialize program parameters
    HMAPaliParams ali_params;
    ApplicationParams app_params;
    RCfile default_rc;
    default_rc >> ali_params >> app_params; 

    // Initialize topfile parameters
    if (!topfile.empty()) {
      RCfile top_rc (topfile);
      top_rc >> ali_params >> app_params;
    }
    
    // Initialize command line overrides
    args >> ali_params >> app_params;
    
    // Finally, check for correct number of regular arguments
    if (args.count()!=2 && args.count()!=3) usage();
    
    // Setup alignment; read in sequences
    cerr << "Reading in query profile" << endl;
    HMAPSequence query (args.getArg(0).str().c_str());

    cerr << "Reading in template profile" << endl;
    Troll::Application trollApp;
    HMAPSequence templ (args.getArg(1).str().c_str());
    
    HMAPaliEval ge (ali_params);

    LogisticNormal ln (query.evd1_field,query.evd2_field,
		       templ.evd1_field,templ.evd2_field);

    DPMatrix<HMAPSequence,HMAPSequence,HMAPaliEval> dpm (query, templ, 
							 ge, fwd,
							 ali_params.align_type);

    clock_t t1 = clock();

    Optimal<HMAPSequence,HMAPSequence,HMAPaliEval> opt(ali_params.align_type);
    AlignmentSet<HMAPSequence,HMAPSequence,HMAPaliEval> alignments (dpm, opt);
    
    cerr << "Added optimal alignment to alignment set." << endl;   

    if (!optflag) {
      if (!ucwflag) {
        cerr << "Now adding constrained suboptimal alignments." << endl;
	SuboptFlags subopt(true,templ.size());
	templ.getDefaultFlags (subopt);
	if (args.count()>2) {
	  ifstream fin (args.getArg(2).str().c_str());
	  fin >> Formats::FastaIn("Flags=suboptimal region",false) >> subopt;
	}
	ConstrainedNearOptimal<HMAPSequence,HMAPSequence,HMAPaliEval> 
	  cno (ali_params,subopt);
	cno.enumerate (dpm,alignments);
      } else {
        cerr << "Now adding unconstrained suboptimal alignments." << endl;
        UnconstrainedNearOptimal<HMAPSequence,HMAPSequence,HMAPaliEval>
          ucw (ali_params);
        ucw.enumerate (dpm,alignments);
      }
    }

    alignments.assignIdentity ();
    alignments.assignSignificance (ln);

    clock_t t2 = clock();
    
    switch (app_params.output_format) {
    case oFASTA:
      cout << Formats::FastaOut(app_params.line_length) << alignments;
      break;
    case oPIR:
      cout << Formats::PIROut(app_params.line_length) << alignments;
      break;
    case oHMAP:
      cout << Formats::HMAPOut(ali_params.submatrix_fn.c_str(),
			       app_params.line_length) << alignments;
      break;
    }

    double algorithm_time = (t2-t1)/(double)CLOCKS_PER_SEC;
    double cpu_time = (t2-t0)/(double)CLOCKS_PER_SEC;

    cerr << endl;
    cerr << "time for alignment was (sec) " << algorithm_time << endl;
    cerr << "total cpu time was (sec) " << cpu_time << endl << endl;

  }
  catch (string e) { cerr << e << endl; exit(-1); }
  
}

void usage() {

  cerr << endl << "Usage: ";
  cerr << "nalign query.prof template.prof [template.flag]" << endl;
  cerr << endl;
  cerr << "   Create optimal and near-optimal alignments using the HMAP method" << endl << endl;
  cerr << "   template.flag  specify regions for suboptimal alignment" << endl;
  cerr << "   -opt           just do an optimal alignment (template.flag is ignored)" << endl;  
  cerr << "   -ucw           do standard waterman suboptimal alignment (template.flag is ignored)" << endl;
  cerr << "   -top <file>    specify a parameter file" << endl;
  //cerr << "   -help          get list of available parameters & default parameters" << endl; 
  cerr << "      It is possible to change these values by using the general format:" << endl;
  cerr << "      --PARAMETER_NAME value" << endl;
  cerr << endl;
  exit(0);
  
}


