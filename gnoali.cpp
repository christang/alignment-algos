
#include <time.h>

#include "application.h"
#include "cw.h"
#include "dpmatrix.h"
#include "fastaio.h"
#include "formats.h"
#include "hmapio.h"
#include "gnoalib.h"
#include "optimal.h"
#include "pirio.h"
#include "sflags.h"

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
    
    // Initialize program parameters
    GnoaliParams ali_params;
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
    parameters.Load(Troll::Application::topology_file);
    SMAPSequence templ (args.getArg(1).str().c_str(),
			app_params.verbosity);
    
    GnoaliEval ge (ali_params);

    LogisticNormal ln (query.evd1_field,query.evd2_field,
		       templ.evd1_field,templ.evd2_field);
    
    DPMatrix<HMAPSequence,SMAPSequence,GnoaliEval> dpm (query, templ, ge, 
							fwd);

    clock_t t1 = clock();

    Optimal<HMAPSequence,SMAPSequence,GnoaliEval> opt;
    AlignmentSet<HMAPSequence,SMAPSequence,GnoaliEval> alignments (dpm, opt);
    
    if (!optflag) {
      SuboptFlags subopt(true,templ.size());
      if (args.count()>2) {
        ifstream fin (args.getArg(2).str().c_str());
        fin >> Formats::FastaIn("Flags=suboptimal region",false) >> subopt;
      }
      ConstrainedNearOptimal<HMAPSequence,SMAPSequence,GnoaliEval> 
	cno (ali_params,subopt);
      cno.enumerate (dpm,alignments);
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
    cerr << "GNOALI GNOAL. GNOA.. GNO... GN.... G....." << endl;

  }
  catch (string e) { cerr << e << endl; exit(-1); }
  
}

void usage() {

  cerr << endl << "Usage: ";
  cerr << "gnoali query.prof template.prof [template.flag]" << endl;
  cerr << endl;
  cerr << "   Create optimal and near-optimal alignments using the GNOALI method" << endl << endl;
  cerr << "   template.flag  specify regions for suboptimal alignment" << endl;
  cerr << "   -opt           just do an optimal alignment (template.flag is ignored)" << endl;  
  cerr << "   -top <file>    specify a parameter file" << endl;
  //cerr << "   -help          get list of available parameters & default parameters" << endl; 
  cerr << "      It is possible to change these values by using the general format:" << endl;
  cerr << "      --PARAMETER_NAME value" << endl;
  cerr << endl;
  exit(0);
  
}


