// this code based on S4_align.cpp
// written by Andy Kuziemko (Nov 13, 2007)

#include "aa_seq.h"
#include "application.h"
#include "ssss.h"
#include "dpmatrix.h"
#include "fastaio.h"
#include "formats.h"
#include "hmapio.h"
#include "akalib.h"
#include "optimal.h"
#include "optimal_rev.h"
#include "pirio.h"
#include "sflags.h"

void usage();

int main ( int argc, const char** argv ) {

  try {

    // Check command line arguments
    if( argc < 7 ) usage();

    Argv args (argc,argv);
    if (args.help()) usage();
    
    // Check command line dashed arguments 
    // (& Remove them from arg list)
    string topfile;
    if (args.getSwitch("-top",false))
      args.getSwitch ("-top", 1) >> topfile;
    
    // Initialize program parameters
    AKaliParams ali_params;
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
    
    // Setup alignment; read in sequences
    cerr << "Reading in query profile" << endl;
    HMAPSequence query (args.getArg(0).str().c_str());

    cerr << "Reading in template profile" << endl;
    Troll::Application trollApp;
    cerr << "Loaded trollApp" << endl;
    parameters.Load(Troll::Application::topology_file);
    Application::check_for_bulges = true;
    cerr << "Loaded topology_file" << endl;
    SMAPSequence templ (args.getArg(1).str().c_str(),
			app_params.verbosity);
    cerr << "Loaded template profile" << endl;

    int num_alis_kept = atoi( args.getArg(2).str().c_str() );
    int num_alis_searched = atoi( args.getArg(3).str().c_str() );
    float min_cov = atof( args.getArg(4).str().c_str() );
    float min_CO = atof( args.getArg(5).str().c_str() );
    int ali_how = atoi( args.getArg(6).str().c_str() );
    float max_avg_shift = atof( args.getArg(7).str().c_str() );
    //    int max_SSE_skip = atoi( args.getArg(8).str().c_str() );
    //    int num_clusters = atoi( args.getArg(9).str().c_str() );

    cerr << "Loaded command-line arguments." << endl;

    AKaliEval akev (ali_params);

    LogisticNormal ln (query.evd1_field,query.evd2_field,
		       templ.evd1_field,templ.evd2_field);

    DPMatrix<HMAPSequence,SMAPSequence,AKaliEval> dpm_fwd( query, templ, akev, 
							   fwd );

    cerr << "Made DPM" << endl;

    Optimal<HMAPSequence,SMAPSequence,AKaliEval> opt;

    cerr << "Found optimal alignment" << endl;

    AlignmentSet<HMAPSequence,SMAPSequence,AKaliEval> alignments( dpm_fwd, opt );

    alignments.clear();

    cerr << "About to instantiate s_four" << endl;

    SSSS<HMAPSequence,SMAPSequence,AKaliEval> s_four ( ali_params, akev, &dpm_fwd,
						       num_alis_kept, num_alis_searched,
						       min_cov, min_CO,
						       ali_how,
						       max_avg_shift );

    cerr << "Constructed s_four" << endl;

    s_four.choose_fragments_for_ali();


  }
  catch (string e) {
    cerr << e << endl; exit(-1);
  }

  cerr << "end of main{}" << endl;  
}

void usage() {

  cerr << endl << "Usage: ";
  cerr << "AK_ali query.prof template.prof num_alis_total max_search min_cov min_CO "
       << "ali_mode max_cluster_size" << endl;
  cerr << endl;
  cerr << "num_alis_total    -> keep this many of the highest-scoring alignments" << endl;
  cerr << "max_search        -> maximum number of alignments to search through" << endl;
  cerr << "min cov           -> minimum coverage for an alignment to be kept" << endl;
  cerr << "min_CO            -> minimum SSE contact order needed to be kept" << endl;
  cerr << "ali_mode      -> 0 = align fragments to full template SSEs" << endl;
  cerr << "                 1 = align to part of template SSES" << endl;
  cerr << "max_cluster_size  -> maximum distance between members of a cluster" << endl;
  //  cerr << "max_SSE_skip      -> maximum deletion size in terms of template SSEs" << endl;
  //  cerr << "num_clusters      -> try to reduce alignments into this many clusters or fewer" << endl;

  cerr << "   -top <file>    specify a parameter file" << endl;
  cerr << "   -help          get list of available parameters & default parameters" << endl; 
  cerr << "      It is possible to change these values by using the general format:" << endl;
  cerr << "      --PARAMETER_NAME value" << endl;
  cerr << endl;
  exit(0);
  
}


