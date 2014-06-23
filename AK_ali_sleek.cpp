// this code based on gnoali.cpp
// written by Andy Kuziemko (Dec 5, 2006)

#include "aa_seq.h"
#include "application.h"
#include "ssss_sleek.h"
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
    if( argc < 10 ) usage();

    Argv args (argc,argv);
    if (args.help()) usage();
    
    // Check command line dashed arguments 
    // (& Remove them from arg list)
    string topfile;
    if (args.getSwitch("-top",false))
      args.getSwitch ("-top", 1) >> topfile;

    //    bool optflag = args.getSwitch("-opt",true);
    
    // Initialize program parameters
    //    GnoaliParams ali_params;
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
    
    // Finally, check for correct number of regular arguments
    //    if (args.count()!=2 && args.count()!=3) usage();
    
    // Setup alignment; read in sequences
    cerr << "Reading in query profile" << endl;
    HMAPSequence query (args.getArg(0).str().c_str());

    cerr << "Reading in template profile" << endl;
    Troll::Application trollApp;
    parameters.Load(Troll::Application::topology_file);
    SMAPSequence templ (args.getArg(1).str().c_str(),
			app_params.verbosity);

    cerr << "loaded template" << endl;

    // get parameters
    int mode = atoi( args.getArg(2).str().c_str() );
    int num_alis_kept = atoi( args.getArg(3).str().c_str() );
    float min_cov = atof( args.getArg(4).str().c_str() );
    float min_CO = atof( args.getArg(5).str().c_str() );
    int max_SSE_skip = atoi( args.getArg(6).str().c_str() );
    int num_clusters = atoi( args.getArg(7).str().c_str() );
    int ali_how = atoi( args.getArg(8).str().c_str() );
    float max_avg_shift = atof( args.getArg(9).str().c_str() );
    ifstream ali_file1( args.getArg(10).str().c_str() );

    // load native sequence ali
    AASequence* templ1 = new AASequence;
    AASequence* query1 = new AASequence;

    FastaRead ali1( ali_file1 );

    ali1.readInto( *templ1 );
    ali1.readInto( *query1 );

    templ1->remove_double_start_ends();
    query1->remove_double_start_ends();

    //    AlignedPairList<AASequence,AASequence> native_ali( *templ1, *query1 );

    AKaliEval akev (ali_params);

    LogisticNormal ln (query.evd1_field,query.evd2_field,
		       templ.evd1_field,templ.evd2_field);

    DPMatrix<HMAPSequence,SMAPSequence,AKaliEval> dpm_fwd( query, templ, akev, 
							   fwd );
    Optimal<HMAPSequence,SMAPSequence,AKaliEval> opt;

    AlignmentSet<HMAPSequence,SMAPSequence,AKaliEval> alignments( dpm_fwd, opt );
    alignments.clear();

    cerr << "about to start SSSS" << endl;

    SSSS<HMAPSequence,SMAPSequence,AKaliEval> s_four ( ali_params, akev,
						       mode,
						       num_alis_kept,
						       min_cov,
						       min_CO,
						       max_SSE_skip,			       
						       num_clusters,
						       ali_how,
						       max_avg_shift );

    s_four.enumerate( dpm_fwd, alignments );

    /*
    alignments.assignIdentity ();
    alignments.assignSignificance (ln);

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
    */

  }
  catch (string e) { cerr << e << endl; exit(-1); }
  
}

void usage() {

  cerr << endl << "Usage: ";
  cerr << "AK_ali query.prof template.prof num_alis_total max_search min_cov min_CO "
       << "max_SSE_skip num_clusters ali_mode max_cluster_size str_ali.fa" << endl;
  cerr << endl;
  cerr << "num_alis_total    -> keep this many of the highest-scoring alignments" << endl;
  cerr << "max_search        -> maximum number of alignments to search through" << endl;
  cerr << "min cov           -> minimum coverage for an alignment to be kept" << endl;
  cerr << "min_CO            -> minimum SSE contact order needed to be kept" << endl;
  cerr << "max_SSE_skip      -> maximum deletion size in terms of template SSEs" << endl;
  cerr << "num_clusters      -> try to reduce alignments into this many clusters or fewer" << endl;
  cerr << "ali_mode      -> 0 = align fragments to full template SSEs" << endl;
  cerr << "                 1 = align to part of template SSES" << endl;
  cerr << "max_cluster_size  -> maximum distance between members of a cluster" << endl;

  cerr << "   -top <file>    specify a parameter file" << endl;
  cerr << "   -help          get list of available parameters & default parameters" << endl; 
  cerr << "      It is possible to change these values by using the general format:" << endl;
  cerr << "      --PARAMETER_NAME value" << endl;
  cerr << endl;
  exit(0);
  
}


