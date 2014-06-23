// this code based on gnoali.cpp
// written by Andy Kuziemko (Dec 5, 2006)

#include "aa_seq.h"
#include "application.h"
#include "ssss.h"
#include "dpmatrix.h"
#include "fastaio.h"
#include "formats.h"
#include "hmapio.h"
#include "hmap2_eval.h"
#include "optimal.h"
#include "optimal_rev.h"
#include "pirio.h"
#include "sflags.h"

void usage();

int main ( int argc, const char** argv ) {

  try {

    // Check command line arguments
    if( argc < 3 ) usage();

    Argv args (argc,argv);
    if (args.help()) usage();
    
    // Check command line dashed arguments 
    // (& Remove them from arg list)
    string topfile;
    if (args.getSwitch("-top",false))
      args.getSwitch ("-top", 1) >> topfile;
    
    // Initialize program parameters
    Gn2Params ali_params;
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
    HMAPSequence query (args.getArg(1).str().c_str());

    cerr << "Reading in template profile" << endl;
    Troll::Application trollApp;
    cerr << "Loaded trollApp" << endl;
    parameters.Load(Troll::Application::topology_file);
    Application::check_for_bulges = true;
    cerr << "Loaded topology_file" << endl;
    SMAPSequence templ (args.getArg(0).str().c_str(),
			app_params.verbosity);
    cerr << "Loaded template profile" << endl;

    // setup paramteres with default values
    int num_alis_returned(1000), num_alis_searched(1000000);
    float min_cov(0.4f), min_CO(0.8f);
    int max_in_betw_shift(-1);
    int ali_mode(1);
    float max_cluster_size(0.0f);
    int tracking_mode(0);
    string native_ali_fn( "" );

    // check for paramters, otherwise set to default values
    if( args.find( "max_returned" ) ) {
      num_alis_returned = atoi( args.getValue( "max_returned" ).str().c_str() );
    }

    if( args.find( "max_searched" ) ) {
      num_alis_searched = atoi( args.getValue( "max_searched" ).str().c_str() );
    }

    if( args.find( "min_cov" ) ) {
      min_cov = atof( args.getValue( "min_cov" ).str().c_str() );
    }

    if( args.find( "min_CO" ) ) {
      min_CO = atof( args.getValue( "min_CO" ).str().c_str() );
    }

    if( args.find( "max_in_betw_shift" ) ) {
      max_in_betw_shift = atoi( args.getValue( "max_in_betw_shift" ).str().c_str() );
    }

    if( args.find( "ali_mode" ) ) {
      ali_mode = atoi( args.getValue( "ali_mode" ).str().c_str() );
    }

    if( args.find( "max_cluster_size" ) ) {
      max_cluster_size = atof( args.getValue( "max_cluster_size" ).str().c_str() );
    }

    if( args.find( "str_ali" ) ) {
      native_ali_fn = args.getValue( "str_ali" ).str();
      tracking_mode = 1;
    }

    cerr << "Loaded command-line arguments." << endl;

    Hmap2Eval akev (ali_params);

    LogisticNormal ln (query.evd1_field,query.evd2_field,
		       templ.evd1_field,templ.evd2_field);

    DPMatrix<HMAPSequence,SMAPSequence,Hmap2Eval> dpm_fwd( query, templ, akev, 
							   fwd );

    cerr << "Made DPM" << endl;

    Optimal<HMAPSequence,SMAPSequence,Hmap2Eval> opt;

    cerr << "Found optimal alignment" << endl;

    AlignmentSet<HMAPSequence,SMAPSequence,Hmap2Eval> alignments( dpm_fwd, opt );

    alignments.clear();

    cerr << "About to instantiate s_four" << endl;

    SSSS<HMAPSequence,SMAPSequence,Hmap2Eval> s_four ( ali_params, akev, &dpm_fwd,
						       num_alis_returned, num_alis_searched,
						       min_cov, min_CO,
						       max_in_betw_shift, ali_mode,
						       max_cluster_size,
						       tracking_mode, native_ali_fn );

    cerr << "Constructed s_four" << endl;

    s_four.enumerate( dpm_fwd, alignments );

    cerr << "Done enumerating suboptimal alignments" << endl;

    cerr << "back from enumerate" << endl;

  }
  catch (string e) {
    cerr << e << endl; exit(-1);
  }

  cerr << "end of main{}" << endl;  
}

void usage() {

  cerr << endl << "Usage: ";
  cerr << "AK_ali template.prof query.prof --max_returned 500 --min_CO 0.5" << endl;
  cerr << endl;
  cerr << "Default values in parentheses if applicable." << endl;
  cerr << endl;
  cerr << "Required:" << endl;
  cerr << "   templ.prof    File name of the template profile" << endl;
  cerr << "   query.prof    File name of the query profile" << endl;
  cerr << endl;
  cerr << "Optional:" << endl;
  cerr << "   --max_returned N ..... will keep the highest-scoring N alignments (1K)" << endl;
  cerr << "   --max_searched S ..... will examine approximately S alignments (1M)" << endl;
  cerr << "   --min_cov ............ minimum coverage of a returned alignment (0.4)" << endl;
  cerr << "   --min_CO ............. minimum contact order of a returned alignment (0.8)" << endl;
  cerr << "   --max_in_betw_shift .. neighborhood width when searching for in-between fragment (-1)" << endl;
  cerr << "   --ali_mode ........... 0 = align fragments to full template SSEs" << endl;
  cerr << "                          1 = align to part of template SSES (default)" << endl;
  cerr << "   --max_cluster_size ... maximum distance between members of a cluster" << endl;
  cerr << endl;
  cerr << "   -top <file>    specify a parameter file" << endl;
  cerr << "   -help          get list of available parameters & default parameters" << endl; 
  cerr << "      It is possible to change these values by using the general format:" << endl;
  cerr << "      --PARAMETER_NAME value" << endl;
  cerr << endl;
  exit(0);
}


