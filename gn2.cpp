// Gn2.  Created by Chris Tang/Honig Lab.  
// Deposited to cvs on 11/19/07.
// Copyright 2007.  All rights reserved.

#include <time.h>

#include "application.h"
#include "cw.h"
#include "crcw.h"
#include "dpmatrix.h"
#include "fastaio.h"
#include "formats.h"
#include "hmapio.h"
#include "gn2_eval.h"
#include "kscw.h"
#include "optimal.h"
#include "pirio.h"
#include "sflags.h"
#include "ucw.h"

void usage();
void smooth_subopt_regions(SuboptFlags&);
void make_subopt_regions(SuboptFlags&,unsigned int);

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
    bool kscwflag = args.getSwitch("-kscw",true);
    bool crcwflag = args.getSwitch("-crcw",true);
    bool showrounds = args.getSwitch("-showrounds",true);

    // Initialize program parameters
    Gn2Params ali_params;                     // New for gn2.
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
    cerr << "Reading in query profile...  ";
    HMAPSequence query (args.getArg(0).str().c_str());
    cerr << "length " << query.seq_length << endl;

    cerr << "Reading in template profile...  ";
    Troll::Application trollApp;
    parameters.Load(Troll::Application::topology_file);
    SMAPSequence templ (args.getArg(1).str().c_str(),
			app_params.verbosity, true);   
                     // The true flag makes gn2 run faster.
    cerr << "length " << templ.seq_length << endl;

    Gn2Eval ge (ali_params);                 // New for gn2

    //A new statistical significance calculator will need to be added in 
    // the future.
    //LogisticNormal ln (query.evd1_field,query.evd2_field,
    //		       templ.evd1_field,templ.evd2_field);
    
    DPMatrix<HMAPSequence,SMAPSequence,Gn2Eval> dpm (query, templ, ge, 
						     fwd);

    clock_t t1 = clock();

    Optimal<HMAPSequence,SMAPSequence,Gn2Eval> opt;
    AlignmentSet<HMAPSequence,SMAPSequence,Gn2Eval> alignments (dpm, opt);

    cerr << "Added optimal alignment to alignment set." << endl;   

    if (!optflag) {
      if (ucwflag) {
        cerr << "Now adding unconstrained suboptimal alignments." << endl;
        UnconstrainedNearOptimal<HMAPSequence,SMAPSequence,Gn2Eval>
          ucw (ali_params);
        ucw.enumerate (dpm,alignments);
      } else if (kscwflag) {
	cerr << "Now adding constrained suboptimal alignments, "
	     << "with branching limited by k-sort." << endl;
        SuboptFlags subopt(true,templ.size());
        templ.getDefaultFlags (subopt);
        if (args.count()>2) {
          ifstream fin (args.getArg(2).str().c_str());
          fin >> Formats::FastaIn("Flags=suboptimal region",false) >> subopt;
        }
	KSConstrainedNearOptimal<HMAPSequence,SMAPSequence,Gn2Eval>
	  kscno (ali_params,subopt);
	kscno.enumerate (dpm,alignments);
      } else if (crcwflag) {
	
	// -- Should be a parameter to specify --
	unsigned int regions = 10;
	// -- -------------------------------- --

        SuboptFlags subopt(true,templ.size());
        templ.getDefaultFlags (subopt);
        if (args.count()>2) {
	  cerr << "Reading suboptimal regions from file." << endl;
          ifstream fin (args.getArg(2).str().c_str());
          fin >> Formats::FastaIn("Flags=suboptimal region",false) >> subopt;
        } else {
	  if (regions==0) {
	    cerr << "Using suboptimal regions defined by template. (Removing SSEs of length < 2)." << endl;
	    smooth_subopt_regions (subopt);
	  } else {
	    cerr << "Generating " << regions << " evenly-divided suboptimal regions." << endl;
	    make_subopt_regions (subopt,regions);
	  }
	}
	cerr << Formats::FastaOut() << subopt << endl;
        cerr << "Now adding constrained suboptimal alignments, "
	     << "with reduced overlap redundancy." << endl;

	// Needs k_limit, user_limit, max_overlap, rounds; now in 
        CRConstrainedNearOptimal<HMAPSequence,SMAPSequence,Gn2Eval>
          crcno (ali_params,subopt);

	int user_n = ali_params.number_suboptimal;
	ali_params.number_suboptimal = ali_params.subopt_per_round;

	AlignmentSet<HMAPSequence,SMAPSequence,Gn2Eval> ali_rounds (dpm, opt);
	  
	for (unsigned int i=1; i<=ali_params.rounds; ++i) {

	  crcno.enumerate (dpm,ali_rounds);
	  if (ali_rounds.size() < 1) {  // Looks like it converged?
	    break; 
	  }

	  templ.updateCore (ali_rounds, 0.33f);
	  dpm.reevaluate ();
	  
	  cerr <<"ROUND "<<i<<" ("<<ali_rounds.size();
	  cerr <<" alignments, opt="<<ali_rounds.front().score<<", k_limit="
               <<ali_params.k_limit<<", sort_limit="<<ali_params.sort_limit<<")"
               << endl;

	  if (showrounds) 
	    switch (app_params.output_format) {
	    case oFASTA:
	      cout << Formats::FastaOut(app_params.line_length) << ali_rounds;
	      break;
	    case oPIR:
	      cout << Formats::PIROut(app_params.line_length) << ali_rounds;
	      break;
	    case oHMAP:
	      cout << Formats::HMAPOut(ali_params.submatrix_fn.c_str(),
				       app_params.line_length) << ali_rounds;
	      break;
	    }
	  
	  ali_rounds.clear();

	}
	
	cerr << "FINAL ROUND" << endl;

        ali_params.max_overlap = ali_params.final_overlap;
	ali_params.number_suboptimal = user_n;
	
	if (ali_params.number_suboptimal == 0) {
	  cerr << "Generating new optimal alignment after rounds" << endl;
	  alignments.clear();
	  opt.enumerate (dpm, alignments);
	} else if (ali_params.number_suboptimal == 1) {
	  cerr << "Adding optimal alignment after rounds" << endl;
	  opt.enumerate (dpm, alignments);
	} else 
	  crcno.enumerate (dpm,alignments);

      } else {
        cerr << "Now adding constrained suboptimal alignments." << endl;
        SuboptFlags subopt(true,templ.size());
        templ.getDefaultFlags (subopt);
        if (args.count()>2) {
          ifstream fin (args.getArg(2).str().c_str());
          fin >> Formats::FastaIn("Flags=suboptimal region",false) >> subopt;
        }
        ConstrainedNearOptimal<HMAPSequence,SMAPSequence,Gn2Eval> 
  	  cno (ali_params,subopt);
        cno.enumerate (dpm,alignments);
      }
    }

    alignments.assignIdentity ();
    //alignments.assignSignificance (ln);

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
  cerr << "gn2 query.prof template.prof [template.flag]" << endl;
  cerr << endl;
  cerr << "   Create optimal and near-optimal alignments using the GN2 method" << endl << endl;
  cerr << "   template.flag  specify regions for suboptimal alignment" << endl;
  cerr << "   -opt           just do an optimal alignment (-ucw,-crcw,-kscw & template.flag are ignored)" << endl; 
  cerr << "   -ucw           do standard waterman suboptimal alignment (template.flag is ignored)" << endl;
  cerr << "   -crcw or -kscw experimental redundancy reduction algorithms" << endl;
  cerr << "   -top <file>    specify a parameter file" << endl;
  //cerr << "   -help          get list of available parameters & default parameters" << endl; 
  cerr << "      It is possible to change these values by using the general format:" << endl;
  cerr << "      --PARAMETER_NAME value" << endl;
  cerr << endl;
  exit(0);
  
}

void smooth_subopt_regions (SuboptFlags& sf) 
{
  // Modifies sflags to remove runs of 1s of length 1.
  // Turns them into state before the run.
  for (unsigned int i=1; i<sf.size()-1; ++i)
    if (sf[i] && !sf[i-1] && !sf[i+1]) sf.Set(i,false);
}

void make_subopt_regions (SuboptFlags& sf, unsigned int regs) 
{
  // Evenly divides 'sf' into 'regs' number of regions
  float len = (float)sf.size()/(float)regs;

  bool flag = true;
  float place = len;
  for (unsigned int i=0; i<sf.size(); ++i) {
    sf.Set (i,flag);
    if (i>place) {
      flag = !flag;
      place += len;
    }
  }
  sf.Set (sf.size()-1,true);
}
