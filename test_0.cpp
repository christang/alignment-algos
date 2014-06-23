#include <iostream>

#include "argv.h"
#include "rcfile.h"
#include "alib.h"
#include "application.h"

using namespace std;

int main (int argc, const char** argv) {

  try {
    RCfile rc;
    Argv z(argc,argv);
    AliParams a;
    rc >> a;
    z  >> a;
    
    cout << a.gap_init_penalty << endl;
    cout << a.gap_extn_penalty << endl;

    cout << "C0 " << z.count() << endl;

    char r;
    z.getSwitch ("-a",1) >> r;
    cout << r << endl;

    cout << "C1 " << z.count() << endl;

    ApplicationParams b;
    rc >> b;
    z  >> b;

    cout << "LEN=" << b.line_length << endl;
      

  } catch (string e) {
    cerr << e << endl;
  }
  
}
