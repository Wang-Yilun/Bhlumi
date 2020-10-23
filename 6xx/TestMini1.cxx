using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TH1.h"

// OUR headers
#include "TMCgenBHL6.h"

int main()
{
ofstream   OutFile("pro.output",ios::out);  // Logfile output

/////////////////////////////////////////////////////////////
// Random number object
TRandom3 *RN_gen = new TRandom3();       // Central r.n.gen.
long    iniseed = 54217137;
RN_gen->SetSeed(iniseed);

/////////////////////////////////////////////////////////////
// MC generator object
TMCgenBHL6 *BHLgen = new TMCgenBHL6("MCgen");
//BHLgen->ls();

/////////////////////////////////////////////////////////////
TH1D *h_NORMA = new TH1D("BHLgen_NORMA","Normalization histo",10000,0,10000);
BHLgen->Initialize(  RN_gen, &OutFile, h_NORMA);

/////////////////////////////////////////////////////////////
// Small loop over MC events
for(int iev=1; iev<=1000; iev++) {
   BHLgen->Generate();
   if(iev <= 10) BHLgen->Print1();
   cout<<" iev ="<< iev<<endl;
}

/////////////////////////////////////////////////////////////
// final printout from BHLUMI4 goes to pro.output
BHLgen->Finalize();

cout << "  |--------------------| "<<endl<<flush;
cout << "  |  TestMini1 Ended   | "<<endl<<flush;
cout << "  |--------------------| "<<endl<<flush;
return 0;
}
