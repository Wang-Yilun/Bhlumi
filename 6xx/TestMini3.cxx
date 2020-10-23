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

TFile HistoFile("histo.root","UPDATE","histograms");
TFile MCgenFile("mcgen.root","UPDATE","Generators");
int main()
{
ofstream   OutFile("pro.output",ios::out);  // Logfile output

HistoFile.cd();
/////////////////////////////////////////////////////////////
// Normalization histogram (normaly in TRobol)
TH1D *h_NORMA = new TH1D("BHLgen_NORMA","Normalization histo",10000,0,10000);

TRandom3  *RN_gen = (TRandom3*)MCgenFile.Get("RN_gen");  // read r.n. generator
//
TMCgenBHL6 *BHLgen = (TMCgenBHL6*)MCgenFile.Get("MCgen");
//BHLgen->ls();

BHLgen->Initialize(  RN_gen, &OutFile, h_NORMA);

/////////////////////////////////////////////////////////////
// Small loop over MC events
for(int iev=1; iev<=1000; iev++) {
   BHLgen->Generate();
   if(iev <= 10) BHLgen->Print1();
   if(iev<10) cout<<" iev ="<< iev<<endl;
}

/////////////////////////////////////////////////////////////
// final printout from BHLUMI4 goes to pro.output
BHLgen->Finalize();

cout << "  |--------------------| "<<endl<<flush;
cout << "  |  TestMini1 Ended   | "<<endl<<flush;
cout << "  |--------------------| "<<endl<<flush;
return 0;
}
