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
//#include "TSystem.h"

// OUR headers
#include "TMCgenBHL6.h"

int main()
{
TFile *MCgenFile = new TFile("mcgen.root","UPDATE","Generators");
ofstream   OutFile("pro.output",ios::out);  // Logfile output

MCgenFile->cd();
/////////////////////////////////////////////////////////////
// Random number object
TRandom3 *RN_gen = new TRandom3();       // Central r.n.gen.
long    iniseed = 54217137;
RN_gen->SetSeed(iniseed);
RN_gen->Write("RN_gen",TObject::kOverwrite);

/////////////////////////////////////////////////////////////
// MC generator object
TMCgenBHL6 *BHLgen = new TMCgenBHL6("MCgen");

BHLgen->Write("MCgen",TObject::kOverwrite);

cout<<"------------------------MCgenFile.GetListOfKeys-----------------------------"<<endl;
MCgenFile->GetListOfKeys()->Print();
cout<<"-------------------------MCgenFile.ShowStreamerInfo-------------------------"<<endl;
MCgenFile->ShowStreamerInfo();
cout<<"-------------------------MCgenFile.ls---------------------------------------"<<endl;
MCgenFile->ls();

MCgenFile->Write();
MCgenFile->Close();

cout << "  |--------------------| "<<endl<<flush;
cout << "  |  TestMini2 Ended   | "<<endl<<flush;
cout << "  |--------------------| "<<endl<<flush;
return 0;
}
