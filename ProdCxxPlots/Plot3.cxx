//    make Plot3-run
//    make Zstudy-pdf

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

// ROOT headers
#include <math.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLatex.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TMarker.h"
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"

#include "TGraphPainter.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TColor.h"

#include "TMCgenBHL6.h"
#include "TRobolProd.h"

#include "HSTplot.h"
#include "HisNorm.h"
///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
//
TFile *gHstFile88 = new TFile("../ProdCxx/88GeV/histo.root");  // File with histograms
TFile *gGenFile88 = new TFile("../ProdCxx/88GeV/mcgen.root");  // File with Generator objects
//
TFile *gHstFile89 = new TFile("../ProdCxx/89GeV/histo.root");  // File with histograms
TFile *gGenFile89 = new TFile("../ProdCxx/89GeV/mcgen.root");  // File with Generator objects
//
TFile *gHstFile90 = new TFile("../ProdCxx/90GeV/histo.root");  // File with histograms
TFile *gGenFile90 = new TFile("../ProdCxx/90GeV/mcgen.root");  // File with Generator objects
//
TFile *gHstFile91 = new TFile("../ProdCxx/91GeV/histo.root");  // File with histograms
TFile *gGenFile91 = new TFile("../ProdCxx/91GeV/mcgen.root");  // File with Generator objects
//
TFile *gHstFile92 = new TFile("../ProdCxx/92GeV/histo.root");  // File with histograms
TFile *gGenFile92 = new TFile("../ProdCxx/92GeV/mcgen.root");  // File with Generator objects
//
TFile *gHstFile93 = new TFile("../ProdCxx/93GeV/histo.root");  // File with histograms
TFile *gGenFile93 = new TFile("../ProdCxx/93GeV/mcgen.root");  // File with Generator objects
//
TFile *gHstFile94 = new TFile("../ProdCxx/94GeV/histo.root");  // File with histograms
TFile *gGenFile94 = new TFile("../ProdCxx/94GeV/mcgen.root");  // File with Generator objects

TH1D *gst_ENE;  //list of energies

double gEne[7] = { 88.471e0, 89.444e0, 90.216e0, 91.227e0, 91.959e0, 93.000e0, 93.710e0};
double gEne1[7]= { 88.371e0, 89.344e0, 90.116e0, 91.127e0, 91.859e0, 92.900e0, 93.610e0};
double gEne2[7]= { 88.571e0, 89.544e0, 90.316e0, 91.327e0, 92.059e0, 93.000e0, 93.810e0};

//
// Dump file for temporarary objects
TFile DiskFileB("Plot1.root","RECREATE","histograms");

ofstream   OutFile("Bhl6.output",ios::out);  // Logfile output
//=============================================================================
//Double_t sqr( const Double_t x ){ return x*x;};
//
float  gXcanv = 50, gYcanv = 50;
TMCgenBHL6 *gBHLgen;
double gCMSENE;

TH1D *hst_Zexp_Al_LCAL, *hst_Zne1_Al_LCAL, *hst_Zne2_Al_LCAL;
TH1D *hst_Zexp_SICAL92, *hst_Zne1_SICAL92, *hst_Zne2_SICAL92;
TH1D *hst_Zexp_Op_FCAL, *hst_Zne1_Op_FCAL, *hst_Zne2_Op_FCAL;
TH1D *hst_Zexp_L3_BGO1, *hst_Zne1_L3_BGO1, *hst_Zne2_L3_BGO1; // L3 BGO
TH1D *hst_Zexp_DEL_SAT, *hst_Zne1_DEL_SAT, *hst_Zne2_DEL_SAT; // DELHI SAT
TH1D *hst_Zexp_SICAL94, *hst_Zne1_SICAL94, *hst_Zne2_SICAL94;

HSTplot LibPLT("LibPLT");
//=============================================================================


///////////////////////////////////////////////////////////////////////////////////
void HistRemake(TFile *HstFile ){
//
cout<<"----------------------------- HistRemake ------------------------------------"<<endl;

TH1D *HST_BHL_NORMA = (TH1D*)HstFile->Get("HST_BHL_NORMA");
///////////////////////////////////////////////////////////////////////////
// renormalizing 1-dim histos
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_ZZ0") );
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_ZZ1") );
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_ZZ2") );
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_ZZ3") );
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_ZZ4") );

    cout<<"----------------------------- HistRemake end ----------------------------------"<<endl;
}//HistRemake


///////////////////////////////////////////////////////////////////////////////////
void TabTrig(TFile *HstFile, TString TeXfile)
{
//------------------------------------------------------------------------
  cout<<" ========================= TabTrig start=========================== "<<endl;
  //
  TH1D *hst_trig_ZZ0  = (TH1D*)HstFile->Get("hst_trig_ZZ0");
  TH1D *hst_trig_ZZ1  = (TH1D*)HstFile->Get("hst_trig_ZZ1");
  TH1D *hst_trig_ZZ2  = (TH1D*)HstFile->Get("hst_trig_ZZ2");
  TH1D *hst_trig_ZZ3  = (TH1D*)HstFile->Get("hst_trig_ZZ3");
  TH1D *hst_trig_ZZ4  = (TH1D*)HstFile->Get("hst_trig_ZZ4");

//  for(int i=0; i<=100; i++ ) hst_trig_ZZ1->SetBinError(i,0); // errors omitted


// Column captions
  int nPlt=5;   // No of histos/column +1
  Char_t *Capt[nPlt+1];
  for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
//
  strcpy(Capt[0],"{\\color{blue} No.}");
  strcpy(Capt[1],"{\\color{blue} Z off }");
  strcpy(Capt[2],"{\\color{blue} Z on }");
  strcpy(Capt[3],"{\\color{blue} Z exp. wt12}");
  strcpy(Capt[4],"{\\color{blue} Z alf1 wt82}");
  strcpy(Capt[5],"{\\color{blue} Z alf1 wt92}");
//  strcpy(Capt[6],"{\\color{blue} new-old }");

// pointers to histograms
  TH1D *iHst[nPlt+1];

// multicolumn caption
  Char_t Mcapt[132];
///************************************
  // Latex source file
  FILE *DiskFileTeX = fopen(TeXfile+".txp","w");
//************************************
// Initialization of the latex source file
  LibPLT.PlInit(DiskFileTeX, 2);

  LibPLT.Setfmt0("$  %10.0f $");
//  LibPLT.Setfmt2("$  %11.5f \\pm %9.5f $");  // next columns

  iHst[1]= hst_trig_ZZ0;  // VP off
  iHst[2]= hst_trig_ZZ1;  // VPold on
  iHst[3]= hst_trig_ZZ2;  // VPold_err
  iHst[4]= hst_trig_ZZ3;  // VPnew
  iHst[5]= hst_trig_ZZ4;  // VPnew err

  strcpy(Mcapt,"{\\color{red} ALEPH LCAL 1990-92 (BHLUMI 2.x) [nb] \\hfill }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "B",1,1, 1); // for 100 bins

  strcpy(Mcapt,"{\\color{red} OPAL FCAL 1990-92 (BHLUMI 4.x)  [nb] \\hfill }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",2,2, 1); // for 100 bins

  strcpy(Mcapt,"{\\color{red} L3 BGO 1990-92 (BHLUMI 2.x)     [nb] \\hfill } ");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",3,3, 1); // for 100 bins

  strcpy(Mcapt,"{\\color{red} DELPHI SAT 1990-92 (BHLUMI 2.x) [nb] \\hfill } ");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",4,4, 1); // for 100 bins

  strcpy(Mcapt,"{\\color{red} ALEPH SICAL 1992 (BHLUMI 2.x)   [nb] \\hfill } ");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",5,5, 1); // for 100 bins

  strcpy(Mcapt,"{\\color{red} ALEPH SICAL 1993, 94, 95 (BHLUMI 4.x) [nb] \\hfill }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "t",11,13, 1); // for 100 bins

  strcpy(Mcapt,"{\\color{red} OPAL OSIW 1993-95 (BHLUMI 4.x)   [nb] \\hfill }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",14,14, 1); // for 100 bins

  strcpy(Mcapt,"{\\color{red} L3 BGO 1993-95 (BHLUMI 4.x)      [nb] \\hfill } ");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",15,15, 1); // for 100 bins

  strcpy(Mcapt,"{\\color{red} DELPHI STIC 1994-95 (BHLUMI 4.x) [nb] \\hfill }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",16,16, 1); // for 100 bins

  strcpy(Mcapt,"{\\color{red} DELPHI SAT     1993 (BHLUMI 4.x) [nb] \\hfill }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "E",17,17, 1); // for 100 bins

// finalizing latex source file
  LibPLT.PlEnd(DiskFileTeX);
//************************************
  fclose(DiskFileTeX);
//************************************
  cout<<" ========================= TabTrig end =========================== "<<endl;
}// TabTrig

///////////////////////////////////////////////////////////////////////////////////
void TabZet()
{
//************************************
  cout<<" ========================= TabZet start ========================= "<<endl;
//
 TH1D *hst_ZZ1_88  = (TH1D*)gHstFile88->Get("hst_trig_ZZ1");
 TH1D *hst_ZZ1_89  = (TH1D*)gHstFile89->Get("hst_trig_ZZ1");
 TH1D *hst_ZZ1_90  = (TH1D*)gHstFile90->Get("hst_trig_ZZ1");
 TH1D *hst_ZZ1_91  = (TH1D*)gHstFile91->Get("hst_trig_ZZ1");
 TH1D *hst_ZZ1_92  = (TH1D*)gHstFile92->Get("hst_trig_ZZ1");
 TH1D *hst_ZZ1_93  = (TH1D*)gHstFile93->Get("hst_trig_ZZ1");
 TH1D *hst_ZZ1_94  = (TH1D*)gHstFile94->Get("hst_trig_ZZ1");
//
 TH1D *hst_ZZ2_88  = (TH1D*)gHstFile88->Get("hst_trig_ZZ2");
 TH1D *hst_ZZ2_89  = (TH1D*)gHstFile89->Get("hst_trig_ZZ2");
 TH1D *hst_ZZ2_90  = (TH1D*)gHstFile90->Get("hst_trig_ZZ2");
 TH1D *hst_ZZ2_91  = (TH1D*)gHstFile91->Get("hst_trig_ZZ2");
 TH1D *hst_ZZ2_92  = (TH1D*)gHstFile92->Get("hst_trig_ZZ2");
 TH1D *hst_ZZ2_93  = (TH1D*)gHstFile93->Get("hst_trig_ZZ2");
 TH1D *hst_ZZ2_94  = (TH1D*)gHstFile94->Get("hst_trig_ZZ2");
 //
 TH1D *hst_ZZ3_88  = (TH1D*)gHstFile88->Get("hst_trig_ZZ3");
 TH1D *hst_ZZ3_89  = (TH1D*)gHstFile89->Get("hst_trig_ZZ3");
 TH1D *hst_ZZ3_90  = (TH1D*)gHstFile90->Get("hst_trig_ZZ3");
 TH1D *hst_ZZ3_91  = (TH1D*)gHstFile91->Get("hst_trig_ZZ3");
 TH1D *hst_ZZ3_92  = (TH1D*)gHstFile92->Get("hst_trig_ZZ3");
 TH1D *hst_ZZ3_93  = (TH1D*)gHstFile93->Get("hst_trig_ZZ3");
 TH1D *hst_ZZ3_94  = (TH1D*)gHstFile94->Get("hst_trig_ZZ3");
 //
 TH1D *hst_ZZ4_88  = (TH1D*)gHstFile88->Get("hst_trig_ZZ4");
 TH1D *hst_ZZ4_89  = (TH1D*)gHstFile89->Get("hst_trig_ZZ4");
 TH1D *hst_ZZ4_90  = (TH1D*)gHstFile90->Get("hst_trig_ZZ4");
 TH1D *hst_ZZ4_91  = (TH1D*)gHstFile91->Get("hst_trig_ZZ4");
 TH1D *hst_ZZ4_92  = (TH1D*)gHstFile92->Get("hst_trig_ZZ4");
 TH1D *hst_ZZ4_93  = (TH1D*)gHstFile93->Get("hst_trig_ZZ4");
 TH1D *hst_ZZ4_94  = (TH1D*)gHstFile94->Get("hst_trig_ZZ4");

 hst_Zexp_Al_LCAL  = new TH1D("hst_Zexp_Al_LCAL", "Z expon. ",  7, 0.0, 7.0); hst_Zexp_Al_LCAL->Sumw2();
 hst_Zne1_Al_LCAL  = new TH1D("hst_Zne1_Al_LCAL", "Z O(alf1)",  7, 0.0, 7.0); hst_Zne1_Al_LCAL->Sumw2();
 hst_Zne2_Al_LCAL  = new TH1D("hst_Zne2_Al_LCAL", "Z O(alf1)",  7, 0.0, 7.0); hst_Zne2_Al_LCAL->Sumw2();

 hst_Zexp_Op_FCAL  = new TH1D("hst_Zexp_Op_FCAL", "Z expon. ",  7, 0.0, 7.0); hst_Zexp_Op_FCAL->Sumw2();
 hst_Zne1_Op_FCAL  = new TH1D("hst_Zne1_Op_FCAL", "Z O(alf1)",  7, 0.0, 7.0); hst_Zne1_Op_FCAL->Sumw2();
 hst_Zne2_Op_FCAL  = new TH1D("hst_Zne2_Op_FCAL", "Z O(alf1)",  7, 0.0, 7.0); hst_Zne2_Op_FCAL->Sumw2();

 hst_Zexp_L3_BGO1  = new TH1D("hst_Zexp_L3_BGO1", "Z expon. ",  7, 0.0, 7.0); hst_Zexp_L3_BGO1->Sumw2();
 hst_Zne1_L3_BGO1  = new TH1D("hst_Zne1_L3_BGO1", "Z O(alf1)",  7, 0.0, 7.0); hst_Zne1_L3_BGO1->Sumw2();
 hst_Zne2_L3_BGO1  = new TH1D("hst_Zne2_L3_BGO1", "Z O(alf1)",  7, 0.0, 7.0); hst_Zne2_L3_BGO1->Sumw2();

 hst_Zexp_DEL_SAT  = new TH1D("hst_Zexp_DEL_SAT", "Z expon. ",  7, 0.0, 7.0); hst_Zexp_DEL_SAT->Sumw2();
 hst_Zne1_DEL_SAT  = new TH1D("hst_Zne1_DEL_SAT", "Z O(alf1)",  7, 0.0, 7.0); hst_Zne1_DEL_SAT->Sumw2();
 hst_Zne2_DEL_SAT  = new TH1D("hst_Zne2_DEL_SAT", "Z O(alf1)",  7, 0.0, 7.0); hst_Zne2_DEL_SAT->Sumw2();

 hst_Zexp_SICAL92  = new TH1D("hst_Zexp_SICAL92", "Z expon. ", 7, 0.0, 7.0); hst_Zexp_SICAL92->Sumw2();
 hst_Zne1_SICAL92  = new TH1D("hst_Zne1_SICAL92", "Z O(alf1)", 7, 0.0, 7.0); hst_Zne1_SICAL92->Sumw2();
 hst_Zne2_SICAL92  = new TH1D("hst_Zne2_SICAL92", "Z O(alf1)", 7, 0.0, 7.0); hst_Zne2_SICAL92->Sumw2();

 hst_Zexp_SICAL94  = new TH1D("hst_Zexp_SICAL94", "Z expon. ", 7, 0.0, 7.0); hst_Zexp_SICAL94->Sumw2();
 hst_Zne1_SICAL94  = new TH1D("hst_Zne1_SICAL94", "Z O(alf1)", 7, 0.0, 7.0); hst_Zne1_SICAL94->Sumw2();
 hst_Zne2_SICAL94  = new TH1D("hst_Zne2_SICAL94", "Z O(alf1)", 7, 0.0, 7.0); hst_Zne2_SICAL94->Sumw2();
 //
 TH1D *hst2,*hst3 ,*hst4, *hst_ZZ1, *hst_ZZ2, *hst_ZZ4, *hst_ZZ3;
 int IBexp=1; double xsec;
 for( int IBexp=1; IBexp<= 20; IBexp++ ) {
 if(         IBexp ==1 ) { hst2=hst_Zexp_Al_LCAL; hst3=hst_Zne1_Al_LCAL; hst4=hst_Zne2_Al_LCAL; // Aleph LCAL
 } else if ( IBexp ==2 ) { hst2=hst_Zexp_Op_FCAL; hst3=hst_Zne1_Op_FCAL; hst4=hst_Zne2_Op_FCAL; // OPAL FCAL
 } else if ( IBexp ==3 ) { hst2=hst_Zexp_L3_BGO1; hst3=hst_Zne1_L3_BGO1; hst4=hst_Zne2_L3_BGO1; // L3 BGO
 } else if ( IBexp ==4 ) { hst2=hst_Zexp_DEL_SAT; hst3=hst_Zne1_DEL_SAT; hst4=hst_Zne2_DEL_SAT; // DELHI SAT
 } else if ( IBexp ==5 ) { hst2=hst_Zexp_SICAL92; hst3=hst_Zne1_SICAL92; hst4=hst_Zne2_SICAL92; // SICAL 92
 } else if ( IBexp ==12) { hst2=hst_Zexp_SICAL94; hst3=hst_Zne1_SICAL94; hst4=hst_Zne2_SICAL94; // SICAL 94
                  } else { hst2 = NULL; } // if (IBexp)
 if ( hst2 != NULL) {
 for( int j=1; j<= 7; j++ ) {
  if(         j ==1 ) { hst_ZZ1 = hst_ZZ1_88; hst_ZZ2 = hst_ZZ2_88; hst_ZZ4=hst_ZZ4_88; hst_ZZ3=hst_ZZ3_88;
  } else if ( j ==2 ) { hst_ZZ1 = hst_ZZ1_89; hst_ZZ2 = hst_ZZ2_89; hst_ZZ4=hst_ZZ4_89; hst_ZZ3=hst_ZZ3_89;
  } else if ( j ==3 ) { hst_ZZ1 = hst_ZZ1_90; hst_ZZ2 = hst_ZZ2_90; hst_ZZ4=hst_ZZ4_90; hst_ZZ3=hst_ZZ3_90;
  } else if ( j ==4 ) { hst_ZZ1 = hst_ZZ1_91; hst_ZZ2 = hst_ZZ2_91; hst_ZZ4=hst_ZZ4_91; hst_ZZ3=hst_ZZ3_91;
  } else if ( j ==5 ) { hst_ZZ1 = hst_ZZ1_92; hst_ZZ2 = hst_ZZ2_92; hst_ZZ4=hst_ZZ4_92; hst_ZZ3=hst_ZZ3_92;
  } else if ( j ==6 ) { hst_ZZ1 = hst_ZZ1_93; hst_ZZ2 = hst_ZZ2_93; hst_ZZ4=hst_ZZ4_93; hst_ZZ3=hst_ZZ3_93;
  } else if ( j ==7 ) { hst_ZZ1 = hst_ZZ1_94; hst_ZZ2 = hst_ZZ2_94; hst_ZZ4=hst_ZZ4_94; hst_ZZ3=hst_ZZ3_94;
  }// if
  xsec = hst_ZZ1->GetBinContent(IBexp);
  hst2->SetBinContent(j, hst_ZZ2->GetBinContent(IBexp)/xsec*10000 );
  hst2->SetBinError(  j, hst_ZZ2->GetBinError(IBexp)  /xsec*10000 );
  hst3->SetBinContent(j, hst_ZZ3->GetBinContent(IBexp)/xsec*10000 );
  hst3->SetBinError(  j, hst_ZZ3->GetBinError(IBexp)  /xsec*10000 );
  hst4->SetBinContent(j, hst_ZZ4->GetBinContent(IBexp)/xsec*10000 );
  hst4->SetBinError(  j, hst_ZZ4->GetBinError(IBexp)  /xsec*10000 );
 } // for j
 } // if(hst)
 } // for i

//////////////////////////////////////////////////////
// Column captions
 int nPlt=4;   // No of histos/column +1
 Char_t *Capt[nPlt+1]; // column captions
 for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
// pointers to histograms
 TH1D *iHst[nPlt+1];
// multicolumn caption
 Char_t Mcapt[132];
// format for columns
 char *fmt0 = "$  %10.2f $";               // 0-th column
 char *fmt1 = "& $ %10.4f \\pm %8.4f $ ";  // 1-st column
 char *fmt2 = "& $ %10.4f \\pm %8.4f $ ";  // next columns
 char *fmt3 = "& $ %10.4f $ ";             // error=0 case
//////////////////////////////////////////////////////

//**************************************************
  TString TeXfile = "TabZet1";  // Latex source file
  FILE *DiskFileTeX = fopen(TeXfile+".txp","w");
  LibPLT.PlInit(DiskFileTeX, 2); // Initialization of the latex source file
  fmt0 = "  $ %10.0f $";             // 0-th column
  fmt1 = "& $ %10.2f \\pm %8.2f $ "; // 1-st column
  fmt2 = "& $ %10.2f \\pm %8.2f $ "; // next columns
  fmt3 = "& $ %10.3f $ ";            // error=0 case
  LibPLT.Setfmt0(fmt0);
  LibPLT.Setfmt1(fmt1);
  LibPLT.Setfmt2(fmt2);
  LibPLT.Setfmt3(fmt3);
  iHst[1]= gst_ENE;     // beam energy
  iHst[2]= hst_Zexp_Al_LCAL;  // Z expon.
  iHst[3]= hst_Zne1_Al_LCAL;  // Z alf1
  iHst[4]= hst_Zne2_Al_LCAL;  // Z alf1
  strcpy(Capt[0],"{\\color{blue} No.}");
  strcpy(Capt[1],"{\\color{blue} $s^{1/2}$  [GeV]}");
  strcpy(Capt[2],"{\\color{blue} Z expon. }");
  strcpy(Capt[3],"{\\color{blue} Z ${\\cal O}(\\alpha^1)$ }");
  strcpy(Capt[4],"{\\color{blue} Z ${\\cal O}(\\alpha^1)$ }");
  strcpy(Mcapt,  "{\\color{red} ALEP LCAL: relat. Z contrib. $\\times~ 10^4 $ }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, " ",1,7, 1); // k1,k2,dk
  LibPLT.PlEnd(DiskFileTeX);  // finalizing latex source file
  fclose(DiskFileTeX);
//**************************************************

//**************************************************
  TeXfile = "TabZet2";  // Latex source file
  DiskFileTeX = fopen(TeXfile+".txp","w");
  LibPLT.PlInit(DiskFileTeX, 2); // Initialization of the latex source file
  fmt0 = "  $ %10.0f $";             // 0-th column
  fmt1 = "& $ %10.2f \\pm %8.2f $ "; // 1-st column
  fmt2 = "& $ %10.2f \\pm %8.2f $ "; // next columns
  fmt3 = "& $ %10.3f $ ";            // error=0 case
  LibPLT.Setfmt0(fmt0);
  LibPLT.Setfmt1(fmt1);
  LibPLT.Setfmt2(fmt2);
  LibPLT.Setfmt3(fmt3);
  iHst[1]= gst_ENE;     // beam energy
  iHst[2]= hst_Zexp_SICAL94;  // Z expon.
  iHst[3]= hst_Zne1_SICAL94;  // Z alf1
  iHst[4]= hst_Zne2_SICAL94;  // Z alf1
  strcpy(Capt[0],"{\\color{blue} No.}");
  strcpy(Capt[1],"{\\color{blue} $s^{1/2}$  [GeV]}");
  strcpy(Capt[2],"{\\color{blue} Z expon. }");
  strcpy(Capt[3],"{\\color{blue} Z ${\\cal O}(\\alpha^1)$ }");
  strcpy(Capt[4],"{\\color{blue} Z ${\\cal O}(\\alpha^1)$ }");
  strcpy(Mcapt,  "{\\color{red} ALEP SICAL: relat. Z contrib. $\\times~ 10^4 $ }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, " ",1,7, 1); // k1,k2,dk
  LibPLT.PlEnd(DiskFileTeX);  // finalizing latex source file
  fclose(DiskFileTeX);
  //**************************************************

  cout<<" ========================= TabZet end =========================== "<<endl;
}// TabZet



///////////////////////////////////////////////////////////////////////////////////
void TabZdel()
{
//************************************
  cout<<" ========================= TabZdel start ========================= "<<endl;
//
 TH1D *hst_Zdel2_Al_LCAL = HstDiff("hst_Zdel2_Al_LCAL", hst_Zexp_Al_LCAL, hst_Zne2_Al_LCAL, kRed );
 TH1D *hst_Zdel2_Op_FCAL = HstDiff("hst_Zdel2_Op_FCAL", hst_Zexp_Op_FCAL, hst_Zne2_Op_FCAL, kRed );
 TH1D *hst_Zdel2_L3_BGO1 = HstDiff("hst_Zdel2_L3_BGO1", hst_Zexp_L3_BGO1, hst_Zne2_L3_BGO1, kRed );
 TH1D *hst_Zdel2_DEL_SAT = HstDiff("hst_Zdel2_DEL_SAT", hst_Zexp_DEL_SAT, hst_Zne2_DEL_SAT, kRed );
 TH1D *hst_Zdel2_SICAL92 = HstDiff("hst_Zdel2_SICAL92", hst_Zexp_SICAL92, hst_Zne2_SICAL92, kRed );

//////////////////////////////////////////////////////
// Column captions
 int nPlt=6;   // No of histos/column +1
 Char_t *Capt[nPlt+1]; // column captions
 for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
// pointers to histograms
 TH1D *iHst[nPlt+1];
// multicolumn caption
 Char_t Mcapt[132];
// format for columns
 char *fmt0 = "$  %10.2f $";               // 0-th column
 char *fmt1 = "& $ %10.4f \\pm %8.4f $ ";  // 1-st column
 char *fmt2 = "& $ %10.4f \\pm %8.4f $ ";  // next columns
 char *fmt3 = "& $ %10.4f $ ";             // error=0 case
//////////////////////////////////////////////////////

//**************************************************
  TString TeXfile = "TabZdel";  // Latex source file
  FILE *DiskFileTeX = fopen(TeXfile+".txp","w");
  LibPLT.PlInit(DiskFileTeX, 2); // Initialization of the latex source file
  fmt0 = "  $ %10.0f $";             // 0-th column
  fmt1 = "& $ %10.2f \\pm %8.2f $ "; // 1-st column
  fmt2 = "& $ %10.2f \\pm %8.2f $ "; // next columns
  fmt3 = "& $ %10.3f $ ";            // error=0 case
  LibPLT.Setfmt0(fmt0);
  LibPLT.Setfmt1(fmt1);
  LibPLT.Setfmt2(fmt2);
  LibPLT.Setfmt3(fmt3);
  iHst[1]= gst_ENE;     // beam energy
  iHst[2]= hst_Zdel2_Al_LCAL;  //
  iHst[3]= hst_Zdel2_Op_FCAL;  //
  iHst[4]= hst_Zdel2_L3_BGO1;  //
  iHst[5]= hst_Zdel2_DEL_SAT;  //
  iHst[6]= hst_Zdel2_SICAL92;  //
  strcpy(Capt[0],"{\\color{blue} No.}");
  strcpy(Capt[1],"{\\color{blue} $s^{1/2}$  [GeV]}");
  strcpy(Capt[2],"{\\color{blue} ALEPH LCAL }");
  strcpy(Capt[3],"{\\color{blue} OPAL FCAL }");
  strcpy(Capt[4],"{\\color{blue} L3 BGO1 }");
  strcpy(Capt[5],"{\\color{blue} DELPHI SAT }");
  strcpy(Capt[6],"{\\color{blue} SICAL 92 }");
  strcpy(Mcapt,  "{\\color{red} Z contrib: $[ {\\cal O}(\\alpha^1)_{exp}-{\\cal O}(\\alpha^1) ] \\times~ 10^4 $ }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, " ",1,7, 1); // k1,k2,dk
  LibPLT.PlEnd(DiskFileTeX);  // finalizing latex source file
  fclose(DiskFileTeX);
//**************************************************

  cout<<" ========================= TabZdel end =========================== "<<endl;
}// TabZdel



///////////////////////////////////////////////////////////////////////////////////////
void Plt_LCAL1()
{
//------------------------------------------------------------------------
  cout<<" ========================= Plt_LCAL1 =========================== "<<endl;
//
 double del_Zexp_Al_LCAL[7], der_Zexp_Al_LCAL[7];
 double del_Zne1_Al_LCAL[7], der_Zne1_Al_LCAL[7];
 double del_Zne2_Al_LCAL[7], der_Zne2_Al_LCAL[7];
 double del_Zde1_Al_LCAL[7], der_Zde1_Al_LCAL[7];
 double del_Zde2_Al_LCAL[7], der_Zde2_Al_LCAL[7];

 TH1D *hst_Zde1_Al_LCAL = HstDiff("hst_Zde1_Al_LCAL", hst_Zexp_Al_LCAL, hst_Zne1_Al_LCAL, kBlack );
 TH1D *hst_Zde2_Al_LCAL = HstDiff("hst_Zde2_Al_LCAL", hst_Zexp_Al_LCAL, hst_Zne2_Al_LCAL, kRed );
 //
 double del_Zexp_SICAL94[7], der_Zexp_SICAL94[7];
 for(int i=1; i<=7; i++) {
  del_Zexp_Al_LCAL[i-1] = hst_Zexp_Al_LCAL->GetBinContent(i);
  der_Zexp_Al_LCAL[i-1] = hst_Zexp_Al_LCAL->GetBinError(i);
  del_Zne1_Al_LCAL[i-1] = hst_Zne1_Al_LCAL->GetBinContent(i);
  der_Zne1_Al_LCAL[i-1] = hst_Zne1_Al_LCAL->GetBinError(i);
  del_Zne2_Al_LCAL[i-1] = hst_Zne2_Al_LCAL->GetBinContent(i);
  der_Zne2_Al_LCAL[i-1] = hst_Zne2_Al_LCAL->GetBinError(i);
  //
  del_Zde1_Al_LCAL[i-1] = hst_Zde1_Al_LCAL->GetBinContent(i);
  der_Zde1_Al_LCAL[i-1] = hst_Zde1_Al_LCAL->GetBinError(i);
  del_Zde2_Al_LCAL[i-1] = hst_Zde2_Al_LCAL->GetBinContent(i);
  der_Zde2_Al_LCAL[i-1] = hst_Zde2_Al_LCAL->GetBinError(i);
 }
  double zero[7] = {0, 0, 0, 0, 0, 0, 0};
//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  //CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.030);
  CaptT->SetTextColor(kRed);
  //
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cPlt_LCAL1 = new TCanvas("cPlt_LCAL1","cPlt_LCAL1", gXcanv,  gYcanv,   600,  600);
  //                               Name    Title           xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  double x1=88, x2=94, y1=  -150.0, y2= 150.0;
  cPlt_LCAL1->DrawFrame(x1,y1, x2, y2); // two corners x1,y1, x2,y2
  gXcanv += 50; gYcanv += 50;
  cPlt_LCAL1->SetFillColor(10);
//////////////////////////////////////////////
  cPlt_LCAL1->cd();
//
  TGraphErrors *graph1 = new TGraphErrors(7, gEne, del_Zexp_Al_LCAL, zero, der_Zexp_Al_LCAL);
  graph1->SetMarkerStyle(25);
  graph1->SetMarkerColor(kRed);
  graph1->SetLineColor(kRed);
  graph1->Draw("CP");
  CaptT->SetTextColor(kRed);
  CaptT->DrawLatex(88.5, 100,"ALEPH LCAL expon.");
//
  TGraphErrors *graph2 = new TGraphErrors(7, gEne, del_Zne1_Al_LCAL, zero, der_Zne1_Al_LCAL);
  graph2->SetMarkerStyle(25);
  graph2->SetMarkerColor(kBlue);
  graph2->SetLineColor(kBlue);
  graph2->Draw("CP");
  CaptT->SetTextColor(kBlue);
  CaptT->DrawLatex(91.5,100,"ALEPH LCAL O(alf1) [82]");
//
  TGraphErrors *graph3 = new TGraphErrors(7, gEne, del_Zne2_Al_LCAL, zero, der_Zne2_Al_LCAL);
  graph3->SetMarkerStyle(25);
  graph3->SetMarkerColor(kGreen);
  graph3->SetLineColor(kGreen);
  graph3->Draw("CP");
  CaptT->SetTextColor(kGreen);
  CaptT->DrawLatex(91.5, 80,"ALEPH LCAL O(alf1) [92]");
  //////////////////////////////////////////////
  cPlt_LCAL1->cd();
  cPlt_LCAL1->SaveAs("cPlt_LCAL1.pdf");

  //
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cPlt_LCAL2 = new TCanvas("cPlt_LCAL2","cPlt_LCAL2", gXcanv,  gYcanv,   600,  600);
  //                               Name    Title           xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  x1=88, x2=94, y1=  -15.0, y2= 15.0;
  cPlt_LCAL2->DrawFrame(x1,y1, x2, y2); // two corners x1,y1, x2,y2
  gXcanv += 50; gYcanv += 50;
  cPlt_LCAL2->SetFillColor(10);
//////////////////////////////////////////////
  cPlt_LCAL2->cd();
//
  CaptT->DrawLatex(88.5, 12," Inset in Fig. 4 in PLB353(1995)349");
  TGraphErrors *graph10 = new TGraphErrors(7, gEne, del_Zde1_Al_LCAL, zero, der_Zde1_Al_LCAL);
  graph10->SetMarkerStyle(25);
  graph10->SetMarkerColor(kRed);
  graph10->SetLineColor(kRed);
  graph10->Draw("CP");
  CaptT->SetTextColor(kRed);
  CaptT->DrawLatex(88.5, 10,"ALEPH LCAL O(alf1)exp-O(alf1) [82]");
//
  TGraphErrors *graph12 = new TGraphErrors(7, gEne, del_Zde2_Al_LCAL, zero, der_Zde2_Al_LCAL);
  graph12->SetMarkerStyle(25);
  graph12->SetMarkerColor(kBlue);
  graph12->SetLineColor(kBlue);
  graph12->Draw("CP");
  CaptT->SetTextColor(kBlue);
  CaptT->DrawLatex(88.5, 8,"ALEPH LCAL O(alf1)exp-O(alf1) [92]");
  //////////////////////////////////////////////
  cPlt_LCAL2->cd();
  cPlt_LCAL2->SaveAs("cPlt_LCAL2.pdf");

}//Plt_LCAL1


///////////////////////////////////////////////////////////////////////////////////////
void Plt_SICAL1()
{
//------------------------------------------------------------------------
  cout<<" ========================= Plt_SICAL1 =========================== "<<endl;
//
 double del_Zexp_SICAL94[7], der_Zexp_SICAL94[7];
 double del_Zne1_SICAL94[7], der_Zne1_SICAL94[7];
 double del_Zne2_SICAL94[7], der_Zne2_SICAL94[7];
 double del_Zde1_SICAL94[7], der_Zde1_SICAL94[7];
 double del_Zde2_SICAL94[7], der_Zde2_SICAL94[7];

 TH1D *hst_Zde1_SICAL94 = HstDiff("hst_Zde1_SICAL94", hst_Zexp_SICAL94, hst_Zne1_SICAL94, kBlack );
 TH1D *hst_Zde2_SICAL94 = HstDiff("hst_Zde2_SICAL94", hst_Zexp_SICAL94, hst_Zne2_SICAL94, kRed );
 //
 for(int i=1; i<=7; i++) {
  del_Zexp_SICAL94[i-1] = hst_Zexp_SICAL94->GetBinContent(i);
  der_Zexp_SICAL94[i-1] = hst_Zexp_SICAL94->GetBinError(i);
  del_Zne1_SICAL94[i-1] = hst_Zne1_SICAL94->GetBinContent(i);
  der_Zne1_SICAL94[i-1] = hst_Zne1_SICAL94->GetBinError(i);
  del_Zne2_SICAL94[i-1] = hst_Zne2_SICAL94->GetBinContent(i);
  der_Zne2_SICAL94[i-1] = hst_Zne2_SICAL94->GetBinError(i);
  //
  del_Zde1_SICAL94[i-1] = hst_Zde1_SICAL94->GetBinContent(i);
  der_Zde1_SICAL94[i-1] = hst_Zde1_SICAL94->GetBinError(i);
  del_Zde2_SICAL94[i-1] = hst_Zde2_SICAL94->GetBinContent(i);
  der_Zde2_SICAL94[i-1] = hst_Zde2_SICAL94->GetBinError(i);
 }
  double zero[7] = {0, 0, 0, 0, 0, 0, 0};
//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  //CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.030);
  CaptT->SetTextColor(kRed);
  //
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cPlt_SICAL1 = new TCanvas("cPlt_SICAL1","cPlt_SICAL1", gXcanv,  gYcanv,   600,  600);
  //                               Name    Title           xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  double x1=88, x2=94, y1=  -30.0, y2= 30.0;
  cPlt_SICAL1->DrawFrame(x1,y1, x2, y2); // two corners x1,y1, x2,y2
  gXcanv += 50; gYcanv += 50;
  cPlt_SICAL1->SetFillColor(10);
//////////////////////////////////////////////
  cPlt_SICAL1->cd();
//
  TGraphErrors *graph1 = new TGraphErrors(7, gEne, del_Zexp_SICAL94, zero, der_Zexp_SICAL94);
  graph1->SetMarkerStyle(25);
  graph1->SetMarkerColor(kRed);
  graph1->SetLineColor(kRed);
  graph1->Draw("CP");
  CaptT->SetTextColor(kRed);
  CaptT->DrawLatex(88.5, 25,"ALEPH SICAL expon.");
//
  TGraphErrors *graph2 = new TGraphErrors(7, gEne, del_Zne1_SICAL94, zero, der_Zne1_SICAL94);
  graph2->SetMarkerStyle(25);
  graph2->SetMarkerColor(kBlue);
  graph2->SetLineColor(kBlue);
  graph2->Draw("CP");
  CaptT->SetTextColor(kBlue);
  CaptT->DrawLatex(91.5, 25,"ALEPH SICAL O(alf1) [82]");
//
  TGraphErrors *graph3 = new TGraphErrors(7, gEne, del_Zne2_SICAL94, zero, der_Zne2_SICAL94);
  graph3->SetMarkerStyle(25);
  graph3->SetMarkerColor(kGreen);
  graph3->SetLineColor(kGreen);
  graph3->Draw("CP");
  CaptT->SetTextColor(kGreen);
  CaptT->DrawLatex(91.5, 20,"ALEPH SICAL O(alf1) [92]");
  //////////////////////////////////////////////
  cPlt_SICAL1->cd();
  cPlt_SICAL1->SaveAs("cPlt_SICAL1.pdf");

  //
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cPlt_SICAL2 = new TCanvas("cPlt_SICAL2","cPlt_SICAL2", gXcanv,  gYcanv,   600,  600);
  //                               Name    Title           xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  x1=88, x2=94, y1=  -4.0, y2= 2.0;
  cPlt_SICAL2->DrawFrame(x1,y1, x2, y2); // two corners x1,y1, x2,y2
  gXcanv += 50; gYcanv += 50;
  cPlt_SICAL2->SetFillColor(10);
//////////////////////////////////////////////
  cPlt_SICAL2->cd();
//
  CaptT->DrawLatex(88.5, 1.8," Inset in Fig. 5 in PLB353(1995)349");
  TGraphErrors *graph10 = new TGraphErrors(7, gEne, del_Zde1_SICAL94, zero, der_Zde1_SICAL94);
  graph10->SetMarkerStyle(25);
  graph10->SetMarkerColor(kRed);
  graph10->SetLineColor(kRed);
  graph10->Draw("CP");
  CaptT->SetTextColor(kRed);
  CaptT->DrawLatex(88.5, 1.4,"ALEPH LCAL O(alf1)exp-O(alf1) [82]");
//
  TGraphErrors *graph12 = new TGraphErrors(7, gEne, del_Zde2_SICAL94, zero, der_Zde2_SICAL94);
  graph12->SetMarkerStyle(25);
  graph12->SetMarkerColor(kBlue);
  graph12->SetLineColor(kBlue);
  graph12->Draw("CP");
  CaptT->SetTextColor(kBlue);
  CaptT->DrawLatex(88.5, 1.0,"ALEPH LCAL O(alf1)exp-O(alf1) [92]");
  //////////////////////////////////////////////
  cPlt_SICAL2->cd();
  cPlt_SICAL2->SaveAs("cPlt_SICAL2.pdf");

}//Plt_SICAL1



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
// libProd.so should be loaded prior to opening GenFile
  gSystem->Load("../ProdCxx/.libs/libProd.so"); // seems to be needed !!!
  DiskFileB.cd();
//
  LibPLT.Initialize(); // latex tables
//
  TH1D *HST_BHL_NORMA = (TH1D*)gHstFile91->Get("HST_BHL_NORMA");
  gBHLgen= (TMCgenBHL6*)gGenFile91->Get("MCgen");
  gBHLgen->ReInitialize(&OutFile, HST_BHL_NORMA);
  gCMSENE  = gBHLgen->m_CMSENE;
  cout<< " Plot2:  gCMSENE="<<gCMSENE<<endl;

  gst_ENE   = new TH1D("gst_ENE",   "Energy list    ", 7, 0.0, 7.0);
  gst_ENE->Sumw2();
  // list of energies
  for( int ib=1; ib<=7; ib++) gst_ENE->SetBinContent(ib, gEne[ib-1]);

  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++

  HistRemake(gHstFile88);     // Renormalization of MC histograms
  HistRemake(gHstFile89);     // Renormalization of MC histograms
  HistRemake(gHstFile90);     // Renormalization of MC histograms
  HistRemake(gHstFile91);     // Renormalization of MC histograms
  HistRemake(gHstFile92);     // Renormalization of MC histograms
  HistRemake(gHstFile93);     // Renormalization of MC histograms
  HistRemake(gHstFile94);     // Renormalization of MC histograms

  //========== Tables ==========

  TabTrig(gHstFile88, "TabZ_88GeV");
  TabTrig(gHstFile89, "TabZ_89GeV");
  TabTrig(gHstFile90, "TabZ_90GeV");
  TabTrig(gHstFile91, "TabZ_91GeV");
  TabTrig(gHstFile92, "TabZ_92GeV");
  TabTrig(gHstFile93, "TabZ_93GeV");
  TabTrig(gHstFile94, "TabZ_94GeV");

  TabZet();
  TabZdel();

 //========== Figures =========
  //Graph1();

  Plt_LCAL1();
  Plt_SICAL1();
  //++++++++++++++++++++++++++++++++++++++++

//  gHstFile91->ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  OutFile.close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}



