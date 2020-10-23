//    make Plot2-run
//    make NewVPs-pdf

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

TH1D *hst_dVP_Al_LCAL, *hst_dVP_Op_FCAL, *hst_dVP_L3_BGO1,  *hst_dVP_DEL_SAT;
TH1D *hst_dVP_Op_OSIW, *hst_dVP_L3_BGO2, *hst_dVP_DEL_STIC, *hst_dVP_DSAT_93;
TH1D *hst_dVP_SICAL92, *hst_dVP_SICAL93, *hst_dVP_SICAL94,  *hst_dVP_SICAL95;

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


HSTplot LibPLT("LibPLT");

///////////////////////////////////////////////////////////////////////////////////
void PlotSame(TH1D *HST, double &ycapt, Int_t kolor, TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, opis);
}// PlotSame


///////////////////////////////////////////////////////////////////////////////////
void PlotSame2(TH1D *HST, double &ycapt, Int_t kolor, double xx,  TString label,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
  double xcapt = 0.40;
  CaptT->DrawLatex(xcapt,ycapt, opis);
  CaptT->DrawLatex(xcapt-0.05,ycapt, label);
  //
  TLatex *CaptS = new TLatex();
  CaptS->SetTextSize(0.040);
  CaptS->SetTextAlign(21);
  CaptS->SetTextColor(kolor);
  int ib = HST->FindBin(xx);
  double yy= HST->GetBinContent(ib);
  CaptS->DrawLatex(xx,yy,label);
}// PlotSame2

///////////////////////////////////////////////////////////////////////////////////
void HistRemake(TFile *HstFile ){
//
cout<<"----------------------------- HistRemake ------------------------------------"<<endl;

TH1D *HST_BHL_NORMA = (TH1D*)HstFile->Get("HST_BHL_NORMA");
///////////////////////////////////////////////////////////////////////////
// renormalizing 1-dim histos
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_VP0") );
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_VP1") );
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_VP2") );
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_VP3") );
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_VP4") );
    HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_trig_VP5") );

    cout<<"----------------------------- HistRemake end ----------------------------------"<<endl;
}//HistRemake


///////////////////////////////////////////////////////////////////////////////////
void TabTrig(TFile *HstFile, TString TeXfile)
{
//------------------------------------------------------------------------
  cout<<" ========================= TabTrig start=========================== "<<endl;
  //
  TH1D *hst_trig_VP0  = (TH1D*)HstFile->Get("hst_trig_VP0");
  TH1D *hst_trig_VP1  = (TH1D*)HstFile->Get("hst_trig_VP1");
  TH1D *hst_trig_VP2  = (TH1D*)HstFile->Get("hst_trig_VP2");
  TH1D *hst_trig_VP3  = (TH1D*)HstFile->Get("hst_trig_VP3");
  TH1D *hst_trig_VP4  = (TH1D*)HstFile->Get("hst_trig_VP4");
  TH1D *hst_trig_VP5  = (TH1D*)HstFile->Get("hst_trig_VP5");

  for(int i=0; i<=100; i++ ) hst_trig_VP4->SetBinError(i,0); // errors omitted
  for(int i=0; i<=100; i++ ) hst_trig_VP2->SetBinError(i,0); // errors omitted
  for(int i=0; i<=100; i++ ) hst_trig_VP5->SetBinError(i,0); // errors omitted


// Column captions
  int nPlt=6;   // No of histos/column +1
  Char_t *Capt[nPlt+1];
  for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
//
  strcpy(Capt[0],"{\\color{blue} No.}");
  strcpy(Capt[1],"{\\color{blue} VP off }");
  strcpy(Capt[2],"{\\color{blue} VPold }");
  strcpy(Capt[3],"{\\color{blue} VPerr }");
  strcpy(Capt[4],"{\\color{blue} VPnew }");
  strcpy(Capt[5],"{\\color{blue} VPerr }");
  strcpy(Capt[6],"{\\color{blue} new-old }");

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

  iHst[1]= hst_trig_VP0;  // VP off
  iHst[2]= hst_trig_VP1;  // VPold on
  iHst[3]= hst_trig_VP2;  // VPold_err
  iHst[4]= hst_trig_VP3;  // VPnew
  iHst[5]= hst_trig_VP4;  // VPnew err
  iHst[6]= hst_trig_VP5;  // VP new-old

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
void TabDiff()
{
//************************************
  cout<<" ========================= TabDiff start ========================= "<<endl;
  //
 TH1D *hst_VP5_88  = (TH1D*)gHstFile88->Get("hst_trig_VP5");
 TH1D *hst_VP5_89  = (TH1D*)gHstFile89->Get("hst_trig_VP5");
 TH1D *hst_VP5_90  = (TH1D*)gHstFile90->Get("hst_trig_VP5");
 TH1D *hst_VP5_91  = (TH1D*)gHstFile91->Get("hst_trig_VP5");
 TH1D *hst_VP5_92  = (TH1D*)gHstFile92->Get("hst_trig_VP5");
 TH1D *hst_VP5_93  = (TH1D*)gHstFile93->Get("hst_trig_VP5");
 TH1D *hst_VP5_94  = (TH1D*)gHstFile94->Get("hst_trig_VP5");
 //
 TH1D *hst_VP4_88  = (TH1D*)gHstFile88->Get("hst_trig_VP4");
 TH1D *hst_VP4_89  = (TH1D*)gHstFile89->Get("hst_trig_VP4");
 TH1D *hst_VP4_90  = (TH1D*)gHstFile90->Get("hst_trig_VP4");
 TH1D *hst_VP4_91  = (TH1D*)gHstFile91->Get("hst_trig_VP4");
 TH1D *hst_VP4_92  = (TH1D*)gHstFile92->Get("hst_trig_VP4");
 TH1D *hst_VP4_93  = (TH1D*)gHstFile93->Get("hst_trig_VP4");
 TH1D *hst_VP4_94  = (TH1D*)gHstFile94->Get("hst_trig_VP4");
 //
 TH1D *hst_VP3_88  = (TH1D*)gHstFile88->Get("hst_trig_VP3");
 TH1D *hst_VP3_89  = (TH1D*)gHstFile89->Get("hst_trig_VP3");
 TH1D *hst_VP3_90  = (TH1D*)gHstFile90->Get("hst_trig_VP3");
 TH1D *hst_VP3_91  = (TH1D*)gHstFile91->Get("hst_trig_VP3");
 TH1D *hst_VP3_92  = (TH1D*)gHstFile92->Get("hst_trig_VP3");
 TH1D *hst_VP3_93  = (TH1D*)gHstFile93->Get("hst_trig_VP3");
 TH1D *hst_VP3_94  = (TH1D*)gHstFile94->Get("hst_trig_VP3");

 hst_dVP_Al_LCAL  = new TH1D("hst_dVP_Al_LCAL", "VP differences ",  7, 0.0, 7.0); hst_dVP_Al_LCAL->Sumw2();
 hst_dVP_Op_FCAL  = new TH1D("hst_dVP_Op_FCAL", "VP differences ",  7, 0.0, 7.0); hst_dVP_Op_FCAL->Sumw2();
 hst_dVP_L3_BGO1  = new TH1D("hst_dVP_L3_BGO1", "VP differences ",  7, 0.0, 7.0); hst_dVP_L3_BGO1->Sumw2();
 hst_dVP_DEL_SAT  = new TH1D("hst_dVP_DEL_SAT", "VP differences ",  7, 0.0, 7.0); hst_dVP_DEL_SAT->Sumw2();
 hst_dVP_DSAT_93  = new TH1D("hst_dVP_DSAT_93", "VP differences ",  7, 0.0, 7.0); hst_dVP_DSAT_93->Sumw2();
 //
 hst_dVP_Op_OSIW  = new TH1D("hst_dVP_Op_OSIW",  "VP differences ", 7, 0.0, 7.0); hst_dVP_Op_OSIW->Sumw2();
 hst_dVP_L3_BGO2  = new TH1D("hst_dVP_L3_BGO2",  "VP differences ", 7, 0.0, 7.0); hst_dVP_L3_BGO2->Sumw2();
 hst_dVP_DEL_STIC = new TH1D("hst_dVP_DEL_STIC", "VP differences ", 7, 0.0, 7.0); hst_dVP_DEL_STIC->Sumw2();
 //
 hst_dVP_SICAL92  = new TH1D("hst_dVP_SICAL92", "VP differences ", 7, 0.0, 7.0); hst_dVP_SICAL92->Sumw2();
 hst_dVP_SICAL93  = new TH1D("hst_dVP_SICAL93", "VP differences ", 7, 0.0, 7.0); hst_dVP_SICAL93->Sumw2();
 hst_dVP_SICAL94  = new TH1D("hst_dVP_SICAL94", "VP differences ", 7, 0.0, 7.0); hst_dVP_SICAL94->Sumw2();
 hst_dVP_SICAL95  = new TH1D("hst_dVP_SICAL95", "VP differences ", 7, 0.0, 7.0); hst_dVP_SICAL95->Sumw2();
 //
 TH1D *hst, *hst_VP3, *hst_VP4, *hst_VP5;
 int IBexp=1; double xsec;
 for( int IBexp=1; IBexp<= 20; IBexp++ ) {
 if(         IBexp ==1 ) { hst= hst_dVP_Al_LCAL; // Aleph LCAL
 } else if ( IBexp ==2 ) { hst= hst_dVP_Op_FCAL; // OPAL FCAL
 } else if ( IBexp ==3 ) { hst= hst_dVP_L3_BGO1; // L3 BGO
 } else if ( IBexp ==4 ) { hst= hst_dVP_DEL_SAT; // DELHI SAT
 } else if ( IBexp ==14) { hst= hst_dVP_Op_OSIW;  // OPAL FCAL
 } else if ( IBexp ==15) { hst= hst_dVP_L3_BGO2;  // L3 BGO
 } else if ( IBexp ==16) { hst= hst_dVP_DEL_STIC; // DELHI STIC
 } else if ( IBexp ==17) { hst= hst_dVP_DSAT_93;  // DELHI SAT 93
 } else if ( IBexp == 5) { hst= hst_dVP_SICAL92;  // SICAL 92
 } else if ( IBexp ==11) { hst= hst_dVP_SICAL93;  // SICAL 93
 } else if ( IBexp ==12) { hst= hst_dVP_SICAL94;  // SICAL 94
 } else if ( IBexp ==13) { hst= hst_dVP_SICAL95;  // SICAL 95
                  } else { hst = NULL; } // if (IBexp)
 if ( hst != NULL) {
 for( int j=1; j<= 7; j++ ) {
  if(         j ==1 ) { hst_VP5 = hst_VP5_88; hst_VP3=hst_VP3_88; hst_VP4=hst_VP4_88;
  } else if ( j ==2 ) { hst_VP5 = hst_VP5_89; hst_VP3=hst_VP3_89; hst_VP4=hst_VP4_89;
  } else if ( j ==3 ) { hst_VP5 = hst_VP5_90; hst_VP3=hst_VP3_90; hst_VP4=hst_VP4_90;
  } else if ( j ==4 ) { hst_VP5 = hst_VP5_91; hst_VP3=hst_VP3_91; hst_VP4=hst_VP4_91;
  } else if ( j ==5 ) { hst_VP5 = hst_VP5_92; hst_VP3=hst_VP3_92; hst_VP4=hst_VP4_92;
  } else if ( j ==6 ) { hst_VP5 = hst_VP5_93; hst_VP3=hst_VP3_93; hst_VP4=hst_VP4_93;
  } else if ( j ==7 ) { hst_VP5 = hst_VP5_94; hst_VP3=hst_VP3_94; hst_VP4=hst_VP4_94;
  }// if
  xsec = hst_VP3->GetBinContent(IBexp);
  hst->SetBinContent(j, hst_VP5->GetBinContent(IBexp)/xsec*10000 );
  hst->SetBinError(  j, hst_VP4->GetBinContent(IBexp)/xsec*10000 );
 } // for j
 } // if(hst)
 } // for i

//////////////////////////////////////////////////////
// Column captions
 int nPlt=5;   // No of histos/column +1
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
  TString TeXfile = "TabDiff";  // Latex source file
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
  iHst[2]= hst_dVP_Al_LCAL;  // VP diff
  iHst[3]= hst_dVP_Op_FCAL;  // VP diff
  iHst[4]= hst_dVP_L3_BGO1;  // VP diff
  iHst[5]= hst_dVP_DEL_SAT;  // VP diff
  strcpy(Capt[0],"{\\color{blue} No.}");
  strcpy(Capt[1],"{\\color{blue} $s^{1/2}$  [GeV]}");
  strcpy(Capt[2],"{\\color{blue} ALEPH LCAL }");
  strcpy(Capt[3],"{\\color{blue} OPAL FCAL }");
  strcpy(Capt[4],"{\\color{blue} L3 BGO1 }");
  strcpy(Capt[5],"{\\color{blue} DELPHI SAT }");
  strcpy(Mcapt,  "{\\color{red} VP: (new-old)/new $\\times~ 10^4 $ }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, " ",1,7, 1); // k1,k2,dk
  LibPLT.PlEnd(DiskFileTeX);  // finalizing latex source file
  fclose(DiskFileTeX);
//**************************************************

//**************************************************
  TeXfile = "TabDiff2";  // Latex source file
  DiskFileTeX = fopen(TeXfile+".txp","w");
  LibPLT.PlInit(DiskFileTeX, 2); // Initialization of the latex source file
  iHst[1]= gst_ENE;     // beam energy
  iHst[2]= hst_dVP_DSAT_93;  // VP diff
  iHst[3]= hst_dVP_Op_OSIW;  // VP diff
  iHst[4]= hst_dVP_L3_BGO2;  // VP diff
  iHst[5]= hst_dVP_DEL_STIC; // VP diff
  strcpy(Capt[0],"{\\color{blue} No.}");
  strcpy(Capt[1],"{\\color{blue} $s^{1/2}$  [GeV]}");
  strcpy(Capt[2],"{\\color{blue} DELPHI SAT 93 }");
  strcpy(Capt[3],"{\\color{blue} OPAL OSIW }");
  strcpy(Capt[4],"{\\color{blue} L3 BGO2 }");
  strcpy(Capt[5],"{\\color{blue} DELPHI STIC }");
  strcpy(Mcapt,  "{\\color{red} VP: (new-old)/new $\\times~ 10^4 $ }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, " ",1,7, 1); // k1,k2,dk
  LibPLT.PlEnd(DiskFileTeX);  // finalizing latex source file
  fclose(DiskFileTeX);
//**************************************************


//**************************************************
  TeXfile = "TabDiff3";  // Latex source file
  DiskFileTeX = fopen(TeXfile+".txp","w");
  LibPLT.PlInit(DiskFileTeX, 2); // Initialization of the latex source file
  iHst[1]= gst_ENE;     // beam energy
  iHst[2]= hst_dVP_SICAL92;  // VP diff
  iHst[3]= hst_dVP_SICAL93;  // VP diff
  iHst[4]= hst_dVP_SICAL94;  // VP diff
  iHst[5]= hst_dVP_SICAL95;  // VP diff
  strcpy(Capt[0],"{\\color{blue} No.}");
  strcpy(Capt[1],"{\\color{blue} $s^{1/2}$  [GeV]}");
  strcpy(Capt[2],"{\\color{blue} ALEPH 92 }");
  strcpy(Capt[3],"{\\color{blue} ALEPH 93 }");
  strcpy(Capt[4],"{\\color{blue} ALEPH 94}");
  strcpy(Capt[5],"{\\color{blue} ALEPH 95 }");
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, " ",1,7, 1); // k1,k2,dk
  LibPLT.PlEnd(DiskFileTeX);  // finalizing latex source file
  fclose(DiskFileTeX);
//**************************************************

  cout<<" ========================= TabDiff end =========================== "<<endl;
}// TabDiff


TH1D *HstTrig(TString title, double yy, double th1w, double th2w, double th1n, double th2n)
// Graphic representation of lumi detector angular range
{
  // makes special histo
  double tfid1=0.018, tfid2=0.145;
  cout<<"Entering MakeCumul for  ";
  int NB=1000;
  TH1D* HST = new TH1D(title,"theta range", NB, tfid1, tfid2);
  HST->Sumw2();
  double theta, err;
  for( int i=1; i<=NB; i++ ) {
        HST->SetBinContent(i,yy);
        theta= tfid1+ (i*(tfid2-tfid1))/NB;
        err =0;
        if( theta >th1w && theta <th2w) err=0.10;
        if( theta >th1n && theta <th2n) err=0.20;
        HST->SetBinError(i,err);
  }
  return HST;
}//HstTrig

///////////////////////////////////////////////////////////////////////////////
void AngRng()
{
//------------------------------------------------------------------------
  cout<<" ========================= AngRng =========================== "<<endl;

  TH1D *h_AlepLcal  = HstTrig("h_AlepLcal",  1.0,   0.04300,  0.12500,   0.05700,   0.10700 );
  TH1D *h_OpalFcal  = HstTrig("h_OpalFcal",  2.0,   0.05500,  0.11500,   0.06500,   0.10500 );
  TH1D *h_L3BGO1    = HstTrig("h_L3BGO1",    3.0,   0.02520,  0.07120,   0.03120,   0.06520 );
  TH1D *h_DelphiSat = HstTrig("h_DelphiSat", 4.0,   0.05272,  0.14182,   0.05602,   0.12862 );
  //
  TH1D *h_Aleph_Sical= HstTrig("h_Aleph_Sical",  5.0,  0.0261250,  0.0558750,   0.0303750,   0.0495000);
  TH1D *h_Opal_OSIW  = HstTrig("h_Opal_OSIW",    6.0,  0.0272325,  0.0556875,   0.0312975,   0.0516225);
  TH1D *h_L3_BGO2    = HstTrig("h_L3_BGO2",      7.0,  0.0270000,  0.0650000,   0.0320000,   0.0540000);
  TH1D *h_Delphi_Stic= HstTrig("h_Delphi_Stic",  8.0,  0.0372000,  0.1268000,   0.0436000,   0.1132000);
  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  //CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  CaptT->SetTextColor(kRed);
  //
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAngRng = new TCanvas("cAngRng","cAngRng", gXcanv,  gYcanv,   800,  600);
  //                               Name    Title           xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAngRng->SetFillColor(10);
  //////////////////////////////////////////////
  cAngRng->cd(1);
  h_AlepLcal->SetStats(0);
  h_AlepLcal->SetTitle(0);
  h_AlepLcal->SetMinimum(0); h_AlepLcal->SetMaximum(9);
  h_AlepLcal->GetXaxis()->SetTitle("#theta [rad]");
  //
  h_AlepLcal->DrawCopy("h");
  h_OpalFcal->DrawCopy("hsame");
  h_L3BGO1->DrawCopy("hsame");
  h_DelphiSat->DrawCopy("hsame");
  //
  h_Aleph_Sical->DrawCopy("hsame");
  h_Opal_OSIW->DrawCopy("hsame");
  h_L3_BGO2->DrawCopy("hsame");
  h_Delphi_Stic->DrawCopy("hsame");

  CaptT->DrawLatex(0.022,1.2,"ALEP LCAL");
  CaptT->DrawLatex(0.022,2.2,"OPAL FCAL");
  CaptT->DrawLatex(0.100,3.2,"L3 BGO1");
  CaptT->DrawLatex(0.022,4.2,"DELHI SAT");
  //
  CaptT->DrawLatex(0.100,5.2,"ALEPH SICAL");
  CaptT->DrawLatex(0.100,6.2,"OPAL OSIW");
  CaptT->DrawLatex(0.100,7.2,"L3 BGO2");
  CaptT->DrawLatex(0.022,8.2,"DELPHI STIC");
  //
  //////////////////////////////////////////////
  //
  cAngRng->cd();
  cAngRng->SaveAs("cAngRng.pdf");
//
}// AngRng

void Graph1()
{
//------------------------------------------------------------------------
// !!!! ROOT tutorial !!!!
// Graphs are drawn via the painter TGraphPainter class.
// This class implements techniques needed to display the various kind of graphs
// i.e.: TGraph, TGraphErrors, TGraphBentErrors and TGraphAsymmErrors.
//
cout<<" ========================= Graph1 =========================== "<<endl;

   TCanvas *c43 = new TCanvas("c43","c43",10,10,800,600);
   c43->DrawFrame(0., -0.5, 6., 2);
   double x[5]    = {1, 2, 3, 4, 5};
   double zero[5] = {0, 0, 0, 0, 0};
   // data set (1) with stat and sys errors
   double py1[5]      = {1.2, 1.15, 1.19, 0.9, 1.4};
   double ey_stat1[5] = {0.2, 0.18, 0.17, 0.2, 0.4};
   double ey_sys1[5]  = {0.5, 0.71, 0.76, 0.5, 0.45};
   // data set (2) with stat and sys errors
   double y2[5]       = {0.25, 0.18, 0.29, 0.2, 0.21};
   double ey_stat2[5] = {0.2, 0.18, 0.17, 0.2, 0.4};
   double ey_sys2[5]  = {0.63, 0.19, 0.7, 0.2, 0.7};
   // Now draw data set (1)
   // We first have to draw it only with the stat errors
   TGraphErrors *graph1 = new TGraphErrors(5, x, py1, zero, ey_stat1);
   graph1->SetMarkerStyle(20);
   graph1->Draw("P");
   // Now we have to somehow depict the sys errors
   TGraphErrors *graph1_sys = new TGraphErrors(5, x, py1, zero, ey_sys1);
   graph1_sys->Draw("[]");
   // Now draw data set (2)
   // We first have to draw it only with the stat errors
   TGraphErrors *graph2 = new TGraphErrors(5, x, y2, zero, ey_stat2);
   graph2->SetMarkerStyle(24);
   graph2->Draw("P");
   // Now we have to somehow depict the sys errors
   TGraphErrors *graph2_sys = new TGraphErrors(5, x, y2, zero, ey_sys2);
   graph2_sys->Draw("[]");

}// Graph1


///////////////////////////////////////////////////////////////////////////////////////
void PltDif1()
{
//------------------------------------------------------------------------
  cout<<" ========================= PltDif1 =========================== "<<endl;
//
//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  //CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.030);
  CaptT->SetTextColor(kRed);
  //
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cPltDif1 = new TCanvas("cPltDif1","cPltDif1", gXcanv,  gYcanv,   800,  600);
  //                               Name    Title           xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  double x1=88, x2=94, y1=  -10.0, y2= 2.0;
  cPltDif1->DrawFrame(x1,y1, x2, y2); // two corners x1,y1, x2,y2
  gXcanv += 50; gYcanv += 50;
  cPltDif1->SetFillColor(10);
//////////////////////////////////////////////

  double del_dVP_Al_LCAL[7], del_dVP_Op_FCAL[7], del_dVP_L3_BGO1[7], del_dVP_DEL_SAT[7];
  double der_dVP_Al_LCAL[7], der_dVP_Op_FCAL[7], der_dVP_L3_BGO1[7], der_dVP_DEL_SAT[7];
  for(int i=1; i<=7; i++) {
	  del_dVP_Al_LCAL[i-1] = hst_dVP_Al_LCAL->GetBinContent(i);
	  del_dVP_Op_FCAL[i-1] = hst_dVP_Op_FCAL->GetBinContent(i);
	  del_dVP_L3_BGO1[i-1] = hst_dVP_L3_BGO1->GetBinContent(i);
	  del_dVP_DEL_SAT[i-1] = hst_dVP_DEL_SAT->GetBinContent(i);
	  //
	  der_dVP_Al_LCAL[i-1] = hst_dVP_Al_LCAL->GetBinError(i);
	  der_dVP_Op_FCAL[i-1] = hst_dVP_Op_FCAL->GetBinError(i);
	  der_dVP_L3_BGO1[i-1] = hst_dVP_L3_BGO1->GetBinError(i);
	  der_dVP_DEL_SAT[i-1] = hst_dVP_DEL_SAT->GetBinError(i);
  }
  double zero[7] = {0, 0, 0, 0, 0, 0, 0};

  cPltDif1->cd();

  TGraphErrors *graph1 = new TGraphErrors(7, gEne, del_dVP_Al_LCAL, zero, der_dVP_Al_LCAL);
  graph1->SetMarkerStyle(25);
  graph1->SetMarkerColor(kRed);
  graph1->SetLineColor(kRed);
  graph1->Draw("LP");
  CaptT->SetTextColor(kRed);
  CaptT->DrawLatex(88.5,-3.80,"ALEPH LCAL");
//
  TGraphErrors *graph2 = new TGraphErrors(7, gEne2, del_dVP_Op_FCAL, zero, der_dVP_Op_FCAL);
  graph2->SetMarkerStyle(20);
  graph2->SetLineColor(kBlue);
  graph2->SetMarkerColor(kBlue);
  graph2->Draw("LP");
  CaptT->SetTextColor(kBlue);
  CaptT->DrawLatex(88.5,-6.00,"OPAL FCAL");
//
  TGraphErrors *graph3 = new TGraphErrors(7, gEne, del_dVP_L3_BGO1, zero, der_dVP_L3_BGO1);
  graph3->SetMarkerStyle(23);
  graph3->SetLineColor(kGreen);
  graph3->SetMarkerColor(kGreen);
  graph3->Draw("LP");
  CaptT->SetTextColor(kGreen);
  CaptT->DrawLatex(88.5, 0.50,"L3 BGO1");
//
  TGraphErrors *graph4 = new TGraphErrors(7, gEne1, del_dVP_DEL_SAT, zero, der_dVP_DEL_SAT);
  graph4->SetMarkerStyle(24);
  graph4->SetLineColor(kMagenta);
  graph4->SetMarkerColor(kMagenta);
  graph4->Draw("LP");
  CaptT->SetTextColor(kMagenta);
  CaptT->DrawLatex(88.5,-2.50,"DELPHI SAT");
  //////////////////////////////////////////////
  cPltDif1->cd();
  cPltDif1->SaveAs("cPltDif1.pdf");

}//PltDif1



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

  TabTrig(gHstFile88, "Tab88GeV");
  TabTrig(gHstFile89, "Tab89GeV");
  TabTrig(gHstFile90, "Tab90GeV");
  TabTrig(gHstFile91, "Tab91GeV");
  TabTrig(gHstFile94, "Tab92GeV");
  TabTrig(gHstFile93, "Tab93GeV");
  TabTrig(gHstFile94, "Tab94GeV");

  TabDiff();

 //========== Figures =========
  //Graph1();

  AngRng();
  PltDif1();

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



