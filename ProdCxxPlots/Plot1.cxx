//    make Plot1-run

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

#include "TMCgenBHL6.h"
#include "TRobolProd.h"

#include "HSTplot.h"
#include "HisNorm.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
TString DiskFileA = "../ProdCxx/Prod2/histo.root"; // current
TString DiskFileG = "../ProdCxx/Prod2/mcgen.root"; // current
//
//TString DiskFileA = "../ProdCxx/Prod2/histo.root_VP2Zon_3G"; //
//TString DiskFileA = "../ProdCxx/Prod2/histo.root_VPZoff_5G"; //
//TString DiskFileA = "../ProdCxx/Prod2/histo.root_VPZon_8G"; //

/// obsolete:(
//TString DiskFileG = "../ProdCxx/Prod2/mcgen.root_92GeV";   //

// Dump file for temporarary objects
TFile DiskFileB("Plot1.root","RECREATE","histograms");
// Latex source file
FILE *DiskFileTeX;

ofstream   OutFile("Bhl6.output",ios::out);  // Logfile output
//=============================================================================
///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
//
float  gXcanv = 50, gYcanv = 50;
//Double_t sqr( const Double_t x ){ return x*x;};
TMCgenBHL6 *gBHLgen;
TRobolProd *gRoboT;

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
//////////////////////
double CMSENE= gBHLgen->m_CMSENE;
cout<< "@@@gBHLgen: CMSENE= "<<CMSENE<<endl;
double th1n  = gRoboT->m_th1n_calo2;     //  CALO2
double th2n  = gRoboT->m_th2n_calo2;     //  CALO2
double th1w  = gRoboT->m_th1w_calo2;     //  CALO2
double th2w  = gRoboT->m_th2w_calo2;     //  CALO2
cout<< "++++HistRemake CALO2: th1n= "<<th1n<<" th2n ="<< th2n<<"  th1w="<<th1w<<"  th2w="<<th2w<<endl;
// Born xsections for CALO2 narrwo and wide angular range
double xBornWW = gBHLgen->BornBhl(CMSENE,th1w,th2w);  //  CALO2
double xBornNN = gBHLgen->BornBhl(CMSENE,th1n,th2n);  //  CALO2
cout<< "@@@gBHLgen: xBornNN= "<<xBornNN<<endl;
cout<< "@@@gBHLgen: xBornWW= "<<xBornWW<<endl;
////////////////////////////////////


TH1D *HST_BHL_NORMA = (TH1D*)HstFile->Get("HST_BHL_NORMA");
///////////////////////////////////////////////////////////////////////////
// renormalizing 1-dim histos
//
  HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_v_bare1_ww") );
  HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_v_bare1_nw") );
  HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_v_bare1_nn") );
//
  HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_v_calo2_ww") );
  HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_v_calo2_nw") );
  HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_v_calo2_nn") );
//
  HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_v_sical2_ww") );
  HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_v_sical2_nw") );
  HisNorm1(HST_BHL_NORMA, (TH1D*)HstFile->Get("hst_v_sical2_nn") );
  // cumulative
  TH1D *hst_v_bare1_ww      = (TH1D*)HstFile->Get("hst_v_bare1_ww");
  TH1D *hst_v_bare1_nw      = (TH1D*)HstFile->Get("hst_v_bare1_nw");
  TH1D *hst_v_bare1_nn      = (TH1D*)HstFile->Get("hst_v_bare1_nn");
  //
  TH1D *hst_v_calo2_ww      = (TH1D*)HstFile->Get("hst_v_calo2_ww");
  TH1D *hst_v_calo2_nw      = (TH1D*)HstFile->Get("hst_v_calo2_nw");
  TH1D *hst_v_calo2_nn      = (TH1D*)HstFile->Get("hst_v_calo2_nn");
  //
  TH1D *hst_v_sical2_ww     = (TH1D*)HstFile->Get("hst_v_sical2_ww");
  TH1D *hst_v_sical2_nw     = (TH1D*)HstFile->Get("hst_v_sical2_nw");
  TH1D *hst_v_sical2_nn     = (TH1D*)HstFile->Get("hst_v_sical2_nn");
  //
  TH1D *HST_v_bare1_ww      = HstCumul("HST_v_bare1_ww", hst_v_bare1_ww);
  TH1D *HST_v_bare1_nw      = HstCumul("HST_v_bare1_nw", hst_v_bare1_nw);
  TH1D *HST_v_bare1_nn      = HstCumul("HST_v_bare1_nn", hst_v_bare1_nn);
  //
  TH1D *HST_v_calo2_ww      = HstCumul("HST_v_calo2_ww", hst_v_calo2_ww);
  TH1D *HST_v_calo2_nw      = HstCumul("HST_v_calo2_nw", hst_v_calo2_nw);
  TH1D *HST_v_calo2_nn      = HstCumul("HST_v_calo2_nn", hst_v_calo2_nn);
  //
  TH1D *HST_v_sical2_ww     = HstCumul("HST_v_sical2_ww", hst_v_sical2_ww);
  TH1D *HST_v_sical2_nw     = HstCumul("HST_v_sical2_nw", hst_v_sical2_nw);
  TH1D *HST_v_sical2_nn     = HstCumul("HST_v_sical2_nn", hst_v_sical2_nn);
///////////////////////////////////
// relative correction
  TH1D *hst_BornN = HstConst("hst_BornN",hst_v_calo2_ww, xBornNN);
  TH1D *hst_BornW = HstConst("hst_BornN",hst_v_calo2_ww, xBornWW);
////////////////////////////////////
  TH1D *HST_v_calo2_ww_del = HstDiff("HST_v_calo2_ww_del",HST_v_calo2_ww,hst_BornW, kRed);
  TH1D *HST_v_calo2_nw_del = HstDiff("HST_v_calo2_nw_del",HST_v_calo2_nw,hst_BornN, kBlue);
  TH1D *HST_v_calo2_nn_del = HstDiff("HST_v_calo2_nn_del",HST_v_calo2_nn,hst_BornN, kMagenta);
//
  HST_v_calo2_ww_del->Divide(hst_BornW);
  HST_v_calo2_nw_del->Divide(hst_BornN);
  HST_v_calo2_nn_del->Divide(hst_BornN);
// relative correction
  TH1D *HST_v_sical2_ww_del = HstDiff("HST_v_sical2_ww_del",HST_v_sical2_ww,hst_BornW, kRed);
  TH1D *HST_v_sical2_nw_del = HstDiff("HST_v_sical2_nw_del",HST_v_sical2_nw,hst_BornN, kBlue);
  TH1D *HST_v_sical2_nn_del = HstDiff("HST_v_sical2_nn_del",HST_v_sical2_nn,hst_BornN, kMagenta);
//
  HST_v_sical2_ww_del->Divide(hst_BornW);
  HST_v_sical2_nw_del->Divide(hst_BornN);
  HST_v_sical2_nn_del->Divide(hst_BornN);


  cout<<"----------------------------- HistRemake end ----------------------------------"<<endl;
}//HistRemake



///////////////////////////////////////////////////////////////////////////////////
void FigWT(TFile *HstFile)
{
//------------------------------------------------------------------------
  cout<<" ========================= FigWT =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;

  TH1D *hst_weight       = (TH1D*)HstFile->Get("hst_weight");
  TH1D *hst_weight2      = (TH1D*)HstFile->Get("hst_weight2");
  TH1D *hst_wt2trig      = (TH1D*)HstFile->Get("hst_wt2trig");

  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigWT = new TCanvas("cFigWT","cFigWT", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cFigWT->SetFillColor(10);
  cFigWT->Divide( 2,  0);
  //////////////////////////////////////////////
  cFigWT->cd(1);
  gPad->SetLogy(); // !!!!!!
  hst_weight2->DrawCopy("h");

//  CaptT->DrawLatex(0.12,0.95,"A_{FB}(v_{max}), ????");
  //////////////////////////////////////////////
  cFigWT->cd(2);
  gPad->SetLogy(); // !!!!!!
  hst_wt2trig->DrawCopy("h");
  //
  cFigWT->cd();
//
}// FigWT


///////////////////////////////////////////////////////////////////////////////////
void FigCalo2a(TFile *HstFile)
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCalo2a =========================== "<<endl;
  TH1D *hst_v_calo2_ww      = (TH1D*)HstFile->Get("hst_v_calo2_ww");
  TH1D *hst_v_calo2_nw      = (TH1D*)HstFile->Get("hst_v_calo2_nw");
  TH1D *hst_v_calo2_nn      = (TH1D*)HstFile->Get("hst_v_calo2_nn");

  TH1D *HST_v_calo2_ww      = (TH1D*)DiskFileB.Get("HST_v_calo2_ww");
  TH1D *HST_v_calo2_nw      = (TH1D*)DiskFileB.Get("HST_v_calo2_nw");
  TH1D *HST_v_calo2_nn      = (TH1D*)DiskFileB.Get("HST_v_calo2_nn");
  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
 ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCalo2a = new TCanvas("cFigCalo2a","cFigCalo2a", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cFigCalo2a->SetFillColor(10);
  cFigCalo2a->Divide( 2,  0);
  //////////////////////////////////////////////
  cFigCalo2a->cd(1);
  gPad->SetLogy(); // !!!!!!
  TH1D *HST1 = hst_v_calo2_ww;
//
  HST1->SetStats(0);
  HST1->SetTitle(0);
//  HST1->SetMaximum(  0.0); HST1->SetMinimum( 40.0);
  HST1->GetXaxis()->SetTitle("v=1-z");

  HST1->DrawCopy("h");
  CaptT->DrawLatex(0.06,0.95, "d#sigma/dv [nb]");
  double ycapt = 0.90; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,"CALO2"); ycapt += -0.02;

  PlotSame2(hst_v_calo2_nw, ycapt, kBlue,     0.25, "(a)", "Narrow-Wide ");
  PlotSame2(hst_v_calo2_nn, ycapt, kMagenta,  0.45, "(b)", "Narrow-Narrow ");
  PlotSame2(hst_v_calo2_ww, ycapt, kRed,      0.35, "(c)", "Wide-Wide ");

//////////////////////////////////////////////
  cFigCalo2a->cd(2);
  TH1D *HST2 = HST_v_calo2_ww;
//
  HST2->SetStats(0);
  HST2->SetTitle(0);
//  HST2->SetMaximum(  0.0);
  HST2->SetMinimum( 40.0);
  HST2->GetXaxis()->SetTitle("v_{max}");

  HST2->DrawCopy("h");
  CaptT->DrawLatex(0.06,0.95, "#sigma(v_{max}) [nb]");
  ycapt = 0.40; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,"CALO2 "); ycapt += -0.02;

  PlotSame2(HST_v_calo2_nw, ycapt, kBlue,     0.25, "(a)", "Narrow-Wide ");
  PlotSame2(HST_v_calo2_nn, ycapt, kMagenta,  0.40, "(b)", "Narrow-Narrow ");
  PlotSame2(HST_v_calo2_ww, ycapt, kRed,      0.55, "(c)", "Wide-Wide ");

//
}// FigCalo2a


///////////////////////////////////////////////////////////////////////////////////
void FigCalo2b(TFile *HstFile)
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCalo2b =========================== "<<endl;
  TH1D *HST_v_calo2_ww_del  = (TH1D*)DiskFileB.Get("HST_v_calo2_ww_del");
  TH1D *HST_v_calo2_nw_del  = (TH1D*)DiskFileB.Get("HST_v_calo2_nw_del");
  TH1D *HST_v_calo2_nn_del  = (TH1D*)DiskFileB.Get("HST_v_calo2_nn_del");

  TH1D *HST_v_sical2_ww_del  = (TH1D*)DiskFileB.Get("HST_v_sical2_ww_del");
  TH1D *HST_v_sical2_nw_del  = (TH1D*)DiskFileB.Get("HST_v_sical2_nw_del");
  TH1D *HST_v_sical2_nn_del  = (TH1D*)DiskFileB.Get("HST_v_sical2_nn_del");

  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCalo2b = new TCanvas("cFigCalo2b","cFigCalo2b", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cFigCalo2b->SetFillColor(10);
  cFigCalo2b->Divide( 2,  0);
  //////////////////////////////////////////////
  cFigCalo2b->cd(1);
  //gPad->SetLogy(); // !!!!!!
  TH1D *HST1 = HST_v_calo2_nw_del;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->SetMaximum(  0.0);
  HST1->SetMinimum(-0.15);
  HST1->GetXaxis()->SetTitle("v_{max}");

  HST1->DrawCopy("h");
  CaptT->DrawLatex(0.06,0.95, "#delta(v_{max}) ");
  double ycapt = 0.99; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,"CALO2 [350]"); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,"hep-ph/9602393, Fig. 13"); ycapt += -0.02;

  PlotSame2(HST_v_calo2_nw_del, ycapt, kBlue,     0.015, "(a)", "Narrow-Wide ");
  PlotSame2(HST_v_calo2_nn_del, ycapt, kMagenta,  0.015, "(b)", "Narrow-Narrow ");
  PlotSame2(HST_v_calo2_ww_del, ycapt, kRed,      0.015, "(c)", "Wide-Wide ");

//////////////////////////////////////////////
  cFigCalo2b->cd(2);
  TH1D *HST2 = HST_v_sical2_nw_del;

  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->SetMaximum(  0.0);
  HST2->SetMinimum(-0.15);
  HST2->GetXaxis()->SetTitle("v_{max}");

  HST2->DrawCopy("h");
  CaptT->DrawLatex(0.06,0.90, "#delta(v_{max}) ");
  ycapt = 0.99; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,"SICAL2 [360]"); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,"hep-ph/9602393, Fig. 14"); ycapt += -0.02;

  PlotSame2(HST_v_sical2_nw_del, ycapt, kBlue,     0.015, "(a)", "Narrow-Wide ");
  PlotSame2(HST_v_sical2_nn_del, ycapt, kMagenta,  0.015, "(b)", "Narrow-Narrow ");
  PlotSame2(HST_v_sical2_ww_del, ycapt, kRed,      0.015, "(c)", "Wide-Wide ");

//================================================
  cFigCalo2b->SaveAs("cFigCalo2b.pdf");
//
}// FigCalo2b


///////////////////////////////////////////////////////////////////////////////////
void TabOldBench()
{
//------------------------------------------------------------------------
  cout<<" ========================= TabOldBench start=========================== "<<endl;
//
  TH1D *HST_v_bare1_ww  = (TH1D*)DiskFileB.Get("HST_v_bare1_ww");
  TH1D *HST_v_bare1_nw  = (TH1D*)DiskFileB.Get("HST_v_bare1_nw");
  TH1D *HST_v_bare1_nn  = (TH1D*)DiskFileB.Get("HST_v_bare1_nn");

  TH1D *HST_v_calo2_ww  = (TH1D*)DiskFileB.Get("HST_v_calo2_ww");
  TH1D *HST_v_calo2_nw  = (TH1D*)DiskFileB.Get("HST_v_calo2_nw");
  TH1D *HST_v_calo2_nn  = (TH1D*)DiskFileB.Get("HST_v_calo2_nn");

  TH1D *HST_v_sical2_ww  = (TH1D*)DiskFileB.Get("HST_v_sical2_ww");
  TH1D *HST_v_sical2_nw  = (TH1D*)DiskFileB.Get("HST_v_sical2_nw");
  TH1D *HST_v_sical2_nn  = (TH1D*)DiskFileB.Get("HST_v_sical2_nn");


// Column captions
  int nPlt=3;   // No of histos/column +1
  Char_t *Capt[nPlt+1];
  for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
  strcpy(Capt[0],"{\\color{blue}$1-z_{\\min}$}");
  strcpy(Capt[1],"{\\color{blue} WW }");
  strcpy(Capt[2],"{\\color{blue} NW }");
  strcpy(Capt[3],"{\\color{red}  NN }");

// formats, not used in PlTable2
//  Char_t fmt[3][10];
//  strcpy(fmt[0],"f10.2");
//  strcpy(fmt[1],"f10.4");
//  strcpy(fmt[2],"f8.4");

// pointers to histograms
  TH1D *iHst[nPlt+1];

  // multicolumn caption
  Char_t Mcapt[132];

///************************************
  DiskFileTeX = fopen("TabOldBench.txp","w");
//************************************
// Initialization of the latex source file
  LibPLT.PlInit(DiskFileTeX, 2);
//  int k1,k2,dk; // index loop

  strcpy(Mcapt,"{\\color{red} BARE1 $\\sigma(z_{\\min})$ [nb]}");
  iHst[1]= HST_v_bare1_ww;  // bare1 WW
  iHst[2]= HST_v_bare1_nw;  // bare1 NW
  iHst[3]= HST_v_bare1_nn;  // bare1 NN
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "B",10,90, 20); // for 100 bins

  strcpy(Mcapt,"{\\color{red} CALO2 $\\sigma(z_{\\min})$ [nb]}");
  iHst[1]= HST_v_calo2_ww;  // CALO2 WW
  iHst[2]= HST_v_calo2_nw;  // CALO2 NW
  iHst[3]= HST_v_calo2_nn;  // CALO2 NN
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",10,90, 20); // for 100 bins

  strcpy(Mcapt,"{\\color{red} SICAL2 $\\sigma(z_{\\min})$ [nb]}");
  iHst[1]= HST_v_sical2_ww;  // CALO2 WW
  iHst[2]= HST_v_sical2_nw;  // CALO2 NW
  iHst[3]= HST_v_sical2_nn;  // CALO2 NN
  LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "E",10,90, 20); // for 100 bins

// finalizing latex source file
  LibPLT.PlEnd(DiskFileTeX);
//************************************
  fclose(DiskFileTeX);
//************************************
  cout<<" ========================= TabOldBench end =========================== "<<endl;
}//TabOldBench


///////////////////////////////////////////////////////////////////////////////////
void TabVP1()
{
//------------------------------------------------------------------------
  cout<<" ========================= TabVP1 start=========================== "<<endl;

    //////////////////////
    double CMSENE  = gBHLgen->m_CMSENE;
    double   th1n  = gRoboT->m_th1n_calo2;     //  CALO2
    double   th2n  = gRoboT->m_th2n_calo2;     //  CALO2
    double   th1w  = gRoboT->m_th1w_calo2;     //  CALO2
    double   th2w  = gRoboT->m_th2w_calo2;     //  CALO2

    double TRMIN = -sqr(CMSENE)*(1-cos(th1n))/2;
    double TRMAX = -sqr(CMSENE)*(1-cos(th2n))/2;
    cout<< " TRMAX= "<<TRMAX<<"  TRMIN= "<< TRMIN<<"  GeV^2"<<endl;
    double TMEDI = 2/(1/TRMIN+1/TRMAX);
    cout<< " TMEDI= " << TMEDI <<"  GeV^2"<<endl;


    double tran[4]; tran[1]=TRMIN; tran[2]=TMEDI; tran[3]=TRMAX;
//    tran[1]= -1.87227;
    TH1D* HST_tra = new TH1D("HST_tra" , "transfers", 3, 0.0, 3.0); HST_tra->Sumw2();
    for( int i=1; i<=3; i++ ) {
          HST_tra->SetBinContent(i,tran[i]); HST_tra->SetBinError(i,0.0);}

    double vp1,vp2,vp3,RePiE,dRePiE; int KeyPia;
    double SINW2 = 0.2319e0; // for VP (Jegerlehner), hardwired in BHLUMI4
    //////////////////////////////////////////////////////////
    // KeyPia =0 OFF, it used in semianalytical tests,
    //        =1 Burkhardt et.al. 1989, as in BHLUMI 2.0x
    //        =2 S. Eidelman, F. Jegerlehner, Z. Phys. C (1995)
    //        =3 Burkhardt and Pietrzyk 1995 (Moriond).
    //        =4 F.Jegerlehner 2017
    //SUBROUTINE  vacpol(KeyPia,Q2,SINW2,RePiE,dRePiE)
    TH1D* HST_VP1 = new TH1D("HST_VP1" , "Vac.Pol.", 3, 0.0, 3.0); HST_VP1->Sumw2();
    TH1D* HST_VP2 = new TH1D("HST_VP2" , "Vac.Pol.", 3, 0.0, 3.0); HST_VP2->Sumw2();
    TH1D* HST_VP3 = new TH1D("HST_VP3" , "Vac.Pol.", 3, 0.0, 3.0); HST_VP3->Sumw2();
    TH1D* HST_VP4 = new TH1D("HST_VP4" , "Vac.Pol.", 3, 0.0, 3.0); HST_VP4->Sumw2();
    for( int i=1; i<=3; i++ ) {
   	     KeyPia=1; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_VP1->SetBinContent(i,RePiE); HST_VP1->SetBinError(i,dRePiE);
    	 KeyPia=2; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_VP2->SetBinContent(i,RePiE); HST_VP2->SetBinError(i,dRePiE);
    	 KeyPia=3; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_VP3->SetBinContent(i,RePiE); HST_VP3->SetBinError(i,dRePiE);
    	 KeyPia=4; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_VP4->SetBinContent(i,RePiE); HST_VP4->SetBinError(i,dRePiE);
         cout<<" tran="<<tran[i]<<"  !!!!! RePiE,dRePiE ="<<RePiE<<"   "<<dRePiE<<endl;
   }
    HST_VP1->Scale(100);HST_VP2->Scale(100);HST_VP3->Scale(100);HST_VP4->Scale(100);

//////////////////////////////////////////////////////////////////////
    char *fmt0 = "$  %10.0f $";
    char *fmt1 = "& $ %10.4f \\pm %8.4f $ ";
    fmt1 = "& $ %10.4f $ "; // No error in 1st column
    LibPLT.Setfmt0(fmt0);
    LibPLT.Setfmt1(fmt1);
// Column captions
    int nPlt=5;   // No of histos/column +1
    Char_t *Capt[nPlt+1];
    for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
    strcpy(Capt[0],"  No. ");
    strcpy(Capt[1],"{\\color{blue} $t~[GeV^2]$ }");
    strcpy(Capt[2],"{\\color{blue} Helm89 }");
    strcpy(Capt[3],"{\\color{red}  Fred95 }");
    strcpy(Capt[4],"{\\color{blue} Bolek95 }");
    strcpy(Capt[5],"{\\color{red}  Fred17 }");
//////////////////////////////////////////////////////////////////////
    TH1D *iHst[nPlt+1]; // pointers to histograms
    Char_t Mcapt[132];  // multicolumn caption
///////////////  LaTeX source file
    DiskFileTeX = fopen("TabVP1.txp","w");
    LibPLT.PlInit(DiskFileTeX, 2); // Initialization of the latex source file

    strcpy(Mcapt,"{\\color{red} Vacuum polarization $\\times$ 100 }");
    iHst[1]= HST_tra;
    iHst[2]= HST_VP1;
    iHst[3]= HST_VP2;
    iHst[4]= HST_VP3;
    iHst[5]= HST_VP4;
    LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, " ",1,3, 1); // for 100 bins

// finalizing latex source file
    LibPLT.PlEnd(DiskFileTeX);
    fclose(DiskFileTeX);
//////////////////////////////////////////////////////////////////////
 cout<<" ========================= TabVP1 end =========================== "<<endl;
}//TabVP1



///////////////////////////////////////////////////////////////////////////////////
void TabVP2()
{
//------------------------------------------------------------------------
 cout<<" ======================== TabVP2 start=========================== "<<endl;

 double CMSENE = gBHLgen->m_CMSENE;
///////////////////////////////////////////////////////
// Event selection params from the Robol object
 double TminF = gRoboT->m_TminF; // fidutial
 double TmaxF = gRoboT->m_TmaxF; // fidutial
 cout<< "!!!gRoboT: TminF= "<<TminF<<"  TmaxF= "<< TmaxF<<endl;
 ////////////////////////
 double   th1n  = gRoboT->m_th1n_calo2;     //  CALO2,SICAL2
 double   th2n  = gRoboT->m_th2n_calo2;     //  CALO2,SICAL2
 double   th1w  = gRoboT->m_th1w_calo2;     //  CALO2,SICAL2
 double   th2w  = gRoboT->m_th2w_calo2;     //  CALO2,SICAL2
 cout<< "====HistRemake CALO2: th1n= "<<th1n<<" th2n ="<< th2n<<"  th1w="<<th1w<<"  th2w="<<th2w<<endl;
 double TRMED,SINW2,RePiE,dRePiE;
 double TRMIN = -sqr(CMSENE)*(1-cos(th1n))/2;
 double TRMAX = -sqr(CMSENE)*(1-cos(th2n))/2;
 double TMEDI = 2/(1/TRMIN+1/TRMAX);

// Born xsections for CALO2 narrwow and wide angular range
// Pure t-chanel, no VP no Z !!!
 double xBornWW = gBHLgen->BornBhl(CMSENE,th1w,th2w);  // CALO2,SICAL2
 double xBornNN = gBHLgen->BornBhl(CMSENE,th1n,th2n);  // CALO2,SICAL2
 cout<< "@@@gBHLgen: xBornNN= "<<xBornNN<<endl;
 cout<< "@@@gBHLgen: xBornWW= "<<xBornWW<<endl;
//////////////////////////////////////////////////////////
// KeyPia =0 OFF, it used in semianalytical tests,
//        =1 H. Burkhardt et.al. 1989, as in BHLUMI 2.0x
//        =2 S. Eidelman, F. Jegerlehner, Z. Phys. C (1995)
//        =3 H. Burkhardt and B. Pietrzyk 1995 (Moriond).
//        =4 F.Jegerlehner 2017
 int KeyPia=0,KeyPib=0;
 double xBorn,xBorn2,dBorn,rBorn,xBorn0,dBorn0;
 TH1D* HST_sig1 = new TH1D("HST_sig1" , "Vac.Pol.", 5, 0.0, 5.0); HST_sig1->Sumw2();
 TH1D* HST_sig2 = new TH1D("HST_sig2" , "Vac.Pol.", 5, 0.0, 5.0); HST_sig2->Sumw2();
 for( int i=0; i<=4; i++ ) {
   KeyPia=i; gBHLgen->SetKeyPia(KeyPia);    // Type of VP
   KeyPib=0; gBHLgen->SetKeyPib(KeyPib);    // xBorn with VP not shifted by its error
   xBorn  =  gBHLgen->BornBhl2(th1n,th2n);  // N-N as in CALO2/SICAL2
   KeyPib=1; gBHLgen->SetKeyPib(KeyPib);    // xBorn with VP shifted by its error
   xBorn2 =  gBHLgen->BornBhl2(th1n,th2n);  // N-N as in CALO2/SICAL2
   KeyPib=0; gBHLgen->SetKeyPib(KeyPib);    // xBorn with VP not shifted by its error
   dBorn = xBorn2-xBorn;
   rBorn = dBorn/xBorn *100;   // relative error due to VP times 100
   HST_sig1->SetBinContent(i+1,xBorn); HST_sig1->SetBinError(i+1,dBorn);
   cout<<" KeyPia="<<KeyPia<<"  xBorn="<<xBorn<<"  xBorn2="<<xBorn2<<"  dBorn="<<dBorn<<"  rBorn="<<rBorn<<endl;
   vacpol_(KeyPia,TMEDI,SINW2,RePiE,dRePiE);
   xBorn0 = xBornNN; dBorn0 = 0;
   if( i>0) xBorn0 = xBornNN/(1+2*RePiE);
   if( i>0) dBorn0 = xBornNN*( 1/(1+2*RePiE)- 1/(1+2*RePiE+2*dRePiE));
   HST_sig2->SetBinContent(i+1,xBorn0); HST_sig2->SetBinError(i+1,dBorn0);
}

 //////////////////////////////////////////////////////////////////////
   char *fmt0 = "$  %10.0f $";
   char *fmt1 = "& $ %10.4f \\pm %8.4f $ ";
   LibPLT.Setfmt0(fmt0);
   LibPLT.Setfmt12(fmt1);
 // Column captions
   int nPlt=2;   // No of histos/column +1
   Char_t *Capt[nPlt+1];
   for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
   strcpy(Capt[0],"  No. ");
   strcpy(Capt[1],"{\\color{blue} $ \\int_{t_1}^{t_2} \\frac{dt}{1-2\\Pi(t)} \\; \\frac{d\\sigma_0}{dt} $ }");
   strcpy(Capt[2],"{\\color{blue} $\\frac{1}{1-2\\Pi(t_{med})} \\int_{t_1}^{t_2} d\\sigma_0$ } ");
 //////////////////////////////////////////////////////////////////////
   TH1D *iHst[nPlt+1]; // pointers to histograms
   Char_t Mcapt[132];  // multicolumn caption
 ///////////////  LaTeX source file
   DiskFileTeX = fopen("TabVP2.txp","w");
   LibPLT.PlInit(DiskFileTeX, 2); // Initialization of the latex source file

   iHst[1]= HST_sig1;
   iHst[2]= HST_sig2;
   strcpy(Mcapt,"{\\color{red} Born x-section [nb], no VP, no $Z$, no $\\gamma_s$)}");
   LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "B",1,1, 1); // for 100 bins
   strcpy(Mcapt,"{\\color{red} Born x-section [nb], with VP, no $Z$, no $\\gamma_s$}");
   LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "E",2,5, 1); // for 100 bins

 // finalizing latex source file
   LibPLT.PlEnd(DiskFileTeX);
   fclose(DiskFileTeX);

//************************************
 cout<<" ======================== TabVP2 end ============================ "<<endl;
}//TabVP2



///////////////////////////////////////////////////////////////////////////////////
void TabVP3()
{
//------------------------------------------------------------------------
  cout<<" ========================= TabVP3 start=========================== "<<endl;

    //////////////////////
    double CMSENE  = gBHLgen->m_CMSENE;
    CMSENE  = 91.224;
    double tminF = 0.020, tmaxF=0.130-1e-7;

    double tran[13],thet[13];
    double RePiE,dRePiE; int KeyPia;
    double SINW2 = 0.2319e0; // for VP (Jegerlehner), hardwired in BHLUMI4
    //////////////////////////////////////////////////////////
    // KeyPia =0 OFF, it used in semianalytical tests,
    //        =1 Burkhardt et.al. 1989, as in BHLUMI 2.0x
    //        =2 S. Eidelman, F. Jegerlehner, Z. Phys. C (1995)
    //        =3 Burkhardt and Pietrzyk 1995 (Moriond).
    //        =4 F.Jegerlehner 2017
    //SUBROUTINE  vacpol(KeyPia,Q2,SINW2,RePiE,dRePiE)
    TH1D* HST_THE = new TH1D("HST_THE" , "theta   ", 12, 0, 12); HST_THE->Sumw2();
    TH1D* HST_TRA = new TH1D("HST_TRA" , "transfer", 12, 0, 12); HST_TRA->Sumw2();
    TH1D* HST_vp1 = new TH1D("HST_vp1" , "Vac.Pol.", 12, 0, 12); HST_vp1->Sumw2();
    TH1D* HST_vp2 = new TH1D("HST_vp2" , "Vac.Pol.", 12, 0, 12); HST_vp2->Sumw2();
    TH1D* HST_vp3 = new TH1D("HST_vp3" , "Vac.Pol.", 12, 0, 12); HST_vp3->Sumw2();
    TH1D* HST_vp4 = new TH1D("HST_vp4" , "Vac.Pol.", 12, 0, 12); HST_vp4->Sumw2();
    TH1D* HST_vp5 = new TH1D("HST_vp5" , "Vac.Pol.", 12, 0, 12); HST_vp5->Sumw2();
    TH1D* HST_vp6 = new TH1D("HST_vp6" , "Vac.Pol.", 12, 0, 12); HST_vp6->Sumw2();
    for( int i=1; i<=12; i++ ) {
    	 thet[i]= tminF +(i-1)*(tmaxF-tminF)/11e0;
         HST_THE->SetBinContent(i,thet[i]); HST_vp4->SetBinError(i,0);
         tran[i]= -sqr(CMSENE)*(1-cos(thet[i]))/2;
         HST_TRA->SetBinContent(i,tran[i]); HST_vp4->SetBinError(i,0);
  	     KeyPia=1; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_vp1->SetBinContent(i,RePiE); HST_vp1->SetBinError(i,dRePiE);
    	 KeyPia=2; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_vp2->SetBinContent(i,RePiE); HST_vp2->SetBinError(i,dRePiE);
    	 KeyPia=3; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_vp3->SetBinContent(i,RePiE); HST_vp3->SetBinError(i,dRePiE);
    	 KeyPia=4; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_vp4->SetBinContent(i,RePiE); HST_vp4->SetBinError(i,dRePiE);
        // New Teubner18
    	 KeyPia=5; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_vp5->SetBinContent(i,RePiE); HST_vp5->SetBinError(i,dRePiE);
         // Davier 2020
    	 KeyPia=6; vacpol_(KeyPia,tran[i],SINW2,RePiE,dRePiE);
         HST_vp6->SetBinContent(i,RePiE); HST_vp6->SetBinError(i,dRePiE);
         //
         cout<<"i="<<i<<" thet="<<thet[i]<<" tran="<<tran[i]<<" RePiE,dRePiE ="<<RePiE<<"   +-"<<dRePiE<<endl;
    }
    HST_vp1->Scale(100);HST_vp2->Scale(100);HST_vp3->Scale(100);HST_vp4->Scale(100);HST_vp5->Scale(100);HST_vp6->Scale(100);

//////////////////////////////////////////////////////////////////////
    char *fmt0 = "$  %10.2f $";               // 0-th column
    char *fmt1 = "& $ %10.4f \\pm %8.4f $ ";  // 1-st column
    char *fmt2 = "& $ %10.4f \\pm %8.4f $ ";  // next columns
    char *fmt3 = "& $ %10.4f $ ";             // error=0 case
    fmt0 = "  $ %10.0f $";             // 0-th column
    fmt1 = "& $ %10.3f \\pm %8.3f $ "; // 1-st column
    fmt2 = "& $ %10.3f \\pm %8.3f $ "; // next columns
    fmt3 = "& $ %10.4f $ ";            // error=0 case
    LibPLT.Setfmt0(fmt0);
    LibPLT.Setfmt1(fmt1);
    LibPLT.Setfmt2(fmt2);
    LibPLT.Setfmt3(fmt3);
// Column captions
    int nPlt=8;   // No of histos/column +1
    Char_t *Capt[nPlt+1];
    for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
    strcpy(Capt[0],"  No. ");
    strcpy(Capt[1],"  $\\theta$ (rad)");
    strcpy(Capt[2],"{\\color{blue} $t~[GeV^2]$ }");
    strcpy(Capt[3],"{\\color{blue} Burkhardt89 }");
    strcpy(Capt[4],"{\\color{red}  Jegerlehner95 }");
    strcpy(Capt[5],"{\\color{blue} Burkhardt95 }");
    strcpy(Capt[6],"{\\color{red}  Jegerlehner19 }");
    strcpy(Capt[7],"{\\color{red}  Teubner20 }");
    strcpy(Capt[8],"{\\color{red}  Davier20 }");
//////////////////////////////////////////////////////////////////////
    TH1D *iHst[nPlt+1]; // pointers to histograms
    Char_t Mcapt[132];  // multicolumn caption
///////////////  LaTeX source file
    DiskFileTeX = fopen("TabVP3.txp","w");
    LibPLT.PlInit(DiskFileTeX, 2); // Initialization of the latex source file

    strcpy(Mcapt,"{\\color{red} VP $\\Re \\Pi(t)\\times 100$:~~~ Leptonic VP of Teuner in Teub20, Fred19, Davier20. Top ON. NEW! }");
    iHst[1]= HST_THE;
    iHst[2]= HST_TRA;
    iHst[3]= HST_vp1;
    iHst[4]= HST_vp2;
    iHst[5]= HST_vp3;
    iHst[6]= HST_vp4;
    iHst[7]= HST_vp5;
    iHst[8]= HST_vp6;
    LibPLT.PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, " ",1,12, 1); // for 11 bins

// finalizing latex source file
    LibPLT.PlEnd(DiskFileTeX);
    fclose(DiskFileTeX);
//////////////////////////////////////////////////////////////////////
 cout<<" ========================= TabVP3 end =========================== "<<endl;
}//TabVP3





///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
// libProd.so should be loaded prior to opening GenFile
  gSystem->Load("../ProdCxx/.libs/libProd.so"); // seems to be needed !!!
//
  TFile *GenFile = new TFile(DiskFileG);  // File with Generator objects
  TFile *HstFile = new TFile(DiskFileA);  // File with histograms
  //
  LibPLT.Initialize();
  //
  DiskFileB.cd();

  TH1D *HST_BHL_NORMA = (TH1D*)HstFile->Get("HST_BHL_NORMA");
  gBHLgen= (TMCgenBHL6*)GenFile->Get("MCgen");
  gBHLgen->ReInitialize(&OutFile, HST_BHL_NORMA);
  gRoboT = (TRobolProd*)GenFile->Get("RoboT");


  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++

  HistRemake(HstFile);     // Renormalization of MC histograms

  //========== PLOTTING ==========
  // Template empty canvas  with 2 figures
  //FigTempl();

  FigWT(HstFile);
  FigCalo2a(HstFile);
  FigCalo2b(HstFile);
  //
  TabOldBench();
  TabVP1();
  TabVP2();
  TabVP3();

  //++++++++++++++++++++++++++++++++++++++++

  HstFile->ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  OutFile.close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


