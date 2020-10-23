//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   HSTplot                                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// HSTplot is auxiliary toolbox for histograming, latex tables from histograms etc
/////////////////////////////////////////////////////////////////////////////

#ifndef HSTplot_H
#define HSTplot_H
#include<stdlib.h>
#include<stdio.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"


class HSTplot{
// Interface and extensions to KKsem toolbox
 public:
	//
 	double m_gnanob;
 	double m_pi;
 	double m_ceuler;
 	//
	const char *fmt0;
	const char *fmt1;
	const char *fmt2;
	const char *fmt3;

//
//------ constructors destructors -------
 public:
  HSTplot(){;}
  ~HSTplot(){;}
  HSTplot(const char* Name);
public:
  // Interfaces to KKsem integration routines using Gauss method
  void Initialize();
  void PlInit(FILE *ltx, int lint);
  void PlTable2(int Ncol, TH1D *iHst[], FILE *ltex, Char_t *Capt[], Char_t Mcapt[] , const char *chr1, int k1,int k2,int dk);
  void PlEnd(FILE *ltex);
  void Setfmt0(char *fmt) { fmt0 = fmt; };
  void Setfmt1(char *fmt) { fmt1 = fmt; };
  void Setfmt2(char *fmt) { fmt2 = fmt; };
  void Setfmt3(char *fmt) { fmt3 = fmt; };
  void Setfmt12(char *fmt) { fmt1 = fmt;fmt2 = fmt; };
////////////////////////////////////////////////////////////////////////////
};// HSTplot

#endif
