#ifndef TRobolProd_H
#define TRobolProd_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS TRobolProd                                           //
//                                                                           //
//     Three-stroke engine for analysis: Initialize-Production-Finalize      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

/// OUR headers
#include "TMCgenBHL6.h"
#include "TRig.h"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"


class TRobolProd : public TObject {
///--- data members
/// Part of former Base class MCdev
 public:
  char     f_Name[64];          // Name of a give instance of the class
  TString  f_HistNormName;      // name of the normalization histo
  TRandom  *f_RNgen;            //! Central RN event generator
  TMCgenBHL6  *f_MCgen;         //! MC event generator
  TRig     *f_TRigLib;          //  LEP trigger collection
  TFile    *f_HstFile;          //! ROOT DiskFile with all histos
  TFile    *f_GenFile;          //! ROOT DiskFile with MC generator
  ofstream *f_Out;              //! Central Logfile for messages
  ofstream  f_TraceFile;        //! Special DiskFile for debug
  double    f_NevGen;           //  event serial number
  double    f_count1;           //  auxiliary event counter (debug)
  long      f_isNewRun;         //  "is this new run?"
//////////////////////////////////////////////////////////////
 public:
  long   m_NevGen;             // event serial number
  long   m_count1;             // auxiliary event counter (debug)
  double m_pi;                 ///
// ===============  params of event selection   =======================
  int      m_NEVTOT;   //   = 84000000, total no. of MC events to generate
  double   m_VMAXE;    //  VMAXE   v_max              trigger maximum v
  int      m_NPHI;     //  NPHI    no of phi   sectors, CALO2, SICAL2
  int      m_NTHE;     //  NTHE    no of theta sectors, CALO2, SICAL2
//
  double   m_TminF;    //  TminW   fidutial theta_min wide, CALO2, SICAL2
  double   m_TmaxF;    //  TmaxW   fidutial theta_max wide, CALO2, SICAL2
  double   m_dcone_calo1;    //  CALO1
  double   m_th1n_calo2;     //  CALO2
  double   m_th2n_calo2;     //  CALO2
  double   m_th1w_calo2;     //  CALO2
  double   m_th2w_calo2;     //  CALO2
  double   m_dlth_calo2;     //  CALO2
  double   m_dlph_calo2;     //  CALO2

  double   m_y1A[4];  //  y of BARE1
  double   m_y1B[4];  //  y of CALO1
  double   m_y1C[4];  //  y of CALO2

  double   m_Z1[4];   //  z of SICAL2
  double   m_Z2[4];   //  z of SICAL2
  double   m_Z3[4];   //  z of SICAL2


  // =============== local mirror of BHLUMI event =======================
  TLorentzVector m_pbea1;      //! initial beams
  TLorentzVector m_pbea2;      //! initial beams
  TLorentzVector m_pfer1;      //! final fermions
  TLorentzVector m_pfer2;      //! final fermions
  int            m_Nphot;      //! photon multiplicity
  TLorentzVector m_phot[100];  //! photon 4-momenta
  // ============== Histograms follow =================================
  //
  TH1D   *hst_weight;           //!  No streamer!!!
  TH1D   *hst_weight2;          //!  No streamer!!!
  TH1D   *hst_wt2trig;          //!  No streamer!!!
  TH1D   *hst_Costhe;           //!  No streamer!!!

  TH1D   *hst_v_bare1_ww;       //!
  TH1D   *hst_v_bare1_nw;       //!
  TH1D   *hst_v_bare1_nn;       //!

  TH1D   *hst_v_calo2_ww;       //!
  TH1D   *hst_v_calo2_nw;       //!
  TH1D   *hst_v_calo2_nn;       //!

  TH1D   *hst_v_sical2_ww;       //!
  TH1D   *hst_v_sical2_nw;       //!
  TH1D   *hst_v_sical2_nn;       //!

  // histos of of xsections as a function of the trigges type
  TH1D   *hst_trig_VP0;       //!
  TH1D   *hst_trig_VP1;       //!
  TH1D   *hst_trig_VP2;       //!
  TH1D   *hst_trig_VP3;       //!
  TH1D   *hst_trig_VP4;       //!
  TH1D   *hst_trig_VP5;       //!

  TH1D   *hst_trig_ZZ0;       //!
  TH1D   *hst_trig_ZZ1;       //!
  TH1D   *hst_trig_ZZ2;       //!
  TH1D   *hst_trig_ZZ3;       //!
  TH1D   *hst_trig_ZZ4;       //!
 //
///////////////////////////////////////////
/// mandatory constructors and destructors
  public:
  TRobolProd();                // explicit default constructor for streamer
  TRobolProd(const char*);     // user constructor
  virtual ~TRobolProd();       // explicit destructor
/// mandatory methods
  virtual void Initialize(ofstream*, TFile*, TFile*);
  virtual void Hbooker();
  virtual void Production(double &);
  void FillTrig92( TH1D *hst_trig , double WtMain );
  void FillTrig95( TH1D *hst_trig , double WtMain );
// from MCdev
  TH1D *TH1D_UP(const char*, const char*, int, double, double);
  TH2D *TH2D_UP(const char*, const char*, int, double, double, int, double, double);
  double sqr( const double x ){ return x*x;};
/// for debug
  void StopM(char* message){
    cout <<"++++ TRobol: "<< message << endl; exit(5);}    //Error message

//////////////////////////////////////////////
// Other user methods
  void Finalize();
  void MomPrint( TLorentzVector&);
////////////////////////////////////////////////////////////////////////////
                      ClassDef(TRobolProd,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
