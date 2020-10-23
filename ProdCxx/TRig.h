#ifndef TRig_H
#define TRig_H

///////////////////////////////////////////////////////////
//               Class   TRig                            //
//            wrappers for f77 triggers                  //
///////////////////////////////////////////////////////////

/// C++ headers
using namespace std;
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;


/// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TLorentzVector.h"

//////////////////////////////////////////////////////////////////////////////
//      Routine for CALO2 trigger of WORKSHOP 95 (Fig. 13)
//      DIMENSION y1A(4),y1B(4),y1C(4)
//      CALL TRIG_trical1(th1n,th2n,th1w,th2w,dcone,dlth,dlph,y1A,y1B,y1C)
//      SUBROUTINE trical1(th1n,th2n,th1w,th2w,delcon,dlthe,dlphi,zcalA,zcalB,zcalC)
extern "C"  void trig_trical1_( const double&,  const double&,  const double&,  const double&,
		                        const double&,  const double&,  const double&,
		                        double[], double[], double[]);
//////////////////////////////////////////////////////////////////////////////
//      Routine for SICAL2 trigger of WORKSHOP 95 (Fig. 14)
//      CALL trig_trisic2W(Nphi,Ntheta,Th1f,Th2f,Nseg,Npad,z3A,z3B,z3C)
//      SUBROUTINE TRIG_TRISIC2W (Nphi,Nthe,Tmind,Tmaxd,Nseg,Npad,Z1,Z2,Z3)
extern "C"  void trig_trisic2w_(const int&,  const int&,  const double&,  const double&,
		                        const int&,  const int&,
		                        double[], double[], double[]);

//////////////////////////////////////////////////////////////////
//*                   ALEPH LCAL semi-realistic
//      SUBROUTINE TRIG_TRIGAS2(TH1N,TH2N,TH1W,TH2W,NPHI,XWIDE,XNARR,XMIX1,XMIX2)
extern "C"  void trig_trigas2_( const double&, const double&, const double&, const double&,
                    const int&, const double&, const double&, const double&, const double&);

/////////////////////////////////////////////////////////////////////////////
//  New configurable version of trisic2w
//      SUBROUTINE TRIG_UniSical(lWW,lNN,lNW,lWN)
extern "C"  void trig_unisical_(const int&,  const int&, const int&,  const int& );

/////////////////////////////////////////////////////////////////////////////
//  New configurable version of
//  trigas1 of TH.6118, and trical2 of workshop96
//    SUBROUTINE TRIG_UniElCal(lWW,lNN,lNW,lWN)
extern "C"  void trig_unielcal_(const int&,  const int&, const int&,  const int& );

///////////////////////////////////////////////
//    Setters for configurable triggers      //
///////////////////////////////////////////////

//      SUBROUTINE TRIG_SetGrid(yTMIND,yTMAXD,yNPHI,yNTHE)
extern "C"  void trig_setgrid_(const double&,  const double&, const int&,  const int& );

//      SUBROUTINE trig_setclust(yNpad,yNseg)
extern "C"  void trig_setclust_(const int&,  const int& );

//      SUBROUTINE TRIG_SetAng(yth1w,yth2w, yth1n,yth2n)
extern "C"  void trig_setang_(const double&,  const double&,  const double&,  const double&);

//      SUBROUTINE TRIG_SetAcol(yAcopl,yAcoli)
extern "C"  void trig_setacol_(const double&,  const double& );

//      SUBROUTINE TRIG_SetEcut(yEcutMin,yEcutMax,yEcutSum,yEcutProd)
extern "C"  void trig_setecut_(const double&,  const double& ,  const double& ,  const double& );

//////////////////////////////////////////////////////////////////////////
//________________________________________________________________________
class TRig: public TObject {
public:
 char  f_Name[64];       // Name of a give instance of the class
 float f_Version;        // Actual VERSION of the program
 char  f_Date[40];       // Release DATE of the program
 //-------------- USER Input parameters--------------------------
public:
 double      m_CMSENE;  // not used
//
public:
TRig();                // explicit default constructor for streamer
TRig(const char*);     // user constructor
~TRig();               // explicit destructor
public:
/////////////////////////////////////////////////////////////////////////////
/// methods
void Initialize();

void PrintAng( const char *name, double th1w, double th2w,  double th1n, double th2n);

///////////////////////////////////////////////////////
//  Old semi-realistic triggers, without setters!    //
///////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void TriCal1( double &th1n,   double &th2n,  double &th1w,  double &th2w,
 		      double &dcone,  double &dlth,  double &dlph,
 		      double y1A[],   double y1B[],  double y1C[]) {
//  SUBROUTINE TRIG_trical1(th1n,th2n,th1w,th2w,delcon,dlthe,dlphi,zcalA,zcalB,zcalC)
    trig_trical1_(th1n,th2n,th1w,th2w,dcone,dlth,dlph,y1A,y1B,y1C);
}

void TriSic2w( int &Nphi,  int &Nthe, double &Tmind, double &Tmaxd,
		       double Z1[],   double Z2[],  double Z3[]) {
    int Nseg=1, Npad=1;           // as in demo2.f (workshop96, Fig.16)
//  SUBROUTINE TRIG_TRISIC2W (Nphi,Nthe,Tmind,Tmaxd,Nseg,Npad,Z1,Z2,Z3)
	trig_trisic2w_(Nphi, Nthe, Tmind, Tmaxd, Nseg, Npad, Z1, Z2, Z3);
}

///////////////////////////////////////////////
//   New configurable triggers with setters  //
///////////////////////////////////////////////

void UniSical(int &lWW,  int &lNN, int &lNW,  int &lWN )
//    SUBROUTINE TRIG_UniSical(lWW,lNN,lNW,lWN)
//    Configurable version of sical2, trisic2, trisic2w
  {  trig_unisical_(lWW,lNN,lNW,lWN); }


void UniElcal(int &lWW,  int &lNN, int &lNW,  int &lWN )
//   SUBROUTINE TRIG_UniElCal(lWW,lNN,lNW,lWN)
//   Formely trigas1 of TH.6118, and trical2 of workshop96
{  trig_unielcal_(lWW,lNN,lNW,lWN); }


//////////////////////////////////////////
//   Setters for configurable triggers  //
//////////////////////////////////////////

void SetTrisic2w()
// Setter for UniSical (trig_unisical_)
// workshop96 benchmark: simplified ALEPH 1993 Sical (trisic2w)
{
//    SUBROUTINE TRIG_SetGrid(yTMIND,yTMAXD,yNPHI,yNTHE)
  double TMIND = 0.024, TMAXD= 0.058;
  int    NPHI  = 32, NTHE = 16;
  trig_setgrid_( TMIND, TMAXD, NPHI, NTHE );
//  SUBROUTINE trig_setclust(yNpad,yNseg)
  int Npad= 1, Nseg = 1; // as in demo2.f (workshop96, Fig.16)
  trig_setclust_( Npad,  Nseg );
//  SUBROUTINE TRIG_SetAng(yth1w,yth2w, yth1n,yth2n)
  double th1w,th2w, th1n, th2n;
  double PAD = (TMAXD-TMIND)/NTHE;
  th1w = TMIND + PAD;
  th2w = TMAXD - PAD;
  th1n = TMIND + 2*PAD;
  th2n = TMAXD - 4*PAD;
  trig_setang_(th1w,th2w, th1n,th2n);
//  SUBROUTINE TRIG_SetAcol(yAcopl,yAcoli)
  double Acopl = 1e9, Acoli =1e9; // no cut-offs
  trig_setacol_( Acopl, Acoli);
//  SUBROUTINE TRIG_SetEcut(yEcutMin,yEcutMax,yEcutSum,yEcutProd)
  double EcutMin =0, EcutMax=0, EcutSum=0, EcutProd=0;
  EcutProd=0.5;  // z>0.5
  trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
} // SetTrisic2w

 void SetUniElcal()
// Setter for UniElCal (trig_unielcal_)
// trigas1 of TH.6118 or trical2 of workshop96
{
  double TMIND = 0.043, TMAXD=0.125;  // the same as wide
  int    NPHI  = 6, NTHE = 1;
  trig_setgrid_( TMIND, TMAXD, NPHI, NTHE );
//
  double th1w,th2w, th1n, th2n;
  th1w = 0.043;  //  wide, ALEPH 1990-92
  th2w = 0.125;  //  wide
  th1n = 0.057;  //  narrow
  th2n = 0.107;  //  narrow
  trig_setang_(th1w,th2w, th1n,th2n);
//
  double EcutMin = 0.440;         // minimum energy
  double EcutMax = 0;             // none
  double EcutSum = 0;             // none
  double EcutProd= 0;             // no cut
  trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
} // SetUniElcal


//////////////////////////////////////////////////////////////////
//   Selectors for 2-nd generation lumi detectors (93->95)      //
//////////////////////////////////////////////////////////////////

void SetAleph93()
// Setter for UniSical (trig_unisical_)
//  ALEPH 1993 SICAL detector
{
//  SUBROUTINE TRIG_SetGrid(yTMIND,yTMAXD,yNPHI,yNTHE)
  double TMIND = 0.024, TMAXD= 0.058;
  int    NPHI  = 32, NTHE = 16;
  trig_setgrid_( TMIND, TMAXD, NPHI, NTHE );
//
//  SUBROUTINE trig_setclust(yNpad,yNseg)
  int Npad= 3, Nseg = 1; // as in TRISIC2() in silicon.f
  trig_setclust_( Npad,  Nseg );
//
// SUBROUTINE TRIG_SetAng(yth1w,yth2w, yth1n,yth2n)
  double th1w,th2w, th1n, th2n;
  double PAD = (TMAXD-TMIND)/NTHE;
  th1w = TMIND + PAD;
  th2w = TMAXD - PAD;
  th1n = TMIND + 3*PAD;
  th2n = TMAXD - 4*PAD;
  trig_setang_(th1w,th2w, th1n,th2n);
  ///[[[[
  // PrintAng "SetAleph93", th1w, th2w,  th1n, th2n);
//
//  SUBROUTINE TRIG_SetAcol(yAcopl,yAcoli)
  double Acopl = 0.52359878; // delta phi
  double Acoli =   1e9;      // delta theta, no cut-offs
  trig_setacol_( Acopl, Acoli);
//
//  SUBROUTINE TRIG_SetEcut(yEcutMin,yEcutMax,yEcutSum,yEcutProd)

  double EcutMin = 0.43865902;  // minimum energy
  double EcutMax = 0.43865902 ; // maximum energy
  double EcutSum = 0.60315615 ; // sum of energies
  double EcutProd= 0.000 ;  // no cut
  trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
} // SetAleph93

void TAleph93(int &lWW,  int &lNN, int &lNW,  int &lWN )
//  ALEPH 1994 Sical detector
{
  SetAleph93();
  //
  trig_unisical_(lWW,lNN,lNW,lWN);
}//TAleph93

void TAleph94(int &lWW,  int &lNN, int &lNW,  int &lWN )
//  ALEPH 1994 Sical detector
{
  SetAleph93();
//  SUBROUTINE TRIG_SetEcut(yEcutMin,yEcutMax,yEcutSum,yEcutProd)
  double EcutMin = 0.43865902;  // minimum energy
  double EcutMax = 0.43865902 ; // maximum energy
  double EcutSum = 0.780;        // sum of energies
  double EcutProd= 0.000;  // no cut
  trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
  //
  trig_unisical_(lWW,lNN,lNW,lWN);
} // TAleph94


void TAleph95(int &lWW,  int &lNN, int &lNW,  int &lWN )
//  ALEPH 1995 Sical detector
{
  SetAleph93();
//  SUBROUTINE TRIG_SetEcut(yEcutMin,yEcutMax,yEcutSum,yEcutProd)
  double EcutMin = 0.43865902;  // minimum energy
  double EcutMax = 0.43865902 ; // maximum energy
  double EcutSum = 0.840;  // sum of energies!!!
  double EcutProd= 0.000;  // no cut
  trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
  //
  trig_unisical_(lWW,lNN,lNW,lWN);
} // TAleph95


void TOpal93(int &lWW,  int &lNN, int &lNW,  int &lWN )
// Setter for UniSical (trig_unisical_)
//  OPAL 1993 Silicon detector
{
  double TMIND = 0.02520;
  double TMAXD = 0.05772;
  int    NPHI  = 32, NTHE = 32;
  trig_setgrid_( TMIND, TMAXD, NPHI, NTHE );
//
  int Npad= 16, Nseg = 2; // from TRIOSIW() in silicon.f
  trig_setclust_( Npad,  Nseg );
//
  double th1w,th2w, th1n, th2n;
  double PAD = (TMAXD-TMIND)/NTHE;
  th1w = TMIND + 2*PAD; // from TRIOSIW() in silicon.f
  th2w = TMAXD - 2*PAD;
  th1n = TMIND + 6*PAD;
  th2n = TMAXD - 6*PAD;
  trig_setang_(th1w,th2w, th1n,th2n);
//
  double Acopl = 0.200;       // delda phi
  double Acoli = 0.010161672; // delta theta
  trig_setacol_( Acopl, Acoli);
  ///[[[[
  // PrintAng "TOpal93", th1w, th2w,  th1n, th2n);
//
  double EcutMin = 0.500;  // minimum energy (min(E1,E2)>EcutMin)
  double EcutMax = 0.500;  // maximum energy (max(E1,E2)>EcutMax)
  double EcutSum = 0.750;  // sum of energies (E1+E2>EcutSum)
  double EcutProd= 0.000;  // no cut
  trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
  //
  trig_unisical_(lWW,lNN,lNW,lWN);
} // TOpal93


void TL3bgo2( int &lWW,  int &lNN, int &lNW,  int &lWN )
// Selection emulation using TRIG_UniSical
// L3 BGO detector 1990-95
//-----------------------
//C      CALL TRIBGO2 (32,90,0.025D0,0.070D0,Z1,Z2,Z3)
  {
  double TMIND = 0.025, TMAXD= 0.070;
  int    NPHI  = 32, NTHE = 90;
  trig_setgrid_( TMIND, TMAXD, NPHI, NTHE );
//
  int Npad= 40, Nseg = 5;
  trig_setclust_( Npad,  Nseg );
//
  double th1w,th2w, th1n, th2n;
  double PAD = (TMAXD-TMIND)/NTHE;
  th1w = TMIND +  4*PAD;
  th2w = TMAXD - 10*PAD;
  th1n = TMIND + 14*PAD;
  th2n = TMAXD - 32*PAD;
  trig_setang_(th1w,th2w, th1n,th2n);
  ///[[[[
  // PrintAng "TL3bgo2", th1w, th2w,  th1n, th2n);
//
  double Acopl = 0.3490658/2.0; // delta phi (+-10deg)
  double Acoli =   1e9; // delta theta, no cut-offs
  trig_setacol_( Acopl, Acoli);
//
  double EcutMin = 0.800;  // minimum energy (min(E1,E2)>EcutMin)
  double EcutMax = 0.400;  // maximum energy (max(E1,E2)>EcutMax)
  double EcutSum = 0.000;  // sum, no cut
  double EcutProd= 0.000;  // no cut
// additional cut xEcutMin=0.95, xEcutMax=0.20 not yet included !!!
  trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
  //
  trig_unisical_(lWW,lNN,lNW,lWN);
  } // TL3bgo2()


void TDelphi94(int &lWW,  int &lNN, int &lNW,  int &lWN )
// Selector using UniSical (trig_unisical_)
// DELPHI STIC detector from 1994
//-----------------------
//C DELPHI
//C    CALL TRISTIC2 (16,135,0.030D0,0.138D0,Z1,Z2,Z3)
{
  double TMIND = 0.030, TMAXD= 0.138;
  int    NPHI  = 16, NTHE = 135;
  trig_setgrid_( TMIND, TMAXD, NPHI, NTHE );
//
  int Npad= 38, Nseg = 3;
  trig_setclust_( Npad,  Nseg );
//
  double th1w,th2w, th1n, th2n;
  double PAD = (TMAXD-TMIND)/NTHE;
  th1w = TMIND +  9*PAD;
  th2w = TMAXD - 14*PAD;
  th1n = TMIND + 17*PAD;
  th2n = TMAXD - 31*PAD;
  trig_setang_(th1w,th2w, th1n,th2n);
  ///[[[[
  // PrintAng "TDelphi94", th1w, th2w,  th1n, th2n);
//
  double Acopl = 0.3490658; // delta phi (+-20deg)
  double Acoli =   1e9; // delta theta, no cut-offs
  trig_setacol_( Acopl, Acoli);
//
  double EcutMin = 0.650;   // minimum energy   (min(E1,E2)>EcutMin)
  double EcutMax = 0.650;   // maximum energy   (max(E1,E2)>EcutMax)
  double EcutSum = 0.000;   // sum, no cut
  double EcutProd= 0.000;   // no cut

  trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
  //
  trig_unisical_(lWW,lNN,lNW,lWN);
  } // TDelphi94()



//////////////////////////////////////////////////////////////////
//   Selectors for 1st generation lumi detectors (91->92)       //
//////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////
//*                   ALEPH LCAL (91->92)
//* Narrow range approximately   3.3-6.3 deg.
//* Energy default cut
//*     min(Ecl_1+Ecl_2)/E_beam  > 0.44 and
//*     (Ecl_1+Ecl_2)/(2*E_beam) > 0.6
//****************************************************************
//* Calling sequence:
//*     CALL TRIGAS2(0.057D0,0.107D0,0.043D0,0.125D0,16,XWIDE,XNARR,XMIX1,XMIX2)
//*     IF ( XMIX1 .LT. 0.6D0 ) GOTO 999
//****************************************************************
//      SUBROUTINE TRIG_TRIGAS2(TH1N,TH2N,TH1W,TH2W,NPHI,XWIDE,XNARR,XMIX1,XMIX2)
//C     **********************************************************
void TAlephLcal(int &lWW,  int &lNN, int &lNW,  int &lWN )
{
 double TH1N=0.057e0;
 double TH2N=0.107e0;
 double TH1W=0.043e0;
 double TH2W=0.125e0;
 ///[[[[
 // PrintAng "TAlephLcal", TH1W, TH2W,TH1N,TH2N);
 int NPHI=16;
 double XWIDE,XNARR,XMIX1,XMIX2;
 trig_trigas2_(TH1N,TH2N,TH1W,TH2W,NPHI,XWIDE,XNARR,XMIX1,XMIX2);
 lWW = XWIDE > 0.6e0;
 lNN = XNARR > 0.6e0;
 lNW = XMIX1 > 0.6e0;
 lWN = XMIX2 > 0.6e0;
}//TAlephLcal


/////////////////////////////////////////////////////////////////////////
//*----------    OPAL FCAL (91->92) -------------
//* theta min = 0.040 rad
//* theta max = 0.150 rad
//* Narrow acceptance : 0.065 -> 0.105 rad
//* Wide acceptance :   0.055 -> 0.115 rad
//* E1,E2 > 0.45 Ebeam
//* E1+E2 > 0.67 √s
//* Acoplanarity < 20 degrees
//* Granularity, clustering, Npad = 11, Nseg = 2
//*      CALL TRIFCAL1 (16,22,0.040D0,0.150D0,Z1,Z2,Z3)
//*      IF ( Z1(3) .GE. 5D0 ) GOTO 999

void TOpalFcal(int &lWW,  int &lNN, int &lNW,  int &lWN )
// Selection using UniSical (trig_unisical_)
// OPAL FCAL detector 1990-92 (emulation using TRIG_UniSical)
{
  double TMIND = 0.040, TMAXD= 0.150;
  int    NPHI  = 16, NTHE = 22;
  trig_setgrid_( TMIND, TMAXD, NPHI, NTHE );
//
  int Npad= 11, Nseg = 2;
  trig_setclust_( Npad,  Nseg );
//
  double th1w,th2w, th1n, th2n;
  double PAD = (TMAXD-TMIND)/NTHE;
  th1w = TMIND +  3*PAD;
  th2w = TMAXD -  7*PAD;
  th1n = TMIND +  5*PAD;
  th2n = TMAXD -  9*PAD;
  trig_setang_(th1w,th2w, th1n,th2n);
  ///[[[[
  // PrintAng "TOpalFcal", th1w, th2w,  th1n, th2n);
//
  double Acopl = 0.3490658; // delta phi (+-20deg)
  double Acoli =   1e9;     // delta theta, no cut-offs
  trig_setacol_( Acopl, Acoli);
 //
 //      PARAMETER ( Ecut = 0.45,    Scut = 0.67 )
  double EcutMin = 0.450;  // minimum energy (min(E1,E2)>EcutMin)
  double EcutMax = 0.000;  // maximum energy, no cut
  double EcutSum = 0.670;  // sum, no cut
  double EcutProd= 0.000;  // no cut
  trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
  //
  trig_unisical_(lWW,lNN,lNW,lWN);
} // TOpalFcal()


//C--------   L3 BGO (1) (90->92)  ------------
//* theta min = 0.0252 rad
//* theta max = 0.0712 rad
//* Narrow acceptance : 0.0312 -> 0.0652 rad
//* Wide acceptance :   0.0252 -> 0.0712 rad
//* (no specific cut on the wide side)
//* Emin > 0.4 Ebeam
//* Emax > 0.8 Ebeam
//* Acoplanariity < 10 degrees
//* PARAMETER ( Acopl = 0.3490658D0/2D0 )
//*-------
//* CALL TRIBGO1 (16,23,0.0252D0,0.0712D0,Z1,Z2,Z3)
//* IF ( Z1(3) .GE. 5D0 ) GOTO 999
//*-------
//* Npad = 11, Nseg = 2
//* PARAMETER ( Ecut = 0.,        Scut = 0. )
//* Lecut= DMAX1(ECL1,ECL2)/Ebeam.GE.0.80 .AND.
//*&       DMIN1(ECL1,ECL2)/Ebeam.GE.0.40
//* Lscut= (ECL1+ECL2)/(2D0*Ebeam).GE.Scut
    void TL3bgo1(int &lWW,  int &lNN, int &lNW,  int &lWN )
    // Selection using UniSical (trig_unisical_)
    // L3 BGO detector 1990-92 (emulation using TRIG_UniSical)
    {
      double TMIND = 0.0252, TMAXD= 0.0712;
      int    NPHI  = 16, NTHE = 23;
      trig_setgrid_( TMIND, TMAXD, NPHI, NTHE );
   //
      int Npad= 11, Nseg = 2;
      trig_setclust_( Npad,  Nseg );
   //
      double th1w,th2w, th1n, th2n;
      double PAD = (TMAXD-TMIND)/NTHE;
      th1w = TMIND +  0*PAD;
      th2w = TMAXD -  0*PAD;
      th1n = TMIND +  3*PAD;
      th2n = TMAXD -  3*PAD;
      trig_setang_(th1w,th2w, th1n,th2n);
      ///[[[[
      // PrintAng "TL3bgo1", th1w, th2w,  th1n, th2n);
   //
      double Acopl = 0.3490658/2.0; // delta phi (+-10deg)
      double Acoli =   1e9;         // delta theta, no cut-offs
      trig_setacol_( Acopl, Acoli);
   //
   //  PARAMETER ( Ecut = 0.45,   Scut = 0.67 ) ???
      double EcutMin = 0.400;  // minimum energy (min(E1,E2)>EcutMin)
      double EcutMax = 0.800;  // maximum energy (max(E1,E2)>EcutMax)
      double EcutSum = 0.000;  // sum
      double EcutProd= 0.000;  // no cut
      trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
  //
      trig_unisical_(lWW,lNN,lNW,lWN);
} // TL3bgo1()


////////////////////////////////////////////////////////////
//*-----------    DELPHI SAT (90->93) -----------
//* The calling and testing sequence is also in the file
//* theta min = 0.4282 rad ??? 0.04282 !!!
//* theta max = 0.1452 rad
//* Narrow acceptance : 0.05602 -> 0.12862 rad
//* Wide acceptance :   0.05272 -> 0.14182 rad
//* E1, E2 > 0.65 Ebeam
//* Acoplanarity < 20 degrees
//* Granularity, clustering: Npad = 13, Nseg = 2
//*      CALL TRISAT1 (16,31,0.04282D0,0.14512D0,Z1,Z2,Z3)
//*      IF ( Z1(3) .GE. 5D0 ) GOTO 999

void TDelphiSat(int &lWW,  int &lNN, int &lNW,  int &lWN )
   // Selection using UniSical (trig_unisical_)
   // DELPHI SAT detector 1990-92 (emulation using TRIG_UniSical)
  {
     double TMIND = 0.04282, TMAXD= 0.14512;
     int    NPHI  = 16, NTHE = 31;
     trig_setgrid_( TMIND, TMAXD, NPHI, NTHE );
  //
     int Npad= 13, Nseg = 2;
     trig_setclust_( Npad,  Nseg );
  //
     double th1w,th2w, th1n, th2n;
     double PAD = (TMAXD-TMIND)/NTHE;
     th1w = TMIND +  3*PAD;
     th2w = TMAXD -  1*PAD;
     th1n = TMIND +  4*PAD;
     th2n = TMAXD -  5*PAD;
     trig_setang_(th1w,th2w, th1n,th2n);
     ///[[[[
     // PrintAng "TDelphiSat", th1w, th2w,  th1n, th2n);
    //
     double Acopl = 0.3490658; // delta phi (+-20deg)
     double Acoli =   1e9;     // delta theta, no cut-offs
     trig_setacol_( Acopl, Acoli);
//
//  PARAMETER ( Ecut = 0.65,   Scut = 0.0 )
     double EcutMin = 0.650;   // minimum energy (min(E1,E2)>EcutMin)
     double EcutMax = 0.000;   // maximum,energy no cut
     double EcutSum = 0.000;   // sum, no cut
     double EcutProd= 0.000;   // no cut
     trig_setecut_( EcutMin,EcutMax,EcutSum, EcutProd );
     //
     trig_unisical_(lWW,lNN,lNW,lWN);
    //
  } // TDelhiSat()


 ////////////////////////////////////////////////////////////////////////////
   ClassDef(TRig,1);  //
 };
 ////////////////////////////////////////////////////////////////////////////
 #endif
