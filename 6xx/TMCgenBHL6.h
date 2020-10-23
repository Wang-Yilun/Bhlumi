#ifndef TMCgenBHL6_H
#define TMCgenBHL6_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   TMCgenBHL6                                          //
//                                                                          //
//            Interface (wrapper)  to MC event generator BHLUMI             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

/// C++ headers
using namespace std;
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>


/// ROOT headers
#include "TRandom.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
#include "TLorentzVector.h"

/// My headers
#include"BXFORMAT.h"


//////////////////////////////////////////////////////////////////////////
//      SUBROUTINE BHLUM4(MODE,XPAR,NPAR)
extern "C"  void bhlum4_( const int&, double xpar[], int npar[]);
//      SUBROUTINE BHL_ReModel(WtMain)
extern "C"  void bhl_remodel_(const double& );
//ENTRY      MARINI(IJKLIN, NTOTIN,NTOT2N)
extern "C" void marini_( const int&, const int&, const int&);
//       SUBROUTINE BHCLEN
extern "C" void bhclean();
//      SUBROUTINE BHL_fort_open(nout,fname)
extern "C" void bhl_fort_open_( const int&, const char*, int);
//      SUBROUTINE BHL_fort_close(nout,fname)
extern "C" void bhl_fort_close_(const int&);
//      SUBROUTINE BHL_Get_PrimNorm( nevgen, XsPrim)
extern "C" void bhl_get_primnorm_(const int& , const double& );
//      SUBROUTINE BHL_GetPhoton1(iphot,phot)
extern "C" void bhl_getphoton1_( int&, double []);
//      SUBROUTINE BHL_GetWt(WtMainX,WtCrudX)
extern "C" void bhl_getwt_(const double& , const double&);
//      SUBROUTINE BHL_GetBeams(p1,p2)
extern "C" void bhl_getbeams_(   double [], double []);
//      SUBROUTINE BHL_GetFermions(q1,q2)
extern "C" void bhl_getfermions_(   double [], double []);
//      SUBROUTINE BHL_GetNphot(nphot5)
extern "C" void bhl_getnphot_( const int&);
//      SUBROUTINE BHL_GetWtAlter(j,WtAlter)
extern "C" void bhl_getwtalter_(const int& , const double& );
//      SUBROUTINE BHL_SetKeyPia(KeyPia1)
extern "C" void bhl_setkeypia_( const int&);
//      SUBROUTINE BHL_SetKeyPib(KeyPib2)
extern "C" void bhl_setkeypib_( const int&);
//      SUBROUTINE BHL_SetKeyZet(KeyZet1)
extern "C" void bhl_setkeyzet_( const int&);
//      SUBROUTINE BHL_SetKeyMod(KeyMod2)
extern "C" void bhl_setkeymod_( const int&);

//      SUBROUTINE DUMPS(NOUT)
extern "C" void dumps_( const int&);
//      SUBROUTINE glimit(lenmx)
extern "C" void glimit_( const int&);
//      SUBROUTINE BHL_SetInOut(NINPx,NOUTx)
extern "C" void bhl_setinout_( const int&, const int& );
//      FUNCTION BORNB(CMSENE,THA,THB)
extern "C" double bornb_( const double&, const double&, const double&);
//      FUNCTION BORNG(THA,THB)
extern "C" double borng_( const double&, const double&);
//      SUBROUTINE  vacpol(KeyPia,Q2,SINW2,RePiE,dRePiE)
extern "C" void vacpol_( const int&, const double&, const double&, const double&, const double&);



//////////////////////////////////////////////////////////////////////////
//________________________________________________________________________
class TMCgenBHL6: public TObject {
public:
 char  f_Name[64];       // Name of a give instance of the class
 float f_Version;         // Actual VERSION of the program
 char  f_Date[40];        // Release DATE of the program
 //-------------- USER Input parameters--------------------------
public:
 /// Engines and services
 TRandom  *f_RNgen;            //  External RN event generator
 TH1D     *f_TMCgen_NORMA;     //! special histo keeping overall normalization
  /// data members
public:
 int       f_IsInitialized;    //  prevents repeating initialization
 TFile    *f_GenFile;          //! ROOT DiskFile with MC generator object and data
 TFile    *f_HstFile;          //! ROOT DiskFile with all histos
 ofstream *f_Out;              //! External Logfile for messages
 double    f_NevGen;           //  event serial number
////////////////////////////////////////
 public:
  long      m_NevTot;         //   total numer of events to be generated
  long      m_EvenCounter;    //   event serial counter
  int       m_jmax;           //   lenght of ymar
  double    m_ypar[10001];    //   xpar input params for BHLUMI4
  int       m_npar[10001];    //   npar input params for BHLUMI4

  double    m_XsNormPb;       //   normalization
  double    m_XsErroPb;       //   normalization
  // Seeds for Pseumar r.n. in KKMC
  int       m_ijkl_new;       // = 54217137;
  int       m_ntot_new;       // = 0;
  int       m_ntot2_new;      // = 0;
  long      m_out;            //   output unit number

//------------------------------------------------------------------------
//                    BHLUMI4 input params
//------------------------------------------------------------------------
//  NPAR( 1)  KeyOpt =1000*KeyGen +100*KeyRem +10*KeyWgt +KeyRnd
//            KeyGen =3 for this sub-generator
//            KeyRem =0,1 removal/no-removal switch, both are OK,
//                   =1 no-removal technicaly simpler/safer
//                   =0 OBLIGATORY for KeyZet =1 !!!
//            KeyRnd =1,2 type of random number generator RANMAR,RANECU
//                   =1 better for parallel production on computer farm.
//            KeyWgt =0,1,2 for constant/variable weight WTMOD,
//                   =0, WTMOD =1 useful for apparatus Monte Carlo.
//                   =1, WTMOD variable, option faster/safer, RECOMMENDED
//                   =2, WTMOD variable, events below trmin generated
  int         m_KeyGen;
  int         m_KeyRem;
  int         m_KeyWgt;
  int         m_KeyRnd;
//  NPAR( 2)  KeyRad =1000*KeyZet+100*KeyUpd+10*KeyMod +KeyPia
//            KeyZet =0,1 test switch,
//                   =0   Z contribution OFF
//                   =1   Z contribution ON, DEFAULT!!!
//            KeyUpd =0,1,2 test switch,
//                   =0 normal position DEFAULT!!!
//                   =1 only upper line bremss., =2 only lower line
//            KeyMod =1,2 type of MODEL subrogram and QED matrix element
//                   =1 version compatible with Comp. Phys. Comm. 70 (1992) 305
//                   =2 version 4.x which is now DEFAULT!
//            KeyPia =0,1,2,3 photon vacuum polarization and s-chanel photon
//                   =0 OFF, it used in semianalytical tests,
//                   =1 ON,  Burkhardt et.al. 1989, as in BHLUMI 2.0x
//                   =2 ON,  S. Eidelman, F. Jegerlehner, Z. Phys. C (1995)
//                   =3 ON,  Burkhardt and Pietrzyk 1995 (Moriond).
  int         m_KeyZet;
  int         m_KeyUpd;
  int         m_KeyMod;
  int         m_KeyPia;
//  XPAR( 1)  CMSENE Total center mass energy [GeV]
//  XPAR( 2)  TRMIN Minimum transfer (positive) [GeV**2]
//  XPAR( 3)  TRMAX Maximum transfer (positive) [GeV**2]
//  XPAR( 4)  EPSCM Dimensionless infrared cut on CMS energy of soft
//                   photons, ( E_phot > CMSENE*EPSCM/2 )
  double      m_CMSENE;
  double      m_angmin;
  double      m_angmax;
  double      m_TRMIN;
  double      m_TRMAX;
  double      m_EPSCM;
  double      m_VMAXG;     //  = 0.9999e0, VMAXG   v_max generation
//
 public:
 TMCgenBHL6();                // explicit default constructor for streamer
 TMCgenBHL6(const char*);     // user constructor
 ~TMCgenBHL6();               // explicit destructor
 public:
/////////////////////////////////////////////////////////////////////////////
/// methods obligatory
  void Initialize(TRandom*, ofstream*, TH1D*);
  void ReInitialize(        ofstream*, TH1D*);
  void Finalize();
  void Generate();

  double sqr( const double x ){ return x*x;};
  int Min( int &i, int &j){ if(i<j) return i; else return j;};
/// for debug
  void StopM(char* message){
    cout <<"++++ TMCgen: "<< message << endl; exit(5);}    //Error message
///  getters
  int GetIsNewRun() const
    { if(f_IsInitialized == 0 ) return 1; else return 0;}
///
///////////////////////////////////////////////////////////////////////////////
  void GetPhoton1(const int iphot, TLorentzVector &phot);
  void GetWt( double&,  double&);
  void GetBeams(    TLorentzVector&,  TLorentzVector&);
  void GetFermions( TLorentzVector&,  TLorentzVector&);
  void GetPrimaNorma( int &, double &);
  void GetNphot(  int &);
  double GetWtAlter(const int );
  void Print1();
  void SetInOut( const int &, const int& );

  //      SUBROUTINE MODEL(MODE,WTM)
  void ReModel( const double &WtMain){
     bhl_remodel_( WtMain);
  }

  //      FUNCTION BORNB(CMSENE,THA,THB)  // pure t-chanel
  double BornBhl( double &cmsene,  double thmin,  double thmax)
      { return bornb_(cmsene,thmin,thmax);  };
  // Born with vacuum polarization, requires initializaction of Bhlumi4
  //      FUNCTION BORNG(THA,THB)
  double BornBhl2( double thmin,  double thmax)
      { return borng_(thmin,thmax);  };
  //
  void SetKeyPia( const int & KeyPia)
      {  bhl_setkeypia_(KeyPia);  };
  //
  void SetKeyZet( const int & KeyZet)
      {  bhl_setkeyzet_(KeyZet);  };
  //
  void SetKeyPib( const int & KeyPib)
  // this is for VP error studies (shift of VP value)
      {  bhl_setkeypib_(KeyPib);  };

  //    SUBROUTINE BHL_SetKeyMod(KeyMod1)
  void SetKeyMod( const int & KeyMod)
  // reseting model type KeyMod=1,2 (old.new)
       {  bhl_setkeymod_(KeyMod);  };


////////////////////////////////////////////////////////////////////////////
  ClassDef(TMCgenBHL6,1); // Monte Carlo generator
};
////////////////////////////////////////////////////////////////////////////
#endif
