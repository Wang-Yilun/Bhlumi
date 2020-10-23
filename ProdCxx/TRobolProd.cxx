///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS ROBOL                                                //
//                                                                           //
//     Three-stroke engine for analysis: Initialize-Production-Finalize      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TRobolProd.h"
#include "TRig.h"

#define SP12 setprecision(7) << setw(12)
#define SP18 setprecision(9) << setw(18)
#define SW5  setw(5)
#define SW9  setw(9)
#define SW12 setw(12)
#define SW18 setw(18)
#define SP21 setw(21)<<setprecision(13)
#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(5)


ClassImp(TRobolProd);


TRobolProd::TRobolProd()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "@@@@> TRobolProd DEFAULT Constructor (for ROOT only) "<<endl;
  f_RNgen    = NULL;
  f_MCgen    = NULL;
  f_HstFile  = NULL;
  f_GenFile  = NULL;
  f_Out      = NULL;
  f_TRigLib  = NULL;
}


///_____________________________________________________________
TRobolProd::TRobolProd(const char* Name)
{
//! Constructor to be used by the user!!!
//! Its important role is to define ALL DEFAULTS.
//! to changed by the user before calling TMCgen::Initialize
///
  cout<< "@@@@> TRobolProd::TRobolProd USER Constructor "<<endl;
///
  if(strlen(Name)  >65) StopM("TRobolProd::TRobolProd: +++ stop, Name too long");
  sprintf(f_Name,"%s",Name);         // Class name
/// From base class
  f_HistNormName = "h_TMCgen_NORMA";
  f_RNgen    = NULL;
  f_MCgen    = NULL;
  f_HstFile  = NULL;
  f_GenFile  = NULL;
  f_Out      = NULL;
  f_TRigLib  = NULL;
/////////////////////////////////////////////////////////
  m_NevGen=0;  // event counter
  m_count1=0;
  m_pi      = 3.1415926535897932;
  m_NEVTOT = 84000000;    //  total no. of MC events to generate
// ==================  params of event selection   ==========================
  m_VMAXE  =    0.5e0;    //  VMAXE   v_max              trigger maximum v
  m_NPHI   =       32;    //  no of phi   sect. CALO2 and SICAL2
  m_NTHE   =       16;    //  no of theta sect. CALO2 and SICAL2
//
//  Fidutial angular range [mrad] narrow and wide
  m_TminF  =     .024;    //  fidutial  theta_min CALO2 and SICAL2
  m_TmaxF  =     .058;    //  fidutial  theta_max CALO2 and SICAL2
  m_dcone_calo1  =  0.010e0;    //  radius of the cone  CALO1
}///TRobolFOAM


///______________________________________________________________________________________
TRobolProd::~TRobolProd()
{
  //!Explicit destructor
  cout<< "@@@@> TRobolProd::TRobolProd !!!! DESTRUCTOR !!!! "<<endl;
}///destructor


//______________________________________________________________________________
//////////////////////////////////////////////////////////////
//   Initialize MC generator and analysis programs          //
//////////////////////////////////////////////////////////////
void TRobolProd::Initialize(
        ofstream *OutFile, /// Central log-file for messages
        TFile *GenFile,    /// ROOT disk file for CRNG and MC gen.
        TFile *HstFile)    /// ROOT disk file for histograms
{
  cout<< "****> TRobolProd::Initialize starts"<<endl;
  //////////////////////////////////////////////////////////////
  //   Initialize MC generator and analysis programs          //
  //   This is part formerly in MCdev base class               //
  //////////////////////////////////////////////////////////////
  //
  f_Out        = OutFile;
  f_GenFile    = GenFile;
  f_HstFile    = HstFile;
  //
  f_NevGen=0;
  f_count1=0;
  //
  f_TRigLib  = new TRig("TRigLib");
  f_TRigLib->Initialize();


  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======    TRobolProd::Initialize  ======");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  ///
  /// book histogram keeping track of overall normalization
  f_HstFile->cd(); /// to keep h_TMCgen_NORMA in  f_HstFile!!!
  TH1D *h_TMCgen_NORMA = TH1D_UP(f_HistNormName,"Normalization histo",10000,0,10000);
  //TH1D *h_TMCgen_NORMA = new TH1D("h_TMCgen_NORMA","Normalization histo",10000,0,10000);
  ///
  /// Read prepared object MC of generator from the disk
  f_GenFile->cd();
//  f_MCgen = (TMCgen*)f_GenFile->Get("MCgen");
  f_RNgen =(TRandom*)f_GenFile->Get("RN_gen");  // read r.n. generator
  ///
  f_MCgen = (TMCgenBHL6*)f_GenFile->Get("MCgen");
  /// This is virgin run if MC generator NOT initialized
  f_isNewRun = f_MCgen->GetIsNewRun();
  ///
  /// Initialize/Redress MC generator
  if( f_isNewRun == 1){
    f_MCgen->f_GenFile = f_GenFile; /// Adding access to disk files,
    f_MCgen->f_HstFile = f_HstFile; /// just in case it is needed.
    f_MCgen->Initialize(f_RNgen,f_Out,h_TMCgen_NORMA);
  } else {
	StopM( "Redress not implemented yet" );
    f_MCgen->f_GenFile = f_GenFile;
//    f_MCgen->Redress(   f_RNgen,f_Out,h_TMCgen_NORMA);
    f_MCgen->f_HstFile = f_HstFile;
  }
  /// and write MC generator into disk file.
  f_MCgen->Write("MCgen",TObject::kOverwrite);   // MC generator status
///////////////////////////////
// End of MCdev import
///////////////////////////////////////////////////////////////////
//
// ****************************************************************
// *************** TRICAL1 *****************
// LCAL type trigger with theta-rings and phi-sectors
// Descendant of TRIGAS0, see zcalA,
// zcalA for BERE1, zcalB for CALO1, zcalC for CALO2
// CALO1: Associate photon+electrons wthin delcon cone
// CALO2: Associate photon+electrons wthin dlthe*dlphi plaquette
// ----------------------------------------------------------------------
   double padthe = (m_TmaxF-m_TminF)/m_NTHE;
   double padphi = 2*m_pi/m_NPHI;
   m_th1w_calo2  = m_TminF   +padthe; // CALO2
   m_th2w_calo2  = m_TmaxF   -padthe; // CALO2
   m_th1n_calo2  = m_TminF +2*padthe; // CALO2
   m_th2n_calo2  = m_TmaxF -4*padthe; // CALO2
   m_dlph_calo2  = 1.5e0 *padphi;     // half-size of theta window CALO2
   m_dlth_calo2  = 1.5e0 *padthe;     // half-size of phi   window CALO2
   cout<< "++++TRobolProd::Initialize: th1n= "
		   <<m_th1n_calo2<<" th2n ="<< m_th2n_calo2<<"  th1w="<<m_th1w_calo2<<"  th2w="<<m_th2w_calo2<<endl;

/// Book histograms or read them from the disk
  Hbooker();
////////////////////////
  BX1I(*f_Out," isNewRun", f_isNewRun," is it new MC run?       =");
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"====== END of TRobolProd::Initialize ===");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  cout<< "####> TRobolProd::Initialize FINISHED"<<endl;

}///Initialize

///////////////////////////////////////////////////////////////////////////////
void TRobolProd::Hbooker()
{
  ///
  cout<< "****> TRobolFOAM::Hbooker: histogram booking STARTS"<<endl;
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======    TRobol::Hbooker    ===========");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  f_HstFile->cd();
  //  ************* user histograms  *************
  int nbin=1000;
  hst_weight   = TH1D_UP("hst_weight" , "MC weight         ", 100,   0.000 ,   2.0);
  hst_weight2  = TH1D_UP("hst_weight2", "MC weight no trig.", 250, -25.000 , 100.0);
  hst_wt2trig  = TH1D_UP("hst_wt2trig", "MC weight CALO2   ", 250, -25.000 , 100.0);

  hst_Costhe   = TH1D_UP("hst_Costhe",  "electr. Costhe", 100,   0.000 , 0.100);

  hst_v_bare1_ww = TH1D_UP("hst_v_bare1_ww",   "zmax bare1",  100,  0.0 , 1.0);
  hst_v_bare1_nw = TH1D_UP("hst_v_bare1_nw",   "zmax bare1",  100,  0.0 , 1.0);
  hst_v_bare1_nn = TH1D_UP("hst_v_bare1_nn",   "zmax bare1",  100,  0.0 , 1.0);

  hst_v_calo2_ww = TH1D_UP("hst_v_calo2_ww",   "zmax calo2",  100,  0.0 , 1.0);
  hst_v_calo2_nw = TH1D_UP("hst_v_calo2_nw",   "zmax calo2",  100,  0.0 , 1.0);
  hst_v_calo2_nn = TH1D_UP("hst_v_calo2_nn",   "zmax calo2",  100,  0.0 , 1.0);

  hst_v_sical2_ww = TH1D_UP("hst_v_sical2_ww", "zmax sical2", 100,  0.0 , 1.0);
  hst_v_sical2_nw = TH1D_UP("hst_v_sical2_nw", "zmax sical2", 100,  0.0 , 1.0);
  hst_v_sical2_nn = TH1D_UP("hst_v_sical2_nn", "zmax sical2", 100,  0.0 , 1.0);
//
  hst_trig_VP0 = TH1D_UP("hst_trig_VP0", "xsect[trig] VP on ", 100, 0.0, 100.0);
  hst_trig_VP1 = TH1D_UP("hst_trig_VP1", "xsect[trig] VP off", 100, 0.0, 100.0);
  hst_trig_VP2 = TH1D_UP("hst_trig_VP2", "xsect[trig] VP err", 100, 0.0, 100.0);
  hst_trig_VP3 = TH1D_UP("hst_trig_VP3", "xsect[trig] VP off", 100, 0.0, 100.0);
  hst_trig_VP4 = TH1D_UP("hst_trig_VP4", "xsect[trig] VP err", 100, 0.0, 100.0);
  hst_trig_VP5 = TH1D_UP("hst_trig_VP5", "xsect[trig] VP dif", 100, 0.0, 100.0);
//
  hst_trig_ZZ0 = TH1D_UP("hst_trig_ZZ0", "xsect[trig] Z off ",  100, 0.0, 100.0);
  hst_trig_ZZ1 = TH1D_UP("hst_trig_ZZ1", "xsect[trig] Z on  ",  100, 0.0, 100.0);
  hst_trig_ZZ2 = TH1D_UP("hst_trig_ZZ2", "xs[trig] Zonly exp.", 100, 0.0, 100.0);
  hst_trig_ZZ3 = TH1D_UP("hst_trig_ZZ3", "xs[trig] Zonly O(alf1)", 100, 0.0, 100.0);
  hst_trig_ZZ4 = TH1D_UP("hst_trig_ZZ4", "xs[trig] Zonly O(alf1)", 100, 0.0, 100.0);

 }// Hbooker


///////////////////////////////////////////////////////////////////////////////
void TRobolProd::Production(double &iEvent)
{
///////////////////////////////////////////////////////////////////////////////
//
//   GENERATE AND ANALYZE SINGLE EVENT
//
///////////////////////////////////////////////////////////////////////////////
  // ****************************************************************
  // ************ Generate event and import it here  ****************
  m_NevGen++;
  TMCgenBHL6 *BHL_generator = (TMCgenBHL6*)f_MCgen;
  BHL_generator->Generate(); // done in user class

  /// *************************************************
  double WtMain,WtCrude;
  BHL_generator->GetWt(WtMain,WtCrude);
  BHL_generator->GetBeams(   m_pbea1,m_pbea2);
  BHL_generator->GetFermions(m_pfer1,m_pfer2);

  BHL_generator->GetNphot(m_Nphot);                    // photon multiplicity
  TLorentzVector VSumPhot;    // By default all components are initialized by zero.
  long iphot,iphot1;
  for(iphot=0;iphot<m_Nphot;iphot++){
    BHL_generator->GetPhoton1(iphot+1,m_phot[iphot]);  // photon 4-momenta
    VSumPhot+= m_phot[iphot];
  }
  if(iEvent<=1){
    cout<<"-----------------------------------------------------------  "<<iEvent;
    cout<<"  -----------------------------------------------------------"<<endl;
    cout<<"  WtMain= "<< WtMain  << "  WtCrude= "<< WtCrude<<endl;
    cout<<" m_Nphot= "<< m_Nphot<<endl;
    cout<<"VSumPhot= "; MomPrint( VSumPhot );
    BHL_generator->Print1();
    (*f_Out)<< "  WtMain= "<< WtMain  << "  WtCrude= "<< WtCrude<<endl;
  }
// ****************************************************************
  double s  =(m_pbea1+m_pbea2)*(m_pbea1+m_pbea2);
  double s1 =(m_pfer1+m_pfer2)*(m_pfer1+m_pfer2);
  double CMSene = sqrt(s);
  double Mff    = sqrt(s1);
  double vv     = 1-s1/s;

  double Costh1 = m_pfer1.CosTheta();
  double Costh2 = m_pfer2.CosTheta();


// ******************************************************************
// ///////////////////  histogramming
  hst_weight->Fill(WtMain);
  hst_weight2->Fill(WtMain);
  if( 1-m_y1C[0] < 0.5 ) hst_wt2trig->Fill(WtMain);

  hst_Costhe->Fill(Costh1,WtMain);
  hst_Costhe->Fill(Costh2,WtMain);

////////////////////////////////////////////////////////////////
//  Event selection BARE1
  double Th1n = m_TminF+(m_TmaxF-m_TminF)/16;
  double Th2n = m_TmaxF-(m_TmaxF-m_TminF)/16;
  f_TRigLib->TriCal1(Th1n,Th2n,  m_TminF,m_TmaxF,              // input
    		             m_dcone_calo1,m_dlth_calo2,m_dlph_calo2,  // input
    		             m_y1A,m_y1B,m_y1C);    // output: BARE1, CALO1,CALO2
  hst_v_bare1_ww->Fill(1-m_y1A[0]+1e-5,WtMain);
  hst_v_bare1_nn->Fill(1-m_y1A[1]+1e-5,WtMain);
  hst_v_bare1_nw->Fill(1-m_y1A[2]+1e-5,WtMain/2);
  hst_v_bare1_nw->Fill(1-m_y1A[3]+1e-5,WtMain/2);

////////////////////////////////////////////////////////////////


//  Event selection CALO2
  f_TRigLib->TriCal1(m_th1n_calo2, m_th2n_calo2,m_th1w_calo2,m_th2w_calo2, // input
  		               m_dcone_calo1,m_dlth_calo2,m_dlph_calo2,                // input
  		               m_y1A,m_y1B,m_y1C);    // output: BARE1, CALO1,CALO2
  if(iEvent<=20){
  (*f_Out)<< " NevGen = "<< m_NevGen;
  (*f_Out)<< "  y1C= "<< m_y1C[0]  <<" "<< m_y1C[1] <<" "<< m_y1C[2] <<" "<< m_y1C[3] <<endl;
  }
  hst_v_calo2_ww->Fill(1-m_y1C[0]+1e-5,WtMain);
  hst_v_calo2_nn->Fill(1-m_y1C[1]+1e-5,WtMain);
  hst_v_calo2_nw->Fill(1-m_y1C[2]+1e-5,WtMain/2);
  hst_v_calo2_nw->Fill(1-m_y1C[3]+1e-5,WtMain/2);

////////////////////////////////////////////////////////////////
//  Event selections SICAL2 (simplified for wshop96)
  f_TRigLib->TriSic2w( m_NPHI,  m_NTHE,  m_TminF, m_TmaxF,
   		               m_Z1,    m_Z2,    m_Z3); // output: SICAL, 3 types of energy variable
  hst_v_sical2_ww->Fill(1-m_Z1[0]+1e-5,WtMain);
  hst_v_sical2_nn->Fill(1-m_Z1[1]+1e-5,WtMain);
  hst_v_sical2_nw->Fill(1-m_Z1[2]+1e-5,WtMain/2);
  hst_v_sical2_nw->Fill(1-m_Z1[3]+1e-5,WtMain/2);

///////////////////////////////////////////
// New configurable silicon-type trigger
// Testing BHL_UniSical
/*
  f_TRigLib->SetTrisic2w(); // reset cutoff parameters
  int lWW,lNN,lNW,lWN;
//SUBROUTINE BHL_UniSical(lWW,lNN,lNW,lWN)
  f_TRigLib->UniSical(lWW,lNN,lNW,lWN);

  int kWW=0,kNN=0,kNW=0,kWN=0;
  if( m_Z1[0]>0.5 ) kWW=1;
  if( m_Z1[1]>0.5 ) kNN=1;
  if( m_Z1[2]>0.5 ) kWN=1;
  if( m_Z1[3]>0.5 ) kNW=1;

  if(iEvent<=10000 && (kWW!=lWW || kNN!= lNN || kWN!=lWN || kNW!=lNW) ){
	  cout<<"--------------------------------------------------------------------------"<<endl;
	  cout<< "  iEvent="<< iEvent   << "  WtMain="<< WtMain<<endl;
      cout<< "  kWW,kNN,kNW,kWN="<< kWW << "  "<< kNN << "  "<< kNW << "  "<< kWN << "  "<<endl;
      cout<< "  lWW,lNN,lNW,lWN="<< lWW << "  "<< lNN << "  "<< lNW << "  "<< lWN << "  "<<endl;
   }
*/

// KeyMod=1, Type of ME. Version as in CPC 2.01 (1991,92)
// KeyMod=2, Type of ME. New version (1993, 1994, 1995)
int KeyPia, KeyPib, KeyMod, KeyZet;
BHL_generator->GetWt(WtMain,WtCrude);
double WtMain100, WtMain200, WtMain110, WtMain111, WtMain230, WtMain231, WtMain240, WtMain241;
if( WtCrude !=0 ){
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  KeyZet = 0;   // Z is OFF
  BHL_generator->SetKeyZet(KeyZet);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//======================
  KeyMod = 1;   // 2.x QED m.e. model 1990-92
  KeyPia = 0;   // VP OFF
  KeyPib = 0;   // VP shift by the error value OFF
  BHL_generator->SetKeyMod(KeyMod);   // steering param. in matrix element
  BHL_generator->SetKeyPia(KeyPia);   // steering param. in matrix element
  BHL_generator->SetKeyPib(KeyPib);   // steering param. in matrix element
  BHL_generator->ReModel(WtMain100);  // recalculating matrix element
  //
  KeyMod = 2;   // 4.x QED m.e. model (1993, 1994, 1995)
  KeyPia = 0;   // VP OFF
  KeyPib = 0;   // VP shift by the error value OFF
  BHL_generator->SetKeyMod(KeyMod);   // steering param. in matrix element
  BHL_generator->SetKeyPia(KeyPia);   // steering param. in matrix element
  BHL_generator->SetKeyPib(KeyPib);   // steering param. in matrix element
  BHL_generator->ReModel(WtMain200);  // recalculating matrix element
//======================
// 2.x+VPold89
  KeyMod = 1;   // 2.x QED m.e. model from 1990-92
  KeyPia = 1;   // VP of Burkhardt 1989 !!!
  KeyPib = 0;   // VP shift by the error value OFF
//
  BHL_generator->SetKeyMod(KeyMod);      // steering param. in matrix element
  BHL_generator->SetKeyPia(KeyPia);      // steering param. in matrix element
  BHL_generator->SetKeyPib(KeyPib);      // steering param. in matrix element
  BHL_generator->ReModel(WtMain110);     // recalculating matrix element
  // 2.x+(VPold89+err)
  KeyMod = 1;   // 2.x QED m.e. model from 1990-92
  KeyPia = 1;   // VP of Burkhardt 1989 !!!
  KeyPib = 1;   // VP SHIFTED by the error value ON
  BHL_generator->SetKeyMod(KeyMod);      // steering param. in matrix element
  BHL_generator->SetKeyPia(KeyPia);      // steering param. in matrix element
  BHL_generator->SetKeyPib(KeyPib);      // steering param. in matrix element
  BHL_generator->ReModel(WtMain111);     // recalculating matrix element
//======================
// 4.x+VPold
  KeyMod = 2;   // 4.x QED m.e. model (1993, 1994, 1995)
  KeyPia = 3;   // VP of Burkhardt and Pietrzyk 1995 (Moriond, PLB96)
  KeyPib = 0;   // VP shifted by the error value OFF
  BHL_generator->SetKeyMod(KeyMod);      // steering param. in matrix element
  BHL_generator->SetKeyPia(KeyPia);      // steering param. in matrix element
  BHL_generator->SetKeyPib(KeyPib);      // steering param. in matrix element
  BHL_generator->ReModel(WtMain230);     // recalculating matrix element
  // 4.x+(VPold+err)
  KeyMod = 2;   // 4.x QED m.e. model (1993, 1994, 1995)
  KeyPia = 3;   // VP of Burkhardt and Pietrzyk 1995 (Moriond, PLB96)
  KeyPib = 1;   // VP SHIFTED by the error value ON
  BHL_generator->SetKeyMod(KeyMod);      // steering param. in matrix element
  BHL_generator->SetKeyPia(KeyPia);      // steering param. in matrix element
  BHL_generator->SetKeyPib(KeyPib);      // steering param. in matrix element
  BHL_generator->ReModel(WtMain231);     // recalculating matrix element
 //======================
 // 4.x+VPnew
  KeyMod = 2;   // 4.x QED m.e. model (1993, 1994, 1995)
  KeyPia = 4;   // VP of Jegerlehner 2019
//  KeyPia = 5;   // VP of Teubner 2018
  KeyPib = 0;   // VP shift by the error value OFF
  BHL_generator->SetKeyMod(KeyMod);      // steering param. in matrix element
  BHL_generator->SetKeyPia(KeyPia);      // steering param. in matrix element
  BHL_generator->SetKeyPib(KeyPib);      // steering param. in matrix element
  BHL_generator->ReModel(WtMain240);     // recalculating matrix element
  // 4.x+(VPnew+err)
  KeyMod = 2;   // 4.x QED m.e. model (1993, 1994, 1995)
  KeyPia = 4;   // VP of Jegerlehner 2019
//  KeyPia = 5;   // VP of Teubner 2018
  KeyPib = 1;   // VP SHIFTED by the error value ON
  BHL_generator->SetKeyMod(KeyMod);      // steering param. in matrix element
  BHL_generator->SetKeyPia(KeyPia);      // steering param. in matrix element
  BHL_generator->SetKeyPib(KeyPib);      // steering param. in matrix element
  BHL_generator->ReModel(WtMain241);     // recalculating matrix element
  //====================== columns in the Table ===========================
  // 2.x or 4x with VP OFF, Z off !!!
  FillTrig92(hst_trig_VP0, WtMain100);   // 2.x VPoff, Z off, raws 1-5
  FillTrig95(hst_trig_VP0, WtMain200);   // 4.x VPoff, Z off, raws 11-17
  // 2.x or 4x with VPold
  FillTrig92(hst_trig_VP1, WtMain110);   // 2.x+VPold89, Z off, raws 1-5
  FillTrig95(hst_trig_VP1, WtMain230);   // 4.x+VPold,   Z off, raws 11-17
  // err of VPold in 2x or 4x
  FillTrig92(hst_trig_VP2, WtMain110-WtMain111);  // VPerrOld89 in 2.x
  FillTrig95(hst_trig_VP2, WtMain230-WtMain231);  // VPerrOld in 4.x
  // 4x with VPnew
  FillTrig92(hst_trig_VP3, WtMain240);   // 4.x+VPnew
  FillTrig95(hst_trig_VP3, WtMain240);   // 4.x+VPnew
  // err of VPnew in 4x
  FillTrig92(hst_trig_VP4, WtMain240-WtMain241);  // VPerrNew in 4.x
  FillTrig95(hst_trig_VP4, WtMain240-WtMain241);  // VPerrNew in 4.x
  // ***** Correction to LEPdata: VPold->VPnew (Z off), optionally 2.x->4.x *****
  FillTrig92(hst_trig_VP5, WtMain240-WtMain110);  // (4.x+VPnew)-(2.x+VPold89)
  FillTrig95(hst_trig_VP5, WtMain240-WtMain230);  // (4.x+VPnew)-(4.x+VPold)
//==================================================================
// ************ Study of Z contribution ************
  FillTrig92(hst_trig_ZZ0, WtMain230);   // 4.x+VPold, Z off, raws 1-5
  FillTrig95(hst_trig_ZZ0, WtMain230);   // 4.x+VPold, Z off, raws 11-17
  double WtMain1230, WtMain1230_12,WtMain1230_82, WtMain1230_92;
// 4.x+VPold, Z on
  KeyZet = 1;   // Z is ON
  BHL_generator->SetKeyZet(KeyZet);
  KeyMod = 2;   // 4.x QED m.e. model (1993, 1994, 1995)
  KeyPia = 3;   // VP of Burkhardt and Pietrzyk 1995 (Moriond, PLB96) MANDATORY for Z
  KeyPib = 0;   // VP shifted by the error value OFF
  BHL_generator->SetKeyMod(KeyMod);      // steering param. in matrix element
  BHL_generator->SetKeyPia(KeyPia);      // steering param. in matrix element
  BHL_generator->SetKeyPib(KeyPib);      // steering param. in matrix element
  BHL_generator->ReModel(WtMain1230);    // recalculating matrix element expon.
  WtMain1230_12 = BHL_generator->GetWtAlter(12); //  Z contrib. in standard expon. case
  WtMain1230_82 = BHL_generator->GetWtAlter(82); //  O(alpha^1), with VP. and Z self en.
  WtMain1230_92 = BHL_generator->GetWtAlter(92); //  O(alpha^1), BABAMC emulation?
  //[[[[[[[[[[[[
  //WtMain1230_82 = BHL_generator->GetWtAlter(81); //  O(alpha^1), with VP. and Z self en.
  //WtMain1230_92 = BHL_generator->GetWtAlter(91); //  O(alpha^1), BABAMC emulation?
  //]]]]]]]]]]]]
  FillTrig92(hst_trig_ZZ1, WtMain1230);      // 4.x+VPold, Z ON, raws 1-5
  FillTrig95(hst_trig_ZZ1, WtMain1230);      // 4.x+VPold, Z ON, raws 11-17
  //
  FillTrig92(hst_trig_ZZ2, WtMain1230_12);   // Z only EXPON, raws 1-5
  FillTrig95(hst_trig_ZZ2, WtMain1230_12);   // Z only EXPON, raws 11-17
  //
  FillTrig92(hst_trig_ZZ3, WtMain1230_82);   // Z only O(alf^1), raws 1-5
  FillTrig95(hst_trig_ZZ3, WtMain1230_82);   // Z only O(alf^1), raws 11-17
  //
  FillTrig92(hst_trig_ZZ4, WtMain1230_92);   // Z only O(alf^1), raws 1-5
  FillTrig95(hst_trig_ZZ4, WtMain1230_92);   // Z only O(alf^1), raws 11-17
}// WtCrude

} //End of Production

void TRobolProd::FillTrig92( TH1D   *hst_trig , double WtMain )
{
  int lWW,lNN,lNW,lWN;
 // BHLIMI 2.x
 // first generation lumi detectors 1990-92
 //
  f_TRigLib->TAlephLcal(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill(1-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill(1-1e-5,WtMain/2);
  // OPAL excepted, redone in 2000 using 4.x
  //f_TRigLib->TOpalFcal(lWW,lNN,lNW,lWN);
  //if(lNW) hst_trig->Fill(2-1e-5,WtMain/2);
  //if(lWN) hst_trig->Fill(2-1e-5,WtMain/2);
  //
  f_TRigLib->TL3bgo1(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill(3-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill(3-1e-5,WtMain/2);
  //
  f_TRigLib->TDelphiSat(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill(4-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill(4-1e-5,WtMain/2);
  //
  // SICAL in 92 but with BHLUMI 2.x
  f_TRigLib->TAleph93(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill( 5-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill( 5-1e-5,WtMain/2);
  //
}// FillTrig92


void TRobolProd::FillTrig95( TH1D   *hst_trig , double WtMain )
{
  int lWW,lNN,lNW,lWN;
  // Bhlumi 4.x
  // second generation lumi detectors 1093-95
  //
  f_TRigLib->TAleph93(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill( 11-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill( 11-1e-5,WtMain/2);
  //
  f_TRigLib->TAleph94(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill( 12-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill( 12-1e-5,WtMain/2);
  //
  f_TRigLib->TAleph95(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill( 13-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill( 13-1e-5,WtMain/2);
  //
  f_TRigLib->TOpal93(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill( 14-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill( 14-1e-5,WtMain/2);
  //
  f_TRigLib->TL3bgo2(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill( 15-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill( 15-1e-5,WtMain/2);
  //
  f_TRigLib->TDelphi94(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill( 16-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill( 16-1e-5,WtMain/2);
  //
  // DELPHI SAT with BHL 4.x in 1993
  f_TRigLib->TDelphiSat(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill(17-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill(17-1e-5,WtMain/2);
  //
  // first generation lumi detectors 1990-92
  // OPAL included because was fully updated to 4.x in 2000
  f_TRigLib->TOpalFcal(lWW,lNN,lNW,lWN);
  if(lNW) hst_trig->Fill( 2-1e-5,WtMain/2);
  if(lWN) hst_trig->Fill( 2-1e-5,WtMain/2);
  //
}// FillTrig95


///////////////////////////////////////////////////////////////////////////////
void TRobolProd::Finalize()
{
//   Finalize MC  run, final printouts, cleaning etc., xcheck of normalization
//   Plotting histograms is done independently using root file
  TMCgenBHL6 *BHL_generator = (TMCgenBHL6*)f_MCgen;
//
  double XsNormPb, XsErroPb;
  BHL_generator->Finalize();
  XsNormPb =BHL_generator->m_XsNormPb;
  XsErroPb =BHL_generator->m_XsErroPb;
  cout << " TRobolProd: XsNormPb [pb] = "<<  XsNormPb << "  +-  "<< XsErroPb <<endl;
//  double xSecPb,xErrPb,xSecNb;
//  BHL_generator->GetXsecMC( xSecPb, xErrPb);
//  xSecNb=xSecPb/1000;
//  cout << " KKMC: xSecPb   [pb] = "<<  xSecPb << "  +-  "<< xErrPb <<endl;
//  cout << " KKMC: xSecNb   [nb] = "<<  xSecNb << "  +-  "<< xErrPb/1000 <<endl;
}



//______________________________________________________________________________
TH1D *TRobolProd::TH1D_UP(const char* name, const char* title,
                       int nbins, double xmin, double xmax)
{
  TH1D *h;
  if (f_isNewRun) {
    h = new TH1D(name,title,nbins,xmin,xmax);
    h->Sumw2();
  } else {
    h = (TH1D*)f_HstFile->Get(name);
  }
  return h;
}

//______________________________________________________________________________
TH2D *TRobolProd::TH2D_UP(const char* name, const char* title,
                       int nbinx, double xmin, double xmax,
                       int nbiny, double ymin, double ymax)
{
  TH2D *h;
  if (f_isNewRun) {
    h = new TH2D(name,title,nbinx,xmin,xmax,nbiny,ymin,ymax);
    h->Sumw2();
  } else {
    h = (TH2D*)f_HstFile->Get(name);
  }
  return h;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                             UTILITIES                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


void TRobolProd::MomPrint( TLorentzVector &Vect){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (with endline)
  for ( int k=0; k < 4 ; k++ )   cout << SP18 << Vect[k];
  cout<<endl;
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           End of Class ROBOL                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
