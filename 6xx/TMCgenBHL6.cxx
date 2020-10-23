//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   TMCgenBHL6                                          //
//                                                                          //
//            Interface (wrapper)  to MC event generator BHLUMI             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include "TMCgenBHL6.h"

#define SP21 setw(21)<<setprecision(13)
#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(5)

ClassImp(TMCgenBHL6);

TMCgenBHL6::TMCgenBHL6()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "----> TMCgenBHL6 Default Constructor (for ROOT only) "<<endl;
  cout<< "====> TMCgen::TMCgen DEFAULT Constructor (for ROOT only) in MCdev"<<endl;
  f_Out          = NULL;
  f_GenFile      = NULL;
  f_HstFile      = NULL;
  f_RNgen        = NULL;
  f_TMCgen_NORMA = NULL;
}

///_____________________________________________________________
TMCgenBHL6::TMCgenBHL6(const char* Name)
{
//! all defaults defined here can be changed by the user
//! before calling TMCgenBHL6::Initialize
  sprintf(f_Name,"%s",Name);         // Class name
  sprintf(f_Date,"%s","  Release date:  yyyy.mm.dd   "); // Release date
  f_Version  = 1.00;                                      // Release version
///
  f_RNgen        = NULL;
  f_Out          = NULL;
  f_GenFile      = NULL;
  f_HstFile      = NULL;
  f_TMCgen_NORMA = NULL;
  ///
  f_IsInitialized = 0; /// prevents Initialization when reset to 1
  f_NevGen  =   0;
///////////////////////////////////////////////////
// Physics
  m_NevTot = 0;
  m_EvenCounter = 0;
// Seeds for Pseumar
  m_ijkl_new  = 54217137;
  m_ntot_new  = 0;
  m_ntot2_new = 0;
//---------------------------------------------------------------------
// ---------------- BHLUMI4 input params ------------------------------
// NPAR( 1) = KeyOpt =1000*KeyGen +100*KeyRem +10*KeyWgt +KeyRnd =3121
  m_KeyGen =3;
  m_KeyRem =1;
  m_KeyWgt =2;
  m_KeyRnd =1;
// NPAR( 2) = KeyRad =1000*KeyZet +100*KeyUpd +10*KeyMod +KeyPia = 20
  m_KeyZet =0;
  m_KeyUpd =0;
  m_KeyMod =2;
  m_KeyPia =0;
//  BHLUMI4 input params
  m_CMSENE  =   92.3e0;    //  CMSENE
  m_angmin  =     .022;    //  theta_min [rad]   genaration
  m_angmax  =     .082;    //  theta_max [rad]   generation
  m_EPSCM   =     1e-3;    //  XK0     eps_CMS   generation
  m_VMAXG   = 0.9999e0;    //  VMAXG   v_max     generation
//
  cout<< "----> TMCgenBHL6::TMCgenFOAM USER Constructor "<<endl;
}///TMCgenBHL6

///______________________________________________________________________________________
TMCgenBHL6::~TMCgenBHL6()
{
  //!Explicit destructor
  cout<< "----> TMCgenBHL6::TMCgenBHL6 !!!! DESTRUCTOR !!!! "<<endl;
}///destructor


///______________________________________________________________________________________
void TMCgenBHL6::Initialize(TRandom *RNgen, ofstream *OutFile, TH1D* h_NORMA)
{
  cout<< "----> TMCgenBHL6::Initialize, Entering "<<endl;

  f_Out     = OutFile;
  f_RNgen   = RNgen;        // connect to external RN generator
  f_TMCgen_NORMA = h_NORMA; // special external normalization
  ///----------------------------------------
  /// physics params, redefined defaults
  ///----------------------------------------
  ///=================================================================
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"====== TMCgenBHL6::Initialize     ======");
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,f_Name);
  BX1F(*f_Out,"  Version",f_Version,  f_Date);
  BXTXT(*f_Out,"============== INPUT ===================");
  BXCLO(*f_Out);

  m_NevTot = 0;
  m_EvenCounter = 0;
  m_jmax  = 10000;

// Initialisation input data before generation
  cout << "*******************************" << endl;
  cout << "**   TMCgenBHL6   Initialize  **" << endl;
  cout << "*******************************" << endl;

// Input data BHLUMI4 style, with default values possibly overwritten

  int KeyOpt =1000*m_KeyGen +100*m_KeyRem +10*m_KeyWgt +m_KeyRnd;
  int KeyRad =1000*m_KeyZet +100*m_KeyUpd +10*m_KeyMod +m_KeyPia;
  double TRMIN = sqr(m_CMSENE)*(1-cos(m_angmin))/2;
  double TRMAX = sqr(m_CMSENE)*(1-cos(m_angmax))/2;

  m_npar[1 -1]=KeyOpt; // npar in c++ numbering (from zero)
  m_npar[2 -1]=KeyRad;
  m_ypar[1 -1]=m_CMSENE; // ypar in c++ numbering
  m_ypar[2 -1]=TRMIN;
  m_ypar[3 -1]=TRMAX;
  m_ypar[4 -1]=m_EPSCM;

  cout<<" ------------------------------------"<<endl;
  cout<<" MCgenBHL::Initialize control printout "<<endl;
  cout<<" CMSENE= "<<m_CMSENE<<endl;
  cout<<" KEYOPT= "<<KeyOpt<<endl;
  cout<<" KEYRAD= "<<KeyRad<<endl;
  cout<<" angmin= "<<m_angmin<<endl;
  cout<<" angmax= "<<m_angmax<<endl;
  cout<<" TRMIN=  "<<TRMIN<<endl;
  cout<<" TRMAX=  "<<TRMAX<<endl;
  cout<<" EPSCM=  "<<m_EPSCM<<endl;
  cout<<" VMAXG=  "<<m_VMAXG<<endl;

  cout<<"%%%%% TMCgenBHL6::Initialize: CMSENE="<<m_CMSENE<<endl;

  //=============================================================
  //   opening separate disk file for fortran part of BHLUMI code
  m_out = 6;
  const char *output_file = "./bhl_f77.output";
  int sl2 = strlen(output_file);
  bhl_fort_open_(m_out,output_file,sl2);
  SetInOut( 5, m_out );  // input and output for fortran
  //***********************************************************//
  m_ijkl_new = RNgen->GetSeed();
  cout<<" TMCgen::Initialize: r.n. generator seed  = "<< m_ijkl_new <<endl;
  marini_(m_ijkl_new, m_ntot_new, m_ntot2_new );

        cout<<" TMCgen::Initialize: RANMAR initialized to "<< m_ijkl_new <<endl;
  (*OutFile)<<" TMCgen::Initialize: RANMAR ijkl, ntot, ntot2 ="
		    << m_ijkl_new <<" "<<m_ntot_new<<"  "<< m_ntot2_new <<endl;

  // Set storage for glibk
  glimit_(50000);
  /////////////////////////
  //  CALL BHLUMI(  -1,XPAR,NPAR)
  bhlum4_(-1, m_ypar, m_npar);
  /////////////////////////
  //*******************************************//

  int NevPrim;
  GetPrimaNorma(NevPrim, m_XsNormPb);     //   Primary Xsection for normalization NANOBARNS ????
  cout<<"/////// TMCgen::Initialize: m_XsNormPb="<<m_XsNormPb<<endl;

  for(int j=0;   j<100; j++)  h_NORMA->SetBinContent(j,     m_ypar[j] );    // ypar encoded
  for(int j=100; j<200; j++)  h_NORMA->SetBinContent(j+100, m_npar[j] );    // npar encoded

  cout<<" TMCgenBHL6::Initialize:  ypar() and npar() filled into h_NORMA  "<<endl;

  /////////////////////////////////////////////////////////
  if(f_IsInitialized != 0)
	  cout<< "----> TMCgenBHL6::Initialize, already initialized "<<endl;

}// Initialize


///______________________________________________________________________________________
void TMCgenBHL6::ReInitialize(ofstream *OutFile, TH1D* h_NORMA)
// this is for reinitializing BHLUMI during analysis of MC results
// random number generator initialization not needed
{
  cout<< "----> TMCgenBHL6::Initialize, Entering "<<endl;

  f_Out     = OutFile;
  f_TMCgen_NORMA = h_NORMA; // special external normalization histogram
  ///----------------------------------------
  /// physics params, redefined defaults
  ///----------------------------------------
  ///=================================================================
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"====== TMCgenBHL6::ReInitialize   ======");
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,f_Name);
  BX1F(*f_Out,"  Version",f_Version,  f_Date);
  BXTXT(*f_Out,"============== INPUT ===================");
  BXCLO(*f_Out);
  cout << "*******************************" << endl;
  cout << "**   TMCgenBHL6   ReInitialize  **" << endl;
  cout << "*******************************" << endl;

// Input data BHLUMI4 style, with default values possibly overwritten

  cout<<" ------------------------------------"<<endl;
  cout<<" MCgenBHL::ReInitialize control printout "<<endl;
  cout<<" CMSENE= "<<m_CMSENE<<endl;
  cout<<" angmin= "<<m_angmin<<endl;
  cout<<" angmax= "<<m_angmax<<endl;
  cout<<" TRMIN=  "<<m_TRMIN<<endl;
  cout<<" TRMAX=  "<<m_TRMAX<<endl;
  cout<<" EPSCM=  "<<m_EPSCM<<endl;
  cout<<" VMAXG=  "<<m_VMAXG<<endl;

  //=============================================================
  //   opening separate disk file for fortran part of BHLUMI code
  m_out = 6;
  const char *output_file = "./bhl_f77.output";
  int sl2 = strlen(output_file);
  bhl_fort_open_(m_out,output_file,sl2);
  SetInOut( 5, m_out );  // input and output for fortran
  //***********************************************************//

  // Set storage for glibk
  glimit_(50000);
  /////////////////////////
  //  CALL BHLUMI(  -1,XPAR,NPAR)
  bhlum4_(-1, m_ypar, m_npar);
  /////////////////////////
  //*******************************************//

  int NevPrim;
  GetPrimaNorma(NevPrim, m_XsNormPb);     //   Primary Xsection for normalization NANOBARNS ????
  cout<<"/////// TMCgen::Initialize: m_XsNormPb="<<m_XsNormPb<<endl;

  for(int j=0;   j<100; j++)  m_ypar[j] = h_NORMA->GetBinContent(j  );      // ypar recovered
  for(int j=100; j<200; j++)  m_npar[j] = h_NORMA->GetBinContent(j+100);    // npar recovered

  cout<<" TMCgenBHL6::ReInitialize:  ypar() and npar() recoveredh_NORMA  "<<endl;

  /////////////////////////////////////////////////////////
 // if(f_IsInitialized == 0)
//	  { cout<< "----> TMCgenBHL6::ReInitialize not initialized !!! "<<endl; exit(10);}

}// ReInitialize


///////////////////////////////////////////////////////////////////////////////
void TMCgenBHL6::Generate()
{
  f_NevGen++;
//****************************************************
// ******* invoke BHLUMI4.x generator in f77 *********
  bhlum4_( 0, m_ypar, m_npar);
//****************************************************
// Filling in normalization histogram using data from f77 part
  double XsPrimPb; int NevPrim;
  GetPrimaNorma( NevPrim, XsPrimPb);  //   Primary Xsection for normalization NANOBARNS
  f_TMCgen_NORMA->SetBinContent(0,XsPrimPb*NevPrim);  //
  f_TMCgen_NORMA->SetEntries(NevPrim);

//  cout<<"//////// TMCgenBHL6::Generate: NevPrim, XsPrimPb= "<< NevPrim<<"  "<< XsPrimPb <<endl;

}
///////////////////////////////////////////////////////////////////////////////
void TMCgenBHL6::Finalize()
{

  cout << "*****************************" << endl;
  cout << "**   TMCgenBHL6   Finalize  **" << endl;
  cout << "*****************************" << endl;
  double nevtot = f_TMCgen_NORMA->GetBinContent(2);
  ///
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"====== TMCgenBHL6::Finalize       ======");
  BXTXT(*f_Out,"========================================");
  BX1I(*f_Out,"   nevtot", nevtot, " number of generated events                      =");
  //
  bhlum4_( 2, m_ypar, m_npar);
  //
  //
  bhl_fort_close_(m_out);
}//Finalize

///////////////////////////////////////////////////////////////////////////////
void TMCgenBHL6::GetPrimaNorma(int &NevPrim, double &XsPrim)
{
// get normalization elements NANOBARNS
  int NevPrim1;
  bhl_get_primnorm_( NevPrim1, XsPrim);
  NevPrim=NevPrim1;
}

///////////////////////////////////////////////////////////////////////////////
void TMCgenBHL6::GetPhoton1(const int iphot, TLorentzVector &phot)
{
// get one photon 4-vector from KKMC
  double p1[4];
  int iphot1=iphot;
  bhl_getphoton1_(iphot1, p1);
  phot.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
}

///////////////////////////////////////////////////////////////////////////////
void TMCgenBHL6::GetWt( double &WtMain,  double &WtCrude)
{
// get MC weight of event
  bhl_getwt_(WtMain, WtCrude);
}

///////////////////////////////////////////////////////////////////////////////
void TMCgenBHL6::GetBeams( TLorentzVector &B1,  TLorentzVector &B2)
{
// get 4-momenta of beams
  double p1[4];
  double p2[4];
  bhl_getbeams_(p1,p2);
  B1.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
  B2.SetPxPyPzE(p2[0],p2[1],p2[2],p2[3]);
}

///////////////////////////////////////////////////////////////////////////////
void TMCgenBHL6::GetFermions( TLorentzVector &B1,  TLorentzVector &B2)
{
// get 4-momenta of beams
  double p1[4];
  double p2[4];
  bhl_getfermions_(p1,p2);
  B1.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
  B2.SetPxPyPzE(p2[0],p2[1],p2[2],p2[3]);
}

///////////////////////////////////////////////////////////////////////////////
void TMCgenBHL6::GetNphot( int &Nphot)
{
// get photon multiplicity from KKMC
  bhl_getnphot_( Nphot);
}

///////////////////////////////////////////////////////////////////////////////
double TMCgenBHL6::GetWtAlter(const int id)
{
  double WtAlter;
  int id1 =id;
  bhl_getwtalter_( id1, WtAlter);
  return WtAlter;
}


///////////////////////////////////////////////////////////////////////////////
void TMCgenBHL6::Print1()
{
// print event using BHLumi4
  dumps_(m_out);
}

void TMCgenBHL6::SetInOut( const int &ninp, const int &nout )
{
// Set input/output for fortran
   bhl_setinout_( ninp, nout );
}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//              End of  Class   TMCgenBHL6                                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
