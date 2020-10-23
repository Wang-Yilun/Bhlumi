//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   TRig                                               //
//                                                                          //
//   Interface (wrapper) event selections (triggers) of Lumi LEP detectors  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include "TRig.h"

#define SP21 setw(21)<<setprecision(13)
#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(5)

ClassImp(TRig);

TRig::TRig()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "----> TRig Default Constructor (for ROOT only) "<<endl;
  cout<< "====> TMCgen::TMCgen DEFAULT Constructor (for ROOT only) in MCdev"<<endl;
}

///_____________________________________________________________
TRig::TRig(const char* Name)
{
//! all defaults defined here can be changed by the user
//! before calling TRig::Initialize
  sprintf(f_Name,"%s",Name);         // Class name
  sprintf(f_Date,"%s","  Release date:  yyyy.mm.dd   "); // Release date
  f_Version  = 1.00;                                      // Release version
///

  m_CMSENE  =   0;    //  CMSENE

//
  cout<< "----> TRig::TMCgenFOAM USER Constructor "<<endl;
}///TRig

///______________________________________________________________________________________
TRig::~TRig()
{
  //!Explicit destructor
  cout<< "----> TRig::TRig !!!! DESTRUCTOR !!!! "<<endl;
}///destructor


///______________________________________________________________________________________
void TRig::Initialize()
{
  cout<< ">>>>>> TRig::Initialize, Entering "<<endl;
}// Initialize

///______________________________________________________________________________________
void TRig::PrintAng( const char *name, double th1w, double th2w,  double th1n, double th2n)
{
  cout<< ">>>>>> "<< name << SP15<< th1w<< SP15<<th2w<< SP15<< th1n<< SP15 <<th2n << endl;
}// Initialize


