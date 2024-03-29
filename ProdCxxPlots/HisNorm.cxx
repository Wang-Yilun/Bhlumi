/////////////////////////////////////////////////////////////////////
// Collection of programs for renormalizing histograms
// and projecting 2D scatergrams into 1D histograms
// Also tools for calculating AFB(vmax)
/////////////////////////////////////////////////////////////////////
#include "HisNorm.h"

double sqr( const Double_t x ){ return x*x;};


// This works for 1-dim histograms
void HisNorm0( long   Nevt, double Xsav, TH1 *Hst){
// normalize histogram according to Xsav
  cout<<"HisNorm0: Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<endl;
  //
  int      nbX  = Hst->GetNbinsX();
  //cout<<"nbt = "<<nbt<<endl;
  Double_t Xmax = Hst->GetXaxis()->GetXmax();
  Double_t Xmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbX*Xsav/(Xmax-Xmin)/Nevt;
  cout<<"HisNorm0: Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}

/////////////////////////////////////////////////////////////////////
// This works for 1-dim histograms
void HisNorm1(TH1D *NorHst, TH1 *Hst){
  // normalize histogram in nanobarns
  Long_t   Nevt = NorHst->GetEntries();
  Double_t Xsav = NorHst->GetBinContent(0)/Nevt; // NANOBARNS
  cout<<"HisNorm1: Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<endl;
  //
  int      nbX  = Hst->GetNbinsX();
  //cout<<"nbt = "<<nbt<<endl;
  Double_t Xmax = Hst->GetXaxis()->GetXmax();
  Double_t Xmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbX*Xsav/(Xmax-Xmin)/Nevt;
  cout<<"HisNorm1: Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}


/////////////////////////////////////////////////////////////////////
void HisNorm1M(TH1D *NorHst, TH1D *Hst, Float_t msize, Int_t mcolor, Int_t mark){
  // normalize histogram in nanobarns and sets plotting params
  Long_t   Nevt = NorHst->GetEntries();
  Double_t Xsav = NorHst->GetBinContent(0)/Nevt/Nevt;
  //
  int      nbX  = Hst->GetNbinsX();
  Double_t Xmax = Hst->GetXaxis()->GetXmax();
  Double_t Xmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbX*Xsav/(Xmax-Xmin);
  //
  Hst->Scale(Fact);
  //
  Hst->SetMarkerStyle( mark);
  Hst->SetMarkerColor( mcolor);
  Hst->SetLineColor(   mcolor);
  Hst->SetMarkerSize(  msize);
}

/////////////////////////////////////////////////////////////////////
// This works for 2-dim histograms
void HisNorm2(TH1D *NorHst, TH2 *Hst){
  // normalize histogram in nanobarns
  Long_t   Nevt = NorHst->GetEntries();
  Double_t Xsav = NorHst->GetBinContent(0)/Nevt; // NANOBARNS
  cout<<"HisNorm2: Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<endl;
  //
  int      nbX  = Hst->GetNbinsX();
  int      nbY  = Hst->GetNbinsY();
  //cout<<"nbt = "<<nbt<<endl;
  Double_t Xmax = Hst->GetXaxis()->GetXmax();
  Double_t Xmin = Hst->GetXaxis()->GetXmin();
  Double_t Ymax = Hst->GetYaxis()->GetXmax();
  Double_t Ymin = Hst->GetYaxis()->GetXmin();
  Double_t Fact = Xsav/Nevt;
  Fact *= nbX/(Xmax-Xmin);
  Fact *= nbY/(Ymax-Ymin);
  cout<<"HisNorm2: Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}

TH1D *HstConst(TString title, TH1D *hst1, double yy)
{
  // makes histogram filled in with constant
  cout<<"Entering HstConst ";
  cout<< hst1->GetName() <<endl;
  int      nbX  = hst1->GetNbinsX();
  Double_t Xmax = hst1->GetXaxis()->GetXmax();
  Double_t Xmin = hst1->GetXaxis()->GetXmin();
  //
  TH1D *hOne = (TH1D*)hst1->Clone("hOne"); // allocate hOne
  hOne->Reset();
  for(int iv=0; iv <= nbX; iv++){
    hOne->SetBinContent(iv, yy);
    hOne->SetBinError(  iv, 0.0);
   }
  hOne->SetName(title);
  return hOne;
}//HstOne


///////////////////////////////////////////////////////////////////////////////////
void ProjX1(TH2D *Scat, TH1D *&HstProjX)
{
  // Simple Projection onto X axis taking into account errors!
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  HstProjX = (TH1D*)Scat->ProjectionX("HstNew",1,nbX,"e"); // option "e" doesnt work
  HstProjX->Reset();
  double sum,sum2, dy;
  dy= (Ymax-Ymin)/nbY; // integration over Y
  for(int ix=0; ix <= nbX+1; ix++){
    sum=0.0; sum2=0.0;
    for(int iy=0; iy <= nbY+1; iy++){
      sum  += Scat->GetBinContent(ix,iy);
      sum2 += sqr(Scat->GetBinError(ix,iy));
    }
    HstProjX->SetBinContent(ix, dy*sum);
    HstProjX->SetBinError(  ix, dy*sqrt(sum2));
  }
}// ProjX1

///////////////////////////////////////////////////////////////////////////////////
void ProjY1(TH2D *Scat, TH1D *&HstProjY)
{
  // Simple Projection onto Y axis taking into account errors!
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  HstProjY = (TH1D*)Scat->ProjectionY("HstNew",1,nbX,"e"); // option "e" doesnt work
  HstProjY->Reset();
  double sum,sum2, dx;
  dx= (Xmax-Xmin)/nbX; // integration over X
  for(int iy=0; iy <= nbY+1; iy++){
    sum=0.0; sum2=0.0;
    for(int ix=0; ix <= nbX+1; ix++){
      sum  += Scat->GetBinContent(ix,iy);
      sum2 += sqr(Scat->GetBinError(ix,iy));
    }
    HstProjY->SetBinContent(iy, dx*sum);
    HstProjY->SetBinError(  iy, dx*sqrt(sum2));
  }
}// ProjY1

TH1D *HstProjV(TString title, TH2D *&Scat, int NbMax)
{ // makes cumulative distribution in sigma(Vmax)
  // integrating over cos(theta) bins up to NbMax.
  cout<<"Entering HstProjV for TH2D ";
  cout<< Scat->GetName() <<endl;
  //  Projection onto v axis, suming over Y=cos(theta) up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  TH1D *hxTot = (TH1D*)Scat->ProjectionX(title,1,nbX,"e");
  hxTot->Reset();
  //
  double forw,forw2, back, back2, dx,dy;
  double Forw,Forw2, Back, Back2;
  dx= (Xmax-Xmin)/nbX; // integration over X
  dy= (Ymax-Ymin)/nbY; // integration over Y
  Forw=0.0; Forw2=0.0;
  Back=0.0; Back2=0.0;
  int nbYhalf = nbY/2;
  int nbY2;
  if( (NbMax>0) && (NbMax<nbYhalf) )
    nbY2 = NbMax;
  else
    nbY2 = nbYhalf;
  for(int ix=0; ix <= nbX+1; ix++){
    forw=0.0; forw2=0.0;
    back=0.0; back2=0.0;
    // loop over cos(theta) bins
    for(int iy=1; iy <= nbY2; iy++){
      forw  += Scat->GetBinContent(  ix, nbYhalf+iy);
      forw2 += sqr(Scat->GetBinError(ix, nbYhalf+iy));
      back  += Scat->GetBinContent(  ix, nbYhalf-iy+1);
      back2 += sqr(Scat->GetBinError(ix, nbYhalf-iy+1));
    }// iy
    Forw  += forw;  Forw2 += forw2;
    Back  += back;  Back2 += back2;
    hxTot->SetBinContent(ix, dx*dy*    (Forw +Back));
    hxTot->SetBinError(  ix, dx*dy*sqrt(Forw2+Back2));
  }
  return hxTot;
}// HstProjV


TH1D *HstProjF(TString title, TH2D *&Scat, int NbMax)
{ // makes cumulative distribution in sigma(Vmax)
  // integrating over cos(theta) bins up to NbMax.
  // FORWARD hemisphere only!!!
  cout<<"Entering HstProjV for TH2D ";
  cout<< Scat->GetName() <<endl;
  //  Projection onto v axis, suming over Y=cos(theta) up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  TH1D *hxForw = (TH1D*)Scat->ProjectionX(title,1,nbX,"e");
  hxForw->Reset();
  //
  double forw,forw2, dx,dy;
  double Forw,Forw2;
  dx= (Xmax-Xmin)/nbX; // integration over X
  dy= (Ymax-Ymin)/nbY; // integration over Y
  Forw=0.0; Forw2=0.0;
  int nbYhalf = nbY/2;
  int nbY2;
  if( (NbMax>0) && (NbMax<nbYhalf) )
    nbY2 = NbMax;
  else
    nbY2 = nbYhalf;
  for(int ix=0; ix <= nbX+1; ix++){
    forw=0.0; forw2=0.0;
    // loop over cos(theta) bins
    for(int iy=1; iy <= nbY2; iy++){
      forw  += Scat->GetBinContent(  ix, nbYhalf+iy);
      forw2 += sqr(Scat->GetBinError(ix, nbYhalf+iy));
    }// iy
    Forw  += forw;  Forw2 += forw2;
    hxForw->SetBinContent(ix, dx*dy*    (Forw));
    hxForw->SetBinError(  ix, dx*dy*sqrt(Forw2));
  }
  return hxForw;
}// HstProjV


TH1D *HstProjAv(TString title, TH2D *&Scat, int NbMax)
{ // makes AFB(V), integrating over cos(theta) bins up to NbMax.
  cout<<"Entering HstProjA for  ";
  cout<< Scat->GetName() <<endl;
  //  Projection onto v axis, suming over Y=cos(theta) up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  double deno;
  //
  TH1D *hxAfb = (TH1D*)Scat->ProjectionX(title,1,nbX,"e");
  hxAfb->Reset();
  //
  double forw,forw2, back, back2;
  double Forw,Forw2, Back, Back2;
//  dx= (Xmax-Xmin)/nbX; // integration over X
//  dy= (Ymax-Ymin)/nbY; // integration over Y
  Forw=0.0; Forw2=0.0;
  Back=0.0; Back2=0.0;
  int nbYhalf = nbY/2;
  int nbY2;
  if( (NbMax>0) && (NbMax<nbYhalf) )
    nbY2 = NbMax;
  else
    nbY2 = nbYhalf;
  for(int ix=0; ix <= nbX+1; ix++){
    forw=0.0; forw2=0.0;
    back=0.0; back2=0.0;
    // loop over cos(theta) bins
    for(int iy=1; iy <= nbY2; iy++){
      forw  += Scat->GetBinContent(  ix, nbYhalf+iy);
      forw2 += sqr(Scat->GetBinError(ix, nbYhalf+iy));
      back  += Scat->GetBinContent(  ix, nbYhalf-iy+1);
      back2 += sqr(Scat->GetBinError(ix, nbYhalf-iy+1));
    }// iy
//    Forw  += forw;  Forw2 += forw2;
//    Back  += back;  Back2 += back2;
    Forw  = forw;  Forw2 = forw2;  // no cumulation in v
    Back  = back;  Back2 = back2;  // no cumulation in v
    deno = Forw +Back;
    if( deno == 0.0 ) {
        hxAfb->SetBinContent(ix, 0 );
        hxAfb->SetBinError(  ix, 0 );
    }else{
       hxAfb->SetBinContent(ix,      (Forw -Back)/(Forw +Back) );
       hxAfb->SetBinError(  ix, sqrt(Forw2+Back2)/(Forw +Back) );
    }
    cout<<"  afb="<< (Forw -Back)/(Forw +Back);
  }//ix
  return hxAfb;
}// HstProjAv



TH1D *HstProjA(TString title, TH2D *&Scat, int NbMax)
{ // makes AFB(Vmax), integrating over cos(theta) bins up to NbMax.
  cout<<"Entering HstProjA for  ";
  cout<< Scat->GetName() <<endl;
  //  Projection onto v axis, suming over Y=cos(theta) up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  TH1D *hxAfb = (TH1D*)Scat->ProjectionX(title,1,nbX,"e");
  hxAfb->Reset();
  //
  double forw,forw2, back, back2, dx,dy;
  double Forw,Forw2, Back, Back2;
  dx= (Xmax-Xmin)/nbX; // integration over X
  dy= (Ymax-Ymin)/nbY; // integration over Y
  Forw=0.0; Forw2=0.0;
  Back=0.0; Back2=0.0;
  int nbYhalf = nbY/2;
  int nbY2;
  if( (NbMax>0) && (NbMax<nbYhalf) )
    nbY2 = NbMax;
  else
    nbY2 = nbYhalf;
  for(int ix=0; ix <= nbX+1; ix++){
    forw=0.0; forw2=0.0;
    back=0.0; back2=0.0;
    // loop over cos(theta) bins
    for(int iy=1; iy <= nbY2; iy++){
      forw  += Scat->GetBinContent(  ix, nbYhalf+iy);
      forw2 += sqr(Scat->GetBinError(ix, nbYhalf+iy));
      back  += Scat->GetBinContent(  ix, nbYhalf-iy+1);
      back2 += sqr(Scat->GetBinError(ix, nbYhalf-iy+1));
    }// iy
    Forw  += forw;  Forw2 += forw2;
    Back  += back;  Back2 += back2;
    hxAfb->SetBinContent(ix,      (Forw -Back)/(Forw +Back) );
    hxAfb->SetBinError(  ix, sqrt(Forw2+Back2)/(Forw +Back) );
  }//ix
  return hxAfb;
}// HstProjA



///////////////////////////////////////////////////////////////////////////////////
void ProjV(TH2D *Scat, TH1D *&hxTot, TH1D *&hxAfb, int NbMax)
{
  //  Projection onto v axis, suming over Y=cos(theta) up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  hxTot = (TH1D*)Scat->ProjectionX("HstNew",1,nbX,"e");
  hxTot->Reset();
  hxAfb = (TH1D*)hxTot->Clone("HstNew2");
  //
  double forw,forw2, back, back2, dx,dy;
  double Forw,Forw2, Back, Back2;
  dx= (Xmax-Xmin)/nbX; // integration over X
  dy= (Ymax-Ymin)/nbY; // integration over Y
  Forw=0.0; Forw2=0.0;
  Back=0.0; Back2=0.0;
  int nbYhalf = nbY/2;
  int nbY2;
  if( (NbMax>0) && (NbMax<nbYhalf) )
    nbY2 = NbMax;
  else
    nbY2 = nbYhalf;
  for(int ix=0; ix <= nbX+1; ix++){
    forw=0.0; forw2=0.0;
    back=0.0; back2=0.0;
    // loop over cos(theta) bins
    for(int iy=1; iy <= nbY2; iy++){
      forw  += Scat->GetBinContent(  ix, nbYhalf+iy);
      forw2 += sqr(Scat->GetBinError(ix, nbYhalf+iy));
      back  += Scat->GetBinContent(  ix, nbYhalf-iy+1);
      back2 += sqr(Scat->GetBinError(ix, nbYhalf-iy+1));
    }// iy
    Forw  += forw;  Forw2 += forw2;
    Back  += back;  Back2 += back2;
    //hxTot->SetBinContent(ix, dy*sum);
    //hxTot->SetBinError(  ix, dy*sqrt(sum2));
    hxTot->SetBinContent(ix, dx*dy*    (Forw +Back));
    hxTot->SetBinError(  ix, dx*dy*sqrt(Forw2+Back2));
    hxAfb->SetBinContent(ix, dx*dy*    (Forw -Back));
    hxAfb->SetBinError(  ix, dx*dy*sqrt(Forw2+Back2));
  }
  hxAfb->Divide(hxTot);
}// ProjV

///////////////////////////////////////////////////////////////////////////////////
void ProjC(TH2D *Scat, TH1D *&hTot, TH1D *&hAsy, int NbMax)
{
  // Projection onto c=cos(theta) axis, suming over v up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  hTot = (TH1D*)Scat->ProjectionY("HstNew",1,nbX,"e"); // option "e" doesnt work
  hTot->Reset();
  hAsy = (TH1D*)hTot->Clone("HstNew1");
  TH1D *hSym = (TH1D*)hTot->Clone("HstNew2"); // local histo
  //
  double forw,forw2;
  double back,back2;
  double dx= (Xmax-Xmin)/nbX; // integration over X
  int ivMax;
  //
  if( (NbMax>0) && (NbMax<nbX) )
    ivMax = NbMax;
  else
    ivMax = nbX;
  int nbYhalf = nbY/2;
  for(int ic=1; ic <= nbYhalf; ic++){
    int iforw = nbYhalf+ic;
    int iback = nbYhalf-ic+1;
    forw=0.0;  forw2=0.0;
    back=0.0;  back2=0.0;
    for(int iv=0; iv <= ivMax; iv++){
      forw  += Scat->GetBinContent(  iv,iforw);
      forw2 += sqr(Scat->GetBinError(iv,iforw));
      back  += Scat->GetBinContent(  iv,iback);
      back2 += sqr(Scat->GetBinError(iv,iback));
    }
    // total
    hTot->SetBinContent(iforw,     dx*(forw));
    hTot->SetBinError(  iforw, dx*sqrt(forw2));
    hTot->SetBinContent(iback,     dx*(back));
    hTot->SetBinError(  iback, dx*sqrt(back2));
    // asymetric
    hAsy->SetBinContent(iforw,     dx*(forw-back));
    hAsy->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hAsy->SetBinContent(iback,    -dx*(forw-back));
    hAsy->SetBinError(  iback, dx*sqrt(forw2+back2));
    // symetric
    hSym->SetBinContent(iforw,     dx*(forw+back));
    hSym->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hSym->SetBinContent(iback,     dx*(forw+back));
    hSym->SetBinError(  iback, dx*sqrt(forw2+back2));
  }
  hAsy->Divide(hSym);
  hSym->Delete();
}// ProjC



///////////////////////////////////////////////////////////////////////////////////
TH1D *HstProjCA(TString title, TH2D *Scat, int NbMax)
{
  // Projection onto c=cos(theta) axis, suming over v up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  TH1D *hAsy = (TH1D*)Scat->ProjectionY(title,1,nbX,"e"); // option "e" doesnt work
  hAsy->Reset();
  TH1D *hSym = (TH1D*)hAsy->Clone("HstNew2"); // local histo
  //
  double forw,forw2;
  double back,back2;
  double dx= (Xmax-Xmin)/nbX; // integration over X
  int ivMax;
  //
  if( (NbMax>0) && (NbMax<nbX) )
    ivMax = NbMax;
  else
    ivMax = nbX;
  int nbYhalf = nbY/2;
  for(int ic=1; ic <= nbYhalf; ic++){
    int iforw = nbYhalf+ic;
    int iback = nbYhalf-ic+1;
    forw=0.0;  forw2=0.0;
    back=0.0;  back2=0.0;
    for(int iv=0; iv <= ivMax; iv++){
      forw  += Scat->GetBinContent(  iv,iforw);
      forw2 += sqr(Scat->GetBinError(iv,iforw));
      back  += Scat->GetBinContent(  iv,iback);
      back2 += sqr(Scat->GetBinError(iv,iback));
    }
   // asymetric
    hAsy->SetBinContent(iforw,     dx*(forw-back));
    hAsy->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hAsy->SetBinContent(iback,    -dx*(forw-back));
    hAsy->SetBinError(  iback, dx*sqrt(forw2+back2));
    // symetric
    hSym->SetBinContent(iforw,     dx*(forw+back));
    hSym->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hSym->SetBinContent(iback,     dx*(forw+back));
    hSym->SetBinError(  iback, dx*sqrt(forw2+back2));
  }
  hAsy->Divide(hSym);
  hSym->Delete();
  return hAsy;
}// HstProjCA



///////////////////////////////////////////////////////////////////////////////////
TH1D *HstProjC(TString title, TH2D *Scat, int NbMax)
{
  // Projection onto c=cos(theta) axis, suming over v up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  TH1D *hTot = (TH1D*)Scat->ProjectionY(title,1,nbX,"e"); // option "e" doesnt work
  hTot->Reset();
  //
  double forw,forw2;
  double back,back2;
  double dx= (Xmax-Xmin)/nbX; // integration over X
  int ivMax;
  //
  if( (NbMax>0) && (NbMax<nbX) )
    ivMax = NbMax;
  else
    ivMax = nbX;
  int nbYhalf = nbY/2;
  for(int ic=1; ic <= nbYhalf; ic++){
    int iforw = nbYhalf+ic;
    int iback = nbYhalf-ic+1;
    forw=0.0;  forw2=0.0;
    back=0.0;  back2=0.0;
    for(int iv=0; iv <= ivMax; iv++){
      forw  += Scat->GetBinContent(  iv,iforw);
      forw2 += sqr(Scat->GetBinError(iv,iforw));
      back  += Scat->GetBinContent(  iv,iback);
      back2 += sqr(Scat->GetBinError(iv,iback));
    }
    // total
    hTot->SetBinContent(iforw,     dx*(forw));
    hTot->SetBinError(  iforw, dx*sqrt(forw2));
    hTot->SetBinContent(iback,     dx*(back));
    hTot->SetBinError(  iback, dx*sqrt(back2));
  }
  return hTot;
}// HstProjC




///////////////////////////////////////////////////////////////////////////////////
void MakeCumul(TH1D *hst1, TH1D *&hcum1)
{
  // makes cumulative distribution
  cout<<"Entering MakeCumul for  ";
  cout<< hst1->GetName() <<endl;
  int      nbX  = hst1->GetNbinsX();
  Double_t Xmax = hst1->GetXaxis()->GetXmax();
  Double_t Xmin = hst1->GetXaxis()->GetXmin();
  //
  hcum1 = (TH1D*)hst1->Clone("hcum1"); // allocate hcum1
  hcum1->Reset();
  double sum=0 ,sum2=0;
  double dx= (Xmax-Xmin)/nbX; // integration over X
  for(int iv=0; iv <= nbX; iv++){
    sum   += hst1->GetBinContent(  iv);
    sum2 += sqr(hst1->GetBinError(iv));
    hcum1->SetBinContent(iv,     dx*(sum));
    hcum1->SetBinError(  iv, dx*sqrt(sum2));
   }
}//MakeCumul

TH1D *HstCumul(TString title, TH1D *hst1)
{
  // makes cumulative distribution
  cout<<"Entering MakeCumul for  ";
  cout<< hst1->GetName() <<endl;
  int      nbX  = hst1->GetNbinsX();
  Double_t Xmax = hst1->GetXaxis()->GetXmax();
  Double_t Xmin = hst1->GetXaxis()->GetXmin();
  //
  TH1D *hcum1 = (TH1D*)hst1->Clone("hcum1"); // allocate hcum1
  hcum1->Reset();
  double sum=0 ,sum2=0;
  double dx= (Xmax-Xmin)/nbX; // integration over X
  for(int iv=0; iv <= nbX; iv++){
    sum   += hst1->GetBinContent(  iv);
    sum2 += sqr(hst1->GetBinError(iv));
    hcum1->SetBinContent(iv,     dx*(sum));
    hcum1->SetBinError(  iv, dx*sqrt(sum2));
   }
  hcum1->SetName(title);
  return hcum1;
}//HstRatio


  ///////////////////////////////////////////////////////////////////////////////////
void MakeAFB(TH1D *hAll, TH1D *&hAFB)
{
  // makes Afb(c) out of hAll(c)
  int      nbX  = hAll->GetNbinsX();
  Double_t Xmax = hAll->GetXaxis()->GetXmax();
  Double_t Xmin = hAll->GetXaxis()->GetXmin();
  //
  hAFB = (TH1D*)hAll->Clone("HstNew3"); // allocate hAFB
  hAFB->Reset();
  TH1D *hSym = (TH1D*)hAFB->Clone("HstNew4"); // local temporary histo
  hAFB->Reset();
  //
  double forw,forw2;
  double back,back2;
  double dx= (Xmax-Xmin)/nbX; // integration over X
  int nbXhalf = nbX/2;
  for(int ic=1; ic <= nbXhalf; ic++){
    int iforw = nbXhalf+ic;
    int iback = nbXhalf-ic+1;
    //cout<<"************** "<<iforw<<" "<<iback<<endl;
    forw  = hAll->GetBinContent(  iforw);
    forw2 = sqr(hAll->GetBinError(iforw));
    back  = hAll->GetBinContent(  iback);
    back2 = sqr(hAll->GetBinError(iback));
    // asymetric
    hAFB->SetBinContent(iforw,     dx*(forw-back));
    hAFB->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hAFB->SetBinContent(iback,    -dx*(forw-back));
    hAFB->SetBinError(  iback, dx*sqrt(forw2+back2));
    // symetric
    hSym->SetBinContent(iforw,     dx*(forw+back));
    hSym->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hSym->SetBinContent(iback,     dx*(forw+back));
    hSym->SetBinError(  iback, dx*sqrt(forw2+back2));
    cout<<"MakeAFB: (f-b)/2,f,b==="<<ic<<" "<< (forw-back)/2<<" b="<<back<<" f="<<forw<<endl;
  }
  hAFB->Divide(hSym);
  hSym->Delete();
  //for(int ic=1; ic <= nbX; ic++)
  //  cout<<"***"<<ic<<" "<< hAFB->GetBinContent(ic)<<endl;
}// ProjC


TH1D *HstDiff(TString title, TH1D *HST1, TH1D *HST2, Int_t kolor)
{
TH1D *Hd12 = (TH1D*)HST1->Clone(title);
Hd12->Add(HST1, HST2,    1.0, -1.0);
Hd12->SetLineColor(kolor);
return Hd12;
}//HstDiff

TH1D *HstAddi(TString title, TH1D *HST1, TH1D *HST2, Int_t kolor)
{
TH1D *Hd12 = (TH1D*)HST1->Clone(title);
Hd12->Add(HST1, HST2,    1.0,  1.0);
Hd12->SetLineColor(kolor);
return Hd12;
}//HstAddi

TH1D *HstRatio(TString title, TH1D *HST1, TH1D *HST2, Int_t kolor)
{
TH1D *Hd12 = (TH1D*)HST1->Clone(title);
Hd12->Divide(HST2);
Hd12->SetLineColor(kolor);
return Hd12;
}//HstRatio


TH1D *HstRatioSc(TString title, TH1D *HST1, TH1D *HST2, Double_t fact)
{
TH1D *Hd12 = (TH1D*)HST1->Clone(title);
Hd12->Divide(HST2);
Hd12->Scale(fact);
return Hd12;
}//HstDiff


TH1D *HstTildeAFB(TString title, TH1D *HST1, TH1D *HST2)
{
TH1D *Hst1 = HstCumul(title,HST1);       // total
TH1D *Hst2 = HstCumul("hst_test",HST2);  // weighted with cos(theta)
Hst1->Divide(Hst2);
Hst1->Scale(3.0/2.0);
return Hst1;
}//HstTildeAFB


TH1D *HstAFB(TString title, TH1D *HST1, TH1D *HST2)
{
TH1D *Hst1 = HstCumul(title,HST1);        // forward F
TH1D *Hst2 = HstCumul("hst_test2",HST2);  // total F+B
Hst1->Add(Hst1, Hst2,    2.0, -1.0);      // 2F-(F+B)=F-B
//
Hst1->Divide(Hst2);
return Hst1;
}//HstAFB

TH1D *HstAFB3(TString title, TH1D *HST1, TH1D *HST2, TH1D *HST3)
{
TH1D *Hst1 = HstCumul(title,HST1);        // forward F
TH1D *Hst2 = HstCumul("hst_test2",HST2);  // total F+B
TH1D *Hst3 = HstCumul("hst_test3",HST3);  // Deno
Hst1->Add(Hst1, Hst2,    2.0, -1.0);      // 2F-(F+B)=F-B
//
Hst1->Divide(Hst3); // (F-B)/Deno
return Hst1;
}//HstAFB3


TH1D *HstAFB4(TString title, TH1D *HST21F, TH1D *HST21, TH1D *HST2F, TH1D *HST2)
// Exact AFB for weight differences with in flight cumulating
{
TH1D *Hst21F = HstCumul(title,      HST21F);  // forward dF21
TH1D *Hst21  = HstCumul("hst_test2",HST21);   // total   dF21+dB21
TH1D *Hst2F  = HstCumul("hst_test3",HST2F);   // F2
TH1D *Hst2   = HstCumul("hst_test4",HST2);    // F2+B2
Hst21F->Add(Hst21F, Hst21,    2.0, -1.0);     // 2dF21-(dF21+dB21)=dF21-dB21
Hst21F->Divide(Hst2); // (dF21-dB21)/(F2+B2)
///////////////////////////////////////////////////
// Constructing second term approximately
TH1D *Hafb2  = HstCumul("hst_test5",HST2F);   // F2
Hafb2->Add(Hafb2, Hst2,    2.0, -1.0);        // 2F-(F+B)=F-B
Hafb2->Divide(Hst2);   // AFB2 ready to go
// Second term, approximate, AFB1->AFB2
Hafb2->Multiply(Hst21);;
Hafb2->Divide(Hst2);
// Subtract approximate second term
//Hst21F->Add(Hst21F, Hafb2,    1.0, -1.0);
///////////////////////////////////////////////////
// Constructing second term exactly
TH1D *Hafb1  = HstCumul("hst_test6",HST2F);  // Sig2F
Hafb1->Add(Hafb1, Hst2,    2.0, -1.0);       // numerator of AFB2, 2F-(F+B)
Hafb1->Add(Hafb1, Hst21F,  1.0, -2.0);       // correcting
Hafb1->Add(Hafb1, Hst21,   1.0,  1.0);       // numerator of AFB1
TH1D *Hsig1  = HstCumul("hst_test6",HST2);   // denominator  AFB1
Hsig1->Add(Hsig1, Hst21,   1.0, -1.0);       // denominator  AFB1
Hafb1->Divide(Hsig1);                        // AFB1 completed
// Second term, exact
Hafb1->Multiply(Hst21);;
Hafb1->Divide(Hst2);
// Subtract exact second term
Hst21F->Add(Hst21F, Hafb1,    1.0, -1.0);
///////////////////////////////////////////////////
return Hst21F;
}//HstAFB4


TH1D *HstAFB2cl(TString title, TH1D *HST1F, TH1D *HST1)
{
TH1D *Hst1 = (TH1D*)HST1F->Clone(title);  // total F+B
//
Hst1->Add(Hst1, HST1,    2.0, -1.0);      // 2F-(F+B)=F-B
//
Hst1->Divide(HST1);
return Hst1;
}//HstAFB


TH1D *HstAFB4cl(TString title, TH1D *HST21F, TH1D *HST21, TH1D *HST2F, TH1D *HST2)
// Exact AFB for weight differences with cloning (no cumulation)
{
TH1D *Hst21F = (TH1D*)HST21F->Clone(title);        // forward dF21
TH1D *Hst21  = (TH1D*)HST21->Clone( "hst_Test2" ); // total   dF21+dB21
TH1D *Hst2F  = (TH1D*)HST2F->Clone( "hst_Test3");  // F2
TH1D *Hst2   = (TH1D*)HST2->Clone(  "hst_Test4");  // F2+B2
Hst21F->Add(Hst21F, Hst21,    2.0, -1.0);     // 2dF21-(dF21+dB21)=dF21-dB21
Hst21F->Divide(Hst2); // (dF21-dB21)/(F2+B2)
///////////////////////////////////////////////////
// Constructing second term approximately
//TH1D *Hafb2  = HstCumul("hst_test5",HST2F);   // F2
TH1D *Hafb2  = (TH1D*)HST2F->Clone("hst_Test5"); // F2
Hafb2->Add(Hafb2, Hst2,    2.0, -1.0);        // 2F-(F+B)=F-B
Hafb2->Divide(Hst2);   // AFB2 ready to go
// Second term, approximate, AFB1->AFB2
Hafb2->Multiply(Hst21);;
Hafb2->Divide(Hst2);
// Subtract approximate second term
//Hst21F->Add(Hst21F, Hafb2,    1.0, -1.0);
///////////////////////////////////////////////////
// Constructing second term exactly
//TH1D *Hafb1  = HstCumul("hst_test6",HST2F);  // Sig2F
TH1D *Hafb1  = (TH1D*)HST2F->Clone("hst_Test6");
Hafb1->Add(Hafb1, Hst2,    2.0, -1.0);       // numerator of AFB2, 2F-(F+B)
Hafb1->Add(Hafb1, Hst21F,  1.0, -2.0);       // correcting
Hafb1->Add(Hafb1, Hst21,   1.0,  1.0);       // numerator of AFB1
//TH1D *Hsig1  = HstCumul("hst_test6",HST2);   // denominator  AFB1
TH1D *Hsig1  = (TH1D*)HST2->Clone("hst_test6");
Hsig1->Add(Hsig1, Hst21,   1.0, -1.0);       // denominator  AFB1
Hafb1->Divide(Hsig1);                        // AFB1 completed
// Second term, exact
Hafb1->Multiply(Hst21);;
Hafb1->Divide(Hst2);
// Subtract exact second term
Hst21F->Add(Hst21F, Hafb1,    1.0, -1.0);
///////////////////////////////////////////////////
return Hst21F;
}//HstAFB4cl

