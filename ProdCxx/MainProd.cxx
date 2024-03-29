///////////////////////////////////////////////////////////////////////////////
// This module handles main Monte Carlo loop and writing histograms onto disk.
// It is UNIVERSAL, works for any MC generator for event analysis.
// MC event generator and analysis is encapsulated inside RobolA module.
// Three different disk files are used
///////////////////////////////////////////////////////////////////////////////
// C++ headers
using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"

// OUR headers
#include "TMCgenBHL6.h"
#include "TSemaf.h"
#include "TRobolProd.h"
#include "TSystem.h"

int main()
{
  gSystem->Load("../.libs/libProd.so");
  TFile *HistoFile = new TFile("histo.root","UPDATE","histograms");
  TFile *MCgenFile = new TFile("mcgen.root","UPDATE","Generators");
  cout<<"------------------------------MCgenFile->ls----------------------------------"<<endl;
  MCgenFile->ls();
  //MCgenFile->GetListOfKeys()->Print();
  //MCgenFile->ShowStreamerInfo();
  //
  char chcyc[100];          // Cycle text variable
  Text_t *Tcycle;           // Cycle text variable
  TFile *SemafFile;    // Semaphore and loop params
  SemafFile =TFile::Open("semaf.root","UPDATE","Semaphor");
  ofstream   OutFile("pro.output",ios::out);  // Logfile output
  cout    <<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  cout    << "  | MainProd  Entering | "<<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  OutFile <<endl<<flush;
  OutFile << "  |--------------------| "<<endl<<flush;
  OutFile << "  | MainProd  Entering | "<<endl<<flush;
  OutFile << "  |--------------------| "<<endl<<flush;
  double NevTot =2e14; // Total (limit) of events, redefined below
  double NGroup =2e5;  // no. of events in the group, TO BE TUNED
  long   iLoop;        // external loop event counter
  long   iGroup;       // internal loop event counter
  double iEvent = 0;   // serial number of event in the MC generation
  //////////////////////////////////////////////////////////////////////////
  TString Status = "START";   // to be overwitten
  //------------------------------------------------------------------------
  SemafFile->cd();
  SemafFile->ls();
  TSemaf *Semafor = (TSemaf*)SemafFile->Get("Semafor"); /// read semafor from the disk
  Status   = Semafor->m_flag;
  NevTot   = Semafor->m_nevtot;
  NGroup   = Semafor->m_nevgrp;
  cout<<"  MainPr:  Status  = "<<  Status << endl;
  cout<<"  MainPr:  NevTot  = "<<  NevTot << endl;
  Semafor->m_flag = "CONTINUE";
  Semafor->Write("Semafor",TObject::kOverwrite);      /// write modified semaphore
  SemafFile->Write("",TObject::kOverwrite,0);         /// write modified semaphore
  SemafFile->ls();
  SemafFile->Close();
  //--------------------------------------------------------------------------
  if( Status == "STOP") cout<<"MainPr: +++Warning+++ Status=STOP in semaphore"<<endl;
  //////////////////////////////////////////////////////////////////////////
  TRobolProd *RoboT =(TRobolProd*)MCgenFile->Get("RoboT");   /// read analysis object
  RoboT->Initialize(
        &OutFile,   /// Central log-file for messages
        MCgenFile,  /// ROOT disk file for CRNG and MC gen.
        HistoFile); /// ROOT disk file for histograms
  TMCgenBHL6 *MCgen = (TMCgenBHL6*)RoboT->f_MCgen;
  cout<< "||||||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
  cout<< "|| MainBHL: to be generated "<<NevTot<< " events "<<endl;
  cout<< "||||||||||||||||||||||||||||||||||||||||||||||||||"<<endl;

  //[[[[[[
  // Strange workaround for ubuntu16 linker to see f77 externals !!!!
  MCgen->Print1();
  //]]]]]]

  ///***********************************************************************
  ///          MAIN LOOP OVER MONTE CARLO EVENTS
  for ( iLoop = 0; ; iLoop++)
    {
      for ( iGroup = 0; iGroup< NGroup; iGroup++) 
	{
	  iEvent++;
	  ///***************************
	  RoboT->Production(iEvent); //
	  ///***************************
	  if( iEvent  >= NevTot )  break;  // production stoped
	}; // for iGroup
      //
      cout    <<"iEvent = "<<iEvent<<endl<<flush;
      OutFile <<"iEvent = "<<iEvent<<endl<<flush;
      ///////////////////////////////////////////////////////////////
      // Disk Write for all histos and r.n. generator
      sprintf(chcyc,"*;%i",iLoop); Tcycle = chcyc;
      //cout<<"Main iLoop, Tcycle = "<< iLoop <<"|||"<< Tcycle <<"|||"<<endl;
      HistoFile->cd();
      HistoFile->Write();  //!!!
      HistoFile->Delete(Tcycle);
      HistoFile->Flush();
      //--------------------------------------------------
      MCgenFile->cd();
      MCgen->Write("MCgen",TObject::kOverwrite);  //
      RoboT->Write("RoboT",TObject::kOverwrite);  //
      MCgenFile->Write();  //!!!
      MCgenFile->Flush();
      //[[[[[ DEBUG !!!!
      //cout<<"-------------(i)-----------------MCgenFile->ls----------------------------------"<<endl;
      //MCgenFile->ls();
      //]]]]]
      ///////////////////////////////////////////////////////////////
      // Checks semaphore for each group of NGroup events
      SemafFile =TFile::Open("semaf.root","UPDATE","Semaphor");
      Semafor = (TSemaf*)SemafFile->Get("Semafor");  // read semafor from the disk
      Status   = Semafor->m_flag;                   // get flag
      Semafor->m_nevent = iEvent;                   // record no. of events
      Semafor->Write("Semafor",TObject::kOverwrite);// write modified semaphore
      SemafFile->Write("",TObject::kOverwrite);
      SemafFile->Close();
      //-------------------------------------------------------------
      if( iEvent  >= NevTot )      break;  // production stoped
      if( Status == "STOP" )       break;  // production stoped
      ///////////////////////////////////////////////////////////////
    } // for iLoop
  ///            End of MC LOOP
  ///************************************************************************
  SemafFile =TFile::Open("semaf.root","UPDATE","Semaphor");
  Semafor = (TSemaf*)SemafFile->Get("Semafor");
  Semafor->m_flag = "CONTINUE";  
  Semafor->Write("Semafor",TObject::kOverwrite);// write modified semaphore
  SemafFile->Write("",TObject::kOverwrite);
  SemafFile->Close();
  cout<< "||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
  cout<< "||  MainBHL: Generated "<<iEvent<< " events    "<<endl;
  cout<< "||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
  //
  RoboT->Finalize();
  //////////////////////////////////////////////////////////////////////
  //  LAST WRITEs
  MCgenFile->cd();
  MCgen->Write("MCgen",TObject::kOverwrite);  //
  RoboT->Write("RoboT",TObject::kOverwrite);  //
  //////////////////////////////////////////////////////////////////////
  //      Some TESTS and control printouts
  cout<<"------------------------------HistoFile->ls----------------------------------"<<endl;
  HistoFile->ls();
  //cout<<"------------------------------HistoFile->Map---------------------------------"<<endl;
  //HistoFile->Map();
  //cout<<"-------------------------MCgenFile->ShowStreamerInfo-------------------------"<<endl;
  //MCgenFile->ShowStreamerInfo();
  cout<<"------------------------HistoFile->GetListOfKeys-----------------------------"<<endl;
  HistoFile->GetListOfKeys()->Print();
  cout<<"------------------------MCgenFile->GetListOfKeys-----------------------------"<<endl;
  MCgenFile->GetListOfKeys()->Print();
  cout<<"-------------------------------SemafFile.ls---------------------------------"<<endl;
  SemafFile->ls();
  cout<<"------------------------------MCgenFile->ls----------------------------------"<<endl;
  MCgenFile->ls();
  cout<<"---------------------------------end---------------------------------------"<<endl;
  //////////////////////////////////////////////////////////////////////
  //                 CLOSE ALL FILES
  HistoFile->Close();
  MCgenFile->Close();
  OutFile.close();
  cout    <<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  cout    << "  |  MainBHL   Ended   | "<<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  OutFile <<endl<<flush;
  OutFile << "  |--------------------| "<<endl<<flush;
  OutFile << "  |  MainBHL   Ended   | "<<endl<<flush;
  OutFile << "  |--------------------| "<<endl<<flush;
  return 0;
  }
