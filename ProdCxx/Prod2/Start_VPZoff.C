{
///============================================================================
///
/// This is configuration/initialization script for ../MainKKMC main program
/// To start MC run in the interactive mode just type "make start" 
///
/// This series of calculations integrates a given distribution 
/// by using MC methods
///============================================================================
gROOT->Reset();
cout<<"%%% ================== Start.C ================== %%%%"<<endl;
gSystem->Load("../.libs/libProd.so");
//gSystem->Load("../../6xx/.libs/libBHL6.so");
TFile HistoFile("histo.root","RECREATE","Histograms");
TFile GenFile(  "mcgen.root","RECREATE","r.n.generator, MCgens");
TFile SemFile(  "semaf.root","RECREATE","Semaphore");
///*****************************************************************
////*****************************************************************
///   Create new instance of Semaphore object
///   and fill it with the MC run general parameters
TString semaf   = "START";
double nevtot   = 1e12; // 1000G
//nevtot = 1e6;
double nevgrp   = 5e5; // 500k
///------------------------------------------------------------------
SemFile.cd();
TSemaf *Semafor = new TSemaf(semaf, nevtot, nevgrp);
Semafor->Write("Semafor",TObject::kOverwrite);
SemFile.Write();
SemFile.Close();
///*****************************************************************
///*****************************************************************
///       Create new instance of RN generator and initialize it
GenFile.cd();
TRandom *RN_gen = new TRandom3();       // Central r.n.gen.
long    iniseed = 54217137;
iniseed = 46785717;
RN_gen->SetSeed(iniseed);
RN_gen->Write("RN_gen",TObject::kOverwrite);
///*****************************************************************
///*****************************************************************
///      Create new instance of MC generator
TMCgenBHL6 *MCgen = new TMCgenBHL6("MCgen");
//########### Change some input parameters ###########
MCgen->m_CMSENE  =  92.3e0;   //  CMSENE
MCgen->m_angmin  =  0.022;    //  theta_min [rad]   generation
MCgen->m_angmax  =  0.082;    //  theta_max [rad]   generation
//####################################################
MCgen->ls();
MCgen->Write("MCgen",TObject::kOverwrite);
///*****************************************************************
///*****************************************************************
///       Create new instance of the MC analysis object
TRobolProd *RoboT = new TRobolProd("RoboT");  /// base clase only
// Parameters of the detectors
RoboT->m_TminF  =     .024;    //  fidutial angle theta_min
RoboT->m_TmaxF  =     .058;    //  fidutial angle theta_max
RoboT->f_HistNormName = "HST_BHL_NORMA"; // Name of normalization histo
RoboT->Write("RoboT",TObject::kOverwrite);
///*****************************************************************
GenFile.Write();
cout<<"---------------------------------------------------------"<<endl;
GenFile.Close();
///*****************************************************************
///*****************************************************************
cout << "===========Output written in histo.root===========" << endl;
HistoFile.Write();
HistoFile.ls();
HistoFile.Close();
cout<<"%%% ===============End Start.C ================== %%%%"<<endl;
return 0;
}
