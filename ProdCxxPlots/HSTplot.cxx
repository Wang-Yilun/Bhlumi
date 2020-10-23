//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   HSTplot                                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include "HSTplot.h"

HSTplot::HSTplot(const char* Name)
{
   cout<< "----> HSTplot USER Constructor "<<endl;
   m_pi = 3.141593e0;
}

///////////////////////////////////////////////////////////////////////////////////
void HSTplot::Initialize(){
//------------------------------------------------------------------------
//------------------------------------------------------------------------
  cout<<"====================================================================="<<endl;
  cout<<"================ HSTplot::initialization ============================"<<endl;

  fmt0 = "$  %10.2f $";
  fmt1 = "& $ %10.4f \\pm %8.4f $ ";
  fmt2 = "& $ %10.4f \\pm %8.4f $ ";
  fmt3 = "& $ %10.4f $ ";

}//Initialize


void HSTplot::PlInit(FILE *ltx, int lint)
{
//----------------------------------------------------------------------
// Lint =0     Normal mode, full LaTeX header
// Lint =1     For TeX file is used in \input, no  LaTeX header
// Lint =2     LaTeX header for one-page plot used as input for postscript
// Negative Lint only for debug, big frame around plot is added.
//----------------------------------------------------------------------
//[[[   m_lint=lint;
//[[[if( abs(lint) == 0){
//[[[// Normal mode, no colors!!!
//[[[   fprintf(ltx,"\\documentclass[12pt]{article}\n");
//[[[   fprintf(ltx,"\\textwidth  = 16cm\n");
//[[[   fprintf(ltx,"\\textheight = 18cm\n");
//[[[   fprintf(ltx,"\\begin{document}\n");
//[[[   fprintf(ltx,"  \n");
//[[[} else if( abs(lint) == 1) {
//[[[// For TeX file is used in \input
//[[[   fprintf(ltx,"  \n");
//[[[} else if( abs(lint) == 2){
// For one-page plot being input for postscript
   fprintf(ltx,"\\documentclass[12pt,dvips]{article}\n");
   fprintf(ltx,"\\usepackage{amsmath}\n");
   fprintf(ltx,"\\usepackage{amssymb}\n");
   fprintf(ltx,"\\usepackage{epsfig}\n");
   fprintf(ltx,"\\usepackage{epic}\n");
   fprintf(ltx,"\\usepackage{eepic}\n");
   fprintf(ltx,"\\usepackage{color}\n"); //<-for colors!!!
   fprintf(ltx,"\\begin{document}\n");
   fprintf(ltx,"\\pagestyle{empty}\n");
   fprintf(ltx,"  \n");
//[[[} else {
//[[[   cout<<"+++STOP in GLK_PlInt, wrong lint =" <<lint<< endl;
//[[[}// lint
}//GLK_PlCap

void HSTplot::PlTable2(int Ncol, TH1D *iHst[], FILE *ltex, Char_t *Capt[], Char_t Mcapt[] , const char *chr1, int k1,int k2,int dk)
//* Tables in TeX, up to 9 columns
//* Ncol          = numbers of columns/histograms
//* idl(1:Npl)    = list of histo id's
//* ccapt(1:Npl+1)= list of column-captions above each column
//* mcapt         = multicolumn header, none if mcapt=' ',
//* chr1          = ' ' normal default, = Header+Table+Ending
//*               = 'B' no page eject,  = Header+Table
//*               = 'E' no page eject,  =        Table+Ending
//*               = 'E' no page eject,  =        Table
//* k1,k2,dk      = range of bins is (k1,k2) with increment dk
{
  int Npl=abs(Ncol);
  if( chr1 == " " || chr1 == "B"){
	  //------------------------------!
	  //           Header
	  //------------------------------!
	  fprintf(ltex," \n");
	  fprintf(ltex,"% ========================================\n");
	  fprintf(ltex,"% ============ begin table ===============\n");
	  //
//
	  fprintf(ltex,"\\noindent\n");

	  //------------------------------!
	  // Tabular header
	  //------------------------------!
      fprintf(ltex,"\\begin{tabular}{|");
	  for(int i=0; i<=Npl; i++ ) fprintf(ltex,"|r");
	  fprintf(ltex,"||}\n");
	  fprintf(ltex,"\\hline\\hline\n");

	  //------------------------------!
	  // Captions in columns
	  //------------------------------!
	  fprintf(ltex,"%s  \n", Capt[0]);
	  for(int i=1; i<=Npl; i++ ) fprintf(ltex,"& %s \n", Capt[i]);

	  fprintf(ltex,"\\\\ \\hline\n");
  }

  //------------------------------------------!
  // Optional Multicolumn caption
  //------------------------------------------!
  if(Ncol>0){
     fprintf(ltex,"& \\multicolumn{ %i }{c||}{",Npl);
     fprintf(ltex,"  %s  } \\\\   \\hline\n", Mcapt);
  }//Mcapt

  // X range taken from 1-st histogram
  int      nbX  = (iHst[1])->GetNbinsX();
  Double_t Xmax = (iHst[1])->GetXaxis()->GetXmax();
  Double_t Xmin = (iHst[1])->GetXaxis()->GetXmin();
  double dx = (Xmax-Xmin)/nbX;
  cout<<"  nbX=  " <<nbX<<endl;
  cout<<"  Xmin= " <<Xmin<<endl;
  cout<<"  Xmax= " <<Xmax<<endl;
  // Raws

  double xk, yi, ei;
  for(int k=k1; k<=k2; k+=dk){    // loop over bin (raw) number
	cout<<" k="<<k<<endl;
	xk = Xmin +dx*k; // right edge
    fprintf(ltex,fmt0, xk);
    for( int j=1; j<=Npl; j++ ){
	    yi = (iHst[j])->GetBinContent(k);
	    ei = (iHst[j])->GetBinError(k);
	    if( ei != 0 ){
	      if( j == 1){
		    fprintf(ltex, fmt1, yi, ei); // 1st column
    	  }else{
		    fprintf(ltex, fmt2, yi, ei); // next columns
    	  }
	    } else{
			fprintf(ltex, fmt3, yi);     // zero error
	    }
     }//j
    fprintf(ltex,"\\\\ \n"); // double slash and newline
 }//k
 fprintf(ltex,"\\hline\n");

  //------------------------------!
  // Ending
  //------------------------------!
  if( chr1 == " " || chr1 == "E"){
	  fprintf(ltex,"\\end{tabular}\n");
	  fprintf(ltex,"% ============ end   table ===============\n");
	  fprintf(ltex,"% ========================================\n");
	  fprintf(ltex," \n");
  }//chr1

}//GLK_PlTable2


void HSTplot::PlEnd(FILE *ltex)
{//---------------------------------------------------
// Note that TeX file is used in \input then you may not want
// to have header and \end{document}
   fprintf(ltex,"\\end{document}");
}//HSTplot::PlEnd



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                End of Class HSTplot                                        //
///////////////////////////////////////////////////////////////////////////////

