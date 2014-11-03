//#include "TopPlotter.h"

void RunTopPlots(TString pathtofile, Int_t verbose=0){
  TString outputdir = pathtofile + "/TopPlots/";
  
  cout << "--------------" << endl;
  cout << "OutputDir is:      " << outputdir << endl;
  cout << "Verbose level is:  " << verbose << endl;
  cout << "--------------" << endl;
  
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
  
  cout << "Loading TopPlotter.cc..."<< endl;
  gROOT->LoadMacro("TopPlotter.cc+");
  
  cout << "Creating TopPlotter class..."<<endl;
  TopPlotter *tA = new TopPlotter();

  cout << "Init TopPlotter..."<<endl;
  tA->Init(pathtofile);

  cout << "Set OutputDir and Vervose level..."<<endl;
  tA->SetOutputDir(outputdir);
  tA->SetVerbose(verbose);
  
  tA->Loop();
  
  delete tA;
}
