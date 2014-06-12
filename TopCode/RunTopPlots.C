//#include "TopPlotter.h"

void RunTopPlots(TString pathtofile, Int_t verbose=0){
  TString outputdir = pathtofile + "/TopPlots/";
  
  cout << "--------------" << endl;
  cout << "OutputDir is:      " << outputdir << endl;
  cout << "Verbose level is:  " << verbose << endl;
  cout << "--------------" << endl;
  
  gROOT->LoadMacro("~folgueras/TOP/TopCode/tdrstyle.C");
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
  
  //  cout << "Printing Yields with MC estimations..." << endl;
  //  tA->PrintYieldsWithMC();

  //  cout << "Draw Kinematic Plots with MC estimations..." << endl;
  // tA->DrawKinematicPlotsWithMC(-1, NBTagsNJets, -1);
  // tA->DrawKinematicPlotsWithMC();
  
  cout << "Calculate Data-Driven backgrounds" << endl;
  tA->CalculateNonWZLeptonsBkg();
  tA->CalculateDYBkg();
  
  cout << "Draw Plots for Likelihood..." << endl;
  tA->DrawNbjetsNjets(false);
  tA->DrawNbjetsNjets(true);
  
  cout << "Saving Plots for Likelihood..." << endl;
  tA->SaveHistosForLH(false);
  tA->SaveHistosForLH(true);
  
  
  cout << "Calculating cross section ... " << endl;
  tA->CalculateCrossSection(false);
  tA->CalculateCrossSection(true );
  
  delete tA;
}
