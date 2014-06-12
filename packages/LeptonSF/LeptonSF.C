#include "LeptonSF.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"

#include <iostream>

LeptonSF::LeptonSF():
  fMuonSF(0),  fElecSF(0),
  fDoubleMuSF(0),  fDoubleElSF(0),  fMuEGSF(0) {
  
  TString path = "http://www.hep.uniovi.es/folgueras/TOP/LeptonSF/";
  LoadMuonSF(path + "ReReco_SF_Mu.root");
  LoadElecSF(path + "ReReco_SF_EG.root");
  LoadElMuSF(path + "ReReco_trigger_SF_emu.root");
}
void LeptonSF::LoadMuonSF(TString file){
  fMuonSF     = GetHistogramFromFile(file, "GlobalSF", "fMuonSF");
  fDoubleMuSF = GetHistogramFromFile(file, "scalefactor eta2d with syst", "fDoubleMuSF");
  return;
}
void LeptonSF::LoadElecSF(TString file){
  fElecSF     = GetHistogramFromFile(file, "GlobalSF", "fElecSF");
  fDoubleElSF = GetHistogramFromFile(file, "scalefactor eta2d with syst", "fDoubleElSF");
  return;
}
void LeptonSF::LoadElMuSF(TString file){
  fMuEGSF = GetHistogramFromFile(file, "scalefactor eta2d with syst", "fMuEGSF");
  return;
}
TH2F* LeptonSF::GetHistogramFromFile(TString filename, 
				     TString histoname, 
				     TString newhname) const {
  TFile* file  = TFile::Open(filename);
  if (!file) {
      std::cerr << "ERROR[TopSelector]: Could not load file" << std::endl
		<< "                        " << filename << std::endl;
      return 0;
  }
  TH2F* h = (TH2F*) file->Get(histoname)->Clone(newhname);
  if (!h) {
    std::cerr << "ERROR[TopSelector]: Could not find histogram " 
	      << histoname << std::endl;
    return 0;
    }
  h->SetDirectory(0);
  file->Close();
  
  return h;
}
