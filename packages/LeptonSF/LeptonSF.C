#include "LeptonSF.h"
#include "TFile.h"

#include <iostream>

LeptonSF::LeptonSF(bool loadhistos):
  fTightMuonIDSF(0),
  fTightMuonIsoSF(0),
  fTightMuonSF(0),
  fTightElectronIDSF(0),
  fTightElectronSF(0),
  fDoubleMuSF(0),
  fDoubleElSF(0),  
  fMuEGSF(0) {

  if (loadhistos) {
    // Load Lepton SFs
    LoadTightMuonIDSF();
    LoadTightMuonIsoSF();
    LoadTightMuonSF();
    LoadTightElectronIDSF();
    LoadTightElectronSF();

    // Load Trigger SFs
    LoadDoubleElectronSF();
    LoadDoubleMuonSF();
    LoadElMuSF();
  }
}




TH2D* LeptonSF::LoadTightMuonIDSF(const char* file, 
					 const char* histo) {
  fTightMuonIDSF = GetHistogramFromFileD(file, histo, "fTightMuonIDSF");
  return fTightMuonIDSF;
}

TH2D* LeptonSF::LoadTightMuonIsoSF(const char* file, 
				   const char* histo) {
  fTightMuonIsoSF = GetHistogramFromFileD(file, histo, "fTightMuonISOSF");
  return fTightMuonIsoSF;
}
TH2D* LeptonSF::LoadTightMuonSF(const char* file, 
				const char* histo) {
  fTightMuonSF = GetHistogramFromFileD(file, histo, "fTightMuonSF");
  return fTightMuonSF;
}

TH2D* LeptonSF::LoadTightElectronIDSF(const char* file, 
				      const char* histo) {
  fTightElectronIDSF = GetHistogramFromFileD(file, histo, "fTightElectronIDSF");
  return fTightElectronIDSF;
}
TH2D* LeptonSF::LoadTightElectronSF(const char* file, 
				    const char* histo) {
  fTightElectronSF = GetHistogramFromFileD(file, histo, "fTightElectronSF");
  return fTightElectronSF;
}


TH2F* LeptonSF::LoadDoubleMuonSF(const char* file, 
				 const char* histo){
  fDoubleMuSF = GetHistogramFromFileF(file, histo, "fDoubleMuSF");
  return fDoubleMuSF;
}

TH2F* LeptonSF::LoadDoubleElectronSF(const char* file, 
				     const char* histo){
  fDoubleElSF = GetHistogramFromFileF(file, histo, "fDoubleElSF");
  return fDoubleElSF;
}

TH2F* LeptonSF::LoadElMuSF(const char* file, 
			   const char* histo){
  fMuEGSF = GetHistogramFromFileF(file, histo, "fMuEGSF");
  return fMuEGSF;
}






TH2D* LeptonSF::GetHistogramFromFileD(const char* filename, 
				      const char* histoname, 
				      const char* newhname) {
  TFile* file  = TFile::Open(filename);
  if (!file) {
    std::cerr << "ERROR[LeptonSF]: Could not load file" << std::endl
	      << "                 " << filename << std::endl;
    return 0;
  }
  TH2D* h = (TH2D*) file->Get(histoname)->Clone(newhname);
  if (!h) {
    std::cerr << "ERROR[LeptonSF]: Could not find histogram " 
	      << histoname << std::endl;
    return 0;
  }
  h->SetDirectory(0);
  file->Close();
  
  return h;
}

TH2F* LeptonSF::GetHistogramFromFileF(const char* filename, 
				      const char* histoname, 
				      const char* newhname) const {
  TFile* file  = TFile::Open(filename);
  if (!file) {
      std::cerr << "ERROR[LeptonSF]: Could not load file" << std::endl
		<< "                 " << filename << std::endl;
      return 0;
  }
  TH2F* h = (TH2F*) file->Get(histoname)->Clone(newhname);
  if (!h) {
    std::cerr << "ERROR[LeptonSF]: Could not find histogram " 
	      << histoname << std::endl;
    return 0;
    }
  h->SetDirectory(0);
  file->Close();
  
  return h;
}

TCanvas* LeptonSF::Draw() {
  //gStyle->SetOptStat(0);
  //gStyle->SetPalette(1);
  //gStyle->SetPaintTextFormat("4.3f");
  TCanvas* c = new TCanvas();
  c->Divide(2,1);
  
  c->cd(1);
  fTightElectronSF->SetTitle("Electron Id/Iso SF");
  fTightElectronSF->SetXTitle("#eta");
  fTightElectronSF->SetYTitle("P_{t}");
  fTightElectronSF->GetYaxis()->SetRange(1, 3);
  //fTightElectronSF->GetYaxis()->SetRangeUser(0.,100.);
  fTightElectronSF->Draw("COLZ, TEXTE");
  
  /*c->cd(2);
  fTightMuonSF->SetTitle("Tight Muon ID");
  fTightMuonSF->SetXTitle("P_{t}");
  fTightMuonSF->SetYTitle("#eta");
  fTightMuonSF->GetXaxis()->SetRange(1,8);
  fTightMuonSF->Draw("COLZ, TEXTE");
  */
  
  c->cd(2);
  fMuEGSF->SetTitle("MuEG Trigger SF");
  fMuEGSF->GetXaxis()->SetNdivisions(505);
  fMuEGSF->Draw("COLZ, TEXTE");
  
  c->SaveAs("sfs.png");

  return c;
}
