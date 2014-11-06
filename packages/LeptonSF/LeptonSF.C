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
  TCanvas* c = new TCanvas();
  c->Divide(3,2);
  c->cd(1);
  fTightElectronIDSF->SetTitle("Tight Electron ID");
  fTightElectronIDSF->SetXTitle("#eta");
  fTightElectronIDSF->SetYTitle("P_{t}");
  fTightElectronIDSF->Draw("COL, TEXT");
  c->cd(2);
  fTightMuonIDSF->SetTitle("Tight Muon ID");
  fTightMuonIDSF->SetXTitle("#eta");
  fTightMuonIDSF->SetYTitle("P_{t}");
  fTightMuonIDSF->Draw("COL, TEXT");
  c->cd(3);
  fTightMuonIsoSF->SetTitle("Tight Muon Isolation");
  fTightMuonIsoSF->SetXTitle("#eta");
  fTightMuonIsoSF->SetYTitle("P_{t}");
  fTightMuonIsoSF->Draw("COL, TEXT");

  c->cd(4);
  fDoubleElSF->SetTitle("Double Electron");
  fDoubleElSF->Draw("COL, TEXT");
  c->cd(5);
  fDoubleMuSF->SetTitle("Double Muon");
  fDoubleMuSF->Draw("COL, TEXT");
  c->cd(6);
  fMuEGSF->SetTitle("MuEG Trigger");
  fMuEGSF->Draw("COL, TEXT");

  return c;
}
