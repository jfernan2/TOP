///////////////////////////////////////////////////////////////////////
//
//    FILE: LeptonSF.h
//   CLASS: LeptonSF
// AUTHORS: S. Folgueras, I. Gonzalez
//    DATE: February, 2014
//
// CONTENT: An utility class to handle Lepton and Trigger Scale Factors.
//
///////////////////////////////////////////////////////////////////////
#ifndef LEPTONSF_H
#define LEPTONSF 1

#include "TH2D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"

class LeptonSF {
 public:
  LeptonSF(bool loadhistos = true);
  ~LeptonSF() {}
  
  
  // Muon SFs
  double GetTightMuonIDSF (double pt, double eta) { 
    return fTightMuonIDSF->GetBinContent(fTightMuonIDSF->FindBin(TMath::Abs(eta), pt));
  }
  double GetTightMuonIsoSF(double pt, double eta) { 
    return fTightMuonIsoSF->GetBinContent(fTightMuonIsoSF->FindBin(TMath::Abs(eta), pt));
  }
  double GetTightMuonSF   (double pt, double eta) { 
    return fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(pt,TMath::Abs(eta)));
  }
  
  // Electron SFs
  float GetTightElectronIDSF  (float pt, float eta) const { // binned in eta, pt
    return fTightElectronIDSF->GetBinContent(fTightElectronIDSF->FindBin(TMath::Abs(eta), pt));
  }
  float GetTightElectronSF  (float pt, float eta) const { // binned in eta, pt
    return fTightElectronSF->GetBinContent(fTightElectronSF->FindBin(eta, pt));
  }

  
  // Trigger SFs
  float GetDoubleMuSF(float eta1, float eta2) const { // binned in eta1, eta2
    return fDoubleMuSF->GetBinContent(fDoubleMuSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }
  float GetDoubleElSF(float eta1, float eta2) const { // binned in eta1, eta2
    return fDoubleElSF->GetBinContent(fDoubleElSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }
  float GetMuEGSF    (float eta1, float eta2) const { // binned in eta1, eta2
    return fMuEGSF->GetBinContent(fMuEGSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }
  // Trigger SF errors
  float GetDoubleMuSF_err(float eta1, float eta2) const { // binned in eta1, eta2
    return fDoubleMuSF->GetBinError(fDoubleMuSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }
  float GetDoubleElSF_err(float eta1, float eta2) const { // binned in eta1, eta2
    return fDoubleElSF->GetBinError(fDoubleElSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }
  float GetMuEGSF_err(float eta1, float eta2) const { // binned in eta1, eta2
    return fMuEGSF->GetBinError(fMuEGSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }

  

  
  /// Methods to get the ERROR 
  double GetTightMuonSF_err(double pt, double eta) {
    return fTightMuonSF->GetBinError(fTightMuonSF->FindBin(TMath::Abs(eta), pt));
  }
  double GetTightElectronSF_err(double pt, double eta) {
    return fTightElectronSF->GetBinError(fTightElectronSF->FindBin(eta, pt));
  }
  
  // Methods to load Scale Factors
  // + Lepton SFs
  TH2D* LoadTightMuonIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/MuonIDSF.root",
			   const char* histo = "DATA_over_MC_Tight");
  TH2D* LoadTightMuonIsoSF(const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/MuonISOSF.root",
			   const char* histo = "DATA_over_MC_Tight");
  TH2D* LoadTightMuonSF(const char* file = "http://www.hep.uniovi.es/folgueras/TOP/LeptonSF/LegacySF/muonEffTop2D.root",
			const char* histo = "eff_pt_eta_fullwithtrackerr");
  
//  TH2D* LoadTightElectronIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/ElectronSF.root",
//			       const char* histo = "sfTIGHT");
  TH2D* LoadTightElectronIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/ElectronMVAIDIsoSF.root",
			       const char* histo = "electronsDATAMCratio_FO_ID_ISO");
  
  TH2D* LoadTightElectronSF (const char* file = "http://www.hep.uniovi.es/folgueras/TOP/LeptonSF/LegacySF/ElecSF_tight_inclerrs.root",
			     const char* histo = "GlobalSF");

  // + Trigger SFs
  TH2F* LoadDoubleElectronSF(const char* file = "http://www.hep.uniovi.es/folgueras/TOP/LeptonSF/LegacySF/triggerSummary_ee_tightleps.root",
 			     const char* histo = "scalefactor_eta2d_with_syst");
  TH2F* LoadDoubleMuonSF(const char* file = "http://www.hep.uniovi.es/folgueras/TOP/LeptonSF/LegacySF/triggerSummary_mumu_tightleps.root",
 			 const char* histo = "scalefactor_eta2d_with_syst");
  TH2F* LoadElMuSF(const char* file = "http://www.hep.uniovi.es/folgueras/TOP/LeptonSF/LegacySF/triggerSummary_emu_tightleps.root",
 		   const char* histo = "scalefactor_eta2d_with_syst");


  // Plot all histograms
  TCanvas* Draw();

 protected:
  TH2D* GetHistogramFromFileD(const char* file, const char* histo, 
			      const char* hname);
  TH2F* GetHistogramFromFileF(const char* file, const char* histo, 
			      const char* hname) const;
  
 private:
  // Lepton SFs
  TH2D*  fTightMuonIDSF;     //Muon ID Scale Factors
  TH2D*  fTightMuonIsoSF;    //Muon Isolation Scale Factors
  TH2D*  fTightMuonSF;    //Muon Isolation Scale Factors
  TH2D*  fTightElectronIDSF; //Electron ID Scale Factors
  TH2D*  fTightElectronSF; //Electron ID Scale Factors

  // Trigger SFs
  TH2F *fDoubleMuSF;
  TH2F *fDoubleElSF;
  TH2F *fMuEGSF;
};
#endif
