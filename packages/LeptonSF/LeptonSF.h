#ifndef LEPTONSF_H
#define LEPTONSF 1

#include "TH2F.h"
#include "TString.h"
#include "TMath.h"

class LeptonSF {
 public:
  LeptonSF();
  ~LeptonSF() {}
  
  void LoadMuonSF(TString);
  void LoadElecSF(TString);
  void LoadElMuSF(TString);
  
  float GetMuonIDSF  (float pt, float eta) const { // binned in eta, pt
    return fMuonSF->GetBinContent(fMuonSF->FindBin(eta, pt));
  }
  float GetElecIDSF  (float pt, float eta) const { // binned in eta, pt
    return fElecSF->GetBinContent(fElecSF->FindBin(eta, pt));
  }
  float GetDoubleMuSF(float eta1, float eta2) const { // binned in eta1, eta2
    return fDoubleMuSF->GetBinContent(fDoubleMuSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }
  float GetDoubleElSF(float eta1, float eta2) const { // binned in eta1, eta2
    return fDoubleElSF->GetBinContent(fDoubleElSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));;
  }
  float GetMuEGSF    (float eta1, float eta2) const { // binned in eta1, eta2
    return fMuEGSF->GetBinContent(fMuEGSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));;
  }
 protected:
  TH2F* GetHistogramFromFile(TString, TString, TString) const;
  
 private:
  TH2F *fMuonSF;
  TH2F *fElecSF;

  TH2F *fDoubleMuSF;
  TH2F *fDoubleElSF;
  TH2F *fMuEGSF;
};
#endif
