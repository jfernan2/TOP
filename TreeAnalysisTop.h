#ifndef TREEANALYSIS_TOP_H
#define TREEANALYSIS_TOP_H 1

#include "PAFAnalysis.h"
#include "TCounterUI.h"
#include "PUWeight.h"
//#include "GlobalVariables.h"
#include "LeptonSF.h"
#include "BTagSFUtil.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"

#include <fstream>
#include <iostream>
#include <vector>

enum gChannel{
  channels_begin,
  Muon = channels_begin,
  Elec,
  ElMu,
  gNCHANNELS,
};
TString gChanLabel[gNCHANNELS] = {"Muon","Elec","ElMu"};
enum gFPSwitch{
  SigSup,
  ZDecay,
  Sig
};
enum iCut{
  iDilepton, 
  iZVeto, 
  iMET, 
  i2jets, 
  i1btag, 
  iExact1btag,
  iExact2btag,
  iNCUTS
};
TString sCut[iNCUTS] = {"dilepton", "ZVeto", "MET", "2jets", "1btag","Exact1btag","Exact2btag"};
enum gSystFlag{
  Norm,
  BtagUp,
  BtagDown,
  MisTagUp,
  MisTagDown,
  JESUp,
  JESDown,
  JER,
  LESUp,
  LESDown,
  /*  LepUp,
      LepDown,
      TrigUp,
      TrigDown,
  */
  PUUp,
  PUDown,
  TopPt,
  gNSYST
};
TString SystName[gNSYST] = {
  "Normal",
  "BtagUp",
  "BtagDown",
  "MisTagUp",
  "MisTagDown",
  "JESUp",
  "JESDown",
  "JER",
  "LESUp",
  "LESDown",
  /*  "LepUp",
  "LepDown",
  "TrigUp",
  "TrigDown",*/
//  "METUp",
//  "METDown",
  "PUUp",
  "PUDown",
  "TopPt",
};
enum FakeSource{
  HF_mu,
  Other_mu,
  HF_el,
  Conv_el,
  Other_el,
  RightSign,
  WrongSign,
  gNFAKESOURCE
};
class lepton{
 public:
  lepton(){}
  lepton(const lepton &l): p(l.p), charge(l.charge), type(l.type), index(l.index){ };
  lepton(TLorentzVector vec, int ch, int ty, int ind){
    p = vec;
    charge = ch;
    type = ty;
    index = ind;
  }
  TLorentzVector p;
  int charge;
  int type; // -1(unknown), 0(mu), 1(ele)                                                                                                                                                                                      
  int index;
};

class jet{
 public:
  jet(){};
  jet(TLorentzVector vec, bool btag, int ind){
    p = vec;
    isbtag = btag;
    index = ind;
  };
  TLorentzVector p;
  bool isbtag;
  int index;
};
 
 
const int gNMuFPtBins = 6;
const int gNMuPPtbins = 10;
const int gNMuEtabins = 5;
const int gNElFPtBins = 8;
const int gNElPPtbins = 10;
const int gNElEtabins = 5;
const int gNElCMIdbins = 2;
 
// Muon Binning
double gMuFPtBins[gNMuFPtBins+1] = {20., 25., 30., 35., 40., 50., 60.}; 
double gMuPPtbins[gNMuPPtbins+1] = {20., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100.};
double gMuEtabins[gNMuEtabins+1] = {0., 0.5, 1.0, 1.479, 2.0, 2.5};
 
// Electron Binning //////////////////////////////////////////////////////////////
double gElFPtBins[gNElFPtBins+1]   = {20., 25., 30., 40., 50., 60., 70., 80., 100.};
double gElPPtbins[gNElPPtbins+1]   = {20., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100.}; 
double gElEtabins[gNElEtabins+1]   = {0., 0.5, 1.0, 1.479, 2.0, 2.5};
 
int     getNFPtBins(gChannel chan){ // fake ratios
  if(chan == Muon || chan == ElMu) return gNMuFPtBins;
  if(chan == Elec)                 return gNElFPtBins;
  else return -99;
};
double *getFPtBins (gChannel chan){
  if(chan == Muon || chan == ElMu) return gMuFPtBins;
  else                             return gElFPtBins;
  //  if(chan == Elec)                 return gElFPtBins;
  //  else return *-999;
};
int     getNPPtBins(gChannel chan){
  if(chan == Muon || chan == ElMu) return gNMuPPtbins;
  if(chan == Elec)                 return gNElPPtbins;
  else return -99;
};
double *getPPtBins (gChannel chan){
  if(chan == Muon || chan == ElMu) return gMuPPtbins;
  else                             return gElPPtbins;
  //  if(chan == Elec)                 return gElPPtbins;
  //  else return -99;
};
int     getNEtaBins(gChannel chan){
  if(chan == Muon || chan == ElMu) return gNMuEtabins;
  if(chan == Elec)                 return gNElEtabins;
  else return -99;
};
double *getEtaBins (gChannel chan){
  if(chan == Muon || chan == ElMu) return gMuEtabins;
  else                             return gElEtabins;
  //  if(chan == Elec)                 return gElEtabins;
  //  else return *-99.;
}; 

class TreeAnalysisTop: public PAFAnalysis {
 public:
  TreeAnalysisTop(TTree *tree=0);
  virtual ~TreeAnalysisTop() {}
  
  virtual void Initialise();
  virtual void InitialiseYieldsHistos();
  virtual void InitialiseKinematicHistos();
  virtual void InitialiseDYHistos();
  virtual void InitialiseGenHistos();
  virtual void InitialiseSystematicHistos();
  virtual void InsideLoop();
  virtual void SetDataMembersAtTermination();
  virtual void Summary();

  // Saving Histograms.
  void WriteHistos();
  void WriteValidationsHistos(){};

  // Counters
  //TCounterUI *nEvents;

  // My member functions
  //----------------------------------------------------------------------------
  void  GetParameters();
  Int_t SelectedVertexIndex();
  
  bool PassTriggerMuMu();
  bool PassTriggerEE();
  bool PassTriggerEMu();
  bool PassesZVeto();
  bool PassesNJetsCut();
  bool PassesMETCut();
  bool PassesNBtagCut();
  bool PassesMllVeto();
  bool Passes3rdLeptonVeto();
  bool PassesMuonEta2p1(gChannel);
  bool PassesTopDCut();

  int   getNJets();
  int   getNBTags();
  int   getLeadingJetbTag();
  float getDRClosestJet(TLorentzVector);
  float getDPhiClosestJet(TLorentzVector);
  //  int   getNBTagsMed();
  void  setMET(float);
  float getMET();
  float getMETPhi();
  float getMT(int, gChannel);
  float getHT();
  float getJetPtIndex(unsigned int);
  float getBtagJetPtIndex(unsigned int);
  float getErrPt(float,float);
  float getJERScale(int);
  float getSF(gChannel);
  float getLeptonError(gChannel);
  float getTriggerError(gChannel);
  float getTopPtSF();
  float getTopD();
  float getDeltaPhillJet();

  // Lepton selection methods
  int  getSelectedLeptons();
  bool IsVetoMuon(unsigned int, float ptcut=20.);
  bool IsTightMuon(unsigned int, float ptcut=20.);
  float getMuonIso(int);
  bool IsVetoElectron(unsigned int,float ptcut=20.);
  bool IsMVAIDElectron(unsigned int);
  bool IsTightElectron(unsigned int,float ptcut=20.);
  float getElecIso(int);
  float getEACorrection(float);
  std::vector<lepton> SortLeptonsByPt(std::vector<lepton>&);
  
  int getSelectedJets();
  bool IsGoodJet(unsigned int, float ptcut=30.);
  std::vector<int> CleanedJetIndices(float);
  void SmearJetPts(int);
  void propagateMET(TLorentzVector,TLorentzVector);
  void ScaleMET(int);
  void ScaleLeptons(int);
  
  //void GetGenMuon();
  //void GetGenElec();
  void SelectedGenLepton();
  
  ///////////////////////////////////////////////////////////////////////////// 
  //    Selecting Methods
  ///////////////////////////////////////////////////////////////////////////// 
  int  IsDileptonEvent();
  bool IsMuMuEvent();
  bool IsElMuEvent();
  bool IsElElEvent();
  ///////////////////////////////////////////////////////////////////////////// 
  //    Filling Methods
  ///////////////////////////////////////////////////////////////////////////// 
  void FillYieldsHistograms(gChannel, iCut, gSystFlag);
  void FillYields(gSystFlag sys=Norm);
  void FillDYHistograms();
  void FillKinematicHistos(gChannel,iCut);
  
  ///////////////////////////////////////////////////////////////////////////// 
  //    Set/Reset methods
  ///////////////////////////////////////////////////////////////////////////// 
  void SetOriginalObjects();
  void ResetOriginalObjects();
  void SetEventObjects();
  void ResetHypLeptons();
  
 protected:
  
  // Input parameters
  //----------------------------------------------------------------------------
  TString gSampleName;
  TString gfileSuffix;
  Float_t gWeight;
  Float_t gLumiForPU;
  Float_t gTotalLumi;
  Int_t   gSysSource;
  Int_t   gSysDirection;
  Bool_t  gDoSystStudies;
  Bool_t  gIsData;
  Bool_t  gUseCSVM;
  Bool_t  gDoSF;
  Bool_t  gDoDF;
  Float_t gStopMass;
  Float_t gLspMass;

  PUWeight *fPUWeight;     //The PU weight utility
  PUWeight *fPUWeightUp;   //The PU weight utility
  PUWeight *fPUWeightDown; //The PU weight utility
  BTagSFUtil *fBTagSF;     //The BTag SF utility 
  LeptonSF *fLeptonSF;
  TRandom3 *fRand3;
  
  // EventWeight
  //----------------------------------------------------------------------------
  float EventWeight;
  float PUSF;
  bool  fChargeSwitch;

  //////////////////////////////////////////////////////////////////////////////
  //               Data members
  //////////////////////////////////////////////////////////////////////////////
  // HISTOGRAMS
  
  //++ Yields
  TH1F* fHDummy;
  TH1F* hWeight;
  TH1F* fHyields  [gNCHANNELS][gNSYST];
  TH1F* fHSSyields[gNCHANNELS][gNSYST];
  TH1F* fHTopPtWeight;
  TH1F* fHpdfWeightSum;
  TH1F* fHpdfWeight;
  TH1F* fHLepSys[gNCHANNELS][iNCUTS];
  TH1F* fHTrigSys[gNCHANNELS][iNCUTS];

  TH2F* fHDY_InvMassVsNPV   [gNCHANNELS][iNCUTS];
  TH2F* fHDY_InvMassVsMET   [gNCHANNELS][iNCUTS];
  TH2F* fHDY_InvMassVsNjets [gNCHANNELS][iNCUTS];
  TH2F* fHDY_InvMassVsNbtags[gNCHANNELS][iNCUTS];
  TH1F* fHDY_InvMass        [gNCHANNELS][iNCUTS];
  
  //++ Origin Histos
//  TH2F* fHSSOrigins[gNCHANNELS][iNCUTS];
//  TH2F* fHOrigins[gNCHANNELS][iNCUTS];
  
  //++ Kinematic  
  TH1F* fHMET[gNCHANNELS][iNCUTS];       
  TH1F* fHDiLepPt[gNCHANNELS][iNCUTS];   
  TH1F* fHLep0Pt[gNCHANNELS][iNCUTS];    
  TH1F* fHLep1Pt[gNCHANNELS][iNCUTS];    
  TH1F* fHDelLepPhi[gNCHANNELS][iNCUTS]; 
  TH1F* fHNJets[gNCHANNELS][iNCUTS];     
  TH1F* fHHT[gNCHANNELS][iNCUTS];        
  TH1F* fHNBtagJets[gNCHANNELS][iNCUTS]; 
  TH1F* fHJet0Pt[gNCHANNELS][iNCUTS];    
  TH1F* fHJet1Pt[gNCHANNELS][iNCUTS];    
  TH1F* fHBtagJet0Pt[gNCHANNELS][iNCUTS];
  
  TH1F* fHInvMass[gNCHANNELS][iNCUTS][gNSYST];   
  TH1F* fHSSInvMass[gNCHANNELS][iNCUTS][gNSYST];   
  TH1F* fHNBtagsNJets[gNCHANNELS][iNCUTS][gNSYST]; 
  TH1F* fHSSNBtagsNJets[gNCHANNELS][iNCUTS][gNSYST]; 
  TH1F* fHCSVTag[gNCHANNELS][iNCUTS]; 
  TH1F* fHTopD[gNCHANNELS][iNCUTS];
  TH1F* fHDelPhillJet[gNCHANNELS][iNCUTS];

  TH1F* fHDRLep[gNCHANNELS][iNCUTS];
  TH1F* fHDRLep0Jet[gNCHANNELS][iNCUTS];
  TH1F* fHDPhiLep0Jet[gNCHANNELS][iNCUTS];
  TH1F* fHLep0Iso[gNCHANNELS][iNCUTS];
  TH1F* fHDRLep1Jet[gNCHANNELS][iNCUTS];
  TH1F* fHDPhiLep1Jet[gNCHANNELS][iNCUTS];
  TH1F* fHLep1Iso[gNCHANNELS][iNCUTS];
  
  /// STOP
  //TH1F* fHAbsDelPhiLep[gNCHANNELS][iNCUTS];
  TH1F* fHminDelRJetsLeps[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHSSminDelRJetsLeps[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHdelPhi2LeadJets[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHSSdelPhi2LeadJets[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHAbsDelPhiLeps[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHSSAbsDelPhiLeps[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHStopMass[gNCHANNELS][iNCUTS];
  TH1F* fHChi0Mass[gNCHANNELS][iNCUTS];
  TH2F* fHChi0StopMass[gNCHANNELS][iNCUTS];
  TH1F* fHvertices[gNCHANNELS][iNCUTS];
  
  //++ Gen Info
  TH1F* fHDeltaRLepJet[gNCHANNELS-1];

  lepton fHypLepton1;
  lepton fHypLepton2;
  
  std::vector<Double_t>       Gen_Muon_Charge;
  std::vector<Double_t>       Gen_Elec_Charge;
  std::vector<TLorentzVector> Gen_Muon;
  std::vector<TLorentzVector> Gen_Elec;
 
  std::vector<Int_t>          NGen_Jet;
  std::vector<Int_t>          NGen_b;
  
  std::vector<Double_t>       PtGen_Jet;
  std::vector<Double_t>       PtGen_b;

  Int_t nGenElec;
  Int_t nGenMuon;
  Int_t nGenTau;
  Int_t nGenLepton;
  Int_t nTauElec;
  Int_t nTauMuon;
  Int_t nSGenMuon;
  Int_t nSGenElec;
  
  Int_t nGoodVertex;
  Int_t nBtags;
  Int_t nJets;
  Int_t nMuon;
  Int_t nElec;
  Int_t nLeptons;

  ///// OBJECTS
  std::vector<lepton> Lepton;
  std::vector<jet>    Jet;
  //  std::vector<jet>    Jet15;
  
  std::vector<float> JetEt;
  std::vector<float> JetPhi;
  std::vector<float> MuPx;
  std::vector<float> MuPy;
  std::vector<float> ElPx;
  std::vector<float> ElPy;
  float MET;
  float MET_Phi;
  
  ClassDef(TreeAnalysisTop,0);
};


#endif

/*  LocalWords:  ifndef
 */
