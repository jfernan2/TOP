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
  iTopD,
  iNCUTS
};
TString sCut[iNCUTS] = {"dilepton", "ZVeto", "MET", "2jets", "1btag","TopD"};
enum gSystFlag{
  Norm,
  JESUp,
  JESDown,
  JER,
  BtagUp,
  BtagDown,
  LESUp,
  LESDown,
  LepUp,
  LepDown,
  TrigUp,
  TrigDown,
  METUp,
  METDown,
  PUUp,
  PUDown,
  TopPtUp,
  TopPtDown,
  gNSYST
};
TString SystName[gNSYST] = {
  "Normal",
  "JESUp",
  "JESDown",
  "JER",
  "BtagUp",
  "BtagDown",
  "LESUp",
  "LESDown",
  "LepUp",
  "LepDown",
  "TrigUp",
  "TrigDown",
  "METUp",
  "METDown",
  "PUUp",
  "PUDown",
  "TopPtUp",
  "TopPtDown"
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
  lepton(){};
  lepton(TLorentzVector vec, int ch, int ty, int ind){
    p = vec;
    charge = ch;
    type = ty;
    index = ind;
  };
  TLorentzVector p;
  int charge;
  int type; // -1(unknown), 0(mu), 1(ele)                                                                                                                                                                                      
  int index;
};
 
struct TLRatios{
  TH2F *fntight;
  TH2F *fnloose;
  TH2F *pntight;
  TH2F *pnloose;
   
  TH1F *fntight_nv;
  TH1F *fnloose_nv;
  TH1F *pntight_nv;
  TH1F *pnloose_nv;
   
  TEfficiency *fratio_pt;
  TEfficiency *pratio_pt;
  TEfficiency *fratio_eta;
  TEfficiency *pratio_eta;
  TEfficiency *fratio_nv;
  TEfficiency *pratio_nv;
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
  virtual void InitialiseTLRatios();
  virtual void InitialiseTree();
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
  void WriteTLRatios();
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
  bool PassSingleMuTrigger();
  bool PassSingleElTrigger();
  bool PassesJetPtdPhiCut();
  bool PassesZVeto();
  bool PassesNJetsCut();
  bool PassesMETCut();
  bool PassesNBtagCut();
  bool PassesMllVeto();
  bool Passes3rdLeptonVeto();
  bool PassesMuonEta2p1(gChannel);
  bool PassesTopDCut();
  
  int   getNTightMuons();
  int   getNTightElectrons();
  int   getNMuons();
  int   getNElectrons();
  int   getNJets();
  int   getNBTags();
  int   getLeadingJet();
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
  float getSF(gChannel, int, int);
  float getTopPtSF();
  float getTopD();
  float getDeltaPhillJet();

  // Lepton selection methods
  bool IsGoodMuon(unsigned int,float ptcut=20.);
  bool IsVetoMuon(unsigned int);
  bool IsLooseMuon(unsigned int);
  bool IsTightMuon(unsigned int);
  float getMuonIso(int);
  bool IsGoodElectron(unsigned int,float ptcut=20.);
  bool IsVetoElectron(unsigned int);  
  bool IsLooseElectron(unsigned int);  
  bool IsMVAIDElectron(unsigned int);
  bool IsTightElectron(unsigned int);
  float getElecIso(int);
  float getEACorrection(float);
  std::vector<lepton> SortLeptonsByPt(std::vector<lepton>&);
  
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
  int  HasLooseMuons(int&, int&);
  int  HasLooseMuons();
  int  HasLooseElectrons(int&, int&);
  int  HasLooseElectrons();
  bool IsZMuMuEvent(int&, int&);
  bool IsSigSupMuEvent(int&);
  bool IsZElElEvent(int&, int&);
  bool IsSigSupElEvent(int&);
  bool IsTTbarElMuEvent(int&, int&);
  bool IsTTbarElElEvent(int&, int&);
  bool IsTTbarMuMuEvent(int&, int&);
  int  IsDileptonEvent(int&, int&);
  bool IsMuMuEvent(int&, int&);
  bool IsElMuEvent(int&, int&);
  bool IsElElEvent(int&, int&);
  ///////////////////////////////////////////////////////////////////////////// 
  //    Filling Methods
  ///////////////////////////////////////////////////////////////////////////// 
  //  void FillTree(const TLorentzVector& lepton1,const TLorentzVector& lepton2, UInt_t iChannel);
  void FillAnalysisTree(int);
  void FillYieldsHistograms(gChannel, iCut, gSystFlag);
  void FillYields(gSystFlag sys=Norm);
  void FillTLRatios();
  void FillDYHistograms();
  void FillKinematicHistos(gChannel,iCut);
  
  ///////////////////////////////////////////////////////////////////////////// 
  //    Set/Reset methods
  ///////////////////////////////////////////////////////////////////////////// 
  void SetHypLepton1(int, gChannel);
  void SetHypLepton2(int, gChannel);

  void SetDataMembers();
  void ResetDataMembers();
  void ResetHypLeptons();
  void ResetAnalysisTree();
  
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
  Bool_t  gDoTLRatios;
  Bool_t  gIsData;
  Bool_t  gUseCSVM;

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
  int nGoodVertex;

  // HISTOGRAMS
  TLRatios tlratios[2];
  
  //++ Yields
  TH1F* fHDummy;
  TH1F* fHyields  [gNCHANNELS][gNSYST];
  TH1F* fHSSyields[gNCHANNELS];
  
  TH2F* fHDY_InvMassVsNPV   [gNCHANNELS][iNCUTS];
  TH2F* fHDY_InvMassVsMET   [gNCHANNELS][iNCUTS];
  TH2F* fHDY_InvMassVsNjets [gNCHANNELS][iNCUTS];
  TH2F* fHDY_InvMassVsNbtags[gNCHANNELS][iNCUTS];
  TH1F* fHDY_InvMass        [gNCHANNELS][iNCUTS];

  //++ Kinematic  
  TH1F* fHMET[gNCHANNELS][iNCUTS];       
  TH1F* fHInvMass[gNCHANNELS][iNCUTS];   
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

//SANTI
//SANTI  Double_t                    missingEt;
//SANTI  Double_t                    dileptonInvMass;
//SANTI  Double_t                    var;
//SANTI  Int_t                       iJet1;
//SANTI  Int_t                       iJet2;
//SANTI  Int_t                       iBtagJet1;
//SANTI  Int_t                       iBtagJet2;
//SANTI  UInt_t                      nJets;
//SANTI  UInt_t                      nJetsBtag;
//SANTI  Double_t                    HT;
  UInt_t                      nGenElec;
  UInt_t                      nGenMuon;
  UInt_t                      nGenTau;
  UInt_t                      nGenLepton;
  UInt_t                      nTauElec;
  UInt_t                      nTauMuon;
  UInt_t                      nSGenMuon;
  UInt_t                      nSGenElec;
  
  std::vector<float> JetEt;
  std::vector<float> MuPx;
  std::vector<float> MuPy;
  std::vector<float> ElPx;
  std::vector<float> ElPy;
  float MET;
  float MET_Phi;

  ///////////////////////////////////////
  //    OUTPUT TREE
  ///////////////////////////////////////
  TTree      *AnalysisTree;

  // Branches
  Long64_t    TEvent;
  int         TLumi;
  int         TRun;
  std::string TSName;
  int         TSType;

  float  TWeight;
  int    TChannel;
  int    TSystFlag;
  int    TTLCat;
  int    TNPV;
  float  TMET;
  float  TMET_Phi;
  
  int   TNTMus;
  int   TNTEls;
  int   TNMus;
  int   TNEls;
  float TInvMass;
  float TLep0Pt;
  float TLep0Eta;
  float TLep0Phi;
  float TLep0Ch;
  float TLep0Iso;
  float TLep0Flav;
  float TLep1Pt;
  float TLep1Eta;
  float TLep1Phi;
  float TLep1Ch;
  float TLep1Iso;
  float TLep1Flav;

  int   TNJets;     
  int   TNbJets;     
  int   TNbJetsMed;     
  int   TNJetsBtag; 
  float THT;
  float TBtagJet0;  
  float TBtagJet1;  

  float TJet0Px;    
  float TJet0Py;    
  float TJet0Pz;    
  float TJet0Et;    
  float TJet0E;    
  float TJet1Px;    
  float TJet1Py;    
  float TJet1Pz;    
  float TJet1Et;    
  float TJet1E;    

  float TBtagJet0Px;    
  float TBtagJet0Py;    
  float TBtagJet0Pz;    
  float TBtagJet0Et;    
  float TBtagJet0E;    
  float TBtagJet1Px;    
  float TBtagJet1Py;    
  float TBtagJet1Pz;    
  float TBtagJet1Et;    
  float TBtagJet1E;    

  float TtPx;
  float TtPy;
  float TtPz;
  float TtE;

  float TtbarPx;
  float TtbarPy;
  float TtbarPz;
  float TtbarE;


  ClassDef(TreeAnalysisTop,0);
};


#endif

/*  LocalWords:  ifndef
 */
