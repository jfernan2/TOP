#include "TreeAnalysisTop.h"

#include <fstream>
#include <iostream>
#include <math.h>

#if !defined(__CINT__)
ClassImp(TreeAnalysisTop);
#endif

const float gJetEtCut = 30.;

TreeAnalysisTop::TreeAnalysisTop(TTree* tree) : PAFAnalysis(tree) {}
//------------------------------------------------------------------------------
// Initialise
//------------------------------------------------------------------------------
void TreeAnalysisTop::Initialise() {
#ifdef DEBUG
  cout << "Initialise(): Enter" << endl;
#endif 
  GetParameters();

  fHDummy = CreateH1F("fHDummy","",1,0,1);
  fHDummy->TH1::SetDefaultSumw2();

  InitialiseYieldsHistos();
#ifdef __ISFR
  InitialiseTLRatios();
#endif
  InitialiseKinematicHistos();

#ifdef __ISMC
  InitialiseGenHistos();
#endif  
  fHTopPtWeight = CreateH1F("H_TopPtWeight","TopPt Weight",100, 0, 2);
  
  if (gSampleName == "DoubleMu"        ||       
      gSampleName == "DoubleElectron"  || 
      gSampleName == "MuEG"            ||       
      gSampleName == "DYJets_Madgraph" ||
      gSampleName == "ZJets_Madgraph")  {
    InitialiseDYHistos();
  }
  InitialiseSystematicHistos();
  //  InitialiseTree();

  /***********
  PU Reweight
  ***********/
  fPUWeight     = new PUWeight(gLumiForPU,Summer12_53X,"2012");
#ifdef __ISMC
  fPUWeightUp   = new PUWeight(18494.9,   Summer12_53X,"2012"); //  18494.9  (5% down)
  fPUWeightDown = new PUWeight(20441.7,   Summer12_53X,"2012"); //  20441.7  (5% up  )
#endif

  if (gUseCSVM) fBTagSF   = new BTagSFUtil("CSVM","ABCD");//ReReco
  else          fBTagSF   = new BTagSFUtil("CSVT","ABCD");//ReReco 
  
  fLeptonSF = new LeptonSF();

  fRand3 = new TRandom3(50);

  // No systematics activaded...
  gSysSource = Norm;
  fChargeSwitch = false;

#ifdef DEBUG
  cout << "Initialise(): Exit" << endl;
#endif 
}
void TreeAnalysisTop::InitialiseTLRatios(){
#ifdef DEBUG
  cout << "InitialiseTLRatios(): Enter" << endl;
#endif 
  
  for(size_t l = 0; l < 2; ++l){
    gChannel c = channels_begin;
    if      (l == 0) c = Muon;
    else if (l == 1) c = Elec;
    
    TString rootname = gChanLabel[c];
    tlratios[l].fntight  = CreateH2F(rootname + "_fNTight",  "fNTight",  getNFPtBins(c), getFPtBins(c), getNEtaBins(c), getEtaBins(c)); tlratios[l].fntight ->Sumw2();
    tlratios[l].fnloose  = CreateH2F(rootname + "_fNLoose",  "fNLoose",  getNFPtBins(c), getFPtBins(c), getNEtaBins(c), getEtaBins(c)); tlratios[l].fnloose ->Sumw2();
    tlratios[l].pntight  = CreateH2F(rootname + "_pNTight",  "pNTight",  getNPPtBins(c), getPPtBins(c), getNEtaBins(c), getEtaBins(c)); tlratios[l].pntight ->Sumw2();
    tlratios[l].pnloose  = CreateH2F(rootname + "_pNLoose",  "pNLoose",  getNPPtBins(c), getPPtBins(c), getNEtaBins(c), getEtaBins(c)); tlratios[l].pnloose ->Sumw2();
    
    tlratios[l].fntight_nv  = CreateH1F(rootname + "_fNTight_nv",  "fNTight_nv", 18, 0., 36.); tlratios[l].fntight_nv ->Sumw2();
    tlratios[l].fnloose_nv  = CreateH1F(rootname + "_fNLoose_nv",  "fNLoose_nv", 18, 0., 36.); tlratios[l].fnloose_nv ->Sumw2();
    tlratios[l].pntight_nv  = CreateH1F(rootname + "_pNTight_nv",  "pNTight_nv", 18, 0., 36.); tlratios[l].pntight_nv ->Sumw2();
    tlratios[l].pnloose_nv  = CreateH1F(rootname + "_pNLoose_nv",  "pNLoose_nv", 18, 0., 36.); tlratios[l].pnloose_nv ->Sumw2();
       
    tlratios[l].fratio_pt  = new TEfficiency(rootname + "_fRatio_pt",  "fRatio_pt",  getNFPtBins(c), getFPtBins(c));
    tlratios[l].fratio_eta = new TEfficiency(rootname + "_fRatio_eta", "fRatio_eta", getNEtaBins(c), getEtaBins(c));
    tlratios[l].pratio_pt  = new TEfficiency(rootname + "_pRatio_pt",  "pRatio_pt",  getNPPtBins(c), getPPtBins(c));
    tlratios[l].pratio_eta = new TEfficiency(rootname + "_pRatio_eta", "pRatio_eta", getNEtaBins(c), getEtaBins(c));
    tlratios[l].fratio_nv  = new TEfficiency(rootname + "_fRatio_nv",  "fRatio_nv",  18, 0., 36.);
    tlratios[l].pratio_nv  = new TEfficiency(rootname + "_pRatio_nv",  "pRatio_nv",  18, 0., 36.);

    fOutput->Add(tlratios[l].fratio_pt  );
    fOutput->Add(tlratios[l].fratio_eta );
    fOutput->Add(tlratios[l].pratio_pt  );
    fOutput->Add(tlratios[l].pratio_eta );
    fOutput->Add(tlratios[l].fratio_nv  );
    fOutput->Add(tlratios[l].pratio_nv  );
  }
  
#ifdef DEBUG
  cout << "InitialiseTLRatios(): Exit" << endl;
#endif 
}
void TreeAnalysisTop::InitialiseGenHistos(){
  fHDeltaRLepJet[Muon] = CreateH1F("H_DeltaRLepJet_"+gChanLabel[Muon],"",1000,0.,5.);
  fHDeltaRLepJet[Elec] = CreateH1F("H_DeltaRLepJet_"+gChanLabel[Elec],"",1000,0.,5.);

  
}
void TreeAnalysisTop::InitialiseDYHistos(){
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    for (size_t cut=0; cut<iNCUTS; cut++){
      TString name = "_"+gChanLabel[ch]+"_"+sCut[cut];
      fHDY_InvMassVsNPV   [ch][cut] = CreateH2F("H_DY_InvMassVsNPV"   +name,"",50, -0.5, 49.5, 200, 0, 200);
      fHDY_InvMassVsMET   [ch][cut] = CreateH2F("H_DY_InvMassVsMET"   +name,"",200, 0  , 200 , 200, 0, 200);
      fHDY_InvMassVsNjets [ch][cut] = CreateH2F("H_DY_InvMassVsNjets" +name,"",10, -0.5,  9.5, 200, 0, 200);
      fHDY_InvMassVsNbtags[ch][cut] = CreateH2F("H_DY_InvMassVsNbtags"+name,"",10, -0.5,  9.5, 200, 0, 200);
      fHDY_InvMass        [ch][cut] = CreateH1F("H_DY_InvMass"        +name,"",                200, 0, 200);
    }
  }
}
void TreeAnalysisTop::InitialiseYieldsHistos(){
  
  //++ Yields histograms
  fHyields[Muon][Norm]    = CreateH1F("H_Yields_"+gChanLabel[Muon],"", iNCUTS, -0.5, iNCUTS-0.5); 
  fHyields[Elec][Norm]   = CreateH1F("H_Yields_"+gChanLabel[Elec],"", iNCUTS, -0.5, iNCUTS-0.5);
  fHyields[ElMu][Norm]   = CreateH1F("H_Yields_"+gChanLabel[ElMu],"", iNCUTS, -0.5, iNCUTS-0.5);
  fHSSyields[Muon][Norm] = CreateH1F("H_SSYields_"+gChanLabel[Muon],"", iNCUTS, -0.5, iNCUTS-0.5); 
  fHSSyields[Elec][Norm] = CreateH1F("H_SSYields_"+gChanLabel[Elec],"", iNCUTS, -0.5, iNCUTS-0.5);
  fHSSyields[ElMu][Norm] = CreateH1F("H_SSYields_"+gChanLabel[ElMu],"", iNCUTS, -0.5, iNCUTS-0.5);
  
  if (gDoSystStudies){
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      for (size_t sys=1; sys<gNSYST; sys++){
	fHyields[chan][sys]   = CreateH1F("H_Yields_"+gChanLabel[chan]+"_"+SystName[sys],"",iNCUTS,-0.5,iNCUTS-0.5);
	fHSSyields[chan][sys] = CreateH1F("H_SSYields_"+gChanLabel[chan]+"_"+SystName[sys],"", iNCUTS, -0.5, iNCUTS-0.5);
      }
    }
  }
  
  for (size_t chan=0; chan<gNCHANNELS; chan++){
    for (size_t cut=0; cut<iNCUTS; cut++){
      fHLepSys [chan][cut] = CreateH1F("H_LepSys_" +gChanLabel[chan]+"_"+sCut[cut],"LepSys" , 400, 0, 0.04);
      fHTrigSys[chan][cut] = CreateH1F("H_TrigSys_"+gChanLabel[chan]+"_"+sCut[cut],"TrigSys", 400, 0, 0.04);
    }
  }
}
void TreeAnalysisTop::InitialiseKinematicHistos(){
  //++ Kinematic histograms
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    for (size_t cut=0; cut<iNCUTS; cut++){
      fHMET[ch][cut]         = CreateH1F("H_MET_"        +gChanLabel[ch]+"_"+sCut[cut],"MET"       ,  5000,0,500);
      fHDiLepPt[ch][cut]     = CreateH1F("H_DiLepPt_"    +gChanLabel[ch]+"_"+sCut[cut],"DiLepPt"   , 1800,20,200); 
      fHLep0Pt[ch][cut]      = CreateH1F("H_Lep0Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Lep0Pt"    , 1800,20,200);
      fHLep1Pt[ch][cut]      = CreateH1F("H_Lep1Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Lep1Pt"    , 1800,20,200);
      fHDelLepPhi[ch][cut]   = CreateH1F("H_DelLepPhi_"  +gChanLabel[ch]+"_"+sCut[cut],"DelLepPhi" , 100,-3.2, 3.2);
      fHNJets[ch][cut]       = CreateH1F("H_NJets_"      +gChanLabel[ch]+"_"+sCut[cut],"NJets"     , 10 ,-0.5, 9.5);
      fHHT[ch][cut]          = CreateH1F("H_HT_"         +gChanLabel[ch]+"_"+sCut[cut],"HT"        , 4700,30,500);
      fHNBtagJets[ch][cut]   = CreateH1F("H_NBtagJets_"  +gChanLabel[ch]+"_"+sCut[cut],"NBtagJets" , 10 ,-0.5, 9.5);
      fHJet0Pt[ch][cut]      = CreateH1F("H_Jet0Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Jet0Pt"    , 2700,30,300);
      fHJet1Pt[ch][cut]      = CreateH1F("H_Jet1Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Jet1Pt"    , 2700,30,300);
      fHBtagJet0Pt[ch][cut]  = CreateH1F("H_BtagJet0Pt_" +gChanLabel[ch]+"_"+sCut[cut],"BtagJet0Pt", 2700,30,300);
      
      fHInvMass[ch][cut][0]       = CreateH1F("H_InvMass_"    +gChanLabel[ch]+"_"+sCut[cut],"InvMass"   , 300, 0., 300.);
      fHSSInvMass[ch][cut][0]     = CreateH1F("HSS_InvMass_"  +gChanLabel[ch]+"_"+sCut[cut],"InvMass"   , 300, 0., 300.);
      
      fHNBtagsNJets[ch][cut][0]   = CreateH1F("H_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut]  ,"NBtagsNJets"   ,15 , -0.5, 14.5);
      fHSSNBtagsNJets[ch][cut][0] = CreateH1F("HSS_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut],"SS_NBtagsNJets",15 , -0.5, 14.5);

      // other variables 
      fHCSVTag[ch][cut]      = CreateH1F("H_CSVTag_"+gChanLabel[ch]+"_"+sCut[cut],"NBtagsNJets", 1000, 0.0, 1.0);
      fHTopD[ch][cut]        = CreateH1F("H_TopD_"       +gChanLabel[ch]+"_"+sCut[cut],"TopDiscriminator",1000,0.0,1.0);
      fHDelPhillJet[ch][cut] = CreateH1F("H_DelPhillJet_"+gChanLabel[ch]+"_"+sCut[cut], "DeltaPhi", 1000,0.0, TMath::Pi());

      // Different Top / Z topologies
      fHDRLep[ch][cut]       = CreateH1F("H_DRLep_"       +gChanLabel[ch]+"_"+sCut[cut], "DeltaRLep",       1000,0.0, 5.0);
      fHDRLep0Jet[ch][cut]   = CreateH1F("H_DRLep0Jet_"   +gChanLabel[ch]+"_"+sCut[cut], "DeltaRLep0Jet",   1000,0.0, 5.0);
      fHDPhiLep0Jet[ch][cut] = CreateH1F("H_DPhiLep0Jet_" +gChanLabel[ch]+"_"+sCut[cut], "DeltaPhiLep0Jet", 1000,0.0, TMath::Pi());
      fHLep0Iso[ch][cut]     = CreateH1F("H_Lep0Iso_"     +gChanLabel[ch]+"_"+sCut[cut], "Lep0Iso",         1000,0.0, 0.5);
      fHDRLep1Jet[ch][cut]   = CreateH1F("H_DRLep1Jet_"   +gChanLabel[ch]+"_"+sCut[cut], "DeltaRLep1Jet",   1000,0.0, 5.0);
      fHDPhiLep1Jet[ch][cut] = CreateH1F("H_DPhiLep1Jet_" +gChanLabel[ch]+"_"+sCut[cut], "DeltaPhiLep1Jet", 1000,0.0, TMath::Pi());
      fHLep1Iso[ch][cut]     = CreateH1F("H_Lep1Iso_"     +gChanLabel[ch]+"_"+sCut[cut], "Lep1Iso",         1000,0.0, 0.5);


      // STOP HISTOGRAMS:
      //      fHAbsDelPhiLep[ch][cut] = CreateH1F("H_AbsDelPhiLep_"+gChanLabel[ch]+"_"+sCut[cut],"AbsDelPhiLep" , 66,0, 3.3);
      fHAbsDelPhiLep[ch][cut]= CreateH1F("H_AbsDelPhiLep_" +gChanLabel[ch]+"_"+sCut[cut],"AbsDelPhiLep" , 28,-0.2, 1.2);

#ifdef __ISSTOP
      fHStopMass[ch][cut]     = CreateH1F("H_StopMass_"     +gChanLabel[ch]+"_"+sCut[cut], "StopMass",    500, 0.0, 500);
      fHChi0Mass[ch][cut]     = CreateH1F("H_Chi0Mass_"     +gChanLabel[ch]+"_"+sCut[cut], "Chi0Mass",    500, 0.0, 500);
      fHChi0StopMass[ch][cut] = CreateH2F("H2_Chi0StopMass_"+gChanLabel[ch]+"_"+sCut[cut], "Chi0Mass vs StopMass", 34, 81.25, 506.25, 34, -18.75, 406.25);
#endif
    }
  }
}
void TreeAnalysisTop::InitialiseSystematicHistos(){
  
  TString histoname = "";
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    for (size_t cut=0; cut<iNCUTS; cut++){
      for (size_t sys=1; sys<gNSYST; sys++){
	histoname = "H_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHNBtagsNJets[ch][cut][sys] = CreateH1F(histoname,"NBtagsNJets", 15 , -0.5, 14.5);
	
	histoname = "H_InvMass_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHInvMass[ch][cut][sys]     = CreateH1F(histoname, "InvMass"   , 300, 0., 300.);

	histoname = "HSS_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHSSNBtagsNJets[ch][cut][sys] = CreateH1F(histoname,"SS_NBtagsNJets", 15 , -0.5, 14.5);
	
	histoname = "HSS_InvMass_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHSSInvMass[ch][cut][sys]     = CreateH1F(histoname,"InvMass"   , 300, 0., 300.);
	
      }
    }
  }
}
void TreeAnalysisTop::InitialiseTree(){
#ifdef DEBUG
  cout << "InitialiseTree(): Enter" << endl;
#endif 
  /********************
     Tree Branches 
  ********************/ 
  AnalysisTree = CreateTree("AnalysisTree","TopTree");
  
  // Run variables
  AnalysisTree->Branch("TEvent", &TEvent, "TEvent/L");
  AnalysisTree->Branch("TLumi",  &TLumi,  "TLumi/I");
  AnalysisTree->Branch("TRun",   &TRun,   "TRun/I");
  
  // Sample variables
  AnalysisTree->Branch("TSName", &TSName);
  AnalysisTree->Branch("TSType", &TSType, "TSType/I");
  
  // Per-Event variables
  AnalysisTree->Branch("TSystFlag", &TSystFlag, "TSystFlag/I"); 
  AnalysisTree->Branch("TWeight",   &TWeight,   "TWeight/F");
  AnalysisTree->Branch("TChannel",  &TChannel,  "TChannel/I");
  AnalysisTree->Branch("TTLCat",    &TTLCat,    "TTLCat/I");  //0-TT; 1-TL; 2-LT; 3-LL
  AnalysisTree->Branch("TNPV",      &TNPV,      "TNPV/I");
  AnalysisTree->Branch("TMET",      &TMET,      "TMET/F");
  AnalysisTree->Branch("TMET_Phi",  &TMET_Phi,  "TMET_Phi/F");
   
  // Lepton variables
  AnalysisTree->Branch("TNTMus",   &TNTMus,   "TNTMus/I");
  AnalysisTree->Branch("TNTEls",   &TNTEls,   "TNTEls/I");
  AnalysisTree->Branch("TNMus",    &TNMus,    "TNMus/I");
  AnalysisTree->Branch("TNEls",    &TNEls,    "TNEls/I");
  AnalysisTree->Branch("TInvMass", &TInvMass, "TInvMass/F");
  
  AnalysisTree->Branch("TLep0Pt",  &TLep0Pt,  "TLep0Pt/F");
  AnalysisTree->Branch("TLep0Eta", &TLep0Eta, "TLep0Eta/F");
  AnalysisTree->Branch("TLep0Phi", &TLep0Phi, "TLep0Phi/F");
  AnalysisTree->Branch("TLep0Ch",  &TLep0Ch,  "TLep0Ch/F");
  AnalysisTree->Branch("TLep0Flav",&TLep0Flav,"TLep0Flav/F");
  AnalysisTree->Branch("TLep1Pt",  &TLep1Pt,  "TLep1Pt/F");
  AnalysisTree->Branch("TLep1Eta", &TLep1Eta, "TLep1Eta/F");
  AnalysisTree->Branch("TLep1Phi", &TLep1Phi, "TLep1Phi/F");
  AnalysisTree->Branch("TLep1Ch",  &TLep1Ch,  "TLep1Ch/F");
  AnalysisTree->Branch("TLep1Flav",&TLep0Flav,"TLep1Flav/F");

  // Jet variables
  AnalysisTree->Branch("TNJets",     &TNJets,     "TNJets/F");
  AnalysisTree->Branch("TNbJets",    &TNbJets,    "TNbJets/F");
  AnalysisTree->Branch("TNbJetsMed", &TNbJetsMed, "TNbJetsMed/F");
  AnalysisTree->Branch("TNJetsBtag", &TNJetsBtag, "TNJetBtag/F");
  AnalysisTree->Branch("THT",        &THT,        "THT/F");
  
  AnalysisTree->Branch("TJet0Px", &TJet0Px, "TJet0Px/F");
  AnalysisTree->Branch("TJet0Py", &TJet0Py, "TJet0Py/F");
  AnalysisTree->Branch("TJet0Pz", &TJet0Pz, "TJet0Pz/F");
  AnalysisTree->Branch("TJet0Et", &TJet0Et, "TJet0Et/F");
  AnalysisTree->Branch("TJet0E",  &TJet0E,  "TJet0E/F");
  AnalysisTree->Branch("TJet1Px", &TJet1Px, "TJet1Px/F");
  AnalysisTree->Branch("TJet1Py", &TJet1Py, "TJet1Py/F");
  AnalysisTree->Branch("TJet1Pz", &TJet1Pz, "TJet1Pz/F");
  AnalysisTree->Branch("TJet1Et", &TJet1Et, "TJet1Et/F");
  AnalysisTree->Branch("TJet1E",  &TJet1E,  "TJet1E/F");
  
  AnalysisTree->Branch("TBtagJet0",   &TBtagJet0,   "TBtagJet0/F");
  AnalysisTree->Branch("TBtagJet1",   &TBtagJet1,   "TBtagJet1/F");
  AnalysisTree->Branch("TBtagJet0Px", &TBtagJet0Px, "TBtagJet0Px/F");
  AnalysisTree->Branch("TBtagJet0Py", &TBtagJet0Py, "TBtagJet0Py/F");
  AnalysisTree->Branch("TBtagJet0Pz", &TBtagJet0Pz, "TBtagJet0Pz/F");
  AnalysisTree->Branch("TBtagJet0Et", &TBtagJet0Et, "TBtagJet0Et/F");
  AnalysisTree->Branch("TBtagJet0E",  &TBtagJet0E,  "TBtagJet0E/F");
  AnalysisTree->Branch("TBtagJet1Px", &TBtagJet1Px, "TBtagJet1Px/F");
  AnalysisTree->Branch("TBtagJet1Py", &TBtagJet1Py, "TBtagJet1Py/F");
  AnalysisTree->Branch("TBtagJet1Pz", &TBtagJet1Pz, "TBtagJet1Pz/F");
  AnalysisTree->Branch("TBtagJet1Et", &TBtagJet1Et, "TBtagJet1Et/F");
  AnalysisTree->Branch("TBtagJet1E",  &TBtagJet1E,  "TBtagJet1E/F");

  // Top pT reweight
  AnalysisTree->Branch("TtPx",    &TtPx,    "TtPx/F");
  AnalysisTree->Branch("TtPy",    &TtPy,    "TtPy/F");
  AnalysisTree->Branch("TtPz",    &TtPz,    "TtPz/F");
  AnalysisTree->Branch("TtE" ,    &TtE ,    "TtE/F");
  AnalysisTree->Branch("TtbarPx", &TtbarPx, "TtbarPx/F");
  AnalysisTree->Branch("TtbarPy", &TtbarPy, "TtbarPy/F");
  AnalysisTree->Branch("TtbarPz", &TtbarPz, "TtbarPz/F");
  AnalysisTree->Branch("TtbarE" , &TtbarE , "TtbarE/F");
  
#ifdef DEBUG
  cout << "InitialiseTree(): Exit" << endl;
#endif 
}

//------------------------------------------------------------------------------
// InsideLoop
//------------------------------------------------------------------------------
void TreeAnalysisTop::SetDataMembers(){
  //ResetAnalysisTree();
  ResetHypLeptons();

  fChargeSwitch = false;

  nGenLepton = 0;
  nGenElec   = 0;
  nGenMuon   = 0;
  nGenTau    = 0;
  nTauElec   = 0;
  nTauMuon   = 0;

  PUSF = 1.;
  EventWeight = 1.;
  nGoodVertex = 0;

  // Reset to No Systematic calculation
  gSysSource = Norm;
  
  JetEt.clear();
  MuPx.clear();
  MuPy.clear();
  ElPx.clear();
  ElPy.clear();
  MET  = 0.;
  
  // Save original values for MET, Jets and Leptons
  for (UInt_t i=0; i<T_JetAKCHS_Et->size(); i++){    
    JetEt.push_back(T_JetAKCHS_Et->at(i));   
  }
  for (UInt_t i=0; i<T_Elec_Energy->size(); i++){    
    ElPx.push_back(T_Elec_Px->at(i));
    ElPy.push_back(T_Elec_Py->at(i));
  }
  for (UInt_t i=0; i<T_Muon_Energy->size(); i++){ 
    MuPx.push_back(T_Muon_Px->at(i)); 
    MuPy.push_back(T_Muon_Py->at(i));   
  }
  MET     = T_METPFType0IxyShift_ET;
  MET_Phi = T_METPFType0IxyShift_Phi;
}
void TreeAnalysisTop::ResetDataMembers(){
  
  // Save original values for MET, Jets and Leptons
  for (UInt_t i=0; i<T_JetAKCHS_Et->size(); i++){    
    JetEt[i] = T_JetAKCHS_Et->at(i);
  }
  for (UInt_t i=0; i<T_Elec_Energy->size(); i++){
    ElPx[i] = T_Elec_Px->at(i);
    ElPy[i] = T_Elec_Py->at(i);
  }
  for (UInt_t i=0; i<T_Muon_Energy->size(); i++){ 
    MuPx[i] = T_Muon_Px->at(i);
    MuPy[i] = T_Muon_Py->at(i);
  }
  setMET(T_METPFType0IxyShift_ET);

}
void TreeAnalysisTop::InsideLoop(){
  // Calculate PU Weight
  PUSF = 1.;
#ifdef __ISMC
  PUSF = fPUWeight->GetWeight((int)T_Event_nTruePU);//True      
#endif
  fHDummy->Fill(0.5);


#ifdef __ISSTOP
  if (gSampleName == "T2tt_150to250LSP1to100_LeptonFilter" &&  
      (abs(T_Gen_StopMass->at(0)-T_Gen_StopMass->at(1)) >  0.1  ||
       abs(T_Gen_Chi0Mass->at(0)-T_Gen_Chi0Mass->at(1)) >  0.1  || 
       abs(T_Gen_StopMass->at(0)-175.)                  >  6.25 ||
       abs(T_Gen_Chi0Mass->at(0)-1.)                    >  6.25)
      ) return;
   #endif

  // Init data members
  //----------------------------------------------------------------------------
  SetDataMembers();
  
  // Get number of generated leptons 
  //----------------------------------------------------------------------------
#ifdef __ISMC
  SelectedGenLepton();
  if (gSampleName == "TTJets_MadSpin"         && nGenLepton != 2) return;
  if (gSampleName == "TTJets_matchingup"      && nGenLepton != 2) return;
  if (gSampleName == "TTJets_matchingdown"    && nGenLepton != 2) return;
  if (gSampleName == "TTJets_scaleup"         && nGenLepton != 2) return;
  if (gSampleName == "TTJets_scaledown"       && nGenLepton != 2) return;
  if (gSampleName == "TTJetsFullLeptMGtauola" && nGenLepton != 2) return;
  

  if (gSampleName == "TTJets_MadSpin"         ||
      gSampleName == "TTJets_matchingup"      ||
      gSampleName == "TTJets_matchingdown"    ||
      gSampleName == "TTJets_scaleup"         ||
      gSampleName == "TTJets_scaledown"       ||
      gSampleName == "TTJetsFullLeptMGtauola" ||
      gSampleName == "TTJetsSemiLeptMGtauola") {
    
    Float_t Weight = 1.; 
    TLorentzVector top;
    for (size_t t=0; t<T_Gen_tSt3_Px->size(); t++){
      top.SetPxPyPzE(T_Gen_tSt3_Px->at(t),T_Gen_tSt3_Py->at(t),T_Gen_tSt3_Pz->at(t),T_Gen_tSt3_Energy->at(t));
      Float_t pt    = TMath::Min(top.Pt(), 400.);
      Float_t topSF = TMath::Exp(0.148 - 0.00129 * pt);
      Weight *= topSF;
    }
    fHTopPtWeight->Fill(TMath::Sqrt(Weight));
  }
  // Fill Gen Info 
  //----------------------------------------------------------------------------
  TLorentzVector lep,jet;
  Float_t minDRmu(999.),minDRel(999.);
  for (UInt_t jt=0; jt<T_Gen_bSt3_Energy->size(); jt++){
    jet.SetPxPyPzE(T_Gen_bSt3_Px->at(jt),T_Gen_bSt3_Py->at(jt),T_Gen_bSt3_Pz->at(jt), T_Gen_bSt3_Energy->at(jt));
    
    for (UInt_t mu=0; mu<T_Gen_MuonSt3_Energy->size(); mu++){
      lep.SetPxPyPzE(T_Gen_MuonSt3_Px->at(mu),T_Gen_MuonSt3_Py->at(mu),
		     T_Gen_MuonSt3_Pz->at(mu),T_Gen_MuonSt3_Energy->at(mu));
      
      if (minDRmu > lep.DeltaR(jet))  minDRmu = lep.DeltaR(jet);
    }
    for (UInt_t el=0; el<T_Gen_ElecSt3_Energy->size(); el++){
      lep.SetPxPyPzE(T_Gen_ElecSt3_Px->at(el),T_Gen_ElecSt3_Py->at(el),
		     T_Gen_ElecSt3_Pz->at(el),T_Gen_ElecSt3_Energy->at(el));

      if (minDRel > lep.DeltaR(jet))  minDRel = lep.DeltaR(jet);
    }
  }
  fHDeltaRLepJet[Muon] -> Fill(minDRmu);
  fHDeltaRLepJet[Elec] -> Fill(minDRel);
#endif


  
  // Accept only events with a good vertex
  //----------------------------------------------------------------------------
  if (SelectedVertexIndex() < 0) return;
  
  // Fill MiniTrees for further processing...
  //----------------------------------------------------------------------------
#ifdef __ISFR  
  FillTLRatios();
#endif
  //       FillAnalysisTree(0);
  FillYields();
  
  // Get SS Yields...
  fChargeSwitch = true;
  FillYields(); /// Get SS yields....
  fChargeSwitch = false;
  
  // Fill DY DD histograms
  if (gSampleName == "DoubleMu"        ||       
      gSampleName == "DoubleElectron"  || 
      gSampleName == "MuEG"            ||       
      gSampleName == "DYJets_Madgraph" ||
      gSampleName == "ZJets_Madgraph")  {
    FillDYHistograms();
  }
  
  // Fill Yields for syst. studies (only for MC)
  //----------------------------------------------------------------------------
  if (gIsData)         return;
  if (!gDoSystStudies) return;
  
  ResetDataMembers();
  gSysSource = BtagUp;
  FillYields(BtagUp);
  fChargeSwitch = true;
  FillYields(BtagUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetDataMembers();
  gSysSource = BtagDown;
  FillYields(BtagDown);
  fChargeSwitch = true;
  FillYields(BtagDown); /// Get SS yields....
  fChargeSwitch = false;

  ResetDataMembers();
  gSysSource = MisTagUp;
  FillYields(MisTagUp);
  fChargeSwitch = true;
  FillYields(MisTagUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetDataMembers();
  gSysSource = MisTagDown;
  FillYields(MisTagDown);
  fChargeSwitch = true;
  FillYields(MisTagDown); /// Get SS yields....
  fChargeSwitch = false;
  
  ResetDataMembers();
  SmearJetPts(1);
  gSysSource = JESUp;
  FillYields(JESUp);
  fChargeSwitch = true;
  FillYields(JESUp); /// Get SS yields....
  fChargeSwitch = false;
  
  ResetDataMembers();
  SmearJetPts(2);
  gSysSource = JESDown;
  FillYields(JESDown);
  fChargeSwitch = true;
  FillYields(JESDown); /// Get SS yields....
  fChargeSwitch = false;
 
  ResetDataMembers();
  SmearJetPts(3);
  gSysSource = JER;
  FillYields(JER);
  fChargeSwitch = true;
  FillYields(JER); /// Get SS yields....
  fChargeSwitch = false;
  
  // LEPTON SCALE
  ResetDataMembers();
  ScaleLeptons(1); //up
  gSysSource = LESUp;
  FillYields(LESUp);
  fChargeSwitch = true;
  FillYields(LESUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetDataMembers();
  ScaleLeptons(2); //down
  gSysSource = LESDown;
  FillYields(LESDown);
  fChargeSwitch = true;
  FillYields(LESDown); /// Get SS yields....
  fChargeSwitch = false;
   
  // PILE UP UNCERTAINTY
  ResetDataMembers();
#ifdef __ISMC
  PUSF = fPUWeightUp->GetWeight((int)T_Event_nTruePU);
#endif
  gSysSource = PUUp;
  FillYields(PUUp);
  
  ResetDataMembers();
#ifdef __ISMC
  PUSF = fPUWeightDown->GetWeight((int)T_Event_nTruePU);
#endif
  gSysSource = PUDown;
  FillYields(PUDown);
  
  
  // TOP PT
  ResetDataMembers();
  gSysSource = TopPtUp;
  FillYields(TopPtUp);
  fChargeSwitch = true;
  FillYields(TopPtUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetDataMembers();
  gSysSource = TopPtDown;
  FillYields(TopPtDown);
  fChargeSwitch = true;
  FillYields(TopPtDown); /// Get SS yields....
  fChargeSwitch = false;
  
  //
}// void(InsideLoop)
//------------------------------------------------------------------------------
// SetDataMembersAtTermination
//------------------------------------------------------------------------------
void TreeAnalysisTop::SetDataMembersAtTermination(){
  GetParameters();
  //  WriteHistos();
}
void TreeAnalysisTop::WriteHistos(){
  WriteTLRatios();
}
void TreeAnalysisTop::WriteTLRatios(){
  
}
//------------------------------------------------------------------------------
// Summary
//------------------------------------------------------------------------------
void TreeAnalysisTop::Summary(){
  
  //  AnalysisTree = ((TTree*) FindOutput("AnalysisTree"));
  //  cout << "---------------------------------------------------" <<endl;
  //  cout << "Number of entries in the Tree= " << AnalysisTree->GetEntries() <<endl;  
  //  cout << "---------------------------------------------------" <<endl;
}
//------------------------------------------------------------------------------
// GetParameters
//------------------------------------------------------------------------------
void TreeAnalysisTop::GetParameters()
{
  
  gSampleName   = GetInputParameters()->TheNamedString("sampleName");
  
  GetInputParameters()->TheNamedBool ("IsData",        gIsData);
  GetInputParameters()->TheNamedFloat("weight",        gWeight); // cross section / events in the sample
  GetInputParameters()->TheNamedFloat("LumiForPU",     gLumiForPU);
  GetInputParameters()->TheNamedFloat("TotalLumi",     gTotalLumi);
  GetInputParameters()->TheNamedBool ("DoSystStudies", gDoSystStudies);
  GetInputParameters()->TheNamedBool ("DoFR"         , gDoTLRatios);
  GetInputParameters()->TheNamedBool ("UseCSVM",       gUseCSVM);
  //  GetInputParameters()->TheNamedInt("SystDirection", gSysDirection);
}
//-----------------------------------------------------------------------------------
//  SET / RESET METHODS
//-----------------------------------------------------------------------------------
void TreeAnalysisTop::ResetHypLeptons(){
  TLorentzVector vec(0., 0., 0., 0.);
  fHypLepton1 = lepton(vec, 0, -1, -1);
  fHypLepton2 = lepton(vec, 0, -1, -1);
}
void TreeAnalysisTop::SetHypLepton1(int index, gChannel chan){
  TLorentzVector vec(0., 0., 0., 0.);
  if(chan == Muon){
    vec.SetPxPyPzE(MuPx.at(index),MuPy.at(index),T_Muon_Pz->at(index),T_Muon_Energy->at(index));
    fHypLepton1 = lepton(vec, T_Muon_Charge->at(index), 0, index);
  }
  else if(chan == Elec){
    vec.SetPxPyPzE(ElPx.at(index),ElPy.at(index),T_Elec_Pz->at(index),T_Elec_Energy->at(index));
    fHypLepton1 = lepton(vec, T_Elec_Charge->at(index), 1, index);
  }
  else exit(-1);
}
void TreeAnalysisTop::SetHypLepton2(int index, gChannel chan){
  TLorentzVector vec(0., 0., 0., 0.);
  if(chan == Muon){
    vec.SetPxPyPzE(MuPx.at(index),MuPy.at(index),T_Muon_Pz->at(index),T_Muon_Energy->at(index));
    fHypLepton2 = lepton(vec, T_Muon_Charge->at(index), 0, index);
  }
  else if(chan == Elec){
    vec.SetPxPyPzE(ElPx.at(index),ElPy.at(index),T_Elec_Pz->at(index),T_Elec_Energy->at(index));
    fHypLepton2 = lepton(vec, T_Elec_Charge->at(index), 1, index);
  }
  else exit(-1);
}

void TreeAnalysisTop::ResetAnalysisTree(){
  // Run variables
  TEvent = -999;  
  TLumi  = -999;   
  TRun   = -999;    
  
  // Sample variables
  TSName = "none";
  TSType = -999; 
  
  // Per-Event variables
  TSystFlag = -999; 
  TWeight   = -999;   
  TChannel  = -999;  
  TTLCat    = -999;    
  TNPV      = -999;      
  TMET      = -999;      
  TMET_Phi  = -999;  
   
  // Lepton variables
  TNTMus   = -999;   
  TNTEls   = -999;   
  TNMus    = -999;    
  TNEls    = -999;    
  TInvMass = -999; 
  
  TLep0Pt   = -999;  
  TLep0Eta  = -999; 
  TLep0Phi  = -999; 
  TLep0Ch   = -999;  
  TLep0Flav = -999;
  TLep1Pt   = -999;  
  TLep1Eta  = -999; 
  TLep1Phi  = -999; 
  TLep1Ch   = -999;  
  TLep1Flav = -999;

  // Jet variables
  TNJets     = -999;     
  TNbJets    = -999;    
  TNbJetsMed = -999; 
  TNJetsBtag = -999; 
  THT        = -999;        
  
  TJet0Px = -999; 
  TJet0Py = -999; 
  TJet0Pz = -999; 
  TJet0Et = -999; 
  TJet0E  = -999;  
  TJet1Px = -999; 
  TJet1Py = -999; 
  TJet1Pz = -999; 
  TJet1Et = -999; 
  TJet1E  = -999;  
  
  TBtagJet0   = -999;    
  TBtagJet1   = -999;    
  TBtagJet0Px = -999;   
  TBtagJet0Py = -999;   
  TBtagJet0Pz = -999;   
  TBtagJet0Et = -999;   
  TBtagJet0E  = -999;    
  TBtagJet1Px = -999;   
  TBtagJet1Py = -999;   
  TBtagJet1Pz = -999;   
  TBtagJet1Et = -999;   
  TBtagJet1E  = -999;   

  // Top pT reweight
  TtPx = -999;     
  TtPy = -999;     
  TtPz = -999;     
  TtE  = -999;     
  TtbarPx = -999;  
  TtbarPy = -999;  
  TtbarPz = -999;  
  TtbarE  = -999;   

}

//------------------------------------------------------------------------------
// SelectedVertexIndex
//------------------------------------------------------------------------------
Int_t TreeAnalysisTop::SelectedVertexIndex()
{
  Int_t goodVertexIndex = -999;
  nGoodVertex           =    0;
  
  for (UInt_t iVertex=0; iVertex<T_Vertex_z->size(); iVertex++) {
    
    if (fabs(T_Vertex_z ->at(iVertex)) < 24 &&
	T_Vertex_rho    ->at(iVertex)  <  2 &&
	T_Vertex_ndof   ->at(iVertex)  >  4 &&
	!T_Vertex_isFake->at(iVertex)) {
      
      nGoodVertex++;
      
      if (nGoodVertex == 1) goodVertexIndex = iVertex;
    }
  }
  
  return goodVertexIndex;
}
//------------------------------------------------------------------------------
// TRIGGER INFORMATION
//------------------------------------------------------------------------------
bool TreeAnalysisTop::PassTriggerMuMu()
{
  Bool_t pass = T_passTriggerDoubleMu; 
  if (T_Event_RunNumber==191090 ||  T_Event_RunNumber==193112 || T_Event_RunNumber==193116) pass=false; 
  return pass;
}
bool TreeAnalysisTop::PassTriggerEE()
{
  Bool_t pass = T_passTriggerDoubleEl; 
  if (T_Event_RunNumber==191090 ||  T_Event_RunNumber==193112 || T_Event_RunNumber==193116) pass=false; 
  return pass;
}
bool TreeAnalysisTop::PassTriggerEMu()
{
  Bool_t pass = T_passTriggerElMu; 
  if (T_Event_RunNumber==191090 ||  T_Event_RunNumber==193112 || T_Event_RunNumber==193116) pass=false; 
  return pass;
}
bool TreeAnalysisTop::PassSingleMuTrigger(){
  if (!gIsData) return true;
  bool passMu8(false),passMu17(false);
#ifdef __ISFR
   passMu8  = (T_HLT_Mu8_v16  || 
	       T_HLT_Mu8_v17  || 
	       T_HLT_Mu8_v18);
   
   passMu17 = (T_HLT_Mu17_v3  ||
	       T_HLT_Mu17_v4  ||
	       T_HLT_Mu17_v5);
#endif
  return (passMu8 || passMu17);
}
bool TreeAnalysisTop::PassSingleElTrigger(){
  if (!gIsData) return true;
  bool passEl8(false), passEl17(false);
#ifdef __ISFR
  passEl8  = (T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v ||
	      T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v);
  passEl17 = (T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v ||
	      T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v);
#endif
  return (passEl8 || passEl17);
}
bool TreeAnalysisTop::PassesJetPtdPhiCut(){
  
  vector<int> jetinds;
  for(unsigned int i = 0; i < T_JetAKCHS_Px->size(); ++i) {
    if(IsGoodJet(i, 50) ) {
      jetinds.push_back(i);
    }
  }
  if (jetinds.size() != 1) return false;
  
  TLorentzVector jet(T_JetAKCHS_Px    ->at(jetinds[0]),
		     T_JetAKCHS_Py    ->at(jetinds[0]),
		     T_JetAKCHS_Pz    ->at(jetinds[0]),
		     T_JetAKCHS_Energy->at(jetinds[0]));
  
  
  float dphi = jet.DeltaPhi(fHypLepton1.p);
  if (fabs(dphi) > 2.) return true;

  return false; 
}

//------------------------------------------------------------------------------
// SELECTORS
//------------------------------------------------------------------------------
int TreeAnalysisTop::HasLooseMuons(int &mu1, int &mu2){
  // Returns the number of loose muons and fills their indices in mu1 and mu2
  // Assumes the muons are sorted by pt in the minitree
  
  vector<int> loosemus;
  mu1 = -1;  mu2 = -1;
  for(unsigned int i = 0; i < T_Muon_Energy->size(); ++i) if(IsLooseMuon(i)) loosemus.push_back(i);
  if(loosemus.size() > 0) mu1 = loosemus[0];
  if(loosemus.size() > 1) mu2 = loosemus[1];
  return loosemus.size();
}
int TreeAnalysisTop::HasLooseMuons(){
  int ind1(-1), ind2(-1);
  return HasLooseMuons(ind1, ind2);
}
int TreeAnalysisTop::HasLooseElectrons(int &el1, int &el2){
  // Returns the number of loose electrons and fills their indices in el1 and el2
  // Assumes the electrons are sorted by pt in the minitree
  vector<int> looseels;
  el1 = -1;  el2 = -1;
  for(unsigned int i = 0; i < T_Elec_Energy->size(); ++i) if(IsLooseElectron(i)) looseels.push_back(i);
  if(looseels.size() > 0) el1 = looseels[0];
  if(looseels.size() > 1) el2 = looseels[1];
  return looseels.size();
}
int TreeAnalysisTop::HasLooseElectrons(){
  int ind1(-1), ind2(-1);
  return HasLooseElectrons(ind1, ind2);
}
bool TreeAnalysisTop::IsSigSupMuEvent(int &mu1){
  int mu2(-1);
  if(HasLooseMuons(mu1, mu2) < 1) return false;
  
  SetHypLepton1(mu1, Muon);
  if(!PassesJetPtdPhiCut())  return false;
  if(getMT(mu1,Muon) > 20.)  return false;
  if(getMET()        > 20.)  return false;
  
  int nmus(0);
  for (unsigned int i=0; i< T_Muon_Energy->size(); ++i){
    if (IsLooseMuon(i)) nmus++;
  }
  
  if (nmus > 1)                   return false;
  if (T_Muon_Energy->size() > 1)  return false;
  
  return true;
}
bool TreeAnalysisTop::IsSigSupElEvent(int &el1){
  int el2(-1);
  if(HasLooseElectrons(el1, el2) < 1) return false;
  
  SetHypLepton1(el1, Elec);
  if(!PassesJetPtdPhiCut())  return false;
  if(getMT(el1,Elec) > 20.)  return false;
  if(getMET()        > 20.)  return false;
  
  int nels(0);
  for (unsigned int i=0; i< T_Elec_Energy->size(); ++i){
    if (IsLooseElectron(i)) nels++;
  }
  
  if (nels > 1)                   return false;
  if (T_Elec_Energy->size() > 1)  return false;
  
  return true;
}
bool TreeAnalysisTop::IsZMuMuEvent(int &mu1, int &mu2){
  if(HasLooseMuons(mu1, mu2) < 2)  return false;
  if(T_Muon_Charge->at(mu1) == T_Muon_Charge->at(mu2)) return false; // os
  
  // Z mass window cut
  TLorentzVector p1, p2;
  p1.SetPxPyPzE(MuPx.at(mu1), MuPy.at(mu1), T_Muon_Pz->at(mu1), T_Muon_Energy->at(mu1));
  p2.SetPxPyPzE(MuPx.at(mu2), MuPy.at(mu2), T_Muon_Pz->at(mu2), T_Muon_Energy->at(mu2));
  
  double m = (p1+p2).M();
  if(fabs(91. - m) > 15.) return false;
  
  SetHypLepton1(mu1, Muon);
  SetHypLepton2(mu2, Muon);
  
  if(getMET() < 20.) return false;
  if(getNJets() < 2) return false;
  return true;
}
bool TreeAnalysisTop::IsZElElEvent(int &el1, int &el2){
  if(HasLooseElectrons(el1, el2) < 2)  return false;
  if(T_Elec_Charge[el1] == T_Elec_Charge[el2]) return false; // os
  
  // Z mass window cut
  TLorentzVector p1, p2;
  p1.SetPxPyPzE(ElPx.at(el1), ElPy.at(el1), T_Elec_Pz->at(el1), T_Elec_Energy->at(el1));
  p2.SetPxPyPzE(ElPx.at(el2), ElPy.at(el2), T_Elec_Pz->at(el2), T_Elec_Energy->at(el2));
  
  double m = (p1+p2).M();
  if(fabs(91. - m) > 15.) return false;
  
  SetHypLepton1(el1, Elec);
  SetHypLepton2(el2, Elec);
  
  if(getMET() < 20.) return false;
  if(getNJets() < 2) return false;
  return true;
}
//////////////////////////////////////////////////////////////////
/// Get METHODS
/////////////////////////////////////////////////////////////////
float TreeAnalysisTop::getHT(){
  float ht(0);
  for (unsigned int i=0; i<T_JetAKCHS_Energy->size(); i++) if(IsGoodJet(i,gJetEtCut)) ht+=JetEt.at(i);
  
  return ht;
}
float TreeAnalysisTop::getJetPtIndex(unsigned int ind){
  vector<float> jetpt;

  for (unsigned int i=0; i<T_JetAKCHS_Energy->size(); i++) 
    if(IsGoodJet(i,gJetEtCut)) jetpt.push_back(JetEt.at(i));
  
  if (jetpt.size() <= ind) return -999.;
  
  return jetpt[ind];
}
float TreeAnalysisTop::getBtagJetPtIndex(unsigned int ind){
  vector<float> jetpt;
  
  int btagSys = 0;
    
  for (unsigned int i=0; i<T_JetAKCHS_Energy->size(); i++){ 
    if (!IsGoodJet(i,gJetEtCut))                                   continue;
        
    if(gIsData  && !(fBTagSF->IsTagged(T_JetAKCHS_Tag_CombSVtx->at(i),-999999, JetEt.at(i), 
				       T_JetAKCHS_Eta->at(i), btagSys))) continue;

    // Split btag-efficiency and mistag-rate
    if(!gIsData) {
      if(TMath::Abs(T_JetAKCHS_Parton_Flavour->at(i)) == 5 || TMath::Abs(T_JetAKCHS_Parton_Flavour->at(i)) == 4){
	if (gSysSource == BtagUp)     btagSys =  1;
	if (gSysSource == BtagDown)   btagSys = -1;
	if (gSysSource == MisTagUp)   btagSys =  0;
	if (gSysSource == MisTagDown) btagSys =  0;
      }
      if(TMath::Abs(T_JetAKCHS_Parton_Flavour->at(i)) != 5 || TMath::Abs(T_JetAKCHS_Parton_Flavour->at(i)) != 4){
	if (gSysSource == BtagUp)     btagSys =  0;
	if (gSysSource == BtagDown)   btagSys =  0;
	if (gSysSource == MisTagUp)   btagSys =  1;
	if (gSysSource == MisTagDown) btagSys = -1;
      }
    }

    if(!gIsData && !(fBTagSF->IsTagged(T_JetAKCHS_Tag_CombSVtx->at(i),T_JetAKCHS_Parton_Flavour->at(i), 
				       JetEt.at(i), T_JetAKCHS_Eta->at(i), btagSys))) continue;

    jetpt.push_back(JetEt.at(i));
    
  }

  if (jetpt.size() <= ind) return -999.;
  
  return jetpt[ind];
}
float TreeAnalysisTop::getMT(int ind, gChannel chan){
  // Calculates MT
  
  TLorentzVector pmet, plep;
  if (chan == Muon) plep.SetPxPyPzE(MuPx.at(ind), MuPy.at(ind), T_Muon_Pz->at(ind), T_Muon_Energy->at(ind));
  if (chan == Elec) plep.SetPxPyPzE(ElPx.at(ind), ElPy.at(ind), T_Elec_Pz->at(ind), T_Elec_Energy->at(ind));
  
  pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
  double ETlept = sqrt(plep.M2() + plep.Perp2());
  
  double MT = sqrt( 2*(ETlept*getMET() - plep.Px()*pmet.Px() - plep.Py()*pmet.Py()));
  return MT;
}
float TreeAnalysisTop::getJERScale(int jet){
  float eta = T_JetAKCHS_Eta->at(jet);
  if(     TMath::Abs(eta) < 0.5) return 1.052;
  else if(TMath::Abs(eta) < 1.1) return 1.057;
  else if(TMath::Abs(eta) < 1.7) return 1.096;
  else if(TMath::Abs(eta) < 2.3) return 1.134;
  else                           return 1.288;
}
float TreeAnalysisTop::getErrPt(float Pt, float Eta) {
  float InvPerr2;
  float N(0.), S(0.), C(0.), m(0.);
  
  if(TMath::Abs(Eta) < 0.5 ) {
    N = 3.96859;
    S = 0.18348;
    C = 0.;
    m = 0.62627;
  } else if( TMath::Abs(Eta) < 1.  ) {
    N = 3.55226;
    S = 0.24026;
    C = 0.;
    m = 0.52571;
  } else if( TMath::Abs(Eta) < 1.5  ) {
    N = 4.54826;
    S = 0.22652;
    C = 0.;
    m = 0.58963;
  } else if( TMath::Abs(Eta) < 2.  ) {
    N = 4.62622;
    S = 0.23664;
    C = 0.;
    m = 0.48738;
  } else if( TMath::Abs(Eta) < 2.5  ) {
    N = 2.53324;
    S = 0.34306;
    C = 0.;
    m = 0.28662;
  } else if( TMath::Abs(Eta) < 3.  ) {
    N = -3.33814;
    S = 0.73360;
    C = 0.;
    m = 0.08264;
  } else if( TMath::Abs(Eta) < 5.  ) {
    N = 2.95397;
    S = 0.11619;
    C = 0.;
    m = 0.96086;
  }
  // this is the absolute resolution (squared), not sigma(pt)/pt
  // so have to multiply by pt^2, thats why m+1 instead of m-1
  InvPerr2 =  (N * TMath::Abs(N) ) + (S * S) * pow(Pt, m+1) + (C * C) * Pt * Pt ;
  return sqrt(InvPerr2);
}
float TreeAnalysisTop::getLeptonError(gChannel chan){
  float err1(0.), err2(0.);
  int ind1 = fHypLepton1.index; 
  int ind2 = fHypLepton2.index;
  if (chan==Muon){
    err1 = fLeptonSF->GetTightMuonSF(T_Muon_Pt->at(ind1), T_Muon_Eta->at(ind1));
    err2 = fLeptonSF->GetTightMuonSF(T_Muon_Pt->at(ind2), T_Muon_Eta->at(ind2));
  }
  if (chan==ElMu){
    err1 = fLeptonSF->GetTightMuonSF    (T_Muon_Pt->at(ind1), T_Muon_Eta->at(ind1));
    err2 = fLeptonSF->GetTightElectronSF(T_Elec_Pt->at(ind2), T_Elec_Eta->at(ind2));
  }
  if (chan==Elec){
    err1 = fLeptonSF->GetTightElectronSF(T_Elec_Pt->at(ind1), T_Elec_Eta->at(ind1));
    err2 = fLeptonSF->GetTightElectronSF(T_Elec_Pt->at(ind2), T_Elec_Eta->at(ind2));
  }
  return TMath::Sqrt(err1*err1+err2*err2);
}
//  if (chan==Muon){
//    err1 = 0.0054;
//    err2 = 0.0054;
//  }
//  if (chan==ElMu){
//    err1 = 0.0054;
//    if (TMath::Abs(fHypLepton2.p.Eta()) < 1.5){
//      if      (fHypLepton2.p.Pt() < 30)                             err2 = 0.014;
//      else if (fHypLepton2.p.Pt() >= 30 && fHypLepton2.p.Pt() < 40) err2 = 0.0028;
//      else if (fHypLepton2.p.Pt() >= 40 && fHypLepton2.p.Pt() < 50) err2 = 0.0014;
//      else                                                          err2 = 0.0041;
//    }
//    else {
//      if      (fHypLepton2.p.Pt() < 30)                             err2 = 0.022;
//      else if (fHypLepton2.p.Pt() >= 30 && fHypLepton2.p.Pt() < 40) err2 = 0.0059;
//      else if (fHypLepton2.p.Pt() >= 40 && fHypLepton2.p.Pt() < 50) err2 = 0.0030;
//      else                                                          err2 = 0.0053;
//    }
//  }
//  if (chan==Elec){
//    if (TMath::Abs(fHypLepton1.p.Eta()) < 1.5){
//      if      (fHypLepton1.p.Pt() < 30)                             err1 = 0.014;
//      else if (fHypLepton1.p.Pt() >= 30 && fHypLepton1.p.Pt() < 40) err1 = 0.0028;
//      else if (fHypLepton1.p.Pt() >= 40 && fHypLepton1.p.Pt() < 50) err1 = 0.0014;
//      else                                                          err1 = 0.0041;
//    }
//    else {
//      if      (fHypLepton1.p.Pt() < 30)                             err1 = 0.022;
//      else if (fHypLepton1.p.Pt() >= 30 && fHypLepton1.p.Pt() < 40) err1 = 0.0059;
//      else if (fHypLepton1.p.Pt() >= 40 && fHypLepton1.p.Pt() < 50) err1 = 0.0030;
//      else                                                          err1 = 0.0053;
//    }
//
//    if (TMath::Abs(fHypLepton2.p.Eta()) < 1.5){
//      if      (fHypLepton2.p.Pt() < 30)                             err2 = 0.014;
//      else if (fHypLepton2.p.Pt() >= 30 && fHypLepton2.p.Pt() < 40) err2 = 0.0028;
//      else if (fHypLepton2.p.Pt() >= 40 && fHypLepton2.p.Pt() < 50) err2 = 0.0014;
//      else                                                          err2 = 0.0041;
//    }
//    else {
//      if      (fHypLepton2.p.Pt() < 30)                             err2 = 0.022;
//      else if (fHypLepton2.p.Pt() >= 30 && fHypLepton2.p.Pt() < 40) err2 = 0.0059;
//      else if (fHypLepton2.p.Pt() >= 40 && fHypLepton2.p.Pt() < 50) err2 = 0.0030;
//      else                                                          err2 = 0.0053;
//    }
//  }
float TreeAnalysisTop::getTriggerError(gChannel chan){
  float trig(0.);
  int ind1 = fHypLepton1.index; 
  int ind2 = fHypLepton2.index;
  if (chan==Muon) trig = fLeptonSF->GetDoubleMuSF(T_Muon_Eta->at(ind1),T_Muon_Eta->at(ind2));
  if (chan==ElMu) trig = fLeptonSF->GetMuEGSF    (T_Elec_Eta->at(ind2),T_Muon_Eta->at(ind1));
  if (chan==Elec) trig = fLeptonSF->GetDoubleElSF(T_Elec_Eta->at(ind1),T_Elec_Eta->at(ind2));
  return trig;
}
float TreeAnalysisTop::getSF(gChannel chan, int ind1, int ind2) {
  if (gIsData)              return 1.; //Don't scale data
  if (ind1 < 0 || ind2 < 0) return 1.; //sanity check!!; 
  
  //  float Id1_err(1.), Id2_err(1.), Trig_err(1.);
  float Id   = 1.;
  float Trig = 1.;
  if (chan == Muon){
    Id   = fLeptonSF->GetTightMuonSF(T_Muon_Pt->at(ind1), T_Muon_Eta->at(ind1));
    Id  *= fLeptonSF->GetTightMuonSF(T_Muon_Pt->at(ind2), T_Muon_Eta->at(ind2));
    Trig = fLeptonSF->GetDoubleMuSF (T_Muon_Eta->at(ind1),T_Muon_Eta->at(ind2)) ;
  } 
  else if (chan == Elec){
    Id   = fLeptonSF->GetTightElectronSF(T_Elec_Pt->at(ind1), T_Elec_Eta->at(ind1));
    Id  *= fLeptonSF->GetTightElectronSF(T_Elec_Pt->at(ind2), T_Elec_Eta->at(ind2));
    Trig = fLeptonSF->GetDoubleElSF     (T_Elec_Eta->at(ind1),T_Elec_Eta->at(ind2));
  }
  else if (chan == ElMu){
    Id   = fLeptonSF->GetTightMuonSF    (T_Muon_Pt->at(ind1), T_Muon_Eta->at(ind1));
    Id  *= fLeptonSF->GetTightElectronSF(T_Elec_Pt->at(ind2), T_Elec_Eta->at(ind2));
    Trig = fLeptonSF->GetMuEGSF         (T_Elec_Eta->at(ind2),T_Muon_Eta->at(ind1));
  }
  return (PUSF*Id*Trig);
}
float TreeAnalysisTop::getTopPtSF(){
  // Return SF of the pt pt of the top 
  // Only apply SF if the process is ttbar...
  if(!gSampleName.Contains("TTJets")) return 1.;
  
  TLorentzVector top;
  Float_t topSF = 0.;
  Float_t Weight = 1.; 
#ifdef __ISMC
  if (T_Gen_tSt3_Px->size() != 2) return 1.;
  
  for (size_t t=0; t<T_Gen_tSt3_Px->size(); t++){
    top.SetPxPyPzE(T_Gen_tSt3_Px->at(t),T_Gen_tSt3_Py->at(t),T_Gen_tSt3_Pz->at(t),T_Gen_tSt3_Energy->at(t));
    Float_t pt = TMath::Min(top.Pt(), 400.);
    topSF = TMath::Exp(0.148 - 0.00129 * pt);
    Weight *= topSF;
  }
  Weight = TMath::Sqrt(Weight);
#endif
  
  if (gSysSource == TopPtDown) { return 1.;          }
  if (gSysSource == TopPtUp)   { return 2*Weight-1.; }
  
  return Weight;
}
void TreeAnalysisTop::FillDYHistograms(){
  float Mll = 0.;
  int ind1(-1),ind2(-1);
  if (PassTriggerEMu()  && IsElMuEvent(ind1,ind2)){
    // Define Hypothesis Leptons...
    SetHypLepton1(ind1, Muon);
    SetHypLepton2(ind2, Elec);
    if (IsTightMuon(ind1) && IsTightElectron(ind2)){
      EventWeight = gWeight * getSF(ElMu,ind1,ind2);
      Mll = (fHypLepton1.p+fHypLepton2.p).M();
      
      if (PassesMllVeto() && PassesMuonEta2p1(ElMu) && Passes3rdLeptonVeto()){
	fHDY_InvMassVsNPV   [ElMu][iDilepton]->Fill(nGoodVertex, Mll, EventWeight);
	fHDY_InvMassVsMET   [ElMu][iDilepton]->Fill(getMET()   , Mll, EventWeight);
	fHDY_InvMassVsNjets [ElMu][iDilepton]->Fill(getNJets() , Mll, EventWeight);
	fHDY_InvMassVsNbtags[ElMu][iDilepton]->Fill(getNBTags(), Mll, EventWeight);
	fHDY_InvMass        [ElMu][iDilepton]->Fill(             Mll, EventWeight);
	
	if (PassesNJetsCut()) {
	  fHDY_InvMassVsNPV   [ElMu][i2jets]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [ElMu][i2jets]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [ElMu][i2jets]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[ElMu][i2jets]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [ElMu][i2jets]->Fill(             Mll, EventWeight);

	  if (PassesMETCut())   {
	    fHDY_InvMassVsNPV   [ElMu][iMET]->Fill(nGoodVertex, Mll, EventWeight);
	    fHDY_InvMassVsMET   [ElMu][iMET]->Fill(getMET()   , Mll, EventWeight);
	    fHDY_InvMassVsNjets [ElMu][iMET]->Fill(getNJets() , Mll, EventWeight);
	    fHDY_InvMassVsNbtags[ElMu][iMET]->Fill(getNBTags(), Mll, EventWeight);
	    fHDY_InvMass        [ElMu][iMET]->Fill(             Mll, EventWeight);
	    
	    if (PassesNBtagCut()) {
	      fHDY_InvMassVsNPV   [ElMu][i1btag]->Fill(nGoodVertex, Mll, EventWeight);
	      fHDY_InvMassVsMET   [ElMu][i1btag]->Fill(getMET()   , Mll, EventWeight);
	      fHDY_InvMassVsNjets [ElMu][i1btag]->Fill(getNJets() , Mll, EventWeight);
	      fHDY_InvMassVsNbtags[ElMu][i1btag]->Fill(getNBTags(), Mll, EventWeight);
	      fHDY_InvMass        [ElMu][i1btag]->Fill(             Mll, EventWeight);
	    }
	  }
	}
      }
    }
  }
  
  ResetHypLeptons(); 
  if (PassTriggerMuMu() && IsMuMuEvent(ind1,ind2)){
    SetHypLepton1(ind1, Muon);
    SetHypLepton2(ind2, Muon);
    
    if (IsTightMuon(ind1) && IsTightMuon(ind2)){
      EventWeight = gWeight * getSF(Muon,ind1,ind2);
      Mll = (fHypLepton1.p+fHypLepton2.p).M();
      
      if (PassesMllVeto() && PassesMuonEta2p1(Muon) && Passes3rdLeptonVeto()){
	fHDY_InvMassVsNPV   [Muon][iDilepton]->Fill(nGoodVertex, Mll, EventWeight);
	fHDY_InvMassVsMET   [Muon][iDilepton]->Fill(getMET()   , Mll, EventWeight);
	fHDY_InvMassVsNjets [Muon][iDilepton]->Fill(getNJets() , Mll, EventWeight);
	fHDY_InvMassVsNbtags[Muon][iDilepton]->Fill(getNBTags(), Mll, EventWeight);
	fHDY_InvMass        [Muon][iDilepton]->Fill(             Mll, EventWeight);
	
	if (PassesNJetsCut()) {
	  fHDY_InvMassVsNPV   [Muon][i2jets]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [Muon][i2jets]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [Muon][i2jets]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[Muon][i2jets]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [Muon][i2jets]->Fill(             Mll, EventWeight);
	  
	  if (PassesMETCut())   {
	    fHDY_InvMassVsNPV   [Muon][iMET]->Fill(nGoodVertex, Mll, EventWeight);
	    fHDY_InvMassVsMET   [Muon][iMET]->Fill(getMET()   , Mll, EventWeight);
	    fHDY_InvMassVsNjets [Muon][iMET]->Fill(getNJets() , Mll, EventWeight);
	    fHDY_InvMassVsNbtags[Muon][iMET]->Fill(getNBTags(), Mll, EventWeight);
	    fHDY_InvMass        [Muon][iMET]->Fill(             Mll, EventWeight);
	    
	    if (PassesNBtagCut()) {
	      fHDY_InvMassVsNPV   [Muon][i1btag]->Fill(nGoodVertex, Mll, EventWeight);
	      fHDY_InvMassVsMET   [Muon][i1btag]->Fill(getMET()   , Mll, EventWeight);
	      fHDY_InvMassVsNjets [Muon][i1btag]->Fill(getNJets() , Mll, EventWeight);
	      fHDY_InvMassVsNbtags[Muon][i1btag]->Fill(getNBTags(), Mll, EventWeight);
	      fHDY_InvMass        [Muon][i1btag]->Fill(             Mll, EventWeight);
	    }
	  }
	}
      }
    }
  }

  ResetHypLeptons(); 
  if (PassTriggerEE()   && IsElElEvent(ind1,ind2)){
    SetHypLepton1(ind1, Elec);
    SetHypLepton2(ind2, Elec);
    
    if (IsTightElectron(ind1) && IsTightElectron(ind2)){
      EventWeight = gWeight * getSF(Elec,ind1,ind2);
      Mll = (fHypLepton1.p+fHypLepton2.p).M();
      
      if (PassesMllVeto() && PassesMuonEta2p1(Elec) && Passes3rdLeptonVeto()){
	fHDY_InvMassVsNPV   [Elec][iDilepton]->Fill(nGoodVertex, Mll, EventWeight);
	fHDY_InvMassVsMET   [Elec][iDilepton]->Fill(getMET()   , Mll, EventWeight);
	fHDY_InvMassVsNjets [Elec][iDilepton]->Fill(getNJets() , Mll, EventWeight);
	fHDY_InvMassVsNbtags[Elec][iDilepton]->Fill(getNBTags(), Mll, EventWeight);
	fHDY_InvMass        [Elec][iDilepton]->Fill(             Mll, EventWeight);
	
	if (PassesNJetsCut()) {
	  fHDY_InvMassVsNPV   [Elec][i2jets]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [Elec][i2jets]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [Elec][i2jets]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[Elec][i2jets]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [Elec][i2jets]->Fill(             Mll, EventWeight);

	  if (PassesMETCut())   {
	    fHDY_InvMassVsNPV   [Elec][iMET]->Fill(nGoodVertex, Mll, EventWeight);
	    fHDY_InvMassVsMET   [Elec][iMET]->Fill(getMET()   , Mll, EventWeight);
	    fHDY_InvMassVsNjets [Elec][iMET]->Fill(getNJets() , Mll, EventWeight);
	    fHDY_InvMassVsNbtags[Elec][iMET]->Fill(getNBTags(), Mll, EventWeight);
	    fHDY_InvMass        [Elec][iMET]->Fill(             Mll, EventWeight);
	    
	    if (PassesNBtagCut()) {
	      fHDY_InvMassVsNPV   [Elec][i1btag]->Fill(nGoodVertex, Mll, EventWeight);
	      fHDY_InvMassVsMET   [Elec][i1btag]->Fill(getMET()   , Mll, EventWeight);
	      fHDY_InvMassVsNjets [Elec][i1btag]->Fill(getNJets() , Mll, EventWeight);
	      fHDY_InvMassVsNbtags[Elec][i1btag]->Fill(getNBTags(), Mll, EventWeight);
	      fHDY_InvMass        [Elec][i1btag]->Fill(             Mll, EventWeight);
	    }
	  }
	}
      }
    }
  }
  ResetHypLeptons();
}
void TreeAnalysisTop::FillKinematicHistos(gChannel chan, iCut cut){
  if (gSysSource != Norm)      return;  //only fill histograms for nominal distributions...
  if (fChargeSwitch == true  ) return;
  if (fHypLepton1.index == -1) return;
  if (fHypLepton2.index == -1) return;
  
  //++ met info
  fHMET[chan][cut]        ->Fill(getMET(),                                EventWeight);
  
  //++ lepton info
  //  fHInvMass[chan][cut]    ->Fill((fHypLepton1.p+fHypLepton2.p).M(),       EventWeight);
  fHDiLepPt[chan][cut]    ->Fill((fHypLepton1.p+fHypLepton2.p).Pt(),      EventWeight);
  fHLep0Pt[chan][cut]     ->Fill(fHypLepton1.p.Pt(),                      EventWeight);
  fHLep1Pt[chan][cut]     ->Fill(fHypLepton2.p.Pt(),                      EventWeight);
  fHDelLepPhi[chan][cut]  ->Fill(fHypLepton1.p.DeltaPhi((fHypLepton2.p)), EventWeight);
  			  
  //++ jet info		  
  int njets = getNJets(); 
  fHNJets[chan][cut]      ->Fill(njets,                                   EventWeight);
  fHHT[chan][cut]         ->Fill(getHT(),                                 EventWeight);
  fHNBtagJets[chan][cut]  ->Fill(getNBTags(),                             EventWeight);
  fHJet0Pt[chan][cut]     ->Fill(getJetPtIndex(0),                        EventWeight);
  fHJet1Pt[chan][cut]     ->Fill(getJetPtIndex(1),                        EventWeight);
  fHBtagJet0Pt[chan][cut] ->Fill(getBtagJetPtIndex(0),                    EventWeight);
  
//  if (njets == 1) fHNBtagsNJets[chan][cut][0]->Fill(getNBTags(),             EventWeight);
//  if (njets == 2) fHNBtagsNJets[chan][cut][0]->Fill(getNBTags()+5,           EventWeight);
//  if (njets == 3) fHNBtagsNJets[chan][cut][0]->Fill(getNBTags()+10,          EventWeight);
//  if (njets >= 4) fHNBtagsNJets[chan][cut][0]->Fill(getNBTags()+15,          EventWeight);

//  if (njets > 0)
//  fHCSVTag[chan][cut] ->Fill(T_JetAKCHS_Tag_CombSVtx->at(getLeadingJet()), EventWeight);
  if (njets > 1)
    fHCSVTag[chan][cut] ->Fill(T_JetAKCHS_Tag_CombSVtx->at(getSecondLeadingJet()), EventWeight);

  fHTopD[chan][cut] ->Fill(getTopD(), EventWeight);
  fHDelPhillJet[chan][cut]->Fill(getDeltaPhillJet(), EventWeight);
  
  //// Top/Z diff topology.
  fHDRLep[chan][cut]        ->Fill(fHypLepton1.p.DeltaR(fHypLepton2.p),     EventWeight);
  if (chan == Muon){
    fHLep0Iso[chan][cut]      ->Fill(getMuonIso(fHypLepton1.index) ,EventWeight);
    fHLep1Iso[chan][cut]      ->Fill(getMuonIso(fHypLepton2.index), EventWeight);
  }
  else if (chan==Elec){
    fHLep0Iso[chan][cut]      ->Fill(getElecIso(fHypLepton1.index) ,EventWeight);
    fHLep1Iso[chan][cut]      ->Fill(getElecIso(fHypLepton2.index), EventWeight);
  }
  else if (chan == ElMu){
    fHLep0Iso[chan][cut]      ->Fill(getMuonIso(fHypLepton1.index) ,EventWeight);
    fHLep1Iso[chan][cut]      ->Fill(getElecIso(fHypLepton2.index), EventWeight);
  }
  
  if (njets > 0) {
    fHDRLep0Jet[chan][cut]    ->Fill(getDRClosestJet(fHypLepton1.p),   EventWeight);
    fHDRLep1Jet[chan][cut]    ->Fill(getDRClosestJet(fHypLepton2.p),   EventWeight);
    fHDPhiLep0Jet[chan][cut]  ->Fill(getDPhiClosestJet(fHypLepton1.p), EventWeight);
    fHDPhiLep1Jet[chan][cut]  ->Fill(getDPhiClosestJet(fHypLepton2.p), EventWeight);
  }

  //  fHAbsDelPhiLep[chan][cut]->Fill(abs(fHypLepton1.p.DeltaPhi((fHypLepton2.p))), EventWeight);
  fHAbsDelPhiLep[chan][cut] ->Fill(TMath::Abs(fHypLepton1.p.DeltaPhi((fHypLepton2.p)))/TMath::Pi(), EventWeight);
#ifdef __ISSTOP  
  if(gSampleName == "T2tt_150to250LSP1to100_LeptonFilter"){
    for (size_t t=0; t<T_Gen_StopMass->size(); t++)
      fHStopMass[chan][cut] ->Fill(T_Gen_StopMass->at(t), EventWeight);
    for (size_t t=0; t<T_Gen_Chi0Mass->size(); t++)
      fHChi0Mass[chan][cut] ->Fill(T_Gen_Chi0Mass->at(t), EventWeight);
    
    fHChi0StopMass[chan][cut] ->Fill(T_Gen_StopMass->at(0),
				     T_Gen_Chi0Mass->at(0), EventWeight);
  }
#endif

}
void TreeAnalysisTop::FillYieldsHistograms(gChannel chan, iCut cut, gSystFlag sys){
  if (fChargeSwitch){   fHSSyields[chan][sys]->Fill(cut, EventWeight);  }
  else {                fHyields[chan][sys]  ->Fill(cut, EventWeight);  }
  
  /// FOR SYSTEMATIC STUDIES
  int njets  = 0; njets  = getNJets();
  int nbtags = 0; nbtags = getNBTags();
  
  if (fChargeSwitch) { 
    fHSSInvMass[chan][cut][sys]->Fill((fHypLepton1.p+fHypLepton2.p).M(), EventWeight);
    if (njets == 0) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags,        EventWeight);
    if (njets == 1) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+1,      EventWeight);
    if (njets == 2) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+3,      EventWeight);
    if (njets == 3) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+6,      EventWeight);
    if (njets >= 4) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+10,     EventWeight);
  }
  else {
    fHInvMass[chan][cut][sys]->Fill((fHypLepton1.p+fHypLepton2.p).M(), EventWeight);
    if (njets == 0) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags,        EventWeight);
    if (njets == 1) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+1,      EventWeight);
    if (njets == 2) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+3,      EventWeight);
    if (njets == 3) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+6,      EventWeight);
    if (njets >= 4) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+10,     EventWeight);
  }
  
  if (!gIsData){
    fHLepSys[chan][cut]->Fill(getLeptonError(chan));
    fHTrigSys[chan][cut]->Fill(getTriggerError(chan));
//    // FOR SS ORIGINS
//    if (fChargeSwitch) fHSSOrigins[chan][cut]->Fill();
//    else               fHOrigins[chan][cut]  ->Fill();
    
  }
  return;
}
void TreeAnalysisTop::FillYields(gSystFlag sys){
  ResetHypLeptons();
  
  int ind1(-1),ind2(-1);
  if (PassTriggerEMu()  && IsElMuEvent(ind1,ind2)){
    // Define Hypothesis Leptons...
    SetHypLepton1(ind1, Muon);
    SetHypLepton2(ind2, Elec);
    if (IsTightMuon(ind1) && IsTightElectron(ind2)){
      EventWeight = gWeight * getSF(ElMu,ind1,ind2) * getTopPtSF();
#ifdef __ISSTOP
      if(gSampleName == "T2tt_150to250LSP1to100_LeptonFilter")
	EventWeight = EventWeight * T_Gen_polWeights->at(10);
#endif
     if (PassesMllVeto() && PassesMuonEta2p1(ElMu) && Passes3rdLeptonVeto()){
	FillYieldsHistograms(ElMu, iDilepton, sys);
	if(sys==Norm) FillKinematicHistos(ElMu,iDilepton);

	FillYieldsHistograms(ElMu, iZVeto, sys);      
	if(sys==Norm) FillKinematicHistos(ElMu, iZVeto);
	
	FillYieldsHistograms(ElMu, iMET, sys);      
	if(sys==Norm) FillKinematicHistos(ElMu,iMET);
	
	if (PassesNJetsCut()) {
	  FillYieldsHistograms(ElMu, i2jets, sys);      
	  if(sys==Norm) FillKinematicHistos(ElMu,i2jets);
	  if (PassesNBtagCut()) {
	    FillYieldsHistograms(ElMu, i1btag, sys);      
	    if(sys==Norm) FillKinematicHistos(ElMu,i1btag);
	  }
	}
      }
    }
  }
  
  ResetHypLeptons(); 
  if (PassTriggerMuMu() && IsMuMuEvent(ind1,ind2)){
    SetHypLepton1(ind1, Muon);
    SetHypLepton2(ind2, Muon);
    if (IsTightMuon(ind1) && IsTightMuon(ind2)){
      EventWeight = gWeight * getSF(Muon,ind1,ind2)  * getTopPtSF();
#ifdef __ISSTOP
      if(gSampleName == "T2tt_150to250LSP1to100_LeptonFilter")
	EventWeight = EventWeight * T_Gen_polWeights->at(10);
#endif
      if (PassesMllVeto() && PassesMuonEta2p1(Muon) && Passes3rdLeptonVeto()){
	FillYieldsHistograms(Muon,iDilepton, sys);
	if(sys==Norm) FillKinematicHistos(Muon,iDilepton);
	if (PassesZVeto())    {
	  FillYieldsHistograms(Muon,iZVeto, sys);      
	  if(sys==Norm) FillKinematicHistos(Muon,iZVeto);
	  if (PassesMETCut())   {
	    FillYieldsHistograms(Muon,iMET, sys);      
	    if(sys==Norm) FillKinematicHistos(Muon,iMET);
	    if (PassesNJetsCut()) {
	      FillYieldsHistograms(Muon,i2jets, sys);      
	      if(sys==Norm) FillKinematicHistos(Muon,i2jets);
	      if (PassesNBtagCut()) {
		FillYieldsHistograms(Muon,i1btag, sys);      
		if(sys==Norm) FillKinematicHistos(Muon,i1btag);
	      }
	    }
	  }
	}
      }
    }
  }

  ResetHypLeptons(); 
  if (PassTriggerEE()   && IsElElEvent(ind1,ind2)){
    SetHypLepton1(ind1, Elec);
    SetHypLepton2(ind2, Elec);
    
    if (IsTightElectron(ind1) && IsTightElectron(ind2)){
      EventWeight = gWeight * getSF(Elec,ind1,ind2) * getTopPtSF();     
#ifdef __ISSTOP
      if(gSampleName == "T2tt_150to250LSP1to100_LeptonFilter")
	EventWeight = EventWeight * T_Gen_polWeights->at(10);
#endif
      if (PassesMllVeto() && PassesMuonEta2p1(Elec) && Passes3rdLeptonVeto()){
	FillYieldsHistograms(Elec,iDilepton, sys);
	if(sys==Norm) FillKinematicHistos(Elec,iDilepton);
	if (PassesZVeto())    {
	  FillYieldsHistograms(Elec,iZVeto, sys);      
	  if(sys==Norm) FillKinematicHistos(Elec,iZVeto);
	  if (PassesMETCut())   {
	    FillYieldsHistograms(Elec,iMET, sys);      
	    if(sys==Norm) FillKinematicHistos(Elec,iMET);
	    if (PassesNJetsCut()) {
	      FillYieldsHistograms(Elec,i2jets, sys);      
	      if(sys==Norm) FillKinematicHistos(Elec,i2jets);
	      if (PassesNBtagCut()) {
		FillYieldsHistograms(Elec,i1btag, sys);      
		if(sys==Norm) FillKinematicHistos(Elec,i1btag);
	      }
	    }
	  }	  
	}
      }
    }
  }
  ResetHypLeptons();
}
void TreeAnalysisTop::FillAnalysisTree(int flag){
#ifdef DEBUG2
  cout << "FillAnalysisTree(): Enter" << endl;
#endif
  
  ResetHypLeptons();
  
  TEvent  = T_Event_EventNumber;
  TLumi   = T_Event_LuminosityBlock;
  TRun    = T_Event_RunNumber;
    
  TSName    = gSampleName;
  TSType    = -1;  //getSampleType(sampleName);
  TSystFlag = flag;
  
  // Per Event Information
  TWeight   = EventWeight;
  TNPV      = nGoodVertex;
  TMET      = getMET();
  TMET_Phi  = getMETPhi();
  
  TNTMus    = getNTightMuons();
  TNTEls    = getNTightElectrons();
  TNMus     = getNMuons();
  TNEls     = getNElectrons();
  
  int ind1(-1),ind2(-1);
  if      (PassTriggerEMu()  && IsTTbarElMuEvent(ind1,ind2)){
    TLep0Pt   = fHypLepton1.p.Pt();
    TLep0Eta  = fHypLepton1.p.Eta();
    TLep0Phi  = fHypLepton1.p.Phi();
    TLep0Ch   = fHypLepton1.charge;
    TLep0Flav = fHypLepton1.type;   
    TLep1Pt   = fHypLepton2.p.Pt(); 
    TLep1Eta  = fHypLepton2.p.Eta();
    TLep1Phi  = fHypLepton2.p.Phi();
    TLep1Ch   = fHypLepton2.charge; 
    TLep1Flav = fHypLepton2.type;   
    TInvMass  = (fHypLepton1.p+fHypLepton2.p).M();

    // Event Classification
    if ( IsTightMuon(ind1) &&  IsTightElectron(ind2)) TTLCat = 0;
    if ( IsTightMuon(ind1) && !IsTightElectron(ind2)) TTLCat = 1;
    if (!IsTightMuon(ind1) &&  IsTightElectron(ind2)) TTLCat = 2;
    if (!IsTightMuon(ind1) && !IsTightElectron(ind2)) TTLCat = 3;

    // Jet variables:
    TNJets = getNJets();
  }
  else if (PassTriggerMuMu() && IsTTbarMuMuEvent(ind1,ind2)){
    TLep0Pt   = fHypLepton1.p.Pt();
    TLep0Eta  = fHypLepton1.p.Eta();
    TLep0Phi  = fHypLepton1.p.Phi();
    TLep0Ch   = fHypLepton1.charge;
    TLep0Flav = fHypLepton1.type;   
    TLep1Pt   = fHypLepton2.p.Pt(); 
    TLep1Eta  = fHypLepton2.p.Eta();
    TLep1Phi  = fHypLepton2.p.Phi();
    TLep1Ch   = fHypLepton2.charge; 
    TLep1Flav = fHypLepton2.type;   
    TInvMass  = (fHypLepton1.p+fHypLepton2.p).M();
    
    // Event Classification
    if ( IsTightMuon(ind1) &&  IsTightMuon(ind2)) TTLCat = 0;
    if ( IsTightMuon(ind1) && !IsTightMuon(ind2)) TTLCat = 1;
    if (!IsTightMuon(ind1) &&  IsTightMuon(ind2)) TTLCat = 2;
    if (!IsTightMuon(ind1) && !IsTightMuon(ind2)) TTLCat = 3;
  }
  else if (PassTriggerEE()   && IsTTbarElElEvent(ind1,ind2)){
    TLep0Pt   = fHypLepton1.p.Pt();
    TLep0Eta  = fHypLepton1.p.Eta();
    TLep0Phi  = fHypLepton1.p.Phi();
    TLep0Ch   = fHypLepton1.charge;
    TLep0Flav = fHypLepton1.type;   
    TLep1Pt   = fHypLepton2.p.Pt(); 
    TLep1Eta  = fHypLepton2.p.Eta();
    TLep1Phi  = fHypLepton2.p.Phi();
    TLep1Ch   = fHypLepton2.charge; 
    TLep1Flav = fHypLepton2.type;   
    TInvMass  = (fHypLepton1.p+fHypLepton2.p).M();
    
    // Event Classification
    if ( IsTightElectron(ind1) &&  IsTightElectron(ind2)) TTLCat = 0;
    if ( IsTightElectron(ind1) && !IsTightElectron(ind2)) TTLCat = 1;
    if (!IsTightElectron(ind1) &&  IsTightElectron(ind2)) TTLCat = 2;
    if (!IsTightElectron(ind1) && !IsTightElectron(ind2)) TTLCat = 3;
  }
  
  AnalysisTree->Fill();
#ifdef DEBUG2
  cout << "FillAnalysisTree(): Exit" << endl;
#endif
}
//int TreeAnalysisTop::IsDileptonEvent(int &ind1, int &ind2){
//  int res = IsLLEvent(ind1, &TreeAnalysisTop::IsTightMuon, ind2, &TreeAnalysisTop::IsTightElectron);
//  if (res > 0) return res;
//  return IsLLEvent(ind1, &TreeAnalysisTop::IsLooseMuon, ind2, &TreeAnalysisTop::IsLooseElectron);
//}
bool TreeAnalysisTop::PassesMuonEta2p1(gChannel chan){
  if (fHypLepton1.index == -1) return false;
  if (fHypLepton2.index == -1) return false;
  
  if (chan == Muon){  
    if (TMath::Abs(fHypLepton1.p.Eta()) < 2.1) return true; 
    if (TMath::Abs(fHypLepton2.p.Eta()) < 2.1) return true;
  }
  else if (chan == ElMu){
    if (TMath::Abs(fHypLepton1.p.Eta()) < 2.1) return true;
  }
  else if (chan == Elec){    
    return true;  
  }
  
  return false;
}
bool TreeAnalysisTop::Passes3rdLeptonVeto(){
  // Return false if there are not 2 signal leptons
  if (fHypLepton1.index == -1) return false;
  if (fHypLepton2.index == -1) return false;
  
  return true; // don't apply third lepton veto...
  
  //  Int_t nvetoleptons = 0;
  for(UInt_t i = 0; i < T_Muon_Pt->size(); ++i){
    if (fHypLepton1.index == i && fHypLepton1.type == 0) continue;
    if (fHypLepton2.index == i && fHypLepton2.type == 0) continue;
    if (IsVetoMuon(i)) return false;
    //    nvetoleptons++;
  }
  
  for(UInt_t i = 0; i < T_Elec_Pt->size(); ++i){
    if (fHypLepton1.index == i && fHypLepton1.type == 1) continue;
    if (fHypLepton2.index == i && fHypLepton2.type == 1) continue;
    if (IsVetoElectron(i)) return false;
    //    nvetoleptons++;
  }
  //  if (nvetoleptons > 0) return false;
  
  return true;
}
bool TreeAnalysisTop::PassesMllVeto(){
  // Check consistency.
  if (fHypLepton1.index == -1) return false;
  if (fHypLepton2.index == -1) return false;
  
  float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
  // remove low resonances
  if (InvMass < 20.)            return false; 

  return true;
}
bool TreeAnalysisTop::PassesZVeto(){
  // Check consistency.
  if (fHypLepton1.index == -1) return false;
  if (fHypLepton2.index == -1) return false;
  
  float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
  // remove Z-peak
  if (InvMass > 76. && InvMass < 106.) return false;

  return true;
}
bool TreeAnalysisTop::PassesTopDCut(){
  // Check consistency.
  if (fHypLepton1.index == -1) return false;
  if (fHypLepton2.index == -1) return false;
  
  if (getTopD() < 0.8) return false;
  return true;
}
bool TreeAnalysisTop::PassesNJetsCut(){
  if (getNJets() <= 1) return false;

  return true;
}
bool TreeAnalysisTop::PassesMETCut(){
  if (getMET() < 40.) return false;
  
  return true;
}
bool TreeAnalysisTop::PassesNBtagCut(){
  if (getNBTags() < 1) return false;
  
  return true;
}

bool TreeAnalysisTop::IsTTbarMuMuEvent(int &ind1, int &ind2){
#ifdef DEBUG2
  cout << "IsMuMuEvent()?" << endl;
#endif
  int nmus = HasLooseMuons(ind1, ind2);
  if (nmus < 2) return false;
  
  // pick dilepton pair...
  if (TMath::Abs(IsDileptonEvent(ind1,ind2)) != 1) return false;
  
  // Define Hypothesis Leptons...
  SetHypLepton1(ind1, Muon);
  SetHypLepton2(ind2, Muon);
  
  if (!PassesZVeto())    return false;
  if (!PassesNJetsCut()) return false;
  if (!PassesMETCut())   return false;
  if (!PassesNBtagCut()) return false;
  
#ifdef DEBUG2
  cout << "IsMuMuEvent()?: YES" << endl;
#endif

  return true;
}
bool TreeAnalysisTop::IsTTbarElElEvent(int &ind1, int &ind2){
#ifdef DEBUG2
  cout << "IsElElEvent()?" << endl;
#endif
  int nels = HasLooseElectrons(ind1, ind2);
  if (nels < 2) return false;
  
  // pick dilepton pair...
  if (TMath::Abs(IsDileptonEvent(ind1,ind2)) != 2) return false;
  
  // Define Hypothesis Leptons...
  SetHypLepton1(ind1, Elec);
  SetHypLepton2(ind2, Elec);
  
  if (!PassesZVeto()) return false;
  if (!PassesNJetsCut())       return false;
  if (!PassesMETCut())         return false;
  if (!PassesNBtagCut())       return false;
  
#ifdef DEBUG2
  cout << "IsElElEvent()?: YES" << endl;
#endif
  return true;
}
bool TreeAnalysisTop::IsTTbarElMuEvent(int &ind1, int &ind2){ 
#ifdef DEBUG2
  cout << "IsElMuEvent()?" << endl;
#endif
  int nmus = HasLooseMuons(ind1, ind2);
  int nels = HasLooseElectrons(ind2, ind1);
  if (nmus < 1 || nels < 1) return false;
  
  if (TMath::Abs(IsDileptonEvent(ind1,ind2)) != 3) return false;
  
  SetHypLepton1(ind1, Muon);
  SetHypLepton2(ind2, Elec);

  //  if (!PassesZVeto()) return false;
  if (!PassesNJetsCut())       return false;
  if (!PassesNBtagCut())       return false;
    
#ifdef DEBUG2
  cout << "IsElMuEvent()?: YES" << endl;
#endif
  return true;
}
bool TreeAnalysisTop::IsElMuEvent(int &ind1, int &ind2){
  // if fChargeSwitch, pick SS event.
  if (fChargeSwitch){      return (IsDileptonEvent(ind1,ind2)  == 3);   }

  return (IsDileptonEvent(ind1,ind2) == -3);
}
bool TreeAnalysisTop::IsMuMuEvent(int &ind1, int &ind2){
  // if fChargeSwitch, pick SS event.
  if (fChargeSwitch){  return (IsDileptonEvent(ind1,ind2)  == 1); }
  
  return (IsDileptonEvent(ind1,ind2) == -1);
}
bool TreeAnalysisTop::IsElElEvent(int &ind1, int &ind2){
  // if fChargeSwitch, pick SS event.
  if (fChargeSwitch){    return (IsDileptonEvent(ind1,ind2)  == 2); }
  
  return (IsDileptonEvent(ind1,ind2) == -2);
}
int TreeAnalysisTop::IsDileptonEvent(int &ind1, int &ind2){
  // bool(TreeAnalysisTop::*muonSelector)(unsigned int), 
  // int &ind2, bool(TreeAnalysisTop::*eleSelector)(unsigned int)){
  // Looks for a pair of leptons with given object selectors
  // Return the channel: 0 = none found
  //                1 / -1 = mu+mu+ (SS) / mu-mu+ (OS) pair
  //                2 / -2 = e+e+   (SS) / e-e+   (OS) pair
  //                3 / -3 = mu+e+  (SS) / mu-e+  (OS) pair
  // The indices in the argument given are sorted by pt unless
  
  vector<lepton> tmp_Loose;
  vector<lepton> tmp_Tight;

  // First store all loose leptons in two vectors according to their charges
  TLorentzVector plep;
  for(UInt_t i = 0; i < T_Muon_Pt->size(); ++i){
    if (IsLooseMuon(i) == false) continue;
    plep.SetPxPyPzE(MuPx.at(i), MuPy.at(i), T_Muon_Pz->at(i), T_Muon_Energy->at(i));
    lepton tmpLepton(plep, T_Muon_Charge->at(i), 0, i);
    tmp_Loose.push_back(tmpLepton);
    
    if (IsTightMuon(i) == false) continue;
    tmp_Tight.push_back(tmpLepton);
  }
  for(UInt_t i = 0; i < T_Elec_Pt->size(); ++i){
    if (IsLooseElectron(i) == false) continue;
    plep.SetPxPyPzE(ElPx.at(i), ElPy.at(i), T_Elec_Pz->at(i), T_Elec_Energy->at(i));
    lepton tmpLepton(plep, T_Elec_Charge->at(i), 1, i);
    tmp_Loose.push_back(tmpLepton);
    
    if (IsTightElectron(i) == false) continue;
    tmp_Tight.push_back(tmpLepton);
  }
  
  // Check for at least one loose pair
  if(tmp_Loose.size() < 2) return 0;
  
  /////////////////////////////////////////////////////////////////////////
  // Sort these vectors by  pt, if there are not two tight leptons use two
  // loose ones.
  vector<lepton> Leptons;
  if (tmp_Tight.size() > 1) Leptons = SortLeptonsByPt(tmp_Tight);
  else                      Leptons = SortLeptonsByPt(tmp_Loose); 
  
  // Proceed to select pair with highest pt
  vector<lepton> selectedPair;
  selectedPair.push_back(Leptons[0]);
  selectedPair.push_back(Leptons[1]);
  int select = Leptons[0].charge*Leptons[1].charge;
  
  /////////////////////////////////////////////////////////////////////////
  int result = 0;
  if (selectedPair[0].type == 0 && selectedPair[1].type == 0) result = 1; // mu/mu
  if (selectedPair[0].type == 1 && selectedPair[1].type == 1) result = 2; // el/el
  if (selectedPair[0].type == 0 && selectedPair[1].type == 1) result = 3; // mu/el
  if (selectedPair[0].type == 1 && selectedPair[1].type == 0) result = 3; // mu/el
  
  // Return values, assigning indexes (first is always a muon, second electron):
  if (result == 3) {
    if      (selectedPair[0].type == 0 && selectedPair[1].type == 1) { 
      ind1 = selectedPair[0].index;
      ind2 = selectedPair[1].index;
    }
    else if (selectedPair[0].type == 1 && selectedPair[1].type == 0){
      ind1 = selectedPair[1].index;
      ind2 = selectedPair[0].index;
    }
  }
  else {
    ind1 = selectedPair[0].index;
    ind2 = selectedPair[1].index;
  } 
  
  result *= select; // Add charge to result
  
  return result;
}
void TreeAnalysisTop::FillTLRatios(){
  int looseMuInd(-1);
  if(PassSingleMuTrigger() && IsSigSupMuEvent(looseMuInd)){
    if( IsTightMuon(looseMuInd) ){
      tlratios[0].fntight   ->Fill(T_Muon_Pt->at(looseMuInd), TMath::Abs(T_Muon_Eta->at(looseMuInd)), EventWeight);
      tlratios[0].fntight_nv->Fill(nGoodVertex,                                     EventWeight);
    }
    if( IsLooseMuon(looseMuInd) ){
      tlratios[0].fratio_pt ->Fill(IsTightMuon(looseMuInd), T_Muon_Pt->at(looseMuInd));
      tlratios[0].fratio_eta->Fill(IsTightMuon(looseMuInd), TMath::Abs(T_Muon_Eta->at(looseMuInd)));
      tlratios[0].fratio_nv ->Fill(IsTightMuon(looseMuInd), nGoodVertex);
      
      tlratios[0].fnloose   ->Fill(T_Muon_Pt->at(looseMuInd), TMath::Abs(T_Muon_Eta->at(looseMuInd)), EventWeight);
      tlratios[0].fnloose_nv->Fill(nGoodVertex,                                     EventWeight);
    }
  }
  // ZMuMu Control Region
  int mu1(-1), mu2(-1);
  if(PassTriggerMuMu() && IsZMuMuEvent(mu1, mu2)){
    if( IsTightMuon(mu2) ){
      tlratios[0].pntight   ->Fill(T_Muon_Pt->at(mu2), TMath::Abs(T_Muon_Eta->at(mu2)), EventWeight);
      tlratios[0].pntight_nv->Fill(nGoodVertex,                       EventWeight);
    }
    if( IsLooseMuon(mu2) ){
      tlratios[0].pratio_pt ->Fill(IsTightMuon(mu2), T_Muon_Pt->at(mu2));
      tlratios[0].pratio_eta->Fill(IsTightMuon(mu2), TMath::Abs(T_Muon_Eta->at(mu2)));
      tlratios[0].pratio_nv ->Fill(IsTightMuon(mu2), nGoodVertex);
      
      tlratios[0].pnloose->Fill(T_Muon_Pt->at(mu2), TMath::Abs(T_Muon_Eta->at(mu2)), EventWeight);
      tlratios[0].pnloose_nv->Fill(nGoodVertex,                    EventWeight);
    }
  }
  ResetHypLeptons();
  
  // Electron QCD Control region
  int looseElInd(-1);
  if(PassSingleElTrigger() && IsSigSupElEvent(looseElInd)){
    if( IsTightElectron(looseElInd) ){
      tlratios[1].fntight   ->Fill(T_Elec_Pt->at(looseElInd), TMath::Abs(T_Elec_Eta->at(looseElInd)), EventWeight);
      tlratios[1].fntight_nv->Fill(nGoodVertex,                                         EventWeight);
    }
    if( IsLooseElectron(looseElInd) ){
      tlratios[1].fratio_pt ->Fill(IsTightElectron(looseElInd), T_Elec_Pt->at(looseElInd));
      tlratios[1].fratio_eta->Fill(IsTightElectron(looseElInd), TMath::Abs(T_Elec_Eta->at(looseElInd)));
      tlratios[1].fratio_nv ->Fill(IsTightElectron(looseElInd), nGoodVertex);
      
      tlratios[1].fnloose   ->Fill(T_Elec_Pt->at(looseElInd), TMath::Abs(T_Elec_Eta->at(looseElInd)), EventWeight);
      tlratios[1].fnloose_nv->Fill(nGoodVertex,                                         EventWeight);
    }
  }
  
  int el1(-1), el2(-1);
  if(PassTriggerEE() && IsZElElEvent(el1, el2)){
    if( IsTightElectron(el2) ){
      tlratios[1].pntight   ->Fill(T_Elec_Pt->at(el2), TMath::Abs(T_Elec_Eta->at(el2)), EventWeight);
      tlratios[1].pntight_nv->Fill(nGoodVertex,                       EventWeight);
    }
    if( IsLooseElectron(el2) ){
      tlratios[1].pratio_pt ->Fill(IsTightElectron(el2), T_Elec_Pt->at(el2));
      tlratios[1].pratio_eta->Fill(IsTightElectron(el2), TMath::Abs(T_Elec_Eta->at(el2)));
      tlratios[1].pratio_nv ->Fill(IsTightElectron(el2), nGoodVertex);
      
      tlratios[1].pnloose   ->Fill(T_Elec_Pt->at(el2), TMath::Abs(T_Elec_Eta->at(el2)), EventWeight);
      tlratios[1].pnloose_nv->Fill(nGoodVertex,                       EventWeight);
    }
  }
}
//==============================================================================
// LEPTON SELECTORS
//==============================================================================
bool momentumComparator(lepton i, lepton j){ return (i.p.Pt()>j.p.Pt()); }
vector<lepton> TreeAnalysisTop::SortLeptonsByPt(vector<lepton>& leptons){
  vector<lepton> theLep = leptons;
  sort (theLep.begin(), theLep.end(), momentumComparator);
  return theLep;
}
//------------------------------------------------------------------------------
// Muon Selectors
//------------------------------------------------------------------------------
int  TreeAnalysisTop::getNTightMuons(){
  int nMus = 0;
  for (UInt_t i=0; i<T_Muon_Energy->size(); i++) {
    if (!IsTightMuon(i)) nMus++;
  }
  return nMus;
}
int  TreeAnalysisTop::getNMuons(){
  int nMus = 0;
  for (UInt_t i=0; i<T_Muon_Energy->size(); i++) {
    if (!IsGoodMuon(i,10)) nMus++;
  }
  return nMus;
}
bool TreeAnalysisTop::IsGoodMuon(unsigned int iMuon, float ptcut){
  if (T_Muon_Pt->at(iMuon) < ptcut)      return false;
  if (TMath::Abs(T_Muon_Eta->at(iMuon)) > 2.4) return false;
  
  return true;
}
bool TreeAnalysisTop::IsVetoMuon(unsigned int iMuon){
  if (T_Muon_Pt->at(iMuon) < 20)               return false;
  if (TMath::Abs(T_Muon_Eta->at(iMuon)) > 2.4) return false;
  
  if (T_Muon_IsGlobalMuon->at(iMuon) == 0 && T_Muon_IsTrackerMuonArbitrated->at(iMuon) == 0) return false;
  
  float relIso = (T_Muon_chargedHadronIsoR04->at(iMuon) + max(0.0 , T_Muon_neutralHadronIsoR04->at(iMuon) + T_Muon_photonIsoR04->at(iMuon)- 0.5*T_Muon_sumPUPtR04->at(iMuon)))/T_Muon_Pt->at(iMuon);
  if (relIso > 0.20) return false;

  return true;
}
bool TreeAnalysisTop::IsLooseMuon(unsigned int iMuon){
  if (!IsGoodMuon(iMuon)) return false;
  
  // POG Tight Muons definition				       
  if (!T_Muon_IsGlobalMuon->at(iMuon))                       return false;
  if (T_Muon_NormChi2GTrk->at(iMuon) >= 10.)                 return false;
  if (T_Muon_NValidHitsGTrk->at(iMuon) < 1)                  return false;
  //this is still not the exact same def.		       
  if (T_Muon_NumOfMatchedStations->at(iMuon) <= 1)           return false; 
  //							       
  if (TMath::Abs(T_Muon_IPwrtAveBSInTrack->at(iMuon)) >= 0.2)      return false; 
  if (TMath::Abs(T_Muon_vz->at(iMuon) - T_Vertex_z->at(0)) >= 0.5) return false;
  if (T_Muon_NValidPixelHitsInTrk->at(iMuon) == 0)           return false;
  if (T_Muon_NLayers->at(iMuon) <= 5)                        return false;
  
  float relIso = (T_Muon_chargedHadronIsoR04->at(iMuon) + max(0.0 , T_Muon_neutralHadronIsoR04->at(iMuon) + T_Muon_photonIsoR04->at(iMuon)- 0.5*T_Muon_sumPUPtR04->at(iMuon)))/T_Muon_Pt->at(iMuon);
  
  if (relIso > 1.0) return false;
  
  return true;
}
bool TreeAnalysisTop::IsTightMuon(unsigned int iMuon){
  if (!IsLooseMuon(iMuon)) return false;
  
  float relIso = getMuonIso(iMuon);
  
  if (relIso > 0.12) return false;
  
  return true;
}
float TreeAnalysisTop::getMuonIso(int iMuon){
  if (iMuon < 0) return 9999.;
  if (iMuon >= (int)T_Muon_chargedHadronIsoR04->size()) return 9999.;

  return (T_Muon_chargedHadronIsoR04->at(iMuon) + max(0.0 , T_Muon_neutralHadronIsoR04->at(iMuon) + T_Muon_photonIsoR04->at(iMuon)- 0.5*T_Muon_sumPUPtR04->at(iMuon)))/T_Muon_Pt->at(iMuon);
}
//------------------------------------------------------------------------------
// Electron Selectors
//------------------------------------------------------------------------------
//SANTIint TreeAnalysisTop::GetSelectedElecInd()
//SANTI{
//SANTI  if (nSelElec > 0) {
//SANTI    cout << "[ERROR] GetSelectedElec(): This function has already been called" << endl;
//SANTI    cout << "[ERROR] GetSelectedElec(): Cleaning previous information        " << endl;
//SANTI    S_Elec.clear();
//SANTI  }
//SANTI  
//SANTI  // Loop over all electrons and keep loose electrons
//SANTI  for (UInt_t i=0; i<T_Elec_Energy->size(); i++) {
//SANTI    if (!IsLooseElectron(i)) continue;
//SANTI    S_Elec.push_back(i);
//SANTI  }
//SANTI  
//SANTI  return S_Elec.size();
//SANTI}
int  TreeAnalysisTop::getNTightElectrons(){
  int nEls = 0;
  for (UInt_t i=0; i<T_Elec_Energy->size(); i++) {
    if (IsTightElectron(i)) nEls++;
  }
  return nEls;
}
int  TreeAnalysisTop::getNElectrons(){
  int nEls = 0;
  for (UInt_t i=0; i<T_Elec_Energy->size(); i++) {
    if (IsGoodElectron(i,10)) nEls++;
  }
  return nEls;
}
bool TreeAnalysisTop::IsGoodElectron(unsigned int iElec, float ptcut){
  if (T_Elec_Pt->at(iElec) < ptcut)            return false;
  if (TMath::Abs(T_Elec_Eta->at(iElec)) > 2.5) return false;
  
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if (sceta > 1.4442 && sceta < 1.566) return false;

  return true;
}
bool TreeAnalysisTop::IsLooseElectron(unsigned int iElec){
  if (!IsGoodElectron(iElec)) return false;
  
  float pt  = T_Elec_Pt->at(iElec);
  
  // Require electrons passing Trigger requirements
  bool passTriggerID = false;
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(T_Elec_sigmaIetaIeta->at(iElec)           < 0.014 &&
       T_Elec_HtoE->at(iElec)                    < 0.15  &&
       T_Elec_dr03TkSumPt->at(iElec)/pt          < 0.2   &&
       T_Elec_dr03EcalSumEt->at(iElec)/pt        < 0.2   &&
       T_Elec_dr03HcalSumEt->at(iElec)/pt        < 0.2   &&
       T_Elec_nLost->at(iElec)                    == 0     )
      passTriggerID = true;
  }
  else {
    if(T_Elec_sigmaIetaIeta->at(iElec)           < 0.035 &&
       T_Elec_HtoE->at(iElec)                    < 0.10  &&
       T_Elec_dr03TkSumPt->at(iElec)/pt          < 0.2   &&
       T_Elec_dr03EcalSumEt->at(iElec)/pt        < 0.2   &&
       T_Elec_dr03HcalSumEt->at(iElec)/pt        < 0.2   &&
       T_Elec_nLost->at(iElec)                    == 0     )
      passTriggerID = true;
  }  
  
  return passTriggerID;
}
bool TreeAnalysisTop::IsVetoElectron(unsigned int iElec){
  if (T_Elec_Pt->at(iElec) < 20)               return false;
  if (TMath::Abs(T_Elec_Eta->at(iElec)) > 2.5) return false;
  
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if (sceta > 1.4442 && sceta < 1.566) return false;
  
  // implement veto electron
  bool passVetoID = false;
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.007 &&
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.8   &&
       T_Elec_sigmaIetaIeta->at(iElec)                      < 0.01  &&
       T_Elec_HtoE->at(iElec)                               < 0.15  &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.04  &&
       TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.2)
      passVetoID = true;
  }
  else {
    if(TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.01  &&
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.7   &&
       T_Elec_sigmaIetaIeta->at(iElec)                      < 0.03  &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.04  &&
       TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.2)
      passVetoID = true;
  }
  
  if (passVetoID == false) return false;
  
  // ISOLATION
  float EA     = getEACorrection(T_Elec_Eta->at(iElec));
  float relIso = (T_Elec_chargedHadronIso->at(iElec) + max((float)0.0, T_Elec_neutralHadronIso->at(iElec) + T_Elec_photonIso->at(iElec) - T_Event_RhoIso*EA))/T_Elec_Pt->at(iElec);
  
  if (relIso > 0.15) return false;
  
  return true;
}
bool TreeAnalysisTop::IsMVAIDElectron(unsigned int iElec){
  // MVA ID
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if      (sceta < 0.8                    ) { if (T_Elec_MVA->at(iElec) > 0.94) return true; }
  else if (sceta >= 0.8   && sceta < 1.479) { if (T_Elec_MVA->at(iElec) > 0.85) return true; }
  else if (sceta >= 1.479 && sceta < 2.5  ) { if (T_Elec_MVA->at(iElec) > 0.92) return true; }
  
  return false;
}
bool TreeAnalysisTop::IsTightElectron(unsigned int iElec){
  if (!IsLooseElectron(iElec))               return false;

  if (!IsMVAIDElectron(iElec))               return false;
  if (!T_Elec_passConversionVeto->at(iElec)) return false;
  if ( T_Elec_nHits->at(iElec) >= 1)         return false;
  
  // ISOLATION
  float relIso =  getElecIso(iElec);

  if (relIso > 0.15) return false;
  
  return true;
}
float TreeAnalysisTop::getElecIso(int iElec){
  if (iElec < 0) return 9999.;
  if (iElec >= (int)T_Elec_chargedHadronIso->size()) return 9999.;

  float EA     = getEACorrection(T_Elec_Eta->at(iElec));
  float relIso = (T_Elec_chargedHadronIso->at(iElec) + 
		  max((float)0.0, 
		      T_Elec_neutralHadronIso->at(iElec) + 
		      T_Elec_photonIso->at(iElec) - 
		      T_Event_RhoIso*EA
		      )
		  )/T_Elec_Pt->at(iElec);

  return relIso;
}
float TreeAnalysisTop::getEACorrection(float eta){  // for a 0.3 CONE
  float abseta = TMath::Abs(eta);
  
  if      (abseta < 1.0)                      return 0.13; // +/- 0.001
  else if (abseta >= 1.0   && abseta < 1.479) return 0.14; // +/- 0.002
  else if (abseta >= 1.479 && abseta < 2.0)   return 0.07; // +/- 0.001
  else if (abseta >= 2.0   && abseta < 2.2)   return 0.09; // +/- 0.001
  else if (abseta >= 2.2   && abseta < 2.3)   return 0.11; // +/- 0.002
  else if (abseta >= 2.3   && abseta < 2.4)   return 0.11; // +/- 0.003
  else if (abseta >= 2.4)                     return 0.14; // +/- 0.004
  
  cout << "[ERROR] getEACorrection(): No correction factor found!!" << endl;
  return -999999.;
}
void TreeAnalysisTop::setMET(float newmet){
  MET = newmet;
}
float TreeAnalysisTop::getMET(){
  return MET; 
}
float TreeAnalysisTop::getMETPhi(){
  return MET_Phi;
}
int TreeAnalysisTop::getNJets(){
  int nj(0);
  for (unsigned int i=0; i<T_JetAKCHS_Energy->size(); i++) if(IsGoodJet(i,gJetEtCut)) nj++;
  
  return nj;
}
int TreeAnalysisTop::getLeadingJet(){
  for (unsigned int i=0; i<T_JetAKCHS_Energy->size(); i++) if(IsGoodJet(i,gJetEtCut)) return i;
  return -1;
}
int TreeAnalysisTop::getSecondLeadingJet(){
  for (unsigned int i=0; i<T_JetAKCHS_Energy->size(); i++) 
    if(i!=getLeadingJet() && IsGoodJet(i,gJetEtCut)) return i;
  return -1;
}

float TreeAnalysisTop::getDRClosestJet(TLorentzVector lep){
  TLorentzVector jet;
  float minDR = 9999.;
  for (unsigned int i=0; i<T_JetAKCHS_Energy->size(); i++) {
    if(!IsGoodJet(i,gJetEtCut)) continue;
    jet.SetPxPyPzE(T_JetAKCHS_Px->at(i),T_JetAKCHS_Py->at(i),T_JetAKCHS_Pz->at(i),T_JetAKCHS_Energy->at(i));
    if (minDR > lep.DeltaR(jet)) minDR = lep.DeltaR(jet);
  }
  
  return minDR;
}
float TreeAnalysisTop::getDPhiClosestJet(TLorentzVector lep){
  TLorentzVector jet;
  float minDphi = 9999.;
  for (unsigned int i=0; i<T_JetAKCHS_Energy->size(); i++) {
    if(!IsGoodJet(i,gJetEtCut)) continue;
    jet.SetPxPyPzE(T_JetAKCHS_Px->at(i),T_JetAKCHS_Py->at(i),T_JetAKCHS_Pz->at(i),T_JetAKCHS_Energy->at(i));
    if (minDphi > TMath::Abs(lep.DeltaPhi(jet))) minDphi = TMath::Abs(lep.DeltaPhi(jet));
  }
  
  return minDphi;
}
int TreeAnalysisTop::getLeadingJetbTag(){
  for (unsigned int i=0; i<T_JetAKCHS_Energy->size(); i++) {
    if(!IsGoodJet(i,gJetEtCut)) continue;
    if(T_JetAKCHS_Tag_CombSVtx->at(i) > 0 && T_JetAKCHS_Tag_CombSVtx->at(i) < 1) return i;
  }
  
  return  -1;
}
int TreeAnalysisTop::getNBTags(){
  int ntags(0);
  
  int btagSys = 0;
  if (gSysSource == BtagUp)   btagSys =  1;
  if (gSysSource == BtagDown) btagSys = -1;
  
  for(UInt_t i = 0; i <T_JetAKCHS_Energy->size(); i++) {
    if(!IsGoodJet(i,gJetEtCut)) continue;
    
    if(gIsData) {
      if(fBTagSF->IsTagged(T_JetAKCHS_Tag_CombSVtx->at(i),-999999, JetEt.at(i), 
			   T_JetAKCHS_Eta->at(i), btagSys)) ntags++;
    }
    else {
      if(fBTagSF->IsTagged(T_JetAKCHS_Tag_CombSVtx->at(i),T_JetAKCHS_Parton_Flavour->at(i), 
			   JetEt.at(i), T_JetAKCHS_Eta->at(i), btagSys)) ntags++;
    }
  }
  return ntags;
}
float TreeAnalysisTop::getDeltaPhillJet(){
  if (fHypLepton1.index == -1) return -999.;
  if (fHypLepton2.index == -1) return -999.;

  Int_t ij = getLeadingJetbTag();
  TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;
  TLorentzVector jet;
  if (ij < 0) return -999.;
  jet.SetPxPyPzE(T_JetAKCHS_Px->at(ij),T_JetAKCHS_Py->at(ij),T_JetAKCHS_Pz->at(ij),T_JetAKCHS_Energy->at(ij));
  
  return TMath::Abs(dilep.DeltaPhi(jet));
}
float TreeAnalysisTop::getTopD(){
  if (fHypLepton1.index == -1) return -999;
  if (fHypLepton2.index == -1) return -999;

  // Make Dilepton system
  TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;

  Float_t DeltaPhi(0.),TopD(0.);
  TLorentzVector jet;
  Int_t ij = getLeadingJet();
  if (ij < 0) return -999.;
  jet.SetPxPyPzE(T_JetAKCHS_Px->at(ij),T_JetAKCHS_Py->at(ij),T_JetAKCHS_Pz->at(ij),T_JetAKCHS_Energy->at(ij));
   
  DeltaPhi = TMath::Abs(dilep.DeltaPhi(jet));
  //  DeltaPhi = TMath::Pi() - TMath::Abs(TMath::Abs(dilep.Phi() - jet.Phi()) - TMath::Pi());
  TopD     = 1 - (DeltaPhi/TMath::Pi()) * (1 - T_JetAKCHS_Tag_CombSVtx->at(ij));
  
  return TopD;
}
bool TreeAnalysisTop::IsGoodJet(unsigned int ijet, float ptcut){
  float minDR = 0.4;
   
  if(JetEt.at(ijet) < ptcut)      return false;
  if(TMath::Abs(T_JetAKCHS_Eta->at(ijet)) > 2.5) return false; // btagging only up to 2.4
  
  // JetID 
  if ( !(T_JetAKCHS_nDaughters->at(ijet)        > 1   ) ) return false;
  if ( !(T_JetAKCHS_NeutHadEnergyFrac->at(ijet) < 0.99) ) return false;
  if ( !(T_JetAKCHS_NeutEmEnergyFrac ->at(ijet) < 0.99) ) return false;
  if (TMath::Abs(T_JetAKCHS_Eta->at(ijet)) < 2.5){
    if ( !(T_JetAKCHS_CharEmEnergyFrac->at(ijet)  < 0.99) ) return false;
    if ( !(T_JetAKCHS_CharHadEnergyFrac->at(ijet) > 0.00) ) return false;
    if ( !(T_JetAKCHS_ChargedMultiplicity->at(ijet) > 0 ) ) return false;
  }
  
  // Remove jets close to hypothesis leptons
  TLorentzVector jet(T_JetAKCHS_Px->at(ijet),
		     T_JetAKCHS_Py->at(ijet),
		     T_JetAKCHS_Pz->at(ijet),
		     T_JetAKCHS_Energy->at(ijet));
  
  if(fHypLepton1.index > -1) if(jet.DeltaR(fHypLepton1.p) < minDR) return false;
  if(fHypLepton2.index > -1) if(jet.DeltaR(fHypLepton2.p) < minDR) return false;
  
  // Remove jets close to all tight leptons
  for(unsigned int imu = 0; imu < T_Muon_Energy->size(); ++imu){
    if(!IsTightMuon(imu)) continue;
    TLorentzVector mu(T_Muon_Px    ->at(imu),     
		      T_Muon_Py    ->at(imu),     
		      T_Muon_Pz    ->at(imu),     
		      T_Muon_Energy->at(imu));
    if(jet.DeltaR(mu) > minDR) continue;
    return false;
  }
  for(unsigned int iel = 0; iel < T_Elec_Energy->size(); ++iel){
    if(!IsTightElectron(iel)) continue;
    TLorentzVector el(T_Elec_Px    ->at(iel),     
		      T_Elec_Py    ->at(iel),     
		      T_Elec_Pz    ->at(iel),     
		      T_Elec_Energy->at(iel));
    
    if(jet.DeltaR(el) > minDR) continue;
    return false;
  }
  return true;
}

//------------------------------------------------------------------------------
// SelectedGenLepton
//------------------------------------------------------------------------------
void TreeAnalysisTop::SelectedGenLepton()
{

#ifdef __ISMC
  // Count generated muons and electrons
  //----------------------------------------------------------------------------
  nGenElec = T_Gen_ElecSt3_PID->size();
  nGenMuon = T_Gen_MuonSt3_PID->size();
  nGenTau  = T_Gen_TauSt3_PID->size();
  nGenLepton = nGenElec + nGenMuon + nGenTau;
  nTauElec = 0;
  nTauMuon = 0;

  for (UInt_t i=0; i<T_Gen_TauSt3_PID->size(); i++) {
    if (T_Gen_TauSt3_IsLepDec->at(i)) {
      if (abs(T_Gen_TauSt3_LepDec_PID->at(i)) == 11) nTauElec++;
      if (abs(T_Gen_TauSt3_LepDec_PID->at(i)) == 13) nTauMuon++;
    }
  }

  //  Including taus as Signal
  if ( nGenLepton == 2 ) {
    if ( nGenElec == 2 || (nGenElec == 1 && nTauElec == 1) || nTauElec == 2 ) nGenElec = 2; 
    if ( nGenMuon == 2 || (nGenMuon == 1 && nTauMuon == 1) || nTauMuon == 2 ) nGenMuon = 2; 
    if ( (nGenElec == 1 && nGenMuon == 1) || (nGenElec == 1 && nTauMuon == 1) ||
	 (nGenMuon == 1 && nTauElec == 1) || (nTauElec == 1 && nTauMuon == 1) ) { nGenElec = 1; nGenMuon = 1; }
  } 
#endif
}
//SANTIint TreeAnalysisTop::muIndexToBin(int ind){
//SANTI  // return the bin to fill for each id/type
//SANTI  int id    = abs(T_Muon_PID->at(ind));
//SANTI  int mid   = abs(T_Muon_GMID->at(ind));
//SANTI  int mtype = abs(MuGenMType[ind]);
//SANTI  if(id  != 13)                                 return 1; // mis id
//SANTI  if(mtype == 1)                                return 2; // W/Z skipped in madgraph event
//SANTI  if(mid == 24)                                 return 2; // W
//SANTI  if(mid == 23)                                 return 3; // Z
//SANTI  if(mtype == 2)                                return 4; // tau
//SANTI  if(mtype == 11 || mtype == 12 || mtype == 18) return 5; // light hadrons
//SANTI  if(mtype == 13 || mtype == 19)                return 6; // strange hadrons
//SANTI  if(mtype == 14 || mtype == 16 || mtype == 20) return 7; // charmed hadrons
//SANTI  if(mtype == 15 || mtype == 17 || mtype == 21) return 8; // bottom hadrons
//SANTI  if(mtype == 91 || mtype == 92)                return 9; // pythia strings
//SANTI  return 12;                                              // uid
//SANTI}
//SANTIint TreeAnalysisTop::elIndexToBin(int ind){
//SANTI  // For the origin histograms
//SANTI  // return the bin to fill for each id/type
//SANTI  int id    = abs(ElGenID[ind]);
//SANTI  int type  = abs(ElGenType[ind]);
//SANTI  int mid   = abs(ElGenMID[ind]);
//SANTI  int mtype = abs(ElGenMType[ind]);
//SANTI  if(id  != 11){                                 // mis id
//SANTI    if(type == 0 || type == 2)                 return 1;  // mis-match
//SANTI    if(id == 22)                               return 2;  // gamma
//SANTI    if(type == 11 || type == 12 || type == 13 ||
//SANTI       type == 18 || type == 19)               return 3;  // Hadr. fake
//SANTI    return 12;                                            // uid
//SANTI  }
//SANTI  if(mtype == 1)                                 return 4;  // W/Z skipped in madgraph event
//SANTI  if(mid == 24)                                  return 4;  // W
//SANTI  if(mid == 23)                                  return 5;  // Z
//SANTI  if(mtype == 2)                                 return 6;  // tau
//SANTI  if(mtype == 11 || mtype == 12 || mtype == 18)  return 7;  // light hadrons
//SANTI  if(mtype == 13 || mtype == 19)                 return 8;  // strange hadrons
//SANTI  if(mtype == 91 || mtype == 92)                 return 11; // pythia strings
//SANTI  return 12;                                                // uid
//SANTI}
void TreeAnalysisTop::propagateMET(TLorentzVector nVec, TLorentzVector oVec){
  TLorentzVector met;
  met.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
  // set the pfMET to the old MET minus original vector plus new vector 
  setMET( (met+oVec-nVec).Pt() );
}
std::vector<int> TreeAnalysisTop::CleanedJetIndices(float pt){
  std::vector<int> cleanJetsInd;
  
  for(UInt_t i = 0; i < T_JetAKCHS_Energy->size(); i++){
    if (IsGoodJet(i,pt)) cleanJetsInd.push_back(i);
  }
  return cleanJetsInd;
}
void TreeAnalysisTop::SmearJetPts(int flag){
  // Modify the jet pt for systematics studies
  // Either shifted or smeared
  // propagate to the MET!!
  if(gIsData)   return; // don't smear data
  if(flag == 0) return; // 0 makes no sense
  // select the jets you want to have... 
  
  std::vector<int> cleanJets = CleanedJetIndices(15.);
  TLorentzVector ojets, jets, tmp;                            // 4-vec of old jets, newjets and a tmp-vector
  
  std::vector<int>::const_iterator it = cleanJets.begin();
  
  for( it = cleanJets.begin(); it != cleanJets.end(); ++it) {
    tmp.SetPxPyPzE(T_JetAKCHS_Px->at(*it), T_JetAKCHS_Py->at(*it), 
		   T_JetAKCHS_Pz->at(*it), T_JetAKCHS_Energy->at(*it));            // set temp to the jet
    float phi = tmp.Phi();
    ojets += tmp;                                                                  // add jet to the old jets vector
    if(flag == 1) JetEt.at(*it) *= (1 + T_JetAKCHS_Uncertainty->at(*it)); // vary up for flag 1
    if(flag == 2) JetEt.at(*it) *= (1 - T_JetAKCHS_Uncertainty->at(*it)); // vary down for flag 2;
    if(flag == 3){
      // get the resolution
      float sigmaMC  = getErrPt(JetEt.at(*it), T_JetAKCHS_Eta->at(*it))/JetEt.at(*it);
      float jerScale = getJERScale(*it);                                  // get JER scale factors 
      
      float factor = fRand3->Gaus(1., sqrt(jerScale*jerScale -1.)*sigmaMC );
      JetEt.at(*it) = JetEt.at(*it) * factor;           // smear for flag 3
    }
    // set tmp to the scaled/smeared jet
    tmp.SetPtEtaPhiE(JetEt.at(*it), T_JetAKCHS_Eta->at(*it), phi, T_JetAKCHS_Energy->at(*it)); 
    jets += tmp;                                                    // add scaled/smeared jet to the new jets
  }
  propagateMET(jets, ojets);                                        // propagate this change to the MET
}
void TreeAnalysisTop::ScaleLeptons(int flag){
  // Shift the lepton pts for systematics studies
  if(gIsData) return; // don't smear data
  
  if(flag == 0) return;
  float scale = 0.002;
  TLorentzVector oleps, leps, tmp;
  for(UInt_t i = 0; i < T_Muon_Energy->size(); ++i){
    tmp.SetPxPyPzE(MuPx.at(i), MuPy.at(i), T_Muon_Pz->at(i), T_Muon_Energy->at(i));
    oleps += tmp;
    if(flag == 1) { MuPx.at(i) += scale*MuPx.at(i); MuPy.at(i) += scale*MuPy.at(i); }
    if(flag == 2) { MuPx.at(i) -= scale*MuPx.at(i); MuPy.at(i) -= scale*MuPy.at(i); }
    tmp.SetPxPyPzE(MuPx.at(i), MuPy.at(i), T_Muon_Pz->at(i), T_Muon_Energy->at(i));
    leps += tmp;
  }
  for(UInt_t i = 0; i < T_Elec_Energy->size(); ++i){
    tmp.SetPxPyPzE(ElPx.at(i), ElPy.at(i), T_Elec_Pz->at(i), T_Elec_Energy->at(i));
    oleps += tmp;
    if(flag == 1) { ElPx.at(i) += scale*ElPx.at(i); ElPy.at(i) += scale*ElPy.at(i); }
    if(flag == 2) { ElPx.at(i) -= scale*ElPx.at(i); ElPy.at(i) -= scale*ElPy.at(i); }
    tmp.SetPxPyPzE(ElPx.at(i), ElPy.at(i), T_Elec_Pz->at(i), T_Elec_Energy->at(i));
    leps += tmp;
  }
  propagateMET(leps, oleps);
  return;
}
void TreeAnalysisTop::ScaleMET(int flag){
  // first try on MET uncertainty
  if(gIsData) return; // don't scale data
  
  TLorentzVector umet, jets, leps, tmp;
  jets.SetPtEtaPhiM(0., 0., 0., 0.); // init
  leps.SetPtEtaPhiM(0., 0., 0., 0.); // init
  tmp.SetPtEtaPhiM(0., 0., 0., 0.);  // init
  umet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.); // add met
  
  // subtract uncleaned jets
  for (UInt_t i=0; i<T_JetAKCHS_Energy->size(); i++) {
    if (!IsGoodJet(i, 15.)) continue; // do this on all jets in the event, not only the good jets with pT > 40
    tmp.SetPxPyPzE(T_JetAKCHS_Px->at(i),T_JetAKCHS_Py->at(i),T_JetAKCHS_Pz->at(i),T_JetAKCHS_Energy->at(i));
    umet += tmp;
    jets += tmp;
    tmp.SetPtEtaPhiE(0., 0., 0., 0.);
  }
  // subtract muons
  for (UInt_t i=0; i<T_Muon_Energy->size(); i++) {
    if (!IsTightMuon(i)) continue;
    tmp.SetPxPyPzE(MuPx.at(i),MuPy.at(i),T_Muon_Pz->at(i),T_Muon_Energy->at(i));
    umet += tmp;
    leps += tmp;
    tmp.SetPtEtaPhiM(0., 0., 0., 0.);
  }
  // subtract electrons
  for (UInt_t i=0; i<T_Elec_Energy->size(); i++) {
    if (!IsTightElectron(i)) continue;
    tmp.SetPxPyPzE(ElPx.at(i),ElPy.at(i),T_Elec_Pz->at(i),T_Elec_Energy->at(i));
    umet += tmp;
    leps += tmp;
    tmp.SetPtEtaPhiM(0., 0., 0., 0.);
  }
  // scale the unclustered energy by 10%
  if (flag == 0) tmp.SetPtEtaPhiE(1.1 * umet.Pt(), umet.Eta(), umet.Phi(), umet.E());
  if (flag == 1) tmp.SetPtEtaPhiE(0.9 * umet.Pt(), umet.Eta(), umet.Phi(), umet.E());
  
  // subtract the leptons and jets again
  tmp -= leps;
  tmp -= jets;
  
  setMET(tmp.Pt());
  return;
}
//------------------------------------------------------------------------------
// GetGenINFO
//------------------------------------------------------------------------------
/*int TreeAnalysisTop::GetFakesSource(){
  // First, make sure we do have two leptons...
  if (fHypLepton1.index == -1) return;  
  if (fHypLepton2.index == -1) return;

#ifdef __ISMC
  // First check for two leptons of Status 3
  TLorentVector gen;
  if (T_Gen_MuonSt3_PID->size() == 2) {
    if (T_Gen_MuonSt3_PID->at(0)*T_Gen_MuonSt3_PID->at(1) < 0) return WrongSign;//OS
    if (T_Gen_MuonSt3_PID->at(0)*T_Gen_MuonSt3_PID->at(1) > 0) return RightSign;//SS
  }
  if (T_Gen_ElecSt3_PID->size() == 2) {
    if (T_Gen_ElecSt3_PID->at(0)*T_Gen_ElecSt3_PID->at(1) < 0) return WrongSign;//OS
    if (T_Gen_ElecSt3_PID->at(0)*T_Gen_ElecSt3_PID->at(1) > 0) return RightSign;//SS
  }
  
  
  // NOW 
  for (UInt_t mu=0; mu<T_Gen_MuonSt3_PID->size(); mu++){
    gen.SetPxPyPzE(T_Gen_MuonSt3_Px->at(mu),  T_Gen_MuonSt3_Py->at(mu),		   
		   T_Gen_MuonSt3_Pz->at(mu),  T_Gen_MuonSt3_Energy->at(mu));
    if(fHypLepton1.type == 0 && mu.DeltaR(fHypLepton1.p) < 0.1) { // gen muon matches hyp. lepton
      if (T_Gen_MuonSt3_PID->at(mu))
    }
    if(fHypLepton2.type == 0 && mu.DeltaR(fHypLepton2.p) < 0.1) { // gen muon matches hyp. lepton
      
    }
  }
  for (UInt_t el=0; el<T_Gen_ElecSt3_PID->size(); el++){
  }

  Bool_t IsMu1FromW = false;
  Bool_t IsMu2FromW = false;
  Bool_t IsMuFromHF = false; 
 TLorentVector genmu;
  for (UInt_t mu=0; mu<T_Gen_Muon_PID->size(); mu++){
    genmu.SetPxPyPzE(T_Gen_Muon_Px->at(mu),T_Gen_Muon_Py->at(mu),T_Gen_Muon_Pz->at(mu),T_Gen_Muon_Energy->at(mu));
    if(fHypLepton1.type == 0 && genmu.DeltaR(fHypLepton1.p) < 0.1) { // gen muon matches hyp. lepton
      if (T_Gen_Muon_MSt->at(mu) == 3 && TMath::Abs(T_Gen_Muon_MPID->at(mu)) == 24) IsMu1FromW  = true; 
      if (T_Gen_Muon_MSt->at(mu) == 3 && TMath::Abs(T_Gen_Muon_MPID->at(mu)) == 5 ) ISMu1FromHF = true;
    }
    if(fHypLepton2.type == 0 && genmu.DeltaR(fHypLepton2.p) < 0.1) { // gen muon matches hyp. lepton
      if (T_Gen_Muon_MSt->at(mu) == 3 && TMath::Abs(T_Gen_Muon_MPID->at(mu)) == 24) IsMu2FromW  = true; 
      if (T_Gen_Muon_MSt->at(mu) == 3 && TMath::Abs(T_Gen_Muon_MPID->at(mu)) == 5 ) ISMu2FromHF = true;
    }
  }
  
  if (IsMu1FromW && IsMu1FromW) return; //dileptonic event in generation
  for (UInt_t el=0; el<T_Gen_Elec_PID->size(); el++){
    if(fHypLepton1.type == 1) ;
    if(fHypLepton2.type == 1) ;
  }
  
#endif 
  return 0; // return 0 for data...
  }*/
//SANTIvoid TreeAnalysisTop::GetGenMuon(){
//SANTI  
//SANTI  UInt_t muonSize = 0;
//SANTI  
//SANTI  muonSize = T_Gen_MuonSt3_Energy->size();
//SANTI  
//SANTI  for (UInt_t i=0; i<muonSize; i++) {
//SANTI    
//SANTI    TLorentzVector Muon(T_Gen_MuonSt3_Px    ->at(i),
//SANTI			T_Gen_MuonSt3_Py    ->at(i),
//SANTI			T_Gen_MuonSt3_Pz    ->at(i),
//SANTI			T_Gen_MuonSt3_Energy->at(i));  
//SANTI    
//SANTI    if (Muon.Pt()>20 &&
//SANTI	TMath::Abs(Muon.Eta()) < 2.1) { 
//SANTI      
//SANTI      if(T_Gen_MuonSt3_PID->at(i)== 13) Gen_Muon_Charge.push_back(-1.);
//SANTI      if(T_Gen_MuonSt3_PID->at(i)==-13) Gen_Muon_Charge.push_back( 1.);
//SANTI      
//SANTI      Gen_Muon.push_back(Muon);
//SANTI      
//SANTI    }// if(GenMuon)
//SANTI  }// for(muonSize)
//SANTI
//SANTI  for (UInt_t i=0; i<T_Gen_TauSt3_PID->size(); i++) {
//SANTI    if (T_Gen_TauSt3_IsLepDec->at(i)) {
//SANTI      if (abs(T_Gen_TauSt3_LepDec_PID->at(i)) == 13){
//SANTI
//SANTI    
//SANTI	TLorentzVector MuonTau(T_Gen_TauSt3_LepDec_Px    ->at(i),
//SANTI			       T_Gen_TauSt3_LepDec_Py    ->at(i),
//SANTI			       T_Gen_TauSt3_LepDec_Pz    ->at(i),
//SANTI			       T_Gen_TauSt3_LepDec_Energy->at(i));  
//SANTI    
//SANTI    if (MuonTau.Pt()>20 &&
//SANTI	TMath::Abs(MuonTau.Eta()) < 2.1) {
//SANTI      
//SANTI      if(T_Gen_TauSt3_LepDec_PID->at(i)== 13) Gen_Muon_Charge.push_back(-1.);
//SANTI      if(T_Gen_TauSt3_LepDec_PID->at(i)==-13) Gen_Muon_Charge.push_back( 1.);
//SANTI      
//SANTI      Gen_Muon.push_back(MuonTau);
//SANTI      
//SANTI    }// if(GenMuonTau)
//SANTI    
//SANTI    
//SANTI      }// if(muon)
//SANTI    }// if(LepDecay)
//SANTI  }// for(Taus)
//SANTI  
//SANTI  nSGenMuon = Gen_Muon.size();
//SANTI}
//SANTI
//SANTI
//SANTI//------------------------------------------------------------------------------
//SANTI// GetGenElec
//SANTI//------------------------------------------------------------------------------
//SANTIvoid TreeAnalysisTop::GetGenElec()
//SANTI{
//SANTI  UInt_t elecSize = 0;
//SANTI
//SANTI  elecSize = T_Gen_ElecSt3_Energy->size();
//SANTI
//SANTI  for (UInt_t i=0; i<elecSize; i++) {
//SANTI
//SANTI    TLorentzVector Elec(T_Gen_ElecSt3_Px    ->at(i),
//SANTI			T_Gen_ElecSt3_Py    ->at(i),
//SANTI			T_Gen_ElecSt3_Pz    ->at(i),
//SANTI			T_Gen_ElecSt3_Energy->at(i));
//SANTI    
//SANTI    if (Elec.Pt()>20 &&
//SANTI	TMath::Abs(Elec.Eta()) < 2.5) { 
//SANTI      
//SANTI      if(T_Gen_ElecSt3_PID->at(i)== 11) Gen_Elec_Charge.push_back(-1.);
//SANTI      if(T_Gen_ElecSt3_PID->at(i)==-11) Gen_Elec_Charge.push_back( 1.);
//SANTI      
//SANTI      Gen_Elec.push_back(Elec);
//SANTI    }
//SANTI  }
//SANTI  
//SANTI  for (UInt_t i=0; i<T_Gen_TauSt3_PID->size(); i++) {
//SANTI    if (T_Gen_TauSt3_IsLepDec->at(i)) {
//SANTI      if (abs(T_Gen_TauSt3_LepDec_PID->at(i)) == 11){
//SANTI
//SANTI    
//SANTI	TLorentzVector ElecTau(T_Gen_TauSt3_LepDec_Px    ->at(i),
//SANTI			       T_Gen_TauSt3_LepDec_Py    ->at(i),
//SANTI			       T_Gen_TauSt3_LepDec_Pz    ->at(i),
//SANTI			       T_Gen_TauSt3_LepDec_Energy->at(i));  
//SANTI    
//SANTI    if (ElecTau.Pt()>20 &&
//SANTI	TMath::Abs(ElecTau.Eta()) < 2.5) {
//SANTI      
//SANTI      if(T_Gen_TauSt3_LepDec_PID->at(i)== 11) Gen_Elec_Charge.push_back(-1.);
//SANTI      if(T_Gen_TauSt3_LepDec_PID->at(i)==-11) Gen_Elec_Charge.push_back( 1.);
//SANTI      
//SANTI      Gen_Elec.push_back(ElecTau);
//SANTI      
//SANTI    }// if(GenMuonTau)
//SANTI    
//SANTI    
//SANTI      }// if(muon)
//SANTI    }// if(LepDecay)
//SANTI  }// for(Taus)
//SANTI
//SANTI  nSGenElec = Gen_Elec.size();
//SANTI}



