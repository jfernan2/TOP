#include "TreeAnalysisTop.h"

#include <fstream>
#include <iostream>
#include <math.h>

#if !defined(__CINT__)
ClassImp(TreeAnalysisTop);
#endif

const float gJetEtCut = 30.;
//#define DEBUG
//#define __ISPDF

TreeAnalysisTop::TreeAnalysisTop(TTree* tree) : PAFAnalysis(tree)
{
  fHDummy        = 0;
  hWeight        = 0;
  fHTopPtWeight  = 0;
  fHpdfWeightSum = 0;
  fHpdfWeight    = 0;
}


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
  InitialiseKinematicHistos();

#ifdef __ISMC
  InitialiseGenHistos();
#endif  
  fHTopPtWeight  = CreateH1F("H_TopPtWeight" , "TopPt Weight"   , 100,    0,    2);
  fHpdfWeightSum = CreateH1F("H_pdfWeightSum", "PDF sum Weights",  52, -0.5, 51.5);
  fHpdfWeight    = CreateH1F("H_pdfWeight"   , "PDF Weights"    ,  52, -0.5, 51.5);

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
  hWeight = CreateH1F("hWeight","",200,0,1);
  //++ Yields histograms
  if (gDoSF) {
    fHyields[Muon][Norm]   = CreateH1F("H_Yields_"+gChanLabel[Muon],"", iNCUTS, -0.5, iNCUTS-0.5); 
    fHyields[Elec][Norm]   = CreateH1F("H_Yields_"+gChanLabel[Elec],"", iNCUTS, -0.5, iNCUTS-0.5);
    fHSSyields[Muon][Norm] = CreateH1F("H_SSYields_"+gChanLabel[Muon],"", iNCUTS, -0.5, iNCUTS-0.5); 
    fHSSyields[Elec][Norm] = CreateH1F("H_SSYields_"+gChanLabel[Elec],"", iNCUTS, -0.5, iNCUTS-0.5);
  }
  if (gDoDF) {
    fHyields[ElMu][Norm]   = CreateH1F("H_Yields_"+gChanLabel[ElMu],"", iNCUTS, -0.5, iNCUTS-0.5);
    fHSSyields[ElMu][Norm] = CreateH1F("H_SSYields_"+gChanLabel[ElMu],"", iNCUTS, -0.5, iNCUTS-0.5);
  }
  
  if (gDoSystStudies){
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      if (!gDoSF && chan==Muon) continue;
      if (!gDoSF && chan==Elec) continue;
      if (!gDoDF && chan==ElMu) continue;
      for (size_t sys=1; sys<gNSYST; sys++){
	fHyields[chan][sys]   = CreateH1F("H_Yields_"+gChanLabel[chan]+"_"+SystName[sys],"",iNCUTS,-0.5,iNCUTS-0.5);
	fHSSyields[chan][sys] = CreateH1F("H_SSYields_"+gChanLabel[chan]+"_"+SystName[sys],"", iNCUTS, -0.5, iNCUTS-0.5);
      }
    }
  }
  
  for (size_t chan=0; chan<gNCHANNELS; chan++){
    if (!gDoSF && chan==Muon) continue;
    if (!gDoSF && chan==Elec) continue;
    if (!gDoDF && chan==ElMu) continue;
    for (size_t cut=0; cut<iNCUTS; cut++){
      fHLepSys [chan][cut] = CreateH1F("H_LepSys_" +gChanLabel[chan]+"_"+sCut[cut],"LepSys" , 400, 0, 0.04);
      fHTrigSys[chan][cut] = CreateH1F("H_TrigSys_"+gChanLabel[chan]+"_"+sCut[cut],"TrigSys", 400, 0, 0.04);
    }
  }
}
void TreeAnalysisTop::InitialiseKinematicHistos(){
  //++ Kinematic histograms
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    if (!gDoSF && ch==Muon) continue;
    if (!gDoSF && ch==Elec) continue;
    if (!gDoDF && ch==ElMu) continue;

    for (size_t cut=0; cut<iNCUTS; cut++){
      fHMET[ch][cut]         = CreateH1F("H_MET_"        +gChanLabel[ch]+"_"+sCut[cut],"MET"       ,  5000,0,500);
      fHDiLepPt[ch][cut]     = CreateH1F("H_DiLepPt_"    +gChanLabel[ch]+"_"+sCut[cut],"DiLepPt"   , 1800,20,200); 
      fHLep0Pt[ch][cut]      = CreateH1F("H_Lep0Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Lep0Pt"    , 1800,20,200);
      fHLep1Pt[ch][cut]      = CreateH1F("H_Lep1Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Lep1Pt"    , 1800,20,200);
      fHDelLepPhi[ch][cut]   = CreateH1F("H_DelLepPhi_"  +gChanLabel[ch]+"_"+sCut[cut],"DelLepPhi" , 100,-3.2, 3.2);
      fHNJets[ch][cut]       = CreateH1F("H_NJets_"      +gChanLabel[ch]+"_"+sCut[cut],"NJets"     , 10 ,-0.5, 9.5);
      fHHT[ch][cut]          = CreateH1F("H_HT_"         +gChanLabel[ch]+"_"+sCut[cut],"HT"        , 4700,30,500);
      fHNBtagJets[ch][cut]   = CreateH1F("H_NBtagJets_"  +gChanLabel[ch]+"_"+sCut[cut],"NBtagJets" , 4 ,-0.5, 3.5);
      fHJet0Pt[ch][cut]      = CreateH1F("H_Jet0Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Jet0Pt"    , 2700,30,300);
      fHJet1Pt[ch][cut]      = CreateH1F("H_Jet1Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Jet1Pt"    , 2700,30,300);
      fHBtagJet0Pt[ch][cut]  = CreateH1F("H_BtagJet0Pt_" +gChanLabel[ch]+"_"+sCut[cut],"BtagJet0Pt", 2700,30,300);
      
      fHInvMass[ch][cut][0]       = CreateH1F("H_InvMass_"    +gChanLabel[ch]+"_"+sCut[cut],"InvMass"   , 300, 0., 300.);
      fHSSInvMass[ch][cut][0]     = CreateH1F("HSS_InvMass_"  +gChanLabel[ch]+"_"+sCut[cut],"InvMass"   , 300, 0., 300.);
      
      fHNBtagsNJets[ch][cut][0]   = CreateH1F("H_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut]  ,"NBtagsNJets"   ,15 , -0.5, 14.5);
      fHSSNBtagsNJets[ch][cut][0] = CreateH1F("HSS_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut],"SS_NBtagsNJets",15 , -0.5, 14.5);
      
      fHAbsDelPhiLeps[ch][cut][0]   = CreateH1F("H_AbsDelPhiLeps_"+gChanLabel[ch]+"_"+sCut[cut]  ,"AbsDelPhiLeps"   , 28,-0.2, 1.2);
      fHSSAbsDelPhiLeps[ch][cut][0] = CreateH1F("HSS_AbsDelPhiLeps_"+gChanLabel[ch]+"_"+sCut[cut],"SS_AbsDelPhiLeps", 28,-0.2, 1.2);
      
      fHdelPhi2LeadJets[ch][cut][0]   = CreateH1F("H_delPhi2LeadJets_"+gChanLabel[ch]+"_"+sCut[cut]  ,"delPhi2LeadJets"   , 28,-0.2, 1.2);
      fHSSdelPhi2LeadJets[ch][cut][0] = CreateH1F("HSS_delPhi2LeadJets_"+gChanLabel[ch]+"_"+sCut[cut],"SS_delPhi2LeadJets", 28,-0.2, 1.2);
      
      fHminDelRJetsLeps[ch][cut][0]   = CreateH1F("H_minDelRJetsLeps_"+gChanLabel[ch]+"_"+sCut[cut]  ,"minDelRJetsLeps"   ,   500,0.0, 5.0);
      fHSSminDelRJetsLeps[ch][cut][0] = CreateH1F("HSS_minDelRJetsLeps_"+gChanLabel[ch]+"_"+sCut[cut],"SS_minDelRJetsLeps",   500,0.0, 5.0);

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
      //fHAbsDelPhiLep[ch][cut]= CreateH1F("H_AbsDelPhiLep_" +gChanLabel[ch]+"_"+sCut[cut],"AbsDelPhiLep" , 28,-0.2, 1.2);
      fHvertices[ch][cut] = CreateH1F("H_Vtx_"+gChanLabel[ch]+"_"+sCut[cut],"", 71, -0.5, 70.5); 
	    

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
    if (!gDoSF && ch==Muon) continue;
    if (!gDoSF && ch==Elec) continue;
    if (!gDoDF && ch==ElMu) continue;

    for (size_t cut=0; cut<iNCUTS; cut++){
      for (size_t sys=1; sys<gNSYST; sys++){
	histoname = "H_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHNBtagsNJets[ch][cut][sys] = CreateH1F(histoname,"NBtagsNJets", 15 , -0.5, 14.5);
	
	histoname = "H_InvMass_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHInvMass[ch][cut][sys]     = CreateH1F(histoname, "InvMass"   , 300, 0., 300.);
	
	histoname = "H_AbsDelPhiLeps_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHAbsDelPhiLeps[ch][cut][sys] = CreateH1F(histoname,"AbsDelPhiLeps", 28,-0.2, 1.2);
	
	histoname = "H_delPhi2LeadJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHdelPhi2LeadJets[ch][cut][sys] = CreateH1F(histoname,"delPhi2LeadJets", 28,-0.2, 1.2);
	
	histoname = "H_minDelRJetsLeps_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHminDelRJetsLeps[ch][cut][sys] = CreateH1F(histoname,"minDelRJetsLeps", 500, 0., 5.0);

	histoname = "HSS_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHSSNBtagsNJets[ch][cut][sys] = CreateH1F(histoname,"SS_NBtagsNJets", 15 , -0.5, 14.5);
	
	histoname = "HSS_InvMass_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHSSInvMass[ch][cut][sys]     = CreateH1F(histoname,"SS_InvMass"   , 300, 0., 300.);
	
	histoname = "HSS_AbsDelPhiLeps_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHSSAbsDelPhiLeps[ch][cut][sys] = CreateH1F(histoname,"SS_AbsDelPhiLeps", 28,-0.2, 1.2);
	
	histoname = "HSS_delPhi2LeadJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHSSdelPhi2LeadJets[ch][cut][sys] = CreateH1F(histoname,"SS_delPhi2LeadJets", 28,-0.2, 1.2);
	
	histoname = "HSS_minDelRJetsLeps_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
	fHSSminDelRJetsLeps[ch][cut][sys] = CreateH1F(histoname,"SS_minDelRJetsLeps", 500, 0., 5.0);
	
      }
    }
  }
}
//------------------------------------------------------------------------------
// InsideLoop
//------------------------------------------------------------------------------
void TreeAnalysisTop::SetOriginalObjects(){
  // To be called once per event, saving information in tmp vectors
  // for systematic studies.
  ResetHypLeptons();
  
  gSysSource = Norm;

  // SAVING ORIGINAL VALUES FOR MET, JET, LEPTONS for SYST
  JetEt.clear();
  JetPhi.clear();
  MuPx.clear();
  MuPy.clear();
  ElPx.clear();
  ElPy.clear();
  MET  = 0.;
  MET_Phi = 0.;
  
  // Save original values for MET, Jets and Leptons
  TLorentzVector j;
  for (UInt_t i=0; i<T_JetAKCHS_Et->size(); i++){    
    j.SetPxPyPzE(T_JetAKCHS_Px->at(i),T_JetAKCHS_Py->at(i),T_JetAKCHS_Pz->at(i),T_JetAKCHS_Energy->at(i));
    JetEt.push_back(j.Et());
    JetPhi.push_back(j.Phi());
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
void TreeAnalysisTop::SetEventObjects(){
  //ResetAnalysisTree();
  ResetHypLeptons();
  
  fChargeSwitch = false;
  
  // EVENT WEIGHTS
  PUSF = 1.;
  EventWeight = 1.;

  // USEFUL COUNTERS
  nGenLepton = 0;
  nGenElec   = 0;
  nGenMuon   = 0;
  nGenTau    = 0;
  nTauElec   = 0;
  nTauMuon   = 0;
  
  nGoodVertex = 0;
  nBtags      = 0;
  nJets       = 0;
  nMuon       = 0;
  nElec       = 0;
  nLeptons    = 0;
  
   
  //// READ AND SAVE OBJETS...
  Jet.clear();
  Lepton.clear();
  
  nLeptons = getSelectedLeptons();
  nJets    = getSelectedJets();
  nBtags   = getNBTags();
}
void TreeAnalysisTop::ResetOriginalObjects(){
  
  // Save original values for MET, Jets and Leptons
  TLorentzVector j;
  for (UInt_t i=0; i<T_JetAKCHS_Et->size(); i++){    
    j.SetPxPyPzE(T_JetAKCHS_Px->at(i),T_JetAKCHS_Py->at(i),T_JetAKCHS_Pz->at(i),T_JetAKCHS_Energy->at(i));
    JetEt[i]  = j.Et();
    JetPhi[i] = j.Phi();
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
      (abs(T_Gen_StopMass->at(0)-T_Gen_StopMass->at(1)) >  0.1 || 
       abs(T_Gen_Chi0Mass->at(0)-T_Gen_Chi0Mass->at(1)) >  0.1 || 
       abs(T_Gen_StopMass->at(0)-gStopMass)             > 6.25 || 
       abs(T_Gen_Chi0Mass->at(0)-gLspMass)              > 6.25 )
      ) return;
#endif
  
  // Init data members
  //----------------------------------------------------------------------------
  SetOriginalObjects();
  SetEventObjects();

#ifdef __ISPDF
  for (int pdf=0; pdf<52; pdf++)
    fHpdfWeightSum->Fill(pdf,T_Event_pdfWeight->at(pdf));
#endif

  // Get number of generated leptons 
  //----------------------------------------------------------------------------
#ifdef __ISMC
  SelectedGenLepton();
  
  if (gSampleName == "TTJets_MadSpin"        ){
    
    Float_t Weight = 1.; 
    TLorentzVector top;
    for (size_t t=0; t<T_Gen_tSt3_Px->size(); t++){
      top.SetPxPyPzE(T_Gen_tSt3_Px->at(t),T_Gen_tSt3_Py->at(t),T_Gen_tSt3_Pz->at(t),T_Gen_tSt3_Energy->at(t));
      Float_t pt    = TMath::Min(top.Pt(), 400.);
      //      Float_t topSF = TMath::Exp(0.148 - 0.00129 * pt);
      Float_t topSF = TMath::Exp(0.156 - 0.00137 * pt);
      Weight *= topSF;
    }
    fHTopPtWeight->Fill(TMath::Sqrt(Weight));
  }

  if (gSampleName == "TTJets_MadSpin"         && nGenLepton != 2) return;
  if (gSampleName == "TTJets_matchingup"      && nGenLepton != 2) return;
  if (gSampleName == "TTJets_matchingdown"    && nGenLepton != 2) return;
  if (gSampleName == "TTJets_scaleup"         && nGenLepton != 2) return;
  if (gSampleName == "TTJets_scaledown"       && nGenLepton != 2) return;
  if (gSampleName == "TTJetsFullLeptMGtauola" && nGenLepton != 2) return;

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

  // FOR PDF Uncertainties
  
  // Accept only events with a good vertex
  //----------------------------------------------------------------------------
  if (SelectedVertexIndex() < 0) return;
  
  // Fill Yields...
  //----------------------------------------------------------------------------
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
  
  ////////////////////////////////////////
  // BTAGGING SYSTEMATICS
  ////////////////////////////////////////
  ResetOriginalObjects();
  gSysSource = BtagUp;
  SetEventObjects();
  FillYields(BtagUp);
  fChargeSwitch = true;  
  FillYields(BtagUp); /* Get SS yields....*/  
  fChargeSwitch = false;
  
  ResetOriginalObjects();
  gSysSource = BtagDown;
  SetEventObjects();  
  FillYields(BtagDown);
  fChargeSwitch = true;
  FillYields(BtagDown); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = MisTagUp;
  SetEventObjects();
  FillYields(MisTagUp);
  fChargeSwitch = true;
  FillYields(MisTagUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = MisTagDown;
  SetEventObjects();
  FillYields(MisTagDown);
  fChargeSwitch = true;
  FillYields(MisTagDown); /// Get SS yields....
  fChargeSwitch = false;
  
  ////////////////////////////////////////
  // JES/JER SYSTEMATICS
  ////////////////////////////////////////
  ResetOriginalObjects();
  SmearJetPts(1);
  gSysSource = JESUp;
  SetEventObjects();
  FillYields(JESUp);
  fChargeSwitch = true;
  FillYields(JESUp); /// Get SS yields....
  fChargeSwitch = false;
  
  ResetOriginalObjects();
  SmearJetPts(2);
  gSysSource = JESDown;
  SetEventObjects();
  FillYields(JESDown);
  fChargeSwitch = true;
  FillYields(JESDown); /// Get SS yields....
  fChargeSwitch = false;
 
  ResetOriginalObjects();
  SmearJetPts(3);
  gSysSource = JER;
  SetEventObjects();
  FillYields(JER);
  fChargeSwitch = true;
  FillYields(JER); /// Get SS yields....
  fChargeSwitch = false;

  /////////////////////////////////////////////////////////////
  // LEPTON SCALE SYSTEMATICS
  ////////////////////////////////////////////////////
  ResetOriginalObjects();
  ScaleLeptons(1); //up
  gSysSource = LESUp;
  SetEventObjects();
  FillYields(LESUp);
  fChargeSwitch = true;
  FillYields(LESUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  ScaleLeptons(2); //down
  gSysSource = LESDown;
  SetEventObjects();
  FillYields(LESDown);
  fChargeSwitch = true;
  FillYields(LESDown); /// Get SS yields....
  fChargeSwitch = false;

  /*  ResetOriginalObjects();
  gSysSource = LepUp;
  SetEventObjects();
  FillYields(LepUp);
  fChargeSwitch = true;
  FillYields(LepUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = LepDown;
  SetEventObjects();
  FillYields(LepDown);
  fChargeSwitch = true;
  FillYields(LepDown); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = TrigUp;
  SetEventObjects();
  FillYields(TrigUp);
  fChargeSwitch = true;
  FillYields(TrigUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = TrigDown;
  SetEventObjects();
  FillYields(TrigDown);
  fChargeSwitch = true;
  FillYields(TrigDown); /// Get SS yields....
  fChargeSwitch = false;
  */
  
  // PILE UP UNCERTAINTY
  ResetOriginalObjects();
#ifdef __ISMC
  PUSF = fPUWeightUp->GetWeight((int)T_Event_nTruePU);
#endif
  gSysSource = PUUp;
  SetEventObjects();
  FillYields(PUUp);
  
  ResetOriginalObjects();
#ifdef __ISMC
  PUSF = fPUWeightDown->GetWeight((int)T_Event_nTruePU);
#endif
  gSysSource = PUDown;
  SetEventObjects();
  FillYields(PUDown); 
  
  // TOP PT
  ResetOriginalObjects();
  gSysSource = TopPt;
  SetEventObjects();
  FillYields(TopPt);
  fChargeSwitch = true;
  FillYields(TopPt); /// Get SS yields....
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
  
  GetInputParameters()->TheNamedBool ("IsData"       , gIsData);
  GetInputParameters()->TheNamedFloat("weight"       , gWeight); // cross section / events in the sample
  GetInputParameters()->TheNamedFloat("LumiForPU"    , gLumiForPU);
  GetInputParameters()->TheNamedFloat("TotalLumi"    , gTotalLumi);
  GetInputParameters()->TheNamedBool ("DoSystStudies", gDoSystStudies);
  GetInputParameters()->TheNamedBool ("UseCSVM"      , gUseCSVM);
  GetInputParameters()->TheNamedBool ("DoSF"         , gDoSF);
  GetInputParameters()->TheNamedBool ("DoDF"         , gDoDF);
  GetInputParameters()->TheNamedFloat("stopMass"     , gStopMass);
  GetInputParameters()->TheNamedFloat("lspMass"      , gLspMass );

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
//////////////////////////////////////////////////////////////////
/// Get METHODS
/////////////////////////////////////////////////////////////////
float TreeAnalysisTop::getHT(){
  float ht(0);
  for (unsigned int i=0; i<Jet.size(); i++) ht+=Jet[i].p.Pt();
  return ht;
}
float TreeAnalysisTop::getJetPtIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  return Jet[ind].p.Pt();
}
float TreeAnalysisTop::getBtagJetPtIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  Int_t btagInd = 0;
  if (ind==0) btagInd = getLeadingJetbTag();
  else  return -999.;
  return Jet[btagInd].p.Pt();
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
  if (chan==Muon){
    err1 = fLeptonSF->GetTightMuonSF_err(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
    err2 = fLeptonSF->GetTightMuonSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
  }
  if (chan==ElMu){
    err1 = fLeptonSF->GetTightMuonSF_err    (fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
    err2 = fLeptonSF->GetTightElectronSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
  }
  if (chan==Elec){
    err1 = fLeptonSF->GetTightElectronSF_err(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
    err2 = fLeptonSF->GetTightElectronSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
  }
  return TMath::Sqrt(err1*err1+err2*err2);
}
float TreeAnalysisTop::getTriggerError(gChannel chan){
  float trig(0.);

  if (chan==Muon) trig = fLeptonSF->GetDoubleMuSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
  if (chan==ElMu) trig = fLeptonSF->GetMuEGSF_err    (fHypLepton2.p.Eta(),fHypLepton1.p.Eta());
  if (chan==Elec) trig = fLeptonSF->GetDoubleElSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
  return trig;
}
float TreeAnalysisTop::getSF(gChannel chan) {
  if (gIsData)              return 1.; //Don't scale data
  
  float id1(1.),id2(1.), trig(1.);
  if (chan == Muon){
    id1  = fLeptonSF->GetTightMuonSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
    id2  = fLeptonSF->GetTightMuonSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
    trig = fLeptonSF->GetDoubleMuSF (fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
    //    if (gSysSource == LepDown)
  } 
  else if (chan == Elec){
    id1  = fLeptonSF->GetTightElectronSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta()); 
    id2  = fLeptonSF->GetTightElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta()); 
    trig = fLeptonSF->GetDoubleElSF     (fHypLepton1.p.Eta(),fHypLepton2.p.Eta()); 
  }
  else if (chan == ElMu){
    id1  = fLeptonSF->GetTightMuonSF    (fHypLepton1.p.Pt(), fHypLepton1.p.Eta()); 
    id2  = fLeptonSF->GetTightElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
    trig = fLeptonSF->GetMuEGSF         (fHypLepton2.p.Eta(),fHypLepton1.p.Eta());
  }
  return (PUSF*id1*id2*trig);
  
}


// Return SF of the pt of the top. Only apply SF if the process is ttbar
float TreeAnalysisTop::getTopPtSF()
{
  Float_t Weight = 1.; 

  if (!gSampleName.Contains("TTJets")) return Weight;
  
  if (gSysSource == TopPt) {
#ifdef __ISMC
    if (T_Gen_tSt3_Px->size() != 2) return Weight;
    
    TLorentzVector top;
    Float_t topSF = 0.;

    for (size_t t=0; t<T_Gen_tSt3_Px->size(); t++) {
      top.SetPxPyPzE(T_Gen_tSt3_Px->at(t),T_Gen_tSt3_Py->at(t),T_Gen_tSt3_Pz->at(t),T_Gen_tSt3_Energy->at(t));
      Float_t pt = TMath::Min(top.Pt(), 400.);
      //      topSF = TMath::Exp(0.148 - 0.00129 * pt);
      topSF = TMath::Exp(0.156 - 0.00137 * pt);
      Weight *= topSF;
    }
    Weight = TMath::Sqrt(Weight);
#endif
  }
  
  return Weight;
}


void TreeAnalysisTop::FillDYHistograms(){

  float Mll = 0.;
  if (PassTriggerEMu()  && IsElMuEvent()){
    // Define Hypothesis Leptons...
    EventWeight = gWeight * getSF(ElMu);
    Mll = (fHypLepton1.p+fHypLepton2.p).M();
    
    if (PassesMllVeto() && PassesMuonEta2p1(ElMu) && Passes3rdLeptonVeto()){
      fHDY_InvMassVsNPV   [ElMu][iDilepton]->Fill(nGoodVertex, Mll, EventWeight);
      fHDY_InvMassVsMET   [ElMu][iDilepton]->Fill(getMET()   , Mll, EventWeight);
      fHDY_InvMassVsNjets [ElMu][iDilepton]->Fill(getNJets() , Mll, EventWeight);
      fHDY_InvMassVsNbtags[ElMu][iDilepton]->Fill(getNBTags(), Mll, EventWeight);
      fHDY_InvMass        [ElMu][iDilepton]->Fill(             Mll, EventWeight);
	
      if (PassesMETCut())   {
        fHDY_InvMassVsNPV   [ElMu][iMET]->Fill(nGoodVertex, Mll, EventWeight);
        fHDY_InvMassVsMET   [ElMu][iMET]->Fill(getMET()   , Mll, EventWeight);
        fHDY_InvMassVsNjets [ElMu][iMET]->Fill(getNJets() , Mll, EventWeight);
        fHDY_InvMassVsNbtags[ElMu][iMET]->Fill(getNBTags(), Mll, EventWeight);
        fHDY_InvMass	    [ElMu][iMET]->Fill( 	    Mll, EventWeight);
      
        if (PassesNJetsCut()) {
	  fHDY_InvMassVsNPV   [ElMu][i2jets]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [ElMu][i2jets]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [ElMu][i2jets]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[ElMu][i2jets]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [ElMu][i2jets]->Fill(             Mll, EventWeight);
	  
	  if (PassesNBtagCut()) {
	    fHDY_InvMassVsNPV   [ElMu][i1btag]->Fill(nGoodVertex, Mll, EventWeight);
	    fHDY_InvMassVsMET   [ElMu][i1btag]->Fill(getMET()   , Mll, EventWeight);
	    fHDY_InvMassVsNjets [ElMu][i1btag]->Fill(getNJets() , Mll, EventWeight);
	    fHDY_InvMassVsNbtags[ElMu][i1btag]->Fill(getNBTags(), Mll, EventWeight);
	    fHDY_InvMass        [ElMu][i1btag]->Fill(             Mll, EventWeight);
	  }
	}
      }
      if (getNBTags() == 1){
	fHDY_InvMassVsNPV   [ElMu][iExact1btag]->Fill(nGoodVertex, Mll, EventWeight);
	fHDY_InvMassVsMET   [ElMu][iExact1btag]->Fill(getMET()   , Mll, EventWeight);
	fHDY_InvMassVsNjets [ElMu][iExact1btag]->Fill(getNJets() , Mll, EventWeight);
	fHDY_InvMassVsNbtags[ElMu][iExact1btag]->Fill(getNBTags(), Mll, EventWeight);
	fHDY_InvMass        [ElMu][iExact1btag]->Fill(             Mll, EventWeight);
      }
      if (getNBTags() == 2){
	fHDY_InvMassVsNPV   [ElMu][iExact2btag]->Fill(nGoodVertex, Mll, EventWeight);
	fHDY_InvMassVsMET   [ElMu][iExact2btag]->Fill(getMET()   , Mll, EventWeight);
	fHDY_InvMassVsNjets [ElMu][iExact2btag]->Fill(getNJets() , Mll, EventWeight);
	fHDY_InvMassVsNbtags[ElMu][iExact2btag]->Fill(getNBTags(), Mll, EventWeight);
	fHDY_InvMass        [ElMu][iExact2btag]->Fill(             Mll, EventWeight);
      }
    }
  }
  
  ResetHypLeptons(); 
  if (PassTriggerMuMu() && IsMuMuEvent()){
    
    EventWeight = gWeight * getSF(Muon);
    Mll = (fHypLepton1.p+fHypLepton2.p).M();
    
    if (PassesMllVeto() && PassesMuonEta2p1(Muon) && Passes3rdLeptonVeto()){
      fHDY_InvMassVsNPV   [Muon][iDilepton]->Fill(nGoodVertex, Mll, EventWeight);
      fHDY_InvMassVsMET   [Muon][iDilepton]->Fill(getMET()   , Mll, EventWeight);
      fHDY_InvMassVsNjets [Muon][iDilepton]->Fill(getNJets() , Mll, EventWeight);
      fHDY_InvMassVsNbtags[Muon][iDilepton]->Fill(getNBTags(), Mll, EventWeight);
      fHDY_InvMass        [Muon][iDilepton]->Fill(             Mll, EventWeight);
	
      if (PassesMETCut())   {
        fHDY_InvMassVsNPV   [Muon][iMET]->Fill(nGoodVertex, Mll, EventWeight);
        fHDY_InvMassVsMET   [Muon][iMET]->Fill(getMET()   , Mll, EventWeight);
        fHDY_InvMassVsNjets [Muon][iMET]->Fill(getNJets() , Mll, EventWeight);
        fHDY_InvMassVsNbtags[Muon][iMET]->Fill(getNBTags(), Mll, EventWeight);
        fHDY_InvMass	    [Muon][iMET]->Fill( 	    Mll, EventWeight);
	if (getNBTags() == 1){
	  fHDY_InvMassVsNPV   [Muon][iExact1btag]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [Muon][iExact1btag]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [Muon][iExact1btag]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[Muon][iExact1btag]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [Muon][iExact1btag]->Fill(             Mll, EventWeight);
	}
	if (getNBTags() == 2){
	  fHDY_InvMassVsNPV   [Muon][iExact2btag]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [Muon][iExact2btag]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [Muon][iExact2btag]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[Muon][iExact2btag]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [Muon][iExact2btag]->Fill(             Mll, EventWeight);
	}
        if (PassesNJetsCut()) {
	  fHDY_InvMassVsNPV   [Muon][i2jets]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [Muon][i2jets]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [Muon][i2jets]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[Muon][i2jets]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [Muon][i2jets]->Fill(             Mll, EventWeight);
	  
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

  ResetHypLeptons(); 
  if (PassTriggerEE()   && IsElElEvent()){
    EventWeight = gWeight * getSF(Elec);
    Mll = (fHypLepton1.p+fHypLepton2.p).M();
    
    if (PassesMllVeto() && PassesMuonEta2p1(Elec) && Passes3rdLeptonVeto()){
      fHDY_InvMassVsNPV   [Elec][iDilepton]->Fill(nGoodVertex, Mll, EventWeight);
      fHDY_InvMassVsMET   [Elec][iDilepton]->Fill(getMET()   , Mll, EventWeight);
      fHDY_InvMassVsNjets [Elec][iDilepton]->Fill(getNJets() , Mll, EventWeight);
      fHDY_InvMassVsNbtags[Elec][iDilepton]->Fill(getNBTags(), Mll, EventWeight);
      fHDY_InvMass        [Elec][iDilepton]->Fill(             Mll, EventWeight);
	
      if (PassesMETCut())   {
        fHDY_InvMassVsNPV   [Elec][iMET]->Fill(nGoodVertex, Mll, EventWeight);
        fHDY_InvMassVsMET   [Elec][iMET]->Fill(getMET()   , Mll, EventWeight);
        fHDY_InvMassVsNjets [Elec][iMET]->Fill(getNJets() , Mll, EventWeight);
        fHDY_InvMassVsNbtags[Elec][iMET]->Fill(getNBTags(), Mll, EventWeight);
        fHDY_InvMass	    [Elec][iMET]->Fill( 	    Mll, EventWeight);
	
	if (getNBTags() == 1){
	  fHDY_InvMassVsNPV   [Elec][iExact1btag]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [Elec][iExact1btag]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [Elec][iExact1btag]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[Elec][iExact1btag]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [Elec][iExact1btag]->Fill(             Mll, EventWeight);
	}
	if (getNBTags() == 2){
	  fHDY_InvMassVsNPV   [Elec][iExact2btag]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [Elec][iExact2btag]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [Elec][iExact2btag]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[Elec][iExact2btag]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [Elec][iExact2btag]->Fill(             Mll, EventWeight);
	}
        if (PassesNJetsCut()) {
	  fHDY_InvMassVsNPV   [Elec][i2jets]->Fill(nGoodVertex, Mll, EventWeight);
	  fHDY_InvMassVsMET   [Elec][i2jets]->Fill(getMET()   , Mll, EventWeight);
	  fHDY_InvMassVsNjets [Elec][i2jets]->Fill(getNJets() , Mll, EventWeight);
	  fHDY_InvMassVsNbtags[Elec][i2jets]->Fill(getNBTags(), Mll, EventWeight);
	  fHDY_InvMass        [Elec][i2jets]->Fill(             Mll, EventWeight);
	  
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
  ResetHypLeptons();
}
void TreeAnalysisTop::FillKinematicHistos(gChannel chan, iCut cut){
#ifdef DEBUG
  cout << "Filling KinematicHistos("<<chan<<","<<cut<<")... ";
  cout << fHypLepton1.index << " , " << fHypLepton2.index << endl;
#endif
  
  if (gSysSource != Norm)      return;  //only fill histograms for nominal distributions...
  if (fChargeSwitch == true  ) return;
  //  if (fHypLepton1.index == -1) return;
  //  if (fHypLepton2.index == -1) return;
  
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

  int ib = getLeadingJetbTag();
  if (ib>=0)
    fHCSVTag[chan][cut] ->Fill(T_JetAKCHS_Tag_CombSVtx->at(ib), EventWeight);
//  if (njets > 1)
//    fHCSVTag[chan][cut] ->Fill(T_JetAKCHS_Tag_CombSVtx->at(Jet[1].index), EventWeight);

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
  
  //fHAbsDelPhiLep[chan][cut] ->Fill(TMath::Abs(fHypLepton1.p.DeltaPhi((fHypLepton2.p)))/TMath::Pi(), EventWeight);
  fHvertices[chan][cut]     ->Fill(nGoodVertex, EventWeight);
#ifdef __ISSTOP  
  if(gSampleName == "T2tt_150to250LSP1to100_LeptonFilter"){
    for (size_t t=0; t<T_Gen_StopMass->size(); t++)
      fHStopMass[chan][cut] ->Fill(T_Gen_StopMass->at(t), EventWeight);
    for (size_t t=0; t<T_Gen_Chi0Mass->size(); t++)
      fHChi0Mass[chan][cut] ->Fill(T_Gen_Chi0Mass->at(t), EventWeight);
      fHChi0StopMass[chan][cut] ->Fill(T_Gen_StopMass->at(0), T_Gen_Chi0Mass->at(0), EventWeight);
  }
#endif
#ifdef DEBUG
  cout << " DONE!" << endl;
#endif
}
void TreeAnalysisTop::FillYieldsHistograms(gChannel chan, iCut cut, gSystFlag sys){
#ifdef DEBUG
  cout << "FillYieldsHistograms("<<chan<<","<<cut<<","<<sys<<")...";
#endif
  if (fChargeSwitch){   fHSSyields[chan][sys]->Fill(cut, EventWeight);  }
  else {                fHyields[chan][sys]  ->Fill(cut, EventWeight);  }
  
#ifdef __ISPDF
  if (cut==i1btag && chan==ElMu && sys==Norm) {
    for (int i=0; i<52; i++) {
      fHpdfWeight->Fill(i, EventWeight*T_Event_pdfWeight->at(i));
    }
  }
#endif

  /// FOR SYSTEMATIC STUDIES
  int njets  = 0; njets  = getNJets();
  int nbtags = 0; nbtags = getNBTags();
  
  if (fChargeSwitch) { 
  //cout << "!!!!!!!"<< endl; 
    fHSSInvMass[chan][cut][sys]->Fill((fHypLepton1.p+fHypLepton2.p).M(), EventWeight);
    if (njets == 0) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags,        EventWeight);
    if (njets == 1) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+1,      EventWeight);
    if (njets == 2) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+3,      EventWeight);
    if (njets == 3) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+6,      EventWeight);
    if (njets >= 4) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+10,     EventWeight);
    fHSSAbsDelPhiLeps[chan][cut][sys]->Fill( TMath::Abs(fHypLepton1.p.DeltaPhi(fHypLepton2.p))/TMath::Pi(), EventWeight);
    if(njets >= 2) fHSSdelPhi2LeadJets[chan][cut][sys]->Fill( TMath::Abs(Jet[0].p.DeltaPhi(Jet[1].p))/TMath::Pi(), EventWeight);
    if (njets > 1) {
       double deltaR_temp1 = 999., deltaR_temp2 = 999.;       
       if( fHypLepton1.p.DeltaR(Jet[0].p) < fHypLepton2.p.DeltaR(Jet[0].p) ) {
          deltaR_temp1 = fHypLepton1.p.DeltaR(Jet[0].p);
          deltaR_temp2 = fHypLepton2.p.DeltaR(Jet[1].p);  
       }else{
          deltaR_temp1 = fHypLepton2.p.DeltaR(Jet[0].p);
          deltaR_temp2 = fHypLepton1.p.DeltaR(Jet[1].p);  
       }
       fHSSminDelRJetsLeps[chan][cut][sys]->Fill(TMath::Min(deltaR_temp1, deltaR_temp2), EventWeight);
    }
  }
  else {
    fHInvMass[chan][cut][sys]->Fill((fHypLepton1.p+fHypLepton2.p).M(), EventWeight);
    if (njets == 0) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags,        EventWeight);
    if (njets == 1) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+1,      EventWeight);
    if (njets == 2) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+3,      EventWeight);
    if (njets == 3) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+6,      EventWeight);
    if (njets >= 4) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+10,     EventWeight);
    fHAbsDelPhiLeps[chan][cut][sys]->Fill( TMath::Abs(fHypLepton1.p.DeltaPhi(fHypLepton2.p))/TMath::Pi(), EventWeight);
    if( njets >= 2) fHdelPhi2LeadJets[chan][cut][sys]->Fill( TMath::Abs(Jet[0].p.DeltaPhi(Jet[1].p))/TMath::Pi(), EventWeight);
    if (njets > 1) {
       double deltaR_temp1 = 999., deltaR_temp2 = 999.;       
       if( fHypLepton1.p.DeltaR(Jet[0].p) < fHypLepton2.p.DeltaR(Jet[0].p) ) {
          deltaR_temp1 = fHypLepton1.p.DeltaR(Jet[0].p);
          deltaR_temp2 = fHypLepton2.p.DeltaR(Jet[1].p);  
       }else{
          deltaR_temp1 = fHypLepton2.p.DeltaR(Jet[0].p);
          deltaR_temp2 = fHypLepton1.p.DeltaR(Jet[1].p);  
       }
       fHminDelRJetsLeps[chan][cut][sys]->Fill(TMath::Min(deltaR_temp1, deltaR_temp2), EventWeight);
    }
  }
  
  if (!gIsData){
    fHLepSys[chan][cut] ->Fill(getLeptonError(chan), EventWeight);
    fHTrigSys[chan][cut]->Fill(getTriggerError(chan),EventWeight);
    //    // FOR SS ORIGINS
    //    if (fChargeSwitch) fHSSOrigins[chan][cut]->Fill();
    //    else               fHOrigins[chan][cut]  ->Fill();
  }
#ifdef DEBUG
  cout << " DONE! " << endl;
#endif
  return;
}
void TreeAnalysisTop::FillYields(gSystFlag sys){
#ifdef DEBUG
  cout << "FillYields("<<sys<<")... ";
#endif
  ResetHypLeptons();  
  if (gDoDF && PassTriggerEMu()  && IsElMuEvent()){
    // Define Hypothesis Leptons...
    EventWeight = gWeight * getSF(ElMu) * getTopPtSF();
    hWeight -> Fill(EventWeight,1.);
#ifdef DEBUG
  cout << " pass trigger + emu, ";
#endif

#ifdef __ISSTOP
    if(gSampleName == "T2tt_150to250LSP1to100_LeptonFilter")
      EventWeight = EventWeight * T_Gen_polWeights->at(10);
#endif
#ifdef __ISMCNLO
    // 0.115 = Fraction events with negative weight
    EventWeight = EventWeight * T_Event_weight /(abs(T_Event_weight)*(1.-2.*0.115)); 
#endif

    if (PassesMllVeto() && PassesMuonEta2p1(ElMu) && Passes3rdLeptonVeto()){
#ifdef DEBUG
      cout << " pass mll, ";
#endif
      FillYieldsHistograms(ElMu, iDilepton, sys);
      if(sys==Norm) FillKinematicHistos(ElMu,iDilepton);
      
      FillYieldsHistograms(ElMu, iZVeto, sys);      
      if(sys==Norm) FillKinematicHistos(ElMu, iZVeto);
      
      FillYieldsHistograms(ElMu, iMET, sys);      
      if(sys==Norm) FillKinematicHistos(ElMu,iMET);
      
      if (PassesNJetsCut()) {
#ifdef DEBUG
	cout << " pass njets with njets = "<<getNJets()<<", ";
#endif
	FillYieldsHistograms(ElMu, i2jets, sys);      
	if(sys==Norm) FillKinematicHistos(ElMu,i2jets);
	if (PassesNBtagCut()) {
#ifdef DEBUG
	  cout << " pass nbjets with nbtags = "<<getNBTags()<<", ";
#endif
	  FillYieldsHistograms(ElMu, i1btag, sys);      
	  if(sys==Norm) FillKinematicHistos(ElMu,i1btag);
	}
      }
      if (getNBTags() == 1){
#ifdef DEBUG
	cout << " pass nbjets=1";
#endif
	FillYieldsHistograms(ElMu, iExact1btag, sys);      
	if(sys==Norm) FillKinematicHistos(ElMu,iExact1btag);
      }
      if (getNBTags() == 2){
#ifdef DEBUG
	cout << " pass nbjets=2";
#endif
	FillYieldsHistograms(ElMu, iExact2btag, sys);      
	if(sys==Norm) FillKinematicHistos(ElMu,iExact2btag);
      }
    }
  }
  //  if(!gDoSF) return;
  
  ResetHypLeptons(); 
  if (gDoSF && PassTriggerMuMu() && IsMuMuEvent()){
    EventWeight = gWeight * getSF(Muon)  * getTopPtSF();
#ifdef __ISSTOP
    if(gSampleName == "T2tt_150to250LSP1to100_LeptonFilter")
      EventWeight = EventWeight * T_Gen_polWeights->at(10);
#endif
#ifdef __ISMCNLO
    // 0.115 = Fraction events with negative weight
    EventWeight = EventWeight * T_Event_weight /(abs(T_Event_weight)*(1.-2.*0.115)); 
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
	  
	  if (getNBTags() == 1){
	    FillYieldsHistograms(Muon, iExact1btag, sys);      
	    if(sys==Norm) FillKinematicHistos(Muon,iExact1btag);
	  }
	  if (getNBTags() == 2){
	    FillYieldsHistograms(Muon, iExact2btag, sys);      
	    if(sys==Norm) FillKinematicHistos(Muon,iExact2btag);
	  }
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

  ResetHypLeptons(); 
  if (gDoSF && PassTriggerEE()   && IsElElEvent()){
    EventWeight = gWeight * getSF(Elec) * getTopPtSF();     
#ifdef __ISSTOP
    if(gSampleName == "T2tt_150to250LSP1to100_LeptonFilter")
      EventWeight = EventWeight * T_Gen_polWeights->at(10);
#endif
#ifdef __ISMCNLO
    // 0.115 = Fraction events with negative weight
    EventWeight = EventWeight * T_Event_weight /(abs(T_Event_weight)*(1.-2.*0.115)); 
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
	  
	  if (getNBTags() == 1){
	    FillYieldsHistograms(Elec, iExact1btag, sys);      
	    if(sys==Norm) FillKinematicHistos(Elec,iExact1btag);
	  }
	  if (getNBTags() == 2){
	    FillYieldsHistograms(Elec, iExact2btag, sys);      
	    if(sys==Norm) FillKinematicHistos(Elec,iExact2btag);
	  }
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
  ResetHypLeptons();
#ifdef DEBUG
  cout << " DONE!"<<endl;
#endif
}
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
  return true; // don't apply third lepton veto...
  
  // Return false if there are not 2 signal leptons
  if (fHypLepton1.index == -1) return false;
  if (fHypLepton2.index == -1) return false;  
  
  //  Int_t nvetoleptons = 0;
  for(UInt_t i = 0; i < T_Muon_Pt->size(); ++i){
    if (fHypLepton1.index > -1 && (UInt_t)fHypLepton1.index == i && fHypLepton1.type == 0) continue;
    if (fHypLepton2.index > -1 && (UInt_t)fHypLepton2.index == i && fHypLepton2.type == 0) continue;
    if (IsVetoMuon(i)) return false;
    //    nvetoleptons++;
  }
  
  for(UInt_t i = 0; i < T_Elec_Pt->size(); ++i){
    if (fHypLepton1.index > -1 && (UInt_t)fHypLepton1.index == i && fHypLepton1.type == 1) continue;
    if (fHypLepton2.index > -1 && (UInt_t)fHypLepton2.index == i && fHypLepton2.type == 1) continue;
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
  //  cout << "[DEBUG]: calling getNBTags from PassesNBtagCut:" << getNBTags() << endl;
  if (getNBTags() < 1) return false;
  
  return true;
}

bool TreeAnalysisTop::IsElMuEvent(){
  // if fChargeSwitch, pick SS event.
  if (fChargeSwitch){      return (IsDileptonEvent()  == 3);   }

  return (IsDileptonEvent() == -3);
}
bool TreeAnalysisTop::IsMuMuEvent(){
  // if fChargeSwitch, pick SS event.
  if (fChargeSwitch){  return (IsDileptonEvent()  == 1); }
  
  return (IsDileptonEvent() == -1);
}
bool TreeAnalysisTop::IsElElEvent(){
  // if fChargeSwitch, pick SS event.
  if (fChargeSwitch){    return (IsDileptonEvent()  == 2); }
  
  return (IsDileptonEvent() == -2);
}
int TreeAnalysisTop::IsDileptonEvent(){
  // bool(TreeAnalysisTop::*muonSelector)(unsigned int), 
  // int &ind2, bool(TreeAnalysisTop::*eleSelector)(unsigned int)){
  // Looks for a pair of leptons with given object selectors
  // Return the channel: 0 = none found
  //                1 / -1 = mu+mu+ (SS) / mu-mu+ (OS) pair
  //                2 / -2 = e+e+   (SS) / e-e+   (OS) pair
  //                3 / -3 = mu+e+  (SS) / mu-e+  (OS) pair
#ifdef DEBUG
  cout << "IsDileptonEvent(): NLeptons =" << Lepton.size()<<", ";
#endif
  // Check for at least one tight pair
  if(Lepton.size() < 2) return 0;
  
  // Proceed to select pair with highest pt
  int select = Lepton[0].charge*Lepton[1].charge;
  
  /////////////////////////////////////////////////////////////////////////
  int result = 0;
  if (Lepton[0].type == 0 && Lepton[1].type == 0) result = 1; // mu/mu
  if (Lepton[0].type == 1 && Lepton[1].type == 1) result = 2; // el/el
  if (Lepton[0].type == 0 && Lepton[1].type == 1) result = 3; // mu/el
  if (Lepton[0].type == 1 && Lepton[1].type == 0) result = 3; // mu/el
  
  // Return values, assigning indexes (first is always a muon, second electron):
  if (result == 3) {
    if      (Lepton[0].type == 0 && Lepton[1].type == 1) { 
      fHypLepton1 = lepton(Lepton[0]);
      fHypLepton2 = lepton(Lepton[1]);
    }
    else if (Lepton[0].type == 1 && Lepton[1].type == 0){
      fHypLepton1 = lepton(Lepton[1]);
      fHypLepton2 = lepton(Lepton[0]);
    }
  }
  else {
    fHypLepton1 = lepton(Lepton[0]);
    fHypLepton2 = lepton(Lepton[1]);
  } 
  
  result *= select; // Add charge to result
  
#ifdef DEBUG
  cout << result;
  cout << " DONE!" << endl;
#endif
  return result;
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
int TreeAnalysisTop::getSelectedLeptons(){
  // Loops over the total number of Muons and Electrons and returns
  // the Number of Leptons.
  if (Lepton.size() > 0) {
    cout << "[WARNING]: you have called this function previously... RESETTING..."<<endl;
    Lepton.clear();
  }
  vector<lepton> tmp_lepton;

  nMuon = 0;
  TLorentzVector lep;
  for (UInt_t i=0; i<T_Muon_Pt->size();i++){
    if (IsTightMuon(i) == false) continue;
    lep.SetPxPyPzE(MuPx.at(i), MuPy.at(i), T_Muon_Pz->at(i), T_Muon_Energy->at(i));
    lepton tmpLepton(lep,T_Muon_Charge->at(i), 0, i);
    tmp_lepton.push_back(tmpLepton);
    nMuon++;
  }
  
  nElec = 0;
  for (UInt_t i=0; i<T_Elec_Pt->size();i++){
    if (IsTightElectron(i) == false) continue;
    lep.SetPxPyPzE(ElPx.at(i), ElPy.at(i), T_Elec_Pz->at(i), T_Elec_Energy->at(i));
    lepton tmpLepton(lep,T_Elec_Charge->at(i), 1, i);
    tmp_lepton.push_back(tmpLepton);
    nElec++;
  }
  
  Lepton = SortLeptonsByPt(tmp_lepton);
  return Lepton.size();
}
//------------------------------------------------------------------------------
// Muon Selectors
//------------------------------------------------------------------------------
bool TreeAnalysisTop::IsVetoMuon(unsigned int iMuon, float ptcut){
  TLorentzVector lep;
  lep.SetPxPyPzE(MuPx.at(iMuon), MuPy.at(iMuon), T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  if (lep.Pt() < ptcut)            return false;
  if (TMath::Abs(lep.Eta()) > 2.4) return false;
  
  if (T_Muon_IsGlobalMuon->at(iMuon) == 0 && T_Muon_IsTrackerMuonArbitrated->at(iMuon) == 0) return false;
  
  float relIso = (T_Muon_chargedHadronIsoR04->at(iMuon) + max(0.0 , T_Muon_neutralHadronIsoR04->at(iMuon) + T_Muon_photonIsoR04->at(iMuon)- 0.5*T_Muon_sumPUPtR04->at(iMuon)))/lep.Pt();
  if (relIso > 0.20) return false;

  return true;
}
bool TreeAnalysisTop::IsTightMuon(unsigned int iMuon,float ptcut){
  TLorentzVector lep;
  lep.SetPxPyPzE(MuPx.at(iMuon), MuPy.at(iMuon), T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  if (lep.Pt() < ptcut)            return false;
  if (TMath::Abs(lep.Eta()) > 2.4) return false;
  
  // POG Tight Muons definition				       
  if (!T_Muon_IsGlobalMuon->at(iMuon))                             return false;
  if (T_Muon_NormChi2GTrk->at(iMuon) >= 10.)                       return false;
  if (T_Muon_NValidHitsGTrk->at(iMuon) < 1)                        return false;
  //this is still not the exact same def.		       
  if (T_Muon_NumOfMatchedStations->at(iMuon) <= 1)                 return false; 
  //							       
  if (TMath::Abs(T_Muon_IPwrtAveBSInTrack->at(iMuon)) >= 0.2)      return false; 
  if (TMath::Abs(T_Muon_vz->at(iMuon) - T_Vertex_z->at(0)) >= 0.5) return false;
  if (T_Muon_NValidPixelHitsInTrk->at(iMuon) == 0)                 return false;
  if (T_Muon_NLayers->at(iMuon) <= 5)                              return false;

  float relIso = getMuonIso(iMuon);
  
  if (relIso > 0.12) return false;
  
  return true;
}
float TreeAnalysisTop::getMuonIso(int iMuon){
  if (iMuon < 0) return 9999.;
  if (iMuon >= (int)T_Muon_chargedHadronIsoR04->size()) return 9999.;

  TLorentzVector lep;
  lep.SetPxPyPzE(MuPx.at(iMuon), MuPy.at(iMuon), T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  
  return (T_Muon_chargedHadronIsoR04->at(iMuon) + max(0.0 , T_Muon_neutralHadronIsoR04->at(iMuon) + T_Muon_photonIsoR04->at(iMuon)- 0.5*T_Muon_sumPUPtR04->at(iMuon)))/lep.Pt();
}
//------------------------------------------------------------------------------
// Electron Selectors
//------------------------------------------------------------------------------
bool TreeAnalysisTop::IsVetoElectron(unsigned int iElec,float ptcut){
  TLorentzVector lep;
  lep.SetPxPyPzE(ElPx.at(iElec), ElPy.at(iElec), T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));

  if (lep.Pt() < ptcut)               return false;
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
  /* float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
     if      (sceta < 0.8                    ) { if (T_Elec_MVA->at(iElec) > 0.94) return true; }
     else if (sceta >= 0.8   && sceta < 1.479) { if (T_Elec_MVA->at(iElec) > 0.85) return true; }
     else if (sceta >= 1.479 && sceta < 2.5  ) { if (T_Elec_MVA->at(iElec) > 0.92) return true; }
  */
  if (T_Elec_MVA->at(iElec) > 0.90) return true;
  return false;
}
bool TreeAnalysisTop::IsTightElectron(unsigned int iElec, float ptcut){
  TLorentzVector lep;
  lep.SetPxPyPzE(ElPx.at(iElec), ElPy.at(iElec), T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));
  
  float pt  = lep.Pt();
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if (sceta > 1.4442 && sceta < 1.566) return false;
  if (lep.Pt() < ptcut)                return false;
  if (TMath::Abs(lep.Eta()) > 2.5)     return false;
  
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
  
  if (!passTriggerID)                        return false;
  if (!IsMVAIDElectron(iElec))               return false;
  if (!T_Elec_passConversionVeto->at(iElec)) return false;
  if ( T_Elec_nHits->at(iElec) >= 1)         return false;
  
  // ISOLATION
  float relIso =  getElecIso(iElec);

  if (relIso > 0.1) return false;
  
  return true;
}
float TreeAnalysisTop::getElecIso(int iElec){
  if (iElec < 0) return 9999.;
  if (iElec >= (int)T_Elec_chargedHadronIso->size()) return 9999.;
  
  TLorentzVector lep;
  lep.SetPxPyPzE(ElPx.at(iElec), ElPy.at(iElec), T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));
  float pt     = lep.Pt();
  float EA     = getEACorrection(T_Elec_Eta->at(iElec));
  float relIso = (T_Elec_chargedHadronIso->at(iElec) + 
		  max((float)0.0, 
		      T_Elec_neutralHadronIso->at(iElec) + 
		      T_Elec_photonIso->at(iElec) - 
		      T_Event_RhoIso*EA
		      )
		  )/pt;

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
  return nJets;
}
float TreeAnalysisTop::getDRClosestJet(TLorentzVector lep){
  float minDR = 9999.;
  for (unsigned int i=0; i<Jet.size(); i++) {
    if (minDR > lep.DeltaR(Jet[i].p)) minDR = lep.DeltaR(Jet[i].p);
  }
  
  return minDR;
}
float TreeAnalysisTop::getDPhiClosestJet(TLorentzVector lep){
  float minDphi = 9999.;
  for (unsigned int i=0; i<Jet.size(); i++) {
    if (minDphi > TMath::Abs(lep.DeltaPhi(Jet[i].p))) minDphi = TMath::Abs(lep.DeltaPhi(Jet[i].p));
  }
  
  return minDphi;
}
int TreeAnalysisTop::getLeadingJetbTag(){
  for (unsigned int i=0; i<Jet.size(); i++) {
    if (Jet[i].isbtag) return i;
  }
  
  return  -1;
}
int TreeAnalysisTop::getNBTags(){
  int ntags(0);
  
  for(UInt_t i = 0; i <Jet.size(); i++){
    if (Jet[i].isbtag) ntags++;
  }
  
  return ntags;
}
float TreeAnalysisTop::getDeltaPhillJet(){
  if (fHypLepton1.index == -1) return -999.;
  if (fHypLepton2.index == -1) return -999.;
  Int_t ij = getLeadingJetbTag();
  if (ij < 0) return -999.; 
  TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;
  TLorentzVector jet = Jet[ij].p; 
  return TMath::Abs(dilep.DeltaPhi(jet));
}
float TreeAnalysisTop::getTopD(){
  if (fHypLepton1.index == -1) return -999;
  if (fHypLepton2.index == -1) return -999;

  // Make Dilepton system
  TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;

  Float_t DeltaPhi(0.),TopD(0.);
  TLorentzVector jet;
  if (nJets == 0) return -999.;
   
  DeltaPhi = TMath::Abs(dilep.DeltaPhi(Jet[0].p));
  TopD     = 1 - (DeltaPhi/TMath::Pi()) * (1 - T_JetAKCHS_Tag_CombSVtx->at(Jet[0].index));
  
  return TopD;
}
int TreeAnalysisTop::getSelectedJets(){
  int nj(0);
  if (Jet.size() > 0) {
    cout << "[WARNING]: you have called this function previously, RESETTING..."<<endl;
    Jet.clear();
  }

  int btagSys = 0;
  TLorentzVector jt;
  for (UInt_t i=0; i<T_JetAKCHS_Energy->size(); i++) {
    if(!IsGoodJet(i,gJetEtCut)) continue;
    jt.SetPtEtaPhiE(JetEt.at(i),T_JetAKCHS_Eta->at(i),JetPhi.at(i),T_JetAKCHS_Energy->at(i));
    bool isbtag = false;
    if (gIsData) {
      isbtag = fBTagSF->IsTagged(T_JetAKCHS_Tag_CombSVtx->at(i),-999999,JetEt.at(i),T_JetAKCHS_Eta->at(i),btagSys);
    }
    else {
      if(TMath::Abs(T_JetAKCHS_Parton_Flavour->at(i)) == 5 || TMath::Abs(T_JetAKCHS_Parton_Flavour->at(i)) == 4){
	if (gSysSource == BtagUp)     btagSys =  1;
	if (gSysSource == BtagDown)   btagSys = -1;
	if (gSysSource == MisTagUp)   btagSys =  0;
	if (gSysSource == MisTagDown) btagSys =  0;
      }
      else {
//if(TMath::Abs(T_JetAKCHS_Parton_Flavour->at(i)) != 5 && TMath::Abs(T_JetAKCHS_Parton_Flavour->at(i)) != 4){
	if (gSysSource == BtagUp)     btagSys =  0;
	if (gSysSource == BtagDown)   btagSys =  0;
	if (gSysSource == MisTagUp)   btagSys =  1;
	if (gSysSource == MisTagDown) btagSys = -1;
      }
      isbtag = fBTagSF->IsTagged(T_JetAKCHS_Tag_CombSVtx->at(i), T_JetAKCHS_Parton_Flavour->at(i), 
				 JetEt.at(i),T_JetAKCHS_Eta->at(i), btagSys);
    }
    jet tmpjet(jt, isbtag, i);
    Jet.push_back(tmpjet);
    nj++;
  }
  return nj;
}
bool TreeAnalysisTop::IsGoodJet(unsigned int ijet, float ptcut){
  float minDR = 0.5;
  if(JetEt.at(ijet) < ptcut)                     return false;
  if(TMath::Abs(T_JetAKCHS_Eta->at(ijet)) > 2.4) return false; // btagging only up to 2.4
  
  // JetID 
  if ( !(T_JetAKCHS_nDaughters->at(ijet)        > 1   ) ) return false;
  if ( !(T_JetAKCHS_NeutHadEnergyFrac->at(ijet) < 0.99) ) return false;
  if ( !(T_JetAKCHS_NeutEmEnergyFrac ->at(ijet) < 0.99) ) return false;
  if (TMath::Abs(T_JetAKCHS_Eta->at(ijet)) < 2.5){
    if ( !(T_JetAKCHS_CharEmEnergyFrac->at(ijet)  < 0.99) ) return false;
    if ( !(T_JetAKCHS_CharHadEnergyFrac->at(ijet) > 0.00) ) return false;
    if ( !(T_JetAKCHS_ChargedMultiplicity->at(ijet) > 0 ) ) return false;
  }
  
  // Remove jets close to all selected leptons... 
  TLorentzVector jet;
  jet.SetPtEtaPhiE(JetEt.at(ijet),T_JetAKCHS_Eta->at(ijet),JetPhi.at(ijet),T_JetAKCHS_Energy->at(ijet));
  
  for(unsigned int i = 0; i < Lepton.size(); ++i){
    if(jet.DeltaR(Lepton[i].p) < minDR) return false;
  }
  return true;
}

//------------------------------------------------------------------------------
// SelectedGenLepton
//------------------------------------------------------------------------------
void TreeAnalysisTop::SelectedGenLepton() {

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
    tmp.SetPtEtaPhiE(JetEt.at(*it),  T_JetAKCHS_Eta->at(*it), 
		     JetPhi.at(*it), T_JetAKCHS_Energy->at(*it));         // set temp to the jet
    ojets += tmp;                                                         // add jet to the old jets vector
    if(flag == 1) JetEt.at(*it) *= (1 + T_JetAKCHS_Uncertainty->at(*it)); // vary up for flag 1
    if(flag == 2) JetEt.at(*it) *= (1 - T_JetAKCHS_Uncertainty->at(*it)); // vary down for flag 2;
    if(flag == 3){

      //TVector3 genJet(T_JetAKCHS_GenJet_Px->at(*it),T_JetAKCHS_GenJet_Py->at(*it),T_JetAKCHS_GenJet_Pz->at(*it)); 
      //if (genJet.Pt() < 15) continue; 
      //if (genJet.DeltaR(tmp) < 0.5) continue;
      
      //float jerScale = getJERScale(*it);
      //else {
	// get the resolution
	float sigmaMC  = getErrPt(JetEt.at(*it), T_JetAKCHS_Eta->at(*it))/JetEt.at(*it);
	float jerScale = getJERScale(*it);                                  // get JER scale factors 
      
	float factor = fRand3->Gaus(1., sqrt(jerScale*jerScale -1.)*sigmaMC );
	JetEt.at(*it) = JetEt.at(*it) * factor;           // smear for flag 3
      //}
    }
    // set tmp to the scaled/smeared jet
    tmp.SetPtEtaPhiE(JetEt.at(*it), T_JetAKCHS_Eta->at(*it),JetPhi.at(*it), T_JetAKCHS_Energy->at(*it)); 
    jets += tmp;                                                    // add scaled/smeared jet to the new jets
  }
  propagateMET(jets, ojets);                                        // propagate this change to the MET
}
void TreeAnalysisTop::ScaleLeptons(int flag){
  // Shift the lepton pts for systematics studies
  if(gIsData) return; // don't smear data
  
  if(flag == 0) return;
  float scale = 0.002; // 0.2%
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



