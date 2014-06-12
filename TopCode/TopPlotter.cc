#include "TopPlotter.h"

#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLatex.h"

using namespace std;

bool gUseTTMadSpin = true;
const float fLumiNorm = 19664.225;

TopPlotter::TopPlotter(){}

void TopPlotter::Init(TString pathtofiles){
  cout << "ResetDataMembers" << endl;
  ResetDataMembers();
  //  ttbar_TLWG = 252.89;
  ttbar_TLWG = 245.8;
  
  LoadSamples(pathtofiles);
}
void TopPlotter::LoadSamples(TString pathtofiles){
  TFile *_file ;
  TH1F *hOSyields;
  TH1F *hOSyields_sys;
  TH1F *hSSyields;
  
  cout << "Loading Samples.... " << endl;
  Float_t Weight = fLumiNorm; 
  for (size_t sample=0; sample<gNSAMPLES; sample++){
    TString samplename = pathtofiles + "/Tree_Legacy_"+SampleName[sample]+".root";
    cout << "Reading " + samplename +" ..." << endl;
    _file = new TFile(samplename);
    Bool_t IsData = (sample == DoubleElectron || sample == DoubleMu || sample == MuEG);
    
    S[sample].name = SampleName[sample];
    
    // CALCULATE WEIGHT
    if (sample==TTJets_MadSpin    || 
	sample==TTJets_matchingup || sample==TTJets_matchingdown ||
	sample==TTJets_scaleup    || sample==TTJets_scaledown      )   
      Weight = fLumiNorm * (0.108*9.)*(0.108*9.);
    else 
      Weight = fLumiNorm; 
    
    // Load numbers... 
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      hOSyields = (TH1F*) _file->Get("H_Yields_"  +gChanLabel[chan]);
      hSSyields = (TH1F*) _file->Get("H_SSYields_"+gChanLabel[chan]);
      
      for (size_t cut=0; cut<iNCUTS; cut++){
	if (IsData){
	  S[sample].Yields       [chan][cut] = hOSyields->GetBinContent(cut+1);
	  S[sample].Yields_stat  [chan][cut] = hOSyields->GetBinError(cut+1);
	  S[sample].SSYields     [chan][cut] = hSSyields->GetBinContent(cut+1);
	  S[sample].SSYields_stat[chan][cut] = hSSyields->GetBinError(cut+1);
	}
	else {  
	  S[sample].Yields       [chan][cut] = hOSyields->GetBinContent(cut+1) * Weight;
	  S[sample].Yields_stat  [chan][cut] = hOSyields->GetBinError(cut+1)   * Weight;
	  S[sample].SSYields     [chan][cut] = hSSyields->GetBinContent(cut+1) * Weight;
	  S[sample].SSYields_stat[chan][cut] = hSSyields->GetBinError(cut+1)   * Weight;
	}
      }
      // Load Systematics (ONLY FOR MC...) 
      if (!IsData){
	for (size_t sys=1; sys<gNSYST; sys++){
	  hOSyields_sys = (TH1F*) _file->Get("H_Yields_"  +gChanLabel[chan]+"_"+SystName[sys]);
	  for (size_t cut=0; cut<iNCUTS; cut++)
	    S[sample].Yields_syst[chan][cut][sys] = hOSyields_sys->GetBinContent(cut+1) * Weight;
	}
      }
    }
    
    // Load kinematic histograms of the samples. 
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      
      if (SampleName[sample] == "DoubleMu"        ||       
	  SampleName[sample] == "DoubleElectron"  || 
	  SampleName[sample] == "MuEG"            ||       
	  SampleName[sample] == "DYJets_Madgraph" ||
	  SampleName[sample] == "ZJets_Madgraph")  {
	
	/// LOAD DY DD estimation:
	TString histoname = "H_DY_InvMassVsNjets_" + gChanLabel[chan] + "_dilepton";
	S[sample].MllHistos[chan] = (TH1F*) (GetHisto2D(samplename, histoname))->ProjectionY();
	if (!IsData)                S[sample].MllHistos[chan]->Scale(Weight);
      }

      for (size_t cut=0; cut<iNCUTS; cut++){
	TString histoname = "";
	for (size_t var=0; var<gNVARS; var++){
	  histoname = "H_" + KinVarName[var] + "_" + gChanLabel[chan] + "_" + sCut[cut];
	  S[sample].KinHistos[chan][cut][var] = GetHisto1D(samplename, histoname);
	  if (!IsData)                S[sample].KinHistos[chan][cut][var]->Scale(Weight);
	}
	
	// SYSTEMATIC UNCERTAINTIES...
	for (size_t sys=0; sys<gNSYST; sys++){
	  if (sys==0) histoname = "H_NBtagsNJets_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "H_NBtagsNJets_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  //	  cout << "Reading... " << histoname << endl;
	  S[sample].SysHistos[chan][cut][sys] = GetHisto1D(samplename, histoname);
	  if (!IsData) S[sample].SysHistos[chan][cut][sys]->Scale(Weight);
	  
	  if (sys==0) histoname = "HSS_NBtagsNJets_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "HSS_NBtagsNJets_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  
	  S[sample].SSSysHistos[chan][cut][sys] = GetHisto1D(samplename, histoname);
	  if (!IsData) S[sample].SSSysHistos[chan][cut][sys]->Scale(Weight);
	}
      }
    }
  }
  LoadCategories();
}

void TopPlotter::LoadCategories(){
  
  // Now Load everything (Yields, SSYields and KinHistos) to the different bkg. categories...
  for (size_t cut=0; cut<iNCUTS; cut++){
    Data.name = "Data";
    
    Data .Yields[Muon][cut] = S[DoubleMu]      .Yields[Muon][cut];
    Data .Yields[Elec][cut] = S[DoubleElectron].Yields[Elec][cut];
    Data .Yields[ElMu][cut] = S[MuEG]          .Yields[ElMu][cut];
    
    Data .Yields_stat[Muon][cut] = S[DoubleMu]      .Yields_stat[Muon][cut];
    Data .Yields_stat[Elec][cut] = S[DoubleElectron].Yields_stat[Elec][cut];
    Data .Yields_stat[ElMu][cut] = S[MuEG]          .Yields_stat[ElMu][cut];
    
    Data .SSYields[Muon][cut] = S[DoubleMu]      .SSYields[Muon][cut];
    Data .SSYields[Elec][cut] = S[DoubleElectron].SSYields[Elec][cut];
    Data .SSYields[ElMu][cut] = S[MuEG]          .SSYields[ElMu][cut];
   
    Data .SSYields_stat[Muon][cut] = S[DoubleMu]      .SSYields_stat[Muon][cut];
    Data .SSYields_stat[Elec][cut] = S[DoubleElectron].SSYields_stat[Elec][cut];
    Data .SSYields_stat[ElMu][cut] = S[MuEG]          .SSYields_stat[ElMu][cut];

    for (size_t chan=0; chan<gNCHANNELS; chan++){
      TTbar.name = "ttbar";
      if (gUseTTMadSpin) TTbar.Yields[chan][cut] = S[TTJets_MadSpin]        .Yields[chan][cut]; 
      else               TTbar.Yields[chan][cut] = S[TTJetsFullLeptMGtauola].Yields[chan][cut]; 
      
      STop.name = "stop";
      STop .Yields[chan][cut]  = S[TbarWDilep].Yields[chan][cut];
      STop .Yields[chan][cut] += S[TWDilep]   .Yields[chan][cut];
      
      DY.name = "dy";
      DY   .Yields[chan][cut]  = S[ZJets_Madgraph].Yields[chan][cut];
      DY   .Yields[chan][cut] += S[DYJets_Madgraph].Yields[chan][cut];
      
      VV.name = "vv";
      VV   .Yields[chan][cut]  = S[WZ]                .Yields[chan][cut];
      VV   .Yields[chan][cut] += S[ZZ]                .Yields[chan][cut];
      VV   .Yields[chan][cut] += S[WWTo2L2Nu_Madgraph].Yields[chan][cut];
      
      Rare.name = "rare";
      Rare .Yields[chan][cut]  = S[TTWJets] .Yields[chan][cut];
      Rare .Yields[chan][cut] += S[TTZJets] .Yields[chan][cut];
      Rare .Yields[chan][cut] += S[TTGJets] .Yields[chan][cut];
      Rare .Yields[chan][cut] += S[TTWWJets].Yields[chan][cut];
      Rare .Yields[chan][cut] += S[WWWJets] .Yields[chan][cut];
      Rare .Yields[chan][cut] += S[WWZJets] .Yields[chan][cut];
      Rare .Yields[chan][cut] += S[WZZJets] .Yields[chan][cut];
      Rare .Yields[chan][cut] += S[ZZZJets] .Yields[chan][cut];

      Fake.name = "fake";
      Fake .Yields[chan][cut]  = S[TTJetsSemiLeptMGtauola].Yields[chan][cut];
      Fake .Yields[chan][cut] += S[Wbb_Madgraph]          .Yields[chan][cut];
      Fake .Yields[chan][cut] += S[WgammaToLNuG]          .Yields[chan][cut];

      Total.name = "total";
      Total.Yields[chan][cut]  = STop.Yields[chan][cut];
      Total.Yields[chan][cut] += DY  .Yields[chan][cut];
      Total.Yields[chan][cut] += VV  .Yields[chan][cut];
      Total.Yields[chan][cut] += Rare.Yields[chan][cut];
      Total.Yields[chan][cut] += Fake.Yields[chan][cut];

      // SYSTEMATIC ERROR
      for (size_t sys=0; sys<gNSYST; sys++){
	if (gUseTTMadSpin) {
	  TTbar.Yields_syst[chan][cut][sys]=S[TTJets_MadSpin].        Yields_syst[chan][cut][sys];
	}
	else               {
	  TTbar.Yields_syst[chan][cut][sys]=S[TTJetsFullLeptMGtauola].Yields_syst[chan][cut][sys];
	}
	
	STop .Yields_syst[chan][cut][sys]  = S[TbarWDilep].Yields_syst[chan][cut][sys];
	STop .Yields_syst[chan][cut][sys] += S[TWDilep]   .Yields_syst[chan][cut][sys];
	
	DY   .Yields_syst[chan][cut][sys]  = S[ZJets_Madgraph] .Yields_syst[chan][cut][sys];
	DY   .Yields_syst[chan][cut][sys] += S[DYJets_Madgraph].Yields_syst[chan][cut][sys];
	DY   .Yields_syst[chan][cut][sys]  = S[ZJets_Madgraph] .Yields_syst[chan][cut][sys];
	DY   .Yields_syst[chan][cut][sys] += S[DYJets_Madgraph].Yields_syst[chan][cut][sys];

	
	VV   .Yields_syst[chan][cut][sys]  = S[WZ]                .Yields_syst[chan][cut][sys];
	VV   .Yields_syst[chan][cut][sys] += S[ZZ]                .Yields_syst[chan][cut][sys];
	VV   .Yields_syst[chan][cut][sys] += S[WWTo2L2Nu_Madgraph].Yields_syst[chan][cut][sys];

	Rare .Yields_syst[chan][cut][sys]  = S[TTWJets] .Yields_syst[chan][cut][sys];
	Rare .Yields_syst[chan][cut][sys] += S[TTZJets] .Yields_syst[chan][cut][sys];
	Rare .Yields_syst[chan][cut][sys] += S[TTGJets] .Yields_syst[chan][cut][sys];
	Rare .Yields_syst[chan][cut][sys] += S[TTWWJets].Yields_syst[chan][cut][sys];
	Rare .Yields_syst[chan][cut][sys] += S[WWWJets] .Yields_syst[chan][cut][sys];
	Rare .Yields_syst[chan][cut][sys] += S[WWZJets] .Yields_syst[chan][cut][sys];
	Rare .Yields_syst[chan][cut][sys] += S[WZZJets] .Yields_syst[chan ][cut][sys];
	Rare .Yields_syst[chan][cut][sys] += S[ZZZJets] .Yields_syst[chan][cut][sys];

	Fake .Yields_syst[chan][cut][sys]  = S[TTJetsSemiLeptMGtauola].Yields_syst[chan][cut][sys];
	Fake .Yields_syst[chan][cut][sys] += S[Wbb_Madgraph]          .Yields_syst[chan][cut][sys];
	Fake .Yields_syst[chan][cut][sys] += S[WgammaToLNuG]          .Yields_syst[chan][cut][sys];

      }
      
      // STATISTICAL ERROR
      if (gUseTTMadSpin) TTbar.Yields_stat[chan][cut]  = S[TTJets_MadSpin].        Yields_stat[chan][cut] * S[TTJets_MadSpin].        Yields_stat[chan][cut];
      else               TTbar.Yields_stat[chan][cut]  = S[TTJetsFullLeptMGtauola].Yields_stat[chan][cut] * S[TTJetsFullLeptMGtauola].Yields_stat[chan][cut];
      
      STop .Yields_stat[chan][cut]  = S[TbarWDilep].Yields_stat[chan][cut] * S[TbarWDilep].Yields_stat[chan][cut];
      STop .Yields_stat[chan][cut] += S[TWDilep]   .Yields_stat[chan][cut] * S[TWDilep]   .Yields_stat[chan][cut];
      
      DY   .Yields_stat[chan][cut]  = S[ZJets_Madgraph] .Yields_stat[chan][cut] * S[ZJets_Madgraph] .Yields_stat[chan][cut];
      DY   .Yields_stat[chan][cut] += S[DYJets_Madgraph].Yields_stat[chan][cut] * S[DYJets_Madgraph].Yields_stat[chan][cut];
      
      VV   .Yields_stat[chan][cut]  = S[WZ]                .Yields_stat[chan][cut] * S[WZ]                .Yields_stat[chan][cut];
      VV   .Yields_stat[chan][cut] += S[ZZ]                .Yields_stat[chan][cut] * S[ZZ]                .Yields_stat[chan][cut];
      VV   .Yields_stat[chan][cut] += S[WWTo2L2Nu_Madgraph].Yields_stat[chan][cut] * S[WWTo2L2Nu_Madgraph].Yields_stat[chan][cut];
	
      Rare .Yields_stat[chan][cut]  = S[TTWJets] .Yields_stat[chan][cut] * S[TTWJets] .Yields_stat[chan][cut];
      Rare .Yields_stat[chan][cut] += S[TTZJets] .Yields_stat[chan][cut] * S[TTZJets] .Yields_stat[chan][cut];
      Rare .Yields_stat[chan][cut] += S[TTGJets] .Yields_stat[chan][cut] * S[TTGJets] .Yields_stat[chan][cut];
      Rare .Yields_stat[chan][cut] += S[TTWWJets].Yields_stat[chan][cut] * S[TTWWJets].Yields_stat[chan][cut];
      Rare .Yields_stat[chan][cut] += S[WWWJets] .Yields_stat[chan][cut] * S[WWWJets] .Yields_stat[chan][cut];
      Rare .Yields_stat[chan][cut] += S[WWZJets] .Yields_stat[chan][cut] * S[WWZJets] .Yields_stat[chan][cut];
      Rare .Yields_stat[chan][cut] += S[WZZJets] .Yields_stat[chan][cut] * S[WZZJets] .Yields_stat[chan][cut];
      Rare .Yields_stat[chan][cut] += S[ZZZJets] .Yields_stat[chan][cut] * S[ZZZJets] .Yields_stat[chan][cut];
	
      Fake .Yields_stat[chan][cut]  = S[TTJetsSemiLeptMGtauola].Yields_stat[chan][cut] * S[TTJetsSemiLeptMGtauola].Yields_stat[chan][cut];
      Fake .Yields_stat[chan][cut] += S[Wbb_Madgraph]          .Yields_stat[chan][cut] * S[Wbb_Madgraph]          .Yields_stat[chan][cut];
      Fake .Yields_stat[chan][cut] += S[WgammaToLNuG]          .Yields_stat[chan][cut] * S[WgammaToLNuG]          .Yields_stat[chan][cut];
      
      TTbar.Yields_stat[chan][cut] = TMath::Sqrt(TTbar.Yields_stat[chan][cut]);
      STop .Yields_stat[chan][cut] = TMath::Sqrt(STop .Yields_stat[chan][cut]);
      DY   .Yields_stat[chan][cut] = TMath::Sqrt(DY   .Yields_stat[chan][cut]);
      VV   .Yields_stat[chan][cut] = TMath::Sqrt(VV   .Yields_stat[chan][cut]);
      Rare .Yields_stat[chan][cut] = TMath::Sqrt(Rare .Yields_stat[chan][cut]);
      Fake .Yields_stat[chan][cut] = TMath::Sqrt(Fake .Yields_stat[chan][cut]);
      
      Total.Yields_stat[chan][cut]  = Fake .Yields_stat[chan][cut] * Fake .Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut] += Rare .Yields_stat[chan][cut] * Rare .Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut] += VV   .Yields_stat[chan][cut] * VV   .Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut] += DY   .Yields_stat[chan][cut] * DY   .Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut] += STop .Yields_stat[chan][cut] * STop .Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut] += TTbar.Yields_stat[chan][cut] * TTbar.Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut]  = TMath::Sqrt(Total.Yields_stat[chan][cut]);
      
      // SS Yields
      if (gUseTTMadSpin) TTbar.SSYields[chan][cut]  = S[TTJets_MadSpin].        SSYields[chan][cut];
      else               TTbar.SSYields[chan][cut]  = S[TTJetsFullLeptMGtauola].SSYields[chan][cut];
      
      STop .SSYields[chan][cut]  = S[TbarWDilep].SSYields[chan][cut];
      STop .SSYields[chan][cut] += S[TWDilep]   .SSYields[chan][cut];
            
      DY   .SSYields[chan][cut]  = S[ZJets_Madgraph].SSYields[chan][cut];
      DY   .SSYields[chan][cut] += S[DYJets_Madgraph].SSYields[chan][cut];
      
      VV   .SSYields[chan][cut]  = S[WZ]                .SSYields[chan][cut];
      VV   .SSYields[chan][cut] += S[ZZ]                .SSYields[chan][cut];
      VV   .SSYields[chan][cut] += S[WWTo2L2Nu_Madgraph].SSYields[chan][cut];
      
      Rare .SSYields[chan][cut]  = S[TTWJets] .SSYields[chan][cut];
      Rare .SSYields[chan][cut] += S[TTZJets] .SSYields[chan][cut];
      Rare .SSYields[chan][cut] += S[TTGJets] .SSYields[chan][cut];
      Rare .SSYields[chan][cut] += S[TTWWJets].SSYields[chan][cut];
      Rare .SSYields[chan][cut] += S[WWWJets] .SSYields[chan][cut];
      Rare .SSYields[chan][cut] += S[WWZJets] .SSYields[chan][cut];
      Rare .SSYields[chan][cut] += S[WZZJets] .SSYields[chan][cut];
      Rare .SSYields[chan][cut] += S[ZZZJets] .SSYields[chan][cut];

      Fake .SSYields[chan][cut]  = S[TTJetsSemiLeptMGtauola].SSYields[chan][cut];
      Fake .SSYields[chan][cut] += S[Wbb_Madgraph]          .SSYields[chan][cut];
      Fake .SSYields[chan][cut] += S[WgammaToLNuG]          .SSYields[chan][cut];

      Total.SSYields[chan][cut]  = STop .SSYields[chan][cut];
      Total.SSYields[chan][cut] += DY   .SSYields[chan][cut];
      Total.SSYields[chan][cut] += VV   .SSYields[chan][cut];
      Total.SSYields[chan][cut] += Rare .SSYields[chan][cut];
      Total.SSYields[chan][cut] += TTbar.SSYields[chan][cut];

      if (gUseTTMadSpin) TTbar.SSYields_stat[chan][cut]  = S[TTJets_MadSpin].        SSYields_stat[chan][cut] * S[TTJets_MadSpin].        SSYields_stat[chan][cut];
      else               TTbar.SSYields_stat[chan][cut]  = S[TTJetsFullLeptMGtauola].SSYields_stat[chan][cut] * S[TTJetsFullLeptMGtauola].SSYields_stat[chan][cut];
      
      STop .SSYields_stat[chan][cut]  = S[TbarWDilep].SSYields_stat[chan][cut] * S[TbarWDilep].SSYields_stat[chan][cut];
      STop .SSYields_stat[chan][cut] += S[TWDilep]   .SSYields_stat[chan][cut] * S[TWDilep]   .SSYields_stat[chan][cut];
      
      DY   .SSYields_stat[chan][cut]  = S[ZJets_Madgraph] .SSYields_stat[chan][cut] * S[ZJets_Madgraph] .SSYields_stat[chan][cut];
      DY   .SSYields_stat[chan][cut] += S[DYJets_Madgraph].SSYields_stat[chan][cut] * S[DYJets_Madgraph].SSYields_stat[chan][cut];
      
      VV   .SSYields_stat[chan][cut]  = S[WZ]                .SSYields_stat[chan][cut] * S[WZ]                .SSYields_stat[chan][cut];
      VV   .SSYields_stat[chan][cut] += S[ZZ]                .SSYields_stat[chan][cut] * S[ZZ]                .SSYields_stat[chan][cut];
      VV   .SSYields_stat[chan][cut] += S[WWTo2L2Nu_Madgraph].SSYields_stat[chan][cut] * S[WWTo2L2Nu_Madgraph].SSYields_stat[chan][cut];
	
      Rare .SSYields_stat[chan][cut]  = S[TTWJets] .SSYields_stat[chan][cut] * S[TTWJets] .SSYields_stat[chan][cut];
      Rare .SSYields_stat[chan][cut] += S[TTZJets] .SSYields_stat[chan][cut] * S[TTZJets] .SSYields_stat[chan][cut];
      Rare .SSYields_stat[chan][cut] += S[TTGJets] .SSYields_stat[chan][cut] * S[TTGJets] .SSYields_stat[chan][cut];
      Rare .SSYields_stat[chan][cut] += S[TTWWJets].SSYields_stat[chan][cut] * S[TTWWJets].SSYields_stat[chan][cut];
      Rare .SSYields_stat[chan][cut] += S[WWWJets] .SSYields_stat[chan][cut] * S[WWWJets] .SSYields_stat[chan][cut];
      Rare .SSYields_stat[chan][cut] += S[WWZJets] .SSYields_stat[chan][cut] * S[WWZJets] .SSYields_stat[chan][cut];
      Rare .SSYields_stat[chan][cut] += S[WZZJets] .SSYields_stat[chan][cut] * S[WZZJets] .SSYields_stat[chan][cut];
      Rare .SSYields_stat[chan][cut] += S[ZZZJets] .SSYields_stat[chan][cut] * S[ZZZJets] .SSYields_stat[chan][cut];
	
      Fake .SSYields_stat[chan][cut]  = S[TTJetsSemiLeptMGtauola].SSYields_stat[chan][cut] * S[TTJetsSemiLeptMGtauola].SSYields_stat[chan][cut];
      Fake .SSYields_stat[chan][cut] += S[Wbb_Madgraph]          .SSYields_stat[chan][cut] * S[Wbb_Madgraph]          .SSYields_stat[chan][cut];
      Fake .SSYields_stat[chan][cut] += S[WgammaToLNuG]          .SSYields_stat[chan][cut] * S[WgammaToLNuG]          .SSYields_stat[chan][cut];
      
      TTbar.SSYields_stat[chan][cut] = TMath::Sqrt(TTbar.SSYields_stat[chan][cut]);
      STop .SSYields_stat[chan][cut] = TMath::Sqrt(STop .SSYields_stat[chan][cut]);
      DY   .SSYields_stat[chan][cut] = TMath::Sqrt(DY   .SSYields_stat[chan][cut]);
      VV   .SSYields_stat[chan][cut] = TMath::Sqrt(VV   .SSYields_stat[chan][cut]);
      Rare .SSYields_stat[chan][cut] = TMath::Sqrt(Rare .SSYields_stat[chan][cut]);
      Fake .SSYields_stat[chan][cut] = TMath::Sqrt(Fake .SSYields_stat[chan][cut]);

      Total.SSYields_stat[chan][cut]  = Rare .SSYields_stat[chan][cut] * Rare .SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut] += VV   .SSYields_stat[chan][cut] * VV   .SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut] += DY   .SSYields_stat[chan][cut] * DY   .SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut] += STop .SSYields_stat[chan][cut] * STop .SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut] += TTbar.SSYields_stat[chan][cut] * TTbar.SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut]  = TMath::Sqrt(Total.SSYields_stat[chan][cut]);
      
    }
  }

  /// FOR DY ESTIMATION:
  Data.MllHistos[Muon] = (TH1F*)S[DoubleMu]      .MllHistos[Muon]->Clone();
  Data.MllHistos[Elec] = (TH1F*)S[DoubleElectron].MllHistos[Elec]->Clone();
  Data.MllHistos[ElMu] = (TH1F*)S[MuEG]          .MllHistos[ElMu]->Clone();
  for (size_t chan=0; chan<gNCHANNELS; chan++){
    DY.MllHistos[chan] = (TH1F*)S[DYJets_Madgraph].MllHistos[chan]->Clone();
    DY.MllHistos[chan] ->Add(   S[ZJets_Madgraph] .MllHistos[chan]);
  }

  for (size_t chan=0; chan<gNCHANNELS+1; chan++){
    for (size_t cut=0; cut<iNCUTS; cut++){
      for (size_t var=0; var<gNVARS; var++){
	if (chan == gNCHANNELS){
	  Data .KinHistos[chan][cut][var] = (TH1F*)S[DoubleMu]              .KinHistos[Muon][cut][var]->Clone();
	  Data .KinHistos[chan][cut][var] ->Add(   S[DoubleElectron]        .KinHistos[Elec][cut][var]);
	  if (gUseTTMadSpin) { 
	    TTbar.KinHistos[chan][cut][var] = (TH1F*)S[TTJets_MadSpin]        .KinHistos[Muon][cut][var]->Clone();
	    TTbar.KinHistos[chan][cut][var] ->Add(   S[TTJets_MadSpin]        .KinHistos[Elec][cut][var]);
	  }
	  else {
	    TTbar.KinHistos[chan][cut][var] = (TH1F*)S[TTJetsFullLeptMGtauola].KinHistos[Muon][cut][var]->Clone();
	    TTbar.KinHistos[chan][cut][var] ->Add(   S[TTJetsFullLeptMGtauola].KinHistos[Elec][cut][var]);
	  }
	  
	  STop .KinHistos[chan][cut][var] = (TH1F*)S[TbarWDilep]            .KinHistos[Muon][cut][var]->Clone();
	  STop .KinHistos[chan][cut][var] ->Add(   S[TbarWDilep]            .KinHistos[Elec][cut][var]);
	  STop .KinHistos[chan][cut][var] ->Add(   S[TWDilep]               .KinHistos[Muon][cut][var]);
	  STop .KinHistos[chan][cut][var] ->Add(   S[TWDilep]               .KinHistos[Elec][cut][var]);
	  DY   .KinHistos[chan][cut][var] = (TH1F*)S[DYJets_Madgraph]       .KinHistos[Muon][cut][var]->Clone();
	  DY   .KinHistos[chan][cut][var] ->Add(   S[ZJets_Madgraph]        .KinHistos[Muon][cut][var]);
	  DY   .KinHistos[chan][cut][var] ->Add(   S[DYJets_Madgraph]       .KinHistos[Elec][cut][var]);
	  DY   .KinHistos[chan][cut][var] ->Add(   S[ZJets_Madgraph]        .KinHistos[Elec][cut][var]);
	  VV   .KinHistos[chan][cut][var] = (TH1F*)S[WWTo2L2Nu_Madgraph]    .KinHistos[Muon][cut][var]->Clone();
	  VV   .KinHistos[chan][cut][var] ->Add(   S[WZ]                    .KinHistos[Muon][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[ZZ]                    .KinHistos[Muon][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[WWTo2L2Nu_Madgraph]    .KinHistos[Elec][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[WZ]                    .KinHistos[Elec][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[ZZ]                    .KinHistos[Elec][cut][var]);
	  Fake .KinHistos[chan][cut][var] = (TH1F*)S[TTJetsSemiLeptMGtauola].KinHistos[Muon][cut][var]->Clone();  
	  Fake .KinHistos[chan][cut][var] ->Add(   S[WgammaToLNuG]          .KinHistos[Muon][cut][var]);
	  Fake .KinHistos[chan][cut][var] ->Add(   S[Wbb_Madgraph]          .KinHistos[Muon][cut][var]);
	  Fake .KinHistos[chan][cut][var] ->Add(   S[TTJetsSemiLeptMGtauola].KinHistos[Elec][cut][var]);
	  Fake .KinHistos[chan][cut][var] ->Add(   S[WgammaToLNuG]          .KinHistos[Elec][cut][var]);
	  Fake .KinHistos[chan][cut][var] ->Add(   S[Wbb_Madgraph]          .KinHistos[Elec][cut][var]);
	  Rare .KinHistos[chan][cut][var] = (TH1F*)S[TTGJets]               .KinHistos[Muon][cut][var]->Clone();  
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTWJets]               .KinHistos[Muon][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTWWJets]              .KinHistos[Muon][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTZJets]               .KinHistos[Muon][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[WWWJets]               .KinHistos[Muon][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[WWZJets]               .KinHistos[Muon][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[WZZJets]               .KinHistos[Muon][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[ZZZJets]               .KinHistos[Muon][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTGJets]               .KinHistos[Elec][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTWJets]               .KinHistos[Elec][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTWWJets]              .KinHistos[Elec][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTZJets]               .KinHistos[Elec][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[WWWJets]               .KinHistos[Elec][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[WWZJets]               .KinHistos[Elec][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[WZZJets]               .KinHistos[Elec][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[ZZZJets]               .KinHistos[Elec][cut][var]);
	}
	else {
	  if (chan == Muon){
	    Data .KinHistos[chan][cut][var] = (TH1F*)S[DoubleMu]            .KinHistos[chan][cut][var]->Clone();
	  }
	  else if (chan == Elec){
	    Data .KinHistos[chan][cut][var] = (TH1F*)S[DoubleElectron]      .KinHistos[chan][cut][var]->Clone();
	  }
	  else if (chan == ElMu){
	    Data .KinHistos[chan][cut][var] = (TH1F*)S[MuEG]                .KinHistos[chan][cut][var]->Clone();
	  }
	  if (gUseTTMadSpin) 
	    TTbar.KinHistos[chan][cut][var] = (TH1F*)S[TTJets_MadSpin]        .KinHistos[chan][cut][var]->Clone();
	  else 
	    TTbar.KinHistos[chan][cut][var] = (TH1F*)S[TTJetsFullLeptMGtauola].KinHistos[chan][cut][var]->Clone();
	  
	  STop .KinHistos[chan][cut][var] = (TH1F*)S[TbarWDilep]            .KinHistos[chan][cut][var]->Clone();
	  STop .KinHistos[chan][cut][var] ->Add(   S[TWDilep]               .KinHistos[chan][cut][var]);
	  DY   .KinHistos[chan][cut][var] = (TH1F*)S[DYJets_Madgraph]       .KinHistos[chan][cut][var]->Clone();
	  DY   .KinHistos[chan][cut][var] ->Add(   S[ZJets_Madgraph]        .KinHistos[chan][cut][var]);
	  VV   .KinHistos[chan][cut][var] = (TH1F*)S[WWTo2L2Nu_Madgraph]    .KinHistos[chan][cut][var]->Clone();
	  VV   .KinHistos[chan][cut][var] ->Add(   S[WZ]                    .KinHistos[chan][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[ZZ]                    .KinHistos[chan][cut][var]);
	  Fake .KinHistos[chan][cut][var] = (TH1F*)S[TTJetsSemiLeptMGtauola].KinHistos[chan][cut][var]->Clone();  
	  Fake .KinHistos[chan][cut][var] ->Add(   S[WgammaToLNuG]          .KinHistos[chan][cut][var]);
	  Fake .KinHistos[chan][cut][var] ->Add(   S[Wbb_Madgraph]          .KinHistos[chan][cut][var]);
	  Rare .KinHistos[chan][cut][var] = (TH1F*)S[TTGJets]               .KinHistos[chan][cut][var]->Clone();  
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTWJets]               .KinHistos[chan][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTWWJets]              .KinHistos[chan][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTZJets]               .KinHistos[chan][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[WWWJets]               .KinHistos[chan][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[WWZJets]               .KinHistos[chan][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[WZZJets]               .KinHistos[chan][cut][var]);
	  Rare .KinHistos[chan][cut][var] ->Add(   S[ZZZJets]               .KinHistos[chan][cut][var]);

	  
	}
	SetupDraw(Data .KinHistos[chan][cut][var], kBlack   , var); 
	SetupDraw(TTbar.KinHistos[chan][cut][var], kRed+1   , var);  
	SetupDraw(STop .KinHistos[chan][cut][var], kPink-3  , var);  
	SetupDraw(DY   .KinHistos[chan][cut][var], kAzure-2 , var);  
	SetupDraw(VV   .KinHistos[chan][cut][var], kOrange-3, var);  
	SetupDraw(Rare .KinHistos[chan][cut][var], kYellow  , var);  
	SetupDraw(Fake .KinHistos[chan][cut][var], kGreen-3 , var);  
      }
    }
  }

  // SYSTEMATIC ERRORS HISTOS
  for (size_t cut=0; cut<iNCUTS; cut++){
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      for (size_t sys=0; sys<gNSYST; sys++){
	if (chan == Muon){
	  Data .SysHistos[chan][cut][sys]   = (TH1F*)S[DoubleMu]        .SysHistos[chan][cut][0]->Clone();
	  Data .SSSysHistos[chan][cut][sys] = (TH1F*)S[DoubleMu]      .SSSysHistos[chan][cut][0]->Clone();
	}
	else if (chan == Elec){
	  Data .SysHistos[chan][cut][sys]   = (TH1F*)S[DoubleElectron]  .SysHistos[chan][cut][0]->Clone();
	  Data .SSSysHistos[chan][cut][sys] = (TH1F*)S[DoubleElectron].SSSysHistos[chan][cut][0]->Clone();
	}
	else if (chan == ElMu){
	  Data .SysHistos[chan][cut][sys]   = (TH1F*)S[MuEG]            .SysHistos[chan][cut][0]->Clone();
	  Data .SSSysHistos[chan][cut][sys] = (TH1F*)S[MuEG]          .SSSysHistos[chan][cut][0]->Clone();
	}
	
	if (gUseTTMadSpin) {
	  TTbar.SysHistos[chan][cut][sys]   = (TH1F*)S[TTJets_MadSpin]        .SysHistos[chan][cut][sys]->Clone();
	  TTbar.SSSysHistos[chan][cut][sys] = (TH1F*)S[TTJets_MadSpin]      .SSSysHistos[chan][cut][sys]->Clone();
	}
	else {
	  TTbar.SysHistos[chan][cut][sys]   = (TH1F*)S[TTJetsFullLeptMGtauola].SysHistos[chan][cut][sys]->Clone();
	  TTbar.SSSysHistos[chan][cut][sys] = (TH1F*)S[TTJetsFullLeptMGtauola].SSSysHistos[chan][cut][sys]->Clone();
	}
	STop .SysHistos[chan][cut][sys] = (TH1F*)S[TbarWDilep]            .SysHistos[chan][cut][sys]->Clone();
	STop .SysHistos[chan][cut][sys] ->Add(   S[TWDilep]               .SysHistos[chan][cut][sys]);
	DY   .SysHistos[chan][cut][sys] = (TH1F*)S[DYJets_Madgraph]       .SysHistos[chan][cut][sys]->Clone();
	DY   .SysHistos[chan][cut][sys] ->Add(   S[ZJets_Madgraph]        .SysHistos[chan][cut][sys]);
	VV   .SysHistos[chan][cut][sys] = (TH1F*)S[WWTo2L2Nu_Madgraph]    .SysHistos[chan][cut][sys]->Clone();
	VV   .SysHistos[chan][cut][sys] ->Add(   S[WZ]                    .SysHistos[chan][cut][sys]);
	VV   .SysHistos[chan][cut][sys] ->Add(   S[ZZ]                    .SysHistos[chan][cut][sys]);
	Fake .SysHistos[chan][cut][sys] = (TH1F*)S[TTJetsSemiLeptMGtauola].SysHistos[chan][cut][sys]->Clone();  
	Fake .SysHistos[chan][cut][sys] ->Add(   S[WgammaToLNuG]          .SysHistos[chan][cut][sys]);
	Fake .SysHistos[chan][cut][sys] ->Add(   S[Wbb_Madgraph]          .SysHistos[chan][cut][sys]);
	Rare .SysHistos[chan][cut][sys] = (TH1F*)S[TTGJets]               .SysHistos[chan][cut][sys]->Clone();  
	Rare .SysHistos[chan][cut][sys] ->Add(   S[TTWJets]               .SysHistos[chan][cut][sys]);
	Rare .SysHistos[chan][cut][sys] ->Add(   S[TTWWJets]              .SysHistos[chan][cut][sys]);
	Rare .SysHistos[chan][cut][sys] ->Add(   S[TTZJets]               .SysHistos[chan][cut][sys]);
	Rare .SysHistos[chan][cut][sys] ->Add(   S[WWWJets]               .SysHistos[chan][cut][sys]);
	Rare .SysHistos[chan][cut][sys] ->Add(   S[WWZJets]               .SysHistos[chan][cut][sys]);
	Rare .SysHistos[chan][cut][sys] ->Add(   S[WZZJets]               .SysHistos[chan][cut][sys]);
	Rare .SysHistos[chan][cut][sys] ->Add(   S[ZZZJets]               .SysHistos[chan][cut][sys]);
	
	SetupDraw(Data .SysHistos[chan][cut][sys], kBlack   , NBTagsNJets); 
	SetupDraw(TTbar.SysHistos[chan][cut][sys], kRed+1   , NBTagsNJets);  
	SetupDraw(STop .SysHistos[chan][cut][sys], kPink-3  , NBTagsNJets);  
	SetupDraw(DY   .SysHistos[chan][cut][sys], kAzure-2 , NBTagsNJets);  
	SetupDraw(VV   .SysHistos[chan][cut][sys], kOrange-3, NBTagsNJets);  
	SetupDraw(Rare .SysHistos[chan][cut][sys], kYellow  , NBTagsNJets);  
	SetupDraw(Fake .SysHistos[chan][cut][sys], kGreen-3 , NBTagsNJets);  

	STop .SSSysHistos[chan][cut][sys] = (TH1F*)S[TbarWDilep]            .SSSysHistos[chan][cut][sys]->Clone();
	STop .SSSysHistos[chan][cut][sys] ->Add(   S[TWDilep]               .SSSysHistos[chan][cut][sys]);
	DY   .SSSysHistos[chan][cut][sys] = (TH1F*)S[DYJets_Madgraph]       .SSSysHistos[chan][cut][sys]->Clone();
	DY   .SSSysHistos[chan][cut][sys] ->Add(   S[ZJets_Madgraph]        .SSSysHistos[chan][cut][sys]);
	VV   .SSSysHistos[chan][cut][sys] = (TH1F*)S[WWTo2L2Nu_Madgraph]    .SSSysHistos[chan][cut][sys]->Clone();
	VV   .SSSysHistos[chan][cut][sys] ->Add(   S[WZ]                    .SSSysHistos[chan][cut][sys]);
	VV   .SSSysHistos[chan][cut][sys] ->Add(   S[ZZ]                    .SSSysHistos[chan][cut][sys]);
	Fake .SSSysHistos[chan][cut][sys] = (TH1F*)S[TTJetsSemiLeptMGtauola].SSSysHistos[chan][cut][sys]->Clone();  
	Fake .SSSysHistos[chan][cut][sys] ->Add(   S[WgammaToLNuG]          .SSSysHistos[chan][cut][sys]);
	Fake .SSSysHistos[chan][cut][sys] ->Add(   S[Wbb_Madgraph]          .SSSysHistos[chan][cut][sys]);
	Rare .SSSysHistos[chan][cut][sys] = (TH1F*)S[TTGJets]               .SSSysHistos[chan][cut][sys]->Clone();  
	Rare .SSSysHistos[chan][cut][sys] ->Add(   S[TTWJets]               .SSSysHistos[chan][cut][sys]);
	Rare .SSSysHistos[chan][cut][sys] ->Add(   S[TTWWJets]              .SSSysHistos[chan][cut][sys]);
	Rare .SSSysHistos[chan][cut][sys] ->Add(   S[TTZJets]               .SSSysHistos[chan][cut][sys]);
	Rare .SSSysHistos[chan][cut][sys] ->Add(   S[WWWJets]               .SSSysHistos[chan][cut][sys]);
	Rare .SSSysHistos[chan][cut][sys] ->Add(   S[WWZJets]               .SSSysHistos[chan][cut][sys]);
	Rare .SSSysHistos[chan][cut][sys] ->Add(   S[WZZJets]               .SSSysHistos[chan][cut][sys]);
	Rare .SSSysHistos[chan][cut][sys] ->Add(   S[ZZZJets]               .SSSysHistos[chan][cut][sys]);
	
	SetupDraw(Data .SSSysHistos[chan][cut][sys], kBlack   , NBTagsNJets); 
	SetupDraw(TTbar.SSSysHistos[chan][cut][sys], kRed+1   , NBTagsNJets);  
	SetupDraw(STop .SSSysHistos[chan][cut][sys], kPink-3  , NBTagsNJets);  
	SetupDraw(DY   .SSSysHistos[chan][cut][sys], kAzure-2 , NBTagsNJets);  
	SetupDraw(VV   .SSSysHistos[chan][cut][sys], kOrange-3, NBTagsNJets);  
	SetupDraw(Rare .SSSysHistos[chan][cut][sys], kYellow  , NBTagsNJets);  
	SetupDraw(Fake .SSSysHistos[chan][cut][sys], kGreen-3 , NBTagsNJets);  
      }
      if (gUseTTMadSpin) {
	TTbar.SysHistos[chan][cut][Q2ScaleUp   ] = (TH1F*)S[TTJets_scaleup]     .SysHistos[chan][cut][0]->Clone();
	TTbar.SysHistos[chan][cut][Q2ScaleDown ] = (TH1F*)S[TTJets_scaledown]   .SysHistos[chan][cut][0]->Clone();
	TTbar.SysHistos[chan][cut][MatchingUp  ] = (TH1F*)S[TTJets_matchingup]  .SysHistos[chan][cut][0]->Clone();
	TTbar.SysHistos[chan][cut][MatchingDown] = (TH1F*)S[TTJets_matchingdown].SysHistos[chan][cut][0]->Clone();
	
	SetupDraw(TTbar.SysHistos[chan][cut][Q2ScaleUp   ], kRed+1   , NBTagsNJets);  
	SetupDraw(TTbar.SysHistos[chan][cut][Q2ScaleDown ], kRed+1   , NBTagsNJets);  
	SetupDraw(TTbar.SysHistos[chan][cut][MatchingUp  ], kRed+1   , NBTagsNJets);  
	SetupDraw(TTbar.SysHistos[chan][cut][MatchingDown], kRed+1   , NBTagsNJets);  	

	//// SS 
	TTbar.SSSysHistos[chan][cut][Q2ScaleUp   ]=(TH1F*)S[TTJets_scaleup]     .SSSysHistos[chan][cut][0]->Clone();
	TTbar.SSSysHistos[chan][cut][Q2ScaleDown ]=(TH1F*)S[TTJets_scaledown]   .SSSysHistos[chan][cut][0]->Clone();
	TTbar.SSSysHistos[chan][cut][MatchingUp  ]=(TH1F*)S[TTJets_matchingup]  .SSSysHistos[chan][cut][0]->Clone();
	TTbar.SSSysHistos[chan][cut][MatchingDown]=(TH1F*)S[TTJets_matchingdown].SSSysHistos[chan][cut][0]->Clone();
	
	SetupDraw(TTbar.SSSysHistos[chan][cut][Q2ScaleUp   ], kRed+1   , NBTagsNJets);  
	SetupDraw(TTbar.SSSysHistos[chan][cut][Q2ScaleDown ], kRed+1   , NBTagsNJets);  
	SetupDraw(TTbar.SSSysHistos[chan][cut][MatchingUp  ], kRed+1   , NBTagsNJets);  
	SetupDraw(TTbar.SSSysHistos[chan][cut][MatchingDown], kRed+1   , NBTagsNJets);  	
      }
    }
  }
}
void TopPlotter::ResetDataMembers(){
  
  // Resetting samples holder...
  for (size_t sample=0; sample<gNSAMPLES; sample++){
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      for (size_t cut=0; cut<iNCUTS; cut++){
	S[sample].Yields     [chan][cut] = 0.;
	S[sample].Yields_stat[chan][cut] = 0.;
	S[sample].SSYields   [chan][cut] = 0.;
	S[sample].SSYields_stat[chan][cut] = 0.;

	for (size_t sys=0; sys<gNSYST; sys++){
	  S[sample].Yields_syst[chan][cut][sys] = 0.;
	}
      }
      for (size_t sys=0; sys<gNSYSTERRTypes; sys++){
	S[sample].SystError[chan][sys] = 0.;
      }
    }
  }
  
  // Resetting background / signal categories...
  for (size_t chan=0; chan<gNCHANNELS; chan++){
    for (size_t cut=0; cut<iNCUTS; cut++){
      Data .  Yields[chan][cut] = 0.;    Data .  Yields_stat[chan][cut] = 0.;
      TTbar.  Yields[chan][cut] = 0.;	 TTbar.  Yields_stat[chan][cut] = 0.;
      STop .  Yields[chan][cut] = 0.;	 STop .  Yields_stat[chan][cut] = 0.;
      DY   .  Yields[chan][cut] = 0.;	 DY   .  Yields_stat[chan][cut] = 0.;
      VV   .  Yields[chan][cut] = 0.;	 VV   .  Yields_stat[chan][cut] = 0.;
      Rare .  Yields[chan][cut] = 0.;	 Rare .  Yields_stat[chan][cut] = 0.;
      Fake .  Yields[chan][cut] = 0.;	 Fake .  Yields_stat[chan][cut] = 0.;
      Total.  Yields[chan][cut] = 0.;	 Total.  Yields_stat[chan][cut] = 0.;

      for (size_t sys=0; sys<gNSYST; sys++){
	Data .  Yields_syst[chan][cut][sys] = 0.;    
	TTbar.  Yields_syst[chan][cut][sys] = 0.;    
	STop .  Yields_syst[chan][cut][sys] = 0.;    
	DY   .  Yields_syst[chan][cut][sys] = 0.;    
	VV   .  Yields_syst[chan][cut][sys] = 0.;    
	Rare .  Yields_syst[chan][cut][sys] = 0.;    
	Fake .  Yields_syst[chan][cut][sys] = 0.;    
	Total.  Yields_syst[chan][cut][sys] = 0.;    

	DD_DY  .Yields_syst[chan][cut][sys] = 0.;    
	DD_NonW.Yields_syst[chan][cut][sys] = 0.;    
	
      }

      Data .SSYields[chan][cut] = 0.; Data .SSYields_stat[chan][cut] = 0.;
      TTbar.SSYields[chan][cut] = 0.; TTbar.SSYields_stat[chan][cut] = 0.;
      STop .SSYields[chan][cut] = 0.; STop .SSYields_stat[chan][cut] = 0.;
      DY   .SSYields[chan][cut] = 0.; DY   .SSYields_stat[chan][cut] = 0.;
      VV   .SSYields[chan][cut] = 0.; VV   .SSYields_stat[chan][cut] = 0.;
      Rare .SSYields[chan][cut] = 0.; Rare .SSYields_stat[chan][cut] = 0.;
      Fake .SSYields[chan][cut] = 0.; Fake .SSYields_stat[chan][cut] = 0.;
      Total.SSYields[chan][cut] = 0.; Total.SSYields_stat[chan][cut] = 0.;

      // Data Driven
      DD_DY  .  Yields[chan][cut] = 0.;   DD_DY  .Yields_stat[chan][cut] = 0.;
      DD_NonW.  Yields[chan][cut] = 0.;	  DD_NonW.Yields_stat[chan][cut] = 0.;
      DD_DY  .SSYields[chan][cut] = 0.;
      DD_NonW.SSYields[chan][cut] = 0.;
    }
    for (size_t sys=0; sys<gNSYSTERRTypes; sys++){
      Data .SystError[chan][sys] = 0.;
      TTbar.SystError[chan][sys] = 0.;
      STop .SystError[chan][sys] = 0.;
      DY   .SystError[chan][sys] = 0.;
      VV   .SystError[chan][sys] = 0.;
      Rare .SystError[chan][sys] = 0.;
      Fake .SystError[chan][sys] = 0.;
      Total.SystError[chan][sys] = 0.;
    }
  }
  
  // Reset XSection...
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    ttbar.xsec     [ch] = 0.;
    ttbar.xsec_syst[ch] = 0.;
    ttbar.xsec_stat[ch] = 0.;
    ttbar.xsec_lumi[ch] = 0.;
    ttbar.acc      [ch] = 0.;
    ttbar.acc_stat [ch] = 0.;
    ttbar.acc_syst [ch] = 0.;

    ttbar.err_VV   [ch] = 0.;
    ttbar.err_STop [ch] = 0.;
    ttbar.err_Fake [ch] = 0.;
    ttbar.err_Rare [ch] = 0.;
    ttbar.err_IDIso[ch] = 0.;
    ttbar.err_Trig [ch] = 0.;
    ttbar.err_LES  [ch] = 0.;
    ttbar.err_JES  [ch] = 0.;
    ttbar.err_JER  [ch] = 0.;
    ttbar.err_Btag [ch] = 0.;
    ttbar.err_PU   [ch] = 0.;
    ttbar.err_Q2   [ch] = 0.;
    ttbar.err_Match[ch] = 0.;
  }
}
void TopPlotter::ResetSystematicErrors(){
  for (size_t chan=0; chan<gNCHANNELS; chan++){
    for (size_t sys=0; sys<gNSYSTERRTypes; sys++){
      //      Data .SystError[chan][sys] = 0.;
      TTbar.SystError[chan][sys] = 0.;
      STop .SystError[chan][sys] = 0.;
      DY   .SystError[chan][sys] = 0.;
      VV   .SystError[chan][sys] = 0.;
      Rare .SystError[chan][sys] = 0.;
      Fake .SystError[chan][sys] = 0.;
      Total.SystError[chan][sys] = 0.;
    }
    
    for (size_t cut=0; cut<iNCUTS; cut++){
      //      Data .  Yields_syst[chan][cut][0] = 0.;    
      TTbar.  Yields_syst[chan][cut][0] = 0.;    
      STop .  Yields_syst[chan][cut][0] = 0.;    
      DY   .  Yields_syst[chan][cut][0] = 0.;    
      VV   .  Yields_syst[chan][cut][0] = 0.;    
      Rare .  Yields_syst[chan][cut][0] = 0.;    
      Fake .  Yields_syst[chan][cut][0] = 0.;    
      Total.  Yields_syst[chan][cut][0] = 0.;    
    }
  }
  // Reset XSection...
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    ttbar.xsec     [ch] = 0.;
    ttbar.xsec_syst[ch] = 0.;
    ttbar.xsec_stat[ch] = 0.;
    ttbar.xsec_lumi[ch] = 0.;
    ttbar.acc      [ch] = 0.;
    ttbar.acc_stat [ch] = 0.;
    ttbar.acc_syst [ch] = 0.;

    ttbar.err_VV   [ch] = 0.;
    ttbar.err_STop [ch] = 0.;
    ttbar.err_Fake [ch] = 0.;
    ttbar.err_Rare [ch] = 0.;
    ttbar.err_IDIso[ch] = 0.;
    ttbar.err_Trig [ch] = 0.;
    ttbar.err_LES  [ch] = 0.;
    ttbar.err_JES  [ch] = 0.;
    ttbar.err_JER  [ch] = 0.;
    ttbar.err_Btag [ch] = 0.;
    ttbar.err_PU   [ch] = 0.;
    ttbar.err_Q2   [ch] = 0.;
    ttbar.err_Match[ch] = 0.;
  }

}
void TopPlotter::PrintSystematicErrors(){
  fOutputSubDir = "XSection/";
  
  TString filename = "";
  if (gUseTTMadSpin) filename = fOutputDir + fOutputSubDir + Form("SystErrors_MadSpin.txt");
  else               filename = fOutputDir + fOutputSubDir + Form("SystErrors_Madgraph.txt");
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  fOUTSTREAM2.open(filename, ios::trunc);

  fOUTSTREAM2 << "===========================================================================" << endl;
  fOUTSTREAM2 << "                      Systematic Uncertainties in pb (%)                   " << endl;
  fOUTSTREAM2 << " Source                |      El/El      |      Mu/Mu     |     El/Mu      " << endl;
  fOUTSTREAM2 << "---------------------------------------------------------------------------" << endl;
  fOUTSTREAM2 << Form(" VV                    | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_VV[Elec],  100 * ttbar.err_VV[Elec] / ttbar.xsec[Elec],
		      ttbar.err_VV[Muon],  100 * ttbar.err_VV[Muon] / ttbar.xsec[Muon],
		      ttbar.err_VV[ElMu],  100 * ttbar.err_VV[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" STop                  | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_STop[Elec],  100 * ttbar.err_STop[Elec] / ttbar.xsec[Elec],
		      ttbar.err_STop[Muon],  100 * ttbar.err_STop[Muon] / ttbar.xsec[Muon],
		      ttbar.err_STop[ElMu],  100 * ttbar.err_STop[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Fake                  | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Fake[Elec],  100 * ttbar.err_Fake[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Fake[Muon],  100 * ttbar.err_Fake[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Fake[ElMu],  100 * ttbar.err_Fake[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Rare                  | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Rare[Elec],  100 * ttbar.err_Rare[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Rare[Muon],  100 * ttbar.err_Rare[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Rare[ElMu],  100 * ttbar.err_Rare[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Lepton Efficiencies   | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_IDIso[Elec],  100 * ttbar.err_IDIso[Elec] / ttbar.xsec[Elec],
		      ttbar.err_IDIso[Muon],  100 * ttbar.err_IDIso[Muon] / ttbar.xsec[Muon],
		      ttbar.err_IDIso[ElMu],  100 * ttbar.err_IDIso[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Trigger Efficiencies  | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Trig[Elec],  100 * ttbar.err_Trig[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Trig[Muon],  100 * ttbar.err_Trig[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Trig[ElMu],  100 * ttbar.err_Trig[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Lepton Energy Scale   | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_LES[Elec],  100 * ttbar.err_LES[Elec] / ttbar.xsec[Elec],
		      ttbar.err_LES[Muon],  100 * ttbar.err_LES[Muon] / ttbar.xsec[Muon],
		      ttbar.err_LES[ElMu],  100 * ttbar.err_LES[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Jet Energy Scale      | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_JES[Elec],  100 * ttbar.err_JES[Elec] / ttbar.xsec[Elec],
		      ttbar.err_JES[Muon],  100 * ttbar.err_JES[Muon] / ttbar.xsec[Muon],
		      ttbar.err_JES[ElMu],  100 * ttbar.err_JES[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Jet Energy Resolution | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_JER[Elec],  100 * ttbar.err_JER[Elec] / ttbar.xsec[Elec],
		      ttbar.err_JER[Muon],  100 * ttbar.err_JER[Muon] / ttbar.xsec[Muon],
		      ttbar.err_JER[ElMu],  100 * ttbar.err_JER[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" b-tagging             | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Btag[Elec],  100 * ttbar.err_Btag[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Btag[Muon],  100 * ttbar.err_Btag[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Btag[ElMu],  100 * ttbar.err_Btag[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Pile Up               | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_PU[Elec],  100 * ttbar.err_PU[Elec] / ttbar.xsec[Elec],
		      ttbar.err_PU[Muon],  100 * ttbar.err_PU[Muon] / ttbar.xsec[Muon],
		      ttbar.err_PU[ElMu],  100 * ttbar.err_PU[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" QCD scale             | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Q2[Elec],  100 * ttbar.err_Q2[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Q2[Muon],  100 * ttbar.err_Q2[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Q2[ElMu],  100 * ttbar.err_Q2[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Matching partons      | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Match[Elec],  100 * ttbar.err_Match[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Match[Muon],  100 * ttbar.err_Match[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Match[ElMu],  100 * ttbar.err_Match[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2.close();
  
  gSystem->Exec("cat "+filename);
}
void TopPlotter::PrintYieldsWithDD(){
  if (gUseTTMadSpin) fOutputSubDir = "Yields_DD_Madspin/";
  else               fOutputSubDir = "Yields_DD_Madgraph/"; 

  TString yieldsfilename = "";
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  for (size_t cut=0; cut<iNCUTS; cut++){
    yieldsfilename = fOutputDir + fOutputSubDir + "Yields_MC_"+sCut[cut]+".txt";
    fOUTSTREAM.open(yieldsfilename, ios::trunc);

    fOUTSTREAM << "/////////////////////////////////////////////////////////////////////////////" << endl;
    fOUTSTREAM << " Producing MC yields for cut " << sCut[cut] << endl;
    fOUTSTREAM << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;
    
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << "                YIELDS  |    El/El   |    Mu/Mu   |    El/Mu   ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    
    float sum_mm(0.), sum_em(0.), sum_ee(0.);
    for (size_t sample=TTJets_MadSpin; sample<gNSAMPLES; sample++){
      if ( gUseTTMadSpin && sample == TTJetsFullLeptMGtauola) continue;
      if (!gUseTTMadSpin && sample == TTJets_MadSpin        ) continue;
      
      float temp_mm  = S[sample].Yields[Muon][cut]; sum_mm  += temp_mm ;
      float temp_em  = S[sample].Yields[ElMu][cut]; sum_em  += temp_em ;
      float temp_ee  = S[sample].Yields[Elec][cut]; sum_ee  += temp_ee ;
      
      fOUTSTREAM << Form("%23s |  %8.2f  |  %8.2f  |  %8.2f  ||\n", SampleName[sample].Data(),
			 temp_ee , temp_mm , temp_em );
    }
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << Form("%23s |  %8.2f  |  %8.2f  |  %8.2f  ||\n", "MC Sum", 
		       sum_ee , sum_mm , sum_em );
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << Form("%23s |  %5.0f     |  %5.0f     |  %5.0f     ||\n", "Data", 
		       S[DoubleElectron].Yields[Elec][cut], S[DoubleMu].Yields[Muon][cut], S[MuEG].Yields[ElMu][cut]); 
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    
    // Now POST a Summary...
    fOUTSTREAM << "=============================================================================" << endl;
    fOUTSTREAM << "                            SUMMARY (stat error only):                       " << endl; 
    fOUTSTREAM << "=============================================================================" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << "          Source  |       El/El      |       Mu/Mu      |       El/Mu      ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Drell-Yan        | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||", 
		       DD_DY.Yields[Elec][cut], DD_DY.Yields_stat[Elec][cut], 
		       DD_DY.Yields[Muon][cut], DD_DY.Yields_stat[Muon][cut], 
		       DD_DY.Yields[ElMu][cut], DD_DY.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Non W/Z leptons  | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       DD_NonW.Yields[Elec][cut], DD_NonW.Yields_stat[Elec][cut], 
		       DD_NonW.Yields[Muon][cut], DD_NonW.Yields_stat[Muon][cut], 
		       DD_NonW.Yields[ElMu][cut], DD_NonW.Yields_stat[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << Form(" Single top quark | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       STop.Yields[Elec][cut], STop.Yields_stat[Elec][cut],
		       STop.Yields[Muon][cut], STop.Yields_stat[Muon][cut], 
		       STop.Yields[ElMu][cut], STop.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Dibosons         | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       VV.Yields[Elec][cut], VV.Yields_stat[Elec][cut], 
		       VV.Yields[Muon][cut], VV.Yields_stat[Muon][cut], 
		       VV.Yields[ElMu][cut], VV.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Rare             | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Rare.Yields[Elec][cut], Rare.Yields_stat[Elec][cut], 
		       Rare.Yields[Muon][cut], Rare.Yields_stat[Muon][cut], 
		       Rare.Yields[ElMu][cut], Rare.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Total background | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Total.Yields[Elec][cut], Total.Yields_stat[Elec][cut], 
		       Total.Yields[Muon][cut], Total.Yields_stat[Muon][cut], 
		       Total.Yields[ElMu][cut], Total.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" TTbar dilepton   | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       TTbar.Yields[Elec][cut], TTbar.Yields_stat[Elec][cut], 
		       TTbar.Yields[Muon][cut], TTbar.Yields_stat[Muon][cut], 
		       TTbar.Yields[ElMu][cut], TTbar.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Data             | %5.0f            | %5.0f            | %5.0f            ||", 
		       Data.Yields[Elec][cut], 
		       Data.Yields[Muon][cut], 
		       Data.Yields[ElMu][cut] )<< endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    fOUTSTREAM.close();
  }
}
void TopPlotter::PrintYieldsWithMC(){
  if (gUseTTMadSpin) fOutputSubDir = "Yields_MC_Madspin/";
  else               fOutputSubDir = "Yields_MC_Madgraph/"; 

  TString yieldsfilename = "";
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  for (size_t cut=0; cut<iNCUTS; cut++){
    yieldsfilename = fOutputDir + fOutputSubDir + "Yields_MC_"+sCut[cut]+".txt";
    fOUTSTREAM.open(yieldsfilename, ios::trunc);

    fOUTSTREAM << "/////////////////////////////////////////////////////////////////////////////" << endl;
    fOUTSTREAM << " Producing MC yields for cut " << sCut[cut] << endl;
    fOUTSTREAM << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;
    
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << "                YIELDS  |    El/El   |    Mu/Mu   |    El/Mu   ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    
    float sum_mm(0.), sum_em(0.), sum_ee(0.);
    for (size_t sample=TTJets_MadSpin; sample<gNSAMPLES; sample++){
      if ( gUseTTMadSpin && sample == TTJetsFullLeptMGtauola) continue;
      if (!gUseTTMadSpin && sample == TTJets_MadSpin        ) continue;
      
      float temp_mm  = S[sample].Yields[Muon][cut]; sum_mm  += temp_mm ;
      float temp_em  = S[sample].Yields[ElMu][cut]; sum_em  += temp_em ;
      float temp_ee  = S[sample].Yields[Elec][cut]; sum_ee  += temp_ee ;
      
      fOUTSTREAM << Form("%23s |  %8.2f  |  %8.2f  |  %8.2f  ||\n", SampleName[sample].Data(),
			 temp_ee , temp_mm , temp_em );
    }
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << Form("%23s |  %8.2f  |  %8.2f  |  %8.2f  ||\n", "MC Sum", 
		       sum_ee , sum_mm , sum_em );
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << Form("%23s |  %5.0f     |  %5.0f     |  %5.0f     ||\n", "Data", 
		       S[DoubleElectron].Yields[Elec][cut], S[DoubleMu].Yields[Muon][cut], S[MuEG].Yields[ElMu][cut]); 
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    
    // Now POST a Summary...
    fOUTSTREAM << "=============================================================================" << endl;
    fOUTSTREAM << "                            SUMMARY (stat error only):                       " << endl; 
    fOUTSTREAM << "=============================================================================" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << "          Source  |       El/El      |       Mu/Mu      |       El/Mu      ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Drell-Yan        | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||", 
		       DY.Yields[Elec][cut], DY.Yields_stat[Elec][cut], 
		       DY.Yields[Muon][cut], DY.Yields_stat[Muon][cut], 
		       DY.Yields[ElMu][cut], DY.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Non W/Z leptons  | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Fake.Yields[Elec][cut], Fake.Yields_stat[Elec][cut], 
		       Fake.Yields[Muon][cut], Fake.Yields_stat[Muon][cut], 
		       Fake.Yields[ElMu][cut], Fake.Yields_stat[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << Form(" Single top quark | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       STop.Yields[Elec][cut], STop.Yields_stat[Elec][cut],
		       STop.Yields[Muon][cut], STop.Yields_stat[Muon][cut], 
		       STop.Yields[ElMu][cut], STop.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Dibosons         | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       VV.Yields[Elec][cut], VV.Yields_stat[Elec][cut], 
		       VV.Yields[Muon][cut], VV.Yields_stat[Muon][cut], 
		       VV.Yields[ElMu][cut], VV.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Rare             | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Rare.Yields[Elec][cut], Rare.Yields_stat[Elec][cut], 
		       Rare.Yields[Muon][cut], Rare.Yields_stat[Muon][cut], 
		       Rare.Yields[ElMu][cut], Rare.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Total background | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Total.Yields[Elec][cut], Total.Yields_stat[Elec][cut], 
		       Total.Yields[Muon][cut], Total.Yields_stat[Muon][cut], 
		       Total.Yields[ElMu][cut], Total.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" TTbar dilepton   | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       TTbar.Yields[Elec][cut], TTbar.Yields_stat[Elec][cut], 
		       TTbar.Yields[Muon][cut], TTbar.Yields_stat[Muon][cut], 
		       TTbar.Yields[ElMu][cut], TTbar.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Data             | %5.0f            | %5.0f            | %5.0f            ||", 
		       Data.Yields[Elec][cut], 
		       Data.Yields[Muon][cut], 
		       Data.Yields[ElMu][cut] )<< endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    fOUTSTREAM.close();
  }
}

void TopPlotter::DrawKinematicPlotsWithMC(Int_t onechan, Int_t onevar, Int_t onecut){
  if (gUseTTMadSpin) fOutputSubDir = "KinematicHistos_MC_MadSpin/";
  else               fOutputSubDir = "KinematicHistos_MC_Madgraph/";

  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);
  TString figname = "";

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1000);
  gStyle->SetOptStat("emruo");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  gROOT->LoadMacro("~folgueras/TOP/TopCode/tdrstyle.h"); 
  setTDRStyle();
  
  TLegend *leg = new TLegend(0.73,0.58,0.90,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  //leg->SetLineWidth(4);
  leg->SetTextFont(62); // Events in the leg!
  leg->SetTextSize(0.04);
  
  leg->AddEntry(Data .KinHistos[0][0][0], "Data", "PL");
  leg->AddEntry(Rare .KinHistos[0][0][0], "t#bar{t}V, VVV", "F");
  leg->AddEntry(VV   .KinHistos[0][0][0], "VV", "F");
  leg->AddEntry(Fake .KinHistos[0][0][0], "Non W/Z", "F");
  leg->AddEntry(STop .KinHistos[0][0][0], "Single Top", "F");
  leg->AddEntry(DY   .KinHistos[0][0][0], "Z/\\gamma* l^{+}l^{-}", "F");
  leg->AddEntry(TTbar.KinHistos[0][0][0],"t#bar{t}", "F");
  
  TCanvas *c1 = new TCanvas("c1","Plot");
  c1->Divide(1,2);
   
  //Plot Pad
  TPad *plot = (TPad*)c1->GetPad(1); 
  plot->SetPad(0.01, 0.23, 0.99, 0.99);
  plot->SetTopMargin(0.1);
  plot->SetRightMargin(0.04);
  
  TPad *ratio = (TPad*)c1->GetPad(2);
  ratio->SetPad(0.01, 0.02, 0.99, 0.3);
  ratio->SetGridx();
  ratio->SetGridy();
  ratio->SetTopMargin(0.05);
  ratio->SetBottomMargin(0.4);
  ratio->SetRightMargin(0.04);
  
  THStack* MC[gNCHANNELS+1][iNCUTS][gNVARS];
  for (Int_t var=0; var<gNVARS; var++){
    if (onevar > -1 && onevar != var) continue;  //print only one 
    for (Int_t ch=0; ch<gNCHANNELS+1; ch++){
      if (onechan > -1 && onechan != ch) continue; //print only one
      for (Int_t cut=0; cut<iNCUTS; cut++){
	if (onecut > -1 && onecut != cut) continue; //print only one
	
//	TTbar.KinHistos[ch][cut][var]->Scale(fLumiNorm);
//	STop .KinHistos[ch][cut][var]->Scale(fLumiNorm);
//	DY   .KinHistos[ch][cut][var]->Scale(fLumiNorm);
//	VV   .KinHistos[ch][cut][var]->Scale(fLumiNorm);
//	Rare .KinHistos[ch][cut][var]->Scale(fLumiNorm);
//	Fake .KinHistos[ch][cut][var]->Scale(fLumiNorm);

	MC[ch][cut][var] = new THStack("MC_"+gChanLabel[ch]+sCut[cut]+KinVarName[var],"");
	
	Total.KinHistos[ch][cut][var] = (TH1F*)TTbar.KinHistos[ch][cut][var]->Clone(); 
	Total.KinHistos[ch][cut][var]->Add(    STop .KinHistos[ch][cut][var]);
	Total.KinHistos[ch][cut][var]->Add(    DY   .KinHistos[ch][cut][var]);
	Total.KinHistos[ch][cut][var]->Add(    VV   .KinHistos[ch][cut][var]);
	Total.KinHistos[ch][cut][var]->Add(    Rare .KinHistos[ch][cut][var]);
	Total.KinHistos[ch][cut][var]->Add(    Fake .KinHistos[ch][cut][var]);
	
	MC[ch][cut][var]->Add(TTbar.KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(STop .KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(DY   .KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(VV   .KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(Rare .KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(Fake .KinHistos[ch][cut][var]);  
	
	plot->cd();
	
	MC[ch][cut][var]->Draw("hist");
	
	MC[ch][cut][var]->GetYaxis()->SetTitle("Events");
	MC[ch][cut][var]->GetYaxis()->SetTitleOffset(1.2);
	MC[ch][cut][var]->GetYaxis()->SetTitleSize(0.07);
	MC[ch][cut][var]->GetYaxis()->SetLabelSize(0.055);
	MC[ch][cut][var]->GetYaxis()->SetNdivisions(607);
	
	MC[ch][cut][var]->GetXaxis()->SetLabelSize(0.0);
	MC[ch][cut][var]->GetXaxis()->SetTitle("");
	
	//	Data.KinHistos[ch][cut][var]->Sumw2();
	Data.KinHistos[ch][cut][var]->SetMarkerStyle(20);
	Data.KinHistos[ch][cut][var]->Draw("SAME");
	
	leg->Draw("SAME");
	
	float maxh = Data.KinHistos[ch][cut][var]->GetMaximum();
	if (maxh < MC[ch][cut][var]->GetMaximum()) maxh = MC[ch][cut][var]->GetMaximum();
	MC[ch][cut][var]->SetMaximum(1.7*maxh);
	//	MC[ch][cut][var]->SetMinimum(50);
	
	// set logY
	if (cut == iDilepton) plot->SetLogy();
	if (cut == iZVeto   ) plot->SetLogy();
	
	ratio->cd();
	
	TH1F *H_Ratio;
	H_Ratio = (TH1F*)Data.KinHistos[ch][cut][var]->Clone();
	H_Ratio->Divide(Total.KinHistos[ch][cut][var]);
	
	H_Ratio->SetFillStyle(1001);
	H_Ratio->SetLineWidth(1);
	H_Ratio->SetFillColor(  kGray+1);
	H_Ratio->SetLineColor(  kGray+1);
	H_Ratio->SetMarkerColor(kGray+1);
	H_Ratio->SetMarkerStyle(20);
	//	  H_Ratio->SetMarkerColor(1);
	//    H_Ratio->SetLineColor(0);
	H_Ratio->SetMaximum(2);
	H_Ratio->SetMinimum(0);
	H_Ratio->SetTitle("");
	
	H_Ratio->GetYaxis()->SetTitle("Obs/Exp");
	H_Ratio->GetYaxis()->CenterTitle();
	H_Ratio->GetYaxis()->SetTitleOffset(0.45);
	H_Ratio->GetYaxis()->SetTitleSize(0.16);
	H_Ratio->GetYaxis()->SetLabelSize(0.15);
	H_Ratio->GetYaxis()->SetNdivisions(402);
	H_Ratio->GetXaxis()->SetTitleOffset(1.1);
	H_Ratio->GetXaxis()->SetLabelSize(0.15);
	H_Ratio->GetXaxis()->SetTitleSize(0.16);
	
	H_Ratio->SetMinimum(0.6);
	H_Ratio->SetMaximum(1.4);
	
	H_Ratio->DrawCopy("E2");
	TGraphErrors *thegraphRatio = new TGraphErrors(H_Ratio);
	thegraphRatio->SetFillStyle(3001);
	thegraphRatio->SetFillColor(kGray+1);
	//   
	thegraphRatio->Draw("E2 SAME");
	
	plot->cd();
	DrawTopLine(ch);
	
	TString channel = "";
	if (ch==Muon)      channel = "_MM";
	else if (ch==Elec) channel = "_EE";
	else if (ch==ElMu) channel = "_DF";
	else               channel = "_SF";
    
	figname = fOutputDir+fOutputSubDir + KinVarName[var] + channel + "_" + sCut[cut];
	c1->SaveAs(figname + ".pdf" );
	c1->SaveAs(figname + ".png" );
	c1->SaveAs(figname + ".root");
	
	plot->SetLogy(0);
      }
    }
  }
}
void TopPlotter::DrawNbjetsNjets(bool DD){
  if (DD) {
    if (gUseTTMadSpin) fOutputSubDir = "NbjetsNjets_DD_MadSpin/";
    else               fOutputSubDir = "NbjetsNjets_DD_Madgraph/";
  }
  else {
    if (gUseTTMadSpin) fOutputSubDir = "NbjetsNjets_MC_MadSpin/";
    else               fOutputSubDir = "NbjetsNjets_MC_Madgraph/";
  }
  
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);
  TString figname = "";

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1000);
  gStyle->SetOptStat("emruo");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  gROOT->LoadMacro("~folgueras/TOP/TopCode/tdrstyle.h"); 
  setTDRStyle();
  
  TLegend *leg = new TLegend(0.73,0.58,0.90,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  //leg->SetLineWidth(4);
  leg->SetTextFont(62); // Events in the leg!
  leg->SetTextSize(0.04);
  
  leg->AddEntry(Data .KinHistos[0][0][0], "Data", "PL");
  leg->AddEntry(Rare .KinHistos[0][0][0], "t#bar{t}V, VVV", "F");
  leg->AddEntry(VV   .KinHistos[0][0][0], "VV", "F");
  leg->AddEntry(Fake .KinHistos[0][0][0], "Non W/Z", "F");
  leg->AddEntry(STop .KinHistos[0][0][0], "Single Top", "F");
  leg->AddEntry(DY   .KinHistos[0][0][0], "Z/\\gamma* l^{+}l^{-}", "F");
  leg->AddEntry(TTbar.KinHistos[0][0][0],"t#bar{t}", "F");
  
  TCanvas *c1 = new TCanvas("c1","Plot");
  c1->Divide(1,2);
   
  //Plot Pad
  TPad *plot = (TPad*)c1->GetPad(1); 
  plot->SetPad(0.01, 0.23, 0.99, 0.99);
  plot->SetTopMargin(0.1);
  plot->SetRightMargin(0.04);
  
  TPad *ratio = (TPad*)c1->GetPad(2);
  ratio->SetPad(0.01, 0.02, 0.99, 0.3);
  ratio->SetGridx();
  ratio->SetGridy();
  ratio->SetTopMargin(0.05);
  ratio->SetBottomMargin(0.4);
  ratio->SetRightMargin(0.04);

  THStack* MC[gNCHANNELS][iNCUTS][gNALLSYST];
  for (Int_t ch=0; ch<gNCHANNELS; ch++){
    if (ch!=ElMu) continue;
    for (Int_t cut=0; cut<iNCUTS; cut++){
      if (cut != iDilepton) continue;
      for (size_t sys=0; sys<gNALLSYST; sys++){
	Int_t syst(sys);
	if (sys >= gNSYST) syst = 0;
	else               syst = sys;
	
	Total.SysHistos[ch][cut][sys] = (TH1F*)TTbar.SysHistos[ch][cut][sys ]->Clone();
	Total.SysHistos[ch][cut][sys]->Add(    STop .SysHistos[ch][cut][syst]);
	Total.SysHistos[ch][cut][sys]->Add(  DY   .SysHistos[ch][cut][syst]);
	if (DD) {
	  //	  Total.SysHistos[ch][cut][sys]->Add(  DY   .SysHistos[ch][cut][syst], DY_SF[ch]);
	  Total.SysHistos[ch][cut][sys]->Add(DD_NonW.SysHistos[ch][cut][syst]);
 	}
	else    {
	  Total.SysHistos[ch][cut][sys]->Add(  Fake .SysHistos[ch][cut][syst]);
	}
	Total.SysHistos[ch][cut][sys]->Add(    VV   .SysHistos[ch][cut][syst]);
	Total.SysHistos[ch][cut][sys]->Add(    Rare .SysHistos[ch][cut][syst]);
	
	MC[ch][cut][sys] = new THStack("MC_"+gChanLabel[ch]+sCut[cut]+SystName[sys],"");
	MC[ch][cut][sys]->Add(TTbar.SysHistos[ch][cut][sys ]);  
	MC[ch][cut][sys]->Add(STop .SysHistos[ch][cut][syst]);  
	MC[ch][cut][sys]->Add(DY     .SysHistos[ch][cut][syst]);
	if (DD) { 
	  //	  DY.SysHistos[ch][cut][syst]->Scale(DY_SF[ch]);
	  MC[ch][cut][sys]->Add(DD_NonW.SysHistos[ch][cut][syst]);
	}
	else    { 
	  //	  MC[ch][cut][sys]->Add(DY  .SysHistos[ch][cut][syst]);  
	  MC[ch][cut][sys]->Add(Fake.SysHistos[ch][cut][syst]);
	}
	MC[ch][cut][sys]->Add(VV   .SysHistos[ch][cut][syst]);  
	MC[ch][cut][sys]->Add(Rare .SysHistos[ch][cut][syst]);  

	plot->cd();

	MC[ch][cut][sys]->Draw("HIST");
	
	MC[ch][cut][sys]->GetYaxis()->SetTitle("Events");
	MC[ch][cut][sys]->GetYaxis()->SetTitleOffset(1.4);
	MC[ch][cut][sys]->GetYaxis()->SetTitleSize(0.05);
	MC[ch][cut][sys]->GetYaxis()->SetLabelSize(0.055);
	MC[ch][cut][sys]->GetYaxis()->SetNdivisions(607);
	
	MC[ch][cut][sys]->GetXaxis()->SetLabelSize(0.0);
	MC[ch][cut][sys]->GetXaxis()->SetTitle("");
	
	Data.KinHistos[ch][cut][syst]->Sumw2();
	Data.SysHistos[ch][cut][syst]->SetMarkerStyle(20);
	Data.SysHistos[ch][cut][syst]->Draw("SAME");
	
	leg->Draw("SAME");
	
	float maxh = Data.SysHistos[ch][cut][syst]->GetMaximum();
	if (maxh < MC[ch][cut][sys]->GetMaximum()) maxh = MC[ch][cut][sys]->GetMaximum();
	MC[ch][cut][sys]->SetMaximum(1.7*maxh);

	ratio->cd();
	
	TH1F *H_Ratio;
	H_Ratio = (TH1F*)Data.SysHistos[ch][cut][syst]->Clone();
	H_Ratio->Divide(Total.SysHistos[ch][cut][sys]);
	
	H_Ratio->SetFillStyle(1001);
	H_Ratio->SetLineWidth(1);
	H_Ratio->SetFillColor(  kGray+1);
	H_Ratio->SetLineColor(  kGray+1);
	H_Ratio->SetMarkerColor(kGray+1);
	H_Ratio->SetMarkerStyle(20);
	//	  H_Ratio->SetMarkerColor(1);
	//    H_Ratio->SetLineColor(0);
	H_Ratio->SetMaximum(2);
	H_Ratio->SetMinimum(0);
	H_Ratio->SetTitle("");
	
	H_Ratio->GetYaxis()->SetTitle("Obs/Exp");
	H_Ratio->GetYaxis()->CenterTitle();
	H_Ratio->GetYaxis()->SetTitleOffset(0.45);
	H_Ratio->GetYaxis()->SetTitleSize(0.16);
	H_Ratio->GetYaxis()->SetLabelSize(0.15);
	H_Ratio->GetYaxis()->SetNdivisions(402);
	H_Ratio->GetXaxis()->SetTitleOffset(1.1);
	H_Ratio->GetXaxis()->SetLabelSize(0.15);
	H_Ratio->GetXaxis()->SetTitleSize(0.16);
	
	H_Ratio->SetMinimum(0.6);
	H_Ratio->SetMaximum(1.4);
	
	H_Ratio->DrawCopy("E2");
	TGraphErrors *thegraphRatio = new TGraphErrors(H_Ratio);
	thegraphRatio->SetFillStyle(3001);
	thegraphRatio->SetFillColor(kGray+1);
	//   
	thegraphRatio->Draw("E2 SAME");
	
	plot->cd();
	DrawTopLine(ch);
   
	TString channel = "";
	if      (ch==Muon) channel = "_MM";
	else if (ch==Elec) channel = "_EE";
	else if (ch==ElMu) channel = "_DF";
	else               channel = "_SF";
	
	figname = fOutputDir+fOutputSubDir + "NbtagsNjets" + channel + "_" + sCut[cut] + "_" + SystName[sys];
	c1->SaveAs(figname + ".pdf" );
	c1->SaveAs(figname + ".png" );
	c1->SaveAs(figname + ".root");
      }
    }
  }
}
void TopPlotter::SaveHistosForLH(bool DD){
  TString filename = fOutputDir+"HistosForHL_MC.root";
  if (DD) filename = fOutputDir+"HistosForHL_DD.root";
  
  TFile *hfile = new TFile(filename,"RECREATE");
  
  TH1F *fHData;
  TH1F *fHTTbar[gNALLSYST];
  TH1F *fHStop [gNALLSYST];
  TH1F *fHDY   [gNALLSYST];
  TH1F *fHVV   [gNALLSYST];
  TH1F *fHRare [gNALLSYST];
  TH1F *fHFake [gNALLSYST];
  
  fHData = (TH1F*) Data.SysHistos[ElMu][iDilepton][0]->Clone("NJetsNBjets__DATA");
  fHData->Write();
  for (size_t sys=0; sys<gNALLSYST; sys++){
    Int_t syst = sys;
    if (sys < gNSYST) syst = sys;
    else              syst = 0; 
    
    fHTTbar[sys ] = (TH1F*)TTbar.SysHistos[ElMu][iDilepton][sys ]->Clone("NJetsNBjets__"+TTbar.name + sysname[sys]);
    fHStop [syst] = (TH1F*)STop .SysHistos[ElMu][iDilepton][syst]->Clone("NJetsNBjets__"+STop .name + sysname[sys]);
    fHDY   [syst] = (TH1F*)DY   .SysHistos[ElMu][iDilepton][syst]->Clone("NJetsNBjets__"+DY   .name + sysname[sys]);
    fHVV   [syst] = (TH1F*)VV   .SysHistos[ElMu][iDilepton][syst]->Clone("NJetsNBjets__"+VV   .name + sysname[sys]);
    if (DD) {
      fHFake[syst] = (TH1F*)DD_NonW.SysHistos[ElMu][iDilepton][syst]->Clone("NJetsNBjets__"+Fake.name+sysname[sys]);
    } else {
      fHFake[syst] = (TH1F*)  Fake.SysHistos[ElMu][iDilepton][syst]->Clone("NJetsNBjets__"+Fake.name+sysname[sys]);
    }
    fHRare [syst] = (TH1F*)Rare .SysHistos[ElMu][iDilepton][syst]->Clone("NJetsNBjets__"+Rare .name + sysname[sys]);
    
    fHTTbar[sys ]->Write();
    fHStop [syst]->Write();
    fHDY   [syst]->Write();
    fHVV   [syst]->Write();
    fHRare [syst]->Write();
    fHFake [syst]->Write();
  }
  hfile->Close();
}
void TopPlotter::CalculateSystematicErrors(Categories &C, Int_t cut){
  
  // SF ID and Iso systematics.. (Still hard-coded).
  if (C.name != "Total"){
    C.SystError[Muon][SFIDISO]  = (2.39/100.0)/(0.999); // * 3);
    C.SystError[Muon][SFTrig]   = (1.25/100.0)/(0.967); // * 2);

    C.SystError[Elec][SFIDISO]  = (2.34/100.0)/(0.980); // * 3);
    C.SystError[Elec][SFTrig]   = (1.65/100.0)/(0.974); // * 2);

    C.SystError[ElMu][SFIDISO]  = (1.67/100.0)/(0.990); // * 3);
    C.SystError[ElMu][SFTrig]   = (1.44/100.0)/(0.953); // * 2);
  }
  
  // Rest 
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    if (C.name == "Total"){
      if (STop.Yields_syst[ch][cut][0] == 0) { 
	cout << "[ERROR]: Please Make sure you're calculating Systematics correctly" << endl;
	return;
      }
      C.Yields_syst[ch][cut][0]  = STop   .Yields_syst[ch][cut][0] * STop   .Yields_syst[ch][cut][0];
      C.Yields_syst[ch][cut][0] += VV     .Yields_syst[ch][cut][0] * VV     .Yields_syst[ch][cut][0];
      C.Yields_syst[ch][cut][0] += DY     .Yields_syst[ch][cut][0] * DY     .Yields_syst[ch][cut][0];
      C.Yields_syst[ch][cut][0] += Rare   .Yields_syst[ch][cut][0] * Rare   .Yields_syst[ch][cut][0];
      C.Yields_syst[ch][cut][0] += Fake   .Yields_syst[ch][cut][0] * Fake   .Yields_syst[ch][cut][0];
      C.Yields_syst[ch][cut][0] += DD_NonW.Yields_syst[ch][cut][0] * DD_NonW.Yields_syst[ch][cut][0];
      
      C.Yields_syst[ch][cut][0]  = TMath::Sqrt(C.Yields_syst[ch][cut][0]);
    }
    else {
      if (C.name == "TTbar"){
	C.SystError[ch][Q2]       = TMath::Max(TMath::Abs(S[TTJets_scaleup].Yields[ch][cut]   - 
							  S[TTJets_MadSpin].Yields[ch][cut]),
					       TMath::Abs(S[TTJets_scaledown].Yields[ch][cut] - 
							  S[TTJets_MadSpin].Yields[ch][cut])
					       ) / (2 * S[TTJets_MadSpin].Yields[ch][cut]);

//SANTI  	cout << gChanLabel[ch] + " (SCALE)  = ";
//SANTI  	cout << 100 * (S[TTJets_scaleup].Yields[ch][iDilepton]   - S[TTJets_MadSpin].Yields[ch][iDilepton]) / S[TTJets_MadSpin].Yields[ch][iDilepton] << " (up) "   ;
//SANTI  	cout << 100 * (S[TTJets_scaledown].Yields[ch][iDilepton] - S[TTJets_MadSpin].Yields[ch][iDilepton]) / S[TTJets_MadSpin].Yields[ch][iDilepton] << " (down) " ;
//SANTI  	cout << endl;
//SANTI  	
//SANTI  	cout << gChanLabel[ch] + " (MATCHING) = ";
//SANTI  	cout << 100 * (S[TTJets_matchingup].Yields[ch][iDilepton]   - S[TTJets_MadSpin].Yields[ch][iDilepton]) / S[TTJets_MadSpin].Yields[ch][iDilepton] << " (up) "   ;
//SANTI  	cout << 100 * (S[TTJets_matchingdown].Yields[ch][iDilepton] - S[TTJets_MadSpin].Yields[ch][iDilepton]) / S[TTJets_MadSpin].Yields[ch][iDilepton] << " (down) " ;
//SANTI  	cout << endl;

	C.SystError[ch][Matching] = TMath::Max(TMath::Abs(S[TTJets_matchingup].Yields[ch][cut]   - 
							  S[TTJets_MadSpin].Yields[ch][cut]),
					       TMath::Abs(S[TTJets_matchingdown].Yields[ch][cut] - 
							  S[TTJets_MadSpin].Yields[ch][cut])
					       ) / (2 * S[TTJets_MadSpin].Yields[ch][cut]);
      }
      else {
	C.SystError[ch][PU]       = 0.;
	C.SystError[ch][BR]       = 0.;
	C.SystError[ch][Hadro]    = 0.;
      }
      
      C.SystError[ch][les]      = TMath::Max(TMath::Abs(C.Yields_syst[ch][cut][LESUp] - C.Yields[ch][cut]), 
					     TMath::Abs(C.Yields_syst[ch][cut][LESDown] - C.Yields[ch][cut])
					     ) / C.Yields[ch][cut];
      C.SystError[ch][jes]      = TMath::Max(TMath::Abs(C.Yields_syst[ch][cut][JESUp] - C.Yields[ch][cut]), 
					     TMath::Abs(C.Yields_syst[ch][cut][JESDown] - C.Yields[ch][cut])
					     ) / C.Yields[ch][cut];
      C.SystError[ch][jer]      = (TMath::Abs(C.Yields_syst[ch][cut][JER] - C.Yields[ch][cut])) / C.Yields[ch][cut];
      
      C.SystError[ch][btag]     = TMath::Max(TMath::Abs(C.Yields_syst[ch][cut][BtagUp] - C.Yields[ch][cut]), 
					     TMath::Abs(C.Yields_syst[ch][cut][BtagDown] - C.Yields[ch][cut])
					     ) / C.Yields[ch][cut];
    }
    
    if (C.name != "Total"){
      for (size_t sys=0; sys<gNSYSTERRTypes; sys++){ 
	C.Yields_syst[ch][cut][0] += C.SystError[ch][sys] *  C.SystError[ch][sys];
      }
      C.Yields_syst[ch][cut][0] = C.Yields[ch][cut] * TMath::Sqrt(C.Yields_syst[ch][cut][0]);
    }
  }
  return;
}
void TopPlotter::CalculateDYBkg(){
  fOutputSubDir = "DataDriven/";
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  TString yieldsfilename = "";
  Float_t R    [gNCHANNELS];
  Float_t N_in [gNCHANNELS];
  Float_t N_out[gNCHANNELS];
  Float_t k_ll [gNCHANNELS];
  Float_t SF   [gNCHANNELS];
    
  // Calculate R:
  Int_t low_in = DY.MllHistos[Muon]->FindBin(76.);
  Int_t  up_in = DY.MllHistos[Muon]->FindBin(106.);
  
  Float_t nout_ee(0.),nin_ee(0.);
  Float_t nout_mm(0.),nin_mm(0.);
  nin_mm  = DY.MllHistos[Muon]->Integral(low_in, up_in);
  nin_ee  = DY.MllHistos[Elec]->Integral(low_in, up_in);

  nout_mm = DY.MllHistos[Muon]->Integral()-nin_mm;
  nout_ee = DY.MllHistos[Elec]->Integral()-nin_ee;
  cout << "nin  (MC): "<< nin_mm  << " - " << nin_ee  << endl;
  cout << "nout (MC): "<< nout_mm << " - " << nout_ee << endl;
  
  R    [Muon] = nout_mm / nin_mm;
  R    [Elec] = nout_ee / nin_ee;
  cout << "R : " << R[Muon] << " , " << R[Elec] << endl;

  N_in [Muon] = Data.MllHistos[Muon]->Integral(low_in, up_in); 
  N_in [Elec] = Data.MllHistos[Elec]->Integral(low_in, up_in);
  N_in [ElMu] = Data.MllHistos[ElMu]->Integral(low_in, up_in);
  cout << "N_in: " << N_in[Muon] << ", " << N_in[Elec] << " , " << N_in[ElMu] << endl;
  
  k_ll [Muon] = 1; //TMath::Sqrt(nin_mm / nin_ee  );
  k_ll [Elec] = 1; //TMath::Sqrt(nin_ee / nin_mm  );
  cout << "k_ll: " << k_ll[Muon] << ", " << k_ll[Elec] << endl;
  
  // CALCULATE DD ESTIMATES... 
  N_out[Muon] = R[Muon] * (N_in[Muon] - 0.5 * N_in[ElMu] * k_ll[Muon]);
  N_out[Elec] = R[Elec] * (N_in[Elec] - 0.5 * N_in[ElMu] * k_ll[Elec]);
  
  cout << "N_out: "<< N_out[Muon] << ", " << N_out[Elec] << endl;
  cout << "nout / N_out: "<< (float) nout_mm / N_out[Muon] << ", " << (float) nout_ee / N_out[Elec] << endl;
  
  SF   [Muon] = nout_mm / N_out[Muon];
  SF   [Elec] = nout_ee / N_out[Elec];
  SF   [ElMu] = TMath::Sqrt(SF[Elec] * SF[Muon]);
  cout << "SF: "<< SF[Muon] << ", " << SF[Elec] << ", " << SF[ElMu] << endl;  

  iCut cut = i1btag;
  yieldsfilename = fOutputDir + fOutputSubDir + "DY.txt";
  fOUTSTREAM.open(yieldsfilename, ios::trunc);
  fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
  fOUTSTREAM << "                  |       El/El      |       Mu/Mu      |       El/Mu      ||" << endl;
  fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form(" Drell-Yan        | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		     DY.Yields[Elec][cut], DY.Yields_stat[Elec][cut], 
		     DY.Yields[Muon][cut], DY.Yields_stat[Muon][cut], 
		     DY.Yields[ElMu][cut], DY.Yields_stat[ElMu][cut]) << endl;
  fOUTSTREAM << Form(" Drell-Yan (DD)   | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		     DY.Yields[Elec][cut] * SF[Elec], DY.Yields_stat[Elec][cut] * SF[Elec], 
		     DY.Yields[Muon][cut] * SF[Muon], DY.Yields_stat[Muon][cut] * SF[Muon], 
		     DY.Yields[ElMu][cut] * SF[ElMu], DY.Yields_stat[ElMu][cut] * SF[ElMu]) << endl;
  fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form(" R(out/in)        | %5.3f           | %5.3f           |                 ||",
		     R[Elec], R[Muon]) << endl;
  fOUTSTREAM << Form(" SF (DD/MC)       | %5.3f           | %5.3f           | %5.3f           ||",
		     SF[Elec], SF[Muon], SF[ElMu]) << endl;
  fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
  fOUTSTREAM << endl;
  
  fOUTSTREAM << endl;
  fOUTSTREAM.close();
  
  DY_SF[Elec] = SF[Elec];
  DY_SF[Muon] = SF[Muon];
  DY_SF[ElMu] = SF[ElMu];
}
void TopPlotter::CalculateNonWZLeptonsBkg(){
  fOutputSubDir = "DataDriven/";
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  TString yieldsfilename = "";
  Float_t R    [gNCHANNELS][iNCUTS];
  Float_t R_err[gNCHANNELS][iNCUTS];
  for (size_t cut=0; cut<iNCUTS; cut++){
    yieldsfilename = fOutputDir + fOutputSubDir + "NonWZ_"+sCut[cut]+".txt";
    fOUTSTREAM.open(yieldsfilename, ios::trunc);
    
    // Calculate R:
    for (size_t ch=0; ch<gNCHANNELS; ch++){
      R    [ch][cut] = Fake.Yields[ch][cut] / Fake.SSYields[ch][cut];
      R_err[ch][cut] = Fake.Yields_stat[ch][cut]*Fake.Yields_stat[ch][cut] / (Fake.SSYields[ch][cut]*Fake.SSYields[ch][cut]) + Fake.SSYields_stat[ch][cut]*Fake.SSYields_stat[ch][cut]*(R[ch][cut] / Fake.SSYields[ch][cut])*(R[ch][cut] / Fake.SSYields[ch][cut]);
      R_err[ch][cut] = TMath::Sqrt(R_err[ch][cut]);

      DD_NonW.Yields[ch][cut]      = R[ch][cut] * Data.SSYields[ch][cut]-Total.SSYields[ch][cut];
      DD_NonW.Yields_stat[ch][cut] = TMath::Sqrt(TMath::Power(R[ch][cut]*Total.SSYields_stat[ch][cut],2) + TMath::Power(R_err[ch][cut]*(Data.SSYields[ch][cut]-Total.SSYields[Elec][cut]),2));
      DD_NonW.Yields_syst[ch][cut][0] = TMath::Abs(1 - R[ch][cut])*DD_NonW.Yields[ch][cut];
    }
    
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << "          Source  |       El/El      |       Mu/Mu      |       El/Mu      ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" TTbar dilepton   | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       TTbar.SSYields[Elec][cut], TTbar.SSYields_stat[Elec][cut], 
		       TTbar.SSYields[Muon][cut], TTbar.SSYields_stat[Muon][cut], 
		       TTbar.SSYields[ElMu][cut], TTbar.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Drell-Yan        | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||", 
		       DY.SSYields[Elec][cut], DY.SSYields_stat[Elec][cut], 
		       DY.SSYields[Muon][cut], DY.SSYields_stat[Muon][cut], 
		       DY.SSYields[ElMu][cut], DY.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Single top quark | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       STop.SSYields[Elec][cut], STop.SSYields_stat[Elec][cut],
		       STop.SSYields[Muon][cut], STop.SSYields_stat[Muon][cut], 
		       STop.SSYields[ElMu][cut], STop.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Dibosons         | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       VV.SSYields[Elec][cut], VV.SSYields_stat[Elec][cut], 
		       VV.SSYields[Muon][cut], VV.SSYields_stat[Muon][cut], 
		       VV.SSYields[ElMu][cut], VV.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Rare             | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Rare.SSYields[Elec][cut], Rare.SSYields_stat[Elec][cut], 
		       Rare.SSYields[Muon][cut], Rare.SSYields_stat[Muon][cut], 
		       Rare.SSYields[ElMu][cut], Rare.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Total background | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Total.SSYields[Elec][cut], Total.SSYields_stat[Elec][cut], 
		       Total.SSYields[Muon][cut], Total.SSYields_stat[Muon][cut], 
		       Total.SSYields[ElMu][cut], Total.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Non W/Z lep (SS) | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Fake.SSYields[Elec][cut], Fake.SSYields_stat[Elec][cut], 
		       Fake.SSYields[Muon][cut], Fake.SSYields_stat[Muon][cut], 
		       Fake.SSYields[ElMu][cut], Fake.SSYields_stat[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << Form(" Non W/Z lep (OS) | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Fake.Yields[Elec][cut], Fake.Yields_stat[Elec][cut], 
		       Fake.Yields[Muon][cut], Fake.Yields_stat[Muon][cut], 
		       Fake.Yields[ElMu][cut], Fake.Yields_stat[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << Form(" R (OS/SS)        | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       R[Elec][cut], R_err[Elec][cut], 
		       R[Muon][cut], R_err[Muon][cut], 
		       R[ElMu][cut], R_err[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Data             | %5.0f            | %5.0f            | %5.0f            ||", 
		       Data.SSYields[Elec][cut], 
		       Data.SSYields[Muon][cut], 
		       Data.SSYields[ElMu][cut] )<< endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" SS prediction    | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Data.SSYields[Elec][cut]-Total.SSYields[Elec][cut], Total.SSYields_stat[Elec][cut],
		       Data.SSYields[Muon][cut]-Total.SSYields[Muon][cut], Total.SSYields_stat[Muon][cut],
		       Data.SSYields[ElMu][cut]-Total.SSYields[ElMu][cut], Total.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" SS x R(OSSS)     | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       DD_NonW.Yields[Elec][cut], DD_NonW.Yields_stat[Elec][cut],
		       DD_NonW.Yields[Muon][cut], DD_NonW.Yields_stat[Muon][cut],
		       DD_NonW.Yields[ElMu][cut], DD_NonW.Yields_stat[ElMu][cut])
      //SANTI		     R[Elec][cut] * Data.SSYields[Elec][cut]-Total.SSYields[Elec][cut], 
      //SANTI		     TMath::Sqrt(TMath::Power(R[Elec][cut]*Total.SSYields_stat[Elec][cut],2) + 
      //SANTI		     TMath::Power(R_err[Elec][cut]*(Data.SSYields[Elec][cut]-Total.SSYields[Elec][cut]),2)),
      //SANTI		     R[Muon][cut] * Data.SSYields[Muon][cut]-Total.SSYields[Muon][cut], 
      //SANTI		     TMath::Sqrt(TMath::Power(R[Muon][cut]*Total.SSYields_stat[Muon][cut],2) + 
      //SANTI		     TMath::Power(R_err[Muon][cut]*(Data.SSYields[Muon][cut]-Total.SSYields[Muon][cut]),2)),
      //SANTI		     R[ElMu][cut] * Data.SSYields[ElMu][cut]-Total.SSYields[ElMu][cut], 
      //SANTI		     TMath::Sqrt(TMath::Power(R[ElMu][cut]*Total.SSYields_stat[ElMu][cut],2) + 
      //SANTI		     TMath::Power(R_err[ElMu][cut]*(Data.SSYields[ElMu][cut]-Total.SSYields[ElMu][cut]),2)))
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    
    fOUTSTREAM << endl;
    fOUTSTREAM.close();
  }
  /////////////////////////////////////////////////
  ////  TEMPLATES 
  /////////////////////////////////////////////////
  for (size_t cut=0; cut<iNCUTS; cut++){
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      for (size_t sys=0; sys<gNSYST; sys++){
	
	DD_NonW.SysHistos[chan][cut][sys] = (TH1F*) Fake.SSSysHistos[chan][cut][sys]->Clone();
	DD_NonW.SysHistos[chan][cut][sys] ->Add(    Fake.SSSysHistos[chan][cut][sys], -1);
	DD_NonW.SysHistos[chan][cut][sys] ->Add(    Data.SSSysHistos[chan][cut][sys],  1);
	DD_NonW.SysHistos[chan][cut][sys] ->Add(    STop.SSSysHistos[chan][cut][sys], -1);
	DD_NonW.SysHistos[chan][cut][sys] ->Add(    VV  .SSSysHistos[chan][cut][sys], -1);
	DD_NonW.SysHistos[chan][cut][sys] ->Add(    DY  .SSSysHistos[chan][cut][sys], -1);
	DD_NonW.SysHistos[chan][cut][sys] ->Add(    Rare.SSSysHistos[chan][cut][sys], -1);
	DD_NonW.SysHistos[chan][cut][sys] ->Add(    STop.SSSysHistos[chan][cut][sys], -1);
      }
    }
  }
}
void TopPlotter::CalculateCrossSection(Bool_t DD){
  // CALCULATE XSECTION without DD calculation (default);
  fOutputSubDir = "XSection/";
  TString filename = "";
  iCut cut = i1btag; // //iZVeto; //i2jets; //i1btag;
  const char* scut = sCut[cut].Data();
  if (DD) {
    if (gUseTTMadSpin) 
      filename = fOutputDir+fOutputSubDir+Form("Cross_Section_%4.1f_MadSpin_%s_DD.txt",fLumiNorm/1000.,scut);
    else 
      filename = fOutputDir+fOutputSubDir+Form("Cross_Section_%4.1f_Madgraph_%s_DD.txt",fLumiNorm/1000.,scut);
  }
  else {
    if (gUseTTMadSpin) 
      filename = fOutputDir+fOutputSubDir+Form("Cross_Section_%4.1f_MadSpin_%s_MC.txt",fLumiNorm/1000.,scut);
    else 
      filename = fOutputDir+fOutputSubDir+Form("Cross_Section_%4.1f_Madgraph_%s_MC.txt",fLumiNorm/1000.,scut);
  }
  
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);
  
   
  // Get Systematic Errors
  ResetSystematicErrors();

  cout << "Calculating Systematic Errors..." << endl;
  CalculateSystematicErrors(TTbar, cut);
  CalculateSystematicErrors(STop , cut);
  CalculateSystematicErrors(VV   , cut);
  CalculateSystematicErrors(DY   , cut);
  CalculateSystematicErrors(Rare , cut);
  if (!DD) CalculateSystematicErrors(Fake , cut);
  else {
    for (Int_t ch = 0; ch<gNCHANNELS; ch++){
      Total.Yields[ch][cut] -= Fake   .Yields[ch][cut];
      Total.Yields[ch][cut] += DD_NonW.Yields[ch][cut];
    }
  }
  CalculateSystematicErrors(Total, cut);
  
  // Get ttbar XSection
  for (Int_t ch = 0; ch<gNCHANNELS; ch++){
    ttbar.xsec     [ch] = ttbar_TLWG * (Data.Yields[ch][cut] - Total.Yields[ch][cut]) / TTbar.Yields[ch][cut]; 
    
    // statistical error
    ttbar.xsec_stat[ch] = ttbar.xsec[ch] * 
      TMath::Sqrt(Data.Yields[ch][cut]) / (Data.Yields[ch][cut] - Total.Yields[ch][cut]);
    
    // systematic error
    ttbar.xsec_syst[ch] = ttbar.xsec[ch] * 
      TMath::Sqrt((Total.Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
		  (Total.Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) + 
		  (TTbar.Yields_syst[ch][cut][0] / TTbar.Yields[ch][cut]) * 
		  (TTbar.Yields_syst[ch][cut][0] / TTbar.Yields[ch][cut]));
    ttbar.xsec_lumi[ch] = ttbar.xsec[ch]  * 0.026; // Lumi uncertainty 
    
    // Acceptances...
    ttbar.acc      [ch] = TTbar.Yields     [ch][cut]    / (fLumiNorm * ttbar_TLWG);
    ttbar.acc_stat [ch] = TTbar.Yields_stat[ch][cut]    / (fLumiNorm * ttbar_TLWG);
    ttbar.acc_syst [ch] = TTbar.Yields_syst[ch][cut][0] / (fLumiNorm * ttbar_TLWG);
    
    // Loading systematics...
    ttbar.err_VV   [ch] = ttbar.xsec[ch] * 
      TMath::Sqrt((VV   .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
     		  (VV   .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    ttbar.err_STop [ch] = ttbar.xsec[ch] * 
      TMath::Sqrt((STop .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
     		  (STop .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    if (DD) {
      ttbar.err_Fake [ch] = ttbar.xsec[ch] * 
	TMath::Sqrt((DD_NonW.Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
		    (DD_NonW.Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    }
    else {
      ttbar.err_Fake [ch] = ttbar.xsec[ch] * 
	TMath::Sqrt((Fake .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
		    (Fake .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    }
    ttbar.err_Rare [ch] = ttbar.xsec[ch] * 
      TMath::Sqrt((Rare .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
		  (Rare .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    
    ttbar.err_IDIso[ch] = ttbar.xsec[ch] * TTbar.SystError[ch][SFIDISO];
    ttbar.err_Trig [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][SFTrig];
    ttbar.err_LES  [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][les];
    ttbar.err_JES  [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][jes];
    ttbar.err_JER  [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][jer];
    ttbar.err_Btag [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][btag];
    ttbar.err_PU   [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][PU];
    ttbar.err_Q2   [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][Q2];
    ttbar.err_Match[ch] = ttbar.xsec[ch] * TTbar.SystError[ch][Matching];
  }
  
  // print systematic errors
  //  PrintSystematicErrors();
  
  fOUTSTREAM.open(filename, ios::trunc);

  fOUTSTREAM << "=======================================================================" << endl;
  fOUTSTREAM << "                                   Systematic Uncertainties (pb)       " << endl;
  fOUTSTREAM << " Source                |     El/El     |     Mu/Mu     |    El/Mu      " << endl;
  fOUTSTREAM << "-----------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form(" VV                    | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%)  ", 
		      ttbar.err_VV[Elec],  100 * ttbar.err_VV[Elec] / ttbar.xsec[Elec],
		      ttbar.err_VV[Muon],  100 * ttbar.err_VV[Muon] / ttbar.xsec[Muon],
		      ttbar.err_VV[ElMu],  100 * ttbar.err_VV[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" STop                  | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_STop[Elec],  100 * ttbar.err_STop[Elec] / ttbar.xsec[Elec],
		      ttbar.err_STop[Muon],  100 * ttbar.err_STop[Muon] / ttbar.xsec[Muon],
		      ttbar.err_STop[ElMu],  100 * ttbar.err_STop[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Fake                  | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_Fake[Elec],  100 * ttbar.err_Fake[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Fake[Muon],  100 * ttbar.err_Fake[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Fake[ElMu],  100 * ttbar.err_Fake[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Rare                  | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_Rare[Elec],  100 * ttbar.err_Rare[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Rare[Muon],  100 * ttbar.err_Rare[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Rare[ElMu],  100 * ttbar.err_Rare[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Lepton Efficiencies   | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_IDIso[Elec],  100 * ttbar.err_IDIso[Elec] / ttbar.xsec[Elec],
		      ttbar.err_IDIso[Muon],  100 * ttbar.err_IDIso[Muon] / ttbar.xsec[Muon],
		      ttbar.err_IDIso[ElMu],  100 * ttbar.err_IDIso[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Trigger Efficiencies  | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_Trig[Elec],  100 * ttbar.err_Trig[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Trig[Muon],  100 * ttbar.err_Trig[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Trig[ElMu],  100 * ttbar.err_Trig[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Lepton Energy Scale   | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_LES[Elec],  100 * ttbar.err_LES[Elec] / ttbar.xsec[Elec],
		      ttbar.err_LES[Muon],  100 * ttbar.err_LES[Muon] / ttbar.xsec[Muon],
		      ttbar.err_LES[ElMu],  100 * ttbar.err_LES[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Jet Energy Scale      | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_JES[Elec],  100 * ttbar.err_JES[Elec] / ttbar.xsec[Elec],
		      ttbar.err_JES[Muon],  100 * ttbar.err_JES[Muon] / ttbar.xsec[Muon],
		      ttbar.err_JES[ElMu],  100 * ttbar.err_JES[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Jet Energy Resolution | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_JER[Elec],  100 * ttbar.err_JER[Elec] / ttbar.xsec[Elec],
		      ttbar.err_JER[Muon],  100 * ttbar.err_JER[Muon] / ttbar.xsec[Muon],
		      ttbar.err_JER[ElMu],  100 * ttbar.err_JER[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" b-tagging             | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_Btag[Elec],  100 * ttbar.err_Btag[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Btag[Muon],  100 * ttbar.err_Btag[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Btag[ElMu],  100 * ttbar.err_Btag[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Pile Up               | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_PU[Elec],  100 * ttbar.err_PU[Elec] / ttbar.xsec[Elec],
		      ttbar.err_PU[Muon],  100 * ttbar.err_PU[Muon] / ttbar.xsec[Muon],
		      ttbar.err_PU[ElMu],  100 * ttbar.err_PU[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" QCD scale             | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_Q2[Elec],  100 * ttbar.err_Q2[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Q2[Muon],  100 * ttbar.err_Q2[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Q2[ElMu],  100 * ttbar.err_Q2[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Matching partons      | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) | %4.2f (%4.2f %%) ", 
		      ttbar.err_Match[Elec],  100 * ttbar.err_Match[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Match[Muon],  100 * ttbar.err_Match[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Match[ElMu],  100 * ttbar.err_Match[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << "=======================================================================" << endl;
  fOUTSTREAM << endl;
  fOUTSTREAM << endl;

  fOUTSTREAM << "============================================================================================================" << endl;
  fOUTSTREAM << "          Source  |            El/El            |            Mu/Mu            |            El/Mu            " << endl;
  fOUTSTREAM << "------------------------------------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form(" Drell-Yan        | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f ", 
		     DY.Yields[Elec][cut], DY.Yields_stat[Elec][cut], DY.Yields_syst[Elec][cut][0], 
		     DY.Yields[Muon][cut], DY.Yields_stat[Muon][cut], DY.Yields_syst[Muon][cut][0], 
		     DY.Yields[ElMu][cut], DY.Yields_stat[ElMu][cut], DY.Yields_syst[ElMu][cut][0])
	     << endl;
  if (DD) {
    fOUTSTREAM << Form(" Non W/Z leptons  | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f ",  
		      DD_NonW.Yields[Elec][cut], DD_NonW.Yields_stat[Elec][cut], DD_NonW.Yields_syst[Elec][cut][0], 
		      DD_NonW.Yields[Muon][cut], DD_NonW.Yields_stat[Muon][cut], DD_NonW.Yields_syst[Muon][cut][0], 
		      DD_NonW.Yields[ElMu][cut], DD_NonW.Yields_stat[ElMu][cut], DD_NonW.Yields_syst[ElMu][cut][0])
	       << endl;	
  }
  else {
    fOUTSTREAM << Form(" Non W/Z leptons  | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f ",  
		       Fake.Yields[Elec][cut], Fake.Yields_stat[Elec][cut], Fake.Yields_syst[Elec][cut][0], 
		       Fake.Yields[Muon][cut], Fake.Yields_stat[Muon][cut], Fake.Yields_syst[Muon][cut][0], 
		       Fake.Yields[ElMu][cut], Fake.Yields_stat[ElMu][cut], Fake.Yields_syst[ElMu][cut][0])
	       << endl;	
  }						  
  fOUTSTREAM << Form(" Single top quark | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f ",  
		     STop.Yields[Elec][cut], STop.Yields_stat[Elec][cut], STop.Yields_syst[Elec][cut][0],
		     STop.Yields[Muon][cut], STop.Yields_stat[Muon][cut], STop.Yields_syst[Muon][cut][0], 
		     STop.Yields[ElMu][cut], STop.Yields_stat[ElMu][cut], STop.Yields_syst[ElMu][cut][0])
	     << endl;
  fOUTSTREAM << Form(" Dibosons         | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f ",  
		     VV.Yields[Elec][cut], VV.Yields_stat[Elec][cut], VV.Yields_syst[Elec][cut][0], 
		     VV.Yields[Muon][cut], VV.Yields_stat[Muon][cut], VV.Yields_syst[Muon][cut][0], 
		     VV.Yields[ElMu][cut], VV.Yields_stat[ElMu][cut], VV.Yields_syst[ElMu][cut][0])
	     << endl;
  fOUTSTREAM << Form(" Rare             | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f ",  
		     Rare.Yields[Elec][cut], Rare.Yields_stat[Elec][cut], Rare.Yields_syst[Elec][cut][0], 
		     Rare.Yields[Muon][cut], Rare.Yields_stat[Muon][cut], Rare.Yields_syst[Muon][cut][0], 
		     Rare.Yields[ElMu][cut], Rare.Yields_stat[ElMu][cut], Rare.Yields_syst[ElMu][cut][0])
	     << endl;
  fOUTSTREAM << "------------------------------------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form(" Total background | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f ",  
		     Total.Yields[Elec][cut], Total.Yields_stat[Elec][cut], Total.Yields_syst[Elec][cut][0], 
		     Total.Yields[Muon][cut], Total.Yields_stat[Muon][cut], Total.Yields_syst[Muon][cut][0], 
		     Total.Yields[ElMu][cut], Total.Yields_stat[ElMu][cut], Total.Yields_syst[ElMu][cut][0])
	     << endl;
  fOUTSTREAM << Form(" TTbar dilepton   | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f | %7.1f +/- %4.1f +/- %6.1f ",  
		     TTbar.Yields[Elec][cut], TTbar.Yields_stat[Elec][cut], TTbar.Yields_syst[Elec][cut][0], 
		     TTbar.Yields[Muon][cut], TTbar.Yields_stat[Muon][cut], TTbar.Yields_syst[Muon][cut][0], 
		     TTbar.Yields[ElMu][cut], TTbar.Yields_stat[ElMu][cut], TTbar.Yields_syst[ElMu][cut][0])
	     << endl;
  fOUTSTREAM << "------------------------------------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form(" Data             | %5.0f                       | %5.0f                       | %5.0f                      ", 
		     Data.Yields[Elec][cut], 
		     Data.Yields[Muon][cut], 
		     Data.Yields[ElMu][cut] )<< endl;
  fOUTSTREAM << "============================================================================================================" << endl;
  fOUTSTREAM << endl;
  fOUTSTREAM << endl;
  fOUTSTREAM << "===================================================================================" << endl;
  fOUTSTREAM << "                           ttbar acceptance x eff. x BR (%)                        " << endl;
  fOUTSTREAM << "            El/El         |           Mu/Mu          |            El/Mu            " << endl;
  fOUTSTREAM << "-----------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form("    %7.4f +/- %5.4f    |    %7.3f +/- %5.4f    |    %7.3f +/- %5.4f ",
		     100*ttbar.acc[Elec], 100*TMath::Sqrt(ttbar.acc_stat[Elec]*ttbar.acc_stat[Elec]+ttbar.acc_syst[Elec]*ttbar.acc_syst[Elec]),
		     100*ttbar.acc[Muon], 100*TMath::Sqrt(ttbar.acc_stat[Muon]*ttbar.acc_stat[Muon]+ttbar.acc_syst[Muon]*ttbar.acc_syst[Muon]),
		     100*ttbar.acc[ElMu], 100*TMath::Sqrt(ttbar.acc_stat[ElMu]*ttbar.acc_stat[ElMu]+ttbar.acc_syst[ElMu]*ttbar.acc_syst[ElMu])) 
	     << endl;
  fOUTSTREAM << "===================================================================================" << endl;
  
  fOUTSTREAM << endl;
  fOUTSTREAM << "==============================================================================================================" << endl;
  fOUTSTREAM << "                                         Cross-section measurement (pb)                                       " << endl;
  fOUTSTREAM << "                El/El               |               Mu/Mu                |              El/Mu                 " << endl;
  fOUTSTREAM << "--------------------------------------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form("   %4.1f +/- %3.1f +/- %4.2f +/- %3.1f   |   %4.1f +/- %3.1f +/- %4.2f +/- %3.1f   |   %4.1f +/- %3.1f +/- %4.2f +/- %3.1f ",
		     ttbar.xsec[Elec], ttbar.xsec_stat[Elec], ttbar.xsec_syst[Elec], ttbar.xsec_lumi[Elec],
		     ttbar.xsec[Muon], ttbar.xsec_stat[Muon], ttbar.xsec_syst[Muon], ttbar.xsec_lumi[Muon],
		     ttbar.xsec[ElMu], ttbar.xsec_stat[ElMu], ttbar.xsec_syst[ElMu], ttbar.xsec_lumi[ElMu]) 
	     << endl;
  fOUTSTREAM << "==============================================================================================================" << endl;
  fOUTSTREAM << endl;
  fOUTSTREAM.close();
  
  gSystem->Exec("cat "+filename);
}

void TopPlotter::DrawTopLine(Int_t chan, Float_t y){

  TString htitleCMSChannel;
  if(chan==Muon)      htitleCMSChannel="#mu#mu";
  else if(chan==Elec) htitleCMSChannel="ee";
  else if(chan==ElMu) htitleCMSChannel="e#mu";
  else                htitleCMSChannel="ee#mu#mu";
  
  TString titlelabel = Form("CMS Preliminary, #sqrt{s}=8TeV, %4.1f fb^{-1}",fLumiNorm/1000.);
  TLatex *title  = new TLatex(-20.,50.,titlelabel);
  title->SetNDC();
  title->SetTextAlign(12);
  title->SetX(0.16);
  title->SetY(y);
  title->SetTextFont(42);
  title->SetTextSize(0.04);
  title->SetTextSizePixels(22);
  title->Draw("SAME");
  
  TLatex *chtitle  = new TLatex(-20.,50.,htitleCMSChannel);
  chtitle->SetNDC();
  chtitle->SetTextAlign(12);
  chtitle->SetX(0.85);
  chtitle->SetY(y);
  chtitle->SetTextFont(42);
  chtitle->SetTextSize(0.04);
  chtitle->SetTextSizePixels(22);
  chtitle->Draw("SAME");
}
TH1F* TopPlotter::GetHisto1D(TString filename, TString histoname) {
  
  TFile* file  = TFile::Open(filename);
  if (!file) {
    std::cerr << "ERROR: Could not load file" << std::endl
	      << "                        " << filename << std::endl;
    return 0;
  }
  TH1F* h = (TH1F*) file->Get(histoname)->Clone(histoname);
  if (!h) {
    std::cerr << "ERROR[TopSelector]: Could not find histogram " 
	      << histoname << std::endl;
    return 0;
  }
  h->SetDirectory(0);
  file->Close();
  
  return h;
};
TH2F* TopPlotter::GetHisto2D(TString filename, TString histoname) {
  
  TFile* file  = TFile::Open(filename);
  if (!file) {
    std::cerr << "ERROR: Could not load file" << std::endl
	      << "                        " << filename << std::endl;
    return 0;
  }
  TH2F* h = (TH2F*) file->Get(histoname)->Clone(histoname);
  if (!h) {
    std::cerr << "ERROR[TopSelector]: Could not find histogram " 
	      << histoname << std::endl;
    return 0;
  }
  h->SetDirectory(0);
  file->Close();
  
  return h;
};
void TopPlotter::SetupDraw(TH1F* h, Int_t color, Int_t var){
  h->Rebin(GetRebin(h,var));
  
  TString xaxis = KinAxisLabel[var];
  
  h->SetTitle(h->GetTitle());
  h->GetXaxis()->SetTitle(xaxis);
  h->SetLineColor(1);
  h->SetFillColor(color);
  h->SetFillStyle(1001);
  
  h->SetBinContent(h->GetNbinsX(),(h->GetBinContent(h->GetNbinsX()+1)+h->GetBinContent(h->GetNbinsX())));
  h->SetBinContent(h->GetNbinsX()+1,0);
  
  if (var == NBTagsNJets) {  //change bin labels
    h->GetXaxis()->SetBinLabel( 1, "0");
    h->GetXaxis()->SetBinLabel( 2, "0");
    h->GetXaxis()->SetBinLabel( 3, "1");
    h->GetXaxis()->SetBinLabel( 4, "0");
    h->GetXaxis()->SetBinLabel( 5, "1");
    h->GetXaxis()->SetBinLabel( 6, "2");
    h->GetXaxis()->SetBinLabel( 7, "0");
    h->GetXaxis()->SetBinLabel( 8, "1");
    h->GetXaxis()->SetBinLabel( 9, "2");
    h->GetXaxis()->SetBinLabel(10, "3");
    h->GetXaxis()->SetBinLabel(11, "0");
    h->GetXaxis()->SetBinLabel(12, "1");
    h->GetXaxis()->SetBinLabel(13, "2");
    h->GetXaxis()->SetBinLabel(14, "3");
    h->GetXaxis()->SetBinLabel(15, "4");
    
    float bincontent(0.);
    for (int bin=15; bin<h->GetNbinsX()+1; bin++) bincontent += h->GetBinContent(bin);
    h->SetBinContent(15, bincontent);
    
    h->GetXaxis()->SetRangeUser(-0.5,14.5);
  }
  return;
}
Int_t TopPlotter::GetRebin(TH1F* h, Int_t var){
  TString histo = h->GetTitle();
  
  Int_t nbins   = h->GetNbinsX();
  float xmax  = h->GetXaxis()->GetBinUpEdge(nbins);
  float xmin  = h->GetXaxis()->GetBinLowEdge(1);

  Int_t rebin = nbins/(Int_t)(xmax-xmin); // normalize to 1GeV?
  if (var==MET        ) return 5 * rebin;
  if (var==InvMass    ) return 5 * rebin;
  if (var==Lep0Pt     ) return 5 * rebin;
  if (var==Lep1Pt     ) return 5 * rebin;
  if (var==DelLepPhi  ) return 1;
  if (var==NJets      ) return rebin;
  if (var==NBtagJets  ) return rebin;
  if (var==Jet0Pt     ) return 5 * rebin;
  if (var==Jet1Pt     ) return 5 * rebin;
  		      
  if (var==CSVTag     ) return 25;
  if (var==TopD       ) return 25;
  if (var==DelPhillJet) return 10;
  
  return rebin;
}
