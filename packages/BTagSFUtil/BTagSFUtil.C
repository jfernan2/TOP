#include "BTagSFUtil.h"
//#include "BTagSFPayloadMoriond13.C"
#include "BTagSFPayloadWinter13.C"
//#include "BTagSFPayloadWinter13NoTTbar.C"
#include "BTagEfficienciesTTbarSummer12.C"
//#include "FastSimCorrectionFactorsSummer12.C"

BTagSFUtil::BTagSFUtil(TString BTagAlgorithm, TString DataPeriod, int Seed) {
  
  rand_ = new TRandom3(Seed);

  GetBTagPayload(BTagAlgorithm, DataPeriod);

  if (TaggerCut<0.) 
    cout << "BTagSFUtil: " << BTagAlgorithm << " not a supported b-tagging algorithm" << endl;
  
  loglevel = 1;
}

BTagSFUtil::~BTagSFUtil() {

  delete rand_;

}

float BTagSFUtil::ScaleFactorB(float JetPt, int SystematicFlag) {
   
  if (JetPt<BTagPtBinEdge[0]) { cout << "SFb is not available for jet pt<20 GeV" << endl; return -1.; }

  float x = JetPt;
  if (JetPt>BTagPtBinEdge[nBTagPtBins]) x = BTagPtBinEdge[nBTagPtBins];
 
  int JetPtBin = -1;
  for (int ptbtagbin = 0; ptbtagbin<nBTagPtBins; ptbtagbin++) 
    if (x>=BTagPtBinEdge[ptbtagbin]) JetPtBin++;

  if (JetPt>x) SystematicFlag *= 2;

  float SFb = funSFb->Eval(x) + SystematicFlag*SFb_error[JetPtBin];

  return SFb;

}

float BTagSFUtil::ScaleFactorLight(float JetPt, float JetEta, int SystematicFlag) {
 
  if (JetPt<BTagPtBinEdge[0]) { cout << "SFlight is not available for jet pt<20 GeV" << endl; return -1.; }
 
  int JetEtaBin = -1;
  for (int etabtagbin = 0; etabtagbin<nBTagEtaBins; etabtagbin++) 
    if (fabs(JetEta)>=BTagEtaBinEdge[etabtagbin]) JetEtaBin++;

  float MaxJetPtLight = funSFlight[JetEtaBin][SystematicFlag+1]->GetMaximumX();

  float x = JetPt;
  if (JetPt>MaxJetPtLight) x = MaxJetPtLight;

  return funSFlight[JetEtaBin][SystematicFlag+1]->Eval(x);

}

float BTagSFUtil::ScaleFactorJet(int JetFlavor, float JetPt, float JetEta, int SystematicFlag) {
  
  float SF = -1.;
 
  if (JetPt<BTagPtBinEdge[0]) { 
    if (loglevel > 2) cout << "SF is not available for jet pt<20 GeV" << endl; 
    return -1.; 
  }
  if (fabs(JetEta)>2.4)       { 
    if (loglevel > 2) cout << "SF is not available for jet |eta|>2.4" << endl; 
    return -1.; 
  }
  
  if (abs(JetFlavor)==5) SF = ScaleFactorB(JetPt, SystematicFlag);
  else if (abs(JetFlavor)==4) SF = ScaleFactorB(JetPt, 2*SystematicFlag);
  else SF = ScaleFactorLight(JetPt, JetEta, SystematicFlag);

  if (SF == -1.) cout << "Jet parameter out of BTV prescriptions" << endl;

  return SF;

}

float BTagSFUtil::JetTagEfficiency(int JetFlavor, float JetPt, float JetEta) {

  if (abs(JetFlavor)==5) return TagEfficiencyB(JetPt, JetEta);
  else if (abs(JetFlavor)==4) return TagEfficiencyC(JetPt, JetEta);
  else return TagEfficiencyLight(JetPt, JetEta);

}

bool BTagSFUtil::IsTagged(float JetDiscriminant, int JetFlavor, float JetPt, float JetEta, int SystematicFlag) {
  
  bool isBTagged = JetDiscriminant>TaggerCut;

  if (JetFlavor==-999999) return isBTagged; // Data: no correction needed

  bool newBTag = isBTagged;

  float Btag_SF = ScaleFactorJet(JetFlavor, JetPt, JetEta, SystematicFlag);

  if (Btag_SF == 1) return newBTag; //no correction needed 

  //throw die
  float coin = rand_->Uniform(1.);    
 
  if(Btag_SF > 1){  // use this if SF>1

    if( !isBTagged ) {

      float Btag_eff = JetTagEfficiency(JetFlavor, JetPt, fabs(JetEta));

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );

      //upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }

  }else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {newBTag = false;}

  }

  return newBTag;

}


