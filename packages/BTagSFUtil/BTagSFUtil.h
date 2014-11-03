/*************************************************************

Class Usage:

This class should only be used for upgrading and downgrading 
if a single operating point is used in an analysis. 

Based on the BTagSFUtil from Michael Segala

*************************************************************/

#include <Riostream.h>
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"

class BTagSFUtil{

 public:
    
  BTagSFUtil(TString BTagAlgorithm, TString DataPeriod = "ABCD", int Seed = 0);
  ~BTagSFUtil();
    
  
  bool IsTagged(float JetDiscriminant, int JetFlavor, float JetPt, float JetEta, int SystematicFlag);

 private:

  void GetBTagPayload(TString BTagAlgorithm, TString DataPeriod);

  float ScaleFactorB(float JetPt, int SystematicFlag);
  float ScaleFactorLight(float JetPt, float JetEta, int SystematicFlag);
  float ScaleFactorJet(int JetFlavor, float JetPt, float JetEta, int SystematicFlag);

  float JetTagEfficiency(int JetFlavor, float JetPt, float JetEta);
  float TagEfficiencyB(float JetPt, float JetEta);
  float TagEfficiencyC(float JetPt, float JetEta);
  float TagEfficiencyLight(float JetPt, float JetEta);
  
  TRandom3* rand_;

  TString TaggerName;
  float TaggerCut;

  TF1 *funSFb, *funSFlight[4][3];

  int nBTagPtBins;
  float BTagPtBinEdge[50];
  float SFb_error[50];

  int nBTagEtaBins;
  float BTagEtaBinEdge[50];
  
  int loglevel;
  
};

