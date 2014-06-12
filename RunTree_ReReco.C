////////////////////////////////////////////////////////////////////////////////
//
//    FILE: RunProof.C
// AUTHORS: I. Gonzalez Caballero, A.Y. Rodriguez Marrero
//    DATE: 2010
//
// CONTENT: Main macro to run over MiniTrees or TESCO Trees using PROOF
//          in PROOF-Lite, PROOF-Cluster or Sequential mode
//
////////////////////////////////////////////////////////////////////////////////

TProof* proof = 0;

void RunTree_ReReco(TString  sampleName     = "TTbar_Madgraph",
		    Int_t    nSlots         =  1,
		    Bool_t   DoSystStudies  =  false,
		    Bool_t   DoFR           =  false,
		    Long64_t nEvents        = -1,
		    TString  suffix         = ""
		    ){
  
  gROOT->LoadMacro("$PAFPATH/PAF.C");
  
  // VARIABLES TO BE USED AS PARAMETERS...
  Float_t G_Total_Lumi    = 19664.225;   // Total Luminosity
  Float_t G_Event_Weight  = 1.0;         // Event Weight
  Bool_t  G_IsData        = false;       // 1 for data, 0 for MC
  Float_t G_LumiForPUData = 19468.3;     // luminosity in http://www.hep.uniovi.es/jfernan/PUhistos
  
  cout << "Params: " << endl;
  cout << "sampleName      " << sampleName      << endl;
  cout << "lumiForNorm     " << G_Total_Lumi    << endl;
  cout << "EventWeight     " << G_Event_Weight  << endl;
  cout << "lumiForPUdata   " << G_LumiForPUData << endl;
  cout << "DoSystStudies   " << DoSystStudies   << endl;
  cout << "DoFR            " << DoFR            << endl;
  cout << "nEvents         " << nEvents         << endl;
  
  // PROOF settings - see scripts/PAFOptions.h
  //----------------------------------------------------------------------------

  // PROOF mode
  //----------------------------------------------------------------------------
  if      (nSlots == 1) gPAFOptions->SetPAFMode(kSequential); // NO PROOF
  else if (nSlots < 8 ) gPAFOptions->SetPAFMode(kLite);       // PROOF Lite
  else                  gPAFOptions->SetPAFMode(kPoD);        // PoD
  
  gPAFOptions->SetNSlots(nSlots);

  //SANTI 
  //  gPAFOptions->proofMode = kSequential;       // No PROOF
  //  gPAFOptions->proofMode = kLite;             // PROOF Lite
  //  gPAFOptions->proofMode = kPoD;              // PoD

  //gPAFOptions->proofMode = kCluster;            // PROOF Cluster
  //  gPAFOptions->NSlots = 50;                   // Number of slots

  // Optional parameters for PROOF Cluster mode
  // gPAFOptions->proofServer = "proof.ifca.es";  // PROOF server
  // gPAFOptions->proofServerPort = 1093;         // PROOF port
  // gPAFOptions->maxSlavesPerNode = 9999;        // Max number of slaves / node
  
  
  // PROOF start
  //----------------------------------------------------------------------------
 cout << ">> Starting PROOF..." << endl;
 proof = InitProof();
 if (!proof && gPAFOptions->GetPAFMode() != kSequential) {
   cerr << "ERROR: I could not initialise a PROOF session!" << endl;
   return;
 }
 if (gPAFOptions->GetPAFMode() != kSequential) 
   gPAFOptions->GetPROOFSession()->SetLogLevel(2, TProofDebug::kOutput); 
  
  // Tree type
  //----------------------------------------------------------------------------
  gPAFOptions->SetTreeType(kMiniTrees);

  // Base path to input files
  //----------------------------------------------------------------------------
  TString dataPath = "/pool/ciencias/MC_Summer12_53X/Legacy/";
  
  ///////////////////////////////
  // INPUT DATA SAMPLE
  //
  TString userhome = "/mnt_pool/fanae105/user/folgueras/";
  gROOT->LoadMacro(userhome+"/Utils/DatasetManager/DatasetManager.C+");
  
  cout << ">> Setting datasets..." << endl;
  DatasetManager* dm = new DatasetManager("Legacy_Summer12_53X");
  dm->RedownloadFiles();
  
  // Deal with data samples
  if (DoFR){
    cout << "   + FR trees..." << endl;
    if (sampleName == "DoubleElectron" ||
	sampleName == "DoubleMu" ||
	sampleName == "MuEG"){
      
      TString datasuffix[] = {
	"A_876",
	"B_4412",
	"C_7016",
	"D_7360"
      };
      const unsigned int nDataSamples = 4;
      for(unsigned int i = 0; i < nDataSamples; i++) {
	TString asample = Form("Tree_%s*%s",sampleName.Data(), datasuffix[i].Data());
	cout << "   + Looking for " << asample << " trees..." << endl;
	gPAFOptions->AddDataFiles(dm->GetRealDataFiles("MC_Summer12_53X/Legacy/FR/",asample));
      }
      G_Event_Weight = 1.;
      G_IsData = true;
    }
    else {
      dm->LoadDataset(sampleName);
      G_Event_Weight = dm->GetCrossSection() / dm->GetEventsInTheSample();
      G_IsData = false;
      
      TString asample = Form("Tree_%s*",sampleName.Data());
      gPAFOptions->AddDataFiles(dm->GetRealDataFiles("MC_Summer12_53X/Legacy/FR/",asample));
      
    }
  }
  else if ((sampleName == "DoubleElectron" ||
	    sampleName == "DoubleMu" ||
	    sampleName == "MuEG") && !DoFR) {
    cout << "   + Data..." << endl;
    
    TString datasuffix[] = {
      "A_876",
      "B_4412",
      "C_7016",
      "D_7360"
    };
    const unsigned int nDataSamples = 4;
    for(unsigned int i = 0; i < nDataSamples; i++) {
      TString asample = Form("Tree_%s*%s",sampleName.Data(), datasuffix[i].Data());
      cout << "   + Looking for " << asample << " trees..." << endl;
      gPAFOptions->AddDataFiles(dm->GetRealDataFiles("MC_Summer12_53X/Legacy/",asample));
    }
    G_Event_Weight = 1.;
    G_IsData = true;
  }
  // Deal with MC samples
  else {
    dm->LoadDataset(sampleName);

    gPAFOptions->AddDataFiles(dm->GetFiles());
    
    G_Event_Weight = dm->GetCrossSection() / dm->GetEventsInTheSample();
    G_IsData = false;
    
    
    cout << endl;
    cout << "      x-section = " << dm->GetCrossSection()      << endl;
    cout << "        nevents = " << dm->GetEventsInTheSample() << endl;
    cout << " base file name = " << dm->GetBaseFileName()      << endl;
    cout << "         weight = " << G_Event_Weight             << endl;
    cout << endl;
  }
    
  // Output file name
  //----------------------------------------------------------------------------
  Bool_t G_Use_CSVM = false;
  TString outputDir = "/mnt_pool/fanae105/user/folgueras/TOP/TopTrees/Jun11_Jet30_LH_withFakesAndDY_CSVT/";
  
  //CSVM_METType0I_3rdLepV_FullSyst_Mar20_jet25/";
  gSystem->mkdir(outputDir, kTRUE);
  
  std::ostringstream oss;      
  oss << G_Total_Lumi;
  
  TString LumiString = oss.str();
  
  TString outputFile = outputDir
    + "/"
    + "Tree_Legacy_"
//    + LumiString
//    + "pb-1_"
    + sampleName
    + suffix
    + ".root";
  
  gPAFOptions->SetOutputFile(outputFile);
  //
  //SANTI  gPAFOptions->outputFile = outputFile;

  // Parameters for the analysis
  //----------------------------------------------------------------------------
  if (!G_IsData) gSystem->AddIncludePath("-D__ISMC");
  if (DoFR)      gSystem->AddIncludePath("-D__ISFR");

  // See packages/InputParameters/InputParameters.h for information on how
  // to use this class.

  gPAFOptions->inputParameters = new InputParameters();
  
  //  gPAFOptions->inputParameters->SetNamedInt("sys_source",     sys_source);
  //  gPAFOptions->inputParameters->SetNamedInt("sys_direction",  sys_direction);
   
  gPAFOptions->inputParameters->SetNamedString("sampleName",    sampleName.Data());
  gPAFOptions->inputParameters->SetNamedBool  ("IsData",        G_IsData         );
  gPAFOptions->inputParameters->SetNamedBool  ("UseCSVM",       G_Use_CSVM       );
  gPAFOptions->inputParameters->SetNamedFloat ("weight",        G_Event_Weight   );
  gPAFOptions->inputParameters->SetNamedFloat ("LumiForPU",     G_LumiForPUData  );
  gPAFOptions->inputParameters->SetNamedFloat ("TotalLumi",     G_Total_Lumi     );
  gPAFOptions->inputParameters->SetNamedBool  ("DoSystStudies", DoSystStudies    );
  gPAFOptions->inputParameters->SetNamedBool  ("DoFR"         , DoFR             );
  
  // Number of events (Long64_t)
  //----------------------------------------------------------------------------
  gPAFOptions->SetNEvents(nEvents);
  // First event (Long64_t)
  //----------------------------------------------------------------------------
  gPAFOptions->SetFirstEvent(0);
  
  // Name of analysis class
  //----------------------------------------------------------------------------
  // If 0 the default name schema will be used, i.e. depending on the value
  // of gPAFOptions->treeType: MyAnalysisTESCO or MyAnalsyisMiniTrees
  gPAFOptions->SetAnalysis("TreeAnalysisTop");

  gPAFOptions->AddPackage("PUWeight");
  gPAFOptions->AddPackage("BTagSFUtil");
  gPAFOptions->AddPackage("LeptonSF");
  //  gPAFOptions->AddPackage("GlobalVariables");
  
  // Control output and checks
  //----------------------------------------------------------------------------
  // + If true (default) PAF checks for new version in CVS every time
  // gPAFOptions->checkVersion = true;
  // + If true (default) the output file is reopened so the objects in the
  //   file can be interactively accessed. The object in the output are also
  //   listed
  gPAFOptions->ReopenOutputFile(false);
  
  // Run the analysis
  //----------------------------------------------------------------------------
  //  gPAFOptions->SetMergeThroughFile();
  if (!RunAnalysis())
    cerr << "ERROR: There was a problem running the analysis!" << endl;
}
