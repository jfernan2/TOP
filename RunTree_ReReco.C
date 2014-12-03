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
		    Long64_t nEvents        = -1,
		    Float_t  stopMass       = 0.0
		    ){
  
  gROOT->LoadMacro("$PAFPATH/PAF.C");
  
  // VARIABLES TO BE USED AS PARAMETERS...
  Float_t G_Total_Lumi    = 19664.225;   // Total Luminosity
  Float_t G_Event_Weight  = 1.0;         // Event Weight
  Bool_t  G_IsData        = false;       // 1 for data, 0 for MC
  Float_t G_LumiForPUData = 19468.3;     // luminosity in http://www.hep.uniovi.es/jfernan/PUhistos
  Bool_t  DoSF            = true;
  Bool_t  DoDF            = true;
  
  cout << "Params: " << endl;
  cout << "sampleName      " << sampleName      << endl;
  cout << "lumiForNorm     " << G_Total_Lumi    << endl;
  cout << "EventWeight     " << G_Event_Weight  << endl;
  cout << "lumiForPUdata   " << G_LumiForPUData << endl;
  cout << "DoSystStudies   " << DoSystStudies   << endl;
  cout << "nEvents         " << nEvents         << endl;
  cout << "massStop        " << stopMass        << endl;
      
  // PROOF settings - see scripts/PAFOptions.h
  //----------------------------------------------------------------------------

  // PROOF mode
  //----------------------------------------------------------------------------
  if      (nSlots == 1) gPAFOptions->SetPAFMode(kSequential); // NO PROOF
  else if (nSlots <= 8) gPAFOptions->SetPAFMode(kLite);       // PROOF Lite
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
  if ((sampleName == "DoubleElectron" ||
       sampleName == "DoubleMu" ||
       sampleName == "MuEG")) {
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
  else if (sampleName == "T2tt_150to250LSP1to100_LeptonFilter") {

    if(stopMass==162.5) G_Event_Weight = (5./9.) * 46.2   /462897.;  // cross section in pb, 5/9 for BR
    if(stopMass==175.0) G_Event_Weight = (5./9.) * 36.7994/446023.;  
    if(stopMass==187.5) G_Event_Weight = (5./9.) * 25.93  /440137.;
    if(stopMass==200.0) G_Event_Weight = (5./9.) * 18.5245/343270.;  
    G_IsData = false;
    
    cout << endl;
    cout << "      x-section = " << 36.7994        << endl;
    cout << "        nevents = " << 446023.        << endl;
    cout << " base file name = " << sampleName     << endl;
    cout << "         weight = " << G_Event_Weight << endl;
    cout << endl;

    TString asample = Form("Tree_T2tt_150to250LSP1to100_LeptonFilter_*");
    cout << "   + Looking for " << asample << " trees..." << endl;
    gPAFOptions->AddDataFiles(dm->GetRealDataFiles("MC_Summer12_53X/Legacy/",asample));
  }
  else {   // Deal with MC samples
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
  Bool_t G_Use_CSVM = true;
  TString outputDir = "/mnt_pool/fanae105/user/folgueras/TOP/TopTrees/Dec03_Jet30_Lep20_CSVM/";
  
  if (DoSF && DoDF) outputDir += "";
  else if (DoSF)    outputDir += "SF/";
  else if (DoDF)    outputDir += "DF/";
  else             {cout << "ERROR, please indicate SF or DF" << endl;  return; }

  gSystem->mkdir(outputDir, kTRUE);

  std::ostringstream oss;      
  oss << G_Total_Lumi;
  
  TString mstop="";
  if(sampleName == "T2tt_150to250LSP1to100_LeptonFilter")  mstop = Form("_%f",stopMass);
  TString LumiString = oss.str();
  TString outputFile = outputDir
    + "/"
    + "Tree_Legacy_"
    + sampleName
    + mstop
    + ".root";
  
  gPAFOptions->SetOutputFile(outputFile);

  // Parameters for the analysis
  //----------------------------------------------------------------------------
  if (!G_IsData) { 
    if(gPAFOptions->GetPAFMode() != kSequential) 
      proof->Exec("gSystem->AddIncludePath(\"-D__ISMC\");"); 
    else      
      gSystem->AddIncludePath("-D__ISMC"); 
  }
  if (sampleName == "T2tt_150to250LSP1to100_LeptonFilter"){
    //cout << "this is a stop sample!!! " << endl;
    if(gPAFOptions->GetPAFMode() != kSequential)
      proof->Exec("gSystem->AddIncludePath(\"-D__ISSTOP\");");
    else
      gSystem->AddIncludePath("-D__ISSTOP");
  }
  if (sampleName == "TTbar_MCatNLO"){
    //cout << "this is a MC@NLO sample!!! " << endl;
    if(gPAFOptions->GetPAFMode() != kSequential)  
      proof->Exec("gSystem->AddIncludePath(\"-D__ISMCNLO\");");
    else
      gSystem->AddIncludePath("-D__ISMCNLO");
  }  

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
  gPAFOptions->inputParameters->SetNamedBool  ("DoSF"         , DoSF             );
  gPAFOptions->inputParameters->SetNamedBool  ("DoDF"         , DoDF             );
  gPAFOptions->inputParameters->SetNamedFloat ("stopMass"     , stopMass         );

  // Number of events (Long64_t)
  //----------------------------------------------------------------------------
  if (nSlots==1) gPAFOptions->SetNEvents(2000);
  else           gPAFOptions->SetNEvents(nEvents);

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
