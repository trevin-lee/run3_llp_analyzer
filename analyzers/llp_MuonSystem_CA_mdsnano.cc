//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/Selector.hh"
#include "llp_MuonSystem_CA_mdsnano.h"
#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"

#include "CACluster.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include <iostream>
#include <fstream>
#include <random>
//C++ includes
#include "assert.h"
//ROOT includes
#include "TH1F.h"

//using namespace fastjet;
using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;

struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  // bool passLooseId;
  // bool passMediumId;
  bool passId;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;
};


struct jets
{
  TLorentzVector jet;
  float time;
  bool passId;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool isCSVL;
  int ecalNRechits;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  bool jetPassMuFrac;
  float jetPtJESUp;
  float jetPtJESDown;
  float jetEJESUp;
  float jetEJESDown;
  float JecUnc;

  float electronEnergyFraction;
  float neutralEmEnergyFraction;
  float chargedHadronEnergyFraction;
  float neutralHadronEnergyFraction;
  float muonEnergyFraction;

};

//lepton highest pt comparator
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;


void llp_MuonSystem_CA_mdsnano::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";
  
  //DEBUGGING: Create debug log file
  ofstream debugLog("debug_analyzer.log");
  debugLog << "=== Debug Log Started ===" << endl;
  debugLog << "IsData = " << isData << endl;
  debugLog << "options = " << options << endl;
  debugLog << "outputfilename = " << outputfilename << endl;
  debugLog << "analysisTag = " << analysisTag << endl;
  debugLog.flush();

  bool signalScan = int(options/10) == 1;
  int option = options%10;



  if( isData )
  {
    debugLog << "[INFO]: running on data with option: " << option << endl;
  }
  else
  {
    debugLog << "[INFO]: running on MC with option: " << option << endl;
  }
  if( signalScan )
  {
    debugLog << "[INFO]: running with Signal scan" << endl;
  }
  else
  {
    debugLog << "[INFO]: running without Signal scan " << option << endl;
  }

  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    analysisTag = "Summer24";
    debugLog << "Set analysisTag to default: " << analysisTag << endl;
    debugLog.flush();
  } else {
    debugLog << "Using provided analysisTag: " << analysisTag << endl;
    debugLog.flush();
  }

  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  debugLog << "Setting up output file..." << endl;
  debugLog.flush();
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");
  debugLog << "Output file created: " << outfilename << endl;
  debugLog.flush();

  debugLog << "Creating TreeMuonSystemCombination..." << endl;
  debugLog.flush();
  TreeMuonSystemCombination *MuonSystem = new TreeMuonSystemCombination;
  debugLog << "Calling CreateTree()..." << endl;
  debugLog.flush();
  MuonSystem->CreateTree();
  debugLog << "Setting AutoFlush..." << endl;
  debugLog.flush();
  MuonSystem->tree_->SetAutoFlush(0);
  debugLog << "Calling InitTree()..." << endl;
  debugLog.flush();
  MuonSystem->InitTree();
  debugLog << "Tree setup completed successfully." << endl;
  debugLog.flush();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*>Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  map<pair<int,int>, TH1F*> accep2D;
  map<pair<int,int>, TH1F*> accep_met2D;
  map<pair<int,int>, TH1F*> Total2D;

  //histogram containing total number of processed events (for normalization)
  debugLog << "Creating histograms..." << endl;
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_csccsc = new TH1F("accep_csccsc", "accep_csccsc", 1, 1, 2);
  TH1F *accep_cscdt = new TH1F("accep_cscdt", "accep_cscdt", 1, 1, 2);
  TH1F *accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

  TH1F *Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
  TH1F *NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
  TH1F *Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
  TH1F *Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
  TH1F *NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);
  debugLog << "Histograms created successfully." << endl;

  //JetDefinition jet_def( antikt_algorithm, .4 );
  //fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4);

  //vector<fastjet::PseudoJet> input_particles;

  debugLog << "Setting up CMSSW path..." << endl;
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  debugLog << "CMSSW_BASE environment variable: " << (cmsswPath ? cmsswPath : "NULL") << endl;
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1
  debugLog << "CMSSW path setup complete. Path: " << pathname << endl;
  debugLog << "option = " << option << endl;
  debugLog.flush();

  //--------------------------------
  //Initialize helper
  //--------------------------------
  debugLog << "Initializing RazorHelper with analysisTag='" << analysisTag << "' and isData=" << isData << endl;
  debugLog.flush();
  RazorHelper *helper = 0;
  
  // DEBUGGING: Try to skip RazorHelper if it's causing issues
  bool skipRazorHelper = false;  // Set to true to bypass RazorHelper for debugging
  
  if (!skipRazorHelper) {
    try {
      debugLog << "About to call RazorHelper constructor..." << endl;
      debugLog.flush();
      helper = new RazorHelper(analysisTag, isData);
      debugLog << "RazorHelper constructor completed successfully." << endl;
      debugLog.flush();
    } catch (const std::exception& e) {
      debugLog << "ERROR: RazorHelper constructor threw exception: " << e.what() << endl;
      debugLog.flush();
      skipRazorHelper = true;
      debugLog << "Continuing without RazorHelper for debugging..." << endl;
    } catch (...) {
      debugLog << "ERROR: RazorHelper constructor threw unknown exception" << endl;
      debugLog.flush();
      skipRazorHelper = true;
      debugLog << "Continuing without RazorHelper for debugging..." << endl;
    }
  } else {
    debugLog << "DEBUGGING: Skipping RazorHelper initialization" << endl;
    debugLog.flush();
  }
  
  if (helper) {
    debugLog << "RazorHelper initialized successfully." << endl;
  } else {
    debugLog << "WARNING: Running without RazorHelper (some features may not work)" << endl;
  }
  debugLog.flush();

  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  debugLog << "Checking fChain..." << endl;
  debugLog.flush();
  if (fChain == 0) {
    debugLog << "ERROR: fChain is null!" << endl;
    debugLog.flush();
    debugLog.close();
    return;
  }
  debugLog << "fChain is valid, getting entries..." << endl;
  debugLog.flush();
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  
  //DEBUGGING: Check if we actually have events to process
  debugLog << "Total events in chain: " << fChain->GetEntries() << endl;
  debugLog.flush();
  if (fChain->GetEntries() == 0) {
    debugLog << "ERROR: No events found in input files!" << endl;
    cout << "ERROR: No events found in input files!" << endl;
    debugLog.close();
    return;
  }
  
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  debugLog << "Starting event loop..." << endl;
  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {
    //begin event
    if(jentry % 1000 == 0)
    {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      debugLog << "Processing entry " << jentry << " (every 1000th event)" << endl;
      start = clock();
    }
    
    //DEBUGGING: Print every event for first 10 events
    if (jentry < 10) {
      debugLog << "=== DEBUGGING: Processing event " << jentry << " ===" << endl;
    }
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) {
      debugLog << "ERROR: LoadTree failed for entry " << jentry << endl;
      break;
    }
    nb = fChain->GetEntry(jentry); 
    if (nb <= 0) {
      debugLog << "ERROR: GetEntry failed for entry " << jentry << ", returned " << nb << endl;
      continue;
    }
    nbytes += nb;

    //fill normalization histogram
    MuonSystem->InitVariables();
    
    //DEBUGGING: Print basic event info for first 10 events
    if (jentry < 10) {
      debugLog << "Event " << jentry << " basic info:" << endl;
      debugLog << "  runNum = " << run << endl;
      debugLog << "  eventNum = " << event << endl;
      debugLog << "  ncscRechits = " << ncscRechits << endl;
      debugLog << "  ndtRecHits = " << ndtRecHits << endl;
      debugLog << "  nMuon = " << nMuon << endl;
      debugLog << "  nJet = " << nJet << endl;
      debugLog.flush();
    }
    
    auto nCscSeg = ncscSegments;
    auto nDtSeg = ndtSegments;
    auto nDtRechits = ndtRecHits;
    auto nRpc = nrpcRecHits;
    auto cscRechitsPhi = cscRechits_Phi;
    auto cscRechitsEta = cscRechits_Eta;
    auto cscRechitsX = cscRechits_X;
    auto cscRechitsY = cscRechits_Y;
    auto cscRechitsZ = cscRechits_Z;
    auto cscRechitsTpeak = cscRechits_Tpeak;
    auto cscRechitsTwire = cscRechits_Twire;
    auto cscRechitsChamber = cscRechits_Chamber;
    auto cscRechitsStation = cscRechits_Station;

    auto dtRechitCorrectX = dtRecHits_X;
    auto dtRechitCorrectY = dtRecHits_Y;
    auto dtRechitCorrectZ = dtRecHits_Z;
    auto dtRechitCorrectEta = dtRecHits_Eta;
    auto dtRechitCorrectPhi = dtRecHits_Phi;
    auto dtRechitStation = dtRecHits_Station;
    auto dtRechitWheel = dtRecHits_Wheel;
    auto dtRechitSuperLayer = dtRecHits_SuperLayer;

    auto dtSegStation = dtSegments_Station;
    auto dtSegPhi = dtSegments_Phi;
    auto dtSegEta = dtSegments_Eta;
    auto rpcPhi = rpcRecHits_Phi;
    auto rpcEta = rpcRecHits_Eta;
    auto rpcX = rpcRecHits_X;
    auto rpcY = rpcRecHits_Y;
    auto rpcZ = rpcRecHits_Z;
    
    auto rpcBx = rpcRecHits_Bx;
    auto rpcRegion = rpcRecHits_Region;
    auto rpcRing = rpcRecHits_Ring;

    auto eleCharge = Electron_charge;       // eleCharge
    auto eleEta = Electron_eta;             // eleEta
    auto elePhi = Electron_phi;             // elePhi
    auto elePt = Electron_pt;               // elePt
    auto ele_dZ = Electron_dz;              // ele_dZ
    bool ele_passCutBasedIDTight[10] = {0}; // ele_passCutBasedIDTight
    bool ele_passCutBasedIDVeto[10] = {0};  // ele_passCutBasedIDVeto

    for (int i = 0; i < nElectron; i++)
    {
      ele_passCutBasedIDTight[i] = Electron_cutBased[i] >= 4;
      ele_passCutBasedIDVeto[i] = Electron_cutBased[i] >= 1;
    }
    auto eventNum = event;                                    // eventNum
    auto fixedGridRhoFastjetAll = Rho_fixedGridRhoFastjetAll; // fixedGridRhoFastjetAll

    Float_t jetE[150] = {0}; // jetE

    for (int i = 0; i < nJet; ++i)
    {
      auto eta = Jet_eta[i];
      auto pt = Jet_pt[i];
      auto pz = pt * TMath::SinH(eta);
      auto mass = Jet_mass[i];
      jetE[i] = TMath::Sqrt(mass * mass + pt * pt + pz * pz);
    }
    auto jetEta = Jet_eta;         // jetEta
    bool jetPassIDTightLepVeto[150] = {0}; // jetPassIDLoose
    bool jetPassIDTight[150] = {0}; // jetPassIDTight
    for (int i = 0; i < nJet; ++i)
    {
      if (helper) {
        jetPassIDTight[i] = helper->jetTightLepVeto(analysisTag, false, Jet_neHEF[i], Jet_neEmEF[i], Jet_chEmEF[i], Jet_muEF[i], Jet_chHEF[i], Jet_chMultiplicity[i], Jet_neMultiplicity[i], Jet_eta[i], Jet_jetId[i]);
        jetPassIDTightLepVeto[i] = helper->jetTightLepVeto(analysisTag, true, Jet_neHEF[i], Jet_neEmEF[i], Jet_chEmEF[i], Jet_muEF[i], Jet_chHEF[i], Jet_chMultiplicity[i], Jet_neMultiplicity[i], Jet_eta[i], Jet_jetId[i]);
      } else {
        // Default values if helper is not available
        jetPassIDTight[i] = true;
        jetPassIDTightLepVeto[i] = true;
      }
      // cout<<jetPassIDTight[i]<<","<<jetPassIDTightLepVeto[i]<<", "<< Jet_neHEF[i]<<", "<<Jet_neEmEF[i]<<", "<<Jet_chEmEF[i]<<", "<<Jet_muEF[i]<<", "<<Jet_chHEF[i]<<", "<<Jet_chMultiplicity[i]<<", "<<Jet_neMultiplicity[i]<<", "<< Jet_eta[i]<<", "<< Jet_jetId[i]<<endl;
    }
    auto jetPhi = Jet_phi;          // jetPhi
    auto jetPt = Jet_pt;            // jetPt
    auto lumiNum = luminosityBlock; // lumiNum
    auto metType1Phi = MET_phi;     // metType1Phi
    auto metType1Pt = MET_pt;       // metType1Pt

    if (analysisTag == "Summer24")
    {
      metType1Phi = PFMET_phi;     // metType1Phi
      metType1Pt = PFMET_pt;       // metType1Pt
    }
    
    auto muonCharge = Muon_charge;  // muonCharge
    Float_t muonE[50] = {0};        // muonE

    for (int i = 0; i < nMuon; ++i)
    {
      auto eta = Muon_eta[i];
      auto pt = Muon_pt[i];
      auto pz = pt * TMath::SinH(eta);
      auto mass = MU_MASS;
      muonE[i] = TMath::Sqrt(mass * mass + pt * pt + pz * pz);
    }
    auto muonEta = Muon_eta;                    // muonEta
    auto muonIsLoose = Muon_looseId;            // muonIsLoose
    auto muonIsTight = Muon_tightId;            // muonIsTight
    auto muonPhi = Muon_phi;                    // muonPhi
    auto muonPt = Muon_pt;                      // muonPt
    auto muon_pfRelIso04_all = Muon_pfRelIso04_all; // muon_chargedIso
    auto muon_dZ = Muon_dz;                     // muon_dZ
    auto muon_isGlobal = Muon_isGlobal;         // muon_isGlobal
    auto nElectrons = nElectron;                // nElectrons
    auto nJets = nJet;                          // nJets
    auto nMuons = nMuon;                        // nMuons
    int nPV = PV_npvs;                          // nPV
    auto runNum = run;                          // runNum

    auto *lheComments = (string *)"123";
    auto genWeight = Generator_weight; // genWeight


    if (!isData && signalScan)
    {

      string mh_substring = lheComments->substr(lheComments->find("MH-")+3);
      int mh = stoi(mh_substring.substr(0,mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-")+3);
      int mx = stoi(mx_substring.substr(0,mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-")+6);
      int ctau = stoi(ctau_substring.substr(0,ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;

      // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
      // cout<<*lheComments<<endl;

      pair<int,int> signalPair = make_pair(mx, ctau);

      if (Files2D.count(signalPair) == 0){ //create file and tree
        //format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end()-5, thisFileName.end());
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
        Trees2D[signalPair] =  MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1,0.5,1.5);
        Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1,0.5,1.5);
        accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1,0.5,1.5);
        accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1,0.5,1.5);

        cout << "Created new output file " << thisFileName << endl;
      }
          //Fill NEvents hist
      NEvents2D[signalPair]->Fill(1.0, genWeight);


    }
    //event info
    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else
    {
      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
    }
    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    //DEBUGGING: Comment out restrictive data run filter
    if (isData && runNum < 360019)continue;
    
    //DEBUGGING: Print basic event info and rechit counts
    if (jentry < 10 || jentry % 1000 == 0) {
      debugLog << "Event " << jentry << " (eventNum=" << eventNum;
      if (isData) debugLog << ", run=" << runNum;
      debugLog << "): CSC rechits=" << ncscRechits << ", DT rechits=" << nDtRechits << endl;
    }
    bool llp_flag = false;
    if (!isData)
    {
      //for DS model
      MuonSystem->nGenParticles = 0;
      for(int i = 0;i < nGenPart; i++)
      {
        if (abs(GenPart_pdgId[i])>=999999)
        {
          // cout<<eventNum<<","<<gParticleId[i]<<","<<MuonSystem->nGLLP<<endl;
          MuonSystem->gParticleEta[MuonSystem->nGenParticles] = GenPart_eta[i];
          MuonSystem->gParticlePhi[MuonSystem->nGenParticles] = GenPart_phi[i];
          MuonSystem->gParticlePt[MuonSystem->nGenParticles] = GenPart_pt[i];
          MuonSystem->gParticleId[MuonSystem->nGenParticles] = GenPart_pdgId[i];
          float decay_x;
          float decay_y;
          float decay_z;
          bool foundDaughter = false;
          for (int j=0; j < nGenPart; j++)
          {
            if(GenPart_genPartIdxMother[j]==i)
            {
              decay_x = GenPart_vx[j];
              decay_y = GenPart_vx[j];
              decay_z = GenPart_vx[j];
              foundDaughter = true;
              break;
            }
          }
          if (!foundDaughter)cout<<"didn't find HV LLP daughter "<<eventNum<<endl;
          MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles] = sqrt(decay_x*decay_x+decay_y*decay_y);
          MuonSystem->gParticle_decay_vertex_x[MuonSystem->nGenParticles] = decay_x;
          MuonSystem->gParticle_decay_vertex_y[MuonSystem->nGenParticles] = decay_y;
          MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles] = decay_z;


          TLorentzVector particle = makeTLorentzVectorPtEtaPhiM( GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i] );

          float gParticle_decay_vertex = sqrt(pow(MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles], 2) + pow(MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles],2));

          
          MuonSystem->gParticle_ctau[MuonSystem->nGenParticles] = gParticle_decay_vertex/(particle.Beta()*particle.Gamma());
          MuonSystem->gParticle_beta[MuonSystem->nGenParticles] = particle.Beta();
            // if (abs(MuonSystem->gLLP_eta[MuonSystem->nGLLP]) < 2.4
            //   && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])>400
            //   && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
            // if (abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])< 661.0
            //   && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 800
            //    && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;
            MuonSystem->nGenParticles++;
          }

      }


      //for Twin higgs model
       MuonSystem->nGLLP = 0;
      for(int i = 0;i < nGenPart; i++)
      {
        if (abs(GenPart_pdgId[i])==9000006)
        {
          // cout<<eventNum<<","<<gParticleId[i]<<","<<MuonSystem->nGLLP<<endl;
          MuonSystem->gLLP_eta[MuonSystem->nGLLP] = GenPart_eta[i];
          MuonSystem->gLLP_phi[MuonSystem->nGLLP] = GenPart_phi[i];
          MuonSystem->gLLP_pt[MuonSystem->nGLLP] = GenPart_pt[i];
          // MuonSystem->gLLP_e[MuonSystem->nGenParticles] = GenPart_pdgId[i];
          float decay_x;
          float decay_y;
          float decay_z;
          bool foundDaughter = false;
          for (int j=0; j < nGenPart; j++)
          {
            if(GenPart_genPartIdxMother[j]==i)
            {
              decay_x = GenPart_vx[j];
              decay_y = GenPart_vy[j];
              decay_z = GenPart_vz[j];
              foundDaughter = true;
              break;
            }
          }
          if (!foundDaughter)cout<<"didn't find LLP daughter "<<eventNum<<endl;
          MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(decay_x*decay_x+decay_y*decay_y);
          MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = decay_x;
          MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = decay_y;
          MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = decay_z;


          TLorentzVector LLP = makeTLorentzVectorPtEtaPhiM( GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i] );

          float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP], 2) + pow(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP],2));

          
          MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex/(LLP.Beta()*LLP.Gamma());
          MuonSystem->gLLP_beta[MuonSystem->nGLLP] = LLP.Beta();
            if (abs(MuonSystem->gLLP_eta[MuonSystem->nGLLP]) < 2.4
              && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])>400
              && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
            if (abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])< 661.0
              && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 800
               && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;
            MuonSystem->nGLLP++;
          }

      }



    //    for(int i = 0; i < 2;i++)
    //    {
	  //      MuonSystem->gLLP_eta[MuonSystem->nGLLP] = gLLP_eta[i];
    //      MuonSystem->gLLP_phi[MuonSystem->nGLLP] = gLLP_phi[i];
    //      MuonSystem->gLLP_e[MuonSystem->nGLLP] = gLLP_e[i];
    //      MuonSystem->gLLP_pt[MuonSystem->nGLLP] = gLLP_pt[i];

    //      MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
    //      MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = gLLP_decay_vertex_x[i];
    //      MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = gLLP_decay_vertex_y[i];
    //      MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = gLLP_decay_vertex_z[i];
    //      float beta = gLLP_beta[i];
    //      float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i],2));
    //      float gamma = 1.0/sqrt(1-beta*beta);
    //      MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex/(beta * gamma);
    //      MuonSystem->gLLP_beta[MuonSystem->nGLLP] = gLLP_beta[i];


    //        if (abs(MuonSystem->gLLP_eta[i]) < 2.4
    //          && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>400
    //          && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
    //        if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 661.0
    //          && MuonSystem->gLLP_decay_vertex_r[i] < 800
    //           && MuonSystem->gLLP_decay_vertex_r[i] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;

	  //  if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 1500 && 
		// 	   MuonSystem->gLLP_decay_vertex_r[i] < 1000 &&
		// 	   MuonSystem->gLLP_decay_vertex_r[i]>100 )llp_flag = true;
    //         MuonSystem->nGLLP++;
    //    }


      MuonSystem->npu = Pileup_nTrueInt;
      if (helper) {
        MuonSystem->pileupWeight = helper->getPileupWeight(Pileup_nTrueInt);
        MuonSystem->pileupWeightUp = helper->getPileupWeightUp(Pileup_nTrueInt) / MuonSystem->pileupWeight;
        MuonSystem->pileupWeightDown = helper->getPileupWeightDown(Pileup_nTrueInt) / MuonSystem->pileupWeight;
      } else {
        // Default pileup weights if helper is not available
        MuonSystem->pileupWeight = 1.0;
        MuonSystem->pileupWeightUp = 1.0;
        MuonSystem->pileupWeightDown = 1.0;
      }

      for (unsigned int i = 0; i < 9; i++)
      {
       MuonSystem->LHEScaleWeight[i]= LHEScaleWeight[i];
      }
     
    //  if (!llp_flag)continue;
    }//end of isData

      //get NPU
      MuonSystem->npv = nPV;
      MuonSystem->rho = fixedGridRhoFastjetAll;
      MuonSystem->met = metType1Pt;
      MuonSystem->metPhi = metType1Phi;

      MuonSystem->Puppimet = PuppiMET_pt;
      MuonSystem->PuppimetPhi = PuppiMET_phi;

      if(signalScan && !isData)Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
      if(signalScan && !isData)
      {
        accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
      }
      else if (!isData)
      {
        if (MuonSystem->gLLP_csc[0] && MuonSystem->gLLP_csc[1]) accep_csccsc->Fill(1.0, genWeight*MuonSystem->pileupWeight);
        if ((MuonSystem->gLLP_dt[0] && MuonSystem->gLLP_csc[1]) || (MuonSystem->gLLP_dt[1] && MuonSystem->gLLP_csc[0])) accep_cscdt->Fill(1.0, genWeight*MuonSystem->pileupWeight);
      }

      // noise filters

      MuonSystem->Flag_goodVertices = Flag_goodVertices;
      MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
      MuonSystem->Flag_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
      MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
      MuonSystem->Flag_BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
      MuonSystem->Flag_hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;
      MuonSystem->Flag_eeBadScFilter = Flag_eeBadScFilter;
      MuonSystem->Flag_all = Flag_eeBadScFilter 
                              && Flag_hfNoisyHitsFilter 
                              && Flag_BadPFMuonDzFilter 
                              && Flag_BadPFMuonFilter 
                              && Flag_EcalDeadCellTriggerPrimitiveFilter
                              && Flag_globalSuperTightHalo2016Filter 
                              && Flag_goodVertices;
                              
      if (analysisTag == "Summer24") MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;

      // Flag_ecalBadCalibFilter for nanoAOD: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#ECal_BadCalibration_Filter_Flag
      if (analysisTag == "Summer24") MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
      else{
        MuonSystem->Flag_ecalBadCalibFilter = true;
        if (isData && runNum >= 362433 && runNum<=367144)
        {
          if (PuppiMET_pt > 100)
          {
            for(int i = 0; i < nJets; i++)
            {
              if (jetPt[i]<50) continue;
              if (!(jetEta[i] <= -0.1 && jetEta[i]>=-0.5 && jetPhi[i] <-1.8 && jetPhi[i]> -2.1)) continue;
              if (!(Jet_neEmEF[i] >0.9 || Jet_chEmEF[i]>0.9)) continue;
              if (deltaPhi(PuppiMET_phi, jetPhi[i])<2.9) continue;
              Flag_ecalBadCalibFilter = false;
            }
          }
        }
      }

      // jet veto map, following selections here: https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
      MuonSystem->jetVeto = true;
      for(int i = 0; i < nJets; i++)
      {
        if (jetPt[i] <= 15) continue;
        if (Jet_neEmEF[i] + Jet_chEmEF[i] >= 0.9) continue;
        if (!jetPassIDTight[i]) continue;
        //remove overlaps
        bool overlap = false;
        for(int j = 0; j < nMuons; j++)
        {
          if (!Muon_isPFcand[j])continue;
          if (RazorAnalyzerMerged::deltaR(jetEta[i],jetPhi[i],muonEta[j], muonPhi[j]) < 0.2) overlap = true;
        }
        if (overlap) continue;
        if (helper) {
          helper->getJetVetoMap(0,1);
          if (helper->getJetVetoMap(jetEta[i],jetPhi[i])>0.0) MuonSystem->jetVeto = false;
          if (analysisTag == "Summer24" && helper->getJetVetoFpixMap(jetEta[i],jetPhi[i])>0.0) MuonSystem->jetVeto = false;
        }
      } 
      
      MuonSystem->HLT_CSCCSC = HLT_CscCluster_Loose;
      MuonSystem->HLT_CSCDT = HLT_L1CSCShower_DTCluster50;

      //*************************************************************************
      //Start Object Selection
      //*************************************************************************
      std::vector<leptons> Leptons;

      //-------------------------------
      //Muons
      //-------------------------------
      for( int i = 0; i < nMuons; i++ )
      {
        if(!muonIsLoose[i]) continue;
        if(muonPt[i] < 25) continue;
        if(fabs(muonEta[i]) > 2.4) continue;

        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (RazorAnalyzerMerged::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;

        leptons tmpMuon;
        tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
        tmpMuon.pdgId = 13 * -1 * muonCharge[i];
        tmpMuon.dZ = muon_dZ[i];
        tmpMuon.passId = muonIsTight[i];

        tmpMuon.passLooseIso = muon_pfRelIso04_all[i]<0.25;
        tmpMuon.passTightIso = muon_pfRelIso04_all[i]<0.15;
        tmpMuon.passVTightIso = muon_pfRelIso04_all[i]<0.10;
        tmpMuon.passVVTightIso = muon_pfRelIso04_all[i]<0.05;

        tmpMuon.passVetoId = false;
        Leptons.push_back(tmpMuon);
      }

      //-------------------------------
      //Electrons
      //-------------------------------
      for( int i = 0; i < nElectrons; i++ )
      {
        if (!ele_passCutBasedIDVeto[i]) continue;
        if(elePt[i] < 35) continue;
        if(fabs(eleEta[i]) > 2.5) continue;

        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (RazorAnalyzerMerged::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;
        leptons tmpElectron;
        tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);
        tmpElectron.pdgId = 11 * -1 * eleCharge[i];
        tmpElectron.dZ = ele_dZ[i];
        tmpElectron.passId = ele_passCutBasedIDTight[i];
        Leptons.push_back(tmpElectron);
      }

      sort(Leptons.begin(), Leptons.end(), my_largest_pt);


      for ( auto &tmp : Leptons )
      {
        MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
        MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
        MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
        MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
        MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
        MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
        MuonSystem->lepTightId[MuonSystem->nLeptons] = tmp.passId;
        MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
        MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
        MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
        MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
        MuonSystem->nLeptons++;
      }


    //-----------------------------------------------
    //Select Jets
    //-----------------------------------------------

    std::vector<jets> Jets;
    float MetX_JESDown = 0;
    float MetY_JESDown = 0;

    float MetX_JESUp = 0;
    float MetY_JESUp = 0;

    for(int i = 0; i < nJets; i++)
    {
      if (fabs(jetEta[i]) >= 3.0)continue;
      if( jetPt[i] < 20 ) continue;
      if (!jetPassIDTight[i] && !jetPassIDTightLepVeto[i]) continue;
      //------------------------------------------------------------
      //exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;

      for(auto& lep : Leptons)
      {
        double thisDR = RazorAnalyzerMerged::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
        if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

      TLorentzVector thisJet = makeTLorentzVector( jetPt[i], jetEta[i], jetPhi[i], jetE[i] );
      
      jets tmpJet;
      tmpJet.jet    = thisJet;
      tmpJet.passId = jetPassIDTightLepVeto[i];

      // calculate jet energy scale uncertainty
      double unc = 0.0;
      if (helper) {
        unc = helper->getJecUnc( jetPt[i], jetEta[i], runNum ); //use run=999 as default
      }
      tmpJet.jetPtJESUp = jetPt[i]*(1+unc);
      tmpJet.jetPtJESDown = jetPt[i]*(1-unc);
      tmpJet.jetEJESUp = jetE[i]*(1+unc);
      tmpJet.jetEJESDown = jetE[i]*(1-unc);
      tmpJet.JecUnc = unc;
      TLorentzVector thisJetJESUp = makeTLorentzVector(tmpJet.jetPtJESUp, jetEta[i], jetPhi[i], tmpJet.jetEJESUp);
      TLorentzVector thisJetJESDown = makeTLorentzVector(tmpJet.jetPtJESDown, jetEta[i], jetPhi[i], tmpJet.jetEJESDown);

      MetX_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
      MetY_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());

      MetX_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
      MetY_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
      
      Jets.push_back(tmpJet);
    }

      sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

      for ( auto &tmp : Jets )
      {
        if(tmp.jet.Pt()<30)continue;

        MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
        MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
        MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
        MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
        MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;
        MuonSystem->nJets++;
      }

      TLorentzVector met = makeTLorentzVectorPtEtaPhiM(MuonSystem->Puppimet, 0, MuonSystem->PuppimetPhi, 0);

      float MetXJESUp   = met.Px() + MetX_JESUp;
      float MetYJESUp   = met.Py() + MetY_JESUp;
      MuonSystem->PuppimetJESUp    = sqrt( pow(MetXJESUp,2) + pow(MetYJESUp,2) );
      MuonSystem->PuppimetPhiJESUp    = atan(MetYJESUp/MetXJESUp);
      if  (MetXJESUp < 0.0) MuonSystem->PuppimetPhiJESUp = RazorAnalyzerMerged::deltaPhi(TMath::Pi() + MuonSystem->PuppimetPhiJESUp,0.0);

      //JES down
      float MetXJESDown   = met.Px() + MetX_JESDown;
      float MetYJESDown   = met.Py() + MetY_JESDown;
      MuonSystem->PuppimetJESDown    =  sqrt( pow(MetXJESDown,2) + pow(MetYJESDown,2) );
      MuonSystem->PuppimetPhiJESDown    = atan(MetYJESDown/MetXJESDown);
      if  (MetXJESDown < 0.0) MuonSystem->PuppimetPhiJESDown = RazorAnalyzerMerged::deltaPhi(TMath::Pi() + MuonSystem->PuppimetPhiJESDown,0.0);

      MuonSystem->nDTRechits  = 0;

      int nDTRechitsChamberMinus12 = 0;
      int nDTRechitsChamberMinus11 = 0;
      int nDTRechitsChamber10 = 0;
      int nDTRechitsChamberPlus11 = 0;
      int nDTRechitsChamberPlus12 = 0;
      int nDTRechitsChamberMinus22 = 0;
      int nDTRechitsChamberMinus21 = 0;
      int nDTRechitsChamber20 = 0;
      int nDTRechitsChamberPlus21 = 0;
      int nDTRechitsChamberPlus22 = 0;
      int nDTRechitsChamberMinus32 = 0;
      int nDTRechitsChamberMinus31 = 0;
      int nDTRechitsChamber30 = 0;
      int nDTRechitsChamberPlus31 = 0;
      int nDTRechitsChamberPlus32 = 0;
      int nDTRechitsChamberMinus42 = 0;
      int nDTRechitsChamberMinus41 = 0;
      int nDTRechitsChamber40 = 0;
      int nDTRechitsChamberPlus41 = 0;
      int nDTRechitsChamberPlus42 = 0;

      for (int i = 0; i < nDtRechits; i++) {
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus12++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) nDTRechitsChamber10++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus12++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus22++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) nDTRechitsChamber20++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus22++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus32++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) nDTRechitsChamber30++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus32++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus42++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) nDTRechitsChamber40++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus42++;
      }


      if ( nDTRechitsChamberMinus12 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus11 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber10 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus11 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus12 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus22 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus21 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber20 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus21 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus22 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus32 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus31 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber30 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus31 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus32 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus42 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus41 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber40 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus41 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus42 > 50) MuonSystem->nDtRings++;

      vector<Rechits> points;
      points.clear();

      int nCscRechitsChamberPlus11 = 0;
      int nCscRechitsChamberPlus12 = 0;
      int nCscRechitsChamberPlus13 = 0;
      int nCscRechitsChamberPlus21 = 0;
      int nCscRechitsChamberPlus22 = 0;
      int nCscRechitsChamberPlus31 = 0;
      int nCscRechitsChamberPlus32 = 0;
      int nCscRechitsChamberPlus41 = 0;
      int nCscRechitsChamberPlus42 = 0;

      int nCscRechitsChamberMinus11 = 0;
      int nCscRechitsChamberMinus12 = 0;
      int nCscRechitsChamberMinus13 = 0;
      int nCscRechitsChamberMinus21 = 0;
      int nCscRechitsChamberMinus22 = 0;
      int nCscRechitsChamberMinus31 = 0;
      int nCscRechitsChamberMinus32 = 0;
      int nCscRechitsChamberMinus41 = 0;
      int nCscRechitsChamberMinus42 = 0;

      
      for (int i = 0; i < ncscRechits; i++) {
        int layer = 0;
        Rechits point;
        point.phi = cscRechitsPhi[i];
        point.eta = cscRechitsEta[i];
        point.x = cscRechitsX[i];
        point.y = cscRechitsY[i];
        point.z = cscRechitsZ[i];
        point.t = cscRechitsTpeak[i];
        point.twire = cscRechitsTwire[i];
        point.station = cscRechitsStation[i];
        point.chamber = cscRechitsChamber[i];
        point.layer = layer;
        point.superlayer = 0;
        point.wheel = 0;
        point.clusterId = -999;
        points.push_back(point);
  
        if (cscRechitsChamber[i] == 11)  nCscRechitsChamberPlus11++;
        if (cscRechitsChamber[i] == 12)  nCscRechitsChamberPlus12++;
        if (cscRechitsChamber[i] == 13)  nCscRechitsChamberPlus13++;
        if (cscRechitsChamber[i] == 21)  nCscRechitsChamberPlus21++;
        if (cscRechitsChamber[i] == 22)  nCscRechitsChamberPlus22++;
        if (cscRechitsChamber[i] == 31)  nCscRechitsChamberPlus31++;
        if (cscRechitsChamber[i] == 32)  nCscRechitsChamberPlus32++;
        if (cscRechitsChamber[i] == 41)  nCscRechitsChamberPlus41++;
        if (cscRechitsChamber[i] == 42)  nCscRechitsChamberPlus42++;
        if (cscRechitsChamber[i] == -11)  nCscRechitsChamberMinus11++;
        if (cscRechitsChamber[i] == -12)  nCscRechitsChamberMinus12++;
        if (cscRechitsChamber[i] == -13)  nCscRechitsChamberMinus13++;
        if (cscRechitsChamber[i] == -21)  nCscRechitsChamberMinus21++;
        if (cscRechitsChamber[i] == -22)  nCscRechitsChamberMinus22++;
        if (cscRechitsChamber[i] == -31)  nCscRechitsChamberMinus31++;
        if (cscRechitsChamber[i] == -32)  nCscRechitsChamberMinus32++;
        if (cscRechitsChamber[i] == -41)  nCscRechitsChamberMinus41++;
        if (cscRechitsChamber[i] == -42)  nCscRechitsChamberMinus42++;
      }

      MuonSystem->nCscRings = 0;
      if ( nCscRechitsChamberPlus11 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus12 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus13 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus21 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus22 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus31 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus32 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus41 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus42 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus11 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus12 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus13 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus21 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus22 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus31 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus32 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus41 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus42 > 50) MuonSystem->nCscRings++;


      int min_point = 50;  
      float epsilon = 0.2;

      CACluster ds_cscRechit(min_point, epsilon, points);

      ds_cscRechit.run();
      ds_cscRechit.clusterProperties();
      //ds_cscRechit.merge_clusters();
      //ds_cscRechit.clusterProperties();
      ds_cscRechit.sort_clusters();

      // Copy cluster IDs to the fixed-size array
      for (size_t i = 0; i < ds_cscRechit.m_points.size() && i < N_MAX_CSCRECHITS; i++) 
      {
        MuonSystem->cscRechitsClusterId[i] = ds_cscRechit.m_points[i].clusterId;
      }
      
      MuonSystem->nCscRechitClusters = 0;
      MuonSystem->nCscRechitClusters_nocut = 0;

      for (auto &cluster : ds_cscRechit.clusters) {

        MuonSystem->nCscRechitClusters_nocut++;

        if (
          cluster.nCscRechitsChamberPlus11 + 
          cluster.nCscRechitsChamberPlus12 + 
          cluster.nCscRechitsChamberMinus11 + 
          cluster.nCscRechitsChamberMinus12 > 0
        ) continue;

        MuonSystem->nCscRechitClusters++;

        MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] = cluster.x;
        MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] = cluster.y;
        MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] = cluster.z;
        
        MuonSystem->cscRechitClusterTimeWeighted[MuonSystem->nCscRechitClusters] = cluster.tWeighted;
        MuonSystem->cscRechitClusterTimeSpreadWeightedAll[MuonSystem->nCscRechitClusters] = cluster.TSpreadWeightedAll;
        
        MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = cluster.tTotal;
        MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = cluster.TSpread;

        MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] = cluster.eta;
        MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = cluster.phi;

        MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = cluster.nhits;
        MuonSystem->cscRechitClusternXY[MuonSystem->nCscRechitClusters] = cluster.nXY;
        MuonSystem->cscRechitClusternZ[MuonSystem->nCscRechitClusters] = cluster.nZ;
        MuonSystem->cscRechitClusterXSpread[MuonSystem->nCscRechitClusters] = cluster.XSpread;
        MuonSystem->cscRechitClusterYSpread[MuonSystem->nCscRechitClusters] = cluster.YSpread;
        MuonSystem->cscRechitClusterZSpread[MuonSystem->nCscRechitClusters] = cluster.ZSpread;
        MuonSystem->cscRechitClusterXYSpread[MuonSystem->nCscRechitClusters] = cluster.XYSpread;
        MuonSystem->cscRechitClusterRSpread[MuonSystem->nCscRechitClusters] = cluster.RSpread;
        MuonSystem->cscRechitClusterEtaPhiSpread[MuonSystem->nCscRechitClusters] = cluster.EtaPhiSpread;
        MuonSystem->cscRechitClusterEtaSpread[MuonSystem->nCscRechitClusters] = cluster.EtaSpread;
        MuonSystem->cscRechitClusterPhiSpread[MuonSystem->nCscRechitClusters] = cluster.PhiSpread;
        MuonSystem->cscRechitClusterDeltaRSpread[MuonSystem->nCscRechitClusters] = cluster.DeltaRSpread;
        MuonSystem->cscRechitClusterMajorAxis[MuonSystem->nCscRechitClusters] = cluster.MajorAxis;
        MuonSystem->cscRechitClusterMinorAxis[MuonSystem->nCscRechitClusters] = cluster.MinorAxis;
        MuonSystem->cscRechitClusterSkewX[MuonSystem->nCscRechitClusters] = cluster.SkewX;
        MuonSystem->cscRechitClusterSkewY[MuonSystem->nCscRechitClusters] = cluster.SkewY;
        MuonSystem->cscRechitClusterSkewZ[MuonSystem->nCscRechitClusters] = cluster.SkewZ;
        MuonSystem->cscRechitClusterKurtX[MuonSystem->nCscRechitClusters] = cluster.KurtX;
        MuonSystem->cscRechitClusterKurtY[MuonSystem->nCscRechitClusters] = cluster.KurtY;
        MuonSystem->cscRechitClusterKurtZ[MuonSystem->nCscRechitClusters] = cluster.KurtZ;

        MuonSystem->cscRechitClusterNRechitChamberPlus11[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberPlus11;
        MuonSystem->cscRechitClusterNRechitChamberPlus12[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberPlus12;
        MuonSystem->cscRechitClusterNRechitChamberPlus13[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberPlus13;
        MuonSystem->cscRechitClusterNRechitChamberPlus21[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberPlus21;
        MuonSystem->cscRechitClusterNRechitChamberPlus22[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberPlus22;
        MuonSystem->cscRechitClusterNRechitChamberPlus31[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberPlus31;
        MuonSystem->cscRechitClusterNRechitChamberPlus32[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberPlus32;
        MuonSystem->cscRechitClusterNRechitChamberPlus41[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberPlus41;
        MuonSystem->cscRechitClusterNRechitChamberPlus42[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberPlus42;
        MuonSystem->cscRechitClusterNRechitChamberMinus11[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberMinus11;
        MuonSystem->cscRechitClusterNRechitChamberMinus12[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberMinus12;
        MuonSystem->cscRechitClusterNRechitChamberMinus13[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberMinus13;
        MuonSystem->cscRechitClusterNRechitChamberMinus21[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberMinus21;
        MuonSystem->cscRechitClusterNRechitChamberMinus22[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberMinus22;
        MuonSystem->cscRechitClusterNRechitChamberMinus31[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberMinus31;
        MuonSystem->cscRechitClusterNRechitChamberMinus32[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberMinus32;
        MuonSystem->cscRechitClusterNRechitChamberMinus41[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberMinus41;
        MuonSystem->cscRechitClusterNRechitChamberMinus42[MuonSystem->nCscRechitClusters] = cluster.nCscRechitsChamberMinus42;


        // MuonSystem->cscRechitClusterHMTEffPlus11[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(11, cluster.nCscRechitsChamberPlus11);
        // MuonSystem->cscRechitClusterHMTEffPlus12[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(12, cluster.nCscRechitsChamberPlus12);
        // MuonSystem->cscRechitClusterHMTEffPlus13[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(13, cluster.nCscRechitsChamberPlus13);
        // MuonSystem->cscRechitClusterHMTEffPlus21[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(21, cluster.nCscRechitsChamberPlus21);
        // MuonSystem->cscRechitClusterHMTEffPlus22[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(22, cluster.nCscRechitsChamberPlus22);
        // MuonSystem->cscRechitClusterHMTEffPlus31[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(31, cluster.nCscRechitsChamberPlus31);
        // MuonSystem->cscRechitClusterHMTEffPlus32[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(32, cluster.nCscRechitsChamberPlus32);
        // MuonSystem->cscRechitClusterHMTEffPlus41[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(41, cluster.nCscRechitsChamberPlus41);
        // MuonSystem->cscRechitClusterHMTEffPlus42[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(42, cluster.nCscRechitsChamberPlus42);
        // MuonSystem->cscRechitClusterHMTEffMinus11[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(11, cluster.nCscRechitsChamberMinus11);
        // MuonSystem->cscRechitClusterHMTEffMinus12[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(12, cluster.nCscRechitsChamberMinus12);
        // MuonSystem->cscRechitClusterHMTEffMinus13[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(13, cluster.nCscRechitsChamberMinus13);
        // MuonSystem->cscRechitClusterHMTEffMinus21[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(21, cluster.nCscRechitsChamberMinus21);
        // MuonSystem->cscRechitClusterHMTEffMinus22[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(22, cluster.nCscRechitsChamberMinus22);
        // MuonSystem->cscRechitClusterHMTEffMinus31[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(31, cluster.nCscRechitsChamberMinus31);
        // MuonSystem->cscRechitClusterHMTEffMinus32[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(32, cluster.nCscRechitsChamberMinus32);
        // MuonSystem->cscRechitClusterHMTEffMinus41[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(41, cluster.nCscRechitsChamberMinus41);
        // MuonSystem->cscRechitClusterHMTEffMinus42[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(42, cluster.nCscRechitsChamberPlus42);


        // if (helper) {
        //   float efficiencies[] = {
        //     helper->getHMTTriggerEff(13, cluster.nCscRechitsChamberPlus13),
        //     helper->getHMTTriggerEff(21, cluster.nCscRechitsChamberPlus21),
        //     helper->getHMTTriggerEff(22, cluster.nCscRechitsChamberPlus22),
        //     helper->getHMTTriggerEff(31, cluster.nCscRechitsChamberPlus31),
        //     helper->getHMTTriggerEff(32, cluster.nCscRechitsChamberPlus32),
        //     helper->getHMTTriggerEff(41, cluster.nCscRechitsChamberPlus41),
        //     helper->getHMTTriggerEff(42, cluster.nCscRechitsChamberPlus42),
        //     helper->getHMTTriggerEff(13, cluster.nCscRechitsChamberMinus13),
        //     helper->getHMTTriggerEff(21, cluster.nCscRechitsChamberMinus21),
        //     helper->getHMTTriggerEff(22, cluster.nCscRechitsChamberMinus22),
        //     helper->getHMTTriggerEff(31, cluster.nCscRechitsChamberMinus31),
        //     helper->getHMTTriggerEff(32, cluster.nCscRechitsChamberMinus32),
        //     helper->getHMTTriggerEff(41, cluster.nCscRechitsChamberMinus41),
        //     helper->getHMTTriggerEff(42, cluster.nCscRechitsChamberPlus42),
        //   };

        //   MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] = 1.0;

        //   for(const float &eff : efficiencies){
        //       MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] *= (1-eff);
        //   }
        //   MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] = 1-MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters];
        // } else {
        //   // Default efficiency if helper is not available
        //   MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] = 0.5;
        // }

        // MuonSystem->cscRechitClusterHMTEffMinus42[MuonSystem->nCscRechitClusters] = getHMTTriggerEff(11, cluster.nCscRechitsChamberPlus11);

        // cout<<cluster.nCscRechitsChamberPlus42<<", "<<helper->getHMTTriggerEff(42,   cluster.nCscRechitsChamberPlus42)<<endl;;


        MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] = cluster.maxChamber;
        MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = 1.0*cluster.maxChamberRechits/cluster.nhits;
        MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] = cluster.nChamber;
        MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = cluster.maxStation;
        MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = 1.0*cluster.maxStationRechits/cluster.nhits;

        MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = cluster.nStation10;
        MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = cluster.avgStation10;

        //Jet veto/ muon veto
        MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;


        // jet veto
        for(int i = 0; i < nJets; i++)
        {
          if (fabs(jetEta[i]>3.0)) continue;
          if (RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
            MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = jetPt[i];
            MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters]  = jetE[i];
            MuonSystem->cscRechitClusterJetVetoTightId[MuonSystem->nCscRechitClusters]  = jetPassIDTight[i];

            double unc = 0.0;
            if (helper) {
              unc = helper->getJecUnc( jetPt[i], jetEta[i], runNum ); //use run=999 as default
            }
            MuonSystem->cscRechitClusterJetVetoPtJESUp[MuonSystem->nCscRechitClusters]  = jetPt[i]*(1+unc);
            MuonSystem->cscRechitClusterJetVetoPtJESDown[MuonSystem->nCscRechitClusters]  = jetPt[i]*(1-unc);
          }

        }
        float min_deltaR = 15.;
        int index = 999;

        for(int i = 0; i < nMuons; i++)
        {
          if (fabs(muonEta[i]>3.0)) continue;
          // float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
          if (RazorAnalyzerMerged::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
            MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]  = muonPt[i];
            MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters]  = muonE[i];
            MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters]  = muon_isGlobal[i];
            MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters]  = muonIsLoose[i];
          }
        }
        if(!isData)
        {
          // match to gen level LLP
          min_deltaR = 15.;
          index = 999;
          for(int j = 0; j < MuonSystem->nGLLP;j++)
          {
            double current_delta_r = RazorAnalyzerMerged::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
            if (current_delta_r < min_deltaR)
            {
              min_deltaR = current_delta_r;
              index = j;
            }
          }
          if (min_deltaR < 0.4)MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = true;
          else MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = false;


          MuonSystem->cscRechitCluster_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters] = min_deltaR;
          MuonSystem->cscRechitCluster_match_gLLP_index[MuonSystem->nCscRechitClusters] = index;
          MuonSystem->cscRechitCluster_match_gLLP_eta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_eta[index];
          MuonSystem->cscRechitCluster_match_gLLP_phi[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_phi[index];
          MuonSystem->cscRechitCluster_match_gLLP_decay_r[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
          MuonSystem->cscRechitCluster_match_gLLP_decay_z[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
          MuonSystem->cscRechitCluster_match_gLLP_csc[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_csc[index];
          MuonSystem->cscRechitCluster_match_gLLP_dt[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_dt[index];
          MuonSystem->cscRechitCluster_match_gLLP_e[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_e[index];
        }

        //match to MB1 DT segments
        for (int i = 0; i < nDtSeg; i++) 
        {
          if (RazorAnalyzerMerged::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
          {
            MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters] ++;
            if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters] ++;
          }
        }

        //match to RPC hits in RE1/2
        for (int i = 0; i < nRpc; i++) 
        {
          float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
          if (RazorAnalyzerMerged::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
          {
            if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730) MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters] ++;
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661) MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters] ++;
          }
        }

        MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzerMerged::deltaPhi(
          MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->metPhi
        );
        MuonSystem->cscRechitClusterPuppiMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzerMerged::deltaPhi(
          MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->PuppimetPhi
        );
        MuonSystem->nCscRechitClusters++;
      }


      // DT cluster
      points.clear();

      for (int i = 0; i < nDtRechits; i++) {
        Rechits point;

        point.phi = dtRechitCorrectPhi[i];
        point.eta = dtRechitCorrectEta[i];
        point.x = dtRechitCorrectX[i];
        point.y = dtRechitCorrectY[i];
        point.z = dtRechitCorrectZ[i];
        point.t = -999.;
        point.twire = -999.;
        point.station = dtRechitStation[i];
        point.chamber = dtRechitWheel[i];
        point.superlayer = dtRechitSuperLayer[i];
        point.wheel = dtRechitWheel[i];
        point.clusterId = -999;
        points.push_back(point);
      }
      // cout<<"here"<<endl;

      int min_point_dt = 50;  //minimum number of segments to call it a cluster (was 50)
      float epsilon_dt = 0.2; //cluster radius parameter

      CACluster ds_dtRechit(min_point_dt, epsilon_dt, points);
      ds_dtRechit.run();
      ds_dtRechit.clusterProperties();
      //ds_dtRechit.merge_clusters();
      //ds_dtRechit.clusterProperties();
      ds_dtRechit.sort_clusters();


      // Copy cluster IDs to the fixed-size array
      for (size_t i = 0; i < ds_dtRechit.m_points.size() && i < N_MAX_CSCRECHITS; i++) {
        MuonSystem->dtRechitsClusterId[i] = ds_dtRechit.m_points[i].clusterId;
      }

      MuonSystem->nDtRechitClusters = 0;
      MuonSystem->nDtRechitClusters_nocut = 0;

      for (auto &cluster : ds_dtRechit.clusters) {


        bool overlap = false;
        for(int i = 0; i < MuonSystem->nCscRechitClusters; i++)
        {
          if (RazorAnalyzerMerged::deltaR(MuonSystem->cscRechitClusterEta[i],MuonSystem->cscRechitClusterPhi[i],cluster.eta, cluster.phi)<0.4) overlap = true;
        }


        if (overlap) continue;
          MuonSystem->nDtRechitClusters_nocut++;
          
	      if (cluster.nDtRechitsStation1 > 0)continue;
          MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] = cluster.x;
          MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] = cluster.y;
          MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] = cluster.z;

          if (abs(cluster.z) < 126.8) {
            MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
          }
          else if (cluster.z > 126.8 && cluster.z < 395.4) {
            MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
          }
          else if (cluster.z < -126.8 && cluster.z > -395.4) {
            MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
          }
          else if (cluster.z<0) {
            MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
          }
          else {
            MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;
          }

          MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] =cluster.eta;
          MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] =cluster.phi;
          MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = cluster.nhits;
          MuonSystem->dtRechitClusternXY[MuonSystem->nDtRechitClusters] = cluster.nXY;
          MuonSystem->dtRechitClusternZ[MuonSystem->nDtRechitClusters] = cluster.nZ;
          MuonSystem->dtRechitClusterXSpread[MuonSystem->nDtRechitClusters] = cluster.XSpread;
          MuonSystem->dtRechitClusterYSpread[MuonSystem->nDtRechitClusters] = cluster.YSpread;
          MuonSystem->dtRechitClusterZSpread[MuonSystem->nDtRechitClusters] = cluster.ZSpread;
          MuonSystem->dtRechitClusterXYSpread[MuonSystem->nDtRechitClusters] = cluster.XYSpread;
          MuonSystem->dtRechitClusterRSpread[MuonSystem->nDtRechitClusters] = cluster.RSpread;
          MuonSystem->dtRechitClusterEtaPhiSpread[MuonSystem->nDtRechitClusters] = cluster.EtaPhiSpread;
          MuonSystem->dtRechitClusterEtaSpread[MuonSystem->nDtRechitClusters] = cluster.EtaSpread;
          MuonSystem->dtRechitClusterPhiSpread[MuonSystem->nDtRechitClusters] = cluster.PhiSpread;
          MuonSystem->dtRechitClusterDeltaRSpread[MuonSystem->nDtRechitClusters] = cluster.DeltaRSpread;
          MuonSystem->dtRechitClusterMajorAxis[MuonSystem->nDtRechitClusters] = cluster.MajorAxis;
          MuonSystem->dtRechitClusterMinorAxis[MuonSystem->nDtRechitClusters] = cluster.MinorAxis;
          MuonSystem->dtRechitClusterSkewX[MuonSystem->nDtRechitClusters] = cluster.SkewX;
          MuonSystem->dtRechitClusterSkewY[MuonSystem->nDtRechitClusters] = cluster.SkewY;
          MuonSystem->dtRechitClusterSkewZ[MuonSystem->nDtRechitClusters] = cluster.SkewZ;
          MuonSystem->dtRechitClusterKurtX[MuonSystem->nDtRechitClusters] = cluster.KurtX;
          MuonSystem->dtRechitClusterKurtY[MuonSystem->nDtRechitClusters] = cluster.KurtY;
          MuonSystem->dtRechitClusterKurtZ[MuonSystem->nDtRechitClusters] = cluster.KurtZ;


          unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
          default_random_engine generator (seed);

	        // default_random_engine generator;
          uniform_real_distribution<double> distribution(0.0,1.0);
          float prob = 0.03;
          for (int i=0; i<12; ++i) {
            if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters]++;
          }
          for (int i=0; i<12; ++i) {
            if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters]++;
          }
          for (int i=0; i<12; ++i) {
            if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters]++;
          }
          for (int i=0; i<8; ++i) {
            if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters]++;
          }

          MuonSystem->dtRechitClusterNoiseHit[MuonSystem->nDtRechitClusters] = MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters] +
                                                                               MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters] +
                                                                               MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters] +
                                                                               MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters];


	        MuonSystem->dtRechitClusterNHitStation1[MuonSystem->nDtRechitClusters] = cluster.nDtRechitsStation1;
        	MuonSystem->dtRechitClusterNHitStation2[MuonSystem->nDtRechitClusters] = cluster.nDtRechitsStation2;
        	MuonSystem->dtRechitClusterNHitStation3[MuonSystem->nDtRechitClusters] = cluster.nDtRechitsStation3;
        	MuonSystem->dtRechitClusterNHitStation4[MuonSystem->nDtRechitClusters] = cluster.nDtRechitsStation4;

          MuonSystem->dtRechitClusterNHitWheel0[MuonSystem->nDtRechitClusters] = cluster.nDtRechitsWheel0;
          MuonSystem->dtRechitClusterNHitWheel1[MuonSystem->nDtRechitClusters] = cluster.nDtRechitsWheel1;
          MuonSystem->dtRechitClusterNHitWheel2[MuonSystem->nDtRechitClusters] = cluster.nDtRechitsWheel2;

        	MuonSystem->dtRechitClusterMaxChamber[MuonSystem->nDtRechitClusters] = cluster.maxChamber;
        	MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = 1.0*cluster.maxChamberRechits/cluster.nhits;
        	MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] = cluster.nChamber;
        	MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = cluster.maxStation;
        	MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = 1.0*cluster.maxStationRechits/cluster.nhits;
          MuonSystem->dtRechitClusterNStation10[MuonSystem->nDtRechitClusters] = cluster.nStation10;
          MuonSystem->dtRechitClusterAvgStation10[MuonSystem->nDtRechitClusters] = cluster.avgStation10;



          //Jet veto/ muon veto
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;


          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters]  = jetPt[i];
              MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters]  = jetE[i];
              MuonSystem->dtRechitClusterJetVetoTightId[MuonSystem->nDtRechitClusters]  = jetPassIDTight[i];

              double unc = 0.0;
              if (helper) {
                unc = helper->getJecUnc( jetPt[i], jetEta[i], runNum ); //use run=999 as default
              }
              MuonSystem->dtRechitClusterJetVetoPtJESUp[MuonSystem->nDtRechitClusters]  = jetPt[i]*(1+unc);
              MuonSystem->dtRechitClusterJetVetoPtJESDown[MuonSystem->nDtRechitClusters]  = jetPt[i]*(1-unc);
            }
          }

          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            // float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzerMerged::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]  = muonPt[i];
              MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters]  = muonE[i];
              MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters]  = muon_isGlobal[i];
              MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters]  = muonIsLoose[i];
            }
          }

          for (int i = 0; i < nDtSeg; i++) 
          {
            if (RazorAnalyzerMerged::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
              if (dtSegStation[i] == 1) MuonSystem->dtRechitClusterNSegStation1[MuonSystem->nDtRechitClusters]  +=1;
              if (dtSegStation[i] == 2) MuonSystem->dtRechitClusterNSegStation2[MuonSystem->nDtRechitClusters]  +=1;
              if (dtSegStation[i] == 3) MuonSystem->dtRechitClusterNSegStation3[MuonSystem->nDtRechitClusters]  +=1;
              if (dtSegStation[i] == 4) MuonSystem->dtRechitClusterNSegStation4[MuonSystem->nDtRechitClusters]  +=1;
            }
            if (abs(RazorAnalyzerMerged::deltaPhi(dtSegPhi[i],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]))>2) {
              if (dtSegStation[i] == 1) MuonSystem->dtRechitClusterNOppositeSegStation1[MuonSystem->nDtRechitClusters]  +=1;
              if (dtSegStation[i] == 2) MuonSystem->dtRechitClusterNOppositeSegStation2[MuonSystem->nDtRechitClusters]  +=1;
              if (dtSegStation[i] == 3) MuonSystem->dtRechitClusterNOppositeSegStation3[MuonSystem->nDtRechitClusters]  +=1;
              if (dtSegStation[i] == 4) MuonSystem->dtRechitClusterNOppositeSegStation4[MuonSystem->nDtRechitClusters]  +=1;
            }
          }



          // match to gen-level LLP
          float min_deltaR = 15.;
          int index = 999;
          if (!isData)
          {
            for(int j = 0; j < MuonSystem->nGLLP;j++)
            {
              double current_delta_r = RazorAnalyzerMerged::deltaR(cluster.eta, cluster.phi, MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
              if (current_delta_r < min_deltaR)
              {
                min_deltaR = current_delta_r;
                index = j;
              }
            }
            if (min_deltaR < 0.4)MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = true;
            else MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = false;

              MuonSystem->dtRechitCluster_match_gLLP_minDeltaR[MuonSystem->nDtRechitClusters] = min_deltaR;
              MuonSystem->dtRechitCluster_match_gLLP_index[MuonSystem->nDtRechitClusters] = index;
              MuonSystem->dtRechitCluster_match_gLLP_eta[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_eta[index];
              MuonSystem->dtRechitCluster_match_gLLP_phi[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_phi[index];
              MuonSystem->dtRechitCluster_match_gLLP_decay_r[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
              MuonSystem->dtRechitCluster_match_gLLP_decay_z[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
              MuonSystem->dtRechitCluster_match_gLLP_csc[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_csc[index];
              MuonSystem->dtRechitCluster_match_gLLP_dt[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_dt[index];
              MuonSystem->dtRechitCluster_match_gLLP_e[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_e[index];

          }


          //match to MB1 DT segments
          MuonSystem->nCscRechits = ncscRechits;

          for (int i = 0; i < nDtRechits; i++) {
            if (RazorAnalyzerMerged::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters] ++;
            }
            if (RazorAnalyzerMerged::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters] ++;
            }
            if(abs(dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])==1 && dtRechitStation[i] == 1)
            {
              if (abs(RazorAnalyzerMerged::deltaPhi(dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi()/4.0 )
              {
                if (dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1) MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters] ++;
                else MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters] ++;
              }
            }

          }

          std::vector<int> dtRechitCluster_match_rpcBx;

          //match to RPC hits with dPhi<0.5 and same wheel in DT
          for (int i = 0; i < nRpc; i++) {
            float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
            if (rpcRegion[i]!=0) continue;
            if (abs(RazorAnalyzerMerged::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
            {
              if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
              {
                dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
                MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
                if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters] ++;
              }
            }
            if(RazorAnalyzerMerged::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
            {
              if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters] ++;
            }
          }

          int max_occurence = 0;
          int max_bx = -999;

          for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
          {
            int counter = 0;
            for(unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j ++)
            {
              if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l]) counter++;
            }
            if (counter>max_occurence)
            {
              max_occurence = counter;
              max_bx = dtRechitCluster_match_rpcBx[l];
            }
          }
          MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;

          MuonSystem->dtRechitClusterMet_dPhi[MuonSystem->nDtRechitClusters] =  RazorAnalyzerMerged::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhi);
          MuonSystem->dtRechitClusterPuppiMet_dPhi[MuonSystem->nDtRechitClusters] =  RazorAnalyzerMerged::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->PuppimetPhi);

          MuonSystem->nDtRechitClusters++;
        }

      if (MuonSystem->nCscRechitClusters_nocut == 0)continue;
      if (MuonSystem->nCscRechitClusters == 0)continue;

      if (MuonSystem->nDtRechitClusters_nocut + MuonSystem->nCscRechitClusters_nocut < 2) continue;
      if (MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters < 2) continue;
      
      debugLog << "Event " << eventNum << ": CSC clusters (nocut/cut): " 
               << MuonSystem->nCscRechitClusters_nocut << "/" << MuonSystem->nCscRechitClusters
               << ", DT clusters (nocut/cut): " 
               << MuonSystem->nDtRechitClusters_nocut << "/" << MuonSystem->nDtRechitClusters << endl;

      if(!isData && signalScan)
      {
        pair<int,int> smsPair = make_pair(MuonSystem->mX, MuonSystem->ctau);
        Trees2D[smsPair]->Fill();
      }
      else
      {
        debugLog << "Filling tree for event " << eventNum << endl;
        MuonSystem->tree_->Fill();
      }


    }
      if(!isData && signalScan)
      {
        for(auto &filePtr : Files2D)
        {
          cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
          debugLog << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
          filePtr.second->cd();
          Trees2D[filePtr.first]->Write();
          NEvents2D[filePtr.first]->Write("NEvents");
          Total2D[filePtr.first]->Write("Total");
          accep2D[filePtr.first]->Write("acceptance");
          accep_met2D[filePtr.first]->Write("acceptance_met");
          filePtr.second->Close();
        }
        debugLog << "=== Signal scan analysis completed ===" << endl;
        debugLog.close();
      }
      else if (!isData)
      {
        cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
        cout << "Writing output trees..." << endl;
        debugLog << "Filled Total of " << NEvents->GetBinContent(1) << " Events" << endl;
        debugLog << "Writing output trees..." << endl;
        outFile->cd();
        MuonSystem->tree_->Write();
        NEvents->Write();
        accep->Write("acceptance");
        accep_csccsc->Write("acceptance_csccsc");
        accep_cscdt->Write("acceptance_cscdt");
        accep_met->Write("acceptance_met");
        outFile->Close();
        debugLog << "=== MC analysis completed ===" << endl;
        debugLog.close();
      }
      else
      {
        cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
        cout << "Writing output trees..." << endl;
        outFile->cd();
        MuonSystem->tree_->Write();
        Nmet200->Write();
        NmetFilter->Write();
        Nlep0->Write();
        Njet1->Write();
        NcosmicVeto->Write();
        NEvents->Write();
        // outFile->Write();
        outFile->Close();
      }
      
      //DEBUGGING: Close debug log with summary
      debugLog << "=== Data analysis completed ===" << endl;
      debugLog << "Total events processed: " << fChain->GetEntries() << endl;
      debugLog << "Events written to tree: " << NEvents->GetBinContent(1) << endl;
      debugLog.close();
}
