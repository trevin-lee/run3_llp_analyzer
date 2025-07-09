#include "RazorHelper.h"
#include <chrono>
#include <fstream>

// Constructor
RazorHelper::RazorHelper(std::string tag_, bool isData_):
        tag(tag_), isData(isData_){
    
    // Create debug log file for RazorHelper
    std::ofstream razorLog("debug_razorhelper.log");
    razorLog << "=== RazorHelper Debug Log Started ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "RazorHelper initializing with tag " << tag << std::endl;
    razorLog << "RazorHelper initializing with tag " << tag << std::endl;
    razorLog.flush();

    // check that CMSSW is set up
    razorLog << "Loading CMSSW path..." << std::endl;
    razorLog.flush();
    auto cmssw_start = std::chrono::high_resolution_clock::now();
    loadCMSSWPath();
    auto cmssw_end = std::chrono::high_resolution_clock::now();
    auto cmssw_duration = std::chrono::duration_cast<std::chrono::milliseconds>(cmssw_end - cmssw_start);
    razorLog << "CMSSW path loading took " << cmssw_duration.count() << " ms" << std::endl;
    razorLog << "CMSSW path: " << cmsswPath << std::endl;
    razorLog.flush();
    
    if (cmsswPath == "") {
        razorLog << "CMSSW path is empty, loading null tag" << std::endl;
        razorLog.flush();
        loadTag_Null();
        razorLog.close();
        return;
    }

    // tag for Run 3
    razorLog << "Loading tag-specific configurations..." << std::endl;
    razorLog.flush();
    auto tag_start = std::chrono::high_resolution_clock::now();
    
    if (tag == "Summer22") {
        razorLog << "Loading Summer22 tag..." << std::endl;
        razorLog.flush();
        loadTag_Summer22();
    }
    else if (tag == "Summer22EE") {
        razorLog << "Loading Summer22EE tag..." << std::endl;
        razorLog.flush();
        loadTag_Summer22EE();
    }
    else if (tag == "Summer23") {
        razorLog << "Loading Summer23 tag..." << std::endl;
        razorLog.flush();
        loadTag_Summer23();
    }
    else if (tag == "Summer23BPix") {
        razorLog << "Loading Summer23BPix tag..." << std::endl;
        razorLog.flush();
        loadTag_Summer23BPix();
    }
    else if (tag == "Summer24") {
        razorLog << "Loading Summer24 tag..." << std::endl;
        razorLog.flush();
        loadTag_Summer24();
    }
    // tag not found
    else {
        std::cout << "Error in RazorHelper::RazorHelper : specified tag " << tag << " is not supported!" << std::endl;
        razorLog << "Error: specified tag " << tag << " is not supported!" << std::endl;
        razorLog.flush();
        loadTag_Null();
        razorLog.close();
        return;
    }
    
    auto tag_end = std::chrono::high_resolution_clock::now();
    auto tag_duration = std::chrono::duration_cast<std::chrono::milliseconds>(tag_end - tag_start);
    razorLog << "Tag loading took " << tag_duration.count() << " ms" << std::endl;
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - start_time);
    razorLog << "Total RazorHelper initialization took " << total_duration.count() << " ms" << std::endl;
    razorLog << "=== RazorHelper initialization completed ===" << std::endl;
    razorLog.close();
}

// Destructor -- close all of the open TFile objects
RazorHelper::~RazorHelper() {
    if (pileupWeightFile) {
        pileupWeightFile->Close();
        delete pileupWeightFile;
    }

}

// Retrieves CMSSW_BASE and stores in variable cmsswPath
void RazorHelper::loadCMSSWPath() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading CMSSW path..." << std::endl;
    razorLog.flush();
    
    char* cmsswPathChar = getenv("CMSSW_BASE");
    if (cmsswPathChar == NULL) {
        std::cout << "Warning in RazorHelper::loadCMSSWPath : CMSSW_BASE not detected." << std::endl;
        razorLog << "Warning: CMSSW_BASE not detected." << std::endl;
        razorLog.flush();
        cmsswPath = "";
        razorLog.close();
        return;
    }
    cmsswPath = std::string(cmsswPathChar);
    razorLog << "CMSSW_BASE found: " << cmsswPath << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadTag_Null() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading null tag (setting everything to 0)..." << std::endl;
    razorLog.flush();
    
    std::cout << "Warning: initializing all RazorHelper files and histograms to 0" << std::endl;

    // pileup weights
    pileupWeightFile = 0;
    pileupWeightHist = 0;
    pileupWeightSysUpHist = 0;
    pileupWeightSysDownHist = 0;

    razorLog << "Null tag loading completed" << std::endl;
    razorLog.flush();
    razorLog.close();
}


void RazorHelper::loadHMTEfficiency2223() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading HMT L1 efficiency histograms for 2022-2023..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading HMT L1 efficiency histograms" << std::endl;
    
    razorLog << "Opening file: L1_efficiencies_2022_2023_032625-Hists-TEff.root" << std::endl;
    razorLog.flush();
    HMTEffFile = TFile::Open("L1_efficiencies_2022_2023_032625-Hists-TEff.root");
    
    if (!HMTEffFile) {
        razorLog << "ERROR: Failed to open HMT efficiency file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "HMT efficiency file opened successfully" << std::endl;
    razorLog.flush();

    razorLog << "Loading individual ME histograms..." << std::endl;
    razorLog.flush();
    
    HMTEffHist[11] = (TEfficiency*)HMTEffFile->Get("ME11");
    razorLog << "Loaded ME11: " << (HMTEffHist[11] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[12] = (TEfficiency*)HMTEffFile->Get("ME12");
    razorLog << "Loaded ME12: " << (HMTEffHist[12] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[13] = (TEfficiency*)HMTEffFile->Get("ME13");
    razorLog << "Loaded ME13: " << (HMTEffHist[13] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[21] = (TEfficiency*)HMTEffFile->Get("ME21");
    razorLog << "Loaded ME21: " << (HMTEffHist[21] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[22] = (TEfficiency*)HMTEffFile->Get("ME22");
    razorLog << "Loaded ME22: " << (HMTEffHist[22] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[31] = (TEfficiency*)HMTEffFile->Get("ME31");
    razorLog << "Loaded ME31: " << (HMTEffHist[31] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[32] = (TEfficiency*)HMTEffFile->Get("ME32");
    razorLog << "Loaded ME32: " << (HMTEffHist[32] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[41] = (TEfficiency*)HMTEffFile->Get("ME41");
    razorLog << "Loaded ME41: " << (HMTEffHist[41] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[42] = (TEfficiency*)HMTEffFile->Get("ME42");
    razorLog << "Loaded ME42: " << (HMTEffHist[42] ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "HMT efficiency loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadHMTEfficiency24() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading HMT L1 efficiency histograms for 2024..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading HMT L1 efficiency histograms" << std::endl;
    
    razorLog << "Opening file: HMT_Efficiencies_2024.root" << std::endl;
    razorLog.flush();
    HMTEffFile = TFile::Open("HMT_Efficiencies_2024.root");
    
    if (!HMTEffFile) {
        razorLog << "ERROR: Failed to open HMT efficiency file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "HMT efficiency file opened successfully" << std::endl;
    razorLog.flush();

    razorLog << "Loading individual ME histograms..." << std::endl;
    razorLog.flush();
    
    HMTEffHist[11] = (TEfficiency*)HMTEffFile->Get("ME11");
    razorLog << "Loaded ME11: " << (HMTEffHist[11] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[12] = (TEfficiency*)HMTEffFile->Get("ME12");
    razorLog << "Loaded ME12: " << (HMTEffHist[12] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[13] = (TEfficiency*)HMTEffFile->Get("ME13");
    razorLog << "Loaded ME13: " << (HMTEffHist[13] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[21] = (TEfficiency*)HMTEffFile->Get("ME21");
    razorLog << "Loaded ME21: " << (HMTEffHist[21] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[22] = (TEfficiency*)HMTEffFile->Get("ME22");
    razorLog << "Loaded ME22: " << (HMTEffHist[22] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[31] = (TEfficiency*)HMTEffFile->Get("ME31");
    razorLog << "Loaded ME31: " << (HMTEffHist[31] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[32] = (TEfficiency*)HMTEffFile->Get("ME32");
    razorLog << "Loaded ME32: " << (HMTEffHist[32] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[41] = (TEfficiency*)HMTEffFile->Get("ME41");
    razorLog << "Loaded ME41: " << (HMTEffHist[41] ? "SUCCESS" : "FAILED") << std::endl;
    HMTEffHist[42] = (TEfficiency*)HMTEffFile->Get("ME42");
    razorLog << "Loaded ME42: " << (HMTEffHist[42] ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "HMT efficiency loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

////////////////////////////////////////////////
//  Summer 22
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer22() {
  std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
  razorLog << "Loading Summer22 tag components..." << std::endl;
  razorLog.flush();
  auto start_time = std::chrono::high_resolution_clock::now();
  
  auto pileup_start = std::chrono::high_resolution_clock::now();
  loadPileup_Summer22();
  auto pileup_end = std::chrono::high_resolution_clock::now();
  auto pileup_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pileup_end - pileup_start);
  razorLog << "Pileup loading took " << pileup_duration.count() << " ms" << std::endl;
  
  auto hmt_start = std::chrono::high_resolution_clock::now();
  loadHMTEfficiency2223();
  auto hmt_end = std::chrono::high_resolution_clock::now();
  auto hmt_duration = std::chrono::duration_cast<std::chrono::milliseconds>(hmt_end - hmt_start);
  razorLog << "HMT efficiency loading took " << hmt_duration.count() << " ms" << std::endl;
  
  auto jetveto_start = std::chrono::high_resolution_clock::now();
  loadJetVeto_Summer22();
  auto jetveto_end = std::chrono::high_resolution_clock::now();
  auto jetveto_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jetveto_end - jetveto_start);
  razorLog << "Jet veto loading took " << jetveto_duration.count() << " ms" << std::endl;
  
  auto jec_start = std::chrono::high_resolution_clock::now();
  loadJECs();
  auto jec_end = std::chrono::high_resolution_clock::now();
  auto jec_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jec_end - jec_start);
  razorLog << "JEC loading took " << jec_duration.count() << " ms" << std::endl;
  
  auto met_start = std::chrono::high_resolution_clock::now();
  loadMetTrigger_Summer22();
  auto met_end = std::chrono::high_resolution_clock::now();
  auto met_duration = std::chrono::duration_cast<std::chrono::milliseconds>(met_end - met_start);
  razorLog << "MET trigger loading took " << met_duration.count() << " ms" << std::endl;
  
  auto total_end = std::chrono::high_resolution_clock::now();
  auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - start_time);
  razorLog << "Total Summer22 loading took " << total_duration.count() << " ms" << std::endl;
  razorLog.flush();
  razorLog.close();
}


void RazorHelper::loadMetTrigger_Summer22() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading MET trigger histograms for Summer22..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
    
    razorLog << "Opening file: METTriggerEff_Summer22.root" << std::endl;
    razorLog.flush();
    MetTriggerFile = TFile::Open("METTriggerEff_Summer22.root");
    
    if (!MetTriggerFile) {
        razorLog << "ERROR: Failed to open MET trigger file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "MET trigger file opened successfully" << std::endl;
    razorLog.flush();
    
    MetTriggerHist = (TH1F*)MetTriggerFile->Get("eff_exp");
    razorLog << "Loaded eff_exp: " << (MetTriggerHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysUpHist = (TH1F*)MetTriggerFile->Get("eff_high");
    razorLog << "Loaded eff_high: " << (MetTriggerSysUpHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysDownHist = (TH1F*)MetTriggerFile->Get("eff_low");
    razorLog << "Loaded eff_low: " << (MetTriggerSysDownHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "MET trigger loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadPileup_Summer22() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading pileup weights for Summer22..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    
    razorLog << "Opening file: PileupReweight_Summer22.root" << std::endl;
    razorLog.flush();
    pileupWeightFile = TFile::Open("PileupReweight_Summer22.root");
    
    if (!pileupWeightFile) {
        razorLog << "ERROR: Failed to open pileup weight file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "Pileup weight file opened successfully" << std::endl;
    razorLog.flush();
    
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    razorLog << "Loaded npu_nominal: " << (pileupWeightHist ? "SUCCESS" : "FAILED") << std::endl;
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    razorLog << "Loaded npu_up: " << (pileupWeightSysUpHist ? "SUCCESS" : "FAILED") << std::endl;
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
    razorLog << "Loaded npu_down: " << (pileupWeightSysDownHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "Pileup loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadJetVeto_Summer22() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading jet veto map for Summer22..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading jet veto map histograms" << std::endl;
    
    razorLog << "Opening file: Summer22_23Sep2023_RunCD_v1.root" << std::endl;
    razorLog.flush();
    JetVetoFile = TFile::Open("Summer22_23Sep2023_RunCD_v1.root");
    
    if (!JetVetoFile) {
        cout<<"Jet Veto File Not Found"<<endl;
        razorLog << "ERROR: Jet Veto File Not Found!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "Jet veto file opened successfully" << std::endl;
    razorLog.flush();
    
    JetVetoHist = (TH2D*)JetVetoFile->Get("jetvetomap");
    razorLog << "Loaded jetvetomap: " << (JetVetoHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "Jet veto loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}
////////////////////////////////////////////////
//  Summer 22EE
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer22EE() {
  std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
  razorLog << "Loading Summer22EE tag components..." << std::endl;
  razorLog.flush();
  auto start_time = std::chrono::high_resolution_clock::now();
  
  auto pileup_start = std::chrono::high_resolution_clock::now();
  loadPileup_Summer22EE();
  auto pileup_end = std::chrono::high_resolution_clock::now();
  auto pileup_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pileup_end - pileup_start);
  razorLog << "Pileup loading took " << pileup_duration.count() << " ms" << std::endl;
  
  auto hmt_start = std::chrono::high_resolution_clock::now();
  loadHMTEfficiency2223();
  auto hmt_end = std::chrono::high_resolution_clock::now();
  auto hmt_duration = std::chrono::duration_cast<std::chrono::milliseconds>(hmt_end - hmt_start);
  razorLog << "HMT efficiency loading took " << hmt_duration.count() << " ms" << std::endl;
  
  auto jetveto_start = std::chrono::high_resolution_clock::now();
  loadJetVeto_Summer22EE();
  auto jetveto_end = std::chrono::high_resolution_clock::now();
  auto jetveto_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jetveto_end - jetveto_start);
  razorLog << "Jet veto loading took " << jetveto_duration.count() << " ms" << std::endl;
  
  auto jec_start = std::chrono::high_resolution_clock::now();
  loadJECs();
  auto jec_end = std::chrono::high_resolution_clock::now();
  auto jec_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jec_end - jec_start);
  razorLog << "JEC loading took " << jec_duration.count() << " ms" << std::endl;
  
  auto met_start = std::chrono::high_resolution_clock::now();
  loadMetTrigger_Summer22EE();
  auto met_end = std::chrono::high_resolution_clock::now();
  auto met_duration = std::chrono::duration_cast<std::chrono::milliseconds>(met_end - met_start);
  razorLog << "MET trigger loading took " << met_duration.count() << " ms" << std::endl;
  
  auto total_end = std::chrono::high_resolution_clock::now();
  auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - start_time);
  razorLog << "Total Summer22EE loading took " << total_duration.count() << " ms" << std::endl;
  razorLog.flush();
  razorLog.close();
}

void RazorHelper::loadMetTrigger_Summer22EE() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading MET trigger histograms for Summer22EE..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
    
    razorLog << "Opening file: METTriggerEff_Summer22EE.root" << std::endl;
    razorLog.flush();
    MetTriggerFile = TFile::Open("METTriggerEff_Summer22EE.root");
    
    if (!MetTriggerFile) {
        razorLog << "ERROR: Failed to open MET trigger file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "MET trigger file opened successfully" << std::endl;
    razorLog.flush();
    
    MetTriggerHist = (TH1F*)MetTriggerFile->Get("eff_exp");
    razorLog << "Loaded eff_exp: " << (MetTriggerHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysUpHist = (TH1F*)MetTriggerFile->Get("eff_high");
    razorLog << "Loaded eff_high: " << (MetTriggerSysUpHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysDownHist = (TH1F*)MetTriggerFile->Get("eff_low");
    razorLog << "Loaded eff_low: " << (MetTriggerSysDownHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "MET trigger loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}
void RazorHelper::loadPileup_Summer22EE() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("PileupReweight_Summer22EE.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
}

void RazorHelper::loadJetVeto_Summer22EE() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    JetVetoFile = TFile::Open("Summer22EE_23Sep2023_RunEFG_v1.root");
    JetVetoHist = (TH2D*)JetVetoFile->Get("jetvetomap");
}
////////////////////////////////////////////////
//  Summer 23
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer23() {
  std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
  razorLog << "Loading Summer23 tag components..." << std::endl;
  razorLog.flush();
  auto start_time = std::chrono::high_resolution_clock::now();
  
  auto pileup_start = std::chrono::high_resolution_clock::now();
  loadPileup_Summer23();
  auto pileup_end = std::chrono::high_resolution_clock::now();
  auto pileup_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pileup_end - pileup_start);
  razorLog << "Pileup loading took " << pileup_duration.count() << " ms" << std::endl;
  
  auto hmt_start = std::chrono::high_resolution_clock::now();
  loadHMTEfficiency2223();
  auto hmt_end = std::chrono::high_resolution_clock::now();
  auto hmt_duration = std::chrono::duration_cast<std::chrono::milliseconds>(hmt_end - hmt_start);
  razorLog << "HMT efficiency loading took " << hmt_duration.count() << " ms" << std::endl;
  
  auto jetveto_start = std::chrono::high_resolution_clock::now();
  loadJetVeto_Summer23();
  auto jetveto_end = std::chrono::high_resolution_clock::now();
  auto jetveto_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jetveto_end - jetveto_start);
  razorLog << "Jet veto loading took " << jetveto_duration.count() << " ms" << std::endl;
  
  auto jec_start = std::chrono::high_resolution_clock::now();
  loadJECs();
  auto jec_end = std::chrono::high_resolution_clock::now();
  auto jec_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jec_end - jec_start);
  razorLog << "JEC loading took " << jec_duration.count() << " ms" << std::endl;
  
  auto met_start = std::chrono::high_resolution_clock::now();
  loadMetTrigger_Summer23();
  auto met_end = std::chrono::high_resolution_clock::now();
  auto met_duration = std::chrono::duration_cast<std::chrono::milliseconds>(met_end - met_start);
  razorLog << "MET trigger loading took " << met_duration.count() << " ms" << std::endl;
  
  auto total_end = std::chrono::high_resolution_clock::now();
  auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - start_time);
  razorLog << "Total Summer23 loading took " << total_duration.count() << " ms" << std::endl;
  razorLog.flush();
  razorLog.close();
}

void RazorHelper::loadMetTrigger_Summer23() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading MET trigger histograms for Summer23..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
    
    razorLog << "Opening file: METTriggerEff_Summer23.root" << std::endl;
    razorLog.flush();
    MetTriggerFile = TFile::Open("METTriggerEff_Summer23.root");
    
    if (!MetTriggerFile) {
        razorLog << "ERROR: Failed to open MET trigger file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "MET trigger file opened successfully" << std::endl;
    razorLog.flush();
    
    MetTriggerHist = (TH1F*)MetTriggerFile->Get("eff_exp");
    razorLog << "Loaded eff_exp: " << (MetTriggerHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysUpHist = (TH1F*)MetTriggerFile->Get("eff_high");
    razorLog << "Loaded eff_high: " << (MetTriggerSysUpHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysDownHist = (TH1F*)MetTriggerFile->Get("eff_low");
    razorLog << "Loaded eff_low: " << (MetTriggerSysDownHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "MET trigger loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadPileup_Summer23() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading pileup weights for Summer23..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    
    razorLog << "Opening file: PileupReweight_Summer23.root" << std::endl;
    razorLog.flush();
    pileupWeightFile = TFile::Open("PileupReweight_Summer23.root");
    
    if (!pileupWeightFile) {
        razorLog << "ERROR: Failed to open pileup weight file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "Pileup weight file opened successfully" << std::endl;
    razorLog.flush();
    
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    razorLog << "Loaded npu_nominal: " << (pileupWeightHist ? "SUCCESS" : "FAILED") << std::endl;
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    razorLog << "Loaded npu_up: " << (pileupWeightSysUpHist ? "SUCCESS" : "FAILED") << std::endl;
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
    razorLog << "Loaded npu_down: " << (pileupWeightSysDownHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "Pileup loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadJetVeto_Summer23() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading jet veto map for Summer23..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    
    razorLog << "Opening file: Summer23Prompt23_RunC_v1.root" << std::endl;
    razorLog.flush();
    JetVetoFile = TFile::Open("Summer23Prompt23_RunC_v1.root");
    
    if (!JetVetoFile) {
        razorLog << "ERROR: Jet Veto File Not Found!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "Jet veto file opened successfully" << std::endl;
    razorLog.flush();
    
    JetVetoHist = (TH2D*)JetVetoFile->Get("jetvetomap");
    razorLog << "Loaded jetvetomap: " << (JetVetoHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "Jet veto loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}
////////////////////////////////////////////////
//  Summer 23BPix
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer23BPix() {
  std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
  razorLog << "Loading Summer23BPix tag components..." << std::endl;
  razorLog.flush();
  auto start_time = std::chrono::high_resolution_clock::now();
  
  auto pileup_start = std::chrono::high_resolution_clock::now();
  loadPileup_Summer23BPix();
  auto pileup_end = std::chrono::high_resolution_clock::now();
  auto pileup_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pileup_end - pileup_start);
  razorLog << "Pileup loading took " << pileup_duration.count() << " ms" << std::endl;
  
  auto hmt_start = std::chrono::high_resolution_clock::now();
  loadHMTEfficiency2223();
  auto hmt_end = std::chrono::high_resolution_clock::now();
  auto hmt_duration = std::chrono::duration_cast<std::chrono::milliseconds>(hmt_end - hmt_start);
  razorLog << "HMT efficiency loading took " << hmt_duration.count() << " ms" << std::endl;
  
  auto jetveto_start = std::chrono::high_resolution_clock::now();
  loadJetVeto_Summer23BPix();
  auto jetveto_end = std::chrono::high_resolution_clock::now();
  auto jetveto_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jetveto_end - jetveto_start);
  razorLog << "Jet veto loading took " << jetveto_duration.count() << " ms" << std::endl;
  
  auto jec_start = std::chrono::high_resolution_clock::now();
  loadJECs();
  auto jec_end = std::chrono::high_resolution_clock::now();
  auto jec_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jec_end - jec_start);
  razorLog << "JEC loading took " << jec_duration.count() << " ms" << std::endl;
  
  auto met_start = std::chrono::high_resolution_clock::now();
  loadMetTrigger_Summer23BPix();
  auto met_end = std::chrono::high_resolution_clock::now();
  auto met_duration = std::chrono::duration_cast<std::chrono::milliseconds>(met_end - met_start);
  razorLog << "MET trigger loading took " << met_duration.count() << " ms" << std::endl;
  
  auto total_end = std::chrono::high_resolution_clock::now();
  auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - start_time);
  razorLog << "Total Summer23BPix loading took " << total_duration.count() << " ms" << std::endl;
  razorLog.flush();
  razorLog.close();
}

void RazorHelper::loadMetTrigger_Summer23BPix() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading MET trigger histograms for Summer23BPix..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
    
    razorLog << "Opening file: METTriggerEff_Summer23BPix.root" << std::endl;
    razorLog.flush();
    MetTriggerFile = TFile::Open("METTriggerEff_Summer23BPix.root");
    
    if (!MetTriggerFile) {
        razorLog << "ERROR: Failed to open MET trigger file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "MET trigger file opened successfully" << std::endl;
    razorLog.flush();
    
    MetTriggerHist = (TH1F*)MetTriggerFile->Get("eff_exp");
    razorLog << "Loaded eff_exp: " << (MetTriggerHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysUpHist = (TH1F*)MetTriggerFile->Get("eff_high");
    razorLog << "Loaded eff_high: " << (MetTriggerSysUpHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysDownHist = (TH1F*)MetTriggerFile->Get("eff_low");
    razorLog << "Loaded eff_low: " << (MetTriggerSysDownHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "MET trigger loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadPileup_Summer23BPix() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading pileup weights for Summer23BPix..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    
    razorLog << "Opening file: PileupReweight_Summer23BPix.root" << std::endl;
    razorLog.flush();
    pileupWeightFile = TFile::Open("PileupReweight_Summer23BPix.root");
    
    if (!pileupWeightFile) {
        razorLog << "ERROR: Failed to open pileup weight file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "Pileup weight file opened successfully" << std::endl;
    razorLog.flush();
    
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    razorLog << "Loaded npu_nominal: " << (pileupWeightHist ? "SUCCESS" : "FAILED") << std::endl;
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    razorLog << "Loaded npu_up: " << (pileupWeightSysUpHist ? "SUCCESS" : "FAILED") << std::endl;
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
    razorLog << "Loaded npu_down: " << (pileupWeightSysDownHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "Pileup loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadJetVeto_Summer23BPix() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading jet veto map for Summer23BPix..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    
    razorLog << "Opening file: Summer23BPixPrompt23_RunD_v1.root" << std::endl;
    razorLog.flush();
    JetVetoFile = TFile::Open("Summer23BPixPrompt23_RunD_v1.root");
    
    if (!JetVetoFile) {
        razorLog << "ERROR: Jet Veto File Not Found!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "Jet veto file opened successfully" << std::endl;
    razorLog.flush();
    
    JetVetoHist = (TH2D*)JetVetoFile->Get("jetvetomap");
    razorLog << "Loaded jetvetomap: " << (JetVetoHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "Jet veto loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}
////////////////////////////////////////////////
//  Summer 24
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer24() {
  std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
  razorLog << "Loading Summer24 tag components..." << std::endl;
  razorLog.flush();
  auto start_time = std::chrono::high_resolution_clock::now();
  
  auto pileup_start = std::chrono::high_resolution_clock::now();
  loadPileup_Summer24();
  auto pileup_end = std::chrono::high_resolution_clock::now();
  auto pileup_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pileup_end - pileup_start);
  razorLog << "Pileup loading took " << pileup_duration.count() << " ms" << std::endl;
  
  auto hmt_start = std::chrono::high_resolution_clock::now();
  loadHMTEfficiency24();
  auto hmt_end = std::chrono::high_resolution_clock::now();
  auto hmt_duration = std::chrono::duration_cast<std::chrono::milliseconds>(hmt_end - hmt_start);
  razorLog << "HMT efficiency loading took " << hmt_duration.count() << " ms" << std::endl;
  
  auto jetveto_start = std::chrono::high_resolution_clock::now();
  loadJetVeto_Summer24();
  auto jetveto_end = std::chrono::high_resolution_clock::now();
  auto jetveto_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jetveto_end - jetveto_start);
  razorLog << "Jet veto loading took " << jetveto_duration.count() << " ms" << std::endl;
  
  auto jec_start = std::chrono::high_resolution_clock::now();
  loadJECs();
  auto jec_end = std::chrono::high_resolution_clock::now();
  auto jec_duration = std::chrono::duration_cast<std::chrono::milliseconds>(jec_end - jec_start);
  razorLog << "JEC loading took " << jec_duration.count() << " ms" << std::endl;
  
  auto met_start = std::chrono::high_resolution_clock::now();
  loadMetTrigger_Summer24();
  auto met_end = std::chrono::high_resolution_clock::now();
  auto met_duration = std::chrono::duration_cast<std::chrono::milliseconds>(met_end - met_start);
  razorLog << "MET trigger loading took " << met_duration.count() << " ms" << std::endl;
  
  auto total_end = std::chrono::high_resolution_clock::now();
  auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - start_time);
  razorLog << "Total Summer24 loading took " << total_duration.count() << " ms" << std::endl;
  razorLog.flush();
  razorLog.close();
}

void RazorHelper::loadMetTrigger_Summer24() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading MET trigger histograms for Summer24..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
    
    razorLog << "Opening file: METTriggerEff_Summer24.root" << std::endl;
    razorLog.flush();
    MetTriggerFile = TFile::Open("METTriggerEff_Summer24.root");
    
    if (!MetTriggerFile) {
        razorLog << "ERROR: Failed to open MET trigger file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "MET trigger file opened successfully" << std::endl;
    razorLog.flush();
    
    MetTriggerHist = (TH1F*)MetTriggerFile->Get("eff_exp");
    razorLog << "Loaded eff_exp: " << (MetTriggerHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysUpHist = (TH1F*)MetTriggerFile->Get("eff_high");
    razorLog << "Loaded eff_high: " << (MetTriggerSysUpHist ? "SUCCESS" : "FAILED") << std::endl;
    MetTriggerSysDownHist = (TH1F*)MetTriggerFile->Get("eff_low");
    razorLog << "Loaded eff_low: " << (MetTriggerSysDownHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "MET trigger loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadPileup_Summer24() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading pileup weights for Summer24..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    
    razorLog << "Opening file: PileupReweight_Summer24.root" << std::endl;
    razorLog.flush();
    pileupWeightFile = TFile::Open("PileupReweight_Summer24.root");
    
    if (!pileupWeightFile) {
        razorLog << "ERROR: Failed to open pileup weight file!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "Pileup weight file opened successfully" << std::endl;
    razorLog.flush();
    
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    razorLog << "Loaded npu_nominal: " << (pileupWeightHist ? "SUCCESS" : "FAILED") << std::endl;
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    razorLog << "Loaded npu_up: " << (pileupWeightSysUpHist ? "SUCCESS" : "FAILED") << std::endl;
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
    razorLog << "Loaded npu_down: " << (pileupWeightSysDownHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "Pileup loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

void RazorHelper::loadJetVeto_Summer24() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading jet veto map for Summer24..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    
    razorLog << "Opening file: Winter24Prompt24_2024BCDEFGHI.root" << std::endl;
    razorLog.flush();
    JetVetoFile = TFile::Open("Winter24Prompt24_2024BCDEFGHI.root");
    
    if (!JetVetoFile) {
        cout<<"Jet Veto File Not Found"<<endl;
        razorLog << "ERROR: Jet Veto File Not Found!" << std::endl;
        razorLog.flush();
        razorLog.close();
        return;
    }
    razorLog << "Jet veto file opened successfully" << std::endl;
    razorLog.flush();
    
    JetVetoHist = (TH2D*)JetVetoFile->Get("jetvetomap");
    razorLog << "Loaded jetvetomap: " << (JetVetoHist ? "SUCCESS" : "FAILED") << std::endl;
    JetVetoFpixHist = (TH2D*)JetVetoFile->Get("jetvetomap_fpix");
    razorLog << "Loaded jetvetomap_fpix: " << (JetVetoFpixHist ? "SUCCESS" : "FAILED") << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "Jet veto loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}

////////////////////////////////////////////////
//  Utilities
////////////////////////////////////////////////


double RazorHelper::getMetTriggerEff(float met) {
    if (MetTriggerHist) {
        return MetTriggerHist->GetBinContent(MetTriggerHist->GetXaxis()->FindFixBin(met));
    }
    else {
        std::cout << "RazorHelper error: MET trigger eff requested, but no histogram available!" << std::endl;
        return 0;
    }
}


double RazorHelper::getMetTriggerEffUp(float met) {
    if (MetTriggerSysUpHist) {
        return MetTriggerSysUpHist->GetBinContent(MetTriggerSysUpHist->GetXaxis()->FindFixBin(met));
    }
    else {
        std::cout << "RazorHelper error: MET trigger eff requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getMetTriggerEffDown(float met) {
    if (MetTriggerSysDownHist) {
        return MetTriggerSysDownHist->GetBinContent(MetTriggerSysDownHist->GetXaxis()->FindFixBin(met));
    }
    else {
        std::cout << "RazorHelper error: MET trigger eff requested, but no histogram available!" << std::endl;
        return 0;
    }
}



//// GET JET UNC for 22-24
double RazorHelper::getJecUnc( float pt, float eta , int run) {

    int foundIndex = -1;
    for (uint i=0; i<JetCorrectionsIOV.size(); i++) {
      if (run >= JetCorrectionsIOV[i].first && run <= JetCorrectionsIOV[i].second) {
        foundIndex = i;
      }
    }
    if (foundIndex == -1) {
      std::cout << "Warning: run = " << run << " was not found in any valid IOV range. use default index = 0 for Jet energy corrections. \n";
      foundIndex = 0;
    }
  
    jecUnc[foundIndex]->setJetPt(pt);
    jecUnc[foundIndex]->setJetEta(eta);
    return jecUnc[foundIndex]->getUncertainty(true);
  }

void RazorHelper::loadJECs() {
    std::ofstream razorLog("debug_razorhelper.log", std::ios::app);
    razorLog << "Loading jet energy corrections..." << std::endl;
    razorLog.flush();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "RazorHelper: loading jet energy correction constants, using 2018_17SeptEarlyReReco" << std::endl;
    // initialize
    std::string jecPathname = "JEC/";
    razorLog << "JEC pathname: " << jecPathname << std::endl;
    razorLog.flush();

    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    
    if (isData) {
      razorLog << "Loading JEC for data..." << std::endl;
      razorLog.flush();
      std::string jecUncPathA = jecPathname+"/Summer22_22Sep2023_V2_MC_Uncertainty_AK4PFPuppi"; //place holder for now, since we dont evaluate on data
      razorLog << "JEC uncertainty path: " << jecUncPathA << std::endl;
      razorLog.flush();
      JetCorrectionUncertainty *jecUncA = new JetCorrectionUncertainty(jecUncPathA);
      jecUnc.push_back(jecUncA);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 352319, 387121 ));
      razorLog << "JEC for data loaded successfully" << std::endl;
    }
    else {
      razorLog << "Loading JEC for MC..." << std::endl;
      razorLog.flush();
      std::cout << "Loading Jet Energy Corrections: 22-24 \n";
      std::string jecUncPath = jecPathname+"/Summer22_22Sep2023_V2_MC_Uncertainty_AK4PFPuppi.txt"; //same file works for 22-24 2025/04/25
      razorLog << "JEC uncertainty path: " << jecUncPath << std::endl;
      razorLog.flush();
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
      razorLog << "JEC for MC loaded successfully" << std::endl;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    razorLog << "JEC loading took " << duration.count() << " ms" << std::endl;
    razorLog.flush();
    razorLog.close();
}


//implemented based on these: https://gitlab.cern.ch/cms-jetmet/coordination/coordination/-/issues/117
bool RazorHelper::jetTightLepVeto(std::string tag, bool tightVeto, float Jet_neHEF, float Jet_neEmEF, float Jet_chEmEF, float Jet_muEF, float Jet_chHEF, UChar_t Jet_chMultiplicity,UChar_t Jet_neMultiplicity, float Jet_eta, bool Jet_jetId){
	bool Jet_passJetIdTightLepVeto = false;
    bool Jet_passJetIdTight = false;

    if (tag == "Summer24"){
        if (abs(Jet_eta) <= 2.6) Jet_passJetIdTight = (Jet_neHEF < 0.99) && (Jet_neEmEF < 0.9) && (Jet_chMultiplicity+Jet_neMultiplicity > 1) && (Jet_chHEF > 0.01) && (Jet_chMultiplicity > 0);
        else if (abs(Jet_eta) > 2.6 && abs(Jet_eta) <= 2.7)Jet_passJetIdTight = (Jet_neHEF < 0.90) && (Jet_neEmEF < 0.99);
        else if (abs(Jet_eta) > 2.7 && abs(Jet_eta) <= 3.0)Jet_passJetIdTight = (Jet_neHEF < 0.99);
        else if (abs(Jet_eta) > 3.0)Jet_passJetIdTight = (Jet_neMultiplicity >= 2) && (Jet_neEmEF < 0.4);

        if (abs(Jet_eta) <= 2.7) Jet_passJetIdTightLepVeto = Jet_passJetIdTight && (Jet_muEF < 0.8) && (Jet_chEmEF < 0.8);
		else Jet_passJetIdTightLepVeto = Jet_passJetIdTight;
	}
	else{
        if (abs(Jet_eta) <= 2.7) Jet_passJetIdTight = Jet_jetId & (1 << 1);
        else if (abs(Jet_eta) > 2.7 && abs(Jet_eta) <= 3.0) Jet_passJetIdTight = (Jet_jetId & (1 << 1)) && (Jet_neHEF < 0.99);
        else if (abs(Jet_eta) > 3.0) Jet_passJetIdTight = (Jet_jetId & (1 << 1)) && (Jet_neEmEF < 0.4);

        if (abs(Jet_eta) <= 2.7) Jet_passJetIdTightLepVeto = Jet_passJetIdTight && (Jet_muEF < 0.8) && (Jet_chEmEF < 0.8);
        else Jet_passJetIdTightLepVeto = Jet_passJetIdTight;

	}
    if (tightVeto) return Jet_passJetIdTightLepVeto;
    else return Jet_passJetIdTight;

}


double RazorHelper::getHMTTriggerEff(int chamber, int nhits){

    map<int, int> hist_cutoff;
    hist_cutoff[11] = 0;
    hist_cutoff[12] = 0;
    hist_cutoff[13] = 600;
    hist_cutoff[21] = 900;
    hist_cutoff[22] = 800;
    hist_cutoff[31] = 900;
    hist_cutoff[32] = 500;
    hist_cutoff[41] = 900;
    hist_cutoff[42] = 500;

    if (HMTEffHist[chamber]){
        return HMTEffHist[chamber]->GetEfficiency(HMTEffHist[chamber]->GetTotalHistogram()->GetXaxis()->FindFixBin(min(nhits,hist_cutoff[chamber])));
    }
    std::cout << "RazorHelper error: HMT efficiency requested, but no histogram available!" << std::endl;
    return 0;
}



double RazorHelper::getPileupWeight(int NPU) {
    if (pileupWeightHist) {
        return pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getPileupWeightUp(int NPU) {
    if (pileupWeightSysUpHist) {
        return pileupWeightSysUpHist->GetBinContent(pileupWeightSysUpHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: 'up' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getPileupWeightDown(int NPU) {
    if (pileupWeightSysDownHist) {
        return pileupWeightSysDownHist->GetBinContent(pileupWeightSysDownHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: 'down' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}


double RazorHelper::getJetVetoMap(float eta, float phi) {
    if (JetVetoHist) {
        return JetVetoHist->GetBinContent(JetVetoHist->GetXaxis()->FindFixBin(eta), JetVetoHist->GetYaxis()->FindFixBin(phi));
    }
    else {
        std::cout << "RazorHelper error: jet veto map requested, but no histogram available!" << std::endl;
        return 0;
    }
}


double RazorHelper::getJetVetoFpixMap(float eta, float phi) {
    if (JetVetoFpixHist) {
        return JetVetoFpixHist->GetBinContent(JetVetoFpixHist->GetXaxis()->FindFixBin(eta), JetVetoFpixHist->GetYaxis()->FindFixBin(phi));
    }
    else {
        std::cout << "RazorHelper error: jet veto map requested, but no histogram available!" << std::endl;
        return 0;
    }
}
