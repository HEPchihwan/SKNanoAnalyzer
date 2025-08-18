#include "DY.h"
#include <utility>

DY::DY() {}
DY::~DY() {}

void DY::initializeAnalyzer() {
    cout << "[DY::initializeAnalyzer] Starting initialization" << endl;
    
    // Check user flags
    RunSyst = HasFlag("RunSyst");
    
    cout << "[DY::initializeAnalyzer] RunSyst = " << RunSyst << endl;
    
    // Muon IDs and scale factor keys
    MuonIDs.clear();
    MuonIDs.push_back(Muon::MuonID::POG_TIGHT);
    //MuonIDSFKeys = {"NUM_TightID_DEN_TrackerMuons"};
    
    // Jet IDs
    
    
    // Era-dependent trigger settings
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP" || DataEra == "2018") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2017") {
        IsoMuTriggerName = "HLT_IsoMu24"; 
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2022") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2022EE") {
        IsoMuTriggerName = "HLT_Mu15";
        TriggerSafePtCut = 15.;
    } else if (DataEra == "2023") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2023BPix") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.; 
    } else {
        cerr << "[DY::initializeAnalyzer] DataEra is not set properly: " << DataEra << endl;
        exit(EXIT_FAILURE);
    }
    
    cout << "[DY::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
    cout << "[DY::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;
    
    // Initialize corrections
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
    
    // Initialize systematic helper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
    }
    
    cout << "[DY::initializeAnalyzer] Initialization complete" << endl;
}

void DY::executeEvent() {
    // Get all physics objects at the beginning to save CPU time
    AllMuons = GetAllMuons();
    
    
    // Loop over systematic sources
    for (const auto &syst_dummy : *systHelper) {
        executeEventFromParameter();
    }
}

void DY::executeEventFromParameter() {
    const TString this_syst = systHelper->getCurrentSysName();
    
    // Get event information
    Event ev = GetEvent();
    FillHist(this_syst + "/CutFlow", 0.0, 1.0, 10, 0., 10.); // Initial event
    
    RVec<Muon> muons = AllMuons;
    

    // Apply HLT trigger (HLT_IsoMu24)
    if (!ev.PassTrigger(IsoMuTriggerName)) return;
    FillHist(this_syst + "/CutFlow", 1.0, 1.0, 10, 0., 10.); // HLT pass
    
    // Copy physics objects for systematic variations
    
    
    // Select muons
    
    selectedMuons = RemoveOverlap(muons);
    selectedMuons = SelectMuons(selectedMuons);
    if (selectedMuons.size() < 2) return;
    if (selectedMuons.size() > 2) return; // Only consider events with exactly 2 muons
    FillHist(this_syst + "/CutFlow", 2.0, 1.0, 10, 0., 10.); // At least 2 muons selected
    
    // Select best Z pair using beamspot constrained chi2
    pair<Muon, Muon> bestZPair = selectBestZPair(selectedMuons);
    Muon leading_muon = bestZPair.first;
    Muon subleading_muon = bestZPair.second;
    
    // Check if valid pair was found
    if (leading_muon.Pt() == 0 || subleading_muon.Pt() == 0) return;
    FillHist(this_syst + "/CutFlow", 3.0, 1.0, 10, 0., 10.); // Valid Z pair found
    
    // Ensure leading muon has higher pT
    if (subleading_muon.Pt() > leading_muon.Pt()) {
        swap(leading_muon, subleading_muon);
    }
    
    FillHist(this_syst + "/CutFlow", 4.0, 1.0, 10, 0., 10.); // Best Z pair selected
    
    float dilepton_mass = (leading_muon + subleading_muon).M();
    //if (dilepton_mass < 15.0 || dilepton_mass > 3000.0) return; // Dilepton mass cut
    
    
    
    // Event weight calculation
    float weight = 1.0;
    if (!IsDATA) {
        weight *= MCweight();
        weight *= ev.GetTriggerLumi("Full");
    }
    
    // Fill histograms
    FillHist(this_syst + "/DileptonMass", dilepton_mass, weight, 3000, 0., 3000.);
    FillHist(this_syst + "/LeadingMuonPt", leading_muon.Pt(), weight, 500, 0., 500.);
    FillHist(this_syst + "/SubleadingMuonPt", subleading_muon.Pt(), weight, 500, 0., 500.);
}

// Helper function implementations

RVec<Muon> DY::SelectMuons(const RVec<Muon>& muons) {
    RVec<Muon> selected_muons;
    for (const auto& muon : muons) {
        if (muon.Pt() > cuts.muon_pt_lead && 
            abs(muon.Eta()) < cuts.muon_eta &&
            muon.PassID(MuonIDs[0])) {
            selected_muons.push_back(muon);
        }
    }
    return selected_muons;
}



RVec<Muon> DY::RemoveOverlap(const RVec<Muon>& muons) {
    RVec<Muon> cleaned_muons;
    for (const auto& muon : muons) {
        bool overlaps = false;
        for (const auto& other_muon : muons) {
            if (muon.DeltaR(other_muon) < cuts.deltaR_overlap && muon != other_muon) {
                overlaps = true;
                break;
            }
        }
        if (!overlaps) {
            cleaned_muons.push_back(muon);
        }
    }
    return cleaned_muons;
}

// Z vertex quality 평가에 활용
pair<Muon, Muon> DY::selectBestZPair(const RVec<Muon>& muons) {
    double bestChi2Sum = 999999;
    pair<int,int> bestPairIndices = {-1, -1};
    
    for (int i = 0; i < muons.size(); i++) {
        for (int j = i+1; j < muons.size(); j++) {
            if (muons[i].Charge() * muons[j].Charge() < 0) {
                // 두 뮤온의 BS chi2 합
                double chi2Sum = muons[i].bsConstrainedChi2() + 
                                muons[j].bsConstrainedChi2();
                
                if (chi2Sum < bestChi2Sum) {
                    bestChi2Sum = chi2Sum;
                    bestPairIndices = {i, j};
                }
            }
        }
    }
    
    // Return the best muon pair
    if (bestPairIndices.first != -1 && bestPairIndices.second != -1) {
        return {muons[bestPairIndices.first], muons[bestPairIndices.second]};
    } else {
        // Return empty pair if no valid pair found
        return {Muon(), Muon()};
    }
}