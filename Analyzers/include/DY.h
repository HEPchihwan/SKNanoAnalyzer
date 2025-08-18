#ifndef DY_h
#define DY_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class DY : public AnalyzerCore {
public:
    DY();
    ~DY();

    void initializeAnalyzer();
    void executeEvent();
    void executeEventFromParameter();

    // Analysis flags
    bool RunSyst;
    
    // Trigger settings
    TString IsoMuTriggerName;
    float TriggerSafePtCut;
    
    // Object ID settings
    RVec<Muon::MuonID> MuonIDs;
    RVec<TString> MuonIDSFKeys;
    
    
    // Physics objects
    RVec<Muon> AllMuons;
    
    RVec<Muon> selectedMuons;
    float dilepton_mass;
    
    // Analysis cuts
    struct AnalysisCuts {
        float muon_pt_lead = 26.0;
        float muon_pt_sublead = 26.0;
        float muon_eta = 2.4;
        float deltaR_overlap = 0.4;
    } cuts;
    
    // Systematic helper
    unique_ptr<SystematicHelper> systHelper;
    
    // Beamspot constrained variables
    float Muon_bsConstrainedChi2;   // Beamspot constraint를 적용한 χ²
    float Muon_bsConstrainedPt;     // BS constraint 후 재계산된 pT
    float Muon_bsConstrainedPtErr;  // pT 오차
    
    // Helper functions
    RVec<Muon> SelectMuons(const RVec<Muon>& muons);
    RVec<Muon> SelectMuonssublead(const RVec<Muon>& muons);
    RVec<Muon> RemoveOverlap(const RVec<Muon>& muons);
    pair<Muon, Muon> selectBestZPair(const RVec<Muon>& muons);
    
    

};

#endif