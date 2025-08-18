#ifndef TTbar_test_h
#define TTbar_test_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class TTbar_test : public AnalyzerCore {
public:
    TTbar_test();
    ~TTbar_test();

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
    RVec<Jet::JetID> JetIDs;
    
    // Physics objects
    RVec<Muon> AllMuons;
    RVec<Jet> AllJets;
    RVec<Muon> selectedMuons;
    RVec<Jet> selectedJets;
    RVec<Jet> selectedBJets;
    RVec<Jet> selectedLightJets;
    
    // Analysis cuts
    struct AnalysisCuts {
        float muon_pt = 30.0;
        float muon_eta = 2.5;
        float jet_pt = 30.0;
        float jet_eta = 2.5;
        float btag_wp = 0.6734; // ParticleNet medium WP
        float deltaR_overlap = 0.4;
        float met_pt = 30.0;
    } cuts;
    
    // Systematic helper
    unique_ptr<SystematicHelper> systHelper;
    
    // Helper functions
    RVec<Muon> SelectMuons(const RVec<Muon>& muons);
    RVec<Jet> SelectJetsFromTTbar(const RVec<Jet>& jets);
    RVec<Jet> SelectBTaggedJets(const RVec<Jet>& jets);
    RVec<Jet> RemoveOverlapWithMuons(const RVec<Jet>& jets, const RVec<Muon>& muons, float deltaR_cut = 0.4);
    RVec<Jet> RemoveOverlapWithJets(const RVec<Jet>& jets, const RVec<Jet>& bjets, float deltaR_cut = 0.4);
    bool PassEventSelection(const RVec<Muon>& muons, const RVec<Jet>& bjets, const RVec<Jet>& lightjets);
    
    // TTbar system observables calculation
    struct TTbarObservables {
        float mt_ttbar_v1;
        float mt_ttbar_v2;
        float mt_W;
        float m_hadronic;
        float m_visible;
        float HT;
        float MET_over_HT;
        float pt_visible;
        float MET_vis_balance;
    };
    
    TTbarObservables CalculateTTbarObservables(const Muon& lepton, 
                                              const Jet& b1, const Jet& b2,
                                              const Jet& j1, const Jet& j2,
                                              float met_pt, float met_phi);
    
    float CalculateTTbarSystemTransverseMass(const Muon& lepton,
                                           const Jet& b1, const Jet& b2,
                                           const Jet& j1, const Jet& j2,
                                           float met_pt, float met_phi);
    
    float CalculateTTbarSystemTransverseMassV2(const Muon& lepton,
                                             const Jet& b1, const Jet& b2,
                                             const Jet& j1, const Jet& j2,
                                             float met_pt, float met_phi);
};

#endif