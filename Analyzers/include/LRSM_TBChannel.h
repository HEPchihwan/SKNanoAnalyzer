#ifndef LRSM_TBChannel_h
#define LRSM_TBChannel_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class LRSM_TBChannel : public AnalyzerCore {
public:
    LRSM_TBChannel();
    ~LRSM_TBChannel();

    void initializeAnalyzer();
    void executeEvent();
    void executeEventFromParameter();

    // Analysis flags
    bool RunSyst;
    bool RunWRCut;
    
    // Selection parameters
    enum class SelectionCuts {
        NO_WR_CUT = 0,
        WR_CUT_2000 = 2000
    };
    
    SelectionCuts WRCutThreshold;
    
    // Trigger settings
    TString IsoMuTriggerName;
    TString Trigger1;
    TString Trigger2;
    TString Trigger3;
    float TriggerSafePtCut;
    
    // Object ID settings
    RVec<Muon::MuonID> MuonIDs;
    RVec<TString> MuonIDSFKeys;
    RVec<Jet::JetID> JetIDs;
    
    // Physics objects
    RVec<Muon> AllMuons;
    RVec<Jet> AllJets;
    RVec<FatJet> AllFatJets;
    RVec<Muon> muon1;
    RVec<Muon> muon2;
    RVec<Muon> muon_overlap_cleaned;
    // Analysis cuts
    struct AnalysisCuts {
        float muon_pt = 50.0;
        float muon_eta = 2.5;
        float jet_pt = 30.0;
        float jet_eta = 2.5;
        float fatjet_pt = 30.0;
        float fatjet_eta = 2.5;
        float btag_wp = 0.6734; // ParticleNet medium WP
        float toptag_score = 0.9;
        float toptag_mass_low = 120.0;
        float toptag_mass_high = 250.0;
        float deltaR_overlap = 0.4;
        float deltaR_fatjet_overlap = 0.8;
        float dilepton_mass_cut = 50;
    } cuts;
    
    // Weight variables
    float weight_Prefire;
    
    // Systematic helper
    unique_ptr<SystematicHelper> systHelper;
    
    // Helper functions
    RVec<Muon> SelectHighPtMuons(const RVec<Muon>& muons);
    RVec<Jet> SelectBTaggedJets(const RVec<Jet>& jets);
    RVec<FatJet> SelectTopTaggedJets(const RVec<FatJet>& fatjets);
    RVec<Muon> RemoveOverlap(const RVec<Muon>& muons, float deltaR_cut = 0.4);
    RVec<Jet> RemoveOverlapWithMuons(const RVec<Jet>& jets, const RVec<Muon>& muons, float deltaR_cut = 0.4);
    RVec<Jet> RemoveOverlapWithFatJets(const RVec<Jet>& jets, const RVec<FatJet>& fatjets, float deltaR_cut = 0.8);
    RVec<FatJet> RemoveOverlapWithMuonsFatJet(const RVec<FatJet>& fatjets, const RVec<Muon>& muons, float deltaR_cut = 0.8);
    bool PassEventSelection(const RVec<Muon>& muons, const RVec<Jet>& bjets, const RVec<FatJet>& topjets);
    bool PassKinematicCuts(const RVec<Muon>& muons);
    bool PassDileptonMassCut(const RVec<Muon>& muons);
    float CalculateWRMass(const RVec<Muon>& muons, const RVec<Jet>& bjets, const RVec<FatJet>& topjets);
    float CalculateNeutrinoMass(const RVec<Muon>& muons, const RVec<Jet>& bjets, const RVec<FatJet>& topjets);
};

#endif