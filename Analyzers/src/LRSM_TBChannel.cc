#include "LRSM_TBChannel.h"

LRSM_TBChannel::LRSM_TBChannel() {}
LRSM_TBChannel::~LRSM_TBChannel() {}

void LRSM_TBChannel::initializeAnalyzer() {
    cout << "[LRSM_TBChannel::initializeAnalyzer] Starting initialization" << endl;
    
    // Check user flags
    RunSyst = HasFlag("RunSyst");
    RunWRCut = HasFlag("RunWRCut");
    
    cout << "[LRSM_TBChannel::initializeAnalyzer] RunSyst = " << RunSyst << endl;
    cout << "[LRSM_TBChannel::initializeAnalyzer] RunWRCut = " << RunWRCut << endl;
    
    // Set WR mass cut threshold
    if (RunWRCut) {
        WRCutThreshold = SelectionCuts::WR_CUT_2000;
    } else {
        WRCutThreshold = SelectionCuts::NO_WR_CUT;
    }
    
    // Muon IDs and scale factor keys
    MuonIDs.clear();
    // Use more reasonable muon IDs for standard analysis
    // POG_GLOBAL_HIGH_PT is very restrictive (for >200 GeV muons)
    // For LRSM analysis, use Tight ID + isolation
    // Your data has HighPtId=1, so use POG_TRACKER_HIGH_PT instead of POG_GLOBAL_HIGH_PT
    MuonIDs.push_back(Muon::MuonID::POG_GLOBAL_HIGH_PT);  // This matches your data (HighPtId=2)
    MuonIDs.push_back(Muon::MuonID::POG_TKISO_TIGHT);      // TkIsoId=2
    
    // Alternative: Use standard IDs if high-pT selection isn't critical
    // MuonIDs.push_back(Muon::MuonID::POG_TIGHT);
    // MuonIDs.push_back(Muon::MuonID::POG_PFISO_TIGHT);
    

    // Jet IDs
    JetIDs = {Jet::JetID::NOCUT};
    
    // Era-dependent trigger settings
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP" || DataEra == "2018") {
        IsoMuTriggerName = "HLT_IsoMu27";
        TriggerSafePtCut = 29.;
    } else if (DataEra == "2017") {
        IsoMuTriggerName = "HLT_IsoMu27"; 
        TriggerSafePtCut = 29.;
    } else if (DataEra == "2022") {
        Trigger1  = "HLT_Mu50";
        Trigger2  = "HLT_CascadeMu100";
        Trigger3  = "HLT_HighPtTkMu100";
        TriggerSafePtCut = 52.;
    } else if (DataEra == "2022EE") {
        Trigger1  = "HLT_Mu50";
        Trigger2  = "HLT_CascadeMu100";
        Trigger3  = "HLT_HighPtTkMu100";
        TriggerSafePtCut = 52.;
    } else if (DataEra == "2023") {
        Trigger1  = "HLT_Mu50";
        Trigger2  = "HLT_CascadeMu100";
        Trigger3  = "HLT_HighPtTkMu100";
        TriggerSafePtCut = 52.;
    } else if (DataEra == "2023BPix") {
        Trigger1  = "HLT_Mu50";
        Trigger2  = "HLT_CascadeMu100";
        Trigger3  = "HLT_HighPtTkMu100";
        TriggerSafePtCut = 52.;
    } else {
        cerr << "[LRSM_TBChannel::initializeAnalyzer] DataEra is not set properly: " << DataEra << endl;
        exit(EXIT_FAILURE);
    }
    
    //cout << "[LRSM_TBChannel::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
    cout << "[LRSM_TBChannel::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;
    
    // Initialize corrections
    if (IsDATA){
        if (DataEra == "2022") {
            corr_C = new MyCorrection("2022", "C", MCSample, true);
            corr_D = new MyCorrection("2022", "D", MCSample, true);
            corr_sm = new MyCorrection("2022", "SingleMuon", MCSample, true);

            
        } else if (DataEra == "2022EE") {
            corr_E = new MyCorrection("2022EE", "E", MCSample, true);
            corr_F = new MyCorrection("2022EE", "F", MCSample, true);
            corr_G = new MyCorrection("2022EE", "G", MCSample, true);
        }
        
    } else{
        if (DataEra == "2022") {
            corr_C = new MyCorrection("2022", "C", MCSample, false);
            corr_D = new MyCorrection("2022", "D", MCSample, false);
            corr_sm = new MyCorrection("2022", "SingleMuon", MCSample, false);
        } else if (DataEra == "2022EE") {
            corr_E = new MyCorrection("2022EE", "E", MCSample, false);
            corr_F = new MyCorrection("2022EE", "F", MCSample, false);
            corr_G = new MyCorrection("2022EE", "G", MCSample, false);
        }
    }

    



    
    // Initialize systematic helper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
    }
    
    cout << "[LRSM_TBChannel::initializeAnalyzer] Initialization complete" << endl;
}

void LRSM_TBChannel::executeEvent() {
    // Get all physics objects at the beginning to save CPU time
    AllMuons = GetAllMuons();
    AllJets = GetAllJets();
    AllFatJets = GetAllFatJets();
    
    
    
    // Loop over systematic sources
    for (const auto &syst_dummy : *systHelper) {
        executeEventFromParameter();
    }
}

void LRSM_TBChannel::executeEventFromParameter() {
    const TString this_syst = systHelper->getCurrentSysName();
    
    // Get event information
    Event ev = GetEvent();
    FillHist(this_syst + "/sumSign" + this_syst, sumSign, 1 , 10 , 0 , 1e+11 );
    FillHist(this_syst + "/CutFlow", 0.0, 1.0, 10, 0., 10.); // Initial event
    // Apply HLT trigger
    if (!(ev.PassTrigger(Trigger1)||ev.PassTrigger(Trigger2)||ev.PassTrigger(Trigger3))) return;
    
    FillHist(this_syst + "/CutFlow", 1.0, 1.0, 10, 0., 10.); // HLT pass
    
    // Copy physics objects for systematic variations
    RVec<Muon> muons = AllMuons;
    RVec<Jet> jets = AllJets;
    RVec<FatJet> fatjets = AllFatJets;
    
    
    
    
    

    // Muon Id pass 
    bool hasGoodMuon = false;
    
    for (const auto& muon : muons) {
        // Debug output to see actual ID values
        FillHist(this_syst + "/MuonhighPtid",muon.PassID(MuonIDs[0]), 1.0, 10, -5., 5.);
        FillHist(this_syst + "/Muonisoid",muon.PassID(MuonIDs[1]), 1.0, 10, -5., 5.);
        
        // Fill additional histograms to understand what IDs are available
        FillHist(this_syst + "/Muon_HighPtId", (int)muon.HighPtId(), 1.0,  10, -5., 5.);
        FillHist(this_syst + "/Muon_TkIsoId", (int)muon.TkIsoId(), 1.0,  10, -5., 5.);
        FillHist(this_syst + "/Muon_TightId", muon.isPOGTightId(), 1.0, 3, 0., 3.);
        FillHist(this_syst + "/Muon_MediumId", muon.isPOGMediumId(), 1.0, 3, 0., 3.);
        FillHist(this_syst + "/Muon_LooseId", muon.isPOGLooseId(), 1.0, 3, 0., 3.);
        
        if (muon.PassID(MuonIDs[1]) and muon.PassID(MuonIDs[0]) ) {
            hasGoodMuon = true;
            break;
        }
    }   
    if (!hasGoodMuon) return;

    // Apply muon selection
    FillHist(this_syst + "/CutFlow", 2.0, 1.0, 10, 0., 10.); // 2 muons
    muons = RemoveOverlap(muons);
    // Require more than 2 muons
    if (muons.size() < 2) return;
    
    // Sort muons by pT
    sort(muons.begin(), muons.end(), PtComparing);
    
    // Apply kinematic cuts
    if (!PassKinematicCuts(muons)) return;
    
    FillHist(this_syst + "/CutFlow", 3.0, 1.0, 10, 0., 10.); // Kinematic cuts
    
    // Apply dilepton mass cut
    if (!PassDileptonMassCut(muons)) return;
    
    Muon muon1 = muons[0];
    Muon muon2 = muons[1];
    muon_overlap_cleaned = { muon1, muon2 };
    


    FillHist(this_syst + "/CutFlow", 4.0, 1.0, 10, 0., 10.); // Dilepton mass cut
    
    
    
    // Select fat jets and remove overlaps
    // FatJet selection - using basic kinematic cuts for now
    RVec<FatJet> selected_fatjets;
    for (const auto& fj : fatjets) {
        if (fj.Pt() > cuts.fatjet_pt && abs(fj.Eta()) < cuts.fatjet_eta) {
            selected_fatjets.push_back(fj);
        }
    }
    fatjets = selected_fatjets;
    FillHist(this_syst + "/FatJetnum", fatjets.size(), 1.0, 10, 0., 10.);
    fatjets = RemoveOverlapWithMuonsFatJet(fatjets, muon_overlap_cleaned);
    FillHist(this_syst + "/FatJetnum_afterOverlap", fatjets.size(), 1.0, 10, 0., 10.);
    RVec<FatJet> topjets = SelectTopTaggedJets(fatjets);
    
    for (const auto& fatjet : fatjets) {
        // Using basic mass cuts for top tagging - update with actual tagger when available
        float toptag_score1 = fatjet.GetTaggerResult(JetTagging::FatJetTaggingtype::ParticleNetWithMass, JetTagging::FatjetTaggingObject::TvsQCD); // placeholder
        float softdrop_mass1 = fatjet.SDMass();
        FillHist(this_syst + "/FatJet_SoftDropMass", softdrop_mass1, 1.0, 100, 0., 1000.);
        FillHist(this_syst + "/FatJet_TopTagScore", toptag_score1, 1.0, 100, 0., 1.);
        
    }
    for (const auto& topjet : topjets) {
        float toptag_score2 = topjet.GetTaggerResult(JetTagging::FatJetTaggingtype::ParticleNetWithMass, JetTagging::FatjetTaggingObject::TvsQCD);
        float softdrop_mass2 = topjet.SDMass();
        FillHist(this_syst + "/topJet_SoftDropmass", softdrop_mass2, 1.0, 100, 0., 1000.);
        FillHist(this_syst + "/topJet_TopTagScore", toptag_score2, 1.0, 100, 0., 1.);
    }
    if (topjets.size() < 1) return;
    sort(topjets.begin(), topjets.end(), PtComparing);
    RVec<FatJet> leading_topjet = {topjets[0]};
    FillHist(this_syst + "/CutFlow", 5.0, 1.0, 10, 0., 10.);
    // Remove overlap between jets and fat jets
    jets = SelectJets(jets, JetIDs[0], cuts.jet_pt, cuts.jet_eta);
    jets = RemoveOverlapWithMuons(jets, muon_overlap_cleaned);
    jets = RemoveOverlapWithFatJets(jets, leading_topjet);
    RVec<Jet> bjets = SelectBTaggedJets(jets);
    sort(bjets.begin(), bjets.end(), PtComparing);
    if (bjets.size() < 1 ) return;
    
    RVec<Jet> leading_bjet = {bjets[0]};
    FillHist(this_syst + "/CutFlow", 6.0, 1.0, 10, 0., 10.); // b-jet and top-jet
    
    
    
    float weight = 1.0;
    // Calculate invariant masses
    float wr_mass = CalculateWRMass(muon_overlap_cleaned, leading_bjet, leading_topjet);
    float dilepton_mass = (muon_overlap_cleaned[0] + muon_overlap_cleaned[1]).M();
    
    // Apply WR mass cut if requested
    if ( wr_mass <0.00) return;
    
    FillHist(this_syst + "/CutFlow", 7.0, 1.0, 10, 0., 10.); // WR mass cut (if applied)
    

    // correction 

    float corr1_1  = corr_C -> GetMuonScaleSF(muon1, variation::this_syst, muon1.Pt());
    float corr1_2 = corr_C -> GetMuonScaleSF(muon2, variation::this_syst, muon2.Pt());

    float coor_2_1  = coor_C -> GetMuonIDSF( , muon1, variation::this_syst);
    float coor_2_2  = coor_C -> GetMuonIDSF( , muon2, variation::this_syst);




    float corr2  = corr_C -> GetMuonIDSF(MuonIDSFKeys[0], muon, variation::nom);
    float corr3  = corr_C -> GetMuonRECOSF(muon, variation::nom);

    // Event weight calculation
    
    if (!IsDATA) {
        weight *= MCweight();
        //cout << "[LRSM_TBChannel::executeEventFromParameter] MC weight: " << MCweight() << endl;
        weight *= ev.GetTriggerLumi("Full");
        //cout << "[LRSM_TBChannel::executeEventFromParameter] Trigger lumi : " << ev.GetTriggerLumi("Full") << endl;
        //cout << "[LRSM_TBChannel::executeEventFromParameter] Event weight: " << weight << endl;
        
        FillHist(this_syst + "/xsec" + this_syst, xsec , 1, 100 , 0 , 1000 );
        FillHist(this_syst + "/Bjetnum", bjets.size(), weight, 10, 0., 10.);
        FillHist(this_syst + "/Topjetnum", topjets.size(), weight, 10, 0., 10.);
        FillHist(this_syst + "/WRMass_" + this_syst, wr_mass, weight, 2000, 0., 2000.);
        FillHist(this_syst + "/DileptonMass_" + this_syst, dilepton_mass, weight, 5000, 0., 5000.);
        FillHist(this_syst + "/LeadingMuonPt_" + this_syst, muon_overlap_cleaned[0].Pt(), weight, 5000, 0., 5000.);
        FillHist(this_syst + "/SubleadingMuonPt_" + this_syst, muon_overlap_cleaned[1].Pt(), weight, 5000, 0., 5000.);
        FillHist(this_syst + "/LeadingBJetPt_" + this_syst, leading_bjet[0].Pt(), weight, 5000, 0., 5000.);
        FillHist(this_syst + "/LeadingTopJetPt_" + this_syst, leading_topjet[0].Pt(), weight, 5000, 0., 5000.);
        
        // Apply systematic weights
        //unordered_map<std::string, float> weight_map = systHelper->calculateWeight();
        //for (const auto &w : weight_map) {
        //    TString weight_suffix = w.first;
        //    float total_weight = weight * w.second;
            
            // Fill histograms with systematic weights
        //    FillHist(this_syst + "/WRMass_" + weight_suffix, wr_mass, total_weight, 100, 0., 8000.);
            //FillHist(this_syst + "/DileptonMass_" + weight_suffix, dilepton_mass, total_weight, 100, 0., 8000.);
            //FillHist(this_syst + "/LeadingMuonPt_" + weight_suffix, muons[0].Pt(), total_weight, 100, 0., 8000.);
            //FillHist(this_syst + "/SubleadingMuonPt_" + weight_suffix, muons[1].Pt(), total_weight, 100, 0., 8000.);
            //FillHist(this_syst + "/LeadingBJetPt_" + weight_suffix, leading_bjet[0].Pt(), total_weight, 100, 0., 8000.);
            //FillHist(this_syst + "/LeadingTopJetPt_" + weight_suffix, leading_topjet[0].Pt(), total_weight, 100, 0., 8000.);
        //}
    } else {
        // For data, only fill nominal histograms
        FillHist(this_syst + "/Bjetnum", bjets.size(), weight, 10, 0., 10.);
        FillHist(this_syst + "/Topjetnum", topjets.size(), weight, 10, 0., 10.);
        FillHist(this_syst + "/WRMass_" + this_syst, wr_mass, weight, 2000, 0., 2000.);
        FillHist(this_syst + "/DileptonMass_" + this_syst, dilepton_mass, weight, 5000, 0., 5000.);
        FillHist(this_syst + "/LeadingMuonPt_" + this_syst, muon_overlap_cleaned[0].Pt(), weight, 5000, 0., 5000.);
        FillHist(this_syst + "/SubleadingMuonPt_" + this_syst, muon_overlap_cleaned[1].Pt(), weight, 5000, 0., 5000.);
        FillHist(this_syst + "/LeadingBJetPt_" + this_syst, leading_bjet[0].Pt(), weight, 5000, 0., 5000.);
        FillHist(this_syst + "/LeadingTopJetPt_" + this_syst, leading_topjet[0].Pt(), weight, 5000, 0., 5000.);
    }
}

// Helper function implementations



RVec<Jet> LRSM_TBChannel::SelectBTaggedJets(const RVec<Jet>& jets) {
    RVec<Jet> btagged_jets;
    for (const auto& jet : jets) {
        if (jet.GetBTaggerResult(JetTagging::JetFlavTagger::ParticleNet) > cuts.btag_wp) {
            btagged_jets.push_back(jet); 
        }
    }
    return btagged_jets;
}

RVec<FatJet> LRSM_TBChannel::SelectTopTaggedJets(const RVec<FatJet>& fatjets) {
    RVec<FatJet> toptagged_jets;
    for (const auto& fatjet : fatjets) {
        // Using basic mass cuts for top tagging - update with actual tagger when available
        float toptag_score = fatjet.GetTaggerResult(JetTagging::FatJetTaggingtype::ParticleNetWithMass, JetTagging::FatjetTaggingObject::TvsQCD); // placeholder
        float softdrop_mass = fatjet.SDMass();
        
        if (toptag_score > cuts.toptag_score &&
            softdrop_mass > cuts.toptag_mass_low &&
            softdrop_mass < cuts.toptag_mass_high) {
            toptagged_jets.push_back(fatjet);
        }
    }
    return toptagged_jets;
}

RVec<Muon> LRSM_TBChannel::RemoveOverlap(const RVec<Muon>& muons, float deltaR_cut) {
    RVec<Muon> cleaned_muons;
    for (size_t i = 0; i < muons.size(); ++i) {
        bool overlaps = false;
        for (size_t j = i + 1; j < muons.size(); ++j) {
            if (muons[i].DeltaR(muons[j]) < deltaR_cut) {
                overlaps = true;
                break;
            }
        }
        if (!overlaps) {
            cleaned_muons.push_back(muons[i]);
        }
    }
    return cleaned_muons;
}

RVec<Jet> LRSM_TBChannel::RemoveOverlapWithMuons(const RVec<Jet>& jets, const RVec<Muon>& muons, float deltaR_cut) {
    RVec<Jet> cleaned_jets;
    for (const auto& jet : jets) {
        bool overlaps = false;
        for (const auto& muon : muons) {
            if (jet.DeltaR(muon) < deltaR_cut) {
                overlaps = true;
                break;
            }
        }
        if (!overlaps) {
            cleaned_jets.push_back(jet);
        }
    }
    return cleaned_jets;
}

RVec<Jet> LRSM_TBChannel::RemoveOverlapWithFatJets(const RVec<Jet>& jets, const RVec<FatJet>& fatjets, float deltaR_cut) {
    RVec<Jet> cleaned_jets;
    for (const auto& jet : jets) {
        bool overlaps = false;
        for (const auto& fatjet : fatjets) {
            if (jet.DeltaR(fatjet) < deltaR_cut) {
                overlaps = true;
                break;
            }
        }
        if (!overlaps) {
            cleaned_jets.push_back(jet);
        }
    }
    return cleaned_jets;
}

RVec<FatJet> LRSM_TBChannel::RemoveOverlapWithMuonsFatJet(const RVec<FatJet>& fatjets, const RVec<Muon>& muons, float deltaR_cut) {
    RVec<FatJet> cleaned_fatjets;
    for (const auto& fatjet : fatjets) {
        bool overlaps = false;
        for (const auto& muon : muons) {
            if (fatjet.DeltaR(muon) < deltaR_cut) {
                overlaps = true;
                break;
            }
        }
        if (!overlaps) {
            cleaned_fatjets.push_back(fatjet);
        }
    }
    return cleaned_fatjets;
}

bool LRSM_TBChannel::PassKinematicCuts(const RVec<Muon>& muons) {
    if (muons.size() < 2) return false;
    
    // Leading muon pT cut
    if (muons[0].Pt() <= 50 ) return false;
    
    // Eta cuts
    for (const auto& muon : muons) {
        if (fabs(muon.Eta()) >= cuts.muon_eta) return false;
    }
    
    return true;
}

bool LRSM_TBChannel::PassDileptonMassCut(const RVec<Muon>& muons) {
    if (muons.size() < 2) return false;
    
    float dilepton_mass = (muons[0] + muons[1]).M();
    return dilepton_mass > cuts.dilepton_mass_cut;
}

float LRSM_TBChannel::CalculateWRMass(const RVec<Muon>& muons, const RVec<Jet>& bjets, const RVec<FatJet>& topjets) {
    if (muons.size() < 2 || bjets.size() < 1 || topjets.size() < 1) return -1.0;
    
    Particle wr_candidate = muons[0] + muons[1] + bjets[0] + topjets[0];
    return wr_candidate.M();
}

