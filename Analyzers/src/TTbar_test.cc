#include "TTbar_test.h"

TTbar_test::TTbar_test() {}
TTbar_test::~TTbar_test() {}

void TTbar_test::initializeAnalyzer() {
    cout << "[TTbar_test::initializeAnalyzer] Starting initialization" << endl;
    
    // Check user flags
    RunSyst = HasFlag("RunSyst");
    
    cout << "[TTbar_test::initializeAnalyzer] RunSyst = " << RunSyst << endl;
    
    // Muon IDs and scale factor keys
    MuonIDs.clear();
    MuonIDs.push_back(Muon::MuonID::POG_TIGHT);
    //MuonIDSFKeys = {"NUM_TightID_DEN_TrackerMuons"};
    
    // Jet IDs
    JetIDs = {Jet::JetID::TIGHTLEPVETO};
    
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
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2023") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2023BPix") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.; 
    } else {
        cerr << "[TTbar_test::initializeAnalyzer] DataEra is not set properly: " << DataEra << endl;
        exit(EXIT_FAILURE);
    }
    
    cout << "[TTbar_test::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
    cout << "[TTbar_test::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;
    
    // Initialize corrections
    myCorr = new MyCorrection(DataEra, IsDATA ? DataStream : MCSample, IsDATA);
    
    // Initialize systematic helper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample);
    }
    
    cout << "[TTbar_test::initializeAnalyzer] Initialization complete" << endl;
}

void TTbar_test::executeEvent() {
    // Get all physics objects at the beginning to save CPU time
    AllMuons = GetAllMuons();
    AllJets = GetAllJets();
    
    // Loop over systematic sources
    for (const auto &syst_dummy : *systHelper) {
        executeEventFromParameter();
    }
}

void TTbar_test::executeEventFromParameter() {
    const TString this_syst = systHelper->getCurrentSysName();
    
    // Get event information
    Event ev = GetEvent();
    FillHist(this_syst + "/CutFlow", 0.0, 1.0, 10, 0., 10.); // Initial event
    
    RVec<Muon> muons = AllMuons;
    RVec<Jet> jets = AllJets;

    // Apply HLT trigger (HLT_IsoMu24)
    if (!ev.PassTrigger(IsoMuTriggerName)) return;
    FillHist(this_syst + "/CutFlow", 1.0, 1.0, 10, 0., 10.); // HLT pass
    
    // Copy physics objects for systematic variations
    
    
    // Select muons
    selectedMuons = SelectMuons(muons);
    if (selectedMuons.size() < 1) return;
    
    // Sort muons by pT and take the leading one
    sort(selectedMuons.begin(), selectedMuons.end(), PtComparing);
    Muon leading_muon = selectedMuons[0];
    FillHist(this_syst + "/CutFlow", 2.0, 1.0, 10, 0., 10.); // Leading muon selected
    
    // Select jets and remove overlap with muons
    selectedJets = SelectJetsFromTTbar(jets);
    selectedJets = RemoveOverlapWithMuons(selectedJets, selectedMuons);
    //selectedJets = RemoveOverlapWithJets(selectedJets,selectedJets );
    // Select b-tagged jets
    selectedBJets = SelectBTaggedJets(selectedJets);
    if (selectedBJets.size() < 2) return;
    sort(selectedBJets.begin(), selectedBJets.end(), PtComparing);
    FillHist(this_syst + "/CutFlow", 3.0, 1.0, 10, 0., 10.); // 2 b-jets selected
    
    // Select light jets (non-b-tagged jets)
    selectedLightJets.clear();
    //selectedLightJets = RemoveOverlapWithJets(selectedLightJets, selectedLightJets);
    selectedLightJets = RemoveOverlapWithJets(selectedJets, selectedBJets);
    
    if (selectedLightJets.size() < 2) return;
    sort(selectedLightJets.begin(), selectedLightJets.end(), PtComparing);
    FillHist(this_syst + "/CutFlow", 4.0, 1.0, 10, 0., 10.); // 2 light jets selected
    
    // Get MET information
    float met_pt = ev.GetMETVector(Event::MET_Type::PUPPI,MyCorrection::variation::nom,Event::MET_Syst::CENTRAL).Pt();
    float met_phi = ev.GetMETVector(Event::MET_Type::PUPPI,MyCorrection::variation::nom,Event::MET_Syst::CENTRAL).Phi();
    
    if (met_pt <= cuts.met_pt) return;
    FillHist(this_syst + "/CutFlow", 5.0, 1.0, 10, 0., 10.); // MET cut
    
    // Calculate TTbar observables
    TTbarObservables obs = CalculateTTbarObservables(leading_muon,
                                                    selectedBJets[0], selectedBJets[1],
                                                    selectedLightJets[0], selectedLightJets[1],
                                                    met_pt, met_phi);
    
    // Event weight calculation
    float weight = 1.0;
    if (!IsDATA) {
        weight *= MCweight();
        weight *= ev.GetTriggerLumi("Full");
    }
    
    // Fill histograms
    FillHist(this_syst + "/TTbarTransverseMass_v1", obs.mt_ttbar_v1, weight, 100, 0., 2000.);
    FillHist(this_syst + "/TTbarTransverseMass_v2", obs.mt_ttbar_v2, weight, 100, 0., 2000.);
    FillHist(this_syst + "/WTransverseMass", obs.mt_W, weight, 100, 0., 200.);
    FillHist(this_syst + "/HadronicMass", obs.m_hadronic, weight, 100, 0., 500.);
    FillHist(this_syst + "/VisibleMass", obs.m_visible, weight, 100, 0., 2000.);
    FillHist(this_syst + "/HT", obs.HT, weight, 100, 0., 1500.);
    FillHist(this_syst + "/MET_over_HT", obs.MET_over_HT, weight, 100, 0., 1.);
    FillHist(this_syst + "/VisiblePt", obs.pt_visible, weight, 100, 0., 500.);
    FillHist(this_syst + "/MET_VisBalance", obs.MET_vis_balance, weight, 100, -1., 1.);
    
    // Fill individual object histograms
    FillHist(this_syst + "/LeadingMuonPt", leading_muon.Pt(), weight, 100, 0., 500.);
    FillHist(this_syst + "/LeadingBJetPt", selectedBJets[0].Pt(), weight, 100, 0., 500.);
    FillHist(this_syst + "/SubleadingBJetPt", selectedBJets[1].Pt(), weight, 100, 0., 500.);
    FillHist(this_syst + "/LeadingLightJetPt", selectedLightJets[0].Pt(), weight, 100, 0., 500.);
    FillHist(this_syst + "/SubleadingLightJetPt", selectedLightJets[1].Pt(), weight, 100, 0., 500.);
    FillHist(this_syst + "/MET", met_pt, weight, 100, 0., 500.);
}

// Helper function implementations

RVec<Muon> TTbar_test::SelectMuons(const RVec<Muon>& muons) {
    RVec<Muon> selected_muons;
    for (const auto& muon : muons) {
        if (muon.Pt() > cuts.muon_pt && 
            abs(muon.Eta()) < cuts.muon_eta &&
            muon.PassID(MuonIDs[0])) {
            selected_muons.push_back(muon);
        }
    }
    return selected_muons;
}

RVec<Jet> TTbar_test::SelectJetsFromTTbar(const RVec<Jet>& jets) {
    RVec<Jet> selected_jets;
    for (const auto& jet : jets) {
        if (jet.Pt() > cuts.jet_pt && abs(jet.Eta()) < cuts.jet_eta 
        ) {
            selected_jets.push_back(jet);
        }
    }
    return selected_jets;
}

RVec<Jet> TTbar_test::SelectBTaggedJets(const RVec<Jet>& jets) {
    RVec<Jet> btagged_jets;
    for (const auto& jet : jets) {
        if (jet.GetBTaggerResult(JetTagging::JetFlavTagger::ParticleNet) > cuts.btag_wp) {
            btagged_jets.push_back(jet); 
        }
    }
    return btagged_jets;
}

RVec<Jet> TTbar_test::RemoveOverlapWithMuons(const RVec<Jet>& jets, const RVec<Muon>& muons, float deltaR_cut) {
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

RVec<Jet> TTbar_test::RemoveOverlapWithJets(const RVec<Jet>& jets, const RVec<Jet>& bjets, float deltaR_cut) {
    RVec<Jet> cleaned_jets;
    for (const auto& jet : jets) {
        bool overlaps = false;
        for (const auto& bjet : bjets) {
            if (jet.DeltaR(bjet) < deltaR_cut) {
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

bool TTbar_test::PassEventSelection(const RVec<Muon>& muons, const RVec<Jet>& bjets, const RVec<Jet>& lightjets) {
    return (muons.size() >= 1 && bjets.size() >= 2 && lightjets.size() >= 2);
}

TTbar_test::TTbarObservables TTbar_test::CalculateTTbarObservables(const Muon& lepton,
                                                                   const Jet& b1, const Jet& b2,
                                                                   const Jet& j1, const Jet& j2,
                                                                   float met_pt, float met_phi) {
    TTbarObservables obs;
    
    // Calculate both versions of ttbar transverse mass
    obs.mt_ttbar_v1 = CalculateTTbarSystemTransverseMass(lepton, b1, b2, j1, j2, met_pt, met_phi);
    obs.mt_ttbar_v2 = CalculateTTbarSystemTransverseMassV2(lepton, b1, b2, j1, j2, met_pt, met_phi);
    
    // W boson transverse mass
    obs.mt_W = sqrt(2.0 * lepton.Pt() * met_pt * (1.0 - cos(lepton.Phi() - met_phi)));
    
    // Hadronic system mass (assuming b1 + j1 + j2 or b2 + j1 + j2)
    Particle hadronic_system = b1 + j1 + j2;
    obs.m_hadronic = hadronic_system.M();
    
    // Total visible mass
    Particle visible_system = lepton + b1 + b2 + j1 + j2;
    obs.m_visible = visible_system.M();
    
    // Total HT (scalar sum of transverse momenta)
    obs.HT = lepton.Pt() + b1.Pt() + b2.Pt() + j1.Pt() + j2.Pt() + met_pt;
    
    // Missing HT ratio
    float visible_ht = lepton.Pt() + b1.Pt() + b2.Pt() + j1.Pt() + j2.Pt();
    obs.MET_over_HT = met_pt / visible_ht;
    
    // Total visible pT
    float vis_px = lepton.Px() + b1.Px() + b2.Px() + j1.Px() + j2.Px();
    float vis_py = lepton.Py() + b1.Py() + b2.Py() + j1.Py() + j2.Py();
    obs.pt_visible = sqrt(vis_px*vis_px + vis_py*vis_py);
    
    // MET vs visible pT balance
    obs.MET_vis_balance = (met_pt - obs.pt_visible) / (met_pt + obs.pt_visible);
    
    return obs;
}

float TTbar_test::CalculateTTbarSystemTransverseMass(const Muon& lepton,
                                                    const Jet& b1, const Jet& b2,
                                                    const Jet& j1, const Jet& j2,
                                                    float met_pt, float met_phi) {
    // Calculate total visible transverse momentum components
    float total_vis_px = lepton.Px() + b1.Px() + b2.Px() + j1.Px() + j2.Px();
    float total_vis_py = lepton.Py() + b1.Py() + b2.Py() + j1.Py() + j2.Py();
    
    // Add MET components
    float met_px = met_pt * cos(met_phi);
    float met_py = met_pt * sin(met_phi);
    
    float total_px = total_vis_px + met_px;
    float total_py = total_vis_py + met_py;
    
    // Calculate total transverse momentum magnitude
    float total_pt = sqrt(total_px*total_px + total_py*total_py);
    
    // Calculate total visible energy and longitudinal momentum
    float total_vis_energy = lepton.E() + b1.E() + b2.E() + j1.E() + j2.E();
    float total_vis_pz = lepton.Pz() + b1.Pz() + b2.Pz() + j1.Pz() + j2.Pz();
    
    float total_vis_pt = sqrt(total_vis_px*total_vis_px + total_vis_py*total_vis_py);
    
    // Transverse mass formula for the complete system
    float mt_ttbar_squared = total_vis_energy*total_vis_energy - total_vis_pz*total_vis_pz + 
                            total_pt*total_pt - total_vis_pt*total_vis_pt;
    
    // Handle negative values
    return sqrt(max(mt_ttbar_squared, 0.0f));
}

float TTbar_test::CalculateTTbarSystemTransverseMassV2(const Muon& lepton,
                                                      const Jet& b1, const Jet& b2,
                                                      const Jet& j1, const Jet& j2,
                                                      float met_pt, float met_phi) {
    // Calculate transverse energies for all visible particles
    float et_lepton = lepton.Et();
    float et_b1 = b1.Et();
    float et_b2 = b2.Et();
    float et_j1 = j1.Et();
    float et_j2 = j2.Et();
    float et_met = met_pt; // MET is already transverse, assume massless
    
    // Sum all transverse energies
    float et_total = et_lepton + et_b1 + et_b2 + et_j1 + et_j2 + et_met;
    
    // Calculate total transverse momentum vector
    float px_total = lepton.Px() + b1.Px() + b2.Px() + j1.Px() + j2.Px() + met_pt * cos(met_phi);
    float py_total = lepton.Py() + b1.Py() + b2.Py() + j1.Py() + j2.Py() + met_pt * sin(met_phi);
    
    float pt_total = sqrt(px_total*px_total + py_total*py_total);
    
    // Transverse mass: mt = sqrt(Et_total^2 - pt_total^2)
    return sqrt(et_total*et_total - pt_total*pt_total);
}

