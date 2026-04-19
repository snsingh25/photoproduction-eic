// =============================================================================
// Thick and Thin Jet Classification Analysis for HERA Photoproduction Events
// Enhanced version with jet constituent particle tracking
// =============================================================================
// Description: Reconstructs jets from photoproduction events and classifies them
//              as "thick" (gluon-like) or "thin" (quark-like) based on integrated
//              jet shape variables, following HERA methodology.
//              Also saves all constituent particles for each jet.
//
// Input: ROOT file with photoproduction events from evtgen.cc
// Output: ROOT file with thick and thin jets in separate directories
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TObjArray.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <memory>
#include <cmath>

using namespace std;
using namespace fastjet;

// =============================================================================
// ANALYSIS CUTS AND PARAMETERS
// =============================================================================

// Jet clustering parameters
const double R = 1.0;                    // Jet radius
const double etMin = 17.0;                // Minimum ET for jets
const double etaMin = -4.0;               // Minimum eta
const double etaMax = 4.0;                // Maximum eta

// Particle selection cuts
const double particle_pt_min = 0.1;      // Minimum particle pT (GeV)
const double particle_eta_max = 5.0;     // Maximum particle |eta|

// Jet shape classification cuts (following HERA methodology)
const double thin_jet_cut = 0.8;         // Ψ(r = 0.3) > 0.8 for thin jets (quark-like)
const double thick_jet_cut = 0.6;        // Ψ(r = 0.3) < 0.6 for thick jets (gluon-like)
const double jet_shape_radius = 0.3;     // Radius for integrated jet shape calculation

// Event progress reporting
const int report_every = 10000;          // Report progress every N events

// Output file name
const string output_filename = "thick_thin_jets_with_constituents.root";

// =============================================================================
// INTEGRATED JET SHAPE CALCULATION FUNCTION
// =============================================================================

double calculateIntegratedJetShape(const PseudoJet& jet, double radius, 
                                 const ClusterSequence& cs) {
    // Get jet constituents using the cluster sequence
    vector<PseudoJet> constituents = cs.constituents(jet);
    
    if (constituents.empty()) return 0.0;
    
    double jet_eta = jet.eta();
    double jet_phi = jet.phi();
    double total_et = jet.Et();
    
    double et_within_radius = 0.0;
    
    // Calculate ET within the specified radius
    for (const auto& constituent : constituents) {
        double delta_eta = constituent.eta() - jet_eta;
        double delta_phi = constituent.phi() - jet_phi;
        
        // Handle phi wrap-around
        while (delta_phi > M_PI) delta_phi -= 2*M_PI;
        while (delta_phi < -M_PI) delta_phi += 2*M_PI;
        
        double delta_r = sqrt(delta_eta*delta_eta + delta_phi*delta_phi);
        
        if (delta_r <= radius) {
            et_within_radius += constituent.Et();
        }
    }
    
    // Return integrated jet shape: fraction of ET within radius
    return (total_et > 0) ? et_within_radius / total_et : 0.0;
}

// =============================================================================
// JET RECONSTRUCTION FUNCTION WITH CLUSTER SEQUENCE PRESERVATION
// =============================================================================

struct JetResult {
    vector<PseudoJet> jets;
    shared_ptr<ClusterSequence> cs;
};

JetResult reconstructJetsWithCS(const vector<Float_t>& px, const vector<Float_t>& py, 
                               const vector<Float_t>& pz, const vector<Float_t>& energy,
                               const vector<Float_t>& eta, const vector<Int_t>& pdgId) {
    
    JetResult result;
    
    // Check input vector sizes
    if (px.empty() || px.size() != py.size() || px.size() != pz.size() || 
        px.size() != energy.size() || px.size() != eta.size() || px.size() != pdgId.size()) {
        return result;
    }
    
    // Create input particles for FastJet
    vector<PseudoJet> input_particles;
    
    for (size_t i = 0; i < px.size(); ++i) {
        // Apply particle-level cuts
        double pt = sqrt(px[i]*px[i] + py[i]*py[i]);
        if (pt < particle_pt_min) continue;
        if (abs(eta[i]) > particle_eta_max) continue;
        
        // Check for valid energy and momentum
        if (energy[i] <= 0) continue;
        
        // Create PseudoJet
        PseudoJet particle(px[i], py[i], pz[i], energy[i]);
        particle.set_user_index(i);  // Store original index
        input_particles.push_back(particle);
    }
    
    if (input_particles.size() < 2) {
        return result;  // Return empty result if too few particles
    }
    
    // Define longitudinally invariant kT algorithm (matching HERA methodology)
    JetDefinition jet_def(kt_algorithm, R);
    
    // Create cluster sequence and store it
    result.cs = make_shared<ClusterSequence>(input_particles, jet_def);
    
    // Get jets and apply cuts
    vector<PseudoJet> jets = result.cs->inclusive_jets();
    
    for (const auto& jet : jets) {
        double et = jet.Et();
        double jet_eta = jet.eta();
        
        // Apply jet-level cuts
        if (et > etMin && jet_eta > etaMin && jet_eta < etaMax) {
            result.jets.push_back(jet);
        }
    }
    
    // Sort jets by ET (descending)
    sort(result.jets.begin(), result.jets.end(), 
         [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
    
    return result;
}

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================

int main() {
    
    cout << "=============================================================================\n";
    cout << "THICK AND THIN JET CLASSIFICATION ANALYSIS FOR HERA PHOTOPRODUCTION\n";
    cout << "WITH JET CONSTITUENT TRACKING\n";
    cout << "=============================================================================\n";
    cout << "Analysis Parameters:\n";
    cout << "  Jet Algorithm:        Longitudinally Invariant kT\n";
    cout << "  Jet Radius:           " << R << "\n";
    cout << "  Min Jet ET:           " << etMin << " GeV\n";
    cout << "  Jet Eta Range:        [" << etaMin << ", " << etaMax << "]\n";
    cout << "  Min Particle pT:      " << particle_pt_min << " GeV\n";
    cout << "  Jet Shape Radius:     " << jet_shape_radius << "\n";
    cout << "  Thin Jet Cut:         Ψ(r=" << jet_shape_radius << ") > " << thin_jet_cut << "\n";
    cout << "  Thick Jet Cut:        Ψ(r=" << jet_shape_radius << ") < " << thick_jet_cut << "\n";
    cout << "=============================================================================\n\n";
    
    // Open input file
    const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt7GeV/hera1M/hera.root";
    TFile* input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file || input_file->IsZombie()) {
        cerr << "ERROR: Cannot open input file: " << input_filename << endl;
        return 1;
    }
    
    cout << "Successfully opened input file: " << input_filename << "\n";
    
    // Get the Combined_Events tree (contains all events)
    TTree* input_tree = (TTree*)input_file->Get("Combined_Events/Combined_Events");
    if (!input_tree) {
        cerr << "ERROR: Cannot find Combined_Events tree in input file" << endl;
        cerr << "Looking for alternative tree locations..." << endl;
        
        // Try other possible locations
        const char* tree_paths[] = {
            "Combined_Events",
            "QQ_Events/QQ_Events", 
            "eventdata",
            "TC"
        };
        
        for (const char* path : tree_paths) {
            input_tree = (TTree*)input_file->Get(path);
            if (input_tree) {
                cout << "Found tree at: " << path << "\n";
                break;
            }
        }
        
        if (!input_tree) {
            cerr << "ERROR: Could not find any suitable event tree" << endl;
            return 1;
        }
    }
    
    Long64_t n_entries = input_tree->GetEntries();
    cout << "Found " << n_entries << " events in tree\n\n";
    
    if (n_entries == 0) {
        cout << "No events found, exiting...\n";
        return 1;
    }
    
    // Set up input branches - using the structure from evtgen.cc
    vector<Float_t> *px = nullptr, *py = nullptr, *pz = nullptr, *energy = nullptr, *eta = nullptr;
    vector<Float_t> *phi = nullptr, *pT = nullptr, *mass = nullptr;
    vector<Int_t> *pdgId = nullptr, *status = nullptr;
    vector<Bool_t> *isCharged = nullptr, *isHadron = nullptr;
    Int_t eventID;
    Float_t crossSection;
    Bool_t isResolved, isDirect;
    
    // Set branch addresses
    input_tree->SetBranchAddress("eventID", &eventID);
    input_tree->SetBranchAddress("px", &px);
    input_tree->SetBranchAddress("py", &py);
    input_tree->SetBranchAddress("pz", &pz);
    input_tree->SetBranchAddress("energy", &energy);
    input_tree->SetBranchAddress("eta", &eta);
    input_tree->SetBranchAddress("phi", &phi);
    input_tree->SetBranchAddress("pT", &pT);
    input_tree->SetBranchAddress("mass", &mass);
    input_tree->SetBranchAddress("pdgId", &pdgId);
    input_tree->SetBranchAddress("status", &status);
    input_tree->SetBranchAddress("isCharged", &isCharged);
    input_tree->SetBranchAddress("isHadron", &isHadron);
    input_tree->SetBranchAddress("crossSection", &crossSection);
    input_tree->SetBranchAddress("isResolved", &isResolved);
    input_tree->SetBranchAddress("isDirect", &isDirect);
    
    cout << "Successfully set up input branches\n";
    
    // Create output file
    TFile* output_file = new TFile(output_filename.c_str(), "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        cerr << "ERROR: Cannot create output file: " << output_filename << endl;
        return 1;
    }
    
    // Create directories for thick and thin jets
    TDirectory* thick_dir = output_file->mkdir("ThickJets");
    TDirectory* thin_dir = output_file->mkdir("ThinJets");
    
    // Create trees for thick jets
    thick_dir->cd();
    TTree* thick_tree = new TTree("thick_jets", "Thick (Gluon-like) Jets with Constituents");
    
    // Variables for thick jets
    Int_t thick_eventID, thick_n_jets;
    Float_t thick_event_weight;
    Bool_t thick_isResolved, thick_isDirect;
    vector<Float_t> thick_jet_mass, thick_jet_pt, thick_jet_px, thick_jet_py, thick_jet_pz;
    vector<Float_t> thick_jet_eta, thick_jet_phi, thick_jet_et, thick_jet_shape;
    vector<Float_t> thick_jet_energy;
    vector<Int_t> thick_jet_nconst;
    
    // NEW: Jet constituent information for thick jets
    vector<vector<Float_t>> thick_const_px, thick_const_py, thick_const_pz, thick_const_energy;
    vector<vector<Float_t>> thick_const_eta, thick_const_phi, thick_const_pt, thick_const_mass;
    vector<vector<Int_t>> thick_const_pdgId, thick_const_status;
    vector<vector<Bool_t>> thick_const_isCharged, thick_const_isHadron;
    vector<vector<Float_t>> thick_const_deltaR;  // Distance from jet axis
    
    thick_tree->Branch("eventID", &thick_eventID);
    thick_tree->Branch("n_jets", &thick_n_jets);
    thick_tree->Branch("event_weight", &thick_event_weight);
    thick_tree->Branch("isResolved", &thick_isResolved);
    thick_tree->Branch("isDirect", &thick_isDirect);
    thick_tree->Branch("jet_mass", &thick_jet_mass);
    thick_tree->Branch("jet_pt", &thick_jet_pt);
    thick_tree->Branch("jet_px", &thick_jet_px);
    thick_tree->Branch("jet_py", &thick_jet_py);
    thick_tree->Branch("jet_pz", &thick_jet_pz);
    thick_tree->Branch("jet_eta", &thick_jet_eta);
    thick_tree->Branch("jet_phi", &thick_jet_phi);
    thick_tree->Branch("jet_et", &thick_jet_et);
    thick_tree->Branch("jet_energy", &thick_jet_energy);
    thick_tree->Branch("jet_shape", &thick_jet_shape);
    thick_tree->Branch("jet_nconst", &thick_jet_nconst);
    
    // NEW: Constituent branches for thick jets
    thick_tree->Branch("const_px", &thick_const_px);
    thick_tree->Branch("const_py", &thick_const_py);
    thick_tree->Branch("const_pz", &thick_const_pz);
    thick_tree->Branch("const_energy", &thick_const_energy);
    thick_tree->Branch("const_eta", &thick_const_eta);
    thick_tree->Branch("const_phi", &thick_const_phi);
    thick_tree->Branch("const_pt", &thick_const_pt);
    thick_tree->Branch("const_mass", &thick_const_mass);
    thick_tree->Branch("const_pdgId", &thick_const_pdgId);
    thick_tree->Branch("const_status", &thick_const_status);
    thick_tree->Branch("const_isCharged", &thick_const_isCharged);
    thick_tree->Branch("const_isHadron", &thick_const_isHadron);
    thick_tree->Branch("const_deltaR", &thick_const_deltaR);
    
    // Create trees for thin jets
    thin_dir->cd();
    TTree* thin_tree = new TTree("thin_jets", "Thin (Quark-like) Jets with Constituents");
    
    // Variables for thin jets
    Int_t thin_eventID, thin_n_jets;
    Float_t thin_event_weight;
    Bool_t thin_isResolved, thin_isDirect;
    vector<Float_t> thin_jet_mass, thin_jet_pt, thin_jet_px, thin_jet_py, thin_jet_pz;
    vector<Float_t> thin_jet_eta, thin_jet_phi, thin_jet_et, thin_jet_shape;
    vector<Float_t> thin_jet_energy;
    vector<Int_t> thin_jet_nconst;
    
    // NEW: Jet constituent information for thin jets
    vector<vector<Float_t>> thin_const_px, thin_const_py, thin_const_pz, thin_const_energy;
    vector<vector<Float_t>> thin_const_eta, thin_const_phi, thin_const_pt, thin_const_mass;
    vector<vector<Int_t>> thin_const_pdgId, thin_const_status;
    vector<vector<Bool_t>> thin_const_isCharged, thin_const_isHadron;
    vector<vector<Float_t>> thin_const_deltaR;  // Distance from jet axis
    
    thin_tree->Branch("eventID", &thin_eventID);
    thin_tree->Branch("n_jets", &thin_n_jets);
    thin_tree->Branch("event_weight", &thin_event_weight);
    thin_tree->Branch("isResolved", &thin_isResolved);
    thin_tree->Branch("isDirect", &thin_isDirect);
    thin_tree->Branch("jet_mass", &thin_jet_mass);
    thin_tree->Branch("jet_pt", &thin_jet_pt);
    thin_tree->Branch("jet_px", &thin_jet_px);
    thin_tree->Branch("jet_py", &thin_jet_py);
    thin_tree->Branch("jet_pz", &thin_jet_pz);
    thin_tree->Branch("jet_eta", &thin_jet_eta);
    thin_tree->Branch("jet_phi", &thin_jet_phi);
    thin_tree->Branch("jet_et", &thin_jet_et);
    thin_tree->Branch("jet_energy", &thin_jet_energy);
    thin_tree->Branch("jet_shape", &thin_jet_shape);
    thin_tree->Branch("jet_nconst", &thin_jet_nconst);
    
    // NEW: Constituent branches for thin jets
    thin_tree->Branch("const_px", &thin_const_px);
    thin_tree->Branch("const_py", &thin_const_py);
    thin_tree->Branch("const_pz", &thin_const_pz);
    thin_tree->Branch("const_energy", &thin_const_energy);
    thin_tree->Branch("const_eta", &thin_const_eta);
    thin_tree->Branch("const_phi", &thin_const_phi);
    thin_tree->Branch("const_pt", &thin_const_pt);
    thin_tree->Branch("const_mass", &thin_const_mass);
    thin_tree->Branch("const_pdgId", &thin_const_pdgId);
    thin_tree->Branch("const_status", &thin_const_status);
    thin_tree->Branch("const_isCharged", &thin_const_isCharged);
    thin_tree->Branch("const_isHadron", &thin_const_isHadron);
    thin_tree->Branch("const_deltaR", &thin_const_deltaR);
    
    // Event counters
    int events_processed = 0;
    int events_with_jets = 0;
    int total_jets = 0;
    int thick_jets_count = 0;
    int thin_jets_count = 0;
    int intermediate_jets_count = 0;
    
    cout << "Starting event processing...\n";
    
    // Process events
    for (Long64_t i = 0; i < n_entries; ++i) {
        
        if (i % report_every == 0) {
            cout << "Processing event " << i << "/" << n_entries 
                 << " (" << (100*i/n_entries) << "%)\r" << flush;
        }
        
        input_tree->GetEntry(i);
        
        // Check if vectors are valid and not empty
        if (!px || !py || !pz || !energy || !eta || !pdgId) {
            continue;
        }
        
        // Check if all vectors have the same size
        if (px->size() != py->size() || px->size() != pz->size() || 
            px->size() != energy->size() || px->size() != eta->size() ||
            px->size() != pdgId->size()) {
            continue;
        }
        
        events_processed++;
        
        // Clear output vectors for jets
        thick_jet_mass.clear(); thick_jet_pt.clear(); thick_jet_px.clear(); 
        thick_jet_py.clear(); thick_jet_pz.clear(); thick_jet_eta.clear(); 
        thick_jet_phi.clear(); thick_jet_et.clear(); thick_jet_shape.clear();
        thick_jet_energy.clear(); thick_jet_nconst.clear();
        
        thin_jet_mass.clear(); thin_jet_pt.clear(); thin_jet_px.clear(); 
        thin_jet_py.clear(); thin_jet_pz.clear(); thin_jet_eta.clear(); 
        thin_jet_phi.clear(); thin_jet_et.clear(); thin_jet_shape.clear();
        thin_jet_energy.clear(); thin_jet_nconst.clear();
        
        // Clear output vectors for constituents
        thick_const_px.clear(); thick_const_py.clear(); thick_const_pz.clear(); thick_const_energy.clear();
        thick_const_eta.clear(); thick_const_phi.clear(); thick_const_pt.clear(); thick_const_mass.clear();
        thick_const_pdgId.clear(); thick_const_status.clear();
        thick_const_isCharged.clear(); thick_const_isHadron.clear(); thick_const_deltaR.clear();
        
        thin_const_px.clear(); thin_const_py.clear(); thin_const_pz.clear(); thin_const_energy.clear();
        thin_const_eta.clear(); thin_const_phi.clear(); thin_const_pt.clear(); thin_const_mass.clear();
        thin_const_pdgId.clear(); thin_const_status.clear();
        thin_const_isCharged.clear(); thin_const_isHadron.clear(); thin_const_deltaR.clear();
        
        // Reconstruct jets with preserved cluster sequence
        JetResult jet_result = reconstructJetsWithCS(*px, *py, *pz, *energy, *eta, *pdgId);
        vector<PseudoJet>& jets = jet_result.jets;
        
        total_jets += jets.size();
        
        if (jets.size() > 0) {
            events_with_jets++;
            
            int thick_jets_this_event = 0;
            int thin_jets_this_event = 0;
            
            // Process each jet and classify as thick or thin
            for (const auto& jet : jets) {
                
                // Calculate integrated jet shape at r = 0.3
                double jet_shape = calculateIntegratedJetShape(jet, jet_shape_radius, *jet_result.cs);
                
                // Get jet constituents
                vector<PseudoJet> constituents = jet_result.cs->constituents(jet);
                int n_constituents = constituents.size();
                
                // Prepare constituent vectors for this jet
                vector<Float_t> const_px_jet, const_py_jet, const_pz_jet, const_energy_jet;
                vector<Float_t> const_eta_jet, const_phi_jet, const_pt_jet, const_mass_jet;
                vector<Int_t> const_pdgId_jet, const_status_jet;
                vector<Bool_t> const_isCharged_jet, const_isHadron_jet;
                vector<Float_t> const_deltaR_jet;
                
                // Fill constituent information
                for (const auto& constituent : constituents) {
                    int orig_index = constituent.user_index();
                    
                    // Store constituent 4-momentum
                    const_px_jet.push_back(constituent.px());
                    const_py_jet.push_back(constituent.py());
                    const_pz_jet.push_back(constituent.pz());
                    const_energy_jet.push_back(constituent.e());
                    const_eta_jet.push_back(constituent.eta());
                    const_phi_jet.push_back(constituent.phi());
                    const_pt_jet.push_back(constituent.pt());
                    const_mass_jet.push_back(constituent.m());
                    
                    // Calculate deltaR from jet axis
                    double delta_eta = constituent.eta() - jet.eta();
                    double delta_phi = constituent.phi() - jet.phi();
                    while (delta_phi > M_PI) delta_phi -= 2*M_PI;
                    while (delta_phi < -M_PI) delta_phi += 2*M_PI;
                    double deltaR = sqrt(delta_eta*delta_eta + delta_phi*delta_phi);
                    const_deltaR_jet.push_back(deltaR);
                    
                    // Store original particle properties if available
                    if (orig_index >= 0 && orig_index < (int)pdgId->size()) {
                        const_pdgId_jet.push_back((*pdgId)[orig_index]);
                        const_status_jet.push_back(status ? (*status)[orig_index] : 1);
                        const_isCharged_jet.push_back(isCharged ? (*isCharged)[orig_index] : false);
                        const_isHadron_jet.push_back(isHadron ? (*isHadron)[orig_index] : false);
                    } else {
                        const_pdgId_jet.push_back(0);
                        const_status_jet.push_back(1);
                        const_isCharged_jet.push_back(false);
                        const_isHadron_jet.push_back(false);
                    }
                }
                
                // Classify jet based on integrated jet shape
                if (jet_shape > thin_jet_cut) {
                    // Thin jet (quark-like)
                    thin_jet_mass.push_back(jet.m());
                    thin_jet_pt.push_back(jet.pt());
                    thin_jet_px.push_back(jet.px());
                    thin_jet_py.push_back(jet.py());
                    thin_jet_pz.push_back(jet.pz());
                    thin_jet_eta.push_back(jet.eta());
                    thin_jet_phi.push_back(jet.phi());
                    thin_jet_et.push_back(jet.Et());
                    thin_jet_energy.push_back(jet.e());
                    thin_jet_shape.push_back(jet_shape);
                    thin_jet_nconst.push_back(n_constituents);
                    
                    // Store constituent information for this thin jet
                    thin_const_px.push_back(const_px_jet);
                    thin_const_py.push_back(const_py_jet);
                    thin_const_pz.push_back(const_pz_jet);
                    thin_const_energy.push_back(const_energy_jet);
                    thin_const_eta.push_back(const_eta_jet);
                    thin_const_phi.push_back(const_phi_jet);
                    thin_const_pt.push_back(const_pt_jet);
                    thin_const_mass.push_back(const_mass_jet);
                    thin_const_pdgId.push_back(const_pdgId_jet);
                    thin_const_status.push_back(const_status_jet);
                    thin_const_isCharged.push_back(const_isCharged_jet);
                    thin_const_isHadron.push_back(const_isHadron_jet);
                    thin_const_deltaR.push_back(const_deltaR_jet);
                    
                    thin_jets_this_event++;
                    thin_jets_count++;
                }
                else if (jet_shape < thick_jet_cut) {
                    // Thick jet (gluon-like)
                    thick_jet_mass.push_back(jet.m());
                    thick_jet_pt.push_back(jet.pt());
                    thick_jet_px.push_back(jet.px());
                    thick_jet_py.push_back(jet.py());
                    thick_jet_pz.push_back(jet.pz());
                    thick_jet_eta.push_back(jet.eta());
                    thick_jet_phi.push_back(jet.phi());
                    thick_jet_et.push_back(jet.Et());
                    thick_jet_energy.push_back(jet.e());
                    thick_jet_shape.push_back(jet_shape);
                    thick_jet_nconst.push_back(n_constituents);
                    
                    // Store constituent information for this thick jet
                    thick_const_px.push_back(const_px_jet);
                    thick_const_py.push_back(const_py_jet);
                    thick_const_pz.push_back(const_pz_jet);
                    thick_const_energy.push_back(const_energy_jet);
                    thick_const_eta.push_back(const_eta_jet);
                    thick_const_phi.push_back(const_phi_jet);
                    thick_const_pt.push_back(const_pt_jet);
                    thick_const_mass.push_back(const_mass_jet);
                    thick_const_pdgId.push_back(const_pdgId_jet);
                    thick_const_status.push_back(const_status_jet);
                    thick_const_isCharged.push_back(const_isCharged_jet);
                    thick_const_isHadron.push_back(const_isHadron_jet);
                    thick_const_deltaR.push_back(const_deltaR_jet);
                    
                    thick_jets_this_event++;
                    thick_jets_count++;
                }
                else {
                    // Intermediate jet (not clearly thick or thin)
                    intermediate_jets_count++;
                }
            }
            
            // Fill thick jets tree if there are thick jets
            if (thick_jets_this_event > 0) {
                thick_eventID = eventID;
                thick_n_jets = thick_jets_this_event;
                thick_event_weight = crossSection;
                thick_isResolved = isResolved;
                thick_isDirect = isDirect;
                thick_tree->Fill();
            }
            
            // Fill thin jets tree if there are thin jets
            if (thin_jets_this_event > 0) {
                thin_eventID = eventID;
                thin_n_jets = thin_jets_this_event;
                thin_event_weight = crossSection;
                thin_isResolved = isResolved;
                thin_isDirect = isDirect;
                thin_tree->Fill();
            }
        }
    }
    
    cout << "\n\n=============================================================================\n";
    cout << "ANALYSIS SUMMARY\n";
    cout << "=============================================================================\n";
    cout << "Events processed: " << events_processed << "\n";
    cout << "Events with jets: " << events_with_jets 
         << " (" << (100.0*events_with_jets/events_processed) << "%)\n";
    cout << "Total jets found: " << total_jets << "\n";
    cout << "Thick jets (gluon-like): " << thick_jets_count 
         << " (" << (100.0*thick_jets_count/total_jets) << "%)\n";
    cout << "Thin jets (quark-like): " << thin_jets_count 
         << " (" << (100.0*thin_jets_count/total_jets) << "%)\n";
    cout << "Intermediate jets: " << intermediate_jets_count 
         << " (" << (100.0*intermediate_jets_count/total_jets) << "%)\n";
    cout << "Average jets/event: " << (events_with_jets > 0 ? (double)total_jets/events_with_jets : 0) << "\n";
    
    // Write trees
    thick_dir->cd();
    thick_tree->Write();
    
    thin_dir->cd();
    thin_tree->Write();
    
    // Close files
    output_file->Write();
    output_file->Close();
    input_file->Close();
    
    cout << "\n=============================================================================\n";
    cout << "ANALYSIS COMPLETED SUCCESSFULLY!\n";
    cout << "=============================================================================\n";
    cout << "Output file: " << output_filename << "\n";
    cout << "\nOutput structure:\n";
    cout << "  /ThickJets/thick_jets - Tree with thick (gluon-like) jets\n";
    cout << "  /ThinJets/thin_jets   - Tree with thin (quark-like) jets\n";
    cout << "\nJet variables saved:\n";
    cout << "  - Kinematics: mass, pt, px, py, pz, eta, phi, et, energy\n";
    cout << "  - Shape: jet_shape (integrated shape at r=0.3)\n";
    cout << "  - Structure: jet_nconst (number of constituents)\n";
    cout << "  - Event info: eventID, event_weight, isResolved, isDirect\n";
    cout << "\nNEW: Constituent particle information saved:\n";
    cout << "  - const_px, const_py, const_pz, const_energy (4-momentum)\n";
    cout << "  - const_eta, const_phi, const_pt, const_mass (kinematics)\n";
    cout << "  - const_pdgId, const_status (particle identification)\n";
    cout << "  - const_isCharged, const_isHadron (particle properties)\n";
    cout << "  - const_deltaR (distance from jet axis)\n";
    cout << "\nData structure:\n";
    cout << "  - Each jet variable is a vector with n_jets entries\n";
    cout << "  - Each constituent variable is a vector<vector<T>> where:\n";
    cout << "    * Outer index = jet number (0 to n_jets-1)\n";
    cout << "    * Inner index = constituent number within that jet\n";
    cout << "\nExample usage in analysis:\n";
    cout << "  for (int ijet = 0; ijet < n_jets; ijet++) {\n";
    cout << "    for (int iconst = 0; iconst < const_px[ijet].size(); iconst++) {\n";
    cout << "      // Access constituent iconst of jet ijet:\n";
    cout << "      // const_px[ijet][iconst], const_py[ijet][iconst], etc.\n";
    cout << "    }\n";
    cout << "  }\n";
    cout << "\nClassification criteria:\n";
    cout << "  - Thin jets (quark-like): Ψ(r=0.3) > 0.8\n";
    cout << "  - Thick jets (gluon-like): Ψ(r=0.3) < 0.6\n";
    cout << "  - Intermediate jets: 0.6 ≤ Ψ(r=0.3) ≤ 0.8\n";
    cout << "=============================================================================\n";
    
    return 0;
}