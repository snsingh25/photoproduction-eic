// =============================================================================
// Jet Reconstruction Analysis for HERA Photoproduction Events
// =============================================================================
// Description: Reconstructs jets from different subprocess categories using 
//              anti-kT algorithm and saves basic jet properties with dijet cuts.
//              Also computes softdrop multiplicity (n_sd) per jet using
//              modified mass drop tagger (zcut=0.1, beta=0) via Cambridge/Aachen
//              reclustering of constituents.
//
// Input: ROOT file with photoproduction events from improved event generator
// Output: ROOT file with reconstructed dijet events for each process category
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
#include <cmath>
#include <fstream>

using namespace std;
using namespace fastjet;

// =============================================================================
// DUAL OUTPUT CLASS FOR CONSOLE AND FILE LOGGING
// =============================================================================
class DualOutput {
private:
    std::ofstream logFile;
    bool fileOpen;
    
public:
    DualOutput(const std::string& filename) : fileOpen(false) {
        logFile.open(filename.c_str());
        fileOpen = logFile.is_open();
        if (!fileOpen) {
            std::cerr << "Warning: Could not open log file: " << filename << std::endl;
        }
    }
    
    ~DualOutput() {
        if (fileOpen) {
            logFile.close();
        }
    }
    
    template<typename T>
    DualOutput& operator<<(const T& data) {
        std::cout << data;  // Output to console
        if (fileOpen) {
            logFile << data;  // Output to file
            logFile.flush();  // Ensure immediate write
        }
        return *this;
    }
    
    // Handle special cases like std::endl, std::flush
    DualOutput& operator<<(std::ostream& (*manip)(std::ostream&)) {
        std::cout << manip;
        if (fileOpen) {
            logFile << manip;
            logFile.flush();
        }
        return *this;
    }
};

// =============================================================================
// ANALYSIS CONFIGURATION
// =============================================================================

// Analysis mode selection
const bool DIJET_ONLY = false;    // Set to true for dijet analysis, false for all jets

// Jet clustering parameters 
const double R = 1.0;             // Jet radius
const double etMin = 17.0;         // Minimum ET for any jet (HERA = 17, EIC = 10), ALL jets must pass this Base Cut
const double etaMin = -4.0;        // Minimum eta for any jet
const double etaMax = 4.0;         // Maximum eta for any jet

// Dijet-specific cuts (only used when DIJET_ONLY = true)
const double leading_jet_et_min = 10.0;      // Minimum ET for leading jet (GeV)
const double subleading_jet_et_min = 7.0;    // Minimum ET for subleading jet (GeV)
const double leading_jet_eta_max = 4.5;      // Maximum |eta| for leading jet
const double subleading_jet_eta_max = 4.5;   // Maximum |eta| for subleading jet

// Additional dijet cuts
// const double dijet_mass_min = 0.0;          // Minimum dijet invariant mass (GeV)
// const double delta_phi_min = 2.5;            // Minimum azimuthal separation (radians)
// const double delta_eta_max = 4.0;            // Maximum rapidity separation

// Particle selection cuts
const double particle_pt_min = 0.1;     // Minimum particle pT (GeV)
const double particle_eta_max = 5.0;    // Maximum particle |eta|

// Softdrop multiplicity parameters (modified mass drop tagger)
// See Frye, Larkoski, Moult, Thaler, JHEP 07 (2017) 064 [arXiv:1704.06266]
const double sd_zcut = 0.1;       // Energy-sharing cut
const double sd_beta = 0.0;       // Angular exponent (0 = modified mass drop)
const double sd_R0   = 1.0;       // Normalization radius (match jet R)

// n_subjets parameters (kT exclusive reclustering of constituents)
// Same ycut as existing nsubjets.cc analysis, for head-to-head comparison.
const double kt_ycut = 0.0005;    // y_cut for exclusive_jets_ycut

// Event progress reporting
const int report_every = 100000;    // Report progres/s every N events

// Input Events — default; override via argv[1]. Points at the repo's data/.
const string default_input_filename =
    "/Users/siddharthsingh/Analysis/repos/photoproduction-eic/data/"
    "allevents_pt7GeV/hera300_pT7/hera300_pT7.root";
// Filled in main() once argv is parsed.
string input_filename;

// =============================================================================
// OUTPUT FILENAME GENERATION
// =============================================================================
string generateOutputFilename(const string& input_filename, double R, double etMin) {
    string baseFileName = input_filename.substr(input_filename.find_last_of("/") + 1);
    baseFileName = baseFileName.substr(0, baseFileName.find_last_of("."));  // Remove .root extension
    string mode_str, leadjet_str, subleadjet_str;
    if (DIJET_ONLY) {
        mode_str = "dijets";
        leadjet_str = to_string((int)leading_jet_et_min);
        subleadjet_str = to_string((int)subleading_jet_et_min);
        return mode_str + "_" + baseFileName + "_R" + to_string((int)(R * 10)) + "_EtMin" + to_string((int)etMin) + "_" + leadjet_str + "_" + subleadjet_str + ".root";
    } else {
        mode_str = "alljets";
        return mode_str + "_" + baseFileName + "_R" + to_string((int)(R * 10)) + "_EtMin" + to_string((int)etMin) + ".root";
    }
}

string generateLogFilename(const string& root_filename) {
    string log_filename = root_filename;
    // Replace .root with .log
    size_t pos = log_filename.find_last_of(".");
    if (pos != string::npos) {
        log_filename = log_filename.substr(0, pos) + ".log";
    } else {
        log_filename += ".log";
    }
    return log_filename;
}

// output_filename / log_filename are derived from input_filename in main().
string output_filename;
string log_filename;

// =============================================================================
// JET RECONSTRUCTION FUNCTION
// =============================================================================
vector<PseudoJet> reconstructJets(const vector<float>& px, const vector<float>& py,
                                  const vector<float>& pz, const vector<float>& energy,
                                  const vector<float>& eta,
                                  vector<int>& jet_nsd_out,
                                  vector<int>& jet_nsubjets_out,
                                  vector<float>& jet_psi03_out,
                                  vector<int>& jet_nsd_beta1_out,
                                  vector<int>& jet_nsd_loose_out,
                                  JetAlgorithm jet_algo_in = antikt_algorithm) {

    // Always start with an empty output so early returns leave it in a sane state
    jet_nsd_out.clear();
    jet_nsubjets_out.clear();
    jet_psi03_out.clear();
    jet_nsd_beta1_out.clear();
    jet_nsd_loose_out.clear();
    
    // Check input vector sizes
    if (px.empty() || px.size() != py.size() || px.size() != pz.size() || 
        px.size() != energy.size() || px.size() != eta.size()) {
        return vector<PseudoJet>();
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
        return vector<PseudoJet>();  
    }
    
    // Define the clustering algorithm (anti-kT by default; kT available
    // via CLI for algorithmic-comparison studies).
    JetDefinition jet_def(jet_algo_in, R);
    
    // Run clustering
    ClusterSequence cs(input_particles, jet_def);
    
    // Get jets and apply cuts
    vector<PseudoJet> jets = cs.inclusive_jets();
    vector<PseudoJet> selected_jets;
    
    for (const auto& jet : jets) {
        double et = jet.Et();
        double jet_eta = jet.eta();
        
        // Apply basic jet-level cuts
        if (et > etMin && jet_eta > etaMin && jet_eta < etaMax) {
            selected_jets.push_back(jet);
        }
    }
    
    // Sort jets by ET (descending)
    sort(selected_jets.begin(), selected_jets.end(), 
         [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
    
    // -------------------------------------------------------------------------
    // Per-jet substructure observables: n_sd (softdrop multiplicity) and
    // n_subjets (kT exclusive with ycut). Both must be computed here, inside
    // reconstructJets(), because the ClusterSequence `cs` is still alive —
    // which means jet.constituents() is valid. Once this function returns,
    // cs dies and any structural access on the returned jets would crash.
    // We keep only the two scalars per jet, independent of cs.
    // -------------------------------------------------------------------------
    jet_nsd_out.reserve(selected_jets.size());
    jet_nsubjets_out.reserve(selected_jets.size());
    jet_psi03_out.reserve(selected_jets.size());
    jet_nsd_beta1_out.reserve(selected_jets.size());
    jet_nsd_loose_out.reserve(selected_jets.size());

    // Integrated jet shape psi(r) at r = 0.3 (matches the convention used
    // by src/jetshapes/integrated/jetrecoint.cc and the paper's Fig 9).
    const double psi_r = 0.3;

    for (const auto& jet : selected_jets) {
        vector<PseudoJet> constituents = jet.constituents();

        // --- psi(r=0.3): fraction of jet ET within Delta-R < 0.3 ---------
        double et_in_r = 0.0;
        const double et_jet = jet.Et();
        if (et_jet > 0.0) {
            for (const auto& c : constituents) {
                if (c.delta_R(jet) <= psi_r) {
                    et_in_r += c.Et();
                }
            }
            jet_psi03_out.push_back(static_cast<float>(et_in_r / et_jet));
        } else {
            jet_psi03_out.push_back(0.0f);
        }

        // --- n_sd via Cambridge/Aachen reclustering + tree walk ------------
        // Compute three configurations in the same walk so a sensitivity
        // study is free:
        //   (z_cut, beta) = (0.1,  0.0)  -- paper baseline (mMDT-like)
        //                 = (0.1,  1.0)  -- original angle-dependent SD
        //                 = (0.05, 0.0)  -- loose z_cut variant
        int n_sd         = 0;   // (0.1, 0.0) — baseline
        int n_sd_beta1   = 0;   // (0.1, 1.0)
        int n_sd_loose   = 0;   // (0.05, 0.0)
        if (constituents.size() >= 2) {
            JetDefinition ca_def(cambridge_algorithm, sd_R0);
            ClusterSequence cs_ca(constituents, ca_def);
            vector<PseudoJet> ca_jets = sorted_by_pt(cs_ca.inclusive_jets());

            if (!ca_jets.empty()) {
                PseudoJet current = ca_jets[0];
                PseudoJet p1, p2;
                while (current.has_parents(p1, p2)) {
                    if (p1.pt() < p2.pt()) std::swap(p1, p2);   // p1 = harder branch
                    double dR = p1.delta_R(p2);
                    double pt_sum = p1.pt() + p2.pt();
                    if (pt_sum > 0.0) {
                        double z = p2.pt() / pt_sum;
                        // Baseline (beta = 0, so the angular factor is 1):
                        if (z > sd_zcut)                              ++n_sd;
                        // Angular-dependent (beta = 1):
                        if (z > sd_zcut * (dR / sd_R0))               ++n_sd_beta1;
                        // Loose z_cut (beta = 0):
                        if (z > 0.05)                                 ++n_sd_loose;
                    }
                    current = p1;   // follow the harder branch inward
                }
            }
        }
        jet_nsd_out.push_back(n_sd);
        jet_nsd_beta1_out.push_back(n_sd_beta1);
        jet_nsd_loose_out.push_back(n_sd_loose);

        // --- n_subjets via kT exclusive with ycut --------------------------
        // Matches the recipe in subprocsubjets/nsubjets.cc so results are
        // directly comparable. Jets with <2 constituents get n_subjets = 0.
        int n_subjets = 0;
        if (constituents.size() >= 2) {
            JetDefinition kt_def(kt_algorithm, R);
            ClusterSequence cs_kt(constituents, kt_def);
            n_subjets = cs_kt.exclusive_jets_ycut(kt_ycut).size();
        }
        jet_nsubjets_out.push_back(n_subjets);

        // cs_ca, cs_kt go out of scope here — the integers are already saved
    }
    // -------------------------------------------------------------------------
    
    return selected_jets;
}

// =============================================================================
// DIJET SELECTION FUNCTION
// =============================================================================
bool passesDijetCuts(const vector<PseudoJet>& jets) {
    // Need at least 2 jets
    if (jets.size() < 2) return false;
    
    PseudoJet leading_jet = jets[0];
    PseudoJet subleading_jet = jets[1];
    
    // Leading jet cuts
    if (leading_jet.Et() < leading_jet_et_min) return false;
    if (abs(leading_jet.eta()) > leading_jet_eta_max) return false;
    
    // Subleading jet cuts
    if (subleading_jet.Et() < subleading_jet_et_min) return false;
    if (abs(subleading_jet.eta()) > subleading_jet_eta_max) return false;
    
    // Calculate dijet properties
    TLorentzVector p4_leading, p4_subleading;
    p4_leading.SetPxPyPzE(leading_jet.px(), leading_jet.py(), leading_jet.pz(), leading_jet.E());
    p4_subleading.SetPxPyPzE(subleading_jet.px(), subleading_jet.py(), subleading_jet.pz(), subleading_jet.E());
    
    // // Dijet invariant mass cut
    // double dijet_mass = (p4_leading + p4_subleading).M();
    // if (dijet_mass < dijet_mass_min) return false;
    
    // // Azimuthal separation cut (back-to-back topology)
    // double delta_phi = abs(leading_jet.phi() - subleading_jet.phi());
    // if (delta_phi > M_PI) delta_phi = 2*M_PI - delta_phi;  // Wrap to [0, π]
    // if (delta_phi < delta_phi_min) return false;
    
    // // Rapidity separation cut
    // double delta_eta = abs(leading_jet.eta() - subleading_jet.eta());
    // if (delta_eta > delta_eta_max) return false;
    
    return true;
}

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================
int main(int argc, char** argv) {

    // Input file: CLI arg 1, else default. Output/log names derived from it.
    //   argv[1] : input ROOT file path
    //   argv[2] : jet algorithm   ("antikt" [default] or "kt")
    input_filename = (argc > 1) ? argv[1] : default_input_filename;

    JetAlgorithm jet_algo = antikt_algorithm;
    string       jet_algo_str = "anti-kT";
    if (argc > 2) {
        string arg2 = argv[2];
        if (arg2 == "kt" || arg2 == "kT") {
            jet_algo = kt_algorithm;
            jet_algo_str = "kT";
        } else if (arg2 == "antikt" || arg2 == "anti-kt" || arg2 == "anti-kT") {
            jet_algo = antikt_algorithm;
            jet_algo_str = "anti-kT";
        } else {
            cerr << "Unknown jet algorithm '" << arg2
                 << "'. Options: antikt (default), kt." << endl;
            return 1;
        }
    }

    output_filename = generateOutputFilename(input_filename, R, etMin);
    log_filename = generateLogFilename(output_filename);

    // Create dual output for console and file logging
    DualOutput dout(log_filename);

    // Open input file
    TFile* input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file || input_file->IsZombie()) {
        cerr << "ERROR: Cannot open input file: " << input_filename << endl;
        return 1;
    }
    
    dout << "=============================================================================\n";
    dout << "JET RECONSTRUCTION ANALYSIS FOR HERA PHOTOPRODUCTION EVENTS\n";
    dout << "=============================================================================\n";
    dout << "Analysis Parameters:\n";
    dout << "  Analysis Mode:    " << (DIJET_ONLY ? "DIJET ONLY" : "ALL JETS") << "\n";
    dout << "  Jet Algorithm:    " << jet_algo_str << "\n";
    dout << "  Jet Radius:       " << R << "\n";
    dout << "  Min Jet ET:       " << etMin << " GeV\n";
    dout << "  Jet Eta Range:    [" << etaMin << ", " << etaMax << "]\n";
    dout << "  Min Particle pT:  " << particle_pt_min << " GeV\n";
    dout << "  Softdrop (zcut, beta, R0): (" << sd_zcut << ", " << sd_beta << ", " << sd_R0 << ")\n";
    dout << "  n_subjets (kT ycut):       " << kt_ycut << "\n";
    
    if (DIJET_ONLY) {
        dout << "\nDijet Selection Cuts:\n";
        dout << "  Leading Jet ET:   > " << leading_jet_et_min << " GeV\n";
        dout << "  Subleading Jet ET:> " << subleading_jet_et_min << " GeV\n";
        dout << "  Leading Jet |η|:  < " << leading_jet_eta_max << "\n";
        dout << "  Subleading |η|:   < " << subleading_jet_eta_max << "\n";
        // dout << "  Dijet Mass:       > " << dijet_mass_min << " GeV\n";
        // dout << "  Δφ(j1,j2):        > " << delta_phi_min << " rad\n";
        // dout << "  Δη(j1,j2):        < " << delta_eta_max << "\n";
    }
    dout << "=============================================================================\n\n";
    
    dout << "Successfully opened input file: " << input_filename << "\n";
    dout << "Log file: " << log_filename << "\n\n";
    
    // Create output file
    TFile* output_file = new TFile(output_filename.c_str(), "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        cerr << "ERROR: Cannot create output file: " << output_filename << endl;
        return 1;
    }
    
    // Define event categories to process
    vector<string> categories = {"QQ_Events", "GG_Events", "GQ_Events", 
                                "Combined_Events", "Resolved_Events", "Direct_Events"};
    
    // Overall statistics
    int total_events_processed = 0;
    int total_events_selected = 0;
    
    // Process each category
    for (const string& category : categories) {
        
        dout << "Processing " << category << "...\n";
        
        // Get tree from input file
        TDirectory* input_dir = (TDirectory*)input_file->Get(category.c_str());
        if (!input_dir) {
            dout << "  Warning: Directory " << category << " not found, skipping...\n";
            continue;
        }
        
        TTree* input_tree = (TTree*)input_dir->Get(category.c_str());
        if (!input_tree) {
            dout << "  Warning: Tree " << category << " not found, skipping...\n";
            continue;
        }
        
        Long64_t n_entries = input_tree->GetEntries();
        dout << "  Found " << n_entries << " events\n";
        
        if (n_entries == 0) {
            dout << "  No events found, skipping...\n";
            continue;
        }
        
        // Set up input branches - check if branches exist first
        vector<float> *px = nullptr, *py = nullptr, *pz = nullptr, *energy = nullptr, *eta = nullptr;

        if (!input_tree->GetBranch("px") || !input_tree->GetBranch("py") ||
            !input_tree->GetBranch("pz") || !input_tree->GetBranch("energy") ||
            !input_tree->GetBranch("eta")) {
            dout << "  Error: Required branches not found in tree, skipping...\n";
            continue;
        }

        input_tree->SetBranchAddress("px", &px);
        input_tree->SetBranchAddress("py", &py);
        input_tree->SetBranchAddress("pz", &pz);
        input_tree->SetBranchAddress("energy", &energy);
        input_tree->SetBranchAddress("eta", &eta);

        // Optional input branches: outgoing hard-process partons
        // (status==23) for jet-to-parton DeltaR matching. Present only
        // in validation samples produced after commit e465c51.
        vector<int>   *parton_pdgId = nullptr;
        vector<float> *parton_eta_in = nullptr, *parton_phi_in = nullptr;
        const bool have_partons =
            input_tree->GetBranch("parton_pdgId") &&
            input_tree->GetBranch("parton_eta") &&
            input_tree->GetBranch("parton_phi");
        if (have_partons) {
            input_tree->SetBranchAddress("parton_pdgId", &parton_pdgId);
            input_tree->SetBranchAddress("parton_eta",   &parton_eta_in);
            input_tree->SetBranchAddress("parton_phi",   &parton_phi_in);
        }
        
        // Create output directory and tree
        output_file->cd();
        TDirectory* output_dir = output_file->mkdir(category.c_str());
        output_dir->cd();
        
        TTree* output_tree = new TTree(("jets_" + category).c_str(), 
                                      ("Reconstructed Jets for " + category).c_str());
        
        // Output variables - basic jet info
        Int_t eventID, n_jets;
        vector<float> jet_et, jet_px, jet_py, jet_pz, jet_eta, jet_phi, jet_mass;
        vector<int>   jet_nsd;          // softdrop multiplicity per saved jet  (z_cut=0.1, beta=0)
        vector<int>   jet_nsd_beta1;    // softdrop multiplicity (z_cut=0.1, beta=1)
        vector<int>   jet_nsd_loose;    // softdrop multiplicity (z_cut=0.05, beta=0)
        vector<int>   jet_nsubjets;     // kT exclusive subjet count per saved jet
        vector<float> jet_psi03;        // integrated jet shape psi(r=0.3)
        vector<float> jet_parton_dR;    // DeltaR to nearest hard parton (−1 if no partons)
        vector<int>   jet_parton_pdgId; // PDG id of nearest hard parton (0 if no partons)
        
        // Dijet-specific variables
        Bool_t passes_dijet_cuts;
        Float_t leading_jet_et, subleading_jet_et;
        Float_t leading_jet_eta, subleading_jet_eta;
        Float_t leading_jet_phi, subleading_jet_phi;
        Float_t dijet_mass, dijet_delta_phi, dijet_delta_eta;
        Float_t dijet_pt, dijet_y;
        
        // Set up output branches - basic jet info
        output_tree->Branch("eventID", &eventID);
        output_tree->Branch("n_jets", &n_jets);
        output_tree->Branch("jet_et", &jet_et);
        output_tree->Branch("jet_px", &jet_px);
        output_tree->Branch("jet_py", &jet_py);
        output_tree->Branch("jet_pz", &jet_pz);
        output_tree->Branch("jet_eta", &jet_eta);
        output_tree->Branch("jet_phi", &jet_phi);
        output_tree->Branch("jet_mass", &jet_mass);
        output_tree->Branch("jet_nsd", &jet_nsd);
        output_tree->Branch("jet_nsd_beta1", &jet_nsd_beta1);
        output_tree->Branch("jet_nsd_loose", &jet_nsd_loose);
        output_tree->Branch("jet_nsubjets", &jet_nsubjets);
        output_tree->Branch("jet_psi03", &jet_psi03);
        output_tree->Branch("jet_parton_dR", &jet_parton_dR);
        output_tree->Branch("jet_parton_pdgId", &jet_parton_pdgId);
        
        // Dijet-specific branches
        output_tree->Branch("passes_dijet_cuts", &passes_dijet_cuts);
        output_tree->Branch("leading_jet_et", &leading_jet_et);
        output_tree->Branch("subleading_jet_et", &subleading_jet_et);
        output_tree->Branch("leading_jet_eta", &leading_jet_eta);
        output_tree->Branch("subleading_jet_eta", &subleading_jet_eta);
        output_tree->Branch("leading_jet_phi", &leading_jet_phi);
        output_tree->Branch("subleading_jet_phi", &subleading_jet_phi);
        output_tree->Branch("dijet_mass", &dijet_mass);
        output_tree->Branch("dijet_delta_phi", &dijet_delta_phi);
        output_tree->Branch("dijet_delta_eta", &dijet_delta_eta);
        output_tree->Branch("dijet_pt", &dijet_pt);
        output_tree->Branch("dijet_y", &dijet_y);
        
        // Event counters
        int events_processed = 0;
        int events_with_jets = 0;
        int events_selected = 0;
        int total_jets = 0;
        
        // Process events
        for (Long64_t i = 0; i < n_entries; ++i) {
            
            if (i % report_every == 0) {
                dout << "    Processing event " << i << "/" << n_entries 
                     << " (" << (100*i/n_entries) << "%)\r" << flush;
            }
            
            input_tree->GetEntry(i);
            
            // Check if vectors are valid and not empty
            if (!px || !py || !pz || !energy || !eta) {
                continue;
            }
            
            // Check if all vectors have the same size
            if (px->size() != py->size() || px->size() != pz->size() || 
                px->size() != energy->size() || px->size() != eta->size()) {
                continue;
            }
            
            events_processed++;
            
            // Clear output vectors
            jet_et.clear();
            jet_px.clear();
            jet_py.clear();
            jet_pz.clear();
            jet_eta.clear();
            jet_phi.clear();
            jet_mass.clear();
            jet_nsd.clear();
            jet_nsd_beta1.clear();
            jet_nsd_loose.clear();
            jet_nsubjets.clear();
            jet_psi03.clear();
            jet_parton_dR.clear();
            jet_parton_pdgId.clear();

            // Reconstruct jets (and compute per-jet n_sd, n_subjets, psi(r) inline).
            // n_sd_beta1 and n_sd_loose are cheap add-ons computed in the same
            // C/A walk; used by the softdrop-parameter-sensitivity study.
            vector<int>   jets_nsd_tmp, jets_nsubjets_tmp;
            vector<int>   jets_nsd_beta1_tmp, jets_nsd_loose_tmp;
            vector<float> jets_psi03_tmp;
            vector<PseudoJet> jets = reconstructJets(*px, *py, *pz, *energy, *eta,
                                                    jets_nsd_tmp, jets_nsubjets_tmp,
                                                    jets_psi03_tmp,
                                                    jets_nsd_beta1_tmp,
                                                    jets_nsd_loose_tmp,
                                                    jet_algo);
            
            // Basic event info
            eventID = i;
            n_jets = jets.size();
            
            // Initialize dijet variables
            passes_dijet_cuts = false;
            leading_jet_et = subleading_jet_et = -999.0;
            leading_jet_eta = subleading_jet_eta = -999.0;
            leading_jet_phi = subleading_jet_phi = -999.0;
            dijet_mass = dijet_delta_phi = dijet_delta_eta = -999.0;
            dijet_pt = dijet_y = -999.0;
            
            // Fill jet info and check dijet criteria
            bool fill_event = false;
            
            if (n_jets >= 2) {
                events_with_jets++;
                
                // Check dijet cuts
                passes_dijet_cuts = passesDijetCuts(jets);
                
                // Fill leading and subleading jet info
                leading_jet_et = jets[0].Et();
                subleading_jet_et = jets[1].Et();
                leading_jet_eta = jets[0].eta();
                subleading_jet_eta = jets[1].eta();
                leading_jet_phi = jets[0].phi();
                subleading_jet_phi = jets[1].phi();
                
                // Calculate dijet system properties
                TLorentzVector p4_leading, p4_subleading;
                p4_leading.SetPxPyPzE(jets[0].px(), jets[0].py(), jets[0].pz(), jets[0].E());
                p4_subleading.SetPxPyPzE(jets[1].px(), jets[1].py(), jets[1].pz(), jets[1].E());
                TLorentzVector dijet_system = p4_leading + p4_subleading;
                
                dijet_mass = dijet_system.M();
                dijet_pt = dijet_system.Pt();
                dijet_y = dijet_system.Rapidity();
                
                dijet_delta_phi = abs(jets[0].phi() - jets[1].phi());
                if (dijet_delta_phi > M_PI) dijet_delta_phi = 2*M_PI - dijet_delta_phi;
                dijet_delta_eta = abs(jets[0].eta() - jets[1].eta());
                
                // Decide whether to fill this event
                if (DIJET_ONLY) {
                    fill_event = passes_dijet_cuts;
                } else {
                    fill_event = true;  // Fill all events with ≥2 jets for "all jets" mode
                }
            }
            
            // Helper: signed delta-phi wrapped into (-pi, pi], used for
            // jet-parton matching.
            auto wrapped_dphi = [](double a, double b) {
                double d = a - b;
                while (d > M_PI) d -= 2.0 * M_PI;
                while (d <= -M_PI) d += 2.0 * M_PI;
                return d;
            };

            // Helper: nearest hard-parton match for a given jet. Returns
            // (dR, pdgId). Falls back to (-1, 0) if this input tree has
            // no parton branches or they are empty/null for this event.
            auto match_parton = [&](const PseudoJet& jet) {
                if (!have_partons || parton_eta_in == nullptr ||
                    parton_phi_in == nullptr || parton_pdgId == nullptr ||
                    parton_eta_in->empty()) {
                    return std::make_pair(-1.0f, 0);
                }
                double best_dR = 1e9;
                int    best_id = 0;
                double je = jet.eta();
                double jp = jet.phi();
                for (size_t k = 0; k < parton_eta_in->size(); ++k) {
                    double deta = je - (*parton_eta_in)[k];
                    double dphi = wrapped_dphi(jp, (*parton_phi_in)[k]);
                    double dR = sqrt(deta * deta + dphi * dphi);
                    if (dR < best_dR) {
                        best_dR = dR;
                        best_id = (*parton_pdgId)[k];
                    }
                }
                return std::make_pair((float)best_dR, best_id);
            };

            // Fill jet vectors based on mode
            if (fill_event) {
                events_selected++;

                if (DIJET_ONLY) {
                    // Only save leading and subleading jets
                    total_jets += 2;
                    for (int j = 0; j < 2; ++j) {
                        jet_et.push_back(jets[j].Et());
                        jet_px.push_back(jets[j].px());
                        jet_py.push_back(jets[j].py());
                        jet_pz.push_back(jets[j].pz());
                        jet_eta.push_back(jets[j].eta());
                        jet_phi.push_back(jets[j].phi());
                        jet_mass.push_back(jets[j].m());
                        jet_nsd.push_back(jets_nsd_tmp[j]);
                        jet_nsd_beta1.push_back(jets_nsd_beta1_tmp[j]);
                        jet_nsd_loose.push_back(jets_nsd_loose_tmp[j]);
                        jet_nsubjets.push_back(jets_nsubjets_tmp[j]);
                        jet_psi03.push_back(jets_psi03_tmp[j]);
                        auto pm = match_parton(jets[j]);
                        jet_parton_dR.push_back(pm.first);
                        jet_parton_pdgId.push_back(pm.second);
                    }
                } else {
                    // Save all jets
                    total_jets += n_jets;
                    for (size_t j = 0; j < jets.size(); ++j) {
                        const auto& jet = jets[j];
                        jet_et.push_back(jet.Et());
                        jet_px.push_back(jet.px());
                        jet_py.push_back(jet.py());
                        jet_pz.push_back(jet.pz());
                        jet_eta.push_back(jet.eta());
                        jet_phi.push_back(jet.phi());
                        jet_mass.push_back(jet.m());
                        jet_nsd.push_back(jets_nsd_tmp[j]);
                        jet_nsd_beta1.push_back(jets_nsd_beta1_tmp[j]);
                        jet_nsd_loose.push_back(jets_nsd_loose_tmp[j]);
                        jet_nsubjets.push_back(jets_nsubjets_tmp[j]);
                        jet_psi03.push_back(jets_psi03_tmp[j]);
                        auto pm = match_parton(jet);
                        jet_parton_dR.push_back(pm.first);
                        jet_parton_pdgId.push_back(pm.second);
                    }
                }

                output_tree->Fill();
            }
        }
        
        dout << "\n  Events processed: " << events_processed << "\n";
        dout << "  Events with ≥2 jets: " << events_with_jets 
             << " (" << (100.0*events_with_jets/events_processed) << "%)\n";
        
        if (DIJET_ONLY) {
            dout << "  Events passing dijet cuts: " << events_selected 
                 << " (" << (100.0*events_selected/events_processed) << "%)\n";
        } else {
            dout << "  Events selected: " << events_selected 
                 << " (" << (100.0*events_selected/events_processed) << "%)\n";
        }
        
        dout << "  Total jets saved: " << total_jets << "\n";
        dout << "  Average jets/selected event: " 
             << (events_selected > 0 ? (double)total_jets/events_selected : 0) << "\n\n";
        
        // Update global counters
        total_events_processed += events_processed;
        total_events_selected += events_selected;
        
        // Write tree
        output_tree->Write();
        
        // Create basic histograms
        TH1F* h_n_jets = new TH1F(("h_n_jets_" + category).c_str(), 
                                 ("Number of Jets per Event - " + category + ";N_{jets};Events").c_str(), 
                                 10, 0, 10);
        TH1F* h_jet_et = new TH1F(("h_jet_et_" + category).c_str(), 
                                 ("Jet E_{T} Distribution - " + category + ";E_{T} [GeV];Jets").c_str(), 
                                 100, 0, 200);
        TH1F* h_jet_eta = new TH1F(("h_jet_eta_" + category).c_str(), 
                                  ("Jet #eta Distribution - " + category + ";#eta;Jets").c_str(), 
                                  80, -4, 4);
        
        // Dijet-specific histograms
        TH1F* h_dijet_mass = new TH1F(("h_dijet_mass_" + category).c_str(),
                                     ("Dijet Mass - " + category + ";M_{jj} [GeV];Events").c_str(),
                                     100, 0, 300);
        TH1F* h_dijet_delta_phi = new TH1F(("h_dijet_delta_phi_" + category).c_str(),
                                          ("Dijet #Delta#phi - " + category + ";#Delta#phi [rad];Events").c_str(),
                                          50, 0, M_PI);
        
        // Fill histograms by reading back the tree
        vector<float> *hist_jet_et = nullptr, *hist_jet_eta = nullptr;
        Int_t hist_n_jets;
        Bool_t hist_passes_cuts;
        Float_t hist_dijet_mass, hist_dijet_dphi;
        
        output_tree->SetBranchAddress("n_jets", &hist_n_jets);
        output_tree->SetBranchAddress("jet_et", &hist_jet_et);
        output_tree->SetBranchAddress("jet_eta", &hist_jet_eta);
        output_tree->SetBranchAddress("passes_dijet_cuts", &hist_passes_cuts);
        output_tree->SetBranchAddress("dijet_mass", &hist_dijet_mass);
        output_tree->SetBranchAddress("dijet_delta_phi", &hist_dijet_dphi);
        
        for (Long64_t i = 0; i < output_tree->GetEntries(); ++i) {
            output_tree->GetEntry(i);
            h_n_jets->Fill(hist_n_jets);
            
            if (hist_jet_et && hist_jet_eta) {
                for (size_t j = 0; j < hist_jet_et->size(); ++j) {
                    h_jet_et->Fill((*hist_jet_et)[j]);
                    h_jet_eta->Fill((*hist_jet_eta)[j]);
                }
            }
            
            if (hist_passes_cuts && hist_dijet_mass > 0) {
                h_dijet_mass->Fill(hist_dijet_mass);
                h_dijet_delta_phi->Fill(hist_dijet_dphi);
            }
        }
        
        h_n_jets->Write();
        h_jet_et->Write();
        h_jet_eta->Write();
        h_dijet_mass->Write();
        h_dijet_delta_phi->Write();
    }
    
    // Write and close files
    output_file->Write();
    output_file->Close();
    input_file->Close();
    
    dout << "=============================================================================\n";
    dout << "ANALYSIS COMPLETED SUCCESSFULLY!\n";
    dout << "=============================================================================\n";
    dout << "Output file: " << output_filename << "\n";
    dout << "Log file: " << log_filename << "\n";
    dout << "\nOverall Statistics:\n";
    dout << "  Total events processed: " << total_events_processed << "\n";
    dout << "  Total events selected: " << total_events_selected 
         << " (" << (100.0*total_events_selected/total_events_processed) << "%)\n";
    
    dout << "\nFor each category, the following are saved:\n";
    dout << "  - Jet trees with basic properties (ET, px, py, pz, eta, phi, mass, nsd, nsubjets)\n";
    dout << "  - Dijet-specific variables (masses, angles, cuts flags)\n";
    dout << "  - Basic distribution histograms\n";
    dout << "\nJets are reconstructed using anti-kT algorithm with:\n";
    dout << "  - R = " << R << "\n";
    dout << "  - ET > " << etMin << " GeV\n";
    dout << "  - " << etaMin << " < eta < " << etaMax << "\n";
    dout << "\nSoftdrop multiplicity n_sd (modified mass drop):\n";
    dout << "  - z_cut = " << sd_zcut << "\n";
    dout << "  - beta  = " << sd_beta << "\n";
    dout << "  - R0    = " << sd_R0 << "\n";
    dout << "\nn_subjets (kT exclusive):\n";
    dout << "  - y_cut = " << kt_ycut << "\n";
    
    if (DIJET_ONLY) {
        dout << "\nDijet selection applied with:\n";
        dout << "  - Leading jet ET > " << leading_jet_et_min << " GeV\n";
        dout << "  - Subleading jet ET > " << subleading_jet_et_min << " GeV\n";
        // dout << "  - Dijet mass > " << dijet_mass_min << " GeV\n";
    }
    dout << "=============================================================================\n";
    
    return 0;
}