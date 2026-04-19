// =============================================================================
// Event Shape Classification Analysis for EIC Photoproduction Events
// =============================================================================
// Description: Classifies dijet events as GG-like, QQ-like, or GQ-like based 
//              on jet shape combinations and analyzes relative fractions
//
// Input: ROOT file with Combined_Events from photoproduction
// Output: Console log with event shape fractions by eta range
// =============================================================================

#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <map>
#include <memory>

using namespace std;
using namespace fastjet;

// =============================================================================
// ANALYSIS CONFIGURATION (MATCHING jetreco.cc)
// =============================================================================

// Analysis mode selection (MUST be true for event shape analysis)
const bool DIJET_ONLY = true;    // Event shapes require exactly 2 jets

// Jet clustering parameters
const double R = 1.0;             // Jet radius
const double etMin = 0.0;         // Minimum ET for any jet
const double etaMin = -4.0;       // Minimum eta for any jet
const double etaMax = 4.0;        // Maximum eta for any jet

// Dijet-specific cuts
const double leading_jet_et_min = 10.0;      // Minimum ET for leading jet (GeV)
const double subleading_jet_et_min = 7.0;    // Minimum ET for subleading jet (GeV)
// const double leading_jet_eta_max = 4.5;      // Maximum |eta| for leading jet
// const double subleading_jet_eta_max = 4.5;   // Maximum |eta| for subleading jet

// Particle selection cuts
// const double particle_pt_min = 0.1;     // Minimum particle pT (GeV)
// const double particle_eta_max = 5.0;    // Maximum particle |eta|

// Thick/Thin jet classification cuts (from HERA studies)
const double thick_jet_cut = 0.6;       // Ψ(r=0.3) < 0.6 for thick jets (gluon-like)
const double thin_jet_cut = 0.8;        // Ψ(r=0.3) > 0.8 for thin jets (quark-like)

// Progress reporting
const int report_every = 100000;

// Input file
// const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt5GeV/eic141_pt5/eic141.root";
const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt7GeV/hera300_pt7/hera300_pt7.root";

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
// OUTPUT FILENAME GENERATION
// =============================================================================
string generateLogFilename(const string& input_filename) {
    string baseFileName = input_filename.substr(input_filename.find_last_of("/") + 1);
    baseFileName = baseFileName.substr(0, baseFileName.find_last_of("."));  // Remove .root extension
    return "event_shape_analysis_" + baseFileName + ".log";
}

const string log_filename = generateLogFilename(input_filename);

// =============================================================================
// INTEGRATED JET SHAPE CALCULATION (FROM PROJECT DATABASE)
// =============================================================================
double calculateIntegratedJetShape(const PseudoJet& jet, double r, const ClusterSequence& cs) {
    
    vector<PseudoJet> constituents = jet.constituents();
    double energy_within_r = 0.0;
    
    // Loop over constituents and select ones within radius r
    for (const auto& constituent : constituents) {
        double dphi = jet.phi() - constituent.phi();
        double drap = jet.eta() - constituent.eta();
        double const_r = sqrt(dphi*dphi + drap*drap);
        
        // Add ET of constituents within radius r
        if (const_r <= r) {
            energy_within_r += constituent.Et();
        }
    }
    
    // Calculate integrated jet shape Ψ(r)
    if (jet.Et() > 0.0) {
        return energy_within_r / jet.Et();
    }
    
    return 0.0;
}

// =============================================================================
// JET RECONSTRUCTION FUNCTION (FROM jetreco.cc)
// =============================================================================
struct JetReconstructionResult {
    vector<PseudoJet> jets;
    shared_ptr<ClusterSequence> cs;
};

JetReconstructionResult reconstructJets(const vector<float>& px, const vector<float>& py, 
                                       const vector<float>& pz, const vector<float>& energy,
                                       const vector<float>& eta) {
    
    JetReconstructionResult result;
    
    // Check input vector sizes
    if (px.empty() || px.size() != py.size() || px.size() != pz.size() || 
        px.size() != energy.size() || px.size() != eta.size()) {
        return result;
    }
    
    // Create input particles for FastJet
    vector<PseudoJet> input_particles;
    
    for (size_t i = 0; i < px.size(); ++i) {
        // Apply particle-level cuts
        double pt = sqrt(px[i]*px[i] + py[i]*py[i]);
        // if (pt < particle_pt_min) continue;
        // if (abs(eta[i]) > particle_eta_max) continue;
        if (energy[i] <= 0) continue;
        
        // Create PseudoJet
        PseudoJet particle(px[i], py[i], pz[i], energy[i]);
        particle.set_user_index(i);
        input_particles.push_back(particle);
    }
    
    if (input_particles.size() < 2) {
        return result;
    }
    
    // Define anti-kT algorithm
    JetDefinition jet_def(antikt_algorithm, R);
    
    // Run clustering - keep ClusterSequence alive
    result.cs = make_shared<ClusterSequence>(input_particles, jet_def);
    
    // Get jets and apply cuts
    vector<PseudoJet> jets = result.cs->inclusive_jets();
    
    for (const auto& jet : jets) {
        double et = jet.Et();
        double jet_eta = jet.eta();
        
        // Apply basic jet-level cuts
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
// DIJET SELECTION FUNCTION (FROM jetreco.cc)
// =============================================================================
bool passesDijetCuts(const vector<PseudoJet>& jets) {
    // Need exactly 2 jets for event shape analysis
    if (jets.size() < 2) return false;
    
    PseudoJet leading_jet = jets[0];
    PseudoJet subleading_jet = jets[1];
    
    // Leading jet cuts
    if (leading_jet.Et() < leading_jet_et_min) return false;
    // if (abs(leading_jet.eta()) > leading_jet_eta_max) return false;
    
    // Subleading jet cuts
    if (subleading_jet.Et() < subleading_jet_et_min) return false;
    // if (abs(subleading_jet.eta()) > subleading_jet_eta_max) return false;
    
    return true;
}

// =============================================================================
// JET CLASSIFICATION FUNCTION
// =============================================================================
string classifyJet(const PseudoJet& jet, const ClusterSequence& cs) {
    
    // Calculate integrated jet shape at r = 0.3 (standard cut from HERA studies)
    double psi_03 = calculateIntegratedJetShape(jet, 0.3, cs);
    
    // Apply thick/thin classification cuts
    if (psi_03 < thick_jet_cut) {
        return "thick";  // Gluon-like (broad energy distribution)
    } else if (psi_03 > thin_jet_cut) {
        return "thin";   // Quark-like (narrow energy distribution)
    } else {
        return "intermediate";  // Intermediate jets
    }
}

// =============================================================================
// EVENT SHAPE CLASSIFICATION FUNCTION
// =============================================================================
string classifyEventShape(const vector<PseudoJet>& jets, const ClusterSequence& cs) {
    
    // Must have exactly 2 jets for event shape analysis
    if (jets.size() != 2) return "invalid";
    
    string jet1_type = classifyJet(jets[0], cs);  // Leading jet
    string jet2_type = classifyJet(jets[1], cs);  // Subleading jet
    
    // Classify event based on jet shape combination
    if (jet1_type == "thick" && jet2_type == "thick") {
        return "GG-like";  // 2 Thick Jets (spherical)
    } else if (jet1_type == "thin" && jet2_type == "thin") {
        return "QQ-like";  // 2 Thin Jets (pencil-like)
    } else if ((jet1_type == "thick" && jet2_type == "thin") || 
               (jet1_type == "thin" && jet2_type == "thick")) {
        return "GQ-like";  // 1 Thick + 1 Thin Jet (mixed)
    } else {
        return "intermediate";  // Contains intermediate jets
    }
}

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================
int main() {

    // Create dual output for console and file logging
    DualOutput dout(log_filename);

    dout << "=============================================================================\n";
    dout << "EVENT SHAPE CLASSIFICATION ANALYSIS\n";
    dout << "=============================================================================\n";
    dout << "Event Shape Classification:\n";
    dout << "  GG-like Events: 2 Thick Jets (spherical, gluon-like)\n";
    dout << "  QQ-like Events: 2 Thin Jets (pencil-like, quark-like)\n";
    dout << "  GQ-like Events: 1 Thick + 1 Thin Jet (mixed)\n";
    dout << "\nJet Classification Criteria:\n";
    dout << "  Thick Jets: Ψ(r=0.3) < " << thick_jet_cut << "\n";
    dout << "  Thin Jets:  Ψ(r=0.3) > " << thin_jet_cut << "\n";
    dout << "\nAnalysis Parameters:\n";
    dout << "  Analysis Mode:    DIJET ONLY (required for event shapes)\n";
    dout << "  Jet Algorithm:    Anti-kT (R = " << R << ")\n";
    dout << "  Leading Jet ET:   > " << leading_jet_et_min << " GeV\n";
    dout << "  Subleading ET:    > " << subleading_jet_et_min << " GeV\n";
    dout << "=============================================================================\n\n";

    // Open input file
    TFile* input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file || input_file->IsZombie()) {
        dout << "ERROR: Cannot open input file: " << input_filename << endl;
        return 1;
    }
    
    dout << "Processing file: " << input_filename << "\n";
    dout << "Log file: " << log_filename << "\n\n";
    
    // Process Combined_Events (contains all subprocess types)
    string category = "Combined_Events";
    
    dout << "Processing " << category << "...\n";
    
    // Get tree from input file
    TDirectory* input_dir = (TDirectory*)input_file->Get(category.c_str());
    if (!input_dir) {
        dout << "ERROR: Directory " << category << " not found!" << endl;
        return 1;
    }
    
    TTree* input_tree = (TTree*)input_dir->Get(category.c_str());
    if (!input_tree) {
        dout << "ERROR: Tree " << category << " not found!" << endl;
        return 1;
    }
    
    Long64_t n_entries = input_tree->GetEntries();
    dout << "Found " << n_entries << " events\n";
    
    // Set up input branches
    vector<float> *px = nullptr, *py = nullptr, *pz = nullptr, *energy = nullptr, *eta = nullptr;
    
    if (!input_tree->GetBranch("px") || !input_tree->GetBranch("py") || 
        !input_tree->GetBranch("pz") || !input_tree->GetBranch("energy") || 
        !input_tree->GetBranch("eta")) {
        dout << "ERROR: Required branches not found in tree!" << endl;
        return 1;
    }
    
    input_tree->SetBranchAddress("px", &px);
    input_tree->SetBranchAddress("py", &py);
    input_tree->SetBranchAddress("pz", &pz);
    input_tree->SetBranchAddress("energy", &energy);
    input_tree->SetBranchAddress("eta", &eta);
    
    // Define eta ranges for counting: -1 to 0, 0 to 1, 1 to 2, 2 to 3, 3 to 4
    vector<pair<double, double>> eta_ranges = {
        {-1.0, 0.0}, {0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}
    };
    
    vector<string> eta_range_labels = {
        "[-1,0)", "[0,1)", "[1,2)", "[2,3)", "[3,4)"
    };
    
    // Storage for event counts: [eta_range_index] = count
    vector<int> gg_like_counts(eta_ranges.size(), 0);
    vector<int> qq_like_counts(eta_ranges.size(), 0);
    vector<int> gq_like_counts(eta_ranges.size(), 0);
    vector<int> intermediate_counts(eta_ranges.size(), 0);
    vector<int> total_counts(eta_ranges.size(), 0);
    
    // Event counters
    int events_processed = 0;
    int events_with_jets = 0;
    int events_selected = 0;
    int total_gg_like = 0;
    int total_qq_like = 0;
    int total_gq_like = 0;
    int total_intermediate = 0;
    
    // Process events
    for (Long64_t i = 0; i < n_entries; ++i) {
        
        if (i % report_every == 0) {
            dout << "Processing event " << i << "/" << n_entries 
                 << " (" << (100*i/n_entries) << "%)\r" << flush;
        }
        
        input_tree->GetEntry(i);
        
        // Check if vectors are valid
        if (!px || !py || !pz || !energy || !eta) continue;
        if (px->size() != py->size() || px->size() != pz->size() || 
            px->size() != energy->size() || px->size() != eta->size()) continue;
        
        events_processed++;
        
        // Reconstruct jets
        JetReconstructionResult jet_result = reconstructJets(*px, *py, *pz, *energy, *eta);
        vector<PseudoJet>& jets = jet_result.jets;
        
        // For event shape analysis: require exactly 2 jets and dijet cuts
        if (jets.size() < 2 || !passesDijetCuts(jets)) continue;
        
        // Keep only the leading 2 jets for dijet analysis
        jets.resize(2);
        
        events_with_jets++;
        events_selected++;
        
        // Classify event shape
        string event_shape = classifyEventShape(jets, *jet_result.cs);
        
        // Count by event shape type
        if (event_shape == "GG-like") {
            total_gg_like++;
        } else if (event_shape == "QQ-like") {
            total_qq_like++;
        } else if (event_shape == "GQ-like") {
            total_gq_like++;
        } else {
            total_intermediate++;
        }
        
        // Determine eta range based on dijet system
        // Use the average eta of the two jets to define event eta
        double dijet_eta = (jets[0].eta() + jets[1].eta()) / 2.0;
        
        // Find which eta range this event belongs to
        for (size_t range_idx = 0; range_idx < eta_ranges.size(); ++range_idx) {
            if (dijet_eta >= eta_ranges[range_idx].first && 
                dijet_eta < eta_ranges[range_idx].second) {
                
                if (event_shape == "GG-like") {
                    gg_like_counts[range_idx]++;
                } else if (event_shape == "QQ-like") {
                    qq_like_counts[range_idx]++;
                } else if (event_shape == "GQ-like") {
                    gq_like_counts[range_idx]++;
                } else {
                    intermediate_counts[range_idx]++;
                }
                total_counts[range_idx]++;
                break;
            }
        }
    }
    
    dout << "\n\nProcessing completed!\n";
    
    // Print results
    dout << "\n=============================================================================\n";
    dout << "EVENT SHAPE CLASSIFICATION RESULTS\n";
    dout << "=============================================================================\n";
    
    dout << "\nOverall Statistics:\n";
    dout << "  Events processed: " << events_processed << "\n";
    dout << "  Events selected: " << events_selected 
         << " (" << (100.0*events_selected/events_processed) << "%)\n";
    dout << "  Total GG-like events: " << total_gg_like << "\n";
    dout << "  Total QQ-like events: " << total_qq_like << "\n";
    dout << "  Total GQ-like events: " << total_gq_like << "\n";
    dout << "  Intermediate events: " << total_intermediate << "\n";
    
    // Print event counts by eta range
    dout << "\nEvent Counts by Eta Range:\n";
    dout << "--------------------------\n";
    dout << setw(12) << "Eta Range" << setw(12) << "GG-like" << setw(12) << "QQ-like" << setw(12) << "GQ-like" << setw(15) << "Intermediate" << setw(12) << "Total" << "\n";
    dout << string(80, '-') << "\n";
    
    for (size_t range_idx = 0; range_idx < eta_ranges.size(); ++range_idx) {
        int total_classified = gg_like_counts[range_idx] + qq_like_counts[range_idx] + 
                              gq_like_counts[range_idx] + intermediate_counts[range_idx];
        dout << setw(12) << eta_range_labels[range_idx]
             << setw(12) << gg_like_counts[range_idx]
             << setw(12) << qq_like_counts[range_idx]
             << setw(12) << gq_like_counts[range_idx]
             << setw(15) << intermediate_counts[range_idx]
             << setw(12) << total_classified << "\n";
    }
    
    // Calculate and print relative fractions
    dout << "\nRelative Event Shape Fractions by Eta Range:\n";
    dout << "--------------------------------------------\n";
    dout << setw(12) << "Eta Range" << setw(15) << "GG-like_frac" << setw(15) << "QQ-like_frac" << setw(15) << "GQ-like_frac" << "\n";
    dout << string(70, '-') << "\n";
    
    for (size_t range_idx = 0; range_idx < eta_ranges.size(); ++range_idx) {
        if (total_counts[range_idx] > 0) {
            double gg_fraction = (double)gg_like_counts[range_idx] / total_counts[range_idx];
            double qq_fraction = (double)qq_like_counts[range_idx] / total_counts[range_idx];
            double gq_fraction = (double)gq_like_counts[range_idx] / total_counts[range_idx];
            
            dout << setw(12) << eta_range_labels[range_idx]
                 << setw(15) << setprecision(4) << fixed << gg_fraction
                 << setw(15) << setprecision(4) << fixed << qq_fraction
                 << setw(15) << setprecision(4) << fixed << gq_fraction << "\n";
        } else {
            dout << setw(12) << eta_range_labels[range_idx]
                 << setw(15) << "0.0000"
                 << setw(15) << "0.0000"
                 << setw(15) << "0.0000" << "\n";
        }
    }
    
    dout << "\n=============================================================================\n";
    dout << "ANALYSIS COMPLETED SUCCESSFULLY!\n";
    dout << "=============================================================================\n";
    
    input_file->Close();
    return 0;
}