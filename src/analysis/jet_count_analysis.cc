// =============================================================================
// Jet Counting Analysis for EIC Photoproduction Events
// =============================================================================
// Description: Counts jets in different eta ranges and categorizes by 
//              subprocess type (QQ, GG, GQ) with relative fractions
//
// Input: ROOT file with photoproduction events
// Output: Console and log file with jet counts and relative fractions
// =============================================================================

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
// ANALYSIS CONFIGURATION (FOLLOWING jetreco.cc APPROACH)
// =============================================================================

// Analysis mode selection (following jetreco.cc)
const bool DIJET_ONLY = true;    // Set to true for dijet analysis, false for all jets

// Jet clustering parameters (from jetreco.cc)
const double R = 1.0;             // Jet radius
const double etMin = 0.0;         // Minimum ET for any jet (base cut)
const double etaMin = -4.0;       // Minimum eta for any jet
const double etaMax = 4.0;        // Maximum eta for any jet

// Dijet-specific cuts (from jetreco.cc)
const double leading_jet_et_min = 10.0;      // Minimum ET for leading jet (GeV)
const double subleading_jet_et_min = 7.0;    // Minimum ET for subleading jet (GeV)
const double leading_jet_eta_max = 4.5;      // Maximum |eta| for leading jet
const double subleading_jet_eta_max = 4.5;   // Maximum |eta| for subleading jet

// Particle selection cuts (from jetreco.cc)
const double particle_pt_min = 0.1;     // Minimum particle pT (GeV)
const double particle_eta_max = 5.0;    // Maximum particle |eta|

// Progress reporting
const int report_every = 100000;    // Report progress every N events

// Input file (from your specification)
const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt5GeV/eic141_pt5/eic141.root";

// =============================================================================
// OUTPUT FILENAME GENERATION
// =============================================================================
string generateLogFilename(const string& input_filename) {
    string baseFileName = input_filename.substr(input_filename.find_last_of("/") + 1);
    baseFileName = baseFileName.substr(0, baseFileName.find_last_of("."));  // Remove .root extension
    return "jet_count_analysis_" + baseFileName + ".log";
}

const string log_filename = generateLogFilename(input_filename);

// =============================================================================
// JET RECONSTRUCTION FUNCTION (FROM jetreco.cc)
// =============================================================================
vector<PseudoJet> reconstructJets(const vector<float>& px, const vector<float>& py, 
                                  const vector<float>& pz, const vector<float>& energy,
                                  const vector<float>& eta) {
    
    // Check input vector sizes
    if (px.empty() || px.size() != py.size() || px.size() != pz.size() || 
        px.size() != energy.size() || px.size() != eta.size()) {
        return vector<PseudoJet>();
    }
    
    // Create input particles for FastJet
    vector<PseudoJet> input_particles;
    
    for (size_t i = 0; i < px.size(); ++i) {
        // Apply particle-level cuts (from jetreco.cc)
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
    
    // Define anti-kT algorithm (from jetreco.cc)
    JetDefinition jet_def(antikt_algorithm, R);
    
    // Run clustering
    ClusterSequence cs(input_particles, jet_def);
    
    // Get jets and apply cuts
    vector<PseudoJet> jets = cs.inclusive_jets();
    vector<PseudoJet> selected_jets;
    
    for (const auto& jet : jets) {
        double et = jet.Et();
        double jet_eta = jet.eta();
        
        // Apply basic jet-level cuts (from jetreco.cc)
        if (et > etMin && jet_eta > etaMin && jet_eta < etaMax) {
            selected_jets.push_back(jet);
        }
    }
    
    // Sort jets by ET (descending) (from jetreco.cc)
    sort(selected_jets.begin(), selected_jets.end(), 
         [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
    
    return selected_jets;
}

// =============================================================================
// DIJET SELECTION FUNCTION (FROM jetreco.cc)
// =============================================================================
bool passesDijetCuts(const vector<PseudoJet>& jets) {
    // Need at least 2 jets
    if (jets.size() < 2) return false;
    
    PseudoJet leading_jet = jets[0];
    PseudoJet subleading_jet = jets[1];
    
    // Leading jet cuts (from jetreco.cc)
    if (leading_jet.Et() < leading_jet_et_min) return false;
    if (abs(leading_jet.eta()) > leading_jet_eta_max) return false;
    
    // Subleading jet cuts (from jetreco.cc)
    if (subleading_jet.Et() < subleading_jet_et_min) return false;
    if (abs(subleading_jet.eta()) > subleading_jet_eta_max) return false;
    
    return true;
}

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================
int main() {

    // Create dual output for console and file logging
    DualOutput dout(log_filename);

    // Open input file
    TFile* input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file || input_file->IsZombie()) {
        cerr << "ERROR: Cannot open input file: " << input_filename << endl;
        return 1;
    }
    
    dout << "=============================================================================\n";
    dout << "JET COUNTING ANALYSIS FOR EIC PHOTOPRODUCTION EVENTS\n";
    dout << "=============================================================================\n";
    dout << "Analysis Parameters:\n";
    dout << "  Analysis Mode:    " << (DIJET_ONLY ? "DIJET ONLY" : "ALL JETS") << "\n";
    dout << "  Jet Algorithm:    Anti-kT\n";
    dout << "  Jet Radius:       " << R << "\n";
    dout << "  Min Jet ET:       " << etMin << " GeV\n";
    dout << "  Jet Eta Range:    [" << etaMin << ", " << etaMax << "]\n";
    dout << "  Min Particle pT:  " << particle_pt_min << " GeV\n";
    
    if (DIJET_ONLY) {
        dout << "\nDijet Selection Cuts:\n";
        dout << "  Leading Jet ET:   > " << leading_jet_et_min << " GeV\n";
        dout << "  Subleading Jet ET:> " << subleading_jet_et_min << " GeV\n";
        dout << "  Leading Jet |η|:  < " << leading_jet_eta_max << "\n";
        dout << "  Subleading |η|:   < " << subleading_jet_eta_max << "\n";
    }
    dout << "=============================================================================\n\n";
    
    dout << "Successfully opened input file: " << input_filename << "\n";
    dout << "Log file: " << log_filename << "\n\n";
    
    // Define event categories to process (focusing on subprocess types)
    vector<string> categories = {"QQ_Events", "GG_Events", "GQ_Events"};
    
    // Define eta ranges for counting: -1 to 0, 0 to 1, 1 to 2, 2 to 3, 3 to 4
    vector<pair<double, double>> eta_ranges = {
        {-1.0, 0.0}, {0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}
    };
    
    // Storage for jet counts: [category][eta_range_index] = count
    map<string, vector<int>> jet_counts;
    map<string, vector<int>> dijet_counts;
    
    // Initialize storage
    for (const string& category : categories) {
        jet_counts[category].resize(eta_ranges.size(), 0);
        dijet_counts[category].resize(eta_ranges.size(), 0);
    }
    
    // Overall statistics
    int total_events_processed = 0;
    int total_dijets_selected = 0;
    
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
        
        // Set up input branches
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
        
        // Event counters
        int events_processed = 0;
        int events_with_jets = 0;
        int events_selected = 0;
        
        // Process events
        for (Long64_t i = 0; i < n_entries; ++i) {
            
            if (i % report_every == 0) {
                dout << "    Processing event " << i << "/" << n_entries 
                     << " (" << (100*i/n_entries) << "%)\r" << flush;
            }
            
            input_tree->GetEntry(i);
            
            // Check if vectors are valid
            if (!px || !py || !pz || !energy || !eta) continue;
            if (px->size() != py->size() || px->size() != pz->size() || 
                px->size() != energy->size() || px->size() != eta->size()) continue;
            
            events_processed++;
            
            // Reconstruct jets (following jetreco.cc approach)
            vector<PseudoJet> jets = reconstructJets(*px, *py, *pz, *energy, *eta);
            
            // Check event selection based on analysis mode
            bool fill_event = false;
            
            if (DIJET_ONLY) {
                // For dijet analysis: require exactly 2 jets and dijet cuts
                if (jets.size() >= 2 && passesDijetCuts(jets)) {
                    // Keep only the leading 2 jets for dijet analysis
                    jets.resize(2);
                    fill_event = true;
                }
            } else {
                // For all-jets analysis: require at least 1 jet
                if (!jets.empty()) {
                    fill_event = true;
                }
            }
            
            if (!fill_event) continue;
            
            events_with_jets++;
            events_selected++;
            
            // Count jets in each eta range
            for (const auto& jet : jets) {
                double jet_eta = jet.eta();
                
                // Find which eta range this jet belongs to
                for (size_t range_idx = 0; range_idx < eta_ranges.size(); ++range_idx) {
                    if (jet_eta >= eta_ranges[range_idx].first && 
                        jet_eta < eta_ranges[range_idx].second) {
                        jet_counts[category][range_idx]++;
                        break;
                    }
                }
            }
            
            // For dijet analysis, count the dijet event in each eta range
            if (DIJET_ONLY && jets.size() == 2) {
                for (const auto& jet : jets) {
                    double jet_eta = jet.eta();
                    
                    // Find which eta range this jet belongs to
                    for (size_t range_idx = 0; range_idx < eta_ranges.size(); ++range_idx) {
                        if (jet_eta >= eta_ranges[range_idx].first && 
                            jet_eta < eta_ranges[range_idx].second) {
                            dijet_counts[category][range_idx]++;
                            break;
                        }
                    }
                }
            }
        }
        
        dout << "\n  Events processed: " << events_processed << "\n";
        dout << "  Events with jets: " << events_with_jets 
             << " (" << (100.0*events_with_jets/events_processed) << "%)\n";
        dout << "  Events selected: " << events_selected 
             << " (" << (100.0*events_selected/events_processed) << "%)\n\n";
        
        // Update global counters
        total_events_processed += events_processed;
        total_dijets_selected += events_selected;
    }
    
    // Calculate and print results
    dout << "=============================================================================\n";
    dout << "JET COUNTING RESULTS\n";
    dout << "=============================================================================\n";
    
    // Print jet counts by eta range
    dout << "\nJet Counts by Eta Range:\n";
    dout << "------------------------\n";
    dout << setw(12) << "Eta Range";
    for (const string& category : categories) {
        dout << setw(12) << category;
    }
    dout << setw(12) << "Total" << "\n";
    
    for (size_t range_idx = 0; range_idx < eta_ranges.size(); ++range_idx) {
        dout << setw(12) << ("[" + to_string(eta_ranges[range_idx].first) + 
                             "," + to_string(eta_ranges[range_idx].second) + ")");
        
        int total_for_range = 0;
        for (const string& category : categories) {
            dout << setw(12) << jet_counts[category][range_idx];
            total_for_range += jet_counts[category][range_idx];
        }
        dout << setw(12) << total_for_range << "\n";
    }
    
    // Calculate and print relative fractions
    dout << "\nRelative Fractions by Eta Range:\n";
    dout << "--------------------------------\n";
    dout << setw(12) << "Eta Range";
    for (const string& category : categories) {
        dout << setw(12) << (category + "_frac");
    }
    dout << "\n";
    
    for (size_t range_idx = 0; range_idx < eta_ranges.size(); ++range_idx) {
        dout << setw(12) << ("[" + to_string(eta_ranges[range_idx].first) + 
                             "," + to_string(eta_ranges[range_idx].second) + ")");
        
        // Calculate total for this eta range
        int total_for_range = 0;
        for (const string& category : categories) {
            total_for_range += jet_counts[category][range_idx];
        }
        
        // Print fractions
        for (const string& category : categories) {
            if (total_for_range > 0) {
                double fraction = (double)jet_counts[category][range_idx] / total_for_range;
                dout << setw(12) << setprecision(4) << fixed << fraction;
            } else {
                dout << setw(12) << "0.0000";
            }
        }
        dout << "\n";
    }
    
    if (DIJET_ONLY) {
        dout << "\nDijet Event Counts by Eta Range (jets in each event):\n";
        dout << "----------------------------------------------------\n";
        dout << setw(12) << "Eta Range";
        for (const string& category : categories) {
            dout << setw(12) << category;
        }
        dout << setw(12) << "Total" << "\n";
        
        for (size_t range_idx = 0; range_idx < eta_ranges.size(); ++range_idx) {
            dout << setw(12) << ("[" + to_string(eta_ranges[range_idx].first) + 
                                 "," + to_string(eta_ranges[range_idx].second) + ")");
            
            int total_for_range = 0;
            for (const string& category : categories) {
                dout << setw(12) << dijet_counts[category][range_idx];
                total_for_range += dijet_counts[category][range_idx];
            }
            dout << setw(12) << total_for_range << "\n";
        }
    }
    
    dout << "\n=============================================================================\n";
    dout << "OVERALL STATISTICS\n";
    dout << "=============================================================================\n";
    dout << "Total events processed: " << total_events_processed << "\n";
    dout << "Total events selected: " << total_dijets_selected 
         << " (" << (100.0*total_dijets_selected/total_events_processed) << "%)\n";
    
    dout << "\nAnalysis completed successfully!\n";
    dout << "Results saved to: " << log_filename << "\n";
    dout << "=============================================================================\n";
    
    input_file->Close();
    return 0;
}