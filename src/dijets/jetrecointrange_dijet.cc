// =============================================================================
//     Integrated Jet Shape Analysis for Multiple Eta Ranges (Table Format)
// =============================================================================
// Description: Calculates mean integrated jet shape <psi(r)> for r = 0.1 to 1.0
//              across multiple eta bins, outputting results in a table format
//              similar to published ZEUS/H1 data.
//
//              Supports both INCLUSIVE and DIJET analysis modes.
//
// Input:  ROOT file with particle-level photoproduction events
// Output: Console/log table with <psi(r)> ± stat for each eta bin
//
// Author: Analysis code for HERA/EIC jet shape studies
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <map>
#include <fstream>
#include <sstream>

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
        std::cout << data;
        if (fileOpen) {
            logFile << data;
            logFile.flush();
        }
        return *this;
    }
    
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

// Input file - CHANGE THIS FOR DIFFERENT DATASETS
const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/data/allevents_pt7GeV/hera300_pT7/hera300_pT7.root"; // "/Users/siddharthsingh/Analysis/ph-new/evt/hera300_pTHat7_pT0Ref4.0_Q2max1_n1000k.root";

// ===== ANALYSIS MODE =====
// Set to true for DIJET analysis, false for INCLUSIVE analysis
const bool DIJET_MODE = true;

// Jet clustering parameters 
const double R = 1.0;                // Jet radius

// ----- INCLUSIVE MODE CUTS -----
const double jet_etMin_inclusive = 17.0;  // Minimum ET for inclusive jets (GeV)

// ----- DIJET MODE CUTS -----
const double leading_jet_etMin = 10.0;     // Minimum ET for leading jet (GeV)
const double subleading_jet_etMin = 7.0;   // Minimum ET for subleading jet (GeV)

// Inelasticity cuts (ZEUS kinematic selection)
const double y_min = 0.2;            // Minimum inelasticity
const double y_max = 0.85;           // Maximum inelasticity
const bool APPLY_INELASTICITY_CUT = false;  // Set to true to enable

// Eta bins to analyze (matching published table)
struct EtaBin {
    double etaMin;
    double etaMax;
    string label;
};

const vector<EtaBin> eta_bins = {
    {-1.0, 0.0, "-1 < eta < 0"},
    {0.0, 1.0, "0 < eta < 1"},
    {1.0, 1.5, "1 < eta < 1.5"},
    {1.5, 2.5, "1.5 < eta < 2.5"}
};

// Integrated jet shape parameters
const double r_min = 0.1;            // Minimum r value
const double r_max = 1.0;            // Maximum r value  
const double r_step = 0.1;           // Step size for r

// Categories to process
const vector<string> categories = {"QQ_Events", "GG_Events", "GQ_Events", "Combined_Events"};

// =============================================================================
// INTEGRATED JET SHAPE CALCULATION
// =============================================================================
// Calculates psi(r) = E_T(within r) / E_T(total) for a single jet
// Uses FastJet's delta_R for proper phi wrapping

double calculatePsiForJet(const PseudoJet& jet, double r_cone) {
    if (jet.Et() <= 0.0) return -1.0;
    
    double et_within_r = 0.0;
    vector<PseudoJet> constituents = jet.constituents();
    
    for (const auto& constituent : constituents) {
        // FastJet's delta_R handles phi wrapping correctly
        double dR = constituent.delta_R(jet);
        if (dR <= r_cone) {
            et_within_r += constituent.Et();
        }
    }
    
    return et_within_r / jet.Et();
}

// =============================================================================
// STATISTICS STORAGE STRUCTURE
// =============================================================================
struct JetShapeStats {
    vector<double> psi_values;  // All psi values for this (category, eta_bin, r)
    
    double getMean() const {
        if (psi_values.empty()) return 0.0;
        double sum = 0.0;
        for (double v : psi_values) sum += v;
        return sum / psi_values.size();
    }
    
    double getStatError() const {
        if (psi_values.size() < 2) return 0.0;
        double mean = getMean();
        double sum_sq = 0.0;
        for (double v : psi_values) {
            sum_sq += (v - mean) * (v - mean);
        }
        double variance = sum_sq / psi_values.size();
        return sqrt(variance / psi_values.size());  // Standard error of mean
    }
    
    int getCount() const {
        return psi_values.size();
    }
};

// =============================================================================
// LOG FILENAME GENERATION
// =============================================================================
string generateLogFilename(const string& input_filename) {
    string baseFileName = input_filename.substr(input_filename.find_last_of("/") + 1);
    baseFileName = baseFileName.substr(0, baseFileName.find_last_of("."));
    
    string mode_tag = DIJET_MODE ? "_DIJET" : "_INCLUSIVE";
    string et_tag;
    
    if (DIJET_MODE) {
        et_tag = "_Et" + to_string((int)leading_jet_etMin) + "_" + to_string((int)subleading_jet_etMin);
    } else {
        et_tag = "_EtMin" + to_string((int)jet_etMin_inclusive);
    }
    
    return "intjetshape_table_" + baseFileName + mode_tag + et_tag + ".log";
}

// =============================================================================
// MAIN ANALYSIS
// =============================================================================
int main() {
    
    gROOT->SetBatch(kTRUE);
    
    const string log_filename = generateLogFilename(input_filename);
    DualOutput dout(log_filename);
    
    // Print header
    dout << "=============================================================================\n";
    dout << "INTEGRATED JET SHAPE TABLE ANALYSIS\n";
    dout << "=============================================================================\n";
    dout << "Input File:     " << input_filename << "\n";
    dout << "Jet Algorithm:  Anti-kT, R = " << R << "\n";
    dout << "Analysis Mode:  " << (DIJET_MODE ? "DIJET" : "INCLUSIVE") << "\n";
    
    if (DIJET_MODE) {
        dout << "Leading Jet:    E_T > " << leading_jet_etMin << " GeV\n";
        dout << "Subleading Jet: E_T > " << subleading_jet_etMin << " GeV\n";
    } else {
        dout << "Jet ET Cut:     > " << jet_etMin_inclusive << " GeV\n";
    }
    
    dout << "Inelasticity:   " << (APPLY_INELASTICITY_CUT ? 
                                   to_string(y_min) + " < y < " + to_string(y_max) : "No cut") << "\n";
    dout << "r Range:        " << r_min << " to " << r_max << " (step " << r_step << ")\n";
    dout << "Eta Bins:       " << eta_bins.size() << " regions\n";
    dout << "Log File:       " << log_filename << "\n";
    dout << "=============================================================================\n\n";
    
    // Open input file
    TFile* input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file || input_file->IsZombie()) {
        cerr << "ERROR: Cannot open file: " << input_filename << endl;
        return 1;
    }
    
    // Storage: [category][eta_bin_idx][r_idx] -> JetShapeStats
    int n_r_points = int((r_max - r_min) / r_step + 0.5) + 1;
    int n_eta_bins = eta_bins.size();
    
    map<string, vector<vector<JetShapeStats>>> all_stats;
    for (const string& cat : categories) {
        all_stats[cat].resize(n_eta_bins);
        for (int e = 0; e < n_eta_bins; e++) {
            all_stats[cat][e].resize(n_r_points);
        }
    }
    
    // Also track total jet counts per category per eta bin
    map<string, vector<int>> jet_counts;
    for (const string& cat : categories) {
        jet_counts[cat].resize(n_eta_bins, 0);
    }

    // Track events and dijets
    int total_events_processed = 0;
    int total_events_with_jets = 0;
    int total_dijet_events = 0;
    
    // Process each category
    for (const string& category : categories) {
        
        dout << "Processing " << category << "...\n";
        
        // Get directory and tree
        TDirectory* dir = (TDirectory*)input_file->Get(category.c_str());
        if (!dir) {
            dout << "  Warning: Directory " << category << " not found, skipping.\n";
            continue;
        }
        
        TTree* tree = (TTree*)dir->Get(category.c_str());
        if (!tree) {
            dout << "  Warning: Tree " << category << " not found, skipping.\n";
            continue;
        }
        
        Long64_t n_entries = tree->GetEntries();
        dout << "  Found " << n_entries << " events\n";
        
        // Setup branches for particle data
        vector<float> *px = nullptr, *py = nullptr, *pz = nullptr, *energy = nullptr, *eta = nullptr;
        tree->SetBranchAddress("px", &px);
        tree->SetBranchAddress("py", &py);
        tree->SetBranchAddress("pz", &pz);
        tree->SetBranchAddress("energy", &energy);
        tree->SetBranchAddress("eta", &eta);
        
        // Setup branch for inelasticity (if it exists)
        Float_t inelasticity = 0.5;  // Default value in valid range
        bool has_inelasticity = false;
        if (tree->GetBranch("inelasticity")) {
            tree->SetBranchAddress("inelasticity", &inelasticity);
            has_inelasticity = true;
        }
        
        int events_processed = 0;
        int events_with_jets = 0;
        int dijet_events = 0;
        int events_cut_by_y = 0;
        
        // Event loop
        for (Long64_t i = 0; i < n_entries; ++i) {
            
            if (i % 100000 == 0 && i > 0) {
                dout << "    Event " << i << "/" << n_entries << "\r" << flush;
            }
            
            tree->GetEntry(i);
            
            // Apply inelasticity cut (ZEUS: 0.2 < y < 0.85)
            if (APPLY_INELASTICITY_CUT && has_inelasticity) {
                if (inelasticity < y_min || inelasticity > y_max) {
                    events_cut_by_y++;
                    continue;
                }
            }
            
            // Validate input
            if (!px || !py || !pz || !energy || !eta) continue;
            if (px->empty()) continue;
            if (px->size() != py->size() || px->size() != pz->size() || 
                px->size() != energy->size() || px->size() != eta->size()) continue;
            
            events_processed++;
            
            // Build particles for FastJet
            vector<PseudoJet> particles;
            for (size_t j = 0; j < px->size(); ++j) {
                if ((*energy)[j] > 0) {
                    PseudoJet p((*px)[j], (*py)[j], (*pz)[j], (*energy)[j]);
                    p.set_user_index(j);
                    particles.push_back(p);
                }
            }
            
            if (particles.size() < 2) continue;
            
            // Cluster jets with anti-kT
            JetDefinition jet_def(antikt_algorithm, R);
            ClusterSequence cs(particles, jet_def);
            vector<PseudoJet> jets = cs.inclusive_jets();
            
            // Sort by ET (descending)
            sort(jets.begin(), jets.end(), 
                 [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
            
            // =====================================================================
            // JET SELECTION BASED ON MODE
            // =====================================================================
            
            vector<PseudoJet> selected_jets;
            
            if (DIJET_MODE) {
                // ----- DIJET MODE -----
                // Require at least 2 jets
                if (jets.size() < 2) continue;
                
                // Check if leading and subleading jets pass ET cuts
                if (jets[0].Et() < leading_jet_etMin) continue;
                if (jets[1].Et() < subleading_jet_etMin) continue;
                
                // Both jets in the dijet are selected for shape analysis
                selected_jets.push_back(jets[0]);
                selected_jets.push_back(jets[1]);
                
                dijet_events++;
                
            } else {
                // ----- INCLUSIVE MODE -----
                // Select all jets passing ET cut
                for (const auto& jet : jets) {
                    if (jet.Et() >= jet_etMin_inclusive) {
                        selected_jets.push_back(jet);
                    }
                }
            }
            
            if (selected_jets.empty()) continue;
            
            events_with_jets++;
            
            // =====================================================================
            // CALCULATE JET SHAPES FOR SELECTED JETS
            // =====================================================================
            
            for (const auto& jet : selected_jets) {
                
                double jet_eta = jet.eta();
                
                // Find which eta bin this jet belongs to
                for (int e = 0; e < n_eta_bins; e++) {
                    if (jet_eta > eta_bins[e].etaMin && jet_eta <= eta_bins[e].etaMax) {
                        
                        jet_counts[category][e]++;
                        
                        // Calculate psi(r) for each r value
                        for (int r_idx = 0; r_idx < n_r_points; r_idx++) {
                            double r = r_min + r_idx * r_step;
                            double psi = calculatePsiForJet(jet, r);
                            
                            if (psi >= 0.0 && psi <= 1.0) {
                                all_stats[category][e][r_idx].psi_values.push_back(psi);
                            }
                        }
                        
                        break;  // Jet belongs to only one eta bin
                    }
                }
            }
        }
        
        total_events_processed += events_processed;
        total_events_with_jets += events_with_jets;
        total_dijet_events += dijet_events;
        
        dout << "\n  Events processed: " << events_processed;
        if (DIJET_MODE) {
            dout << ", dijet events: " << dijet_events;
        }
        dout << ", with selected jets: " << events_with_jets;
        if (APPLY_INELASTICITY_CUT && has_inelasticity) {
            dout << ", cut by inelasticity: " << events_cut_by_y;
        }
        dout << "\n\n";
    }
    
    input_file->Close();
    
    // ==========================================================================
    // OUTPUT RESULTS IN TABLE FORMAT
    // ==========================================================================
    
    dout << "\n";
    dout << "=============================================================================\n";
    dout << "INTEGRATED JET SHAPE RESULTS: <psi(r)> +/- stat\n";
    dout << "Analysis Mode: " << (DIJET_MODE ? "DIJET" : "INCLUSIVE") << "\n";
    
    if (DIJET_MODE) {
        dout << "Jet Selection: Leading E_T > " << leading_jet_etMin 
             << " GeV, Subleading E_T > " << subleading_jet_etMin << " GeV\n";
    } else {
        dout << "Jet Selection: E_T > " << jet_etMin_inclusive << " GeV\n";
    }
    
    dout << "Jet Algorithm: Anti-kT R=" << R << "\n";
    if (APPLY_INELASTICITY_CUT) {
        dout << "Inelasticity Cut: " << y_min << " < y < " << y_max << "\n";
    }
    dout << "=============================================================================\n\n";
    
    // Print table for each category
    for (const string& category : categories) {
        
        dout << ">>> " << category << " <<<\n";
        dout << string(100, '-') << "\n";
        
        // Header row with eta bins
        dout << setw(6) << "r";
        for (int e = 0; e < n_eta_bins; e++) {
            dout << setw(24) << eta_bins[e].label;
        }
        dout << "\n";
        
        // Separator
        dout << setw(6) << "---";
        for (int e = 0; e < n_eta_bins; e++) {
            dout << setw(24) << "--------------------";
        }
        dout << "\n";
        
        // Data rows
        for (int r_idx = 0; r_idx < n_r_points; r_idx++) {
            double r = r_min + r_idx * r_step;
            
            dout << setw(6) << fixed << setprecision(1) << r;
            
            for (int e = 0; e < n_eta_bins; e++) {
                JetShapeStats& stats = all_stats[category][e][r_idx];
                
                if (stats.getCount() > 0) {
                    stringstream ss;
                    ss << fixed << setprecision(4) << stats.getMean() 
                       << " +/- " << setprecision(4) << stats.getStatError();
                    dout << setw(24) << ss.str();
                } else {
                    dout << setw(24) << "N/A";
                }
            }
            dout << "\n";
        }
        
        // Jet counts
        dout << "\nJet counts per eta bin:\n";
        for (int e = 0; e < n_eta_bins; e++) {
            dout << "  " << eta_bins[e].label << ": " << jet_counts[category][e] << " jets\n";
        }
        dout << "\n\n";
    }
    
    // ==========================================================================
    // SUMMARY STATISTICS
    // ==========================================================================
    
    dout << "=============================================================================\n";
    dout << "SUMMARY\n";
    dout << "=============================================================================\n";
    dout << "Total events processed: " << total_events_processed << "\n";
    dout << "Events with selected jets: " << total_events_with_jets << "\n";
    if (DIJET_MODE) {
        dout << "Dijet events: " << total_dijet_events << "\n";
        dout << "Dijet efficiency: " << fixed << setprecision(2) 
             << (100.0 * total_dijet_events / total_events_processed) << "%\n";
    }
    
    // Total jets across all categories (divide by 4 since Combined includes all)
    int total_jets = 0;
    for (int e = 0; e < n_eta_bins; e++) {
        total_jets += jet_counts["Combined_Events"][e];
    }
    dout << "Total jets analyzed: " << total_jets << "\n";
    
    dout << "\n=============================================================================\n";
    dout << "ANALYSIS COMPLETE\n";
    dout << "Log file saved to: " << log_filename << "\n";
    dout << "=============================================================================\n";
    
    return 0;
}