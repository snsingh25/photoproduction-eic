// =============================================================================
// JET SHAPE CALCULATION
// =============================================================================// =============================================================================
// Efficiency and Purity Analysis for Quark/Gluon Jet Classification in Dijets
// =============================================================================
// Description: Analyzes subprocess contributions to thick/thin jet classification
//              using integrated jet shapes in dijet events from HERA data
//              EVENT-LEVEL EFFICIENCY: Both jets must satisfy expected condition
//
// Input: ROOT file with photoproduction events and processType flags
// Output: Detailed efficiency and purity analysis with subprocess breakdown
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1F.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <map>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace fastjet;

// =============================================================================
// ANALYSIS CONFIGURATION
// =============================================================================
// Jet shape cuts for classification (based on HERA studies and your paper)
const double THIN_JET_CUT = 0.8;    // Ψ(r=0.3) > 0.8 for thin jets (quark-like)
const double THICK_JET_CUT = 0.6;   // Ψ(r=0.3) < 0.6 for thick jets (gluon-like)
const double JET_SHAPE_RADIUS = 0.3;

// Dijet analysis parameters
const double R = 1.0;               // Jet radius
const double leading_jet_etMin = 10.0;    // Minimum ET for leading jet (10 GeV)
const double subleading_jet_etMin = 7.0;  // Minimum ET for subleading jet (7 GeV)
const double etaMin = -1.0;         // Minimum eta
const double etaMax = 0.0;          // Maximum eta

// Input file — override by passing a path as argv[1]. Output log name is
// derived from the input filename and written to the current working
// directory (so callers can redirect into a per-sample output folder).
const string default_input_filename =
    "/Users/siddharthsingh/Analysis/repos/photoproduction-eic/data/"
    "allevents_pt7GeV/hera300_pT7/hera300_pT7.root";

// Populated in main() after argv parsing.
string input_filename;
string log_filename;

static string formatFloat2(double value) {
    stringstream ss;
    ss << fixed << setprecision(2) << value;
    return ss.str();
}

static string makeLogFilename(const string& input_path,
                              double lead_et_min, double sublead_et_min,
                              double eta_lo, double eta_hi) {
    string base = input_path.substr(input_path.find_last_of("/") + 1);
    base = base.substr(0, base.find_last_of("."));
    return base +
           "_eT_" + to_string((int)lead_et_min) + "_" + to_string((int)sublead_et_min) +
           "_eta_" + formatFloat2(eta_lo) + "_" + formatFloat2(eta_hi) +
           "_dijet_EventEff.log";
}

// Progress reporting
const int report_every = 100000;    // Report progress every N events

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

// Add a function to check strict dijet selection
// =============================================================================
// DIJET SELECTION FUNCTION
// =============================================================================
bool passesDijetCuts(const vector<PseudoJet>& jets, double lead_et_min, double sublead_et_min, 
                     double eta_min, double eta_max) {
    // Must have exactly 2 jets
    if (jets.size() != 2) return false;
    
    PseudoJet leading_jet = jets[0];
    PseudoJet subleading_jet = jets[1];
    
    // Leading jet cuts
    if (leading_jet.Et() < lead_et_min) return false;
    if (leading_jet.eta() < eta_min || leading_jet.eta() > eta_max) return false;
    
    // Subleading jet cuts
    if (subleading_jet.Et() < sublead_et_min) return false;
    if (subleading_jet.eta() < eta_min || subleading_jet.eta() > eta_max) return false;
    
    return true;
}

double calculateIntegratedJetShape(const vector<PseudoJet>& constituents, 
                                   const PseudoJet& jet, double r) {
    double energy_within_r = 0.0;
    
    for (const auto& constituent : constituents) {
        double dphi = jet.phi() - constituent.phi();
        double drap = jet.rapidity() - constituent.rapidity();
        double const_r = sqrt(dphi*dphi + drap*drap);
        
        if (const_r <= r) {
            energy_within_r += constituent.Et();
        }
    }
    
    return (jet.Et() > 0.0) ? energy_within_r / jet.Et() : 0.0;
}

// =============================================================================
// CLASSIFICATION FUNCTIONS
// =============================================================================

enum JetType {
    UNCLASSIFIED = 0,
    THIN_JET = 1,      // Quark-like
    THICK_JET = 2      // Gluon-like
};

JetType classifyJetByShape(double integrated_shape) {
    if (integrated_shape > THIN_JET_CUT) {
        return THIN_JET;
    } else if (integrated_shape < THICK_JET_CUT) {
        return THICK_JET;
    } else {
        return UNCLASSIFIED;
    }
}

// =============================================================================
// ANALYSIS RESULTS STRUCTURE
// =============================================================================

struct AnalysisResults {
    // Total counts
    int total_events = 0;
    int total_jets = 0;  // Should be 2 * total_events for dijets
    
    // Event-level classification patterns
    map<int, int> events_TT;     // Both jets classified as Thin
    map<int, int> events_TkTk;   // Both jets classified as Thick  
    map<int, int> events_TTk;    // Leading Thin, Subleading Thick [T-Tk]
    map<int, int> events_TkT;    // Leading Thick, Subleading Thin [Tk-T]
    map<int, int> events_other;  // Any other combination (includes unclassified)
    
    // Individual jet classifications by subprocess
    map<int, int> thin_jets_by_subprocess;
    map<int, int> thick_jets_by_subprocess;
    map<int, int> unclassified_jets_by_subprocess;
    
    // Total jets per subprocess (for efficiency calculation)
    map<int, int> total_jets_by_subprocess;
    
    void printResults(DualOutput& dout) {
        dout << "\n=============================================================================\n";
        dout << "EFFICIENCY AND PURITY ANALYSIS RESULTS (EVENT-LEVEL EFFICIENCY)\n";
        dout << "=============================================================================\n";
        
        dout << "Classification Cuts Used:\n";
        dout << "  Thin Jets (Quark-like):  Ψ(r=" << JET_SHAPE_RADIUS << ") > " << THIN_JET_CUT << "\n";
        dout << "  Thick Jets (Gluon-like): Ψ(r=" << JET_SHAPE_RADIUS << ") < " << THICK_JET_CUT << "\n\n";
        
        dout << "Overall Statistics:\n";
        dout << "  Total Events Analyzed: " << total_events << "\n";
        dout << "  Total Jets Analyzed:   " << total_jets << " (should be " << (2*total_events) << ")\n\n";
        
        // Event-level patterns
        dout << "=============================================================================\n";
        dout << "EVENT-LEVEL CLASSIFICATION PATTERNS\n";
        dout << "=============================================================================\n";
        dout << "Format: [Leading-Subleading] classification patterns\n\n";
        
        dout << "Subprocess | Total Events | [T-T] | [Tk-Tk] | [T-Tk] | [Tk-T] | Other | Dominant Pattern\n";
        dout << "-----------|-------------|-------|---------|--------|--------|-------|------------------\n";
        
        for (auto& pair : total_jets_by_subprocess) {
            int subprocess = pair.first;
            int total_events_sub = pair.second / 2;  // Convert back to events
            int tt = events_TT[subprocess];
            int tktk = events_TkTk[subprocess];
            int ttk = events_TTk[subprocess];  // [T-Tk] pattern
            int tkt = events_TkT[subprocess];  // [Tk-T] pattern
            int other = events_other[subprocess];
            
            string name = getSubprocessName(subprocess);
            string dominant = getDominantPattern(tt, tktk, ttk, tkt, other);
            
            dout << setw(10) << name << " | " << setw(11) << total_events_sub << " | "
                 << setw(5) << tt << " | " << setw(7) << tktk << " | " 
                 << setw(6) << ttk << " | " << setw(6) << tkt << " | "
                 << setw(5) << other << " | " << dominant << "\n";
        }
        
        // MODIFIED EFFICIENCY ANALYSIS - EVENT-LEVEL
        dout << "\n=============================================================================\n";
        dout << "EFFICIENCY ANALYSIS (EVENT-LEVEL, BOTH JETS MUST SATISFY CONDITION)\n";
        dout << "=============================================================================\n";
        dout << "Event-level efficiency: Both jets must satisfy the expected condition\n";
        dout << "Expected: QQ→[T-T], GG→[Tk-Tk], GQ→[T-Tk] or [Tk-T] (mixed)\n\n";
        
        dout << "Subprocess | Total Events | Correctly Classified | Efficiency | Performance\n";
        dout << "-----------|-------------|---------------------|------------|-------------\n";
        
        for (auto& pair : total_jets_by_subprocess) {
            int subprocess = pair.first;
            int total_events_sub = pair.second / 2;  // Convert jets back to events
            int correctly_classified = 0;
            string expected_pattern = "";
            
            if (subprocess == 1) {  // QQ subprocess
                correctly_classified = events_TT[subprocess];  // Both jets Thin
                expected_pattern = "[T-T] Both Thin";
            } else if (subprocess == 2) {  // GG subprocess  
                correctly_classified = events_TkTk[subprocess];  // Both jets Thick
                expected_pattern = "[Tk-Tk] Both Thick";
            } else if (subprocess == 3) {  // GQ subprocess
                correctly_classified = events_TTk[subprocess] + events_TkT[subprocess];  // Either mixed pattern
                expected_pattern = "[T-Tk] or [Tk-T] Mixed";
            }
            
            string name = getSubprocessName(subprocess);
            double efficiency = (total_events_sub > 0) ? (double)correctly_classified / total_events_sub : 0.0;
            string performance = getPerformanceRating(efficiency);
            
            dout << setw(10) << name << " | " << setw(11) << total_events_sub << " | "
                 << setw(19) << correctly_classified << " | " << setw(10) << fixed << setprecision(3) 
                 << efficiency << " | " << performance << "\n";
        }
        
        dout << "\nDetailed Expected Performance:\n";
        
        // Calculate individual subprocess efficiencies
        double qq_efficiency = 0.0, gg_efficiency = 0.0, gq_efficiency = 0.0;
        
        for (auto& pair : total_jets_by_subprocess) {
            int subprocess = pair.first;
            int total_events_sub = pair.second / 2;
            
            if (total_events_sub > 0) {
                if (subprocess == 1) {  // QQ
                    qq_efficiency = (double)events_TT[subprocess] / total_events_sub;
                } else if (subprocess == 2) {  // GG
                    gg_efficiency = (double)events_TkTk[subprocess] / total_events_sub;
                } else if (subprocess == 3) {  // GQ
                    gq_efficiency = (double)(events_TTk[subprocess] + events_TkT[subprocess]) / total_events_sub;
                }
            }
        }
        
        dout << "  QQ → [T-T] Both Thin efficiency:  " << setprecision(3) << fixed << qq_efficiency 
             << " (" << (qq_efficiency*100) << "%) " << getPerformanceRating(qq_efficiency) << "\n";
        dout << "  GG → [Tk-Tk] Both Thick efficiency: " << gg_efficiency 
             << " (" << (gg_efficiency*100) << "%) " << getPerformanceRating(gg_efficiency) << "\n";
        dout << "  GQ → [T-Tk] or [Tk-T] Mixed efficiency: " << gq_efficiency 
             << " (" << (gq_efficiency*100) << "%) " << getPerformanceRating(gq_efficiency) << "\n";
        
        // Purity analysis (jet-level for sample composition)
        dout << "\n=============================================================================\n";
        dout << "PURITY ANALYSIS (JET-LEVEL SAMPLE COMPOSITION)\n";
        dout << "=============================================================================\n";
        
        // Calculate total thin and thick jets
        int total_thin_jets = 0, total_thick_jets = 0;
        for (auto& pair : thin_jets_by_subprocess) {
            total_thin_jets += pair.second;
        }
        for (auto& pair : thick_jets_by_subprocess) {
            total_thick_jets += pair.second;
        }
        
        dout << "Total Thin Jets Classified: " << total_thin_jets << "\n";
        dout << "Thin Jet Sample Composition:\n";
        for (auto& pair : thin_jets_by_subprocess) {
            int subprocess = pair.first;
            int count = pair.second;
            string name = getSubprocessName(subprocess);
            double purity = (total_thin_jets > 0) ? (double)count / total_thin_jets : 0.0;
            
            dout << "  " << name << " jets: " << count << " (" 
                 << setprecision(1) << fixed << (purity*100) << "%)\n";
        }
        
        dout << "\nTotal Thick Jets Classified: " << total_thick_jets << "\n";
        dout << "Thick Jet Sample Composition:\n";
        for (auto& pair : thick_jets_by_subprocess) {
            int subprocess = pair.first;
            int count = pair.second;
            string name = getSubprocessName(subprocess);
            double purity = (total_thick_jets > 0) ? (double)count / total_thick_jets : 0.0;
            
            dout << "  " << name << " jets: " << count << " (" 
                 << setprecision(1) << fixed << (purity*100) << "%)\n";
        }
        
        // Physics interpretation
        dout << "\n=============================================================================\n";
        dout << "PHYSICS INTERPRETATION\n";
        dout << "=============================================================================\n";
        
        dout << "Event-Level Classification Performance:\n";
        dout << "  QQ → [T-T] Both Thin:  " << setprecision(3) << fixed << qq_efficiency 
             << " (" << (qq_efficiency*100) << "%) " << getPerformanceRating(qq_efficiency) << "\n";
        dout << "  GG → [Tk-Tk] Both Thick: " << gg_efficiency 
             << " (" << (gg_efficiency*100) << "%) " << getPerformanceRating(gg_efficiency) << "\n";
        dout << "  GQ → [T-Tk] or [Tk-T] Mixed: " << gq_efficiency 
             << " (" << (gq_efficiency*100) << "%) " << getPerformanceRating(gq_efficiency) << "\n";
        
        dout << "\nExpected vs Observed:\n";
        dout << "  - QQ events should produce two quark jets → Both classified as Thin [T-T]\n";
        dout << "  - GG events should produce two gluon jets → Both classified as Thick [Tk-Tk]\n";
        dout << "  - GQ events should produce mixed jets → Either [T-Tk] or [Tk-T] patterns\n";
        
        dout << "\n=============================================================================\n";
    }
    
private:
    string getSubprocessName(int processType) {
        switch(processType) {
            case 1: return "QQ";
            case 2: return "GG";
            case 3: return "GQ";
            default: return "Unknown";
        }
    }
    
    string getDominantPattern(int tt, int tktk, int ttk, int tkt, int other) {
        int max_val = max({tt, tktk, ttk, tkt, other});
        if (max_val == tt) return "[T-T] Thin-dominated";
        if (max_val == tktk) return "[Tk-Tk] Thick-dominated";
        if (max_val == ttk) return "[T-Tk] Mixed (Q-leading)";
        if (max_val == tkt) return "[Tk-T] Mixed (G-leading)";
        return "Other/Unclassified";
    }
    
    string getPerformanceRating(double efficiency) {
        if (efficiency > 0.7) return "✓ Good";
        if (efficiency > 0.5) return "~ Fair";
        return "✗ Poor";
    }
};

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================

int main(int argc, char** argv) {

    input_filename = (argc > 1) ? argv[1] : default_input_filename;
    log_filename = makeLogFilename(input_filename, leading_jet_etMin,
                                    subleading_jet_etMin, etaMin, etaMax);

    // Create dual output for console and file logging
    DualOutput dout(log_filename);
    
    dout << "=============================================================================\n";
    dout << "EVENT-LEVEL EFFICIENCY ANALYSIS FOR DIJET QUARK/GLUON CLASSIFICATION\n";
    dout << "=============================================================================\n";
    dout << "Analysis Parameters:\n";
    dout << "  Analysis Type:        STRICT Dijet Events (exactly 2 jets)\n";
    dout << "  Jet Algorithm:        Anti-kT\n";
    dout << "  Jet Radius:           " << R << "\n";
    dout << "  Leading Jet ET:       > " << leading_jet_etMin << " GeV\n";
    dout << "  Subleading Jet ET:    > " << subleading_jet_etMin << " GeV\n";
    dout << "  Jet Eta Range:        [" << etaMin << ", " << etaMax << "]\n";
    dout << "  Shape Radius:         r = " << JET_SHAPE_RADIUS << "\n";
    dout << "  Thin Jet Cut:         Ψ(r) > " << THIN_JET_CUT << "\n";
    dout << "  Thick Jet Cut:        Ψ(r) < " << THICK_JET_CUT << "\n";
    dout << "  Dijet Requirement:    Exactly 2 jets per event\n";
    dout << "=============================================================================\n\n";
    
    dout << "Input File: " << input_filename << "\n";
    dout << "Log File: " << log_filename << "\n\n";
    
    // Open input file
    TFile* input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file || input_file->IsZombie()) {
        dout << "ERROR: Cannot open input file: " << input_filename << endl;
        return 1;
    }
    
    dout << "Successfully opened input file.\n\n";
    
    // Get Combined_Events directory and tree
    TDirectory* input_dir = (TDirectory*)input_file->Get("Combined_Events");
    if (!input_dir) {
        dout << "ERROR: Combined_Events directory not found!" << endl;
        return 1;
    }
    
    TTree* input_tree = (TTree*)input_dir->Get("Combined_Events");
    if (!input_tree) {
        dout << "ERROR: Combined_Events tree not found!" << endl;
        return 1;
    }
    
    Long64_t n_entries = input_tree->GetEntries();
    dout << "Found " << n_entries << " events in Combined_Events tree.\n\n";
    
    if (n_entries == 0) {
        dout << "No events found, exiting...\n";
        return 1;
    }
    
    // Set up input branches
    vector<float> *px = nullptr, *py = nullptr, *pz = nullptr, *energy = nullptr, *eta = nullptr;
    Int_t processType;
    
    // Check for required branches
    if (!input_tree->GetBranch("px") || !input_tree->GetBranch("py") || 
        !input_tree->GetBranch("pz") || !input_tree->GetBranch("energy") || 
        !input_tree->GetBranch("eta") || !input_tree->GetBranch("processType")) {
        dout << "ERROR: Required branches not found in tree!" << endl;
        return 1;
    }
    
    // Set branch addresses
    input_tree->SetBranchAddress("px", &px);
    input_tree->SetBranchAddress("py", &py);
    input_tree->SetBranchAddress("pz", &pz);
    input_tree->SetBranchAddress("energy", &energy);
    input_tree->SetBranchAddress("eta", &eta);
    input_tree->SetBranchAddress("processType", &processType);
    
    dout << "Successfully set up input branches.\n\n";
    
    // Jet definition for anti-kT algorithm
    JetDefinition jet_def(antikt_algorithm, R);
    
    // Analysis results structure
    AnalysisResults results;
    
    // Counters
    int events_processed = 0;
    int dijet_events = 0;
    
    dout << "Starting event processing...\n\n";
    
    // Event loop
    for (Long64_t entry = 0; entry < n_entries; entry++) {
        input_tree->GetEntry(entry);
        events_processed++;
        
        // Progress reporting
        if (events_processed % report_every == 0) {
            dout << "Processed " << events_processed << " events..." << endl;
        }
        
        // Skip events with unknown process type
        if (processType < 1 || processType > 3) continue;
        
        // Create input particles for FastJet
        vector<PseudoJet> input_particles;
        for (size_t i = 0; i < px->size(); i++) {
            double pt = sqrt((*px)[i]*(*px)[i] + (*py)[i]*(*py)[i]);
            if (pt > 0.1) {  // Minimum pT cut for stability
                PseudoJet particle((*px)[i], (*py)[i], (*pz)[i], (*energy)[i]);
                particle.set_user_index(i);
                input_particles.push_back(particle);
            }
        }
        
        // Cluster jets
        ClusterSequence cs(input_particles, jet_def);
        vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets()); // Minimum pT cut
        
        // Apply loose pre-selection and collect all jets
        vector<PseudoJet> all_jets;
        for (const auto& jet : jets) {
            if (jet.Et() > 1.0 && jet.eta() > etaMin && jet.eta() < etaMax) {
                all_jets.push_back(jet);
            }
        }
        
        // Sort by ET (descending)
        sort(all_jets.begin(), all_jets.end(), 
             [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
        
        // STRICT DIJET SELECTION: Must have exactly 2 jets passing all cuts
        if (!passesDijetCuts(all_jets, leading_jet_etMin, subleading_jet_etMin, etaMin, etaMax)) {
            continue;
        }
        
        // At this point we have exactly 2 jets
        PseudoJet leading_jet = all_jets[0];
        PseudoJet subleading_jet = all_jets[1];
        
        // Event passed all dijet selection criteria
        dijet_events++;
        results.total_events++;
        results.total_jets += 2;  // Count both jets
        
        // Count total jets for this subprocess
        results.total_jets_by_subprocess[processType] += 2;
        
        // Analyze both jets
        vector<PseudoJet> leading_constituents = leading_jet.constituents();
        double leading_shape = calculateIntegratedJetShape(leading_constituents, leading_jet, JET_SHAPE_RADIUS);
        JetType leading_classification = classifyJetByShape(leading_shape);
        
        vector<PseudoJet> subleading_constituents = subleading_jet.constituents();
        double subleading_shape = calculateIntegratedJetShape(subleading_constituents, subleading_jet, JET_SHAPE_RADIUS);
        JetType subleading_classification = classifyJetByShape(subleading_shape);
        
        // Count individual jet classifications
        if (leading_classification == THIN_JET) {
            results.thin_jets_by_subprocess[processType]++;
        } else if (leading_classification == THICK_JET) {
            results.thick_jets_by_subprocess[processType]++;
        } else {
            results.unclassified_jets_by_subprocess[processType]++;
        }
        
        if (subleading_classification == THIN_JET) {
            results.thin_jets_by_subprocess[processType]++;
        } else if (subleading_classification == THICK_JET) {
            results.thick_jets_by_subprocess[processType]++;
        } else {
            results.unclassified_jets_by_subprocess[processType]++;
        }
        
        // Count event-level patterns (distinguish [T-Tk] from [Tk-T])
        if (leading_classification == THIN_JET && subleading_classification == THIN_JET) {
            results.events_TT[processType]++;
        } else if (leading_classification == THICK_JET && subleading_classification == THICK_JET) {
            results.events_TkTk[processType]++;
        } else if (leading_classification == THIN_JET && subleading_classification == THICK_JET) {
            results.events_TTk[processType]++;  // [T-Tk] pattern
        } else if (leading_classification == THICK_JET && subleading_classification == THIN_JET) {
            results.events_TkT[processType]++;  // [Tk-T] pattern
        } else {
            results.events_other[processType]++;
        }
    }
    
    dout << "\n\nEvent processing completed!\n";
    dout << "Events processed: " << events_processed << "\n";
    dout << "Dijet events: " << dijet_events << setprecision(2) << fixed << "(" << (100.0 * dijet_events / events_processed) << ")" << "%\n";
    // dout << "Selection efficiency: " << setprecision(2) << fixed 
    //      << (100.0 * dijet_events / events_processed) << "%\n";
    
    // Print detailed results
    results.printResults(dout);
    
    input_file->Close();
    
    dout << "\nAnalysis completed successfully!\n";
    dout << "Results saved to: " << log_filename << "\n";
    dout << "=============================================================================\n";
    
    return 0;
}