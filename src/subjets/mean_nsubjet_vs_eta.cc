// =============================================================================
// Mean Subjet Multiplicity vs Eta Analysis
// =============================================================================
// Description: Calculates <n_subjet> for QQ and GG jets in different eta bins
//              using exclusive kT clustering with ycut parameter
//
// Input: ROOT file with photoproduction events
// Output: Text file with mean subjet multiplicity per eta bin for plotting
// =============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TDirectoryFile.h"
#include "TROOT.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

using namespace std;
using namespace fastjet;

// =============================================================================
// ANALYSIS CONFIGURATION
// =============================================================================

// Jet clustering parameters
const double R = 1.0;
const double ycut = 0.0005;

// Dijet selection cuts
const double leading_jet_etMin = 10.0;
const double subleading_jet_etMin = 7.0;

// Eta bins for analysis
const int N_ETA_BINS = 4;
const double eta_bin_edges[N_ETA_BINS + 1] = {-1.0, 0.0, 1.0, 1.5, 2};

// Input file
const string input_file = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt7GeV/hera300_pT7/hera300_pT7.root";
// const string input_file = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt5GeV/eic141_pt5/eic141_pT5.root";

// Experiment label (for output)
const string exp_label = "hera300";
const double sqrt_s = 300.0;  // GeV

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

// Find eta bin index (-1 if outside all bins)
int getEtaBin(double eta) {
    for (int i = 0; i < N_ETA_BINS; i++) {
        if (eta >= eta_bin_edges[i] && eta < eta_bin_edges[i + 1]) {
            return i;
        }
    }
    return -1;
}

// Get eta bin center
double getEtaBinCenter(int bin) {
    return (eta_bin_edges[bin] + eta_bin_edges[bin + 1]) / 2.0;
}

// Generate output filename
string generateOutputFilename() {
    // Extract base filename from input path
    string baseFileName = input_file.substr(input_file.find_last_of("/") + 1);
    // Remove .root extension
    baseFileName = baseFileName.substr(0, baseFileName.find_last_of("."));
    
    ostringstream oss;
    oss << baseFileName << "_" 
        << (int)leading_jet_etMin << "_" 
        << (int)subleading_jet_etMin << ".txt";
    return oss.str();
}

// =============================================================================
// MAIN ANALYSIS
// =============================================================================

int main() {
    gROOT->SetBatch(kTRUE);

    cout << "\n========================================" << endl;
    cout << "  Mean Subjet Multiplicity vs Eta" << endl;
    cout << "========================================" << endl;
    cout << "Experiment: " << exp_label << " (sqrt(s) = " << sqrt_s << " GeV)" << endl;
    cout << "Jet Algorithm: anti-kT, R = " << R << endl;
    cout << "Subjet Algorithm: kT, ycut = " << ycut << endl;
    cout << "Leading Jet ET > " << leading_jet_etMin << " GeV" << endl;
    cout << "Subleading Jet ET > " << subleading_jet_etMin << " GeV" << endl;
    cout << "Eta bins: ";
    for (int i = 0; i < N_ETA_BINS; i++) {
        cout << "[" << eta_bin_edges[i] << ", " << eta_bin_edges[i + 1] << ") ";
    }
    cout << "\n========================================\n" << endl;

    // Open ROOT file
    TFile* file = TFile::Open(input_file.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "ERROR: Cannot open file: " << input_file << endl;
        return 1;
    }
    cout << "Successfully opened: " << input_file << "\n" << endl;

    // Access Combined_Events tree
    TDirectoryFile* dir = (TDirectoryFile*)file->Get("Combined_Events");
    if (!dir) {
        cerr << "ERROR: Cannot find Combined_Events directory" << endl;
        return 1;
    }
    
    TTree* tree = (TTree*)dir->Get("Combined_Events");
    if (!tree) {
        cerr << "ERROR: Cannot find Combined_Events tree" << endl;
        return 1;
    }

    // Set up branch addresses
    vector<float>* px = nullptr;
    vector<float>* py = nullptr;
    vector<float>* pz = nullptr;
    vector<float>* energy = nullptr;
    vector<float>* eta = nullptr;
    int processType = 0;

    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("processType", &processType);

    // Storage for subjet counts per eta bin
    // [eta_bin][0=QQ, 1=GG]
    double sum_nsubjet[N_ETA_BINS][2] = {{0}};  // Sum of n_subjet
    int n_jets[N_ETA_BINS][2] = {{0}};           // Number of jets

    // Event counters
    int total_events = 0;
    int dijet_events = 0;
    int qq_events = 0;
    int gg_events = 0;

    Long64_t nEntries = tree->GetEntries();
    cout << "Total events in tree: " << nEntries << endl;

    // =============================================================================
    // EVENT LOOP
    // =============================================================================
    
    for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {
        
        if (iEvent % 100000 == 0 && iEvent > 0) {
            cout << "  Processed " << iEvent << " events..." << endl;
        }

        tree->GetEntry(iEvent);
        total_events++;

        // Check if particles exist
        if (!px || px->empty()) continue;

        // Only process QQ (processType=1) and GG (processType=2)
        if (processType != 1 && processType != 2) continue;
        int proc_idx = (processType == 1) ? 0 : 1;  // 0=QQ, 1=GG

        // Create particles for FastJet
        vector<PseudoJet> particles;
        for (size_t j = 0; j < px->size(); ++j) {
            if ((*energy)[j] <= 0) continue;
            particles.push_back(PseudoJet((*px)[j], (*py)[j], (*pz)[j], (*energy)[j]));
        }

        if (particles.size() < 2) continue;

        // Cluster jets with anti-kT algorithm
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);
        vector<PseudoJet> jets = cs.inclusive_jets();

        // Filter jets by minimum ET (subleading cut as baseline)
        vector<PseudoJet> selected_jets;
        for (const auto& jet : jets) {
            if (jet.Et() > subleading_jet_etMin) {
                selected_jets.push_back(jet);
            }
        }

        // Sort jets by ET (descending)
        sort(selected_jets.begin(), selected_jets.end(),
             [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });

        // Require exactly 2 jets (dijet)
        if (selected_jets.size() != 2) continue;

        // Apply leading and subleading ET cuts
        if (selected_jets[0].Et() < leading_jet_etMin) continue;
        if (selected_jets[1].Et() < subleading_jet_etMin) continue;

        dijet_events++;
        if (processType == 1) qq_events++;
        if (processType == 2) gg_events++;

        // =============================================================================
        // SUBJET ANALYSIS FOR EACH JET
        // =============================================================================
        
        for (const auto& jet : selected_jets) {
            
            // Find eta bin for this jet
            int eta_bin = getEtaBin(jet.eta());
            if (eta_bin < 0) continue;  // Outside eta range

            // Get jet constituents
            vector<PseudoJet> constituents = jet.constituents();
            if (constituents.size() < 2) continue;

            // Recluster constituents with kT algorithm
            JetDefinition kt_def(kt_algorithm, R);
            ClusterSequence cs_kt(constituents, kt_def);

            // Get exclusive subjets with ycut
            vector<PseudoJet> subjets = cs_kt.exclusive_jets_ycut(ycut);
            int n_subjets = subjets.size();

            // Accumulate statistics
            sum_nsubjet[eta_bin][proc_idx] += n_subjets;
            n_jets[eta_bin][proc_idx]++;
        }

    }  // End of event loop

    file->Close();

    // =============================================================================
    // CALCULATE MEAN SUBJET MULTIPLICITIES
    // =============================================================================
    
    cout << "\n========================================" << endl;
    cout << "  Results" << endl;
    cout << "========================================" << endl;
    cout << "Total events processed: " << total_events << endl;
    cout << "Dijet events: " << dijet_events << endl;
    cout << "  QQ events: " << qq_events << endl;
    cout << "  GG events: " << gg_events << endl;
    cout << "\n";

    double mean_nsubjet_QQ[N_ETA_BINS];
    double mean_nsubjet_GG[N_ETA_BINS];

    cout << fixed << setprecision(4);
    cout << "Eta Bin\t\t\t<n_subjet>_QQ\t<n_subjet>_GG\tN_jets_QQ\tN_jets_GG" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    for (int i = 0; i < N_ETA_BINS; i++) {
        mean_nsubjet_QQ[i] = (n_jets[i][0] > 0) ? sum_nsubjet[i][0] / n_jets[i][0] : 0.0;
        mean_nsubjet_GG[i] = (n_jets[i][1] > 0) ? sum_nsubjet[i][1] / n_jets[i][1] : 0.0;

        cout << "[" << eta_bin_edges[i] << ", " << eta_bin_edges[i + 1] << ")\t\t"
             << mean_nsubjet_QQ[i] << "\t\t"
             << mean_nsubjet_GG[i] << "\t\t"
             << n_jets[i][0] << "\t\t"
             << n_jets[i][1] << endl;
    }

    // =============================================================================
    // WRITE OUTPUT FILE
    // =============================================================================
    
    string output_filename = generateOutputFilename();
    ofstream outfile(output_filename);

    outfile << "# Mean Subjet Multiplicity vs Eta" << endl;
    outfile << "# Experiment: " << exp_label << endl;
    outfile << "# sqrt_s = " << sqrt_s << endl;
    outfile << "# R = " << R << endl;
    outfile << "# ycut = " << ycut << endl;
    outfile << "# leading_jet_etMin = " << leading_jet_etMin << endl;
    outfile << "# subleading_jet_etMin = " << subleading_jet_etMin << endl;
    outfile << "# total_events = " << total_events << endl;
    outfile << "# dijet_events = " << dijet_events << endl;
    outfile << "# qq_events = " << qq_events << endl;
    outfile << "# gg_events = " << gg_events << endl;
    outfile << "# Columns: eta_bin_low, eta_bin_high, eta_bin_center, mean_nsubjet_QQ, mean_nsubjet_GG, n_jets_QQ, n_jets_GG" << endl;

    for (int i = 0; i < N_ETA_BINS; i++) {
        outfile << eta_bin_edges[i] << " " 
                << eta_bin_edges[i + 1] << " "
                << getEtaBinCenter(i) << " "
                << mean_nsubjet_QQ[i] << " " 
                << mean_nsubjet_GG[i] << " "
                << n_jets[i][0] << " "
                << n_jets[i][1] << endl;
    }

    outfile.close();
    cout << "\n========================================" << endl;
    cout << "Output saved to: " << output_filename << endl;
    cout << "========================================" << endl;

    return 0;
}