#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TDirectoryFile.h"
#include "TBranch.h"

#include <iostream>
#include <filesystem>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

using namespace std;
using namespace fastjet;

// Analysis parameters
const double R = 1.0;
const double ycut = 0.0005;
const double leading_jet_etMin = 10.0;
const double subleading_jet_etMin = 7.0;

int main(int argc, char** argv) {
    gROOT->SetBatch(kTRUE);

    if (argc < 4) {
        cerr << "Usage: " << argv[0]
             << " <input.root> <etaMin> <etaMax> [out_dir] [alljets|dijets]\n"
             << "  input.root : event-level ROOT (with QQ_Events/GG_Events/GQ_Events trees)\n"
             << "  out_dir    : where to write the .txt (default: current dir)\n"
             << "  selection  : 'alljets' (default) or 'dijets' (require exactly 2 jets)\n";
        return 1;
    }

    const string inputFile = argv[1];
    const double etaMin    = atof(argv[2]);
    const double etaMax    = atof(argv[3]);
    const string outDir    = (argc > 4) ? argv[4] : ".";
    const bool isDijet     = (argc > 5) && string(argv[5]) == "dijets";

    const string expName = std::filesystem::path(inputFile).stem().string();

    // Best-effort sqrt(s) extraction so the plotter can label correctly.
    // Recognises eic64/eic105/eic141/hera300 anywhere in the stem.
    int sqrts_GeV = 0;
    {
        struct Tag { const char* needle; int sqrts; };
        const Tag tags[] = {
            {"hera300", 300}, {"eic141", 141}, {"eic105", 105}, {"eic64", 64}
        };
        for (const auto& t : tags) {
            if (expName.find(t.needle) != string::npos) { sqrts_GeV = t.sqrts; break; }
        }
    }
    
    cout << "\n========================================" << endl;
    cout << "  Subjet Multiplicity Analysis" << endl;
    cout << "========================================" << endl;
    cout << "Experiment : " << expName << endl;
    cout << "Jet Algorithm: anti-kT, R = " << R << endl;
    cout << "Subjet Algorithm: kT, ycut = " << ycut << endl;
    cout << "Leading Jet ET > " << leading_jet_etMin << " GeV" << endl; 
    cout << "Sub-Leading Jet ET > " << subleading_jet_etMin << " GeV" << endl;
    // cout << "All Jets ET > " << jet_etMin << " GeV" << endl;
    cout << "Jet Eta Range: " << etaMin << " < eta < " << etaMax << endl;
    cout << "========================================\n" << endl;

    TFile* file = TFile::Open(inputFile.c_str(), "READ");
    
    if (!file || file->IsZombie()) {
        cerr << "ERROR: Cannot open file: " << inputFile << endl;
        return 1;
    }
    
    cout << "Successfully opened: " << inputFile << "\n" << endl;
    
    // Subprocess names and colors
    vector<string> subprocess_names = {"QQ_Events", "GG_Events", "GQ_Events"};
    vector<string> subprocess_labels = {"Quarks", "Gluons", "Quark-Gluon"};
    
    // Create histograms for statistics collection
    TH1D* h_qq = new TH1D("h_qq", "Subjet Multiplicity", 10, 0, 10);
    TH1D* h_gg = new TH1D("h_gg", "Subjet Multiplicity", 10, 0, 10);
    TH1D* h_gq = new TH1D("h_gq", "Subjet Multiplicity", 10, 0, 10);
    
    vector<TH1D*> histograms = {h_qq, h_gg, h_gq};
    
    // Loop over two subprocesses
    for (int iProc = 0; iProc < 2; iProc++) {
        
        cout << "\n----------------------------------------" << endl;
        cout << "Processing: " << subprocess_names[iProc] << endl;
        cout << "----------------------------------------" << endl;
        
        // Access directory and tree
        TDirectoryFile* dir = (TDirectoryFile*)file->Get(subprocess_names[iProc].c_str());
        if (!dir) {
            cerr << "ERROR: Cannot find directory: " << subprocess_names[iProc] << endl;
            continue;
        }
        
        TTree* tree = (TTree*)dir->Get(subprocess_names[iProc].c_str());
        if (!tree) {
            cerr << "ERROR: Cannot find tree in directory: " << subprocess_names[iProc] << endl;
            continue;
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
        
        int total_events = 0;
        int selected_events = 0;
        int total_jets = 0;
        
        Long64_t nEntries = tree->GetEntries();
        cout << "Total events in tree: " << nEntries << endl;
        
        // Event loop
        for (Long64_t i = 0; i < nEntries; ++i) {

            // if (i >= 10000) break;
            
            // if (i % 50000 == 0 && i > 0) {
            //     cout << "  Processed " << i << " events..." << endl;
            // }
            
            tree->GetEntry(i);
            total_events++;
            
            // Check if particles exist
            if (!px || px->empty()) continue;
            
            // Create particles for FastJet
            vector<PseudoJet> particles;
            for (size_t j = 0; j < px->size(); ++j) {
                if ((*energy)[j] <= 0) continue;
                particles.push_back(PseudoJet((*px)[j], (*py)[j], (*pz)[j], (*energy)[j]));
            }
            
            if (particles.empty()) continue;
            
            // Cluster jets with anti-kT algorithm
            JetDefinition jet_def(antikt_algorithm, R);
            ClusterSequence cs(particles, jet_def);
            vector<PseudoJet> jets = cs.inclusive_jets();
            
            vector<PseudoJet> selected_jets; // Stores all the jets needed

            // Filter out the jets whose minimum Et is that of the subleading jet cut
            for (const auto& jet : jets) {
                if (jet.Et() > subleading_jet_etMin) {
                    selected_jets.push_back(jet);
                }
            }
            
            // Sort jets by ET (descending)
            sort(selected_jets.begin(), selected_jets.end(), 
                 [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
            
            // Require at least one jet
            if (selected_jets.empty()) continue;
            // Require Dijets
            if (isDijet) if (selected_jets.size() != 2) continue; // if (selected_jets.size() != 2) then its a dijet
            
            // Apply leading and subleading ET cuts
            if (selected_jets[0].Et() < leading_jet_etMin) continue;  // Leading > 10 GeV
            if (selected_jets[1].Et() < subleading_jet_etMin) continue;  // Subleading > 7 GeV

            if (selected_jets[0].eta() < etaMin || selected_jets[0].eta() > etaMax) continue;
            if (selected_jets[1].eta() < etaMin || selected_jets[1].eta() > etaMax) continue;
            
            selected_events++;
            
            // Loop over ALL selected jets in this event
            for (const auto& jet : selected_jets) {
                
                // Get jet constituents
                vector<PseudoJet> constituents = jet.constituents();
                
                if (constituents.size() < 2) continue; // Need at least 2 constituents
                
                // Recluster constituents with kT algorithm
                JetDefinition kt_def(kt_algorithm, R);
                ClusterSequence cs_kt(constituents, kt_def);
                
                // Get exclusive subjets with ycut
                vector<PseudoJet> subjets = cs_kt.exclusive_jets_ycut(ycut);
                int n_subjets = subjets.size();
                
                // Fill histogram
                histograms[iProc]->Fill(n_subjets);
                total_jets++;
            }
            
        } // End of event loop
        
        cout << "----------------------------------------" << endl;
        cout << "Summary for " << subprocess_labels[iProc] << ":" << endl;
        cout << "  Total events processed: " << total_events << endl;
        cout << "  Events passing cuts: " << selected_events << endl;
        cout << "  Total jets analyzed: " << total_jets << endl;
        cout << "  Jets in histogram: " << histograms[iProc]->GetEntries() << endl;
        cout << "----------------------------------------" << endl;
        
    } // End of subprocess loop
    
    // Normalize histograms
    cout << "\nNormalizing histograms..." << endl;
    for (int i = 0; i < 3; i++) {
        double integral = histograms[i]->Integral();
        if (integral > 0) {
            histograms[i]->Scale(1.0 / integral);
            cout << subprocess_labels[i] << " normalized (integral was " << integral << ")" << endl;
        }
    }

    // Write output data file with proper decimal formatting
    std::filesystem::create_directories(outDir);
    ostringstream filename_stream;
    filename_stream << outDir << "/subjet_multiplicity_data_" << expName << "_"
                    << (int)leading_jet_etMin << "_"
                    << (int)subleading_jet_etMin << "_"
                    << fixed << setprecision(1) << etaMin << "p"
                    << fixed << setprecision(1) << etaMax << ".txt";

    ofstream outfile(filename_stream.str());

    outfile << "# Subjet Multiplicity Data" << endl;
    outfile << "# Experiment: " << expName << endl;
    outfile << "# CMS_Energy: " << sqrts_GeV << endl;
    outfile << "# ycut = " << ycut << endl;
    outfile << "# etaMin = " << etaMin << endl;
    outfile << "# etaMax = " << etaMax << endl;
    outfile << "# Leading_Jet_ET_Min = " << leading_jet_etMin << endl;
    outfile << "# Subleading_Jet_ET_Min = " << subleading_jet_etMin << endl;
    outfile << "# Selection = " << (isDijet ? "dijets" : "alljets") << endl;
    outfile << "# Columns: bin_center, QQ_normalized, GG_normalized" << endl;
    outfile << "# Mean_QQ: " << h_qq->GetMean() << endl;
    outfile << "# Mean_GG: " << h_gg->GetMean() << endl;
    outfile << "# Entries_QQ: " << h_qq->GetEntries() << endl;
    outfile << "# Entries_GG: " << h_gg->GetEntries() << endl;

    for (int bin = 1; bin <= h_qq->GetNbinsX(); bin++) {
        double bin_center = h_qq->GetBinCenter(bin);
        double qq_val = h_qq->GetBinContent(bin);
        double gg_val = h_gg->GetBinContent(bin);
        outfile << bin_center << " " << qq_val << " " << gg_val << endl;
    }
    outfile.close();
    
    cout << "\n========================================" << endl;
    cout << "  Analysis Complete" << endl;
    cout << "========================================" << endl;
    cout << "Data saved to: " << filename_stream.str() << endl;
    
    cout << "\nHistogram Statistics:" << endl;
    cout << "Jet Eta Range: " << etaMin << " < eta < " << etaMax << endl;
    for (int i = 0; i < 2; i++) {
        cout << subprocess_labels[i] << " (Mean n_subjets: " << histograms[i]->GetMean() << ", Total entries: " << histograms[i]->GetEntries() << ")" << endl;
    }
    cout << "========================================\n" << endl;
    
    // Clean up
    for (auto h : histograms) delete h;
    file->Close();
    
    return 0;
}