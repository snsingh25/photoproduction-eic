// =============================================================================
// Differential Jet Shapes Analysis for HERA Photoproduction Events
// =============================================================================
// Description: Calculates differential jet shapes for r = 0.1 to 1.0
//              for different subprocess categories (QQ, GG, GQ) and 
//              photoproduction types (Combined, Direct, Resolved)
//              with automatic plotting and data storage
//
// Input: ROOT file with photoproduction events
// Output: Console printout of differential jet shape values + ROOT plots
// Author: Siddharth Singh (Enhanced)
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "THStack.h"
#include "TColor.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <map>

using namespace std;
using namespace fastjet;

// =============================================================================
// ANALYSIS CONFIGURATION
// =============================================================================

// **DIJET ANALYSIS SWITCH**
const bool DIJET_ANALYSIS = true;    // Set to true for dijet-only analysis (exactly 2 jets)
                                     // Set to false for inclusive analysis (>=1 jet)

// Jet clustering parameters 
const double R = 0.4;             // Jet radius
const double jet_etMin = 17.0;     // Minimum ET for ALL jets
const double etaMin = -1.0;        // Minimum eta
const double etaMax = 4.0;         // Maximum eta

// Differential jet shape parameters
const double delta_r = 0.10;      // Annulus width (±0.05 around r) it will be divided by 2 in the code
const double r_min = 0.1;         // Minimum r value
const double r_max = 1.0;         // Maximum r value
const double r_step = 0.1;        // Step size for r

// Progress reporting
const int report_every = 1000000;   // Report progress every N events

// Input file path
const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt7GeV/hera1M/hera.root";

// Global data storage for plotting
map<string, vector<double>> global_diff_shapes_data;
const int n_global_r_points = int((r_max - r_min) / r_step) + 1;

// =============================================================================
// DIFFERENTIAL JET SHAPE CALCULATION
// =============================================================================

double calculateDifferentialJetShape(const vector<PseudoJet>& jets, double r) {
    
    // Apply jet requirements based on analysis mode
    if (DIJET_ANALYSIS) {
        // Require exactly 2 jets for dijet analysis
        if (jets.size() != 2) return 0.0;
    } else {
        // Require at least 1 jet for inclusive analysis
        if (jets.size() < 1) return 0.0;
    }
    
    double R1 = r - (delta_r / 2.0);  // Inner radius
    double R2 = r + (delta_r / 2.0);  // Outer radius
    double rho_total = 0.0;
    int n_jets_to_analyze = jets.size();
    
    // For dijet analysis, we analyze both jets
    // For inclusive analysis, we analyze all jets
    for (int i = 0; i < n_jets_to_analyze; i++) {
        
        vector<PseudoJet> constituents = jets[i].constituents();
        double select_const_Et = 0.0;
        
        // Loop over constituents and select ones in annulus [R1, R2]
        for (unsigned int j = 0; j < constituents.size(); j++) {
            
            double dphi = jets[i].phi() - constituents[j].phi();
            double drap = jets[i].rapidity() - constituents[j].rapidity();
            double const_r = sqrt(dphi*dphi + drap*drap);
            
            // Add ET of constituents in the annulus
            if (const_r >= R1 && const_r <= R2) {
                select_const_Et += constituents[j].Et();
            }
        }
        
        // Calculate contribution from this jet
        double jet_contribution = select_const_Et / jets[i].Et();
        rho_total += (jet_contribution / delta_r);
    }
    
    return rho_total;
}

// =============================================================================
// JET RECONSTRUCTION AND SHAPE ANALYSIS FUNCTION
// =============================================================================

vector<double> reconstructJetsAndCalculateShapes(const vector<float>& px, const vector<float>& py, 
                                                 const vector<float>& pz, const vector<float>& energy,
                                                 const vector<float>& eta) {
    
    // Check input vector sizes
    if (px.empty() || px.size() != py.size() || px.size() != pz.size() || 
        px.size() != energy.size() || px.size() != eta.size()) {
        return vector<double>();
    }
    
    // Create input particles for FastJet (include all soft particles)
    vector<PseudoJet> input_particles;
    
    for (size_t i = 0; i < px.size(); ++i) {
        // Only apply minimal cuts - no pT cut to include soft particles
        if (energy[i] <= 0) continue;
        
        // Create PseudoJet
        PseudoJet particle(px[i], py[i], pz[i], energy[i]);
        particle.set_user_index(i);
        input_particles.push_back(particle);
    }
    
    if (input_particles.size() < 2) {
        return vector<double>();
    }
    
    // Define anti-kT algorithm
    JetDefinition jet_def(antikt_algorithm, R);
    
    // Run clustering - ClusterSequence stays in scope for constituent access
    ClusterSequence cs(input_particles, jet_def);
    
    // Get all jets
    vector<PseudoJet> jets = cs.inclusive_jets();
    
    // Sort jets by ET (descending) to identify leading jet
    sort(jets.begin(), jets.end(), 
         [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
    
    // Apply cuts: ALL jets must have ET > 17 GeV and be in eta range
    vector<PseudoJet> selected_jets;
    
    for (const auto& jet : jets) {
        double et = jet.Et();
        double jet_eta = jet.eta();
        
        // Apply same cuts to all jets: ET > 17 GeV and eta in range
        if (et > jet_etMin && jet_eta > etaMin && jet_eta < etaMax) {
            selected_jets.push_back(jet);
        }
    }

    // Apply jet multiplicity requirements based on analysis mode
    if (DIJET_ANALYSIS) {
        // Require exactly 2 jets for dijet analysis
        if (selected_jets.size() != 2) {
            return vector<double>();
        }
    } else {
        // Require at least 1 jet for inclusive analysis
        if (selected_jets.empty()) {
            return vector<double>();
        }
    }
    
    // Calculate differential jet shapes for all r values while ClusterSequence is in scope
    vector<double> shapes;
    int n_r_points = int((r_max - r_min) / r_step) + 1;
    
    for (int r_idx = 0; r_idx < n_r_points; r_idx++) {
        double r = r_min + r_idx * r_step;
        double diff_shape = calculateDifferentialJetShape(selected_jets, r);
        shapes.push_back(diff_shape);
    }
    
    return shapes;
}

// =============================================================================
// DATA STORAGE AND ARRAY GENERATION FUNCTION
// =============================================================================

void storeDataForPlotting(const map<string, vector<vector<double>>>& diff_shapes) {
    
    cout << "\n=============================================================================\n";
    cout << "GENERATING C++ ARRAYS FOR PLOTTING\n";
    cout << "=============================================================================\n\n";
    
    // Category mapping for cleaner output
    map<string, string> category_labels = {
        {"QQ_Events", "quark_values"},
        {"GG_Events", "gluon_values"},
        {"GQ_Events", "gluon_quark_values"},
        {"Combined_Events", "combined_values"},
        {"Resolved_Events", "resolved_values"},
        {"Direct_Events", "direct_values"}
    };
    
    // Generate arrays
    for (const auto& [category, label] : category_labels) {
        if (diff_shapes.find(category) != diff_shapes.end() && 
            !diff_shapes.at(category)[0].empty()) {
            
            cout << "// " << category << "\n";
            cout << "double " << label << "[] = {";
            
            // Calculate averages for each r value
            for (int r_idx = 0; r_idx < n_global_r_points; r_idx++) {
                if (!diff_shapes.at(category)[r_idx].empty()) {
                    double sum = 0.0;
                    for (double val : diff_shapes.at(category)[r_idx]) {
                        sum += val;
                    }
                    double average = sum / diff_shapes.at(category)[r_idx].size();
                    cout << fixed << setprecision(4) << average;
                    
                    // Store for plotting
                    global_diff_shapes_data[category].push_back(average);
                } else {
                    cout << "0.0000";
                    global_diff_shapes_data[category].push_back(0.0);
                }
                
                if (r_idx < n_global_r_points - 1) cout << ", ";
            }
            
            cout << "};\n";
        } else {
            // Create empty array for missing categories
            cout << "// " << category << " (No data)\n";
            cout << "double " << label << "[] = {";
            for (int r_idx = 0; r_idx < n_global_r_points; r_idx++) {
                cout << "0.0000";
                global_diff_shapes_data[category].push_back(0.0);
                if (r_idx < n_global_r_points - 1) cout << ", ";
            }
            cout << "};\n";
        }
    }
    
    cout << "\n=============================================================================\n\n";
}

// =============================================================================
// PLOTTING FUNCTION
// =============================================================================

void createDifferentialJetShapePlot() {
    
    cout << "Creating differential jet shape plots...\n";
    
    // Set global style
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetFrameBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    
    // Set serif fonts (Times-Roman)
    gStyle->SetTextFont(132);     // Times-Roman serif
    gStyle->SetLabelFont(132,"x");
    gStyle->SetLabelFont(132,"y");
    gStyle->SetTitleFont(132,"x");
    gStyle->SetTitleFont(132,"y");
    gStyle->SetLegendFont(132);
    
    // Increase tick sizes
    gStyle->SetTickLength(0.03, "x");
    gStyle->SetTickLength(0.03, "y");
    gStyle->SetLabelSize(0.05, "xyz");
    
    // Define r values and number of bins
    const int nPoints = 10;
    double rBins[nPoints+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    
    // Create main canvas
    TCanvas *c = new TCanvas("c", "Differential Jet Shapes", 800, 600);
    c->SetMargin(0.12, 0.05, 0.12, 0.10);
    
    // Create histograms
    TH1F *hQuark = new TH1F("hQuark", "", nPoints, rBins);
    TH1F *hGluon = new TH1F("hGluon", "", nPoints, rBins);
    TH1F *hGQ = new TH1F("hGQ", "", nPoints, rBins);
    TH1F *hCombined = new TH1F("hCombined", "", nPoints, rBins);
    TH1F *hResolved = new TH1F("hResolved", "", nPoints, rBins);
    TH1F *hDirect = new TH1F("hDirect", "", nPoints, rBins);
    
    // Fill histograms with stored data
    for (int j = 0; j < nPoints; j++) {
        if (global_diff_shapes_data["QQ_Events"].size() > j) {
            hQuark->SetBinContent(j+1, global_diff_shapes_data["QQ_Events"][j]);
        }
        if (global_diff_shapes_data["GG_Events"].size() > j) {
            hGluon->SetBinContent(j+1, global_diff_shapes_data["GG_Events"][j]);
        }
        if (global_diff_shapes_data["GQ_Events"].size() > j) {
            hGQ->SetBinContent(j+1, global_diff_shapes_data["GQ_Events"][j]);
        }
        if (global_diff_shapes_data["Combined_Events"].size() > j) {
            hCombined->SetBinContent(j+1, global_diff_shapes_data["Combined_Events"][j]);
        }
        if (global_diff_shapes_data["Resolved_Events"].size() > j) {
            hResolved->SetBinContent(j+1, global_diff_shapes_data["Resolved_Events"][j]);
        }
        if (global_diff_shapes_data["Direct_Events"].size() > j) {
            hDirect->SetBinContent(j+1, global_diff_shapes_data["Direct_Events"][j]);
        }
    }
    
    // Style the histograms
    hQuark->SetLineColor(TColor::GetColor("#2A9D8F")); // Teal
    hQuark->SetLineWidth(2);
    
    hGluon->SetLineColor(TColor::GetColor("#E76F51")); // Terracotta
    hGluon->SetLineWidth(2);
    
    hGQ->SetLineColor(TColor::GetColor("#264653"));   // Dark Teal
    hGQ->SetLineWidth(2);
    
    hCombined->SetLineColor(TColor::GetColor("#E9C46A")); // Yellow
    hCombined->SetLineWidth(2);
    
    hResolved->SetLineColor(TColor::GetColor("#F4A261")); // Orange
    hResolved->SetLineWidth(2);
    
    hDirect->SetLineColor(TColor::GetColor("#2A9D8F")); // Same as quark
    hDirect->SetLineWidth(2);
    hDirect->SetLineStyle(2); // Dashed
    
    // Set up the frame
    hQuark->SetMinimum(0);
    hQuark->SetMaximum(6.0);
    hQuark->GetXaxis()->SetRangeUser(0, 1.0);
    hQuark->GetXaxis()->SetTitle("r");
    hQuark->GetYaxis()->SetTitle("#rho");
    hQuark->GetXaxis()->SetTitleFont(132);
    hQuark->GetYaxis()->SetTitleFont(132);
    hQuark->GetXaxis()->SetLabelFont(132);
    hQuark->GetYaxis()->SetLabelFont(132);
    
    // Draw histograms
    hQuark->Draw("HIST");
    hGluon->Draw("HIST SAME");
    hGQ->Draw("HIST SAME");
    
    // Create legend
    TLegend *leg = new TLegend(0.55, 0.60, 0.90, 0.85);
    leg->AddEntry(hQuark, "Quark (QQ)", "l");
    leg->AddEntry(hGluon, "Gluon (GG)", "l");
    leg->AddEntry(hGQ, "Quark+Gluon (GQ)", "l");
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(132);
    leg->SetTextSize(0.04);
    leg->Draw();
    
    // Add analysis info
    TLatex *t = new TLatex();
    t->SetTextFont(132);
    t->SetTextSize(0.04);
    t->SetTextAlign(11);
    
    string analysis_mode = DIJET_ANALYSIS ? "Dijet Analysis" : "Inclusive Analysis";
    t->DrawLatex(0.15, 5.5, analysis_mode.c_str());
    t->DrawLatex(0.15, 5.1, Form("E_{T}^{jet} > %.0f GeV", jet_etMin));
    t->DrawLatex(0.15, 4.7, Form("%.1f < #eta_{jet} < %.1f", etaMin, etaMax));
    
    // Save plots
    string output_name = DIJET_ANALYSIS ? "differential_jet_shapes_dijet" : "differential_jet_shapes_inclusive";
    c->SaveAs((output_name + ".pdf").c_str());
    c->SaveAs((output_name + ".png").c_str());
    
    // Save to ROOT file
    TFile *outputFile = new TFile((output_name + ".root").c_str(), "RECREATE");
    c->Write("canvas");
    hQuark->Write();
    hGluon->Write();
    hGQ->Write();
    hCombined->Write();
    hResolved->Write();
    hDirect->Write();
    outputFile->Close();
    
    cout << "Plots saved as " << output_name << ".pdf/.png/.root\n";
    
    delete c;
}

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================

int main() {
    
    cout << "=============================================================================\n";
    cout << "DIFFERENTIAL JET SHAPES ANALYSIS FOR EIC/HERA PHOTOPRODUCTION EVENTS\n";
    cout << "=============================================================================\n";
    cout << "Analysis Mode:          " << (DIJET_ANALYSIS ? "DIJET (exactly 2 jets)" : "INCLUSIVE (>=1 jet)") << "\n";
    cout << "Analysis Parameters:\n";
    cout << "  Jet Algorithm:      Anti-kT\n";
    cout << "  Jet Radius:         " << R << "\n";
    cout << "  ALL Jets ET:        > " << jet_etMin << " GeV\n";
    cout << "  Jet Eta Range:      [" << etaMin << ", " << etaMax << "]\n";
    cout << "  Particles:          All particles (no pT cut - includes soft)\n";
    cout << "  r Range:            [" << r_min << ", " << r_max << "] with step " << r_step << "\n";
    cout << "  Annulus Width:      ±" << (delta_r/2.0) << " around r\n";
    if (DIJET_ANALYSIS) {
        cout << "  Event Selection:    Exactly 2 jets passing cuts\n";
    } else {
        cout << "  Event Selection:    At least 1 jet passing cuts\n";
    }
    cout << "=============================================================================\n\n";
    
    // Open input file
    TFile* input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file || input_file->IsZombie()) {
        cerr << "ERROR: Cannot open input file: " << input_filename << endl;
        return 1;
    }
    
    cout << "Successfully opened input file: " << input_filename << "\n\n";
    
    // Define event categories to process
    vector<string> categories = {"QQ_Events", "GG_Events", "GQ_Events", 
                                "Combined_Events", "Resolved_Events", "Direct_Events"};
    
    // Storage for differential jet shapes
    map<string, vector<vector<double>>> diff_shapes; // [category][r_index][event_values]
    
    // Initialize storage
    int n_r_points = int((r_max - r_min) / r_step) + 1;
    for (const string& category : categories) {
        diff_shapes[category].resize(n_r_points);
    }
    
    // Process each category
    for (const string& category : categories) {
        
        cout << "Processing " << category << "...\n";
        
        // Get tree from input file
        TDirectory* input_dir = (TDirectory*)input_file->Get(category.c_str());
        if (!input_dir) {
            cout << "  Warning: Directory " << category << " not found, skipping...\n";
            continue;
        }
        
        TTree* input_tree = (TTree*)input_dir->Get(category.c_str());
        if (!input_tree) {
            cout << "  Warning: Tree " << category << " not found, skipping...\n";
            continue;
        }
        
        Long64_t n_entries = input_tree->GetEntries();
        cout << "  Found " << n_entries << " events\n";
        
        if (n_entries == 0) {
            cout << "  No events found, skipping...\n";
            continue;
        }
        
        // Set up input branches
        vector<float> *px = nullptr, *py = nullptr, *pz = nullptr, *energy = nullptr, *eta = nullptr;
        
        if (!input_tree->GetBranch("px") || !input_tree->GetBranch("py") || 
            !input_tree->GetBranch("pz") || !input_tree->GetBranch("energy") || 
            !input_tree->GetBranch("eta")) {
            cout << "  Error: Required branches not found in tree, skipping...\n";
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
        int events_with_selected_jets = 0;
        
        // Process events
        for (Long64_t i = 0; i < n_entries; ++i) {
            
            if (i % report_every == 0) {
                cout << "    Processing event " << i << "/" << n_entries 
                     << " (" << (100*i/n_entries) << "%)\r" << flush;
            }
            
            input_tree->GetEntry(i);
            
            // Check if vectors are valid
            if (!px || !py || !pz || !energy || !eta) continue;
            if (px->size() != py->size() || px->size() != pz->size() || 
                px->size() != energy->size() || px->size() != eta->size()) continue;
            
            events_processed++;
            
            // Reconstruct jets and calculate all differential shapes in one go
            vector<double> shapes = reconstructJetsAndCalculateShapes(*px, *py, *pz, *energy, *eta);
            
            if (shapes.empty()) continue;
            events_with_jets++;
            
            if (DIJET_ANALYSIS) {
                events_with_selected_jets++; // For dijet, this is the same as events_with_jets
            } else {
                events_with_selected_jets++; // For inclusive, this is the same as events_with_jets
            }
            
            // Store the shapes for each r value
            for (int r_idx = 0; r_idx < (int)shapes.size() && r_idx < n_r_points; r_idx++) {
                diff_shapes[category][r_idx].push_back(shapes[r_idx]);
            }
        }
        
        cout << "\n  Events processed: " << events_processed << "\n";
        cout << "  Events with jets: " << events_with_jets 
             << " (" << (100.0*events_with_jets/events_processed) << "%)\n";
        
        if (DIJET_ANALYSIS) {
            cout << "  Events with exactly 2 jets: " << events_with_selected_jets 
                 << " (" << (100.0*events_with_selected_jets/events_processed) << "%)\n\n";
        } else {
            cout << "  Events with >=1 jet: " << events_with_selected_jets 
                 << " (" << (100.0*events_with_selected_jets/events_processed) << "%)\n\n";
        }
    }
    
    // Calculate and print results
    cout << "=============================================================================\n";
    cout << "DIFFERENTIAL JET SHAPE RESULTS\n";
    cout << "=============================================================================\n";
    cout << "Analysis Mode: " << (DIJET_ANALYSIS ? "DIJET" : "INCLUSIVE") << "\n";
    if (DIJET_ANALYSIS) {
        cout << "Requirement: Exactly 2 jets with ET > " << jet_etMin << " GeV, " << etaMin << " < eta < " << etaMax << "\n";
    } else {
        cout << "Requirement: At least 1 jet with ET > " << jet_etMin << " GeV, " << etaMin << " < eta < " << etaMax << "\n";
    }
    cout << "ALL Jets: ET > " << jet_etMin << " GeV, " << etaMin << " < eta < " << etaMax << "\n";
    cout << "Particles: All particles included (no pT cut)\n";
    cout << "=============================================================================\n\n";
    
    // Print header
    cout << setw(6) << "r";
    for (const string& category : categories) {
        cout << setw(15) << category;
    }
    cout << "\n";
    
    // Print separator
    cout << setw(6) << "-----";
    for (const string& category : categories) {
        cout << setw(15) << "-------------";
    }
    cout << "\n";
    
    // Print results for each r value
    for (int r_idx = 0; r_idx < n_r_points; r_idx++) {
        double r = r_min + r_idx * r_step;
        cout << setw(6) << setprecision(1) << fixed << r;
        
        for (const string& category : categories) {
            if (!diff_shapes[category][r_idx].empty()) {
                // Calculate average differential jet shape for this r
                double sum = 0.0;
                for (double val : diff_shapes[category][r_idx]) {
                    sum += val;
                }
                double average = sum / diff_shapes[category][r_idx].size();
                cout << setw(15) << setprecision(4) << fixed << average;
            } else {
                cout << setw(15) << "N/A";
            }
        }
        cout << "\n";
    }
    
    cout << "\n=============================================================================\n";
    cout << "Number of events contributing to each category:\n";
    cout << "=============================================================================\n";
    
    for (const string& category : categories) {
        if (!diff_shapes[category][0].empty()) {
            cout << setw(20) << category << ": " << diff_shapes[category][0].size() << " events\n";
        }
    }
    
    // Store data for plotting and generate C++ arrays
    storeDataForPlotting(diff_shapes);
    
    // Create plots
    createDifferentialJetShapePlot();
    
    cout << "\n=============================================================================\n";
    cout << "ANALYSIS COMPLETED SUCCESSFULLY!\n";
    cout << "Analysis Mode: " << (DIJET_ANALYSIS ? "DIJET (exactly 2 jets)" : "INCLUSIVE (>=1 jet)") << "\n";
    cout << "=============================================================================\n";
    
    input_file->Close();
    return 0;
}