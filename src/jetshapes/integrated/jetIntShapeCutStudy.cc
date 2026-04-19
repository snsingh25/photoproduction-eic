// Dijet Jet Shape Analysis with Classification

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream> 
using namespace std;
using namespace fastjet;

const double R = 1.0;
const double JET_SHAPE_RADIUS = 0.3;
const double leading_jet_etMin = 10.0;
const double subleading_jet_etMin = 7.0;
const double etaMin = -1.0;
const double etaMax = 2;
const string input_file = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt7GeV/hera300_pT7/hera300_pT7.root";
// const string input_file = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt5GeV/eic141_pt5/eic141_pT5.root";

// Classification cuts
const double THIN_JET_CUT = 0.8;   // psi > 0.8 for thin jets (quark-like)
const double THICK_JET_CUT = 0.6;  // psi < 0.6 for thick jets (gluon-like)
const double SOFT_JET_CUT = 0.001;  // psi > 0.001 to remove artifacts

// Calculate jet shape (normalized by jet ET)
double calculateJetShape(const vector<PseudoJet>& constituents, const PseudoJet& jet) {
    double energy_within_r = 0.0;
    
    for (const auto& constituent : constituents) {
        double delta_eta = constituent.eta() - jet.eta();
        double delta_phi = constituent.phi() - jet.phi();
        // while (delta_phi > TMath::Pi()) delta_phi -= 2*TMath::Pi();
        // // while (delta_phi < -TMath::Pi()) delta_phi += 2*TMath::Pi();
        double delta_R = sqrt(delta_eta*delta_eta + delta_phi*delta_phi);
        if (delta_R <= JET_SHAPE_RADIUS) energy_within_r += constituent.Et();
    }
    return (jet.Et() > 0.0) ? energy_within_r / jet.Et() : 0.0;
}

int main() {
   
    // Create histograms
    TH1F* h_QQ = new TH1F("h_QQ", "QQ Jets", 50, 0, 1);
    TH1F* h_GG = new TH1F("h_GG", "GG Jets", 50, 0, 1);
    TH1F* h_GQ = new TH1F("h_GQ", "GQ Jets", 50, 0, 1);
    
    h_QQ->SetLineColor(kRed);
    h_QQ->SetFillColor(kRed);
    h_QQ->SetFillStyle(3004);
    h_QQ->SetLineWidth(2);
    
    h_GG->SetLineColor(kBlue);
    h_GG->SetFillColor(kBlue);
    h_GG->SetFillStyle(3005);
    h_GG->SetLineWidth(2);
    
    h_GQ->SetLineColor(kGreen+2);
    h_GQ->SetFillColor(kGreen+2);
    h_GQ->SetFillStyle(3006);
    h_GQ->SetLineWidth(2);
    
    // Classification counters for events (not individual jets)
    int qq_tt_events = 0, qq_tkTk_events = 0, qq_total_events = 0;
    int gg_tt_events = 0, gg_tkTk_events = 0, gg_total_events = 0;
    int gq_tt_events = 0, gq_tkTk_events = 0, gq_total_events = 0;
    int qq_thin_jets = 0, qq_thick_jets = 0, qq_total_clean = 0;
    int gg_thin_jets = 0, gg_thick_jets = 0, gg_total_clean = 0;
    int gq_thin_jets = 0, gq_thick_jets = 0, gq_total_clean = 0;
    int rejected_soft_jets = 0;

    /*
    Event-level counters (dijet events):
        qq_tt_events: QQ events where both jets are thin (Ïˆ > 0.8)
        qq_tkTk_events: QQ events where both jets are thick (Ïˆ < 0.6)
        qq_total_events: Total QQ dijet events analyzed
    Individual jet counters:
        qq_thin_jets: Individual QQ jets classified as thin
        qq_thick_jets: Individual QQ jets classified as thick
        qq_total_clean: Total QQ jets after removing soft jets
    Subprocess categories:
        qq: Quark-quark scattering (processType = 1)
        gg: Gluon-gluon scattering (processType = 2)
        gq: Quark-gluon scattering (processType = 3)

    ___Example: One QQ dijet event with one thin jet and one thick jet contributes___
        0 to qq_tt_events (not both thin)
        0 to qq_tkTk_events (not both thick)
        1 to qq_thin_jets and 1 to qq_thick_jets
    */
    
    // Open input file
    TFile* file = TFile::Open(input_file.c_str());
    if (!file) {
        cout << "ERROR: Cannot open file!" << endl;
        return 1;
    }
    
    TDirectoryFile* dir = (TDirectoryFile*)file->Get("Combined_Events");
    TTree* tree = (TTree*)dir->Get("Combined_Events");
    
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
    int dijet_events = 0;
    int ev = tree->GetEntries();
    
    // Process events
    for (Long64_t i = 0; i < ev; ++i) {
        if (i % 100000 == 0) cout << "Event Processed : " << i << endl;
        
        tree->GetEntry(i);
        total_events++;
        
        if (!px || px->empty()) continue;
        
        // Create particles
        vector<PseudoJet> particles;
        for (size_t j = 0; j < px->size(); ++j) {
            if ((*energy)[j] <= 0) continue;
            particles.push_back(PseudoJet((*px)[j], (*py)[j], (*pz)[j], (*energy)[j]));
        }
        
        if (particles.size() < 2) continue; // Need at least 2 particles to form dijet
        
        // Cluster jets
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);
        vector<PseudoJet> jets = cs.inclusive_jets();
        vector<PseudoJet> sorted_jets_E = sorted_by_E(jets);

        // Select jets
        vector<PseudoJet> selected_jets;
        for (const auto& jet : jets) {
            if (jet.Et() > subleading_jet_etMin && jet.eta() > etaMin && jet.eta() < etaMax) {
                selected_jets.push_back(jet);
            }
        }
        
        if (selected_jets.size() != 2) continue; // Require exactly 2 jets
        if (selected_jets[0].Et() < leading_jet_etMin) continue; // Leading jet cut
        if (selected_jets[1].Et() < subleading_jet_etMin) continue; // Subleading jet cut
        
        dijet_events++;
        
        // Calculate jet shapes and classify
        vector<double> event_psi;
        bool valid_event = true;
        
        for (const auto& jet : selected_jets) {
            vector<PseudoJet> constituents = jet.constituents();
            double psi = calculateJetShape(constituents, jet);
            
            if (psi < 0.0) { // ignore if psi > 1.0, can be an artifact!
                valid_event = false;
                break;
            }
            
            // Remove soft jets (artifacts near 0)
            if (psi < SOFT_JET_CUT) {
                rejected_soft_jets++;
                valid_event = false;
                break;
            }
            
            event_psi.push_back(psi);
            
            // Fill histograms for individual jets
            if (processType == 1) {
                h_QQ->Fill(psi);
                qq_total_clean++;
                if (psi > THIN_JET_CUT) qq_thin_jets++;
                if (psi < THICK_JET_CUT) qq_thick_jets++;
            }
            else if (processType == 2) {
                h_GG->Fill(psi);
                gg_total_clean++;
                if (psi > THIN_JET_CUT) gg_thin_jets++;
                if (psi < THICK_JET_CUT) gg_thick_jets++;
            }
            else if (processType == 3) {
                h_GQ->Fill(psi);
                gq_total_clean++;
                if (psi > THIN_JET_CUT) gq_thin_jets++;
                if (psi < THICK_JET_CUT) gq_thick_jets++;
            }
        }
        
        // Classify events (both jets must be thin for TT, both thick for TkTk)
        if (valid_event && event_psi.size() == 2) {
            bool both_thin = (event_psi[0] > THIN_JET_CUT) && (event_psi[1] > THIN_JET_CUT);
            bool both_thick = (event_psi[0] < THICK_JET_CUT) && (event_psi[1] < THICK_JET_CUT);

            bool thinthick = (event_psi[0] > THIN_JET_CUT) && (event_psi[1] < THICK_JET_CUT);
            bool thickthin = (event_psi[0] < THICK_JET_CUT) && (event_psi[1] > THIN_JET_CUT);
            
            if (processType == 1) {
                qq_total_events++;
                if (both_thin) qq_tt_events++;
                if (both_thick) qq_tkTk_events++;
            }
            else if (processType == 2) {
                gg_total_events++;
                if (both_thin) gg_tt_events++;
                if (both_thick) gg_tkTk_events++;
            }
            else if (processType == 3) {
                gq_total_events++;
                if (both_thin) gq_tt_events++;
                if (both_thick) gq_tkTk_events++;
            }
        }
    }
    
    file->Close();
    
    // Print results
    cout << "\n=============================================================================\n";
    cout << "INTEGRATED JET SHAPE ANALYSIS\n";
    cout << "=============================================================================\n";
    cout << "Total events: " << total_events << endl;
    cout << "Dijet events: " << dijet_events << endl;
    cout << "R : " << R << "; Integrated Jet Radius : " << JET_SHAPE_RADIUS << endl; 
    cout << "Eta range: [" << etaMin << ", " << etaMax << "]" << endl;
    cout << "Leading jet ET > " << leading_jet_etMin << " GeV" << endl;
    cout << "Subleading jet ET > " << subleading_jet_etMin << " GeV" << endl;
    cout << "Rejected soft jets (psi < " << SOFT_JET_CUT << "): " << rejected_soft_jets << endl;
    cout << "\nClean jet statistics (after removing soft jets):" << endl;
    cout << "QQ jets: " << qq_total_clean << " (mean: " << h_QQ->GetMean() << ")" << endl;
    cout << "GG jets: " << gg_total_clean << " (mean: " << h_GG->GetMean() << ")" << endl;
    cout << "GQ jets: " << gq_total_clean << " (mean: " << h_GQ->GetMean() << ")" << endl;
    
    cout << "\n=============================================================================\n";
    cout << "JET CLASSIFICATION RESULTS (TT_QQ_Events/Total_QQ_Events)\n";
    cout << "=============================================================================\n";
    cout << "TT Events (both jets thin, psi > " << THIN_JET_CUT << "):\n";
    cout << "  QQ: " << qq_tt_events << " / " << qq_total_events 
         << " (" << (100.0*qq_tt_events/qq_total_events) << "%)" << endl;
    cout << "  GG: " << gg_tt_events << " / " << gg_total_events 
         << " (" << (100.0*gg_tt_events/gg_total_events) << "%)" << endl;
    cout << "  GQ: " << gq_tt_events << " / " << gq_total_events 
         << " (" << (100.0*gq_tt_events/gq_total_events) << "%)" << endl;
    
    cout << "\nTkTk Events (both jets thick, psi < " << THICK_JET_CUT << "):\n";
    cout << "  QQ: " << qq_tkTk_events << " / " << qq_total_events 
         << " (" << (100.0*qq_tkTk_events/qq_total_events) << "%)" << endl;
    cout << "  GG: " << gg_tkTk_events << " / " << gg_total_events 
         << " (" << (100.0*gg_tkTk_events/gg_total_events) << "%)" << endl;
    cout << "  GQ: " << gq_tkTk_events << " / " << gq_total_events 
         << " (" << (100.0*gq_tkTk_events/gq_total_events) << "%)" << endl;
    
    // Calculate efficiency and purity
    int total_thin = qq_thin_jets + gg_thin_jets + gq_thin_jets;
    int total_thick = qq_thick_jets + gg_thick_jets + gq_thick_jets;
    
    cout << "\n=============================================================================\n";
    cout << "CLASSIFICATION PERFORMANCE (Total Thin or Thick Jets, from all Subprocesses)\n";
    cout << "=============================================================================\n";
    cout << "Total thin jets: " << total_thin << endl;
    cout << "Total thick jets: " << total_thick << endl;
    
    if (total_thin > 0) {
        cout << "\nThin jet purity (QQ_T_Jets/Total_T_Jets)): " << (100.0*qq_thin_jets/total_thin) << "%" << endl;
    }
    if (total_thick > 0) {
        cout << "Thick jet purity (GG_Tk_Jets/Total_Tk_Jets): " << (100.0*gg_thick_jets/total_thick) << "%" << endl;
    }
    
    if (qq_total_clean > 0) {
        cout << "QQ efficiency (thin classification): " << (100.0*qq_thin_jets/qq_total_clean) << "%" << endl;
    }
    if (gg_total_clean > 0) {
        cout << "GG efficiency (thick classification): " << (100.0*gg_thick_jets/gg_total_clean) << "%" << endl;
    }
    
    // NEW: Gluon Efficiency using weighted formula
    cout << "\n=============================================================================\n";
    cout << "GLUON EFFICIENCY (WEIGHTED BY PHYSICS CONTENT)\n";
    cout << "=============================================================================\n";
    cout << "Formula: Gluon Efficiency = (N_GG^thick + 0.5*N_GQ^thick) / (N_QQ^thick + N_GG^thick + 0.5*N_GQ^thick)\n\n";
    
    double numerator_gluon = gg_thick_jets + 0.5 * gq_thick_jets;
    double denominator_thick = qq_thick_jets + gg_thick_jets + 0.5 * gq_thick_jets;
    
    if (denominator_thick > 0.0) {
        double gluon_efficiency = numerator_gluon / denominator_thick;
        cout << "N_GG^thick:    " << gg_thick_jets << " jets\n";
        cout << "N_QQ^thick:    " << qq_thick_jets << " jets\n";
        cout << "N_GQ^thick:    " << gq_thick_jets << " jets\n\n";
        cout << "Numerator:     " << gg_thick_jets << " + 0.5 × " << gq_thick_jets 
             << " = " << numerator_gluon << "\n";
        cout << "Denominator:   " << qq_thick_jets << " + " << gg_thick_jets 
             << " + 0.5 × " << gq_thick_jets << " = " << denominator_thick << "\n\n";
        cout << "Gluon Efficiency = " << fixed << setprecision(4) << gluon_efficiency 
             << " (" << (100.0*gluon_efficiency) << "%)" << endl;
        
        // Complementary metrics
        cout << "\n--- Complementary Metrics (Thick Sample) ---\n";
        double gluon_contamination_in_thick = gg_thick_jets / denominator_thick;
        double quark_contamination_in_thick = qq_thick_jets / denominator_thick;
        double gq_contribution_in_thick = (0.5 * gq_thick_jets) / denominator_thick;
        
        cout << "Gluon Fraction = " << fixed << setprecision(4) 
             << gluon_contamination_in_thick << " (" << (100.0*gluon_contamination_in_thick) << "%)" << endl;
        cout << "Quark Contamination in Thick Sample = " << fixed << setprecision(4) 
             << quark_contamination_in_thick << " (" << (100.0*quark_contamination_in_thick) << "%)" << endl;
        cout << "GQ Mixed Contribution (weighted) = " << fixed << setprecision(4) 
             << gq_contribution_in_thick << " (" << (100.0*gq_contribution_in_thick) << "%)" << endl;
    } else {
        cout << "WARNING: No thick jets found for gluon efficiency calculation!" << endl;
    }
    
    // NEW: Quark Efficiency using weighted formula
    cout << "\n=============================================================================\n";
    cout << "QUARK EFFICIENCY (WEIGHTED BY PHYSICS CONTENT)\n";
    cout << "=============================================================================\n";
    cout << "Formula: Quark Efficiency = (N_QQ^thin + 0.5*N_GQ^thin) / (N_QQ^thin + N_GG^thin + 0.5*N_GQ^thin)\n\n";
    
    double numerator_quark = qq_thin_jets + 0.5 * gq_thin_jets;
    double denominator_thin = qq_thin_jets + gg_thin_jets + 0.5 * gq_thin_jets;
    
    if (denominator_thin > 0.0) {
        double quark_efficiency = numerator_quark / denominator_thin;
        cout << "N_QQ^thin:     " << qq_thin_jets << " jets\n";
        cout << "N_GG^thin:     " << gg_thin_jets << " jets\n";
        cout << "N_GQ^thin:     " << gq_thin_jets << " jets\n\n";
        cout << "Numerator:     " << qq_thin_jets << " + 0.5 × " << gq_thin_jets 
             << " = " << numerator_quark << "\n";
        cout << "Denominator:   " << qq_thin_jets << " + " << gg_thin_jets 
             << " + 0.5 × " << gq_thin_jets << " = " << denominator_thin << "\n\n";
        cout << "Quark Efficiency = " << fixed << setprecision(4) << quark_efficiency 
             << " (" << (100.0*quark_efficiency) << "%)" << endl;
        
        // Complementary metrics
        cout << "\n--- Complementary Metrics (Thin Sample) ---\n";
        double quark_fraction_in_thin = qq_thin_jets / denominator_thin;
        double gluon_contamination_in_thin = gg_thin_jets / denominator_thin;
        double gq_contribution_in_thin = (0.5 * gq_thin_jets) / denominator_thin;
        
        cout << "Quark Fraction = " << fixed << setprecision(4) 
             << quark_fraction_in_thin << " (" << (100.0*quark_fraction_in_thin) << "%)" << endl;
        cout << "Gluon Contamination in Thin Sample = " << fixed << setprecision(4) 
             << gluon_contamination_in_thin << " (" << (100.0*gluon_contamination_in_thin) << "%)" << endl;
        cout << "GQ Mixed Contribution (weighted) = " << fixed << setprecision(4) 
             << gq_contribution_in_thin << " (" << (100.0*gq_contribution_in_thin) << "%)" << endl;
    } else {
        cout << "WARNING: No thin jets found for quark efficiency calculation!" << endl;
    }

    // Write histogram data to text file for Python plotting
    ofstream outfile("int_jet_shape_" + to_string((int)etaMin) + "_" + to_string((int)etaMax) + ".txt");
    outfile << "# Integrated Jet Shape Data" << endl;
    outfile << "# Input: " << input_file << endl;
    outfile << "# R = " << R << endl;
    outfile << "# Integrated radius = " << JET_SHAPE_RADIUS << endl;
    outfile << "# etaMin = " << etaMin << endl;
    outfile << "# etaMax = " << etaMax << endl;
    outfile << "# leading_jet_etMin = " << leading_jet_etMin << endl;
    outfile << "# subleading_jet_etMin = " << subleading_jet_etMin << endl;
    outfile << "# Mean_QQ: " << h_QQ->GetMean() << endl;
    outfile << "# Mean_GG: " << h_GG->GetMean() << endl;
    outfile << "# Mean_GQ: " << h_GQ->GetMean() << endl;
    outfile << "# Entries_QQ: " << qq_total_clean << endl;
    outfile << "# Entries_GG: " << gg_total_clean << endl;
    outfile << "# Entries_GQ: " << gq_total_clean << endl;
    outfile << "# Columns: bin_center, QQ_raw, GG_raw, GQ_raw" << endl;

    for (int bin = 1; bin <= h_QQ->GetNbinsX(); bin++) {
        double bin_center = h_QQ->GetBinCenter(bin);
        double qq_val = h_QQ->GetBinContent(bin);
        double gg_val = h_GG->GetBinContent(bin);
        double gq_val = h_GQ->GetBinContent(bin);
        outfile << bin_center << " " << qq_val << " " << gg_val << " " << gq_val << endl;
    }
    outfile.close();
    // cout << "Data saved to: int_jet_shape_" << to_string((int)etaMin) << "_" << to_string((int)etaMax) << ".txt" << endl;
    cout << "Data saved to: int_jet_shape_" << to_string((int)etaMin) << "_" << to_string((int)etaMax) 
     << "_" << ([&](double v) { std::ostringstream s; s << std::fixed << std::setprecision(1) << v; return s.str(); }(etaMin)) 
     << "_" << ([&](double v) { std::ostringstream s; s << std::fixed << std::setprecision(1) << v; return s.str(); }(etaMax)) 
     << ".txt" << endl;

    // // Make plot
    // gStyle->SetOptStat(0);
    // TCanvas* c1 = new TCanvas("c1", "Jet Shapes (Clean Jets)", 800, 800);
    
    // double max_val = max({h_QQ->GetMaximum(), h_GG->GetMaximum(), h_GQ->GetMaximum()});
    // h_QQ->SetMaximum(max_val * 1.2);
    // h_QQ->SetTitle("Jet Shape #psi(r=0.3); #psi(r=0.3); Entries");
    
    // h_QQ->Draw("HIST");
    // h_GG->Draw("HIST SAME");
    // h_GQ->Draw("HIST SAME");
    
    // TLegend* legend = new TLegend(0.15, 0.70, 0.45, 0.85);
    // legend->AddEntry(h_QQ, Form("QQ (%d)", qq_total_clean), "lf");
    // legend->AddEntry(h_GG, Form("GG (%d)", gg_total_clean), "lf");
    // legend->AddEntry(h_GQ, Form("GQ (%d)", gq_total_clean), "lf");
    // legend->SetBorderSize(1);
    // legend->SetFillColor(kWhite);
    // legend->Draw();
    
    // c1->Print(Form("IntJetShapes0p3_%.1fto%.1f.pdf", etaMin, etaMax));
    // cout << "=============================================================================\n";

    // Cleanup
    delete h_QQ;
    delete h_GG;
    delete h_GQ;
    // delete c1;
    // delete legend;
    
    return 0;
}