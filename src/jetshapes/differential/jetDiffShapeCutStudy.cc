// Differential Jet Shape Analysis with Classification for r=0.3

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
using namespace std;
using namespace fastjet;

const double R = 1.0;
const double DIFF_SHAPE_RADIUS = 0.2;     // Center radius for differential shape
const double DELTA_R = 0.01;              // Annulus width (Δr)
const double leading_jet_etMin = 10.0;
const double subleading_jet_etMin = 7.0;
const double etaMin = -1.0;
const double etaMax = 4.0;
const string input_file = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt7GeV/hera300_pT7/hera300_pT7.root";

// Classification cuts based on differential jet shape Rho(r=0.3)
const double THIN_JET_CUT = 2.5;    // Rho > 2.0 for thin jets (quark-like)
const double THICK_JET_CUT = 2.0;   // Rho < 2.0 for thick jets (gluon-like)
const double SOFT_JET_CUT = 0.001;  // Rho > 0.001 to remove artifacts

// Calculate differential jet shape using the formula:
// Rho(r) = (1/N_jets) * (1/Δr) * Σ[E_T(r-Δr/2, r+Δr/2) / E_T(0, r=1)]
double calculateDifferentialJetShape(const vector<PseudoJet>& constituents, const PseudoJet& jet) {
    double R1 = DIFF_SHAPE_RADIUS - (DELTA_R / 2.0);  // Inner radius
    double R2 = DIFF_SHAPE_RADIUS + (DELTA_R / 2.0);  // Outer radius
    
    double energy_in_annulus = 0.0;
    double total_energy = 0.0;
    
    // Calculate E_T in annulus [R1, R2] and total E_T up to r=1
    for (const auto& constituent : constituents) {
        double delta_eta = constituent.eta() - jet.eta();
        double delta_phi = constituent.phi() - jet.phi();
        
        double delta_R = sqrt(delta_eta*delta_eta + delta_phi*delta_phi);
        
        // Energy in annulus [R1, R2]
        if (delta_R >= R1 && delta_R <= R2) {
            energy_in_annulus += constituent.Et();
        }
        
        // Total energy within r=1 for normalization
        if (delta_R <= 1.0) {
            total_energy += constituent.Et();
        }
    }
    
    if (total_energy > 0.0) {
        return (energy_in_annulus) / total_energy;
    }
    return 0.0;
}

int main() {
   
    // Create histograms for differential jet shapes
    TH1F* h_QQ = new TH1F("h_QQ", "QQ Jets - Differential Shape", 50, 0, 1);
    TH1F* h_GG = new TH1F("h_GG", "GG Jets - Differential Shape", 50, 0, 1);
    TH1F* h_GQ = new TH1F("h_GQ", "GQ Jets - Differential Shape", 50, 0, 1);
    
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
        qq_tt_events: QQ events where both jets are thin (Rho > 2.0)
        qq_tkTk_events: QQ events where both jets are thick (Rho < 2.0)
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
        
        // Sort jets by ET
        sort(jets.begin(), jets.end(), [](const PseudoJet &a, const PseudoJet &b) {
            return a.Et() > b.Et();
        });
        
        // Apply jet cuts and require exactly 2 jets for dijet analysis
        vector<PseudoJet> selected_jets;
        for (const auto& jet : jets) {
            if (jet.Et() > subleading_jet_etMin && 
                jet.eta() >= etaMin && jet.eta() <= etaMax) {
                selected_jets.push_back(jet);
            }
        }
        
        // Require dijet topology
        if (selected_jets.size() < 2) continue;
        if (selected_jets[0].Et() < leading_jet_etMin) continue;
        if (selected_jets[1].Et() < subleading_jet_etMin) continue;
        
        // Keep only the two leading jets
        selected_jets.resize(2);
        dijet_events++;
        
        // Calculate differential jet shapes for both jets
        vector<double> jet_rho_values;
        bool both_jets_valid = true;
        
        for (const auto& jet : selected_jets) {
            vector<PseudoJet> constituents = jet.constituents();
            double rho = calculateDifferentialJetShape(constituents, jet);
            
            // Apply soft jet rejection
            if (rho < SOFT_JET_CUT) {
                rejected_soft_jets++;
                both_jets_valid = false;
                break;
            }
            
            jet_rho_values.push_back(rho);
        }
        
        if (!both_jets_valid) continue;
        
        // Fill histograms and classify jets
        vector<bool> thin_jets, thick_jets;
        
        for (double rho : jet_rho_values) {
            bool is_thin = (rho > THIN_JET_CUT);
            bool is_thick = (rho < THICK_JET_CUT);
            
            thin_jets.push_back(is_thin);
            thick_jets.push_back(is_thick);
            
            // Fill histograms based on process type
            if (processType == 1) {  // QQ
                h_QQ->Fill(rho);
                qq_total_clean++;
                if (is_thin) qq_thin_jets++;
                if (is_thick) qq_thick_jets++;
            }
            else if (processType == 2) {  // GG
                h_GG->Fill(rho);
                gg_total_clean++;
                if (is_thin) gg_thin_jets++;
                if (is_thick) gg_thick_jets++;
            }
            else if (processType == 3) {  // GQ
                h_GQ->Fill(rho);
                gq_total_clean++;
                if (is_thin) gq_thin_jets++;
                if (is_thick) gq_thick_jets++;
            }
        }
        
        // Event-level classification (both jets must satisfy condition)
        bool both_thin = (thin_jets[0] && thin_jets[1]);
        bool both_thick = (thick_jets[0] && thick_jets[1]);
        
        // Count events by process type
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
    
    file->Close();
    
    // Print results
    cout << "\n=============================================================================\n";
    cout << "DIFFERENTIAL JET SHAPE ANALYSIS (r = " << DIFF_SHAPE_RADIUS << ")\n";
    cout << "=============================================================================\n";
    cout << "Total events: " << total_events << endl;
    cout << "Dijet events: " << dijet_events << endl;
    cout << "R : " << R << "; Differential Shape Radius : " << DIFF_SHAPE_RADIUS << endl; 
    cout << "Annulus width (Δr): " << DELTA_R << endl;
    cout << "Eta range: [" << etaMin << ", " << etaMax << "]" << endl;
    cout << "Leading jet ET > " << leading_jet_etMin << " GeV" << endl;
    cout << "Subleading jet ET > " << subleading_jet_etMin << " GeV" << endl;
    cout << "Rejected soft jets (Rho < " << SOFT_JET_CUT << "): " << rejected_soft_jets << endl;
    cout << "\nClean jet statistics (after removing soft jets):" << endl;
    cout << "QQ jets: " << qq_total_clean << " (mean: " << h_QQ->GetMean() << ")" << endl;
    cout << "GG jets: " << gg_total_clean << " (mean: " << h_GG->GetMean() << ")" << endl;
    cout << "GQ jets: " << gq_total_clean << " (mean: " << h_GQ->GetMean() << ")" << endl;
    
    cout << "\n=============================================================================\n";
    cout << "DIFFERENTIAL JET SHAPE CLASSIFICATION RESULTS\n";
    cout << "=============================================================================\n";
    cout << "Classification cuts:" << endl;
    cout << "  Thin jets (quark-like):  Rho > " << THIN_JET_CUT << endl;
    cout << "  Thick jets (gluon-like): Rho < " << THICK_JET_CUT << endl;
    
    cout << "\nIndividual Jet Classification:" << endl;
    cout << "QQ subprocess:" << endl;
    cout << "  Thin jets: " << qq_thin_jets << " / " << qq_total_clean 
         << " (" << (qq_total_clean > 0 ? 100.0*qq_thin_jets/qq_total_clean : 0) << "%)" << endl;
    cout << "  Thick jets: " << qq_thick_jets << " / " << qq_total_clean 
         << " (" << (qq_total_clean > 0 ? 100.0*qq_thick_jets/qq_total_clean : 0) << "%)" << endl;
    
    cout << "GG subprocess:" << endl;
    cout << "  Thin jets: " << gg_thin_jets << " / " << gg_total_clean 
         << " (" << (gg_total_clean > 0 ? 100.0*gg_thin_jets/gg_total_clean : 0) << "%)" << endl;
    cout << "  Thick jets: " << gg_thick_jets << " / " << gg_total_clean 
         << " (" << (gg_total_clean > 0 ? 100.0*gg_thick_jets/gg_total_clean : 0) << "%)" << endl;
    
    cout << "GQ subprocess:" << endl;
    cout << "  Thin jets: " << gq_thin_jets << " / " << gq_total_clean 
         << " (" << (gq_total_clean > 0 ? 100.0*gq_thin_jets/gq_total_clean : 0) << "%)" << endl;
    cout << "  Thick jets: " << gq_thick_jets << " / " << gq_total_clean 
         << " (" << (gq_total_clean > 0 ? 100.0*gq_thick_jets/gq_total_clean : 0) << "%)" << endl;
    
    cout << "\nDijet Event Classification (both jets must satisfy condition):" << endl;
    cout << "QQ events: " << qq_total_events << endl;
    cout << "  Both thin: " << qq_tt_events << " (" << (qq_total_events > 0 ? 100.0*qq_tt_events/qq_total_events : 0) << "%)" << endl;
    cout << "  Both thick: " << qq_tkTk_events << " (" << (qq_total_events > 0 ? 100.0*qq_tkTk_events/qq_total_events : 0) << "%)" << endl;
    
    cout << "GG events: " << gg_total_events << endl;
    cout << "  Both thin: " << gg_tt_events << " (" << (gg_total_events > 0 ? 100.0*gg_tt_events/gg_total_events : 0) << "%)" << endl;
    cout << "  Both thick: " << gg_tkTk_events << " (" << (gg_total_events > 0 ? 100.0*gg_tkTk_events/gg_total_events : 0) << "%)" << endl;
    
    cout << "GQ events: " << gq_total_events << endl;
    cout << "  Both thin: " << gq_tt_events << " (" << (gq_total_events > 0 ? 100.0*gq_tt_events/gq_total_events : 0) << "%)" << endl;
    cout << "  Both thick: " << gq_tkTk_events << " (" << (gq_total_events > 0 ? 100.0*gq_tkTk_events/gq_total_events : 0) << "%)" << endl;
    
    cout << "\n=============================================================================\n";
    cout << "DISCRIMINATION POWER ANALYSIS\n";
    cout << "=============================================================================\n";
    
    // Calculate total thin/thick fractions
    int total_thin = qq_thin_jets + gg_thin_jets + gq_thin_jets;
    int total_thick = qq_thick_jets + gg_thick_jets + gq_thick_jets;
    int total_clean = qq_total_clean + gg_total_clean + gq_total_clean;
    
    cout << "Overall jet classification:" << endl;
    cout << "  Total clean jets: " << total_clean << endl;
    cout << "  Thin jets: " << total_thin << " (" << (total_clean > 0 ? 100.0*total_thin/total_clean : 0) << "%)" << endl;
    cout << "  Thick jets: " << total_thick << " (" << (total_clean > 0 ? 100.0*total_thick/total_clean : 0) << "%)" << endl;
    
    // Calculate purity and efficiency for quark vs gluon discrimination
    double quark_purity_thin = (total_thin > 0) ? 100.0 * qq_thin_jets / total_thin : 0.0;
    double gluon_purity_thick = (total_thick > 0) ? 100.0 * gg_thick_jets / total_thick : 0.0;
    double quark_efficiency = (qq_total_clean > 0) ? 100.0 * qq_thin_jets / qq_total_clean : 0.0;
    double gluon_efficiency = (gg_total_clean > 0) ? 100.0 * gg_thick_jets / gg_total_clean : 0.0;
    
    cout << "\nQuark-Gluon discrimination:" << endl;
    cout << "  Quark purity in thin jets: " << quark_purity_thin << "%" << endl;
    cout << "  Gluon purity in thick jets: " << gluon_purity_thick << "%" << endl;
    cout << "  Quark tagging efficiency: " << quark_efficiency << "%" << endl;
    cout << "  Gluon tagging efficiency: " << gluon_efficiency << "%" << endl;
    
    // Create and save plots
    TCanvas* c1 = new TCanvas("c1", "Differential Jet Shape Distributions", 800, 600);
    
    h_GG->Draw("HIST");
    h_QQ->Draw("HIST SAME");
    h_GQ->Draw("HIST SAME");
    
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h_QQ, "QQ (quark-quark)", "f");
    legend->AddEntry(h_GG, "GG (gluon-gluon)", "f");
    legend->AddEntry(h_GQ, "GQ (quark-gluon)", "f");
    legend->Draw();
    
    c1->SetTitle("Differential Jet Shape Rho(r=0.3)");
    h_GG->GetXaxis()->SetTitle("Rho(r=0.3)");
    h_GG->GetYaxis()->SetTitle("Entries");
    double max_val = max({h_QQ->GetMaximum(), h_GG->GetMaximum(), h_GQ->GetMaximum()});
    h_GG->SetMaximum(max_val * 1.1); 

    c1->Print(Form("DiffJetShapes0p3_%.1fto%.1f.pdf", etaMin, etaMax));
    // c1->SaveAs("Try1.pdf");
    
    cout << "=============================================================================\n";
    
    // Clean up
    delete c1;
    delete legend;
    delete h_QQ;
    delete h_GG;
    delete h_GQ;
    
    return 0;
}