// =============================================================================
// Purity Evolution Plot for HERA Dijet Classification
// =============================================================================
// Description: Creates clean line plots showing purity evolution
//              QQ purity in thin jets and GG purity in thick jets vs eta
//
// Usage: root -l purity_evolution.cc
// Output: purity_evolution.pdf
// =============================================================================

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGaxis.h"

void purity_evolution() {
    
    // Set ROOT style for minimalistic look
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetLabelSize(0.045, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(0.8, "X");
    gStyle->SetTitleOffset(0.8, "Y");
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    
    // ==========================================================================
    // PURITY DATA FROM HERA ANALYSIS RESULTS
    // ==========================================================================
    
    // X-axis points: center of eta ranges
    double eta_centers[4] = {-0.5, 0.5, 1.5, 2.5};  // Centers of -1to0, 0to1, 1to2, 2to3
    
    // QQ Purity in Thin Jets (%)
    double qq_purity_thin[4] = {64.5, 57.2, 51.2, 44.2};
    
    // GG Purity in Thick Jets (%)
    double gg_purity_thick[4] = {8.2, 20.5, 33.2, 38.6};
    
    // ==========================================================================
    // CREATE CANVAS AND GRAPHS
    // ==========================================================================
    
    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetFrameFillColor(0);
    
    // Create graphs
    TGraph *g_qq_purity = new TGraph(4, eta_centers, qq_purity_thin);
    TGraph *g_gg_purity = new TGraph(4, eta_centers, gg_purity_thick);
    
    // ==========================================================================
    // SET COLORS AND STYLES
    // ==========================================================================
    
    // Define colors
    Int_t color_QQ = kRed+1;      // Red
    Int_t color_GG = kGreen+2;    // Green
    
    // Set graph styles for QQ purity in thin jets
    g_qq_purity->SetLineColor(color_QQ);
    g_qq_purity->SetMarkerColor(color_QQ);
    g_qq_purity->SetLineWidth(4);
    g_qq_purity->SetMarkerStyle(20);  // Circle
    g_qq_purity->SetMarkerSize(1.8);
    
    // Set graph styles for GG purity in thick jets
    g_gg_purity->SetLineColor(color_GG);
    g_gg_purity->SetMarkerColor(color_GG);
    g_gg_purity->SetLineWidth(4);
    g_gg_purity->SetMarkerStyle(21);  // Square
    g_gg_purity->SetMarkerSize(1.8);
    
    // ==========================================================================
    // CREATE MULTIGRAPH AND DRAW
    // ==========================================================================
    
    TMultiGraph *mg = new TMultiGraph();
    
    // Add graphs
    mg->Add(g_qq_purity, "lp");
    mg->Add(g_gg_purity, "lp");
    
    // Draw the multigraph
    mg->Draw("A");
    
    // Set axis properties  
    mg->GetXaxis()->SetTitle("Pseudorapidity");
    mg->GetYaxis()->SetTitle("Purity (%)");
    mg->GetXaxis()->SetLimits(-1, 3);  // Set exact limits
    mg->GetYaxis()->SetRangeUser(0, 80);
    mg->GetXaxis()->SetNdivisions(405, kFALSE);  // 4 intervals = 5 tick marks
    mg->GetXaxis()->CenterTitle(kTRUE);
    mg->GetYaxis()->CenterTitle(kTRUE);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetTitleOffset(0.8);
    mg->GetXaxis()->SetLabelOffset(0.005);
    mg->GetYaxis()->SetLabelOffset(0.005);
    
    // ==========================================================================
    // ADD LEGEND
    // ==========================================================================
    
    TLegend *legend = new TLegend(0.68, 0.72, 0.93, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.045);
    legend->SetHeader("Purity Evolution");
    
    legend->AddEntry(g_qq_purity, "QQ in Thin Jets", "lp");
    legend->AddEntry(g_gg_purity, "GG in Thick Jets", "lp");
    
    legend->Draw();
    
    // ==========================================================================
    // SAVE OUTPUT
    // ==========================================================================
    
    c1->SaveAs("purity_evolution.pdf");
    
    // Print summary information
    cout << "=============================================" << endl;
    cout << "PURITY EVOLUTION SUMMARY" << endl;
    cout << "=============================================" << endl;
    cout << "Eta Range    | QQ in Thin (%) | GG in Thick (%)" << endl;
    cout << "-------------|----------------|----------------" << endl;
    
    const char* eta_labels[4] = {"-1 to 0", "0 to 1", "1 to 2", "2 to 3"};
    for (int i = 0; i < 4; i++) {
        printf("%-12s | %13.1f  | %14.1f  \n", 
               eta_labels[i], qq_purity_thin[i], gg_purity_thick[i]);
    }
    cout << "=============================================" << endl;
    cout << "Physics Trends:" << endl;
    cout << "- QQ purity decreases forward: 64.5% -> 44.2%" << endl;
    cout << "- GG purity increases forward: 8.2% -> 38.6%" << endl;
    cout << "=============================================" << endl;
    
    cout << "\nPlot saved as: purity_evolution.pdf" << endl;
    
    // Keep canvas open
    // c1->WaitPrimitive();
}

// Main function for standalone compilation
#ifndef __CINT__
int main() {
    purity_evolution();
    return 0;
}
#endif