// =============================================================================
// Minimalistic Purity Analysis Plots for HERA Dijet Classification
// =============================================================================
// Description: Creates clean line plots showing subprocess composition
//              in thin and thick jet samples across different eta ranges
//
// Usage: root -l purity_plots.cc
// Output: purity_analysis.png, purity_analysis.pdf
// =============================================================================

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGaxis.h"

void purity_plots() {
    
    // Set ROOT style for minimalistic look
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetLabelSize(0.045, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.3, "Y");
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    
    // ==========================================================================
    // DATA FROM HERA ANALYSIS RESULTS
    // ==========================================================================
    
    // X-axis points: center of eta ranges
    double eta_centers[4] = {-0.5, 0.5, 1.5, 2.5};  // Centers of -1to0, 0to1, 1to2, 2to3
    
    // Thin jet composition (percentages)
    double thin_QQ[4] = {64.5, 57.2, 51.2, 44.2};
    double thin_GG[4] = {1.3, 2.5, 4.9, 7.5};
    double thin_GQ[4] = {34.3, 40.3, 44.0, 48.3};
    
    // Thick jet composition (percentages) - COMMENTED OUT FOR NOW
    /*
    double thick_QQ[4] = {42.5, 28.4, 19.0, 14.5};
    double thick_GG[4] = {8.2, 20.5, 33.2, 38.6};
    double thick_GQ[4] = {49.4, 51.2, 47.8, 46.9};
    */
    
    // ==========================================================================
    // CREATE CANVAS AND GRAPHS
    // ==========================================================================
    
    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetFrameFillColor(0);
    
    // Create graphs for thin jets
    TGraph *g_thin_QQ = new TGraph(4, eta_centers, thin_QQ);
    TGraph *g_thin_GG = new TGraph(4, eta_centers, thin_GG);
    TGraph *g_thin_GQ = new TGraph(4, eta_centers, thin_GQ);
    
    // Create graphs for thick jets - COMMENTED OUT
    /*
    TGraph *g_thick_QQ = new TGraph(4, eta_centers, thick_QQ);
    TGraph *g_thick_GG = new TGraph(4, eta_centers, thick_GG);
    TGraph *g_thick_GQ = new TGraph(4, eta_centers, thick_GQ);
    */
    
    // ==========================================================================
    // SET COLORS AND STYLES
    // ==========================================================================
    
    // Define colors
    // Int_t color_QQ = TColor::GetColor("#FF6B6B");  // Red
    // Int_t color_GG = TColor::GetColor("#4ECDC4");  // Teal
    // Int_t color_GQ = TColor::GetColor("#45B7D1");  // Blue
    Int_t color_QQ = kRed+1;      // Red
    Int_t color_GG = kGreen+2;    // Green
    Int_t color_GQ = kBlue+1;     // Blue
    
    // Set graph styles for thin jets
    g_thin_QQ->SetLineColor(color_QQ);
    g_thin_QQ->SetMarkerColor(color_QQ);
    g_thin_QQ->SetLineWidth(4);
    g_thin_QQ->SetMarkerStyle(20);
    g_thin_QQ->SetMarkerSize(1.4);
    
    g_thin_GG->SetLineColor(color_GG);
    g_thin_GG->SetMarkerColor(color_GG);
    g_thin_GG->SetLineWidth(4);
    g_thin_GG->SetMarkerStyle(21);
    g_thin_GG->SetMarkerSize(1.4);
    
    g_thin_GQ->SetLineColor(color_GQ);
    g_thin_GQ->SetMarkerColor(color_GQ);
    g_thin_GQ->SetLineWidth(4);
    g_thin_GQ->SetMarkerStyle(22);
    g_thin_GQ->SetMarkerSize(1.4);
    
    // Set styles for thick jets - COMMENTED OUT
    /*
    g_thick_QQ->SetLineColor(color_QQ);
    g_thick_QQ->SetMarkerColor(color_QQ);
    g_thick_QQ->SetLineWidth(3);
    g_thick_QQ->SetMarkerStyle(24);
    g_thick_QQ->SetMarkerSize(1.2);
    g_thick_QQ->SetLineStyle(2);  // Dashed line
    
    g_thick_GG->SetLineColor(color_GG);
    g_thick_GG->SetMarkerColor(color_GG);
    g_thick_GG->SetLineWidth(3);
    g_thick_GG->SetMarkerStyle(25);
    g_thick_GG->SetMarkerSize(1.2);
    g_thick_GG->SetLineStyle(2);  // Dashed line
    
    g_thick_GQ->SetLineColor(color_GQ);
    g_thick_GQ->SetMarkerColor(color_GQ);
    g_thick_GQ->SetLineWidth(3);
    g_thick_GQ->SetMarkerStyle(26);
    g_thick_GQ->SetMarkerSize(1.2);
    g_thick_GQ->SetLineStyle(2);  // Dashed line
    */
    
    // ==========================================================================
    // CREATE MULTIGRAPH AND DRAW
    // ==========================================================================
    
    TMultiGraph *mg = new TMultiGraph();
    
    // Add thin jet graphs
    mg->Add(g_thin_QQ, "lp");
    mg->Add(g_thin_GG, "lp");
    mg->Add(g_thin_GQ, "lp");
    
    // Add thick jet graphs - COMMENTED OUT
    /*
    mg->Add(g_thick_QQ, "lp");
    mg->Add(g_thick_GG, "lp");
    mg->Add(g_thick_GQ, "lp");
    */
    
    // Draw the multigraph
    mg->Draw("A");
    
    // Set axis properties
    mg->GetXaxis()->SetTitle("Pseudorapidity");
    mg->GetYaxis()->SetTitle("Percentage (%)");
    mg->GetYaxis()->SetTitleOffset(0.9);
    mg->GetXaxis()->SetTitleOffset(0.95);
    mg->GetXaxis()->SetRangeUser(-1.2, 3.2);
    mg->GetYaxis()->SetRangeUser(0, 100);
    
    // Set custom x-axis labels at specific positions
    // mg->GetXaxis()->SetNdivisions(5, 0, 0, kFALSE);  // 5 main divisions, no subdivisions
    // mg->GetXaxis()->CenterTitle(kTRUE);
    // mg->GetYaxis()->CenterTitle(kTRUE);
    // mg->GetXaxis()->SetNdivisions(505, kFALSE);  // 5 major ticks only
    mg->GetXaxis()->CenterTitle(kTRUE);
    mg->GetYaxis()->CenterTitle(kTRUE);
    mg->GetXaxis()->SetLabelOffset(0.01);
    mg->GetYaxis()->SetLabelOffset(0.01);
    
    // Manually set x-axis tick positions
    mg->GetXaxis()->SetLimits(-1.2, 3.2);
    // TGaxis *axis = new TGaxis(-1.2, 0, 3.2, 0, -1, 3, 5, "");
    // axis->SetLabelSize(0.045);
    // axis->Draw();
    
    // ==========================================================================
    // ADD LEGEND (MINIMALISTIC)
    // ==========================================================================
    
    TLegend *legend = new TLegend(0.75, 0.75, 0.93, 0.93);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.045);
    legend->SetHeader("Thin Jets");
    legend->AddEntry(g_thin_QQ, "QQ", "lp");
    legend->AddEntry(g_thin_GG, "GG", "lp");
    legend->AddEntry(g_thin_GQ, "GQ", "lp");
    
    // Add thick jet entries to legend - COMMENTED OUT
    /*
    legend->AddEntry(g_thick_QQ, "QQ (Thick)", "lp");
    legend->AddEntry(g_thick_GG, "GG (Thick)", "lp");
    legend->AddEntry(g_thick_GQ, "GQ (Thick)", "lp");
    */
    
    legend->Draw();
    
    // ==========================================================================
    // SAVE OUTPUTS
    // ==========================================================================
    
    // c1->SaveAs("purity_analysis_thin_jets.png");
    c1->SaveAs("purity_analysis_thin_jets.pdf");
    // c1->SaveAs("purity_analysis_thin_jets.eps");
    
    // Print summary information
    cout << "=============================================" << endl;
    cout << "PURITY ANALYSIS SUMMARY (THIN JETS)" << endl;
    cout << "=============================================" << endl;
    cout << "Eta Range    |  QQ(%)  |  GG(%)  |  GQ(%)  " << endl;
    cout << "-------------|---------|---------|--------- " << endl;
    
    const char* eta_labels[4] = {"-1 to 0", "0 to 1", "1 to 2", "2 to 3"};
    for (int i = 0; i < 4; i++) {
        printf("%-12s | %6.1f  | %6.1f  | %6.1f  \n", 
               eta_labels[i], thin_QQ[i], thin_GG[i], thin_GQ[i]);
    }
    cout << "=============================================" << endl;
    
    // Print thick jet summary - COMMENTED OUT
    /*
    cout << "=============================================" << endl;
    cout << "PURITY ANALYSIS SUMMARY (THICK JETS)" << endl;
    cout << "=============================================" << endl;
    cout << "Eta Range    |  QQ(%)  |  GG(%)  |  GQ(%)  " << endl;
    cout << "-------------|---------|---------|--------- " << endl;
    
    for (int i = 0; i < 4; i++) {
        printf("%-12s | %6.1f  | %6.1f  | %6.1f  \n", 
               eta_labels[i], thick_QQ[i], thick_GG[i], thick_GQ[i]);
    }
    cout << "=============================================" << endl;
    */
    
    cout << "\nPlots saved as:" << endl;
    // cout << "- purity_analysis_thin_jets.png" << endl;
    cout << "- purity_analysis_thin_jets.pdf" << endl;
    // cout << "- purity_analysis_thin_jets.eps" << endl;
    
    // Keep canvas open
    // c1->WaitPrimitive();
}

// Main function for standalone compilation
#ifndef __CINT__
int main() {
    purity_plots_thin();
    return 0;
}
#endif