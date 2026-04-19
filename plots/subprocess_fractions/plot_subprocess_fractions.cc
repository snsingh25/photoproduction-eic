// =============================================================================
// Simple ROOT Jet Subprocess Fractions Plotter
// =============================================================================

#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "THStack.h"
#include "TROOT.h"

#include <vector>

using namespace std;

// Data
const vector<double> eta_centers = {-0.5, 0.5, 1.5, 2.5, 3.5};
const vector<double> qq_fractions = {0.4663, 0.3980, 0.3287, 0.2880, 0.2802};
const vector<double> gg_fractions = {0.0742, 0.1263, 0.1767, 0.2100, 0.2101};
const vector<double> gq_fractions = {0.4595, 0.4758, 0.4947, 0.5019, 0.5097};

// Colors
const int color_QQ = kRed+1;
const int color_GG = kGreen+2;
const int color_GQ = kBlue+1;

void setStyle() {
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetLabelSize(0.03, "XYZ");
    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.3, "Y");
    gStyle->SetMarkerSize(2.0);
    gStyle->SetLineWidth(2);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
}

void plot_jet_fractions() {
    
    setStyle();
    
    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 1000);
    gPad->SetMargin(0.12, 0.05, 0.15, 0.10);
    
    // Create graphs
    TGraph* graph_QQ = new TGraph(eta_centers.size());
    TGraph* graph_GG = new TGraph(eta_centers.size());
    TGraph* graph_GQ = new TGraph(eta_centers.size());
    
    for (size_t i = 0; i < eta_centers.size(); ++i) {
        graph_QQ->SetPoint(i, eta_centers[i], qq_fractions[i]);
        graph_GG->SetPoint(i, eta_centers[i], gg_fractions[i]);
        graph_GQ->SetPoint(i, eta_centers[i], gq_fractions[i]);
    }
    
    // Style graphs
    graph_QQ->SetMarkerStyle(20); graph_QQ->SetMarkerColor(color_QQ); 
    graph_QQ->SetLineColor(color_QQ); graph_QQ->SetLineWidth(4);
    
    graph_GG->SetMarkerStyle(21); graph_GG->SetMarkerColor(color_GG); 
    graph_GG->SetLineColor(color_GG); graph_GG->SetLineWidth(4);
    
    graph_GQ->SetMarkerStyle(22); graph_GQ->SetMarkerColor(color_GQ); 
    graph_GQ->SetLineColor(color_GQ); graph_GQ->SetLineWidth(4);
    
    // Create multigraph
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(graph_QQ);
    mg->Add(graph_GG);
    mg->Add(graph_GQ);
    
    mg->Draw("ALP");
    mg->GetXaxis()->SetTitle("Pseudorapidity Range (#eta)");
    mg->GetYaxis()->SetTitle("Relative Fraction");
    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);

    mg->GetXaxis()->SetLimits(-1.0, 4.0);
    mg->GetXaxis()->SetRangeUser(-1.0, 4.0);
    mg->GetXaxis()->SetNdivisions(505, kFALSE);
    mg->GetYaxis()->SetRangeUser(0.0, 0.60);
    mg->GetYaxis()->SetNdivisions(503, kFALSE);

    mg->GetXaxis()->SetTickLength(0.02);  // Default is 0.03
    mg->GetYaxis()->SetTickLength(0.02);
    gPad->SetTickx(1);  // Ticks on top and bottom
    gPad->SetTicky(1);  // Ticks on left and right
    mg->GetXaxis()->SetLabelSize(0.03);   // Default is 0.03
    mg->GetYaxis()->SetLabelSize(0.03);

    // // Legend
    TLegend* legend = new TLegend(0.75, 0.20, 0.85, 0.35);
    legend->SetTextSize(0.03);
    legend->AddEntry(graph_QQ, "QQ", "lp");
    legend->AddEntry(graph_GG, "GG", "lp");
    legend->AddEntry(graph_GQ, "GQ", "lp");
    legend->Draw();
    
    // Save as PDF
    canvas->SaveAs("jet_subprocess_fractions.pdf");
}