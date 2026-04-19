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
// Relative fractions for thin, thick, and combined jet
// const vector<double> thinthin_event_fraction = {0.3014, 0.2663, 0.2493, 0.2049, 0.1696};
// const vector<double> thickthick_event_fraction = {0.6986, 0.7337, 0.7507, 0.7951, 0.8304};
// const vector<double> thickthin_event_fraction = {1.0000, 1.0000, 1.0000, 1.0000, 1.0000};

// Relative fractions for thin, thick, and combined jet
const vector<double> thinthin_event_fraction = {0.0343,0.0420,0.0382,0.0243,0.0081};
const vector<double> thickthick_event_fraction = {0.2858,0.3273,0.3970,0.4866,0.5161};
const vector<double> thickthin_event_fraction = {0.1851,0.1777,0.1648,0.1348,0.0887};

// Colors
const int color_Thin = kRed+1;
const int color_Thick = kGreen+2;
const int color_Combined = kBlue+1;

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
    TCanvas* canvas = new TCanvas("canvas", "", 800, 800);
    gPad->SetMargin(0.12, 0.05, 0.15, 0.10);
    gPad->SetLogy();
    
    // Create graphs
    TGraph* graph_Thin = new TGraph(eta_centers.size());
    TGraph* graph_Thick = new TGraph(eta_centers.size());
    TGraph* graph_ThickThin = new TGraph(eta_centers.size());
    
    for (size_t i = 0; i < eta_centers.size(); ++i) {
        graph_Thin->SetPoint(i, eta_centers[i], thinthin_event_fraction[i]);
        graph_Thick->SetPoint(i, eta_centers[i], thickthick_event_fraction[i]);
        graph_ThickThin->SetPoint(i, eta_centers[i], thickthin_event_fraction[i]);
    }
    
    // Style graphs
    graph_Thin->SetMarkerStyle(20); graph_Thin->SetMarkerColor(color_Thin); 
    graph_Thin->SetLineColor(color_Thin); graph_Thin->SetLineWidth(4);
    
    graph_Thick->SetMarkerStyle(21); graph_Thick->SetMarkerColor(color_Thick); 
    graph_Thick->SetLineColor(color_Thick); graph_Thick->SetLineWidth(4);
    
    graph_ThickThin->SetMarkerStyle(22); graph_ThickThin->SetMarkerColor(color_Combined); 
    graph_ThickThin->SetLineColor(color_Combined); graph_ThickThin->SetLineWidth(4);
    
    // Create multigraph
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(graph_Thin);
    mg->Add(graph_Thick);
    mg->Add(graph_ThickThin);
    
    mg->Draw("ALP");
    mg->GetXaxis()->SetTitle("Pseudorapidity Range (#eta)");
    mg->GetYaxis()->SetTitle("Relative Event Fraction");
    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetRangeUser(-1, 4);
    mg->GetXaxis()->SetNdivisions(6, 0, 0);
    mg->SetMinimum(0.001);
    mg->SetMaximum(0.9);
    
    // Legend
    TLegend* legend = new TLegend(0.65, 0.15, 0.88, 0.35);
    legend->SetTextSize(0.05);
    legend->AddEntry(graph_Thin, "QQ-Events", "lp");
    legend->AddEntry(graph_Thick, "GG-Events", "lp");
    legend->AddEntry(graph_ThickThin, "QG-Events", "lp");
    legend->Draw();
    
    // Save as PDF
    canvas->SaveAs("jet_fractions_thickthin.pdf");
}