// =============================================================================
// Jet ET Distribution Comparison Script
// =============================================================================
// Description: Plots jet_eta distributions for QQ, GG, GQ, and Combined events
//              from photoproduction analysis on the same canvas for comparison
//
// Author: Analysis script for photoproduction jet studies
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TColor.h"
#include "TPad.h"
#include "TLatex.h"

#include <iostream>
#include <vector>
#include <string>

void plotJetEtaComparison(bool normalize = false) {
    
    // =============================================================================
    // CONFIGURATION
    // =============================================================================
    
    const std::string inputFile = "/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/alljets/eic141_alljet/alljets_eic141_R10_EtMin10.root";
    
    // Subprocess categories to analyze
    std::vector<std::string> categories = {"QQ_Events", "GG_Events", "GQ_Events", "Combined_Events"};
    std::vector<std::string> labels = {"QQ", "GG", "GQ", "Combined"};
    std::vector<int> colors = {kBlue, kRed, kGreen+2, kBlack};
    
    // Histogram settings
    const int nBins = 120;
    const double etaMin = -1.5;
    const double etaMax = 4.0; 
    
    std::cout << "=============================================================================\n";
    std::cout << "JET ET DISTRIBUTION COMPARISON\n";
    // std::cout << "=============================================================================\n";
    // std::cout << "Input file: " << inputFile << "\n";
    // std::cout << "Analyzing categories: ";
    // for (const auto& cat : categories) std::cout << cat << " ";
    // std::cout << "\n=============================================================================\n\n";
    
    // =============================================================================
    // OPEN INPUT FILE
    // =============================================================================
    
    TFile* file = TFile::Open(inputFile.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Cannot open input file: " << inputFile << std::endl;
        return;
    }
    
    std::cout << "Successfully opened input file.\n\n";
    
    // =============================================================================
    // SET ROOT STYLE
    // =============================================================================
    
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(0.8);
    gStyle->SetLineWidth(2);

    gStyle->SetTextFont(42);           // Helvetica
    gStyle->SetLabelFont(42, "XYZ");   // All axes
    gStyle->SetTitleFont(42, "XYZ");   // All titles
    gStyle->SetStatFont(42);
    // Font sizes
    gStyle->SetTextSize(0.04);
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    // Better spacing
    gStyle->SetTitleOffset(1.4, "Y");  // Move Y-title away from axis
    gStyle->SetTitleOffset(1.2, "X");  // Move X-title down
    gStyle->SetLabelOffset(0.01, "XYZ");
    // Tick marks
    gStyle->SetPadTickX(1);         // Tick marks on top
    gStyle->SetPadTickY(1);         // Tick marks on right
    gStyle->SetTickLength(0.02, "XYZ");

    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(kWhite);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetHistFillColor(kWhite);
    
    // =============================================================================
    // CREATE HISTOGRAMS
    // =============================================================================
    
    std::vector<TH1F*> histograms;
    std::vector<int> totalEvents;
    
    for (size_t i = 0; i < categories.size(); ++i) {
        const std::string& category = categories[i];
        
        // std::cout << "Processing " << category << "...\n";
        
        // Get directory and tree
        TDirectory* dir = (TDirectory*)file->Get(category.c_str());
        if (!dir) {
            std::cout << "  Warning: Directory " << category << " not found, skipping...\n";
            histograms.push_back(nullptr);
            totalEvents.push_back(0);
            continue;
        }
        
        TTree* tree = (TTree*)dir->Get(("jets_" + category).c_str());
        if (!tree) {
            std::cout << "  Warning: Tree jets_" << category << " not found, skipping...\n";
            histograms.push_back(nullptr);
            totalEvents.push_back(0);
            continue;
        }
        
        Long64_t nEntries = tree->GetEntries();
        std::cout << "  Found " << nEntries << " events\n";
        
        if (nEntries == 0) {
            std::cout << "  No events found, skipping...\n";
            histograms.push_back(nullptr);
            totalEvents.push_back(0);
            continue;
        }
        
        // Create histogram
        std::string histName = "h_jet_eta_" + category;
        TH1F* hist = new TH1F(histName.c_str(), "", nBins, etaMin, etaMax);
        hist->SetLineColor(colors[i]);
        hist->SetLineWidth(2);
        hist->GetXaxis()->SetTitle("Jet #eta");
        hist->GetYaxis()->SetTitle("Events");
        hist->GetXaxis()->CenterTitle();
        hist->GetYaxis()->CenterTitle();
        hist->GetXaxis()->SetLabelFont(42);
        hist->GetYaxis()->SetLabelFont(42);
        
        // Set up branch reading
        std::vector<float>* jet_eta = nullptr;
        tree->SetBranchAddress("jet_eta", &jet_eta);
        
        int eventsProcessed = 0;
        int totalJets = 0;
        
        // Fill histogram
        for (Long64_t j = 0; j < nEntries; ++j) {
            tree->GetEntry(j);
            
            if (!jet_eta || jet_eta->empty()) continue;
            
            eventsProcessed++;
            
            // Fill histogram with all jets in the event
            for (size_t k = 0; k < jet_eta->size(); ++k) {
                hist->Fill((*jet_eta)[k]);
                totalJets++;
            }
        }
        
        // std::cout << "  Events processed: " << eventsProcessed << "\n";
        // std::cout << "  Total jets filled: " << totalJets << "\n";
        // std::cout << "  Histogram entries: " << hist->GetEntries() << "\n";
        
        // Normalize if requested
        if (normalize && hist->GetEntries() > 0) {
            hist->Scale(1.0 / hist->Integral());
            std::cout << "  Histogram normalized\n";
        }
        // std::cout << "\n";
        
        histograms.push_back(hist);
        totalEvents.push_back(eventsProcessed);
    }
    
    // =============================================================================
    // CREATE CANVAS AND PLOT
    // =============================================================================
    
    TCanvas* canvas = new TCanvas("c_jet_et_comparison", "Jet E_{T} Distribution Comparison", 1200, 1200);
    canvas->SetLeftMargin(0.15);    
    canvas->SetBottomMargin(0.15);  
    canvas->SetTopMargin(0.05);     
    canvas->SetRightMargin(0.05);
    // canvas->SetLogy(1); // Log scale for better visibility
    
    // Find maximum for scaling
    double maxVal = 0;
    for (auto* hist : histograms) {
        if (hist && hist->GetMaximum() > maxVal) {
            maxVal = hist->GetMaximum();
        }
    }
    
    bool first = true;
    TLegend* legend = new TLegend(0.65, 0.65, 0.92, 0.88);
    legend->SetTextSize(0.045);
    legend->SetTextFont(42);
    legend->SetBorderSize(0);        // No border
    legend->SetFillColor(0);         // Transparent
    legend->SetFillStyle(0);         // Transparent
    legend->SetLineWidth(0);         // No line
    // legend->SetHeader("Subprocess", "C"); // Centered header
    
    for (size_t i = 0; i < histograms.size(); ++i) {
        TH1F* hist = histograms[i];
        if (!hist) continue;
        
        if (first) {
            hist->SetMaximum(maxVal * 1.5);
            hist->SetMinimum(0.1); // For log scale
            hist->Draw("HIST");
            first = false;
        } else {
            hist->Draw("HIST SAME");
        }
        
        // Add to legend with event count
        std::string legendEntry = labels[i]; // + " (" + std::to_string(totalEvents[i]) + " events)";
        legend->AddEntry(hist, legendEntry.c_str(), "l");
    }
    
    legend->Draw();
    
    // =============================================================================
    // PRINT STATISTICS
    // =============================================================================
    
    std::cout << "=============================================================================\n";
    std::cout << "COMPARISON STATISTICS\n";
    std::cout << "=============================================================================\n";
    
    for (size_t i = 0; i < categories.size(); ++i) {
        if (histograms[i]) {
            std::cout << labels[i] << ":\n";
            std::cout << "  Events: " << totalEvents[i] << "\n";
            std::cout << "  Total jets: " << (int)histograms[i]->GetEntries() << "\n";
            std::cout << "  Mean Eta: " << histograms[i]->GetMean() << " GeV\n";
            std::cout << "  RMS: " << histograms[i]->GetRMS() << "\n\n";
        }
    }
    
    // =============================================================================
    // SAVE CANVAS
    // =============================================================================
    
    // Extract root file name without path and extension
    std::string rootFileName = inputFile.substr(inputFile.find_last_of("/") + 1);
    rootFileName = rootFileName.substr(0, rootFileName.find_last_of("."));
    
    std::string outputName = "Eta_" + rootFileName + ".pdf";
    canvas->SaveAs(outputName.c_str());
    std::cout << "Canvas saved as: " << outputName << "\n";
    
    // outputName = "Eta_" + rootFileName + ".png";
    // canvas->SaveAs(outputName.c_str());
    // std::cout << "Canvas saved as: " << outputName << "\n";
    

    if (normalize) {
        std::cout << "Histograms are normalized to unit area for shape comparison.\n";
    } else {
        std::cout << "Histograms show absolute event counts.\n";
    }
    
    file->Close();
}

// Main function for standalone execution
int main() {
    plotJetEtaComparison(false); // Default: no normalization
    return 0;
}