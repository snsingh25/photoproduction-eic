#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPad.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TColor.h"
#include <vector>
#include <iostream>
#include <iomanip>

void plot_combined_et_vs_eta() {
    
    // Enable LaTeX fonts and set publication quality style
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    // LaTeX font settings
    gStyle->SetTextFont(42);           // Helvetica (closest to Computer Modern in ROOT)
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetStatFont(42);
    gStyle->SetTextSize(0.04);
    gStyle->SetLabelSize(0.05, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    
    // Use ROOT's default jet colormap (similar to matplotlib's jet)
    gStyle->SetPalette(1);  // Rainbow palette
    gStyle->SetNumberContours(100);
    
    // Clean appearance
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);
    
    // Dataset configuration
    struct DatasetInfo {
        const char* filepath;
        const char* label;
        const char* shortname;
    };
    
    DatasetInfo datasets[4] = {
        {"/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/hera300_dijets_pT7/dijets_hera300_R10_EtMin0_10_7.root", 
         "ep, 300 GeV", "hera300"},
        {"/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic141_dijets_pT5/dijets_eic141_R10_EtMin0_10_7.root", 
         "ep, 141 GeV", "eic141"},
        {"/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic105_dijets_pT5/dijets_eic105_R10_EtMin0_10_7.root", 
         "ep, 105 GeV", "eic105"},
        {"/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic64_dijets_pT5/dijets_eic64_R10_EtMin0_10_7.root", 
         "ep, 64 GeV", "eic64"}
    };
    
    const char* eventTypes[2] = {"QQ_Events", "GG_Events"};
    const char* eventLabels[2] = {"QQ", "GG"};
    
    // Create canvas with precise dimensions
    TCanvas* canvas = new TCanvas("canvas", "Jet E_{T} vs #eta Analysis", 1200, 1200);
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetFrameFillColor(0);
    canvas->SetFrameBorderMode(0);
    
    // Create 4x3 grid - left plots, middle colorbars, right plots
    TPad* pads[8];  // Only 8 plotting pads (skip middle colorbar pads)
    const double plotWidth = 0.32;    // Width for each plot
    const double cbarWidth = 0.001;    // Width for colorbar space
    const double plotHeight = 0.22;
    const double leftMargin = 0.08;
    const double topMargin = 0.04;
    const double vSpacing = 0.0;
    
    
    // Create only the 8 plotting pads (4 rows x 2 columns)
    int padCount = 0;
    for (int row = 0; row < 4; row++) {
        for (int plotCol = 0; plotCol < 2; plotCol++) {  // Only 2 plotting columns
            
            double x1, x2;
            if (plotCol == 0) {        // Left plot
                x1 = leftMargin;
                x2 = x1 + plotWidth;
            } else {                   // Right plot 
                x1 = leftMargin + plotWidth + cbarWidth;
                x2 = x1 + plotWidth;   // Same width as left
            }
            
            double y1 = 1.0 - topMargin - (row + 1) * plotHeight - row * vSpacing;
            double y2 = y1 + plotHeight;
            
            TString padName = Form("pad_%d", padCount);
            pads[padCount] = new TPad(padName, "", x1, y1, x2, y2);
            pads[padCount]->SetFillColor(0);
            pads[padCount]->SetBorderMode(0);
            pads[padCount]->SetFrameBorderMode(0);
            
            // Set margins for individual colorbars
            pads[padCount]->SetLeftMargin(plotCol == 0 ? 0.06 : 0.06);   
            pads[padCount]->SetRightMargin(0.15);  // Space for individual colorbars
            pads[padCount]->SetBottomMargin(row == 3 ? 0.12 : 0.01); 
            pads[padCount]->SetTopMargin(0.01);
            
            // Enable ticks and linear scale
            pads[padCount]->SetTickx(1);
            pads[padCount]->SetTicky(1);
            pads[padCount]->SetLogz(0);  // Linear scale for colors
            
            canvas->cd();
            pads[padCount]->Draw();
            padCount++;
        }
    }
    
    // Storage for histograms
    TH2F* histograms[8];
    for (int i = 0; i < 8; i++) {
        histograms[i] = nullptr;
    }
    
    // Statistics storage
    struct Statistics {
        int n_jets;
        double mean_eta;
        double mean_et;
        double std_eta;
        double std_et;
    };
    Statistics stats[8];
    
    // Process each dataset
    for (int fileIdx = 0; fileIdx < 4; fileIdx++) {
        
        std::cout << "Processing " << datasets[fileIdx].label << "..." << std::endl;
        
        TFile* file = TFile::Open(datasets[fileIdx].filepath, "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Cannot open file " << datasets[fileIdx].filepath << std::endl;
            continue;
        }
        
        for (int eventIdx = 0; eventIdx < 2; eventIdx++) {
            
            int histIdx = fileIdx * 2 + eventIdx;
            
            // Navigate to event directory and get tree
            TDirectory* dir = (TDirectory*)file->Get(eventTypes[eventIdx]);
            if (!dir) {
                std::cerr << "Error: Cannot find directory " << eventTypes[eventIdx] << std::endl;
                histograms[histIdx] = nullptr;
                continue;
            }
            
            TString treeName = TString("jets_") + eventTypes[eventIdx];
            TTree* tree = (TTree*)dir->Get(treeName);
            if (!tree) {
                std::cerr << "Error: Cannot find tree " << treeName << std::endl;
                histograms[histIdx] = nullptr;
                continue;
            }
            
            TString histName = Form("h_%s_%s", datasets[fileIdx].shortname, eventLabels[eventIdx]);
            TString histTitle = Form("%s %s;#eta_{jet};E_{T} [GeV]", 
                                   datasets[fileIdx].label, eventLabels[eventIdx]);
            
            histograms[histIdx] = new TH2F(histName, histTitle,
                                         60, -1.5, 4.0,   // eta range (30 bins)
                                         60, 5.0, 70.0);  // ET range (30 bins)
            histograms[histIdx]->SetDirectory(0);
            
            // Set up branches
            std::vector<float> *jet_et = nullptr;
            std::vector<float> *jet_eta = nullptr;
            
            if (!tree->GetBranch("jet_et") || !tree->GetBranch("jet_eta")) {
                std::cerr << "Error: Required branches not found" << std::endl;
                delete histograms[histIdx];
                histograms[histIdx] = nullptr;
                continue;
            }
            
            tree->SetBranchAddress("jet_et", &jet_et);
            tree->SetBranchAddress("jet_eta", &jet_eta);
            
            Long64_t nEntries = tree->GetEntries();
            std::cout << "  Processing " << nEntries << " events for " << eventLabels[eventIdx] << std::endl;
            
            // Fill histogram and calculate statistics
            int totalJets = 0;
            double sum_eta = 0.0, sum_et = 0.0;
            double sum_eta2 = 0.0, sum_et2 = 0.0;
            
            for (Long64_t i = 0; i < nEntries; i++) {
                
                if (i % 100000 == 0 && i > 0) {
                    std::cout << "    Event " << i << "/" << nEntries 
                             << " (" << (100.0 * i / nEntries) << "%)" << std::endl;
                }
                
                tree->GetEntry(i);
                
                if (!jet_et || !jet_eta || jet_et->size() != jet_eta->size()) continue;
                
                for (size_t j = 0; j < jet_et->size(); j++) {
                    if (std::isfinite((*jet_eta)[j]) && std::isfinite((*jet_et)[j])) {
                        histograms[histIdx]->Fill((*jet_eta)[j], (*jet_et)[j]);
                        
                        // Calculate statistics
                        sum_eta += (*jet_eta)[j];
                        sum_et += (*jet_et)[j];
                        sum_eta2 += (*jet_eta)[j] * (*jet_eta)[j];
                        sum_et2 += (*jet_et)[j] * (*jet_et)[j];
                        totalJets++;
                    }
                }
            }
            
            // Store statistics
            stats[histIdx].n_jets = totalJets;
            if (totalJets > 0) {
                stats[histIdx].mean_eta = sum_eta / totalJets;
                stats[histIdx].mean_et = sum_et / totalJets;
                stats[histIdx].std_eta = sqrt(sum_eta2 / totalJets - stats[histIdx].mean_eta * stats[histIdx].mean_eta);
                stats[histIdx].std_et = sqrt(sum_et2 / totalJets - stats[histIdx].mean_et * stats[histIdx].mean_et);
            }
            
            std::cout << "  Total jets: " << totalJets << std::endl;
        }
        
        file->Close();
        delete file;
    }
    
    // Draw histograms
    for (int histIdx = 0; histIdx < 8; histIdx++) {
        if (!histograms[histIdx]) continue;
        
        pads[histIdx]->cd();
        
        // Set axis ranges
        histograms[histIdx]->GetXaxis()->SetRangeUser(-1.5, 4.0);
        histograms[histIdx]->GetYaxis()->SetRangeUser(5.0, 70.0);
        
        // Set individual color ranges for each pad
        if (histIdx == 0) {        // HERA 300 QQ
            histograms[histIdx]->SetMinimum(0);
            histograms[histIdx]->SetMaximum(1200);
        } else if (histIdx == 1) { // HERA 300 GG
            histograms[histIdx]->SetMinimum(0);
            histograms[histIdx]->SetMaximum(900);
        } else if (histIdx == 2) { // EIC 141 QQ
            histograms[histIdx]->SetMinimum(0);
            histograms[histIdx]->SetMaximum(500);
        } else if (histIdx == 3) { // EIC 141 GG
            histograms[histIdx]->SetMinimum(0);
            histograms[histIdx]->SetMaximum(300);
        } else if (histIdx == 4) { // EIC 105 QQ
            histograms[histIdx]->SetMinimum(0);
            histograms[histIdx]->SetMaximum(500);
        } else if (histIdx == 5) { // EIC 105 GG
            histograms[histIdx]->SetMinimum(0);
            histograms[histIdx]->SetMaximum(180);
        } else if (histIdx == 6) { // EIC 64 QQ
            histograms[histIdx]->SetMinimum(0);
            histograms[histIdx]->SetMaximum(500);
        } else if (histIdx == 7) { // EIC 64 GG
            histograms[histIdx]->SetMinimum(0);
            histograms[histIdx]->SetMaximum(80);
        }
        
        // Configure axis position
        int row = histIdx / 2;
        int col = histIdx % 2;
        
        // Style the axes
        histograms[histIdx]->GetXaxis()->SetTitleFont(42);
        histograms[histIdx]->GetXaxis()->SetLabelFont(42);
        histograms[histIdx]->GetYaxis()->SetTitleFont(42);
        histograms[histIdx]->GetYaxis()->SetLabelFont(42);
        
        // Remove individual axis titles
        histograms[histIdx]->GetXaxis()->SetTitle("");
        histograms[histIdx]->GetYaxis()->SetTitle("");
        
        // Show axis values and increase tick sizes
        histograms[histIdx]->GetXaxis()->SetLabelSize(0.08);
        // histograms[histIdx]->GetYaxis()->SetLabelSize(0.08);
        histograms[histIdx]->GetXaxis()->SetTickLength(0.05);
        histograms[histIdx]->GetYaxis()->SetTickLength(0.03);
        
        // Show Y-axis values only on left column
        if (col == 0) {
            histograms[histIdx]->GetYaxis()->SetLabelSize(0.08);
            histograms[histIdx]->GetYaxis()->SetNdivisions(506);  // 5 major divisions, 6 minor
            // This will show: 10, 20, 30, 40, 50, 60 (skipping 70)
        } else {
            histograms[histIdx]->GetYaxis()->SetLabelSize(0.08);
            histograms[histIdx]->GetYaxis()->SetNdivisions(506);
        }
        
        // Show X-axis values only on bottom row
        if (row == 3) {
            histograms[histIdx]->GetXaxis()->SetLabelSize(0.08);
        } else {
            histograms[histIdx]->GetXaxis()->SetLabelSize(0);
        }
        
        // Colorbar values
        double zmin = histograms[histIdx]->GetMinimum();
        double zmax = histograms[histIdx]->GetMaximum();
        histograms[histIdx]->GetZaxis()->SetRangeUser(zmin + 1, zmax - 1);
        
        // Draw with individual colorbars
        histograms[histIdx]->Draw("COLZ");
        
        // Increase colorbar font size
        TPaletteAxis *palette = (TPaletteAxis*)histograms[histIdx]->GetListOfFunctions()->FindObject("palette");
        if (palette) {
            palette->SetLabelSize(0.06);
        }

        // Add dataset and process labels at top right
        TLatex* topRightLabel = new TLatex();
        topRightLabel->SetNDC();
        topRightLabel->SetTextFont(42);
        topRightLabel->SetTextSize(0.08);
        topRightLabel->SetTextColor(kBlack);
        topRightLabel->SetTextAlign(13);  // Right top alignment
        TLatex* topRightLabel_event = new TLatex();
        topRightLabel_event->SetNDC();
        topRightLabel_event->SetTextFont(42);
        topRightLabel_event->SetTextSize(0.08);
        topRightLabel_event->SetTextColor(kBlack);
        topRightLabel_event->SetTextAlign(13);  // Right top alignment
        // Create the label text
        TString labelText = Form("%s", datasets[histIdx/2].label);
        TString labelText_event = Form("%s", eventLabels[histIdx%2]);
        // topRightLabel->DrawLatex(3.8, 58, labelText);  
        topRightLabel->DrawLatex(0.12, 0.90, labelText);
        topRightLabel_event->DrawLatex(0.12, 0.80, labelText_event);
        
        // Update pad
        pads[histIdx]->Modified();
        pads[histIdx]->Update();
    }
    
    // Add common axis titles
    canvas->cd();
    
    TLatex* commonYTitle = new TLatex();
    commonYTitle->SetNDC();
    commonYTitle->SetTextFont(42);
    commonYTitle->SetTextSize(0.035);
    commonYTitle->SetTextAlign(22);
    commonYTitle->SetTextAngle(90);
    commonYTitle->DrawLatex(0.05, 0.5, "E_{T} [GeV]");  // Moved closer
    
    TLatex* commonXTitle = new TLatex();
    commonXTitle->SetNDC();
    commonXTitle->SetTextFont(42);
    commonXTitle->SetTextSize(0.035);
    commonXTitle->SetTextAlign(22);
    commonXTitle->DrawLatex(0.4, 0.05, "#eta_{jet}");
    
    canvas->Update();
    
    // Save PDF
    canvas->SaveAs("jet_et_vs_eta_root.pdf");
    std::cout << "Saved: jet_et_vs_eta_root.pdf" << std::endl;
    
    // Print statistics
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "ANALYSIS SUMMARY" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    for (int fileIdx = 0; fileIdx < 4; fileIdx++) {
        std::cout << "\n" << datasets[fileIdx].label << ":" << std::endl;
        for (int eventIdx = 0; eventIdx < 2; eventIdx++) {
            int histIdx = fileIdx * 2 + eventIdx;
            if (histograms[histIdx]) {
                std::cout << "  " << eventLabels[eventIdx] << ": " 
                         << stats[histIdx].n_jets << " jets, "
                         << "<eta> = " << std::fixed << std::setprecision(2) 
                         << stats[histIdx].mean_eta << "±" << stats[histIdx].std_eta << ", "
                         << "<E_T> = " << std::setprecision(1) 
                         << stats[histIdx].mean_et << "±" << stats[histIdx].std_et << " GeV" << std::endl;
            } else {
                std::cout << "  " << eventLabels[eventIdx] << ": No data available" << std::endl;
            }
        }
    }
    
    std::cout << "\nAnalysis completed successfully!" << std::endl;
}