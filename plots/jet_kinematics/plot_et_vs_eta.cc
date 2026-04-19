#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include <vector>
#include <iostream>

void plot_et_vs_eta() {
    
    // Set ROOT style for better plots
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetNumberContours(50);
    // Modern serif font styling
    // gStyle->SetTextFont(42);        // Times-Roman serif font
    // gStyle->SetLabelFont(42, "XYZ"); // Apply to all axes
    gStyle->SetTextSizePixels(24);  
    gStyle->SetTitleFont(42, "XYZ"); // Apply to all titles
    gStyle->SetStatFont(42);        // Statistics box font
    
    // Open the ROOT file
    const char* filename = "/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic64_dijets_pT5/dijets_eic64_R10_EtMin0_10_7.root";

    // const char* hera300-filename = "/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic64_dijets_pT7/dijets_hera300_R10_EtMin0_10_7.root";
    // const char* eic141-filename = "/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic141_dijets_pT5/dijets_eic141_R10_EtMin0_10_7.root";
    // const char* eic105-filename = "/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic105_dijets_pT5/dijets_eic105_R10_EtMin0_10_7.root";
    // const char* eic64-filename = "/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic64_dijets_pT5/dijets_eic64_R10_EtMin0_10_7.root";

    TFile* file = TFile::Open(filename, "READ");
    
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    
    // Navigate to the GG_Events directory and get the tree
    // GG_Events
    TDirectory* dir = (TDirectory*)file->Get("GG_Events");
    if (!dir) {
        std::cerr << "Error: Cannot find directory" << std::endl;
        file->Close();
        return;
    }
    // jets_GG_Events
    TTree* tree = (TTree*)dir->Get("jets_GG_Events");
    if (!tree) {
        std::cerr << "Error: Cannot find tree" << std::endl;
        file->Close();
        return;
    }
    
    // Set up branches for reading
    std::vector<float> *jet_et = nullptr;
    std::vector<float> *jet_eta = nullptr;
    Int_t n_jets;
    
    tree->SetBranchAddress("jet_et", &jet_et);
    tree->SetBranchAddress("jet_eta", &jet_eta);
    tree->SetBranchAddress("n_jets", &n_jets);
    
    TH2F* h_et_vs_eta = new TH2F("h_et_vs_eta", 
                                // ";#eta;E_{T} [GeV]",
                                ";;",
                                60, -1.5, 4.0,   // eta range 
                                60, 5.0, 60.0);  // ET range
    
    // Get number of entries
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    
    // Event loop
    int totalJets = 0;
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        // Progress indicator
        if (i % 10000 == 0) {
            std::cout << "Processing event " << i << "/" << nEntries 
                     << " (" << (100.0 * i / nEntries) << "%)" << std::endl;
        }
        
        // Check if vectors are valid
        if (!jet_et || !jet_eta) continue;
        if (jet_et->size() != jet_eta->size()) continue;
        
        // Fill histogram for all jets in the event
        for (size_t j = 0; j < jet_et->size(); j++) {
            h_et_vs_eta->Fill((*jet_eta)[j], (*jet_et)[j]);
            totalJets++;
        }
    }
    
    std::cout << "Total jets processed: " << totalJets << std::endl;
    
    TCanvas* canvas = new TCanvas("canvas", "ET vs Eta Distribution", 1200, 1200);
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.12);  
    canvas->SetBottomMargin(0.12);
    canvas->SetTopMargin(0.08);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    
    // Draw the histogram without color bar
    h_et_vs_eta->Draw("COLZ");
    
    // h_et_vs_eta->GetXaxis()->SetTitle("");
    // h_et_vs_eta->GetYaxis()->SetTitle("");
    // TMathText* xTitle = new TMathText();
    // TMathText* yTitle = new TMathText();
    // xTitle->SetTextSize(0.045);xTitle->SetTextAlign(22);xTitle->DrawMathText(1.25, 7.5, "#eta"); 
    // yTitle->SetTextSize(0.045);yTitle->SetTextAlign(22);yTitle->SetTextAngle(90);yTitle->DrawMathText(-1.2, 17, "E_{T} [GeV]"); 

    // h_et_vs_eta->GetXaxis()->SetTitle("\\eta");
    // h_et_vs_eta->GetYaxis()->SetTitle("E_{T} \\text{ [GeV]}");

    h_et_vs_eta->GetXaxis()->SetTitle("#eta");
    h_et_vs_eta->GetYaxis()->SetTitle("E_{T} [GeV]");

    // Modern font and size styling
    h_et_vs_eta->GetXaxis()->SetTitleFont(42);
    h_et_vs_eta->GetYaxis()->SetTitleFont(42);
    h_et_vs_eta->GetXaxis()->SetLabelFont(42);
    h_et_vs_eta->GetYaxis()->SetLabelFont(42);

    h_et_vs_eta->GetXaxis()->SetTitleSize(0.045);
    h_et_vs_eta->GetYaxis()->SetTitleSize(0.045);
    h_et_vs_eta->GetXaxis()->SetLabelSize(0.040);
    h_et_vs_eta->GetYaxis()->SetLabelSize(0.040);

    h_et_vs_eta->GetXaxis()->SetTitleOffset(1.1);
    h_et_vs_eta->GetYaxis()->SetTitleOffset(1.2);
    // Center the axis titles
    h_et_vs_eta->GetXaxis()->CenterTitle(true);
    h_et_vs_eta->GetYaxis()->CenterTitle(true);

    h_et_vs_eta->GetXaxis()->SetAxisColor(kBlack);
    h_et_vs_eta->GetYaxis()->SetAxisColor(kBlack);
    h_et_vs_eta->GetXaxis()->SetTickLength(0.02);
    h_et_vs_eta->GetYaxis()->SetTickLength(0.02);
    
    // Extract filenames and tree names for dynamic naming
    TString fileNameOnly = gSystem->BaseName(filename);
    fileNameOnly.ReplaceAll(".root", "");  // Remove .root extension
    TString treeName = tree->GetName();
    
    // Create dynamic output filename
    TString outputName = Form("%s_%s_eTvsEta.pdf", fileNameOnly.Data(), treeName.Data());
    
    // Add statistics box
    canvas->Update();
    
    // Save only the high-quality PDF
    canvas->SaveAs(outputName.Data());
    // canvas->SaveAs("plot.tex");
    
    std::cout << "\nPlot saved as: " << outputName.Data() << std::endl;
    
    // Print some statistics
    std::cout << "\nHistogram Statistics:" << std::endl;
    std::cout << "  Entries: " << h_et_vs_eta->GetEntries() << std::endl;
    std::cout << "  Mean Eta: " << h_et_vs_eta->GetMean(1) << std::endl;
    std::cout << "  Mean ET: " << h_et_vs_eta->GetMean(2) << " GeV" << std::endl;
    std::cout << "  RMS Eta: " << h_et_vs_eta->GetRMS(1) << std::endl;
    std::cout << "  RMS ET: " << h_et_vs_eta->GetRMS(2) << " GeV" << std::endl;
    
    // Clean up
    file->Close();
    delete file;
    
    std::cout << "\nAnalysis completed successfully!" << std::endl;
}