#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TString.h"
#include "TBox.h"

void plot_qq_differential_jetshapes() {
    // Set ROOT style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.05);

    // Use proper fonts
    gStyle->SetTextFont(42); // Helvetica
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    
    // Data for r values (x-axis)
    Double_t r_values[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    
    // Function to create and save plot for each eta bin
    auto createPlot = [&](const Double_t qq_64[10], const Double_t qq_105[10], const Double_t qq_141[10], const Double_t qq_300[10],
                         const Double_t gg_64[10], const Double_t gg_105[10], const Double_t gg_141[10], const Double_t gg_300[10],
                         const TString& etaLabel, const TString& fileName) {
        
        // Create TGraphs for QQ
        TGraph *gr_qq_64 = new TGraph(10, r_values, qq_64);
        TGraph *gr_qq_105 = new TGraph(10, r_values, qq_105);
        TGraph *gr_qq_141 = new TGraph(10, r_values, qq_141);
        TGraph *gr_qq_300 = new TGraph(10, r_values, qq_300);
        
        // Create TGraphs for GG
        TGraph *gr_gg_64 = new TGraph(10, r_values, gg_64);
        TGraph *gr_gg_105 = new TGraph(10, r_values, gg_105);
        TGraph *gr_gg_141 = new TGraph(10, r_values, gg_141);
        TGraph *gr_gg_300 = new TGraph(10, r_values, gg_300);
        
        // Create canvas
        TCanvas *c1 = new TCanvas("c1", "QQ Differential Jet Shapes", 800, 800);
        c1->SetLeftMargin(0.12);
        c1->SetBottomMargin(0.12);
        c1->SetRightMargin(0.05);
        c1->SetTopMargin(0.08);
        
        // Set graph styles - QQ in red, GG in blue, same line styles for corresponding energies
        
        // QQ graphs - red color with different line styles
        gr_qq_64->SetLineColor(kRed+1);
        gr_qq_64->SetLineStyle(3);  // dotted
        gr_qq_64->SetLineWidth(2);
        
        gr_qq_105->SetLineColor(kRed+1);
        gr_qq_105->SetLineStyle(2);  // dash-dot
        gr_qq_105->SetLineWidth(2);
        
        gr_qq_141->SetLineColor(kRed+1);
        gr_qq_141->SetLineStyle(10);  // dashed
        gr_qq_141->SetLineWidth(2);
        
        gr_qq_300->SetLineColor(kRed+1);
        gr_qq_300->SetLineStyle(1);  // solid
        gr_qq_300->SetLineWidth(2);
        
        // GG graphs - blue color with same line styles
        gr_gg_64->SetLineColor(kBlue+1);
        gr_gg_64->SetLineStyle(3);  // dotted
        gr_gg_64->SetLineWidth(2);
        
        gr_gg_105->SetLineColor(kBlue+1);
        gr_gg_105->SetLineStyle(2);  // dash-dot
        gr_gg_105->SetLineWidth(2);
        
        gr_gg_141->SetLineColor(kBlue+1);
        gr_gg_141->SetLineStyle(10);  // dashed
        gr_gg_141->SetLineWidth(2);
        
        gr_gg_300->SetLineColor(kBlue+1);
        gr_gg_300->SetLineStyle(1);  // solid
        gr_gg_300->SetLineWidth(2);
        
        // Draw the first graph with axes
        gr_qq_300->Draw("AC");
        
        // Set axis ranges and labels
        gr_qq_300->GetXaxis()->SetRangeUser(0.05, 1.0);
        gr_qq_300->GetYaxis()->SetRangeUser(0.0, 4.0);
        
        gr_qq_300->GetXaxis()->SetTitle("r");
        gr_qq_300->GetYaxis()->SetTitle("#rho(r)");
        
        gr_qq_300->GetXaxis()->SetTitleSize(0.045);
        gr_qq_300->GetYaxis()->SetTitleSize(0.045);
        gr_qq_300->GetXaxis()->SetLabelSize(0.04);
        gr_qq_300->GetYaxis()->SetLabelSize(0.04);
        gr_qq_300->GetXaxis()->SetTitleOffset(1.1);
        gr_qq_300->GetYaxis()->SetTitleOffset(1.2);
        // Center the axis titles
        gr_qq_300->GetXaxis()->CenterTitle();
        gr_qq_300->GetYaxis()->CenterTitle();
        
        // Draw all other graphs
        gr_qq_64->Draw("C SAME");
        gr_qq_105->Draw("C SAME");
        gr_qq_141->Draw("C SAME");
        
        gr_gg_64->Draw("C SAME");
        gr_gg_105->Draw("C SAME");
        gr_gg_141->Draw("C SAME");
        gr_gg_300->Draw("C SAME");
        
        // Create legend
        TLegend *legend = new TLegend(0.65, 0.55, 0.92, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.035);

        // Add colored boxes for QQ and GG
        TBox *qq_box = new TBox();
        qq_box->SetFillColor(kRed+1);
        TBox *gg_box = new TBox();
        gg_box->SetFillColor(kBlue+1);
        
        // Add header for legend
        // legend->SetHeader("QQ (red) / GG (blue)", "C");
        legend->SetHeader("#bf{#color[2]{QQ} / #color[4]{GG}}", "C");
        // legend->AddEntry(gr_qq_64, "ep 64 GeV", "");
        // legend->AddEntry(gr_qq_105, "ep 105 GeV", "");
        // legend->AddEntry(gr_qq_141, "ep 141 GeV", "");
        // legend->AddEntry(gr_qq_300, "ep 300 GeV", "");
        // Create black versions for legend display
        TGraph *gr_legend_64 = (TGraph*)gr_qq_64->Clone();
        gr_legend_64->SetLineColor(kBlack);

        TGraph *gr_legend_105 = (TGraph*)gr_qq_105->Clone();  
        gr_legend_105->SetLineColor(kBlack);

        TGraph *gr_legend_141 = (TGraph*)gr_qq_141->Clone();
        gr_legend_141->SetLineColor(kBlack);

        TGraph *gr_legend_300 = (TGraph*)gr_qq_300->Clone();
        gr_legend_300->SetLineColor(kBlack);

        // Use the black versions in legend
        legend->AddEntry(gr_legend_64, "ep 64 GeV", "l");
        legend->AddEntry(gr_legend_105, "ep 105 GeV", "l");
        legend->AddEntry(gr_legend_141, "ep 141 GeV", "l");
        legend->AddEntry(gr_legend_300, "ep 300 GeV", "l");
        
        legend->Draw();
        
        // // Add text for eta bin
        // TText *eta_text = new TText(0.15, 0.85, etaLabel);
        // eta_text->SetNDC();
        // eta_text->SetTextSize(0.04);
        // eta_text->SetTextFont(42);
        // eta_text->Draw();
        TLatex *eta_text = new TLatex();
        eta_text->SetNDC();
        eta_text->SetTextSize(0.04);
        eta_text->SetTextFont(42); // Helvetica font
        eta_text->DrawLatex(0.15, 0.85, etaLabel.Data());
            
        // Update canvas
        c1->Update();
        // Save the plot
        c1->SaveAs(fileName);
        cout << "Plot created for " << etaLabel << " -> " << fileName << endl;
        // Clean up
        delete c1;
    };
    
    // Define all data arrays for each eta bin
    
    // ETA BIN 1: (-1 < eta < 0)
    Double_t qq_64_eta1[10] = {2.7067, 2.1164, 1.3982, 0.8514, 0.5688, 0.396, 0.3346, 0.2481, 0.1452, 0.0181};
    Double_t qq_105_eta1[10] = {3.2144, 2.4553, 1.5482, 0.989, 0.6732, 0.478, 0.3704, 0.2872, 0.1802, 0.0227};
    Double_t qq_141_eta1[10] = {3.1969, 2.3892, 1.4322, 0.8947, 0.6225, 0.4458, 0.3458, 0.2681, 0.1735, 0.0244};
    Double_t qq_300_eta1[10] = {3.1127, 2.2566, 1.3594, 0.871, 0.6213, 0.4571, 0.3672, 0.3044, 0.2078, 0.0335};
    
    Double_t gg_64_eta1[10] = {1.2393, 1.7212, 1.4479, 1.1181, 0.9598, 0.9094, 0.5361, 0.6709, 0.2544, 0.0162};
    Double_t gg_105_eta1[10] = {1.1629, 1.7121, 1.5188, 1.2898, 0.9642, 0.7612, 0.7216, 0.5244, 0.2693, 0.0227};
    Double_t gg_141_eta1[10] = {1.2213, 1.5389, 1.477, 1.2965, 1.0624, 0.8243, 0.6386, 0.5299, 0.3096, 0.0457};
    Double_t gg_300_eta1[10] = {1.4184, 1.7653, 1.5791, 1.299, 1.0366, 0.8584, 0.7193, 0.5869, 0.3915, 0.064};
    
    // ETA BIN 2: (0 < eta < 1)
    Double_t qq_64_eta2[10] = {3.2144, 2.4553, 1.5482, 0.989, 0.6732, 0.478, 0.3704, 0.2872, 0.1802, 0.0227};
    Double_t qq_105_eta2[10] = {3.0471, 2.4194, 1.4517, 0.9281, 0.6276, 0.4586, 0.3492, 0.2701, 0.1728, 0.0229};
    Double_t qq_141_eta2[10] = {3.0451, 2.3809, 1.4758, 0.9425, 0.673, 0.4836, 0.3814, 0.3, 0.1988, 0.0283};
    Double_t qq_300_eta2[10] = {3.0955, 2.3324, 1.4234, 0.8844, 0.6038, 0.447, 0.3407, 0.2733, 0.1797, 0.0235};
    
    Double_t gg_64_eta2[10] = {1.4013, 1.8083, 1.6386, 1.2993, 1.0333, 0.822, 0.6873, 0.5343, 0.3137, 0.0494};
    Double_t gg_105_eta2[10] = {1.3564, 1.8137, 1.5139, 1.3047, 1.0164, 0.8221, 0.6815, 0.5232, 0.3261, 0.0516};
    Double_t gg_141_eta2[10] = {1.2158, 1.5969, 1.5567, 1.3275, 1.0829, 0.9092, 0.7596, 0.6151, 0.3974, 0.0674};
    Double_t gg_300_eta2[10] = {1.4663, 1.7894, 1.6186, 1.3163, 1.057, 0.8749, 0.7399, 0.6255, 0.4401, 0.0749};
    
    // ETA BIN 3: (1 < eta < 1.5)
    Double_t qq_64_eta3[10] = {3.1969, 2.3892, 1.4322, 0.8947, 0.6225, 0.4458, 0.3458, 0.2681, 0.1735, 0.0244};
    Double_t qq_105_eta3[10] = {3.1965, 2.3338, 1.4193, 0.8811, 0.6035, 0.4455, 0.348, 0.2743, 0.179, 0.026};
    Double_t qq_141_eta3[10] = {3.074, 2.2656, 1.4013, 0.9164, 0.6398, 0.4709, 0.3716, 0.3072, 0.205, 0.0308};
    Double_t qq_300_eta3[10] = {3.0471, 2.4194, 1.3875, 0.9281, 0.6276, 0.4586, 0.3492, 0.2701, 0.1728, 0.0229};
    
    Double_t gg_64_eta3[10] = {1.5157, 1.9342, 1.7226, 1.3443, 1.1092, 0.8719, 0.7479, 0.5836, 0.3877, 0.0596};
    Double_t gg_105_eta3[10] = {1.4394, 1.8063, 1.6565, 1.3035, 1.0745, 0.8891, 0.7483, 0.5888, 0.4094, 0.0607};
    Double_t gg_141_eta3[10] = {1.3305, 1.67, 1.5658, 1.3375, 1.1313, 0.9602, 0.7895, 0.6468, 0.4453, 0.0725};
    Double_t gg_300_eta3[10] = {1.499, 1.8235, 1.6476, 1.3386, 1.0851, 0.9005, 0.7639, 0.6471, 0.4643, 0.0799};
    
    // ETA BIN 4: (1.5 < eta < -2)
    Double_t qq_64_eta4[10] = {3.0729, 2.2427, 1.3847, 0.8767, 0.6124, 0.4351, 0.328, 0.2706, 0.1829, 0.0241};
    Double_t qq_105_eta4[10] = {3.2012, 2.3082, 1.4196, 0.9112, 0.6204, 0.46, 0.3616, 0.2799, 0.1963, 0.0282};
    Double_t qq_141_eta4[10] = {2.149, 2.0844, 1.5643, 1.1801, 0.8991, 0.7025, 0.5879, 0.4761, 0.3349, 0.0515};
    Double_t qq_300_eta4[10] = {3.1127, 2.2566, 1.3594, 0.871, 0.6213, 0.4571, 0.357, 0.3044, 0.2078, 0.0335};
    
    Double_t gg_64_eta4[10] = {1.4914, 1.7605, 1.5187, 1.4257, 1.1041, 0.8632, 0.7323, 0.5641, 0.3963, 0.0613};
    Double_t gg_105_eta4[10] = {1.4394, 1.8063, 1.6565, 1.3035, 1.0745, 0.8891, 0.7483, 0.5888, 0.4094, 0.0607};
    Double_t gg_141_eta4[10] = {1.3708, 1.7405, 1.5921, 1.3568, 1.168, 0.9773, 0.8288, 0.676, 0.4751, 0.0811};
    Double_t gg_300_eta4[10] = {1.4022, 1.7214, 1.5166, 1.225, 0.9666, 0.775, 0.6535, 0.5014, 0.328, 0.0502};
    
    // Create plots for all eta bins
    createPlot(qq_64_eta1, qq_105_eta1, qq_141_eta1, qq_300_eta1,
               gg_64_eta1, gg_105_eta1, gg_141_eta1, gg_300_eta1,
               "-1 < #eta < 0", "qq_gg_differential_jetshapes_eta_minus1_0.pdf");
    
    createPlot(qq_64_eta2, qq_105_eta2, qq_141_eta2, qq_300_eta2,
               gg_64_eta2, gg_105_eta2, gg_141_eta2, gg_300_eta2,
               "0 < #eta < 1", "qq_gg_differential_jetshapes_eta_0_1.pdf");
    
    createPlot(qq_64_eta3, qq_105_eta3, qq_141_eta3, qq_300_eta3,
               gg_64_eta3, gg_105_eta3, gg_141_eta3, gg_300_eta3,
               "1 < #eta < 1.5", "qq_gg_differential_jetshapes_eta_1_1p5.pdf");
    
    createPlot(qq_64_eta4, qq_105_eta4, qq_141_eta4, qq_300_eta4,
               gg_64_eta4, gg_105_eta4, gg_141_eta4, gg_300_eta4,
               "1.5 < #eta < 2", "qq_gg_differential_jetshapes_eta_1p5_minus2.pdf");
    
    cout << "\nAll plots created successfully!" << endl;
    cout << "Files saved:" << endl;
    cout << "  - qq_gg_differential_jetshapes_eta_minus1_0.pdf" << endl;
    cout << "  - qq_gg_differential_jetshapes_eta_0_1.pdf" << endl;
    cout << "  - qq_gg_differential_jetshapes_eta_1_1p5.pdf" << endl;
    cout << "  - qq_gg_differential_jetshapes_eta_1p5_minus2.pdf" << endl;
}