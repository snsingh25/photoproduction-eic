// For all eta ranges EIC, HERA -4 to 4

#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPad.h>
#include <TStyle.h>
#include <TFile.h>
#include <THStack.h>

void plotdiffjetsubplotsall() {
    // Set global style
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetFrameBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    
    // Set non-bold fonts
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42,"x");
    gStyle->SetLabelFont(42,"y");
    gStyle->SetTitleFont(42,"x");
    gStyle->SetTitleFont(42,"y");
    
    // Increase tick sizes
    gStyle->SetTickLength(0.03, "x");
    gStyle->SetTickLength(0.03, "y");
    gStyle->SetLabelSize(0.05, "xyz");
    
    // Define r values and number of bins
    const int nPoints = 10;
    // double rBins[nPoints+1] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
    double rBins[nPoints+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};


    // EIC 64 -4 to 4, 10 GeV Cut
    double quark_values[] = {4.0742, 3.0304, 1.8749, 1.1804, 0.7980, 0.5787, 0.4420, 0.3400, 0.2134, 0.0288};
    double gluon_values[] = {1.8683, 2.2251, 1.8797, 1.5105, 1.1599, 0.9186, 0.7097, 0.5976, 0.3276, 0.0482};
    double gluon_quark_values[] = {3.3582, 2.8521, 1.9543, 1.3377, 0.9575, 0.7191, 0.5425, 0.4236, 0.2703, 0.0368};
    double combined_values[] = {3.6359, 2.9099, 1.9129, 1.2698, 0.8899, 0.6606, 0.5017, 0.3912, 0.2455, 0.0335};
    double resolved_values[] = {3.2562, 2.7334, 1.8658, 1.2795, 0.9139, 0.6850, 0.5252, 0.4176, 0.2573, 0.0360};
    double direct_values[] = {4.0677, 3.1106, 1.9663, 1.2588, 0.8625, 0.6327, 0.4749, 0.3612, 0.2322, 0.0306};

    // EIC 105 -4 to 4, 10 GeV Cut 
    double quark_values_2[] = {4.2309, 3.1495, 1.8968, 1.1836, 0.7812, 0.5813, 0.4403, 0.3446, 0.2254, 0.0309};
    double gluon_values_2[] = {1.9396, 2.3131, 1.9770, 1.6216, 1.2180, 0.9463, 0.7692, 0.5830, 0.3885, 0.0570};
    double gluon_quark_values_2[] = {3.4005, 2.8527, 1.9469, 1.3452, 0.9682, 0.7254, 0.5670, 0.4470, 0.2935, 0.0418};
    double combined_values_2[] = {3.6430, 2.9379, 1.9275, 1.2974, 0.9071, 0.6807, 0.5283, 0.4134, 0.2716, 0.0383};
    double resolved_values_2[] = {3.2679, 2.7660, 1.8852, 1.3104, 0.9391, 0.7112, 0.5635, 0.4386, 0.2888, 0.0416};
    double direct_values_2[] = {4.3067, 3.2420, 2.0023, 1.2744, 0.8504, 0.6267, 0.4660, 0.3689, 0.2412, 0.0325};

    // EIC 141 -4 to 4, 10 GeV Cut
    double quark_values_3[] = {4.3648, 3.2278, 1.9475, 1.2198, 0.8301, 0.6095, 0.4731, 0.3819, 0.2572, 0.0358};
    double gluon_values_3[] = {2.0503, 2.5228, 2.1952, 1.7422, 1.3845, 1.1140, 0.9263, 0.7585, 0.5138, 0.0828};
    double gluon_quark_values_3[] = {3.3891, 2.9446, 2.0608, 1.4483, 1.0586, 0.8130, 0.6512, 0.5243, 0.3566, 0.0528};
    double combined_values_3[] = {3.6331, 3.0118, 2.0302, 1.3890, 1.0028, 0.7645, 0.6100, 0.4929, 0.3340, 0.0493};
    double resolved_values_3[] = {3.2868, 2.8815, 2.0270, 1.4331, 1.0585, 0.8216, 0.6673, 0.5448, 0.3719, 0.0559};
    double direct_values_3[] = {4.4014, 3.3006, 2.0372, 1.2912, 0.8793, 0.6378, 0.4829, 0.3778, 0.2500, 0.0345};
    
    // Hera -4 to 4, 17 GeV Cut
    double quark_values_4[] = {5.2709, 2.7920, 1.5292, 0.9040, 0.6041, 0.4580, 0.3398, 0.2784, 0.2252, 0.0410};
    double gluon_values_4[] = {2.5964, 2.4778, 1.8710, 1.3170, 1.0125, 0.7854, 0.6305, 0.5161, 0.3929, 0.0744};
    double gluon_quark_values_4[] = {4.1289, 2.7390, 1.7630, 1.1515, 0.8064, 0.6157, 0.4709, 0.3915, 0.2924, 0.0550};
    double combined_values_4[] = {4.3347, 2.7223, 1.6925, 1.0839, 0.7610, 0.5815, 0.4450, 0.3674, 0.2817, 0.0526};
    double resolved_values_4[] = {3.9234, 2.6240, 1.6990, 1.1155, 0.7961, 0.6121, 0.4763, 0.3964, 0.3038, 0.0566};
    double direct_values_4[] = {5.4626, 2.9920, 1.6748, 0.9975, 0.6647, 0.4976, 0.3592, 0.2879, 0.2209, 0.0414};
    
    // Create main canvas
    TCanvas *c = new TCanvas("c", "Differential Jet Shapes", 1200, 1200);
    c->SetMargin(0, 0, 0, 0);
    
    // Create four pads
    TPad *pad1 = new TPad("pad1", "pad1", 0.00, 0.50, 0.50, 1.00);
    TPad *pad2 = new TPad("pad2", "pad2", 0.50, 0.50, 1.00, 1.00);
    TPad *pad3 = new TPad("pad3", "pad3", 0.00, 0.00, 0.50, 0.50);
    TPad *pad4 = new TPad("pad4", "pad4", 0.50, 0.00, 1.00, 0.50);
    
    // Set pad margins
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.00);
    pad1->SetTopMargin(0.10);
    pad1->SetBottomMargin(0.00);
    
    pad2->SetLeftMargin(0.00);
    pad2->SetRightMargin(0.05);
    pad2->SetTopMargin(0.10);
    pad2->SetBottomMargin(0.00);
    
    pad3->SetLeftMargin(0.15);
    pad3->SetRightMargin(0.00);
    pad3->SetTopMargin(0.00);
    pad3->SetBottomMargin(0.15);
    
    pad4->SetLeftMargin(0.00);
    pad4->SetRightMargin(0.05);
    pad4->SetTopMargin(0.00);
    pad4->SetBottomMargin(0.15);
    
    pad1->SetFillStyle(0); 
    pad2->SetFillStyle(0);
    pad3->SetFillStyle(0);
    pad4->SetFillStyle(0);
    
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();
    
    // Data arrays and pad arrays
    TPad *pads[4] = {pad1, pad2, pad3, pad4};
    const char* titles[4] = {"-4 < #eta_{jet} < 4", "-4 < #eta_{jet} < 4", 
                            "-4 < #eta_{jet} < 4", "-4 < #eta_{jet} < 4"};
    
    double* combinedData[4] = {combined_values, combined_values_2, combined_values_3, combined_values_4};
    double* quarkData[4] = {quark_values, quark_values_2, quark_values_3, quark_values_4};
    double* gluonData[4] = {gluon_values, gluon_values_2, gluon_values_3, gluon_values_4};
    
    // Store histograms for ROOT file
    TH1F *hCombined[4], *hQuark[4], *hGluon[4];
    
    // Process each subplot
    for (int i = 0; i < 4; i++) {
        pads[i]->cd();
        // No grid
        pads[i]->SetGridx(0);
        pads[i]->SetGridy(0);
        pads[i]->SetFrameBorderSize(0);
        
        // Create histograms
        hCombined[i] = new TH1F(Form("hCombined_%d", i), "", nPoints, rBins);
        hQuark[i] = new TH1F(Form("hQuark_%d", i), "", nPoints, rBins);
        hGluon[i] = new TH1F(Form("hGluon_%d", i), "", nPoints, rBins);
        
        // Fill histograms
        for (int j = 0; j < nPoints; j++) {
            hCombined[i]->SetBinContent(j+1, combinedData[i][j]);
            hQuark[i]->SetBinContent(j+1, quarkData[i][j]);
            hGluon[i]->SetBinContent(j+1, gluonData[i][j]);
        }
        
        // Style for Combined - Blue (to match example)
        hCombined[i]->SetLineWidth(2);
        hCombined[i]->SetLineColor(kBlue);
        
        // Style for Quark - Orange/Red (to match example)
        hQuark[i]->SetLineWidth(2);
        hQuark[i]->SetLineColor(kOrange+7);
        
        // Style for Gluon - Green (to match example)
        hGluon[i]->SetLineWidth(2);
        hGluon[i]->SetLineColor(kGreen+2);
        
        // // Set axis titles and range
        // if (i >= 2) { // Bottom row
        //     hCombined[i]->GetXaxis()->SetTitle("Jet Radius (r)");
        //     hCombined[i]->GetXaxis()->SetLabelSize(0.06);
        //     hCombined[i]->GetXaxis()->SetTitleSize(0.06);
        //     hCombined[i]->GetXaxis()->SetTitleOffset(1.0);
        // } else { // Top row
        //     hCombined[i]->GetXaxis()->SetLabelSize(0);
        // }
        
        // if (i % 2 == 0) { // Left column
        //     hCombined[i]->GetYaxis()->SetTitle("Differential Jet Shape (#rho)");
        //     hCombined[i]->GetYaxis()->SetLabelSize(0.06);
        //     hCombined[i]->GetYaxis()->SetTitleSize(0.06);
        //     hCombined[i]->GetYaxis()->SetTitleOffset(1.2);
        // } else { // Right column
        //     hCombined[i]->GetYaxis()->SetLabelSize(0);
        // }
        
        // Set axis ranges to match example image
        hCombined[i]->SetMinimum(0);
        hCombined[i]->SetMaximum(5.5);
        hCombined[i]->GetXaxis()->SetRangeUser(0, 1.0);
        
        // Draw histograms - using HIST for step-like appearance as in example
        hCombined[i]->Draw("HIST");
        hQuark[i]->Draw("HIST SAME");
        hGluon[i]->Draw("HIST SAME");
        
        // Add energy and eta range labels
        TLatex *t = new TLatex();
        t->SetTextFont(42);
        t->SetTextSize(0.05);
        t->SetTextAlign(31);
        t->DrawLatex(0.966, 1.779, titles[i]);

        if (i == 0){
            t->DrawLatex(0.966, 1.282, Form("E_{T}^{jet} > 10 GeV"));
            t->DrawLatex(0.966, 2.302, "ep (64 GeV)");
        }
        if (i == 1){
            t->DrawLatex(0.966, 1.282, Form("E_{T}^{jet} > 10 GeV"));
            t->DrawLatex(0.966, 2.302, "ep (105 GeV)");
        }
            
        if (i == 2){
            t->DrawLatex(0.966, 1.282, Form("E_{T}^{jet} > 10 GeV"));
            t->DrawLatex(0.966, 2.302, "ep (141 GeV)");
        }
            
        if (i == 3){
            t->DrawLatex(0.966, 1.282, Form("E_{T}^{jet} > 17 GeV"));
            t->DrawLatex(0.966, 2.302, "ep (300 GeV)");
        }
        
        // Add legend only on the first panel
        if (i == 1) {
            // t->DrawLatex(1.065, 2.486, "ep (64 GeV)");
            TLegend *leg = new TLegend(0.58, 0.56, 0.98, 0.76);
            leg->AddEntry(hQuark[i], "Quark", "l");
            leg->AddEntry(hGluon[i], "Gluon", "l");
            leg->AddEntry(hCombined[i], "Combined", "l");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.05);
            leg->SetTextFont(42);
            leg->Draw();
        }
    }
    // Add shared axis labels on the canvas
    c->cd();

    // Vertical label for y-axis
    TLatex yLabel;
    yLabel.SetTextFont(42);
    yLabel.SetTextSize(0.03);
    yLabel.SetTextAlign(22); // Center alignment
    yLabel.SetTextAngle(90);
    yLabel.DrawLatexNDC(0.02, 0.50, "#rho");

    // Horizontal label for x-axis
    TLatex xLabel;
    xLabel.SetTextFont(42);
    xLabel.SetTextSize(0.03);
    xLabel.SetTextAlign(22); // Center alignment
    xLabel.DrawLatexNDC(0.50, 0.02, "r");

    // Save to file
    c->SaveAs("differential_jet_shapes_all.pdf");
    // c->SaveAs("differential_jet_shapes_hera.png");
    
    // Save to ROOT file
    // TFile *outputFile = new TFile("differential_jet_shapes_hera.root", "RECREATE");
    // c->Write("canvas");
    
    for (int i = 0; i < 4; i++) {
        hCombined[i]->Write();
        hQuark[i]->Write();
        hGluon[i]->Write();
    }
    
    // outputFile->Close();
    
    std::cout << "Plot saved to differential_jet_shapes_all.pdf" << std::endl;
}