#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPad.h>
#include <TStyle.h>
#include <TFile.h>

void plotintjetsubplots() {
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
    
    // Define r values and number of points
    const int nPoints = 10;
    double rValues[nPoints] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    // EIC 64 -4 to 4, 10 GeV Cut
    double quark_values[] = {0.2038, 0.4671, 0.6393, 0.7461, 0.8158, 0.8655, 0.9024, 0.9309, 0.9521, 0.9602};
    double gluon_values[] = {0.0913, 0.2701, 0.4369, 0.5765, 0.6899, 0.7747, 0.8412, 0.8970, 0.9366, 0.9513};
    double gluon_quark_values[] = {0.1703, 0.4066, 0.5789, 0.6968, 0.7789, 0.8397, 0.8853, 0.9210, 0.9473, 0.9579};
    double combined_values[] = {0.1828, 0.4296, 0.6016, 0.7151, 0.7927, 0.8492, 0.8916, 0.9247, 0.9492, 0.9587};
    double resolved_values[] = {0.1684, 0.4059, 0.5777, 0.6951, 0.7776, 0.8385, 0.8840, 0.9204, 0.9476, 0.9580};
    double direct_values[] = {0.1993, 0.4565, 0.6288, 0.7379, 0.8098, 0.8614, 0.9001, 0.9295, 0.9510, 0.9596};

    // EIC 105 -4 to 4, 10 GeV Cut 
    double quark_values_2[] = {0.2142, 0.4857, 0.6583, 0.7640, 0.8321, 0.8798, 0.9162, 0.9440, 0.9654, 0.9740};
    double gluon_values_2[] = {0.0928, 0.2769, 0.4493, 0.5954, 0.7082, 0.7940, 0.8635, 0.9183, 0.9586, 0.9755};
    double gluon_quark_values_2[] = {0.1739, 0.4140, 0.5870, 0.7056, 0.7883, 0.8496, 0.8961, 0.9335, 0.9613, 0.9730};
    double combined_values_2[] = {0.1848, 0.4340, 0.6068, 0.7220, 0.8008, 0.8582, 0.9022, 0.9368, 0.9629, 0.9736};
    double resolved_values_2[] = {0.1701, 0.4098, 0.5824, 0.7016, 0.7849, 0.8467, 0.8946, 0.9328, 0.9614, 0.9733};
    double direct_values_2[] = {0.2109, 0.4769, 0.6501, 0.7581, 0.8290, 0.8787, 0.9156, 0.9439, 0.9656, 0.9743};

    
    // EIC 141 -4 to 4, 10 GeV Cut
    double quark_values_3[] = {0.2127, 0.4807, 0.6506, 0.7546, 0.8226, 0.8705, 0.9068, 0.9359, 0.9585, 0.9680};
    double gluon_values_3[] = {0.0865, 0.2609, 0.4332, 0.5724, 0.6829, 0.7712, 0.8436, 0.9034, 0.9504, 0.9704};
    double gluon_quark_values_3[] = {0.1618, 0.3934, 0.5646, 0.6835, 0.7681, 0.8318, 0.8818, 0.9224, 0.9538, 0.9672};
    double combined_values_3[] = {0.1739, 0.4138, 0.5847, 0.6998, 0.7806, 0.8406, 0.8876, 0.9257, 0.9553, 0.9679};
    double resolved_values_3[] = {0.1582, 0.3872, 0.5579, 0.6769, 0.7626, 0.8275, 0.8792, 0.9216, 0.9548, 0.9690};
    double direct_values_3[] = {0.2088, 0.4730, 0.6440, 0.7505, 0.8205, 0.8697, 0.9063, 0.9349, 0.9564, 0.9653};

    
    // Hera -4 to 4, 17 GeV Cut
    double quark_values_4[] = {0.3150, 0.5863, 0.7289, 0.8101, 0.8611, 0.8980, 0.9256, 0.9479, 0.9661, 0.9758};
    double gluon_values_4[] = {0.1328, 0.3452, 0.5191, 0.6463, 0.7407, 0.8110, 0.8698, 0.9156, 0.9534, 0.9736};
    double gluon_quark_values_4[] = {0.2356, 0.4815, 0.6420, 0.7447, 0.8130, 0.8631, 0.9025, 0.9333, 0.9590, 0.9720};
    double combined_values_4[] = {0.2504, 0.5010, 0.6568, 0.7550, 0.8206, 0.8687, 0.9064, 0.9362, 0.9608, 0.9736};
    double resolved_values_4[] = {0.2318, 0.4752, 0.6333, 0.7367, 0.8068, 0.8583, 0.8995, 0.9323, 0.9595, 0.9737};
    double direct_values_4[] = {0.3015, 0.5718, 0.7210, 0.8054, 0.8586, 0.8970, 0.9253, 0.9469, 0.9646, 0.9734};
    
    // Create main canvas
    TCanvas *c = new TCanvas("c", "Integrated Jet Shapes", 1200, 1200);
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
    
    double* combinedData[4] = {gluon_quark_values, gluon_quark_values_2, gluon_quark_values_3, gluon_quark_values_4};
    double* quarkData[4] = {quark_values, quark_values_2, quark_values_3, quark_values_4};
    double* gluonData[4] = {gluon_values, gluon_values_2, gluon_values_3, gluon_values_4};
    
    // Store graphs for ROOT file
    TGraph *gCombined[4], *gQuark[4], *gGluon[4];
    
    // Process each subplot
    for (int i = 0; i < 4; i++) {
        pads[i]->cd();
        pads[i]->SetLogy(1); 
        gStyle->SetOptTitle(0);
        // No grid
        pads[i]->SetGridx(0);
        pads[i]->SetGridy(0);
        pads[i]->SetFrameBorderSize(0);
        
        // Create graphs
        gCombined[i] = new TGraph(nPoints, rValues, combinedData[i]);
        gQuark[i] = new TGraph(nPoints, rValues, quarkData[i]);
        gGluon[i] = new TGraph(nPoints, rValues, gluonData[i]);
        
        // Style for Combined - Blue
        gCombined[i]->SetLineWidth(3);
        gCombined[i]->SetLineColor(kBlue);
        gCombined[i]->SetMarkerStyle(20);
        gCombined[i]->SetMarkerColor(kBlue);
        gCombined[i]->SetMarkerSize(1.0);
        
        // Style for Quark - Orange/Red
        gQuark[i]->SetLineWidth(3);
        gQuark[i]->SetLineColor(kOrange+7);
        gQuark[i]->SetMarkerStyle(21);
        gQuark[i]->SetMarkerColor(kOrange+7);
        gQuark[i]->SetMarkerSize(1.0);
        
        // Style for Gluon - Green
        gGluon[i]->SetLineWidth(3);
        gGluon[i]->SetLineColor(kGreen+2);
        gGluon[i]->SetMarkerStyle(22);
        gGluon[i]->SetMarkerColor(kGreen+2);
        gGluon[i]->SetMarkerSize(1.0);
        
        // Set axis ranges for integrated jet shapes (0 to 1)
        gCombined[i]->SetMinimum(0.06);
        gCombined[i]->SetMaximum(1.1);
        gCombined[i]->GetXaxis()->SetLimits(0.0, 1.1);
        
        // Set axis labels and title sizes
        gCombined[i]->GetXaxis()->SetLabelSize(0.06);
        gCombined[i]->GetYaxis()->SetLabelSize(0.06);
        gCombined[i]->GetXaxis()->SetTitleSize(0.06);
        gCombined[i]->GetYaxis()->SetTitleSize(0.06);
        
        // Draw graphs - using "ALP" for first graph, "LP SAME" for others
        gCombined[i]->Draw("APC");
        gQuark[i]->Draw("PC SAME");
        gGluon[i]->Draw("PC SAME");
        
        // Add energy and eta range labels
        TLatex *t = new TLatex();
        t->SetTextFont(42);
        t->SetTextSize(0.05);
        t->SetTextAlign(31);
        t->DrawLatex(1.05, 0.12, titles[i]);
        
        if (i == 0){
            t->DrawLatex(1.05, 0.15, Form("E_{T}^{jet} > 10 GeV"));
            t->DrawLatex(1.05, 0.075, "ep (64 GeV)");
        }
        if (i == 1){
            t->DrawLatex(1.05, 0.15, Form("E_{T}^{jet} > 10 GeV"));
            t->DrawLatex(1.05, 0.075, "ep (105 GeV)");
        }
            
        if (i == 2){
            t->DrawLatex(1.05, 0.15, Form("E_{T}^{jet} > 10 GeV"));
            t->DrawLatex(1.05, 0.075, "ep (141 GeV)");
        }
            
        if (i == 3){
            t->DrawLatex(1.05, 0.15, Form("E_{T}^{jet} > 17 GeV"));
            t->DrawLatex(1.05, 0.075, "ep (300 GeV)");
        }
            

        // Add legend only on the first panel
        if (i == 1) {
            TLegend *leg = new TLegend(0.287, 0.056, 0.637, 0.256);
            leg->AddEntry(gQuark[i], "QQ", "lp");
            leg->AddEntry(gGluon[i], "GG", "lp");
            leg->AddEntry(gCombined[i], "GQ", "lp");
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
    yLabel.DrawLatexNDC(0.02, 0.50, "#psi(r)");

    // Horizontal label for x-axis
    TLatex xLabel;
    xLabel.SetTextFont(42);
    xLabel.SetTextSize(0.03);
    xLabel.SetTextAlign(22); // Center alignment
    xLabel.DrawLatexNDC(0.50, 0.02, "r");

    // Save to file
    c->SaveAs("integrated_jet_shapes_all.pdf");
    
    // Save to ROOT file (optional)
    TFile *outputFile = new TFile("integrated_jet_shapes_all.root", "RECREATE");
    c->Write("canvas");
    
    for (int i = 0; i < 4; i++) {
        gCombined[i]->Write(Form("gCombined_%d", i));
        gQuark[i]->Write(Form("gQuark_%d", i));
        gGluon[i]->Write(Form("gGluon_%d", i));
    }
    
    outputFile->Close();
    
    std::cout << "Plot saved to integrated_jet_shapes_all.pdf" << std::endl;
}