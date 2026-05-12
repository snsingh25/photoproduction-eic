#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPad.h>
#include <TStyle.h>
#include <TFile.h>
#include <THStack.h>

void plotdiffjetsubplots() {
    // Set global style
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetFrameBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    
    // // Set non-bold fonts
    // gStyle->SetTextFont(42);
    // gStyle->SetLabelFont(42,"x");
    // gStyle->SetLabelFont(42,"y");
    // gStyle->SetTitleFont(42,"x");
    // gStyle->SetTitleFont(42,"y");

    // Set serif fonts (Times-Roman)
    gStyle->SetTextFont(132);     // Times-Roman serif
    gStyle->SetLabelFont(132,"x");
    gStyle->SetLabelFont(132,"y");
    gStyle->SetTitleFont(132,"x");
    gStyle->SetTitleFont(132,"y");
    gStyle->SetLegendFont(132);   // Add this line
    
    // Increase tick sizes
    gStyle->SetTickLength(0.03, "x");
    gStyle->SetTickLength(0.03, "y");
    gStyle->SetLabelSize(0.05, "xyz");
    
    // Define r values and number of bins
    const int nPoints = 10;
    // double rBins[nPoints+1] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
    double rBins[nPoints+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    // // EIC 64 -4 < eta < 4 ET_j > 10GeV
    double quark_values[] = {4.0742, 3.0304, 1.8749, 1.1804, 0.7980, 0.5787, 0.4420, 0.3400, 0.2134, 0.0288};
    double gluon_values[] = {1.8683, 2.2251, 1.8797, 1.5105, 1.1599, 0.9186, 0.7097, 0.5976, 0.3276, 0.0482};
    double gluon_quark_values[] = {3.3582, 2.8521, 1.9543, 1.3377, 0.9575, 0.7191, 0.5425, 0.4236, 0.2703, 0.0368};
    double combined_values[] = {3.6359, 2.9099, 1.9129, 1.2698, 0.8899, 0.6606, 0.5017, 0.3912, 0.2455, 0.0335};
    double resolved_values[] = {3.2562, 2.7334, 1.8658, 1.2795, 0.9139, 0.6850, 0.5252, 0.4176, 0.2573, 0.0360};
    double direct_values[] = {4.0677, 3.1106, 1.9663, 1.2588, 0.8625, 0.6327, 0.4749, 0.3612, 0.2322, 0.0306};
    // // EIC 105 -4 < eta < 4 ET_j > 10GeV
    double quark_values_2[] = {4.2309, 3.1495, 1.8968, 1.1836, 0.7812, 0.5813, 0.4403, 0.3446, 0.2254, 0.0309};
    double gluon_values_2[] = {1.9396, 2.3131, 1.9770, 1.6216, 1.2180, 0.9463, 0.7692, 0.5830, 0.3885, 0.0570};
    double gluon_quark_values_2[] = {3.4005, 2.8527, 1.9469, 1.3452, 0.9682, 0.7254, 0.5670, 0.4470, 0.2935, 0.0418};
    double combined_values_2[] = {3.6430, 2.9379, 1.9275, 1.2974, 0.9071, 0.6807, 0.5283, 0.4134, 0.2716, 0.0383};
    double resolved_values_2[] = {3.2679, 2.7660, 1.8852, 1.3104, 0.9391, 0.7112, 0.5635, 0.4386, 0.2888, 0.0416};
    double direct_values_2[] = {4.3067, 3.2420, 2.0023, 1.2744, 0.8504, 0.6267, 0.4660, 0.3689, 0.2412, 0.0325};
    // // EIC 141 -4 < eta < 4 ET_j > 10GeV
    double quark_values_3[] = {4.3648, 3.2278, 1.9475, 1.2198, 0.8301, 0.6095, 0.4731, 0.3819, 0.2572, 0.0358};
    double gluon_values_3[] = {2.0503, 2.5228, 2.1952, 1.7422, 1.3845, 1.1140, 0.9263, 0.7585, 0.5138, 0.0828};
    double gluon_quark_values_3[] = {3.3891, 2.9446, 2.0608, 1.4483, 1.0586, 0.8130, 0.6512, 0.5243, 0.3566, 0.0528};
    double combined_values_3[] = {3.6331, 3.0118, 2.0302, 1.3890, 1.0028, 0.7645, 0.6100, 0.4929, 0.3340, 0.0493};
    double resolved_values_3[] = {3.2868, 2.8815, 2.0270, 1.4331, 1.0585, 0.8216, 0.6673, 0.5448, 0.3719, 0.0559};
    double direct_values_3[] = {4.4014, 3.3006, 2.0372, 1.2912, 0.8793, 0.6378, 0.4829, 0.3778, 0.2500, 0.0345};
    // // HERA 300 -4 < eta < 4 ET_j > 17GeV
    double quark_values_4[] = {5.2709, 2.7920, 1.5292, 0.9040, 0.6041, 0.4580, 0.3398, 0.2784, 0.2252, 0.0410};
    double gluon_values_4[] = {2.5964, 2.4778, 1.8710, 1.3170, 1.0125, 0.7854, 0.6305, 0.5161, 0.3929, 0.0744};
    double gluon_quark_values_4[] = {4.1289, 2.7390, 1.7630, 1.1515, 0.8064, 0.6157, 0.4709, 0.3915, 0.2924, 0.0550};
    double combined_values_4[] = {4.3347, 2.7223, 1.6925, 1.0839, 0.7610, 0.5815, 0.4450, 0.3674, 0.2817, 0.0526};
    double resolved_values_4[] = {3.9234, 2.6240, 1.6990, 1.1155, 0.7961, 0.6121, 0.4763, 0.3964, 0.3038, 0.0566};
    double direct_values_4[] = {5.4626, 2.9920, 1.6748, 0.9975, 0.6647, 0.4976, 0.3592, 0.2879, 0.2209, 0.0414};

    // // EIC 140 GeV data: ET_j > 17GeV 
    // // -1 to 0
    // double quark_values[] = {2.5949, 2.1248, 1.3888, 0.8835, 0.5965, 0.4390, 0.3344, 0.2520, 0.1504, 0.0199};
    // double gluon_values[] = {1.6266, 1.7094, 1.5375, 1.0247, 1.0366, 0.7110, 0.5932, 0.4829, 0.2439, 0.0518};
    // double gluon_quark_values[] = {2.3634, 2.1100, 1.3915, 0.9234, 0.6640, 0.4841, 0.3758, 0.2839, 0.1641, 0.0219};
    // double combined_values[] = {2.4708, 2.1086, 1.3935, 0.9043, 0.6363, 0.4651, 0.3586, 0.2713, 0.1586, 0.0215};
    // double resolved_values[] = {2.3069, 2.0324, 1.4097, 0.9755, 0.7033, 0.5152, 0.3809, 0.3069, 0.1770, 0.0251};
    // double direct_values[] = {2.6250, 2.1803, 1.3783, 0.8372, 0.5734, 0.4179, 0.3377, 0.2379, 0.1413, 0.0181};
    // // 0 to 1 
    // double quark_values_2[] = {3.2965, 2.4950, 1.5614, 0.9744, 0.6688, 0.4778, 0.3671, 0.2752, 0.1755, 0.0226};
    // double gluon_values_2[] = {1.5318, 1.8849, 1.6213, 1.3365, 1.0034, 0.8205, 0.6185, 0.5090, 0.2906, 0.0395};
    // double gluon_quark_values_2[] = {2.8258, 2.3630, 1.6010, 1.0929, 0.7693, 0.5840, 0.4344, 0.3381, 0.2143, 0.0289};
    // double combined_values_2[] = {2.9982, 2.4066, 1.5829, 1.0461, 0.7307, 0.5428, 0.4097, 0.3150, 0.1989, 0.0263};
    // double resolved_values_2[] = {2.6622, 2.2809, 1.5752, 1.0636, 0.7596, 0.5732, 0.4367, 0.3413, 0.2135, 0.0288};
    // double direct_values_2[] = {3.3524, 2.5392, 1.5909, 1.0276, 0.7003, 0.5107, 0.3812, 0.2872, 0.1834, 0.0237};
    // // 1 to 1.5
    // double quark_values_3[] = {3.3087, 2.3549, 1.4268, 0.9090, 0.5935, 0.4408, 0.3340, 0.2692, 0.1614, 0.0237};
    // double gluon_values_3[] = {1.6867, 1.9927, 1.6734, 1.3889, 1.0440, 0.8043, 0.6255, 0.5399, 0.2790, 0.0433};
    // double gluon_quark_values_3[] = {2.6714, 2.2927, 1.5745, 1.0834, 0.7947, 0.5878, 0.4453, 0.3525, 0.2305, 0.0325};
    // double combined_values_3[] = {2.9347, 2.3098, 1.5080, 1.0129, 0.7090, 0.5266, 0.3997, 0.3206, 0.1995, 0.0287};
    // double resolved_values_3[] = {2.7451, 2.2614, 1.5077, 1.0642, 0.7574, 0.5568, 0.4352, 0.3531, 0.2116, 0.0312};
    // double direct_values_3[] = {3.1415, 2.3625, 1.5083, 0.9570, 0.6561, 0.4936, 0.3609, 0.2851, 0.1862, 0.0260};
    // // 1.5 to 2.0
    // double quark_values_4[] = {3.1835, 2.2921, 1.3523, 0.8728, 0.5699, 0.4308, 0.3095, 0.2507, 0.1657, 0.0227};
    // double gluon_values_4[] = {1.8487, 2.1424, 1.6471, 1.2406, 0.9744, 0.7582, 0.6141, 0.5338, 0.2813, 0.0469};
    // double gluon_quark_values_4[] = {2.5221, 2.1825, 1.5396, 1.0580, 0.7943, 0.5820, 0.4501, 0.3413, 0.2264, 0.0304};
    // double combined_values_4[] = {2.8193, 2.2351, 1.4516, 0.9745, 0.6912, 0.5151, 0.3879, 0.3052, 0.1988, 0.0273};
    // double resolved_values_4[] = {2.7354, 2.1753, 1.4490, 1.0053, 0.7215, 0.5486, 0.4170, 0.3369, 0.2122, 0.0307};
    // double direct_values_4[] = {2.9219, 2.3081, 1.4548, 0.9369, 0.6542, 0.4741, 0.3523, 0.2665, 0.1825, 0.0233};


    // // HERA data: ET_j > 17GeV
    // // -1 to 0
    // double quark_values[] = {3.6176, 1.9996, 1.0736, 0.6431, 0.4157, 0.3105, 0.2438, 0.1935, 0.1404, 0.0225};
    // double gluon_values[] = {1.7372, 1.8840, 1.4762, 1.0965, 0.8233, 0.6661, 0.5107, 0.4032, 0.2900, 0.0496};
    // double gluon_quark_values[] = {3.1398, 1.9317, 1.1541, 0.7473, 0.5178, 0.3894, 0.3109, 0.2451, 0.1793, 0.0311};
    // double combined_values[] = {3.2719, 1.9600, 1.1372, 0.7212, 0.4899, 0.3704, 0.2925, 0.2313, 0.1683, 0.0283};
    // double resolved_values[] = {3.0212, 1.9357, 1.1780, 0.7696, 0.5391, 0.4146, 0.3298, 0.2624, 0.1912, 0.0324};
    // double direct_values[] = {3.6811, 1.9996, 1.0706, 0.6422, 0.4096, 0.2983, 0.2316, 0.1805, 0.1311, 0.0216};
    // // 0 to 1   
    // double quark_values_2[] = {3.9771, 2.1735, 1.1460, 0.6930, 0.4625, 0.3359, 0.2651, 0.2164, 0.1617, 0.0300};
    // double gluon_values_2[] = {1.9279, 2.0059, 1.5531, 1.1426, 0.8753, 0.6899, 0.5493, 0.4554, 0.3365, 0.0640};
    // double gluon_quark_values_2[] = {3.2383, 2.1184, 1.3082, 0.8565, 0.6071, 0.4576, 0.3627, 0.2943, 0.2228, 0.0408};
    // double combined_values_2[] = {3.4005, 2.1290, 1.2681, 0.8202, 0.5766, 0.4326, 0.3427, 0.2796, 0.2099, 0.0389};
    // double resolved_values_2[] = {3.0484, 2.0791, 1.3030, 0.8651, 0.6241, 0.4758, 0.3808, 0.3133, 0.2358, 0.0439};
    // double direct_values_2[] = {4.1354, 2.2330, 1.1953, 0.7266, 0.4774, 0.3425, 0.2632, 0.2091, 0.1558, 0.0284};
    // // 1 to 1.5
    // double quark_values_3[] = {3.8519, 2.0592, 1.0964, 0.6533, 0.4320, 0.3191, 0.2586, 0.2084, 0.1656, 0.0309};
    // double gluon_values_3[] = {2.0687, 2.0465, 1.5607, 1.1331, 0.8458, 0.6574, 0.5559, 0.4649, 0.3704, 0.0723};
    // double gluon_quark_values_3[] = {3.0965, 2.0803, 1.3350, 0.8697, 0.6270, 0.4687, 0.3734, 0.3113, 0.2375, 0.0456};
    // double combined_values_3[] = {3.2478, 2.0681, 1.2745, 0.8224, 0.5820, 0.4368, 0.3538, 0.2924, 0.2276, 0.0435};
    // double resolved_values_3[] = {2.9957, 2.0554, 1.3231, 0.8760, 0.6357, 0.4790, 0.3949, 0.3293, 0.2571, 0.0493};
    // double direct_values_3[] = {3.8640, 2.0992, 1.1558, 0.6913, 0.4508, 0.3337, 0.2533, 0.2023, 0.1552, 0.0294};
    // // 1.5 to 2.0
    // double quark_values_4[] = {3.8815, 2.0337, 1.0959, 0.6561, 0.4460, 0.3279, 0.2598, 0.2200, 0.1684, 0.0329};
    // double gluon_values_4[] = {2.1704, 2.1106, 1.5487, 1.1265, 0.8703, 0.6887, 0.5681, 0.4858, 0.3826, 0.0752};
    // double gluon_quark_values_4[] = {3.0415, 2.1336, 1.3646, 0.9098, 0.6536, 0.4986, 0.3986, 0.3318, 0.2575, 0.0496};
    // double combined_values_4[] = {3.2022, 2.0949, 1.2986, 0.8540, 0.6141, 0.4680, 0.3760, 0.3163, 0.2455, 0.0477};
    // double resolved_values_4[] = {2.9848, 2.0959, 1.3449, 0.9023, 0.6597, 0.5099, 0.4139, 0.3521, 0.2767, 0.0538};
    // double direct_values_4[] = {3.8268, 2.0920, 1.1654, 0.7155, 0.4830, 0.3475, 0.2673, 0.2134, 0.1561, 0.0303};
    
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
    const char* titles[4] = {"-1 < #eta_{jet} < 0", "0 < #eta_{jet} < 1", 
                            "1 < #eta_{jet} < 1.5", "1.5 < #eta_{jet} < 2"};
    
    // double* combinedData[4] = {combined_values, combined_values_2, combined_values_3, combined_values_4};
    double* gluonquarkData[4] = {gluon_quark_values, gluon_quark_values_2, gluon_quark_values_3, gluon_quark_values_4};
    double* quarkData[4] = {quark_values, quark_values_2, quark_values_3, quark_values_4};
    double* gluonData[4] = {gluon_values, gluon_values_2, gluon_values_3, gluon_values_4};
    
    // Store histograms for ROOT file
    TH1F *hGQ[4], *hQuark[4], *hGluon[4]; //*hCombined[4]
    
    // Process each subplot
    for (int i = 0; i < 4; i++) {
        pads[i]->cd();
        // No grid
        pads[i]->SetGridx(0);
        pads[i]->SetGridy(0);
        pads[i]->SetFrameBorderSize(0);
        
        // Create histograms
        // hCombined[i] = new TH1F(Form("hCombined_%d", i), "", nPoints, rBins);
        hGQ[i] = new TH1F(Form("hGQ_%d", i), "", nPoints, rBins);
        hQuark[i] = new TH1F(Form("hQuark_%d", i), "", nPoints, rBins);
        hGluon[i] = new TH1F(Form("hGluon_%d", i), "", nPoints, rBins);
        
        // Fill histograms
        for (int j = 0; j < nPoints; j++) {
            // hCombined[i]->SetBinContent(j+1, combinedData[i][j]);
            hGQ[i]->SetBinContent(j+1, gluonquarkData[i][j]);
            hQuark[i]->SetBinContent(j+1, quarkData[i][j]);
            hGluon[i]->SetBinContent(j+1, gluonData[i][j]);
        }

        // Style for Combined - Blue (to match example)
        hGQ[i]->SetLineWidth(2);
        // hGQ[i]->SetLineColor(kBlue);
        // hGQ[i]->SetLineColor(TColor::GetColor("#74B9FF"));  
        
        // Style for Quark - Orange/Red (to match example)
        hQuark[i]->SetLineWidth(2);
        // hQuark[i]->SetLineColor(kOrange+7);
        // hQuark[i]->SetLineColor(TColor::GetColor("#FFEAA7")); 
        
        // Style for Gluon - Green (to match example)
        hGluon[i]->SetLineWidth(2);
        // hGluon[i]->SetLineColor(kGreen+2);
        // hGluon[i]->SetLineColor(TColor::GetColor("#00B894"));

        hGQ[i]->SetLineColor(TColor::GetColor("#264653"));   // Dark Teal
        hQuark[i]->SetLineColor(TColor::GetColor("#2A9D8F")); // Teal
        hGluon[i]->SetLineColor(TColor::GetColor("#E76F51")); // Terracotta

        // hGQ[i]->SetLineColor(TColor::GetColor("#4A5568"));   // Gray Blue
        // hQuark[i]->SetLineColor(TColor::GetColor("#E53E3E")); // Red
        // hGluon[i]->SetLineColor(TColor::GetColor("#38A169")); // Green

        hGQ[i]->SetMinimum(0);
        hGQ[i]->SetMaximum(5.5);
        hGQ[i]->GetXaxis()->SetRangeUser(0, 1.0);
        
        // Draw histograms - using HIST for step-like appearance as in example
        // hCombined[i]->Draw("HIST");
        hGQ[i]->Draw("HIST");
        hQuark[i]->Draw("HIST SAME");
        hGluon[i]->Draw("HIST SAME");
        
        // Add energy and eta range labels
        TLatex *t = new TLatex();
        // t->SetTextFont(42);
        t->SetTextFont(132);
        t->SetTextSize(0.05);
        t->SetTextAlign(31);
        // t->DrawLatex(0.90, 1.8, Form("E_{T}^{jet} > 10 GeV"));
        // t->DrawLatex(0.90, 1.5, titles[i]);
        
        // Add legend per panel
        if (i == 0){
            t->DrawLatex(0.97, 5.11, "ep (64 GeV)");
            t->DrawLatex(0.97, 1.8, Form("E_{T}^{jet} > 10 GeV"));
            t->DrawLatex(0.97, 1.5, Form("-4 < #eta_{jet} < 4"));
            hGQ[i]->GetYaxis()->ChangeLabel(1, -1, 0);  // Hide the "0" label
        }
        if (i == 1) {
            t->DrawLatex(0.97, 5.11, "ep (105 GeV)");
            t->DrawLatex(0.97, 1.8, Form("E_{T}^{jet} > 10 GeV"));
            t->DrawLatex(0.97, 1.5, Form("-4 < #eta_{jet} < 4"));
            TLegend *leg = new TLegend(0.58, 0.56, 0.98, 0.76);
            leg->AddEntry(hQuark[i], "Quark", "l");
            leg->AddEntry(hGluon[i], "Gluon", "l");
            leg->AddEntry(hGQ[i], "Quark+Gluon", "l");
            // leg->AddEntry(hCombined[i], "Res+Dir", "l");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            // leg->SetTextSize(0.05);
            // leg->SetTextFont(42);
            leg->SetTextSize(0.05);
            leg->SetTextFont(132);    // Change to serif
            leg->SetFillColor(kWhite);
            leg->SetFillStyle(1001);  // Solid fill
            leg->SetLineColor(kBlack);
            leg->SetLineWidth(1);
            leg->Draw();
        }
        if (i == 2){
            t->DrawLatex(0.97, 5.11, "ep (141 GeV)");
            t->DrawLatex(0.97, 1.8, Form("E_{T}^{jet} > 10 GeV"));
            t->DrawLatex(0.97, 1.5, Form("-4 < #eta_{jet} < 4"));
            hGQ[i]->GetXaxis()->SetRangeUser(0, 1);
            hGQ[i]->GetXaxis()->SetNdivisions(505, kFALSE); 
            hGQ[i]->GetXaxis()->ChangeLabel(6, -1, 0); 
        }
        if (i == 3){
            t->DrawLatex(0.97, 5.11, "ep (300 GeV)");
            t->DrawLatex(0.97, 1.8, Form("E_{T}^{jet} > 17 GeV"));
            t->DrawLatex(0.97, 1.5, Form("-4 < #eta_{jet} < 4"));
            hGQ[i]->GetXaxis()->SetRangeUser(0, 1);
            hGQ[i]->GetXaxis()->SetNdivisions(505, kFALSE); 
            hGQ[i]->GetXaxis()->ChangeLabel(1, -1, 0); 
        }
    }
    // Add shared axis labels on the canvas
    c->cd();

    // Vertical label for y-axis
    TLatex yLabel;
    yLabel.SetTextFont(132);
    yLabel.SetTextSize(0.03);
    yLabel.SetTextAlign(132); // Center alignment
    yLabel.SetTextAngle(90);
    yLabel.DrawLatexNDC(0.02, 0.50, "#rho");

    // Horizontal label for x-axis
    TLatex xLabel;
    xLabel.SetTextFont(132);
    xLabel.SetTextSize(0.03);
    xLabel.SetTextAlign(132); // Center alignment
    xLabel.DrawLatexNDC(0.50, 0.02, "r");

    // Save to file
    // c->SaveAs("differential_jet_shapes_eic64.pdf");
    // c->SaveAs("differential_jet_shapes_hera.png");
    c->SaveAs("diff_jet_shapes_sub_eta_all.pdf");
    
    // Save to ROOT file
    // TFile *outputFile = new TFile("differential_jet_shapes_hera.root", "RECREATE");
    // c->Write("canvas");
    
    // for (int i = 0; i < 4; i++) {
    //     // hCombined[i]->Write();
    //     hGQ[i]->Write();
    //     hQuark[i]->Write();
    //     hGluon[i]->Write();
    // }
    
    // outputFile->Close();
    
    std::cout << "PDF Saved" << std::endl;
}