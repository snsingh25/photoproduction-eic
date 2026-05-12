// =============================================================================
//            Integrated Jet Shapes Analysis with PDF Plotting
// =============================================================================
// Modified from jetrecoint.cc to create PDF plots for different eta bins

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLine.h"


// =============================================================================
// DUAL OUTPUT CLASS FOR CONSOLE AND FILE LOGGING
// =============================================================================
class DualOutput {
private:
    std::ofstream logFile;
    bool fileOpen;
    
public:
    DualOutput(const std::string& filename) : fileOpen(false) {
        logFile.open(filename.c_str());
        fileOpen = logFile.is_open();
        if (!fileOpen) {
            std::cerr << "Warning: Could not open log file: " << filename << std::endl;
        }
    }
    
    ~DualOutput() {
        if (fileOpen) {
            logFile.close();
        }
    }
    
    template<typename T>
    DualOutput& operator<<(const T& data) {
        std::cout << data;  // Output to console
        if (fileOpen) {
            logFile << data;  // Output to file
            logFile.flush();  // Ensure immediate write
        }
        return *this;
    }
    
    // Handle special cases like std::endl, std::flush
    DualOutput& operator<<(std::ostream& (*manip)(std::ostream&)) {
        std::cout << manip;
        if (fileOpen) {
            logFile << manip;
            logFile.flush();
        }
        return *this;
    }
};

// =============================================================================
// ANALYSIS CONFIGURATION
// =============================================================================
const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt5GeV/eic141_pt5/eic141_pt5.root";

// Analysis mode selection
const bool DIJET_ONLY = false;
const double R = 1.0;             
const double jet_etMin = 10.0;    

// Define four eta bins
struct EtaBin {
    double etaMin, etaMax;
    string label;
    string filename_suffix;
};

vector<EtaBin> etaBins = {
    {-1.0, 0.0, "-1 < #eta < 0", "eta_minus1_0"},
    {0.0, 1.0, "0 < #eta < 1", "eta_0_1"},
    {1.0, 1.5, "1 < #eta < 1.5", "eta_1_1p5"},
    {1.5, 2.0, "1.5 < #eta < 2", "eta_1p5_2"}
};

// Integrated jet shape parameters
const double r_min = 0.1;         
const double r_max = 1.0;         
const double r_step = 0.1;        

// Progress reporting
const int report_every = 100000;

// Your provided data values for the 4 different PDFs
// ETA BIN 1: (-1 < eta < 0)
double qq_64_eta1[10] = {2.7067, 2.1164, 1.3982, 0.8514, 0.5688, 0.396, 0.3346, 0.2481, 0.1452, 0.0181};
double qq_105_eta1[10] = {3.2144, 2.4553, 1.5482, 0.989, 0.6732, 0.478, 0.3704, 0.2872, 0.1802, 0.0227};
double qq_141_eta1[10] = {3.1969, 2.3892, 1.4322, 0.8947, 0.6225, 0.4458, 0.3458, 0.2681, 0.1735, 0.0244};
double qq_300_eta1[10] = {3.1127, 2.2566, 1.3594, 0.871, 0.6213, 0.4571, 0.3672, 0.3044, 0.2078, 0.0335};

double gg_64_eta1[10] = {1.2393, 1.7212, 1.4479, 1.1181, 0.9598, 0.9094, 0.5361, 0.6709, 0.2544, 0.0162};
double gg_105_eta1[10] = {1.1629, 1.7121, 1.5188, 1.2898, 0.9642, 0.7612, 0.7216, 0.5244, 0.2693, 0.0227};
double gg_141_eta1[10] = {1.2213, 1.5389, 1.477, 1.2965, 1.0624, 0.8243, 0.6386, 0.5299, 0.3096, 0.0457};
double gg_300_eta1[10] = {1.4184, 1.7653, 1.5791, 1.299, 1.0366, 0.8584, 0.7193, 0.5869, 0.3915, 0.064};

// ETA BIN 2: (0 < eta < 1)
double qq_64_eta2[10] = {3.2144, 2.4553, 1.5482, 0.989, 0.6732, 0.478, 0.3704, 0.2872, 0.1802, 0.0227};
double qq_105_eta2[10] = {3.0471, 2.4194, 1.4517, 0.9281, 0.6276, 0.4586, 0.3492, 0.2701, 0.1728, 0.0229};
double qq_141_eta2[10] = {3.0451, 2.3809, 1.4758, 0.9425, 0.673, 0.4836, 0.3814, 0.3, 0.1988, 0.0283};
double qq_300_eta2[10] = {3.0955, 2.3324, 1.4234, 0.8844, 0.6038, 0.447, 0.3407, 0.2733, 0.1797, 0.0235};

double gg_64_eta2[10] = {1.4013, 1.8083, 1.6386, 1.2993, 1.0333, 0.822, 0.6873, 0.5343, 0.3137, 0.0494};
double gg_105_eta2[10] = {1.3564, 1.8137, 1.5139, 1.3047, 1.0164, 0.8221, 0.6815, 0.5232, 0.3261, 0.0516};
double gg_141_eta2[10] = {1.2158, 1.5969, 1.5567, 1.3275, 1.0829, 0.9092, 0.7596, 0.6151, 0.3974, 0.0674};
double gg_300_eta2[10] = {1.4663, 1.7894, 1.6186, 1.3163, 1.057, 0.8749, 0.7399, 0.6255, 0.4401, 0.0749};

// ETA BIN 3: (1 < eta < 1.5)
double qq_64_eta3[10] = {3.1969, 2.3892, 1.4322, 0.8947, 0.6225, 0.4458, 0.3458, 0.2681, 0.1735, 0.0244};
double qq_105_eta3[10] = {3.1965, 2.3338, 1.4193, 0.8811, 0.6035, 0.4455, 0.348, 0.2743, 0.179, 0.026};
double qq_141_eta3[10] = {3.074, 2.2656, 1.4013, 0.9164, 0.6398, 0.4709, 0.3716, 0.3072, 0.205, 0.0308};
double qq_300_eta3[10] = {3.0471, 2.4194, 1.3875, 0.9281, 0.6276, 0.4586, 0.3492, 0.2701, 0.1728, 0.0229};

double gg_64_eta3[10] = {1.5157, 1.9342, 1.7226, 1.3443, 1.1092, 0.8719, 0.7479, 0.5836, 0.3877, 0.0596};
double gg_105_eta3[10] = {1.4394, 1.8063, 1.6565, 1.3035, 1.0745, 0.8891, 0.7483, 0.5888, 0.4094, 0.0607};
double gg_141_eta3[10] = {1.3305, 1.67, 1.5658, 1.3375, 1.1313, 0.9602, 0.7895, 0.6468, 0.4453, 0.0725};
double gg_300_eta3[10] = {1.499, 1.8235, 1.6476, 1.3386, 1.0851, 0.9005, 0.7639, 0.6471, 0.4643, 0.0799};

// ETA BIN 4: (1.5 < eta < 2)
double qq_64_eta4[10] = {3.0729, 2.2427, 1.3847, 0.8767, 0.6124, 0.4351, 0.328, 0.2706, 0.1829, 0.0241};
double qq_105_eta4[10] = {3.2012, 2.3082, 1.4196, 0.9112, 0.6204, 0.46, 0.3616, 0.2799, 0.1963, 0.0282};
double qq_141_eta4[10] = {3.0965, 2.2534, 1.4213, 0.8935, 0.6349, 0.4695, 0.3896, 0.2980, 0.2164, 0.0312};
double qq_300_eta4[10] = {3.1127, 2.2566, 1.3594, 0.871, 0.6213, 0.4571, 0.357, 0.3044, 0.2078, 0.0335};

double gg_64_eta4[10] = {1.4914, 1.7605, 1.5187, 1.4257, 1.1041, 0.8632, 0.7323, 0.5641, 0.3963, 0.0613};
double gg_105_eta4[10] = {1.4394, 1.8063, 1.6565, 1.3035, 1.0745, 0.8891, 0.7483, 0.5888, 0.4094, 0.0607};
double gg_141_eta4[10] = {1.3708, 1.7405, 1.5921, 1.3568, 1.168, 0.9773, 0.8288, 0.676, 0.4751, 0.0811};
double gg_300_eta4[10] = {1.4022, 1.7214, 1.5166, 1.225, 0.9666, 0.775, 0.6535, 0.5014, 0.328, 0.0502};

// =============================================================================
// PLOTTING FUNCTION
// =============================================================================
void createPlotForEtaBin(int etaBinIndex, const string& outputPrefix) {
    
    // r values
    double r_values[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    
    // Get data pointers based on eta bin
    double* qq_64_data, *qq_105_data, *qq_141_data, *qq_300_data;
    double* gg_64_data, *gg_105_data, *gg_141_data, *gg_300_data;
    
    switch(etaBinIndex) {
        case 0:  // eta bin 1
            qq_64_data = qq_64_eta1; qq_105_data = qq_105_eta1; 
            qq_141_data = qq_141_eta1; qq_300_data = qq_300_eta1;
            gg_64_data = gg_64_eta1; gg_105_data = gg_105_eta1; 
            gg_141_data = gg_141_eta1; gg_300_data = gg_300_eta1;
            break;
        case 1:  // eta bin 2
            qq_64_data = qq_64_eta2; qq_105_data = qq_105_eta2; 
            qq_141_data = qq_141_eta2; qq_300_data = qq_300_eta2;
            gg_64_data = gg_64_eta2; gg_105_data = gg_105_eta2; 
            gg_141_data = gg_141_eta2; gg_300_data = gg_300_eta2;
            break;
        case 2:  // eta bin 3
            qq_64_data = qq_64_eta3; qq_105_data = qq_105_eta3; 
            qq_141_data = qq_141_eta3; qq_300_data = qq_300_eta3;
            gg_64_data = gg_64_eta3; gg_105_data = gg_105_eta3; 
            gg_141_data = gg_141_eta3; gg_300_data = gg_300_eta3;
            break;
        case 3:  // eta bin 4
            qq_64_data = qq_64_eta4; qq_105_data = qq_105_eta4; 
            qq_141_data = qq_141_eta4; qq_300_data = qq_300_eta4;
            gg_64_data = gg_64_eta4; gg_105_data = gg_105_eta4; 
            gg_141_data = gg_141_eta4; gg_300_data = gg_300_eta4;
            break;
    }
    
    // Create canvas
    TCanvas* c = new TCanvas("c", "Integrated Jet Shapes", 800, 800);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    
    // Create TGraphs for each energy and process type
    TGraph* g_qq_64 = new TGraph(10, r_values, qq_64_data);
    TGraph* g_qq_105 = new TGraph(10, r_values, qq_105_data);
    TGraph* g_qq_141 = new TGraph(10, r_values, qq_141_data);
    TGraph* g_qq_300 = new TGraph(10, r_values, qq_300_data);
    
    TGraph* g_gg_64 = new TGraph(10, r_values, gg_64_data);
    TGraph* g_gg_105 = new TGraph(10, r_values, gg_105_data);
    TGraph* g_gg_141 = new TGraph(10, r_values, gg_141_data);
    TGraph* g_gg_300 = new TGraph(10, r_values, gg_300_data);
    
    // Set colors and styles for QQ (red family)
    g_qq_64->SetLineColor(kRed-7);     // Light red
    g_qq_105->SetLineColor(kRed);      // Red  
    g_qq_141->SetLineColor(kRed+1);    // Dark red
    g_qq_300->SetLineColor(kRed+2);    // Darker red
    
    // Set colors and styles for GG (blue family)
    g_gg_64->SetLineColor(kBlue-7);    // Light blue
    g_gg_105->SetLineColor(kBlue);     // Blue
    g_gg_141->SetLineColor(kBlue+1);   // Dark blue
    g_gg_300->SetLineColor(kBlue+2);   // Darker blue
    
    // Set line styles
    g_qq_64->SetLineStyle(3);   // Dotted
    g_qq_105->SetLineStyle(2);  // Dashed
    g_qq_141->SetLineStyle(4);  // Dash-dot
    g_qq_300->SetLineStyle(1);  // Solid
    
    g_gg_64->SetLineStyle(3);   // Dotted
    g_gg_105->SetLineStyle(2);  // Dashed  
    g_gg_141->SetLineStyle(4);  // Dash-dot
    g_gg_300->SetLineStyle(1);  // Solid
    
    // Set line width
    g_qq_64->SetLineWidth(2); g_qq_105->SetLineWidth(2);
    g_qq_141->SetLineWidth(2); g_qq_300->SetLineWidth(2);
    g_gg_64->SetLineWidth(2); g_gg_105->SetLineWidth(2);
    g_gg_141->SetLineWidth(2); g_gg_300->SetLineWidth(2);
    
    // Draw first graph to set up axes
    g_qq_300->SetTitle("");
    g_qq_300->GetXaxis()->SetTitle("r");
    g_qq_300->GetYaxis()->SetTitle("#Psi^{int}(r)");
    g_qq_300->GetXaxis()->SetRangeUser(0.0, 1.1);
    g_qq_300->GetYaxis()->SetRangeUser(0.0, 4.0);
    g_qq_300->GetXaxis()->SetTitleSize(0.05);
    g_qq_300->GetYaxis()->SetTitleSize(0.05);
    g_qq_300->GetXaxis()->SetLabelSize(0.04);
    g_qq_300->GetYaxis()->SetLabelSize(0.04);
    g_qq_300->Draw("AL");
    
    // Draw all graphs
    g_qq_64->Draw("L SAME"); g_qq_105->Draw("L SAME"); g_qq_141->Draw("L SAME");
    g_gg_64->Draw("L SAME"); g_gg_105->Draw("L SAME"); 
    g_gg_141->Draw("L SAME"); g_gg_300->Draw("L SAME");
    
    // Create legend
    TLegend* leg = new TLegend(0.65, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(g_qq_64, "64 GeV", "l");
    leg->AddEntry(g_qq_105, "105 GeV", "l");
    leg->AddEntry(g_qq_141, "141 GeV", "l");
    leg->AddEntry(g_qq_300, "300 GeV", "l");
    leg->Draw();
    
    // Add QQ/GG labels
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.05);
    latex->SetTextColor(kRed);
    latex->SetTextFont(62);
    latex->DrawLatex(0.5, 3.5, "QQ");
    
    latex->SetTextColor(kBlue);
    latex->DrawLatex(0.6, 3.5, "/ GG");
    
    // Add eta bin label
    latex->SetTextColor(kBlack);
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.1, 3.7, etaBins[etaBinIndex].label.c_str());
    
    // Save plot
    string filename = outputPrefix + etaBins[etaBinIndex].filename_suffix + ".pdf";
    c->SaveAs(filename.c_str());
    
    cout << "Created plot: " << filename << endl;
    
    delete c; delete leg; delete latex;
    delete g_qq_64; delete g_qq_105; delete g_qq_141; delete g_qq_300;
    delete g_gg_64; delete g_gg_105; delete g_gg_141; delete g_gg_300;
}

// =============================================================================
// MAIN FUNCTION
// =============================================================================
int main() {
    
    cout << "=============================================================================\n";
    cout << "INTEGRATED JET SHAPES PLOTTING FOR DIFFERENT ETA BINS\n";
    cout << "=============================================================================\n";
    
    // Initialize ROOT
    gROOT->SetBatch(kTRUE);  // Run in batch mode (no X11)
    gStyle->SetOptStat(0);   // Turn off statistics box
    
    // Create plots for first two eta bins (first PDF)
    cout << "\nCreating first PDF with eta bins 1 and 2...\n";
    createPlotForEtaBin(0, "qq_gg_integrated_jetshapes_");  // eta bin 1: -1 < eta < 0
    createPlotForEtaBin(1, "qq_gg_integrated_jetshapes_");  // eta bin 2: 0 < eta < 1
    
    // Create plots for last two eta bins (second PDF) 
    cout << "\nCreating second PDF with eta bins 3 and 4...\n";
    createPlotForEtaBin(2, "qq_gg_integrated_jetshapes_");  // eta bin 3: 1 < eta < 1.5
    createPlotForEtaBin(3, "qq_gg_integrated_jetshapes_");  // eta bin 4: 1.5 < eta < 2
    
    cout << "\n=============================================================================\n";
    cout << "All plots created successfully!\n";
    cout << "Files saved:\n";
    cout << "  - qq_gg_integrated_jetshapes_eta_minus1_0.pdf\n";
    cout << "  - qq_gg_integrated_jetshapes_eta_0_1.pdf\n";
    cout << "  - qq_gg_integrated_jetshapes_eta_1_1p5.pdf\n";
    cout << "  - qq_gg_integrated_jetshapes_eta_1p5_2.pdf\n";
    cout << "=============================================================================\n";
    
    return 0;
}