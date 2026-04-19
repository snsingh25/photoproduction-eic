// =============================================================================
// Enhanced Jet Shape Analysis Script - Inspired by Siddharth's Approach
// =============================================================================
// Description: Reads photoproduction data, performs jet clustering, and 
//              calculates differential and integrated jet shapes following
//              the validated approach from original analysis
//
// Input: photoproduction_data.root (from robust event generator)
// Output: Single ROOT file with jet shape data for physics analysis
//
// Authors: Siddharth Singh, Enhanced for EIC studies
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <iomanip>
#include <map>
#include <cmath>
#include <numeric>

using namespace fastjet;
using namespace std;

// =============================================================================
// CONFIGURABLE ANALYSIS PARAMETERS
// =============================================================================

struct AnalysisConfig {
    // Jet clustering
    double jetRadius = 1.0;
    
    // Kinematic cuts (easily changeable)
    double etMin = 17.0;        // Minimum ET cut (GeV) - HERA default
    double etaMin = -4.0;       // Minimum eta 
    double etaMax = 4.0;        // Maximum eta
    
    // Jet shape calculation parameters
    vector<double> radii = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};  // Radial points
    double deltaR = 0.1;        // Radial bin width for differential shapes
    
    // Specific cuts for quark/gluon discrimination (from your analysis)
    double thinJetCut = 0.8;    // ψ(0.3) > 0.8 for thin (quark) jets
    double thickJetCut = 0.6;   // ψ(0.3) < 0.6 for thick (gluon) jets
    double discriminatorRadius = 0.3;  // Radius for q/g discrimination
    
    // Additional eta ranges for studies
    vector<pair<double, double>> etaRanges = {
        {-1.0, 0.0},   // Forward region
        {0.0, 1.0},    // Central region  
        {1.0, 1.5},    // Transition region
        {1.5, 2.0}     // Far forward region
    };
    
    void print() const {
        cout << "Analysis Configuration:\n";
        cout << "  Jet Radius:       " << jetRadius << "\n";
        cout << "  ET Cut:           " << etMin << " GeV\n";
        cout << "  Eta Range:        [" << etaMin << ", " << etaMax << "]\n";
        cout << "  Radial Points:    " << radii.size() << "\n";
        cout << "  Delta R:          " << deltaR << "\n";
        cout << "  Q/G Cuts:         thin > " << thinJetCut << ", thick < " << thickJetCut << "\n";
    }
};

// =============================================================================
// JET SHAPE CALCULATION FUNCTIONS (Based on Siddharth's validated approach)
// =============================================================================

// Structure to hold jet shape results
struct JetShapeResult {
    vector<double> radii;
    vector<double> differentialShape;  // ρ(r)
    vector<double> integratedShape;    // ψ(r)
    vector<double> errors;             // Statistical errors
    int nJets;
    double avgET;
    double avgEta;
    int nThinJets;  // Number of quark-like jets
    int nThickJets; // Number of gluon-like jets
    
    JetShapeResult(int nBins) : radii(nBins), differentialShape(nBins), 
                               integratedShape(nBins), errors(nBins), 
                               nJets(0), avgET(0), avgEta(0), nThinJets(0), nThickJets(0) {}
};

// Function to calculate distance in eta-phi space (following your method)
double deltaR_distance(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    
    // Handle phi wraparound
    while (dphi > M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    
    return sqrt(deta*deta + dphi*dphi);
}

// Function to calculate average (utility function like yours)
double average(const vector<double>& v) {
    if (v.empty()) return 0.0;
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================

int main() {
    
    try {
        
        // =====================================================================
        // CONFIGURATION
        // =====================================================================
        
        AnalysisConfig config;
        
        // Easy parameter modification (following your style)
        config.etMin = 17.0;      // HERA energy cut (change to 10.0 for EIC)
        config.etaMin = -4.0;
        config.etaMax = 4.0;
        
        const string inputFileName = "/Users/siddharthsingh/Analysis/ph-new/evt/photoproduction_data.root";
        const string outputFileName = "photoproduction_jetshapes.root";
        
        cout << "=============================================================================\n";
        cout << "PHOTOPRODUCTION JET SHAPE ANALYSIS (Following Siddharth's Method)\n";
        cout << "=============================================================================\n";
        cout << "Input File:       " << inputFileName << "\n";
        cout << "Output File:      " << outputFileName << "\n";
        config.print();
        cout << "=============================================================================\n\n";

        // =====================================================================
        // OPEN INPUT FILE
        // =====================================================================
        
        unique_ptr<TFile> inputFile(TFile::Open(inputFileName.c_str(), "READ"));
        if (!inputFile || inputFile->IsZombie()) {
            throw runtime_error("Could not open input file: " + inputFileName);
        }
        
        cout << "Successfully opened input file.\n";

        // =====================================================================
        // CREATE OUTPUT FILE WITH DIRECTORY STRUCTURE
        // =====================================================================
        
        unique_ptr<TFile> outputFile(new TFile(outputFileName.c_str(), "RECREATE"));
        if (!outputFile || outputFile->IsZombie()) {
            throw runtime_error("Failed to create output file: " + outputFileName);
        }
        
        // Create directories (following your structure)
        TDirectory* dirQQ = outputFile->mkdir("QQ_JetShapes");
        TDirectory* dirGG = outputFile->mkdir("GG_JetShapes");
        TDirectory* dirGQ = outputFile->mkdir("GQ_JetShapes");
        TDirectory* dirCombined = outputFile->mkdir("Combined_JetShapes");
        TDirectory* dirResolved = outputFile->mkdir("Resolved_JetShapes");
        TDirectory* dirDirect = outputFile->mkdir("Direct_JetShapes");
        TDirectory* dirThin = outputFile->mkdir("Thin_JetShapes");     // Quark-like
        TDirectory* dirThick = outputFile->mkdir("Thick_JetShapes");   // Gluon-like
        TDirectory* dirConfig = outputFile->mkdir("Analysis_Config");

        // =====================================================================
        // CREATE TREES FOR JET SHAPE DATA
        // =====================================================================
        
        struct JetShapeData {
            vector<Float_t> radii;
            vector<Float_t> differentialShape;  // ρ(r)
            vector<Float_t> integratedShape;    // ψ(r)
            vector<Float_t> errors;
            Int_t nJets;
            Int_t nThinJets;   // Quark-like jets
            Int_t nThickJets;  // Gluon-like jets
            Float_t avgET;
            Float_t avgEta;
            TString etaRange;
            Float_t discriminatorValue;  // ψ(0.3) for q/g separation
        };
        
        auto createJetShapeTree = [](TDirectory* dir, const string& name) -> pair<TTree*, JetShapeData*> {
            dir->cd();
            TTree* tree = new TTree(name.c_str(), ("Jet Shapes for " + name).c_str());
            JetShapeData* data = new JetShapeData();
            
            tree->Branch("radii", &data->radii);
            tree->Branch("differentialShape", &data->differentialShape);
            tree->Branch("integratedShape", &data->integratedShape);
            tree->Branch("errors", &data->errors);
            tree->Branch("nJets", &data->nJets);
            tree->Branch("nThinJets", &data->nThinJets);
            tree->Branch("nThickJets", &data->nThickJets);
            tree->Branch("avgET", &data->avgET);
            tree->Branch("avgEta", &data->avgEta);
            tree->Branch("etaRange", &data->etaRange);
            tree->Branch("discriminatorValue", &data->discriminatorValue);
            
            return make_pair(tree, data);
        };
        
        // Create trees using traditional approach (avoiding C++17 structured bindings)
        pair<TTree*, JetShapeData*> qqPair = createJetShapeTree(dirQQ, "QQ_JetShapes");
        TTree* treeQQ = qqPair.first;
        JetShapeData* dataQQ = qqPair.second;
        
        pair<TTree*, JetShapeData*> ggPair = createJetShapeTree(dirGG, "GG_JetShapes");
        TTree* treeGG = ggPair.first;
        JetShapeData* dataGG = ggPair.second;
        
        pair<TTree*, JetShapeData*> gqPair = createJetShapeTree(dirGQ, "GQ_JetShapes");
        TTree* treeGQ = gqPair.first;
        JetShapeData* dataGQ = gqPair.second;
        
        pair<TTree*, JetShapeData*> combinedPair = createJetShapeTree(dirCombined, "Combined_JetShapes");
        TTree* treeCombined = combinedPair.first;
        JetShapeData* dataCombined = combinedPair.second;
        
        pair<TTree*, JetShapeData*> resolvedPair = createJetShapeTree(dirResolved, "Resolved_JetShapes");
        TTree* treeResolved = resolvedPair.first;
        JetShapeData* dataResolved = resolvedPair.second;
        
        pair<TTree*, JetShapeData*> directPair = createJetShapeTree(dirDirect, "Direct_JetShapes");
        TTree* treeDirect = directPair.first;
        JetShapeData* dataDirect = directPair.second;
        
        pair<TTree*, JetShapeData*> thinPair = createJetShapeTree(dirThin, "Thin_JetShapes");
        TTree* treeThin = thinPair.first;
        JetShapeData* dataThin = thinPair.second;
        
        pair<TTree*, JetShapeData*> thickPair = createJetShapeTree(dirThick, "Thick_JetShapes");
        TTree* treeThick = thickPair.first;
        JetShapeData* dataThick = thickPair.second;

        // =====================================================================
        // READ EVENTS AND COLLECT JETS BY TYPE
        // =====================================================================
        
        cout << "Reading events and clustering jets...\n";
        
        TDirectory* combinedDir = inputFile->GetDirectory("Combined");
        if (!combinedDir) {
            throw runtime_error("Could not find Combined directory in input file");
        }
        
        TTree* combinedTree = (TTree*)combinedDir->Get("CombinedEvents");
        if (!combinedTree) {
            throw runtime_error("Could not find CombinedEvents tree");
        }
        
        TTreeReader reader(combinedTree);
        
        // Event-level branches
        TTreeReaderValue<Int_t> eventID(reader, "eventID");
        TTreeReaderValue<Int_t> processCode(reader, "processCode");
        TTreeReaderValue<Int_t> processType(reader, "processType");
        TTreeReaderValue<Bool_t> isResolved(reader, "isResolved");
        TTreeReaderValue<Bool_t> isDirect(reader, "isDirect");
        
        // Particle arrays
        TTreeReaderArray<Float_t> px(reader, "px");
        TTreeReaderArray<Float_t> py(reader, "py");
        TTreeReaderArray<Float_t> pz(reader, "pz");
        TTreeReaderArray<Float_t> energy(reader, "energy");
        TTreeReaderArray<Int_t> pdgId(reader, "pdgId");
        
        JetDefinition jetDef(kt_algorithm, config.jetRadius);
        Selector jetSelector = SelectorEtMin(config.etMin) && SelectorRapRange(config.etaMin, config.etaMax);
        
        int totalEvents = 0;
        int validEvents = 0;
        
        // Accumulators for jet shape calculations
        struct JetShapeAccumulator {
            vector<double> rho_sum;
            vector<double> psi_sum;
            int nJets;
            double totalET;
            double totalEta;
            int nThinJets;
            int nThickJets;
            
            JetShapeAccumulator(int nRadii) : rho_sum(nRadii, 0.0), psi_sum(nRadii, 0.0), 
                                             nJets(0), totalET(0.0), totalEta(0.0), 
                                             nThinJets(0), nThickJets(0) {}
        };
        
        JetShapeAccumulator qqAcc(config.radii.size());
        JetShapeAccumulator ggAcc(config.radii.size());
        JetShapeAccumulator gqAcc(config.radii.size());
        JetShapeAccumulator combinedAcc(config.radii.size());
        JetShapeAccumulator resolvedAcc(config.radii.size());
        JetShapeAccumulator directAcc(config.radii.size());
        JetShapeAccumulator thinAcc(config.radii.size());
        JetShapeAccumulator thickAcc(config.radii.size());
        
        // Event loop with immediate jet shape calculation
        while (reader.Next()) {
            totalEvents++;
            
            if (totalEvents % 10000 == 0) {
                cout << "Processed " << totalEvents << " events\r" << flush;
            }
            
            if (*processCode == 0) continue;
            validEvents++;
            
            // Create particles for jet clustering
            vector<PseudoJet> particles;
            for (int i = 0; i < px.GetSize(); ++i) {
                // Skip neutrinos (following your approach)
                if (abs(pdgId[i]) == 12 || abs(pdgId[i]) == 14 || abs(pdgId[i]) == 16) {
                    continue;
                }
                particles.push_back(PseudoJet(px[i], py[i], pz[i], energy[i]));
            }
            
            if (particles.size() < 2) continue;
            
            // Cluster jets
            ClusterSequence clusterSeq(particles, jetDef);
            vector<PseudoJet> jets = sorted_by_pt(jetSelector(clusterSeq.inclusive_jets()));
            
            if (jets.empty()) continue;
            
            // Calculate discriminator for quark/gluon separation inline
            double discriminator = 0.0;
            if (!jets.empty()) {
                double psi_sum = 0.0;
                for (const auto& jet : jets) {
                    if (jet.Et() < config.etMin) continue;
                    
                    double etWithinRadius = 0.0;
                    vector<PseudoJet> constituents = jet.constituents();
                    for (const auto& constituent : constituents) {
                        double dr = deltaR_distance(jet.eta(), jet.phi(), 
                                                  constituent.eta(), constituent.phi());
                        if (dr < config.discriminatorRadius) {
                            etWithinRadius += constituent.Et();
                        }
                    }
                    if (jet.Et() > 0) {
                        psi_sum += (etWithinRadius / jet.Et());
                    }
                }
                discriminator = psi_sum / jets.size();
            }
            
            // Helper function to accumulate jet shapes
            auto accumulateJetShapes = [&](JetShapeAccumulator& acc, const vector<PseudoJet>& eventJets) {
                if (eventJets.empty()) return;
                
                for (const auto& jet : eventJets) {
                    if (jet.Et() < config.etMin || jet.eta() < config.etaMin || jet.eta() > config.etaMax) continue;
                    
                    acc.nJets++;
                    acc.totalET += jet.Et();
                    acc.totalEta += jet.eta();
                    
                    // Calculate jet shapes at each radius
                    for (size_t r = 0; r < config.radii.size(); ++r) {
                        double radius = config.radii[r];
                        
                        // Calculate differential shape ρ(r)
                        double R1 = radius - (config.deltaR / 2.0);
                        double R2 = radius + (config.deltaR / 2.0);
                        double etInAnnulus = 0.0;
                        
                        // Calculate integrated shape ψ(r)
                        double etWithinRadius = 0.0;
                        
                        vector<PseudoJet> constituents = jet.constituents();
                        for (const auto& constituent : constituents) {
                            double dr = deltaR_distance(jet.eta(), jet.phi(), 
                                                      constituent.eta(), constituent.phi());
                            
                            // For differential shape
                            if (dr > R1 && dr < R2) {
                                etInAnnulus += constituent.Et();
                            }
                            
                            // For integrated shape
                            if (dr < radius) {
                                etWithinRadius += constituent.Et();
                            }
                        }
                        
                        if (jet.Et() > 0) {
                            acc.rho_sum[r] += (etInAnnulus / jet.Et()) / config.deltaR;
                            acc.psi_sum[r] += (etWithinRadius / jet.Et());
                        }
                    }
                    
                    // Quark/gluon classification
                    if (discriminator > config.thinJetCut) {
                        acc.nThinJets++;
                    } else if (discriminator < config.thickJetCut) {
                        acc.nThickJets++;
                    }
                }
            };
            
            // Accumulate for all relevant categories
            accumulateJetShapes(combinedAcc, jets);
            
            if (*processType == 1) accumulateJetShapes(qqAcc, jets);
            else if (*processType == 2) accumulateJetShapes(ggAcc, jets);
            else if (*processType == 3) accumulateJetShapes(gqAcc, jets);
            
            if (*isResolved) accumulateJetShapes(resolvedAcc, jets);
            if (*isDirect) accumulateJetShapes(directAcc, jets);
            
            // Separate thin and thick jets
            vector<PseudoJet> thinJetsEvent, thickJetsEvent;
            if (discriminator > config.thinJetCut) {
                thinJetsEvent = jets;
            } else if (discriminator < config.thickJetCut) {
                thickJetsEvent = jets;
            }
            
            if (!thinJetsEvent.empty()) accumulateJetShapes(thinAcc, thinJetsEvent);
            if (!thickJetsEvent.empty()) accumulateJetShapes(thickAcc, thickJetsEvent);
        }
        
        cout << "\nCompleted event processing.\n";
        cout << "Total events: " << totalEvents << ", Valid events: " << validEvents << "\n";
        cout << "Jets collected:\n";
        cout << "  QQ Jets: " << qqAcc.nJets << "\n";
        cout << "  GG Jets: " << ggAcc.nJets << "\n";
        cout << "  GQ Jets: " << gqAcc.nJets << "\n";
        cout << "  Combined: " << combinedAcc.nJets << "\n";
        cout << "  Resolved: " << resolvedAcc.nJets << "\n";
        cout << "  Direct: " << directAcc.nJets << "\n";
        cout << "  Thin (Quark-like): " << thinAcc.nJets << "\n";
        cout << "  Thick (Gluon-like): " << thickAcc.nJets << "\n\n";

        // =====================================================================
        // CONVERT ACCUMULATORS TO FINAL RESULTS
        // =====================================================================
        
        cout << "Finalizing jet shape calculations...\n";
        
        auto convertAccumulatorToResult = [&](const JetShapeAccumulator& acc) -> JetShapeResult {
            JetShapeResult result(config.radii.size());
            
            if (acc.nJets == 0) return result;
            
            result.nJets = acc.nJets;
            result.nThinJets = acc.nThinJets;
            result.nThickJets = acc.nThickJets;
            result.avgET = acc.totalET / acc.nJets;
            result.avgEta = acc.totalEta / acc.nJets;
            
            for (size_t i = 0; i < config.radii.size(); ++i) {
                result.radii[i] = config.radii[i];
                result.differentialShape[i] = acc.rho_sum[i] / acc.nJets;
                result.integratedShape[i] = acc.psi_sum[i] / acc.nJets;
                result.errors[i] = result.differentialShape[i] / sqrt(max(1, acc.nJets));
            }
            
            return result;
        };
        
        // Convert all accumulators to results
        auto qqShapes = convertAccumulatorToResult(qqAcc);
        auto ggShapes = convertAccumulatorToResult(ggAcc);
        auto gqShapes = convertAccumulatorToResult(gqAcc);
        auto combinedShapes = convertAccumulatorToResult(combinedAcc);
        auto resolvedShapes = convertAccumulatorToResult(resolvedAcc);
        auto directShapes = convertAccumulatorToResult(directAcc);
        auto thinShapes = convertAccumulatorToResult(thinAcc);
        auto thickShapes = convertAccumulatorToResult(thickAcc);
        
        // =====================================================================
        // FILL TREES WITH JET SHAPE DATA
        // =====================================================================
        
        auto fillJetShapeData = [&](JetShapeData* data, const JetShapeResult& result, const string& etaRange) {
            data->radii.clear();
            data->differentialShape.clear();
            data->integratedShape.clear();
            data->errors.clear();
            
            for (size_t i = 0; i < result.radii.size(); ++i) {
                data->radii.push_back(result.radii[i]);
                data->differentialShape.push_back(result.differentialShape[i]);
                data->integratedShape.push_back(result.integratedShape[i]);
                data->errors.push_back(result.errors[i]);
            }
            
            data->nJets = result.nJets;
            data->nThinJets = result.nThinJets;
            data->nThickJets = result.nThickJets;
            data->avgET = result.avgET;
            data->avgEta = result.avgEta;
            data->etaRange = etaRange;
            
            // Calculate discriminator value
            if (!result.integratedShape.empty() && result.radii.size() >= 3) {
                data->discriminatorValue = result.integratedShape[2]; // ψ(0.3)
            } else {
                data->discriminatorValue = 0.0;
            }
        };
        
        // Fill trees
        fillJetShapeData(dataQQ, qqShapes, "overall");
        treeQQ->Fill();
        
        fillJetShapeData(dataGG, ggShapes, "overall");
        treeGG->Fill();
        
        fillJetShapeData(dataGQ, gqShapes, "overall");
        treeGQ->Fill();
        
        fillJetShapeData(dataCombined, combinedShapes, "overall");
        treeCombined->Fill();
        
        fillJetShapeData(dataResolved, resolvedShapes, "overall");
        treeResolved->Fill();
        
        fillJetShapeData(dataDirect, directShapes, "overall");
        treeDirect->Fill();
        
        fillJetShapeData(dataThin, thinShapes, "overall");
        treeThin->Fill();
        
        fillJetShapeData(dataThick, thickShapes, "overall");
        treeThick->Fill();

        // =====================================================================
        // SAVE ANALYSIS CONFIGURATION AND SUMMARY
        // =====================================================================
        
        dirConfig->cd();
        TTree* configTree = new TTree("AnalysisConfig", "Analysis Configuration Parameters");
        
        Double_t jetRadius = config.jetRadius;
        Double_t etMin = config.etMin;
        Double_t etaMin = config.etaMin;
        Double_t etaMax = config.etaMax;
        Double_t deltaR = config.deltaR;
        Double_t thinJetCut = config.thinJetCut;
        Double_t thickJetCut = config.thickJetCut;
        
        configTree->Branch("jetRadius", &jetRadius);
        configTree->Branch("etMin", &etMin);
        configTree->Branch("etaMin", &etaMin);
        configTree->Branch("etaMax", &etaMax);
        configTree->Branch("deltaR", &deltaR);
        configTree->Branch("thinJetCut", &thinJetCut);
        configTree->Branch("thickJetCut", &thickJetCut);
        
        configTree->Fill();
        
        // Create summary histograms
        TH1F* hJetShapeComparison = new TH1F("hJetShapeComparison", 
            "Integrated Jet Shapes at r=0.3;Process Type;#psi(0.3)", 8, 0, 8);
        hJetShapeComparison->GetXaxis()->SetBinLabel(1, "QQ");
        hJetShapeComparison->GetXaxis()->SetBinLabel(2, "GG");
        hJetShapeComparison->GetXaxis()->SetBinLabel(3, "GQ");
        hJetShapeComparison->GetXaxis()->SetBinLabel(4, "Combined");
        hJetShapeComparison->GetXaxis()->SetBinLabel(5, "Resolved");
        hJetShapeComparison->GetXaxis()->SetBinLabel(6, "Direct");
        hJetShapeComparison->GetXaxis()->SetBinLabel(7, "Thin");
        hJetShapeComparison->GetXaxis()->SetBinLabel(8, "Thick");
        
        if (qqShapes.integratedShape.size() >= 3) hJetShapeComparison->SetBinContent(1, qqShapes.integratedShape[2]);
        if (ggShapes.integratedShape.size() >= 3) hJetShapeComparison->SetBinContent(2, ggShapes.integratedShape[2]);
        if (gqShapes.integratedShape.size() >= 3) hJetShapeComparison->SetBinContent(3, gqShapes.integratedShape[2]);
        if (combinedShapes.integratedShape.size() >= 3) hJetShapeComparison->SetBinContent(4, combinedShapes.integratedShape[2]);
        if (resolvedShapes.integratedShape.size() >= 3) hJetShapeComparison->SetBinContent(5, resolvedShapes.integratedShape[2]);
        if (directShapes.integratedShape.size() >= 3) hJetShapeComparison->SetBinContent(6, directShapes.integratedShape[2]);
        if (thinShapes.integratedShape.size() >= 3) hJetShapeComparison->SetBinContent(7, thinShapes.integratedShape[2]);
        if (thickShapes.integratedShape.size() >= 3) hJetShapeComparison->SetBinContent(8, thickShapes.integratedShape[2]);

        // =====================================================================
        // WRITE OUTPUT AND SUMMARY
        // =====================================================================
        
        cout << "Writing output file...\n";
        outputFile->cd();
        outputFile->Write();
        
        cout << "\n=============================================================================\n";
        cout << "JET SHAPE ANALYSIS COMPLETED SUCCESSFULLY\n";
        cout << "=============================================================================\n";
        cout << "Output File: " << outputFileName << "\n";
        cout << "\nJet Shape Results:\n";
        cout << "  QQ_JetShapes       : " << setw(4) << treeQQ->GetEntries() << " entries (avg ψ(0.3) = " 
             << fixed << setprecision(3) << (qqShapes.integratedShape.size() >= 3 ? qqShapes.integratedShape[2] : 0.0) << ")\n";
        cout << "  GG_JetShapes       : " << setw(4) << treeGG->GetEntries() << " entries (avg ψ(0.3) = " 
             << (ggShapes.integratedShape.size() >= 3 ? ggShapes.integratedShape[2] : 0.0) << ")\n";
        cout << "  GQ_JetShapes       : " << setw(4) << treeGQ->GetEntries() << " entries (avg ψ(0.3) = " 
             << (gqShapes.integratedShape.size() >= 3 ? gqShapes.integratedShape[2] : 0.0) << ")\n";
        cout << "  Combined_JetShapes : " << setw(4) << treeCombined->GetEntries() << " entries (avg ψ(0.3) = " 
             << (combinedShapes.integratedShape.size() >= 3 ? combinedShapes.integratedShape[2] : 0.0) << ")\n";
        cout << "  Resolved_JetShapes : " << setw(4) << treeResolved->GetEntries() << " entries (avg ψ(0.3) = " 
             << (resolvedShapes.integratedShape.size() >= 3 ? resolvedShapes.integratedShape[2] : 0.0) << ")\n";
        cout << "  Direct_JetShapes   : " << setw(4) << treeDirect->GetEntries() << " entries (avg ψ(0.3) = " 
             << (directShapes.integratedShape.size() >= 3 ? directShapes.integratedShape[2] : 0.0) << ")\n";
        cout << "  Thin_JetShapes     : " << setw(4) << treeThin->GetEntries() << " entries (quark-like, ψ(0.3) = " 
             << (thinShapes.integratedShape.size() >= 3 ? thinShapes.integratedShape[2] : 0.0) << ")\n";
        cout << "  Thick_JetShapes    : " << setw(4) << treeThick->GetEntries() << " entries (gluon-like, ψ(0.3) = " 
             << (thickShapes.integratedShape.size() >= 3 ? thickShapes.integratedShape[2] : 0.0) << ")\n";
        
        cout << "\nQuark vs Gluon Discrimination Results:\n";
        cout << "  Thin Jet Cut (ψ > " << config.thinJetCut << "): " << thinAcc.nJets << " jets classified as quark-like\n";
        cout << "  Thick Jet Cut (ψ < " << config.thickJetCut << "): " << thickAcc.nJets << " jets classified as gluon-like\n";
        cout << "  Total Jets: " << combinedAcc.nJets << "\n";
        
        // Print physics-relevant comparisons
        cout << "\nPhysics Summary:\n";
        if (qqShapes.integratedShape.size() >= 3 && ggShapes.integratedShape.size() >= 3) {
            double qqPsi = qqShapes.integratedShape[2];
            double ggPsi = ggShapes.integratedShape[2];
            cout << "  QQ vs GG ψ(0.3) difference: " << abs(qqPsi - ggPsi) 
                 << " (QQ=" << qqPsi << ", GG=" << ggPsi << ")\n";
            cout << "  " << (qqPsi > ggPsi ? "QQ jets are thinner (more quark-like)" : "GG jets are thinner") << "\n";
        }
        
        if (resolvedShapes.integratedShape.size() >= 3 && directShapes.integratedShape.size() >= 3) {
            double resPsi = resolvedShapes.integratedShape[2];
            double dirPsi = directShapes.integratedShape[2];
            cout << "  Resolved vs Direct ψ(0.3) difference: " << abs(resPsi - dirPsi)
                 << " (Res=" << resPsi << ", Dir=" << dirPsi << ")\n";
        }
        
        cout << "\nData Structure:\n";
        cout << "  Each tree contains jet shapes at " << config.radii.size() << " radial points\n";
        cout << "  Differential shapes ρ(r) and integrated shapes ψ(r) are stored\n";
        cout << "  Statistical errors and jet classification info included\n";
        cout << "\nReady for jet substructure analysis and quark/gluon discrimination studies!\n";
        cout << "=============================================================================\n";
        
        // =====================================================================
        // PRINT DETAILED JET SHAPE VALUES (Similar to your output style)
        // =====================================================================
        
        cout << "\n=======================================================================================\n";
        cout << "DETAILED JET SHAPE ANALYSIS RESULTS\n";
        cout << "=======================================================================================\n";
        
        // Print differential shapes at specific radius (following your rr_diff approach)
        double r_diff = 0.6;  // Your differential shape radius
        int r_diff_index = -1;
        for (size_t i = 0; i < config.radii.size(); ++i) {
            if (abs(config.radii[i] - r_diff) < 0.01) {
                r_diff_index = i;
                break;
            }
        }
        
        if (r_diff_index >= 0) {
            cout << "Differential Shape ρ(" << r_diff << ") for ET_min(" << config.etMin << " GeV):\n";
            cout << fixed << setprecision(4);
            cout << "  Combined Jets: " << (combinedShapes.differentialShape.size() > r_diff_index ? 
                                            combinedShapes.differentialShape[r_diff_index] : 0.0) << "\n";
            cout << "  Thin Jets (Quark-like): " << (thinShapes.differentialShape.size() > r_diff_index ? 
                                                    thinShapes.differentialShape[r_diff_index] : 0.0) << "\n";
            cout << "  Thick Jets (Gluon-like): " << (thickShapes.differentialShape.size() > r_diff_index ? 
                                                     thickShapes.differentialShape[r_diff_index] : 0.0) << "\n";
        }
        
        // Print integrated shapes at discriminator radius (following your rr_int approach)
        cout << "\nIntegrated Shape ψ(" << config.discriminatorRadius << ") for ET_min(" << config.etMin << " GeV):\n";
        cout << fixed << setprecision(4);
        if (combinedShapes.integratedShape.size() >= 3) {
            cout << "  Combined Jets: " << combinedShapes.integratedShape[2] << "\n";
            cout << "  Resolved Jets: " << (resolvedShapes.integratedShape.size() >= 3 ? resolvedShapes.integratedShape[2] : 0.0) << "\n";
            cout << "  Direct Jets: " << (directShapes.integratedShape.size() >= 3 ? directShapes.integratedShape[2] : 0.0) << "\n";
            cout << "  ------------------------------------------------\n";
            cout << "  Thin Jets (Quark-like): " << (thinShapes.integratedShape.size() >= 3 ? thinShapes.integratedShape[2] : 0.0) << "\n";
            cout << "  Thick Jets (Gluon-like): " << (thickShapes.integratedShape.size() >= 3 ? thickShapes.integratedShape[2] : 0.0) << "\n";
        }
        
        // Print subprocess fractions (following your style)
        cout << "\nSubprocess Analysis:\n";
        cout << fixed << setprecision(1);
        int totalProcessJets = qqAcc.nJets + ggAcc.nJets + gqAcc.nJets;
        if (totalProcessJets > 0) {
            cout << "  QQ Jets: " << qqAcc.nJets << " | Fraction: " << (100.0 * qqAcc.nJets / totalProcessJets) << "%\n";
            cout << "  GG Jets: " << ggAcc.nJets << " | Fraction: " << (100.0 * ggAcc.nJets / totalProcessJets) << "%\n";
            cout << "  GQ Jets: " << gqAcc.nJets << " | Fraction: " << (100.0 * gqAcc.nJets / totalProcessJets) << "%\n";
        }
        
        cout << "\nQuark/Gluon Classification Results:\n";
        int totalClassifiedJets = thinAcc.nJets + thickAcc.nJets;
        if (totalClassifiedJets > 0) {
            cout << "  Thin Jets (Quark-like): " << thinAcc.nJets << " | Fraction: " << (100.0 * thinAcc.nJets / totalClassifiedJets) << "%\n";
            cout << "  Thick Jets (Gluon-like): " << thickAcc.nJets << " | Fraction: " << (100.0 * thickAcc.nJets / totalClassifiedJets) << "%\n";
            cout << "  Unclassified: " << (combinedAcc.nJets - totalClassifiedJets) << " jets\n";
        }
        
        cout << "\nPhotoproduction Process Analysis:\n";
        int totalPhotoJets = resolvedAcc.nJets + directAcc.nJets;
        if (totalPhotoJets > 0) {
            cout << "  Resolved Photoproduction: " << resolvedAcc.nJets << " | Fraction: " << (100.0 * resolvedAcc.nJets / totalPhotoJets) << "%\n";
            cout << "  Direct Photoproduction: " << directAcc.nJets << " | Fraction: " << (100.0 * directAcc.nJets / totalPhotoJets) << "%\n";
        }
        
        cout << "=======================================================================================\n";
        
        // =====================================================================
        // CLEANUP AND FINALIZATION
        // =====================================================================
        
        // Cleanup dynamically allocated data
        delete dataQQ;
        delete dataGG;
        delete dataGQ;
        delete dataCombined;
        delete dataResolved;
        delete dataDirect;
        delete dataThin;
        delete dataThick;
        
        cout << "\nJet shape analysis completed successfully!\n";
        cout << "Output saved to: " << outputFileName << "\n";
        cout << "\nNext steps:\n";
        cout << "1. Use Thin_JetShapes and Thick_JetShapes for quark vs gluon studies\n";
        cout << "2. Compare QQ_JetShapes vs GG_JetShapes for subprocess differences\n";
        cout << "3. Analyze Resolved_JetShapes vs Direct_JetShapes for photoproduction mechanisms\n";
        cout << "4. Plot differential and integrated jet shapes for physics interpretation\n";
        cout << "\nData is ready for integrated jet shape analysis and subjet multiplicity studies!\n";
        
    } catch (const exception& e) {
        cerr << "\n❌ ERROR: " << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "\n❌ UNKNOWN ERROR occurred during jet shape analysis!" << endl;
        return 1;
    }
    
    return 0;
}