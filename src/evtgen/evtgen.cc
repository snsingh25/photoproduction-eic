// =============================================================================
// Photoproduction Event Generation Script (Improved Version)
// =============================================================================
// Description: Generates photoproduction events using PYTHIA8 and stores
//              final state particles categorized by process type.
//              Output filename automatically includes key PYTHIA parameters.
//
// Output: ROOT file with particle data for different process types
// Author: Siddharth Singh
// Version: 2.0 (with automatic naming and improved diagnostics)
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TNamed.h"
#include "Pythia8/Pythia.h"

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>

using namespace Pythia8;
using namespace std;

// =============================================================================
// COLLIDER PRESETS
// =============================================================================
enum ColliderPreset {
    HERA_300,      // HERA: 820 GeV p + 27.5 GeV e = 300 GeV
    EIC_64,        // EIC:  100 GeV p + 10 GeV e   = 63.2 GeV
    EIC_105,       // EIC:  275 GeV p + 10 GeV e   = 104.9 GeV
    EIC_141,       // EIC:  275 GeV p + 18 GeV e   = 140.7 GeV
    CUSTOM         // User-defined energies
};

// =============================================================================
// ANALYSIS CONFIGURATION
// =============================================================================
struct EventConfig {
    // ----- Collider Selection -----
    ColliderPreset preset = HERA_300;  // Change this to switch colliders
    
    // ----- Event Generation -----
    int nEvents = 1000000;             // Number of events to generate
    
    // ----- Beam Configuration (auto-set by preset, or manual for CUSTOM) -----
    double protonEnergy = 820.0;       // GeV
    double electronEnergy = 27.5;      // GeV
    
    // ----- Physics Cuts -----
    double pTHatMin = 7.0;             // Minimum partonic pT (GeV)
    double Q2max = 1.0;                // Maximum Q² for photoproduction (GeV²)
    
    // ----- PYTHIA Tuning Parameters -----
    double pT0Ref = 4.0;               // MPI regularization scale (GeV) - CRITICAL for γp!
    double StringZ_aLund = 0.68;       // Lund fragmentation parameter a
    double StringZ_bLund = 0.98;       // Lund fragmentation parameter b (GeV⁻²)
    double StringPT_sigma = 0.335;     // Fragmentation pT width (GeV)
    
    // ----- Tune Selection -----
    int tunePP = 14;                   // pp tune (14 = Monash 2013)
    int tuneEE = 7;                    // ee tune (7 = Monash 2013)
    
    // ----- Output Options -----
    string outputDir = "./";           // Output directory
    string customTag = "";             // Optional custom tag for filename
    
    // Apply preset configuration
    void applyPreset() {
        switch (preset) {
            case HERA_300:
                protonEnergy = 820.0;
                electronEnergy = 27.5;
                pTHatMin = 7.0;        // For E_T^jet > 17 GeV studies
                pT0Ref = 4.0;          // Tuned for HERA γp
                break;
            case EIC_64:
                protonEnergy = 100.0;
                electronEnergy = 10.0;
                pTHatMin = 5.0;        // Lower for EIC jet studies
                pT0Ref = 3.0;
                break;
            case EIC_105:
                protonEnergy = 275.0;
                electronEnergy = 10.0;
                pTHatMin = 5.0;
                pT0Ref = 3.0;
                break;
            case EIC_141:
                protonEnergy = 275.0;
                electronEnergy = 18.0;
                pTHatMin = 5.0;
                pT0Ref = 3.0;
                break;
            case CUSTOM:
                // Use manually set values
                break;
        }
    }
    
    // Calculate center-of-mass energy
    double sqrtS() const {
        return 2.0 * sqrt(protonEnergy * electronEnergy);
    }
    
    // Get collider name string
    string getColliderName() const {
        switch (preset) {
            case HERA_300: return "hera300";
            case EIC_64:   return "eic64";
            case EIC_105:  return "eic105";
            case EIC_141:  return "eic141";
            case CUSTOM:   return "custom" + to_string((int)sqrtS());
            default:       return "unknown";
        }
    }
    
    // Generate output filename with all relevant parameters
    string generateOutputFilename() const {
        ostringstream oss;
        oss << outputDir;
        oss << getColliderName();
        oss << "_pTHat" << (int)pTHatMin;
        oss << "_pT0Ref" << fixed << setprecision(1) << pT0Ref;
        oss << "_Q2max" << (int)Q2max;
        oss << "_n" << (nEvents / 1000) << "k";
        if (!customTag.empty()) {
            oss << "_" << customTag;
        }
        oss << ".root";
        return oss.str();
    }
    
    // Generate log filename
    string generateLogFilename() const {
        string rootName = generateOutputFilename();
        return rootName.substr(0, rootName.length() - 5) + ".log";
    }
    
    // Print configuration
    void print(ostream& out = cout) const {
        out << "=============================================================================\n";
        out << "EVENT GENERATION CONFIGURATION\n";
        out << "=============================================================================\n";
        out << "\n--- Collider Setup ---\n";
        out << "  Preset:           " << getColliderName() << "\n";
        out << "  Proton Energy:    " << protonEnergy << " GeV\n";
        out << "  Electron Energy:  " << electronEnergy << " GeV\n";
        out << "  √s:               " << fixed << setprecision(1) << sqrtS() << " GeV\n";
        out << "\n--- Event Generation ---\n";
        out << "  Number of Events: " << nEvents << "\n";
        out << "\n--- Physics Cuts ---\n";
        out << "  pTHatMin:         " << pTHatMin << " GeV\n";
        out << "  Q²max:            " << Q2max << " GeV²\n";
        out << "\n--- PYTHIA Tuning ---\n";
        out << "  Tune:pp:          " << tunePP << " (Monash 2013)\n";
        out << "  Tune:ee:          " << tuneEE << " (Monash 2013)\n";
        out << "  pT0Ref:           " << pT0Ref << " GeV (MPI regularization)\n";
        out << "  StringZ:aLund:    " << StringZ_aLund << "\n";
        out << "  StringZ:bLund:    " << StringZ_bLund << " GeV⁻²\n";
        out << "  StringPT:sigma:   " << StringPT_sigma << " GeV\n";
        out << "\n--- Output ---\n";
        out << "  ROOT File:        " << generateOutputFilename() << "\n";
        out << "  Log File:         " << generateLogFilename() << "\n";
        out << "=============================================================================\n\n";
    }
};

// =============================================================================
// DUAL OUTPUT CLASS (Console + Log File)
// =============================================================================
class DualOutput {
private:
    ofstream logFile;
    bool fileOpen;
    
public:
    DualOutput(const string& filename) : fileOpen(false) {
        logFile.open(filename.c_str());
        fileOpen = logFile.is_open();
        if (!fileOpen) {
            cerr << "Warning: Could not open log file: " << filename << endl;
        }
    }
    
    ~DualOutput() {
        if (fileOpen) logFile.close();
    }
    
    template<typename T>
    DualOutput& operator<<(const T& data) {
        cout << data;
        if (fileOpen) { logFile << data; logFile.flush(); }
        return *this;
    }
    
    DualOutput& operator<<(ostream& (*manip)(ostream&)) {
        cout << manip;
        if (fileOpen) { logFile << manip; logFile.flush(); }
        return *this;
    }
};

// =============================================================================
// SUBPROCESS CLASSIFICATION
// =============================================================================
// Classification based on INITIATING parton from the PROTON side
// For jet shape studies, what matters is the color charge of the 
// parton that initiates the hard scatter

int classifySubprocess(int processCode) {
    // -----------------------------------------------------------------
    // Classification is by OUTGOING hard-process parton configuration:
    //   1 (QQ) : both outgoing hard partons are quarks  (|pdg| in 1..6)
    //   2 (GG) : both outgoing hard partons are gluons  (pdg == 21)
    //   3 (GQ) : one outgoing quark + one outgoing gluon
    // This matches the labelling convention the paper uses (and against
    // which ZEUS-era jet-shape measurements are quoted) and is the
    // convention under which the existing event samples in
    // data/allevents_pt7GeV/ were produced (audit: ≥99.99% match,
    // see scripts/audit_truth_labels.py).
    //
    // Process-code list taken from getProcessName() below; cross-
    // reference against the PYTHIA 8 manual for HardQCD and
    // PhotonParton.
    // -----------------------------------------------------------------

    // QQ: both outgoing hard partons are quarks
    if (processCode == 112 ||    // g g → q qbar          (resolved)
        processCode == 114 ||    // q q' → q q'           (resolved)
        processCode == 116 ||    // q qbar → q' qbar'     (resolved)
        processCode == 121 ||    // g g → c cbar          (resolved)
        processCode == 122 ||    // q qbar → c cbar       (resolved)
        processCode == 123 ||    // g g → b bbar          (resolved)
        processCode == 124 ||    // q qbar → b bbar       (resolved)
        processCode == 271 ||    // g γ → q qbar (uds)    (direct BGF)
        processCode == 272 ||    // g γ → c cbar          (direct BGF)
        processCode == 273 ||    // g γ → b bbar          (direct BGF)
        processCode == 281 ||    // γ g → q qbar (uds)    (direct BGF)
        processCode == 282 ||    // γ g → c cbar          (direct BGF)
        processCode == 283) {    // γ g → b bbar          (direct BGF)
        return 1;
    }

    // GG: both outgoing hard partons are gluons
    if (processCode == 111 ||    // g g → g g             (resolved)
        processCode == 115) {    // q qbar → g g          (resolved)
        return 2;
    }

    // GQ: one outgoing quark, one outgoing gluon
    if (processCode == 113 ||    // q g → q g             (resolved)
        processCode == 274 ||    // q γ → q g             (direct QCDC)
        processCode == 284) {    // γ q → q g             (direct QCDC)
        return 3;
    }

    return 0;
}

// Get process name for logging
string getProcessName(int code) {
    switch (code) {
        case 111: return "g g → g g";
        case 112: return "g g → q qbar";
        case 113: return "q g → q g";
        case 114: return "q q' → q q'";
        case 115: return "q qbar → g g";
        case 116: return "q qbar → q' qbar'";
        case 121: return "g g → c cbar";
        case 122: return "q qbar → c cbar";
        case 123: return "g g → b bbar";
        case 124: return "q qbar → b bbar";
        case 271: return "g γ → q qbar";
        case 272: return "g γ → c cbar";
        case 273: return "g γ → b bbar";
        case 274: return "q γ → q g";
        case 281: return "γ g → q qbar";
        case 282: return "γ g → c cbar";
        case 283: return "γ g → b bbar";
        case 284: return "γ q → q g";
        default:  return "Unknown (" + to_string(code) + ")";
    }
}

// =============================================================================
// MAIN FUNCTION
// =============================================================================

int main(int argc, char** argv) {

    try {

        // =====================================================================
        // CONFIGURATION
        // =====================================================================
        //   ./bin/evtgen [preset] [nEvents] [customTag]
        //
        // Defaults: HERA_300, 1,000,000 attempts, no tag (matches the
        // original production). Any arg you leave out keeps its default.

        EventConfig config;
        config.preset   = HERA_300;
        config.nEvents  = 1000000;
        config.customTag = "";

        if (argc > 1) {
            string presetArg = argv[1];
            if      (presetArg == "HERA_300") config.preset = HERA_300;
            else if (presetArg == "EIC_64")   config.preset = EIC_64;
            else if (presetArg == "EIC_105")  config.preset = EIC_105;
            else if (presetArg == "EIC_141")  config.preset = EIC_141;
            else if (presetArg == "CUSTOM")   config.preset = CUSTOM;
            else {
                cerr << "Unknown preset '" << presetArg << "'. "
                     << "Options: HERA_300, EIC_64, EIC_105, EIC_141, CUSTOM.\n";
                return 1;
            }
        }
        if (argc > 2) config.nEvents = stol(argv[2]);
        if (argc > 3) config.customTag = argv[3];

        // --- Apply preset and generate filenames ---
        config.applyPreset();
        
        const string outputFileName = config.generateOutputFilename();
        const string logFileName = config.generateLogFilename();
        
        // Create dual output (console + log file)
        DualOutput dout(logFileName);
        
        // Print timestamp
        time_t now = time(0);
        dout << "Event Generation Started: " << ctime(&now);
        
        // Print configuration
        config.print(cout);
        
        // Also write to log file directly
        ofstream configLog(logFileName, ios::app);
        config.print(configLog);
        configLog.close();

        // =====================================================================
        // PYTHIA INITIALIZATION
        // =====================================================================
        
        Pythia pythia;
        Settings& settings = pythia.settings;

        // Reduce output verbosity
        pythia.readString("Init:showMultipartonInteractions = off");
        pythia.readString("Init:showChangedSettings = on");  // Show what we changed
        pythia.readString("Init:showChangedParticleData = off");
        pythia.readString("Next:numberCount = 100000");
        pythia.readString("Next:numberShowInfo = 0");
        pythia.readString("Next:numberShowProcess = 0");
        pythia.readString("Next:numberShowEvent = 0");

        // ----- Beam Configuration -----
        pythia.readString("Beams:frameType = 2");
        pythia.readString("Beams:idA = 2212");              // Proton
        pythia.readString("Beams:idB = 11");                // Electron
        pythia.readString("Beams:eA = " + to_string(config.protonEnergy));
        pythia.readString("Beams:eB = " + to_string(config.electronEnergy));
        pythia.readString("PDF:beamB2gamma = on");          // Photon from electron

        // ----- Photoproduction Settings -----
        settings.mode("Photon:ProcessType", 0);             // Auto direct/resolved mix
        pythia.readString("Photon:Q2max = " + to_string(config.Q2max));
        pythia.readString("PhaseSpace:pTHatMin = " + to_string(config.pTHatMin));
        
        // ----- MPI Tuning (CRITICAL for γp!) -----
        pythia.readString("MultipartonInteractions:pT0Ref = " + to_string(config.pT0Ref));
        
        // ----- Fragmentation Parameters (optional tuning) -----
        pythia.readString("StringZ:aLund = " + to_string(config.StringZ_aLund));
        pythia.readString("StringZ:bLund = " + to_string(config.StringZ_bLund));
        pythia.readString("StringPT:sigma = " + to_string(config.StringPT_sigma));
        
        // ----- Enable Processes -----
        pythia.readString("HardQCD:all = on");              // Resolved processes
        pythia.readString("PhotonParton:all = on");         // Direct processes

        // Initialize generator
        if (!pythia.init()) {
            throw runtime_error("Failed to initialize PYTHIA");
        }
        
        // Print PYTHIA configuration check
        dout << "\n=== PYTHIA CONFIGURATION VERIFICATION ===\n";
        dout << "Tune:pp = " << pythia.settings.mode("Tune:pp") << "\n";
        dout << "Tune:ee = " << pythia.settings.mode("Tune:ee") << "\n";
        dout << "StringZ:aLund = " << pythia.settings.parm("StringZ:aLund") << "\n";
        dout << "StringZ:bLund = " << pythia.settings.parm("StringZ:bLund") << "\n";
        dout << "StringPT:sigma = " << pythia.settings.parm("StringPT:sigma") << "\n";
        dout << "MultipartonInteractions:pT0Ref = " << pythia.settings.parm("MultipartonInteractions:pT0Ref") << "\n";
        dout << "PhaseSpace:pTHatMin = " << pythia.settings.parm("PhaseSpace:pTHatMin") << "\n";
        dout << "Photon:Q2max = " << pythia.settings.parm("Photon:Q2max") << "\n";
        dout << "==========================================\n\n";

        // =====================================================================
        // CREATE OUTPUT FILE
        // =====================================================================
        
        TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
        if (!outputFile || outputFile->IsZombie()) {
            throw runtime_error("Failed to create output file: " + outputFileName);
        }
        
        // Create directory structure
        TDirectory* dirQQ = outputFile->mkdir("QQ_Events");
        TDirectory* dirGG = outputFile->mkdir("GG_Events");
        TDirectory* dirGQ = outputFile->mkdir("GQ_Events");
        TDirectory* dirCombined = outputFile->mkdir("Combined_Events");
        TDirectory* dirResolved = outputFile->mkdir("Resolved_Events");
        TDirectory* dirDirect = outputFile->mkdir("Direct_Events");
        TDirectory* dirConfig = outputFile->mkdir("Config");

        // =====================================================================
        // CREATE TREES FOR EVENT DATA
        // =====================================================================
        
        struct EventData {
            Int_t eventID, processCode, processType, photonMode;
            Bool_t isResolved, isDirect;
            Float_t crossSection, inelasticity;
            Float_t photonEnergy, scatteredElectronEnergy;
            
            vector<Float_t> px, py, pz, energy;
            vector<Float_t> pT, eta, phi, mass;
            vector<Int_t> pdgId, status;
            vector<Bool_t> isCharged, isHadron;

            // Outgoing hard-process partons (PYTHIA status == 23).
            // Used for jet-to-parton DeltaR matching checks downstream.
            vector<Int_t>   parton_pdgId;
            vector<Float_t> parton_pT, parton_eta, parton_phi;

            Float_t electronPx, electronPy, electronPz, electronE;
            Int_t nParticles;

            void clear() {
                px.clear(); py.clear(); pz.clear(); energy.clear();
                pT.clear(); eta.clear(); phi.clear(); mass.clear();
                pdgId.clear(); status.clear();
                isCharged.clear(); isHadron.clear();
                parton_pdgId.clear();
                parton_pT.clear(); parton_eta.clear(); parton_phi.clear();
            }
        };
        
        // Create data structures
        EventData* dataQQ = new EventData();
        EventData* dataGG = new EventData();
        EventData* dataGQ = new EventData();
        EventData* dataCombined = new EventData();
        EventData* dataResolved = new EventData();
        EventData* dataDirect = new EventData();
        
        // Helper to create tree with branches
        auto createTree = [](TDirectory* dir, const string& name, EventData* data) -> TTree* {
            dir->cd();
            TTree* tree = new TTree(name.c_str(), name.c_str());
            
            tree->Branch("eventID", &data->eventID);
            tree->Branch("processCode", &data->processCode);
            tree->Branch("processType", &data->processType);
            tree->Branch("photonMode", &data->photonMode);
            tree->Branch("isResolved", &data->isResolved);
            tree->Branch("isDirect", &data->isDirect);
            tree->Branch("crossSection", &data->crossSection);
            tree->Branch("inelasticity", &data->inelasticity);
            tree->Branch("nParticles", &data->nParticles);
            
            tree->Branch("px", &data->px);
            tree->Branch("py", &data->py);
            tree->Branch("pz", &data->pz);
            tree->Branch("energy", &data->energy);
            tree->Branch("pT", &data->pT);
            tree->Branch("eta", &data->eta);
            tree->Branch("phi", &data->phi);
            tree->Branch("mass", &data->mass);
            tree->Branch("pdgId", &data->pdgId);
            tree->Branch("status", &data->status);
            tree->Branch("isCharged", &data->isCharged);
            tree->Branch("isHadron", &data->isHadron);

            tree->Branch("parton_pdgId", &data->parton_pdgId);
            tree->Branch("parton_pT",    &data->parton_pT);
            tree->Branch("parton_eta",   &data->parton_eta);
            tree->Branch("parton_phi",   &data->parton_phi);

            return tree;
        };
        
        TTree* treeQQ = createTree(dirQQ, "QQ_Events", dataQQ);
        TTree* treeGG = createTree(dirGG, "GG_Events", dataGG);
        TTree* treeGQ = createTree(dirGQ, "GQ_Events", dataGQ);
        TTree* treeCombined = createTree(dirCombined, "Combined_Events", dataCombined);
        TTree* treeResolved = createTree(dirResolved, "Resolved_Events", dataResolved);
        TTree* treeDirect = createTree(dirDirect, "Direct_Events", dataDirect);

        // =====================================================================
        // EVENT GENERATION LOOP
        // =====================================================================
        
        dout << "Starting event generation...\n";
        dout << "Target: " << config.nEvents << " events\n\n";
        
        int totalEvents = 0, validEvents = 0;
        int qqEvents = 0, ggEvents = 0, gqEvents = 0, unclassifiedEvents = 0;
        int resolvedEvents = 0, directEvents = 0;
        double totalCrossSection = 0.0;
        double qqCrossSection = 0.0, ggCrossSection = 0.0, gqCrossSection = 0.0;
        
        // Process code counter for diagnostics
        map<int, int> processCodeCounts;
        
        int failedEvents = 0;
        
        for (int iEvent = 0; iEvent < config.nEvents; ++iEvent) {
            totalEvents++;
            
            // Progress - show every 10%
            if (iEvent % (config.nEvents/10) == 0) {
                double efficiency = (totalEvents > 0) ? 100.0 * validEvents / totalEvents : 0;
                dout << "  Progress: " << (100 * iEvent / config.nEvents) << "% | "
                     << "Valid: " << validEvents << " | "
                     << "Efficiency: " << fixed << setprecision(1) << efficiency << "%\n";
            }
            
            if (!pythia.next()) {
                failedEvents++;
                continue;
            }
            
            validEvents++;
            
            // Event classification
            int processCode = pythia.info.code();
            int processType = classifySubprocess(processCode);
            int photonMode = pythia.info.photonMode();
            
            // Classify direct vs resolved based on process code
            // Direct processes: photon interacts as point-like (codes 271-284)
            // Resolved processes: photon has partonic structure (codes 111-124)
            bool isDirect = (processCode >= 271 && processCode <= 284);
            bool isResolved = (processCode >= 111 && processCode <= 124);
            
            double weight = pythia.info.sigmaGen();
            
            // Count process codes
            processCodeCounts[processCode]++;
            
            // Update counters
            totalCrossSection += weight;
            if (processType == 1) { qqEvents++; qqCrossSection += weight; }
            else if (processType == 2) { ggEvents++; ggCrossSection += weight; }
            else if (processType == 3) { gqEvents++; gqCrossSection += weight; }
            else { unclassifiedEvents++; }
            
            if (isResolved) resolvedEvents++;
            if (isDirect) directEvents++;
            
            // Fill event data
            auto fillEventData = [&](EventData* data) {
                data->clear();
                data->eventID = validEvents;
                data->processCode = processCode;
                data->processType = processType;
                data->photonMode = photonMode;
                data->isResolved = isResolved;
                data->isDirect = isDirect;
                data->crossSection = weight;
                data->inelasticity = pythia.info.y();
                
                int nPart = 0;
                for (int i = 0; i < pythia.event.size(); ++i) {
                    const Particle& particle = pythia.event[i];

                    // Save the outgoing hard-process partons (status == 23).
                    // PYTHIA stores them with status code 23 BEFORE the
                    // shower turns them into descendants; that makes them
                    // the "LO outgoing parton" reference for jet matching.
                    if (particle.status() == 23) {
                        int absid = std::abs(particle.id());
                        if (absid == 21 || (absid >= 1 && absid <= 6)) {
                            data->parton_pdgId.push_back(particle.id());
                            data->parton_pT.push_back(particle.pT());
                            data->parton_eta.push_back(particle.eta());
                            data->parton_phi.push_back(particle.phi());
                        }
                    }

                    if (particle.isFinal() && particle.isVisible()) {
                        data->px.push_back(particle.px());
                        data->py.push_back(particle.py());
                        data->pz.push_back(particle.pz());
                        data->energy.push_back(particle.e());
                        data->pT.push_back(particle.pT());
                        data->eta.push_back(particle.eta());
                        data->phi.push_back(particle.phi());
                        data->mass.push_back(particle.m());
                        data->pdgId.push_back(particle.id());
                        data->status.push_back(particle.status());
                        data->isCharged.push_back(particle.isCharged());
                        data->isHadron.push_back(particle.isHadron());
                        nPart++;
                    }
                }
                data->nParticles = nPart;
            };
            
            // Fill combined tree
            fillEventData(dataCombined);
            treeCombined->Fill();
            
            // Fill subprocess trees
            if (processType == 1) { fillEventData(dataQQ); treeQQ->Fill(); }
            else if (processType == 2) { fillEventData(dataGG); treeGG->Fill(); }
            else if (processType == 3) { fillEventData(dataGQ); treeGQ->Fill(); }
            
            // Fill photoproduction type trees
            if (isResolved) { fillEventData(dataResolved); treeResolved->Fill(); }
            if (isDirect) { fillEventData(dataDirect); treeDirect->Fill(); }
        }
        
        dout << "\nEvent generation complete!\n\n";

        // =====================================================================
        // SAVE CONFIGURATION TO ROOT FILE
        // =====================================================================
        
        dirConfig->cd();
        
        // Save as TNamed objects (readable with TFile::ls())
        TNamed("ColliderPreset", config.getColliderName().c_str()).Write();
        TNamed("sqrtS", to_string(config.sqrtS()).c_str()).Write();
        TNamed("protonEnergy", to_string(config.protonEnergy).c_str()).Write();
        TNamed("electronEnergy", to_string(config.electronEnergy).c_str()).Write();
        TNamed("pTHatMin", to_string(config.pTHatMin).c_str()).Write();
        TNamed("Q2max", to_string(config.Q2max).c_str()).Write();
        TNamed("pT0Ref", to_string(config.pT0Ref).c_str()).Write();
        TNamed("StringZ_aLund", to_string(config.StringZ_aLund).c_str()).Write();
        TNamed("StringZ_bLund", to_string(config.StringZ_bLund).c_str()).Write();
        TNamed("StringPT_sigma", to_string(config.StringPT_sigma).c_str()).Write();
        TNamed("nEventsGenerated", to_string(validEvents).c_str()).Write();
        TNamed("TunePP", to_string(config.tunePP).c_str()).Write();
        TNamed("TuneEE", to_string(config.tuneEE).c_str()).Write();
        
        // Save as TTree for easy reading
        TTree* configTree = new TTree("GenerationConfig", "PYTHIA Configuration");
        Double_t cfg_sqrtS = config.sqrtS();
        Double_t cfg_protonE = config.protonEnergy;
        Double_t cfg_electronE = config.electronEnergy;
        Double_t cfg_pTHatMin = config.pTHatMin;
        Double_t cfg_Q2max = config.Q2max;
        Double_t cfg_pT0Ref = config.pT0Ref;
        Double_t cfg_aLund = config.StringZ_aLund;
        Double_t cfg_bLund = config.StringZ_bLund;
        Double_t cfg_sigma = config.StringPT_sigma;
        Int_t cfg_nEvents = validEvents;
        Int_t cfg_tunePP = config.tunePP;
        Int_t cfg_tuneEE = config.tuneEE;
        
        configTree->Branch("sqrtS", &cfg_sqrtS);
        configTree->Branch("protonEnergy", &cfg_protonE);
        configTree->Branch("electronEnergy", &cfg_electronE);
        configTree->Branch("pTHatMin", &cfg_pTHatMin);
        configTree->Branch("Q2max", &cfg_Q2max);
        configTree->Branch("pT0Ref", &cfg_pT0Ref);
        configTree->Branch("StringZ_aLund", &cfg_aLund);
        configTree->Branch("StringZ_bLund", &cfg_bLund);
        configTree->Branch("StringPT_sigma", &cfg_sigma);
        configTree->Branch("nEventsGenerated", &cfg_nEvents);
        configTree->Branch("TunePP", &cfg_tunePP);
        configTree->Branch("TuneEE", &cfg_tuneEE);
        configTree->Fill();

        // =====================================================================
        // PRINT STATISTICS
        // =====================================================================
        
        pythia.stat();
        
        dout << "\n=============================================================================\n";
        dout << "EVENT GENERATION SUMMARY\n";
        dout << "=============================================================================\n";
        dout << "Output File: " << outputFileName << "\n";
        dout << "Log File:    " << logFileName << "\n";
        dout << "\n--- Event Counts ---\n";
        dout << "Total Attempts:     " << totalEvents << "\n";
        dout << "Valid Events:       " << validEvents << "\n";
        dout << "Failed Events:      " << failedEvents << "\n";
        dout << "Efficiency:         " << fixed << setprecision(1) << (100.0 * validEvents / totalEvents) << "%\n";
        dout << "\n--- Subprocess Classification ---\n";
        dout << "  QQ (quark-init):  " << setw(8) << qqEvents 
             << " (" << fixed << setprecision(1) << (100.0 * qqEvents / validEvents) << "%)\n";
        dout << "  GG (gluon-init):  " << setw(8) << ggEvents 
             << " (" << (100.0 * ggEvents / validEvents) << "%)\n";
        dout << "  GQ (mixed):       " << setw(8) << gqEvents 
             << " (" << (100.0 * gqEvents / validEvents) << "%)\n";
        dout << "  Unclassified:     " << setw(8) << unclassifiedEvents 
             << " (" << (100.0 * unclassifiedEvents / validEvents) << "%)\n";
        dout << "\n--- Photoproduction Type ---\n";
        dout << "  Resolved:         " << setw(8) << resolvedEvents 
             << " (" << (100.0 * resolvedEvents / validEvents) << "%)\n";
        dout << "  Direct:           " << setw(8) << directEvents 
             << " (" << (100.0 * directEvents / validEvents) << "%)\n";
        
        dout << "\n--- Process Code Breakdown ---\n";
        dout << setw(8) << "Code" << setw(25) << "Process" << setw(10) << "Count" 
             << setw(10) << "%" << setw(10) << "Class\n";
        dout << string(63, '-') << "\n";
        
        for (const auto& pc : processCodeCounts) {
            int code = pc.first;
            int count = pc.second;
            string procName = getProcessName(code);
            int procType = classifySubprocess(code);
            string className = (procType == 1) ? "QQ" : (procType == 2) ? "GG" : (procType == 3) ? "GQ" : "??";
            
            dout << setw(8) << code 
                 << setw(25) << procName 
                 << setw(10) << count 
                 << setw(9) << fixed << setprecision(2) << (100.0 * count / validEvents) << "%"
                 << setw(10) << className << "\n";
        }
        
        dout << "\n--- Trees Saved ---\n";
        dout << "  QQ_Events:       " << setw(8) << treeQQ->GetEntries() << " events\n";
        dout << "  GG_Events:       " << setw(8) << treeGG->GetEntries() << " events\n";
        dout << "  GQ_Events:       " << setw(8) << treeGQ->GetEntries() << " events\n";
        dout << "  Combined_Events: " << setw(8) << treeCombined->GetEntries() << " events\n";
        dout << "  Resolved_Events: " << setw(8) << treeResolved->GetEntries() << " events\n";
        dout << "  Direct_Events:   " << setw(8) << treeDirect->GetEntries() << " events\n";
        dout << "=============================================================================\n";

        // =====================================================================
        // WRITE AND CLOSE
        // =====================================================================
        
        outputFile->cd();
        outputFile->Write();
        outputFile->Close();  // Properly close the file (ROOT handles tree cleanup)
        
        dout << "\nEvent generation completed successfully!\n";
        dout << "Output: " << outputFileName << "\n";
        
        // Cleanup - only delete our data structures, not ROOT objects
        delete dataQQ;
        delete dataGG;
        delete dataGQ;
        delete dataCombined;
        delete dataResolved;
        delete dataDirect;
        
        // Note: Don't delete outputFile - Close() handles cleanup
        // ROOT manages TTree and TDirectory objects internally
        
    } catch (const exception& e) {
        cerr << "\nERROR: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}