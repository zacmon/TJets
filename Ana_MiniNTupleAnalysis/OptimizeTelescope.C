#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>

void OptimizeTelescope() {
    int numRadii = 40;
    std::vector<double> radii;
    double minRadius = 0.1;
    double maxRadius = 1.5;
    double deltaR = (maxRadius - minRadius) / (numRadii - 1);

    for (double r = minRadius; r <= maxRadius; r += deltaR) {
	radii.push_back(r);
    }
    
    TString fileName = "../Ana_EventGeneration/ntuple_tt_50000v32.root";
    TFile* file = new TFile(fileName, "READ");
    TTree* tree = (TTree*)file->Get("JetTree");

    std::vector<std::vector<double>> *vol3 = nullptr;
    std::vector<double> *v32 = nullptr;
    
    tree->SetBranchAddress("TruthRawTrim_T3Volatility", &vol3);
    tree->SetBranchAddress("v32", &v32);
    
    Long64_t numEntries = tree->GetEntries();

    TH1D* h = new TH1D("h", ";v32", 100, 0, 3);
    TH2D* hist = new TH2D("hist", ";r;volatility", 40, minRadius-deltaR, maxRadius+deltaR, 50, 0, 0.3);
    
    for (Long64_t i = 0; i < numEntries; ++i) {
	tree->GetEntry(i);
	for (unsigned int j = 0; j < vol3->size(); ++j) {
	    for (unsigned int k = 0; k < vol3->at(j).size(); ++k) {
		hist->Fill(radii[k], vol3->at(j).at(k));
	    }
	}
	for (unsigned int j = 0; j < v32->size(); ++j) {
	    h->Fill(v32->at(j));
	}
    }

    TCanvas* c = new TCanvas("c", "c");
    h->Draw();
    
}
