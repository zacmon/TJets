#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>

void OptimizeTelescope() {
    int numRadii = 36;
    std::vector<double> radii;
    double minRadius = 0.1;
    double maxRadius = 1;
    double deltaR = (maxRadius - minRadius) / (numRadii);

    for (double r = minRadius; r <= maxRadius; r += deltaR) {
	radii.push_back(r);
    }
    
    TString fileName = "../Ana_EventGeneration/";
    TFile* newFile = new TFile(fileName+"blah.root", "READ");
    TTree* newTree = (TTree*)newFile->Get("JetTree");

    TFile* oldFile = new TFile(fileName+"ntuple_tt_test50000.root", "READ");
    TTree* oldTree = (TTree*)oldFile->Get("JetTree");

    TString var = "TruthRawTrim_T3jet";
    
    newTree->Draw(var+">>h1(100,0,1.0)");
    oldTree->Draw(var+">>h2(100,0,1.0)");
    TH1D *h1 = (TH1D*)gDirectory->Get("h1");
    TH1D *h2 = (TH1D*)gDirectory->Get("h2");
    h1->Divide(h2);
    TCanvas* c = new TCanvas("c", "c");
    h1->Draw();   
}
