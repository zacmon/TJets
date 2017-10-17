#include <iostream>
#include <vector>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>

double getAverage(const std::vector<double>& values) {
    if (values.empty()) throw std::length_error("Asking for average of empty vector.\n");
    double v = 0;
    for (unsigned int i = 0; i < values.size(); i++) {
	v += values.at(i);
    }
    return v/values.size();
}

double getRMS(const std::vector<double>& values) {
    if (values.empty()) throw std::length_error("Asking for rms of empty vector.\n");
    double v = getAverage(values);
    double v2 = 0;
    for(unsigned int i=0; i < values.size(); i++) {
	v2 += (values.at(i) - v) * (values.at(i) - v);
    }
    double sqRMS = v2/values.size();
    return std::sqrt(sqRMS);
}

double getVolatility(const std::vector<double>& values) {
    if (values.empty()) throw std::length_error("Asking for volatility of empty vector.\n");
    return getRMS(values) / getAverage(values);
}

std::vector<double> setVectorWithNElements(const std::vector<double>& v, double numElements) {
    unsigned int currentNumElements = v.size();
    double stepSize = currentNumElements / numElements;
    if (stepSize < 1) return {};

    std::vector<double> newVec;
    double k = 0;
    while (k < currentNumElements - 1) {
	newVec.push_back(v[std::round(k)]);
	k += stepSize;
    }

    return newVec;
}

std::vector<double> getVolForEvent(const std::vector<double>& masses, unsigned int maxMinRIndex ) {
    double numElements = 33;
    
    std::vector<double> volatilities;
    for (unsigned int j = 0; j < maxMinRIndex; ++j) {
	//  Start vector at correct radius size. This insures we don't go beyond maxMinR.
	std::vector<double> trimmedVector(masses.begin() + j, masses.end());
	//  Makes sure each vector has the same amount of elements so volatility comparisons are legitimate.
	std::vector<double> v = setVectorWithNElements(trimmedVector, numElements);
	volatilities.push_back(getVolatility(trimmedVector));
    }
    return volatilities;
}

void OptimizeTParam() {
    std::vector<double> massBounds = {160, 190};
    std::vector<double> pTBounds = {800, 1000};
    
    int numRadii = 100;
    std::vector<double> radii;
    double minRadius = 0.01;
    double maxRadius = 1;
    double deltaR = (maxRadius - minRadius) / (numRadii);

    unsigned int maxMinRIndex = 9999;
    double maxMinR = 0.5;
    
    for (double r = minRadius; r < maxRadius + deltaR; r += deltaR) {
	radii.push_back(r);
	if ((r - maxMinR) < deltaR) maxMinRIndex = radii.size();
    }

    std::vector<TString> fileNames = {"Dijet"};//, "Top", "W"};
    std::vector<TFile*> files;
    std::vector<TTree*> trees;
    for (auto const &name : fileNames) {
    	TFile* file = new TFile("~/optimizeHigh" + name + ".root", "READ");
    	TTree* tree = (TTree*)file->Get("JetTree");
    	files.push_back(file);
    	trees.push_back(tree);
    }

    std::vector<std::vector<double>>* t1Masses = nullptr;
    std::vector<std::vector<double>>* t2Masses = nullptr;
    std::vector<std::vector<double>>* t3Masses = nullptr;

    std::vector<double>* mass = nullptr;
    std::vector<double>* pT = nullptr;
    
    for (auto const &tree : trees){
	Long64_t numEntries = tree->GetEntries();
	
	std::vector<std::vector<double>> t1Volatilities(numEntries);
	std::vector<std::vector<double>> t2Volatilities(numEntries);
	std::vector<std::vector<double>> t3Volatilities(numEntries);

	tree->SetBranchAddress("TruthRawTrim_m", &mass);
	tree->SetBranchAddress("TruthRawTrim_pt", &pT);
	
    	tree->SetBranchAddress("TruthRawTrim_T1masses", &t1Masses);
    	tree->SetBranchAddress("TruthRawTrim_T2masses", &t2Masses);
    	tree->SetBranchAddress("TruthRawTrim_T3masses", &t3Masses);
	
	
	for (Int_t entry = 0; entry < numEntries; ++entry) {
	    tree->GetEntry(entry);

	    for (unsigned int i = 0; i < mass->size(); ++i) {
		if ((mass->at(i) > massBounds[0] && mass->at(i) < massBounds[1]) && (pT->at(i) >  pTBounds[0] && pT->at(i) < pTBounds[1])) {
		    std::cout << "mass: " << mass->at(i) << " pt: " << pT->at(i) << std::endl;
       	    	    t1Volatilities[entry] = getVolForEvent(t1Masses->at(i), maxMinRIndex);
		    t2Volatilities[entry] = getVolForEvent(t2Masses->at(i), maxMinRIndex);
		    t3Volatilities[entry] = getVolForEvent(t3Masses->at(i), maxMinRIndex);
		}
	    }
	}
    }
    //  TODO
    //  Create TTree and TFile to save results
    //  Save mass and pt for mass and pt cuts in OverlayROCS
    //  Save vector of vectors to TTree -> three vector of vectors
    //  Look at interesting areas, say minR=0.01 to minR=0.2, spaced out evenly and see ROC curves by Overlay ROC curve stuff
    //  ID best minR
    //  With new best minR, find best number of radii
}
