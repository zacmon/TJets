//  Telescoping Jets Package
//
//  Alexander Emerman, Yang-Ting Chien, Shih-Chieh Hsu, Zachary Montague
//

#ifndef __TELESCOPINGJETS_HH__
#define __TELESCOPINGJETS_HH__

#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <memory> 

#include <TROOT.h>
#include "TLorentzVector.h"

//#include "fastjet/contrib/AxesFinder.hh"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>

#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "fastjet/contrib/AxesDefinition.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

struct tSub {
    double minAngle = -1;
    double massVariability = -1;
    double pTVariability = -1;
    
    std::vector<std::vector<double>> masses;
    std::vector<std::vector<double>> pTs;
    std::vector<TLorentzVector> tauAxes;

    double targetMass = -1;
    double targetMassVariability = -1;
    double targetPTVariability = -1;
};

double getVariability(const std::vector<double>& values);
double getAverage(const std::vector<double>& values);
double getRMS(const std::vector<double>& values);
unsigned int choose(unsigned int n, unsigned int k);
TLorentzVector convertPseudoJet2TLV(fastjet::PseudoJet pseudoJet);

class TelescopingJets {

public:
    TelescopingJets(const fastjet::PseudoJet& pseudoJet);
    TelescopingJets(const fastjet::PseudoJet& pseudoJet, int axesType);
    TelescopingJets(const fastjet::PseudoJet& pseudoJet, int axesType, int scale);
    ~TelescopingJets();

    std::vector<double> getTelescopingParameterSet(double minParameter, double maxParameter, unsigned int stepSize);  

    double tPruning(double minDCut, double maxDCut, unsigned int numDCuts);
    double tTrimming(double minFCut, double maxFCut, unsigned int numFCuts);
        
    double tReclustering(int algorithm, double minRadius, double maxRadius, unsigned int numRadii);

    tSub tNSubjet(unsigned int numSubjets, double minRadius, double maxRadius, unsigned int numRadii, double targetMass);

    double tNsubjettiness(int numSubjets, double minBeta, double maxBeta, unsigned int numBetas);
    double tNsubjettinessRatio(int nNumerator, int nDemoninator, double minBeta, double maxBeta, unsigned int numBetas);

    double tEnergyCorrelator_C2(double minBeta, double maxBeta, unsigned int numBetas);
    double tEnergyCorrelator_D2(double minBeta, double maxBeta, unsigned int numBetas);
    double tEnergyCorrelator_C3(double minBeta, double maxBeta, unsigned int numBetas);

private:
    const fastjet::PseudoJet& input;
    fastjet::contrib::AxesDefinition* axesDefinition;
    const int stepScale;

    fastjet::contrib::AxesDefinition* getAxesDefinition(int axesType);
    std::vector<TLorentzVector> getTauAxes(unsigned int numSubjets, double beta);
    std::vector<double> getAnglesBetweenTauAxes(unsigned int numSubjets, std::vector<TLorentzVector> pTauAxes);    
    std::vector<std::vector<std::pair<TLorentzVector, double>>> sortConstituents(unsigned int numSubjets, std::vector<TLorentzVector> pTauAxes);

    bool emptyTelescopingMasses(std::vector<double> telescopingMasses, std::vector<std::vector<double>> combinedTelescopingMasses);

    tSub telescopeSubjets(unsigned int numSubjets, tSub result, std::vector<double> subjetRadii, std::vector<std::vector<std::pair<TLorentzVector, double>>> constituents);
    tSub telescopeSubjets(unsigned int numSubjets, tSub result, std::vector<double> subjetRadii, std::vector<std::vector<std::pair<TLorentzVector, double>>> constituents, double targetMass); 

    unsigned int getCandidateIndex(int targetMass, std::vector<std::vector<double>> telescopingMasses);
};

#endif  // __TELESCOPINGJETS_HH__
