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
    double massVolatility = -1;
    double pTVolatility = -1;
    std::vector<double> masses;
    std::vector<double> pTs;

    double targetMass = -1;
    double targetMassVolatility = -1;
    double targetPTVolatility = -1;
};

//  Helper functions for calculation volatilities.
double getVolatility(const std::vector<double>& values);
double getAverage(const std::vector<double>& values);
double getRMS(const std::vector<double>& values);
unsigned long long getChoose(unsigned long long n, unsigned long long k);

class TelescopingJets {

public:
    TelescopingJets(const fastjet::PseudoJet& pseudoJet);
    ~TelescopingJets();


    std::vector<double> getTelescopingParameterSet(double minParameter, double maxParameter, int numParameter, int stepScale);  

    double tPruning(double minDCut, double maxDCut, int numDCuts);
    double tTrimming(double minFCut, double maxFCut, int numFCuts);
        
    double tReclustering(int algorithm, double minRadius, double maxRadius, int numRadii, int stepScale);

    tSub tNSubjet(unsigned int numSubjets, double minRadius, double maxRadius, int numRadii, int stepScale, int axesType, double targetMass);

    double tNsubjettiness(int numSubjets, double minBeta, double maxBeta, int numBetas, int axesType);
    double tNsubjettinessRatio(int nNumerator, int nDemoninator, double minBeta, double maxBeta, int numBetasi, int axesType);

    double tEnergyCorrelator_C2(double minBeta, double maxBeta, int numBetas);
    double tEnergyCorrelator_D2(double minBeta, double maxBeta, int numBetas);
    double tEnergyCorrelator_C3(double minBeta, double maxBeta, int numBetas);

private:
    const fastjet::PseudoJet& input;

    std::vector<TLorentzVector> convertPseudoJet2TLV(std::vector<fastjet::PseudoJet> pseudoJet);

    const fastjet::contrib::AxesDefinition* getAxesDefinition(int axesType);
    std::vector<TLorentzVector> getTauAxes(unsigned int numSubjets, double beta, int axesType);
    std::vector<double> getAnglesBetweenTauAxes(unsigned int numSubjets, std::vector<TLorentzVector> pTauAxes);    
    std::vector<std::vector<std::pair<TLorentzVector, double>>> sortConstituents(unsigned int numSubjets, std::vector<TLorentzVector> pTauAxes);

    bool emptyTelescopingMasses(std::vector<std::vector<double>> telescopingMasses);

    tSub telescopeSubjets(unsigned int numSubjets, tSub result, std::vector<double> subjetRadii, std::vector<std::vector<std::pair<TLorentzVector, double>>> constituents);
    tSub telescopeSubjets(unsigned int numSubjets, tSub result, std::vector<double> subjetRadii, std::vector<std::vector<std::pair<TLorentzVector, double>>> constituents, double targetMass); 

};

#endif  // __TELESCOPINGJETS_HH__
