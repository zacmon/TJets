//  Telescoping Jets Package
//
//  Alexander Emerman, Yang-Ting Chien, Shih-Chieh Hsu, Zachary Montague
//

#ifndef __TELESCOPINGJETS_HH__
#define __TELESCOPINGJETS_HH__

//#include "fastjet/contrib/AxesFinder.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/JetDefintion.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "fastjet/contrib/AxesDefinition.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

#include <cmath>
#include <vector>
#include <stdexcept>

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

struct t2Sub {
    double minAngle = -1;
    double massVolatility = -1;
    double pTVolatility = -1;
    std::vector<double> masses;

    double wMass = -1;
    double wMassVolatility = -1;
};
    
struct t3Sub {
    double minAngle = -1;
    double massVolatility = -1;
    double pTVolatility = -1;
    std::vector<double> masses;

    double midAngle = -1;
    double maxAngle = -1;

    double wMass = -1;
    double wMassVolatility = -1;
    double wPTVolatility = -1;
};

//  Functions to calculate volatility.
double getVolatility(const std::vector<double>& values);
double getAverage(const std::vector<double>& values);
double getRMS(const std::vector<double>& values);

class TelescopingJets {

public:
    TelescopingJets(const fastjet::PseudoJet& input);
    ~TelescopingJets();
    
    double tPruning(double minDCut, double maxDCut, int numDCuts);
    double tTrimming(double minFCut, double maxFCut, int numFCuts);

    double tReclustering(int algorithm, double minRadius, double maxRadius, int numRadii);

    tSub tNSubjet(unsigned int numSubjets, double minRadius, double maxRadius, int numRadii, int stepScale);
    t2Sub t2Subjet(double minRadius, double maxRadius, int numRadii, int stepScale);
    t3Sub t3Subjet(double minRadius, double maxRadius, int numRadii, int stepScale);

    double tNsubjettiness(int numSubjets, double minBeta, double maxBeta, int numBetas);
    double tNsubjettinessRatio(int nNumerator, int nDemoninator, double minBeta, double maxBeta, int numBetas);

    double T_EnergyCorrelator_C2(double minBeta, double maxBeta, int numBetas);
    double T_EnergyCorrelator_D2(double minBeta, double maxBeta, int numBetas);
    double T_EnergyCorrelator_C3(double minBeta, double maxBeta, int numBetas);

private:
    fastjet::PseduoJet& input;

    std::vector<TLorentzVector> convertPsuedoJet2TLV(std::vector<PseudoJet> pseudoJet);
    std::vector<double> getAnglesBetweenTauAxes(std::vector<TLorentzVector> pTauAxes)
    std::vector<std::vector<TLorentzVector>> sortConstituents(std::vector<TLorentzVector> pTauAxes);
    

};

#endif  // __TELESCOPINGJETS_HH__
