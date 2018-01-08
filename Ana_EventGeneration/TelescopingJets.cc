//  Telescoping Jets Package
//
//  Alexander Emerman, Yang-Ting Chien, Shih-Chieh Hsu, Zachary Montague
//

#include "TelescopingJets.hh"

// Useful mathematical functions to find the average, rms and volatility
// of a set of values
double getAverage(const std::vector<double>& values) {
    if (values.empty()) throw std::length_error("Asking for average of empty vector.\n");
    double v = 0;
    for (unsigned int i = 0; i < values.size(); ++i) {
	v += values.at(i);
    }
    return v/values.size();
}
double getRMS(const std::vector<double>& values) {
    if (values.empty()) throw std::length_error("Asking for rms of empty vector.\n");
    double v = getAverage(values);
    double v2 = 0;
    for (unsigned int i = 0; i < values.size(); ++i) {
	v2 += (values.at(i) - v) * (values.at(i) - v);
    }
    double sqRMS = v2 / values.size();
    return std::sqrt(sqRMS);
}
double getVolatility(const std::vector<double>& values) {
    if (values.empty()) throw std::length_error("Asking for volatility of empty vector.\n");
    return getRMS(values) / getAverage(values);
}
   
unsigned int getChoose(unsigned int n, unsigned int k) {
   if (k > n) return 0;

   unsigned int r = 1;
   for (unsigned int d = 1; d <= k; ++d) {
       r *= n--;
       r /= d;
   }
   return r;
}

TelescopingJets::TelescopingJets(const fastjet::PseudoJet& pseudoJet) : input(pseudoJet), axesDefinition(getAxesDefinition(0)), stepScale(0) {
}

TelescopingJets::TelescopingJets(const fastjet::PseudoJet& pseudoJet, int axesType) : input(pseudoJet), axesDefinition(getAxesDefinition(axesType)), stepScale(0) {
}

TelescopingJets::TelescopingJets(const fastjet::PseudoJet& pseudoJet, int axesType, int scale) : input(pseudoJet), axesDefinition(getAxesDefinition(axesType)), stepScale(scale) {
}

TelescopingJets::~TelescopingJets() {
    delete axesDefinition;
}

TLorentzVector TelescopingJets::convertPseudoJet2TLV(fastjet::PseudoJet pseudoJet) {
    TLorentzVector tLorentzVector;
    tLorentzVector.SetPxPyPzE(pseudoJet.px(), pseudoJet.py(), pseudoJet.pz(), pseudoJet.e());
    return tLorentzVector;
}

std::vector<double> TelescopingJets::getTelescopingParameterSet(double minParameter, double maxParameter, unsigned int numSteps) {
    double stepSize = -999;
    //  Set the step for linear (0) or log (1).
    //  TODO 
    //  Add other options?
    if (stepScale == 0) stepSize = (maxParameter - minParameter) / (numSteps);
    else if (stepScale == 1) stepSize = log10(maxParameter / minParameter) / (numSteps);
    
    std::vector<double> parameterSet(numSteps + 1);
   
    //  Add values to the set.
    for (unsigned int i = 0; i <= numSteps; ++i) {
	if (stepScale == 0) parameterSet[i] = minParameter + i * stepSize;
	else if (stepScale == 1) parameterSet[i] = minParameter * pow(10, i * stepSize);
    }

    return parameterSet;
}

///=========================================
/// Telescoping Pruning
///=========================================
double TelescopingJets::tPruning(double minDCut, double maxDCut, unsigned int numDCuts) {
    //  Single choice of zcut but can be telescoped further.
    double zCut = 0.1;
    
    std::vector<double> dCuts = getTelescopingParameterSet(minDCut, maxDCut, numDCuts);
    std::vector<double> telescopingMasses(dCuts.size());
    
    //  Obtain a set of jet masses by stepping through the DCut
    //  pruning parameter.
    for (unsigned int i = 0; i < dCuts.size(); ++i) {
        fastjet::Pruner pruner(fastjet::cambridge_algorithm, zCut, dCuts[i]);
        fastjet::PseudoJet prunedJet = pruner(input);
        double mass = prunedJet.m();
        if (mass > 0.01) telescopingMasses[i] = mass;
    }
    
    if (!telescopingMasses.empty()) return getVolatility(telescopingMasses);
    
    std::cout << "WARNING zero entries for T_Pruning!   minDCut: " << minDCut <<
	"   maxDCut: " << maxDCut <<
	"   numDCuts: " << numDCuts <<
	"   input mass: " << input.m() << std::endl;
    return -1;
}


///=========================================
/// Telescoping Trimming
///=========================================
double TelescopingJets::tTrimming(double minFCut, double maxFCut, unsigned int numFCuts) {
    //  The subjet radius used should depend on the input jet's
    //  pT--higher pT should have a smaller subjet radius.
    double Rfilt = 0.1; 
    if (input.pt() > 500) Rfilt = 0.05;
    
    std::vector<double> fCuts = getTelescopingParameterSet(minFCut, maxFCut, numFCuts);
    std::vector<double> telescopingMasses(fCuts.size());

    //  Obtain a set of jet masses by stepping through the f trimming parameter.
    for (unsigned int i = 0; i < fCuts.size(); ++i) {
        fastjet::Filter trimmer(Rfilt, fastjet::SelectorPtFractionMin(fCuts[i]));
        fastjet::PseudoJet trimmedJet = trimmer(input);
        double mass = trimmedJet.m();
        if (mass > 0.01) telescopingMasses[i] = mass;
    }

    if (!telescopingMasses.empty()) return getVolatility(telescopingMasses);
    
    std::cout << "WARNING zero entries for T_Trimming!   minFCut: "<< minFCut <<
	"   maxFCut: "<< maxFCut <<
	"   numFCuts: "<< numFCuts <<
	"   input mass: " << input.m() << std::endl;
    return -1;
}

///=========================================
/// Telescoping reclustering
/// for kt, set algorithm = 0
/// for antikt, set algorithm = 2
///=========================================
double TelescopingJets::tReclustering(int algorithm, double minRadius, double maxRadius, unsigned numRadii) {
    std::vector<double> subjetRadii = getTelescopingParameterSet(minRadius, maxRadius, numRadii);
    std::vector<double> telescopingMass(subjetRadii.size());

    for (unsigned int i = 0; i < subjetRadii.size(); ++i) {
        fastjet::JetDefinition TjetDefinition(fastjet::JetAlgorithm(algorithm), subjetRadii[i]);
        fastjet::ClusterSequence TClusterSequence(input.constituents(), TjetDefinition);
        std::vector<fastjet::PseudoJet> recoTJets = sorted_by_pt(TClusterSequence.inclusive_jets());
	
        if (recoTJets.empty()) {
            std::cout << "WARNING! Recluster number of subjets is " << recoTJets.size() << std::endl;
            continue;
        }
        if (recoTJets.size() == 1) {
            TLorentzVector Tsubjet1 = convertPseudoJet2TLV(recoTJets[0]); 
            double mass = Tsubjet1.M();
            if (mass > 0.01) telescopingMass[i] = mass;
        }
        else if (recoTJets.size() >= 2) {
            TLorentzVector Tsubjet1 = convertPseudoJet2TLV(recoTJets[0]);
            TLorentzVector Tsubjet2 = convertPseudoJet2TLV(recoTJets[1]);
            double mass = (Tsubjet1 + Tsubjet2).M();
            if (mass > 0.01) telescopingMass[i] = mass;
        }
    }
    
    if (!telescopingMass.empty()) return getVolatility(telescopingMass);
    
    std::cout << "WARNING: zero entries for T_reclustering!   Algorithm" << algorithm <<
	"   minRadius: " << minRadius <<
	"   maxRadius: " << maxRadius <<
	"   numRadii: " << numRadii <<  std::endl;
    return -1;
}

///=========================================
/// Telescoping Subjet
///=========================================
fastjet::contrib::AxesDefinition* TelescopingJets::getAxesDefinition(int axesType) { 
    fastjet::contrib::AxesDefinition* axesDefinition = nullptr;
    if (axesType == 0) axesDefinition = new fastjet::contrib::OnePass_KT_Axes();
    else if (axesType == 1) axesDefinition = new fastjet::contrib::WTA_KT_Axes();
    else std::cout << "Error! Invalid axesType entered. Please enter \"0\" for OnePass_KT_Axes or \"1\" for WTA_KT_Axes." << std::endl;
    return axesDefinition;
}

std::vector<TLorentzVector> TelescopingJets::getTauAxes(unsigned int numSubjets, double beta) {   
    //  Create nSubjettiness object.
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSubjettiness(numSubjets, *axesDefinition, nSubMeasure);

    //  Get tau axes.
    double tauN = nSubjettiness.result(input);
    std::vector<fastjet::PseudoJet> tauAxes = nSubjettiness.currentAxes();

    std::vector<TLorentzVector> pTauAxes;
    std::transform(tauAxes.begin(), tauAxes.end(), std::back_inserter(pTauAxes), 
            std::bind(&TelescopingJets::convertPseudoJet2TLV));
    
    return pTauAxes;
}

std::vector<double> TelescopingJets::getAnglesBetweenTauAxes(unsigned int numSubjets, std::vector<TLorentzVector> pTauAxes) {
    std::vector<double> angles;

    //  Calculate deltaR between each pair of tau axes.
    for (unsigned int i = 0; i < numSubjets; ++i) {
	for (unsigned int j = i + 1; j < numSubjets; ++j) {
	    angles.push_back(pTauAxes[i].DeltaR(pTauAxes[j]));
	}
    }

    return angles;
}

std::vector<std::vector<std::pair<TLorentzVector, double>>> TelescopingJets::sortConstituents(unsigned int numSubjets, std::vector<TLorentzVector> pTauAxes) {
    //  Get constituents of the input jet.
    std::vector<fastjet::PseudoJet> constituents = input.constituents();

    //  Initialize container that will hold constiuents sorted by
    //  their distance to some tau axis and then sorted by increasing distance.
    std::vector<std::vector<std::pair<TLorentzVector, double>>> sortedConstituents(numSubjets);
    
    for (auto const &constituent : constituents) {
	//  Convert PseudoJet four-vector to TLorentz four-vector.
	TLorentzVector pConstituent = convertPseudoJet2TLV(constituent);

	//  Get distances to different tau axes.
	std::vector<double> distanceToTauAxes(numSubjets);
	for (unsigned int i = 0; i < numSubjets; ++i) {
	    distanceToTauAxes[i] = pConstituent.DeltaR(pTauAxes[i]);
	}

	//  Identify tau axes associated with minimum distance.
	auto minDistanceToTauAxesIt = std::min_element(distanceToTauAxes.cbegin(), distanceToTauAxes.cend());
	double minDistanceToTauAxes = *minDistanceToTauAxesIt;
	unsigned int minDistanceToTauAxesIndex = std::distance(distanceToTauAxes.cbegin(), minDistanceToTauAxesIt);
	
	sortedConstituents[minDistanceToTauAxesIndex].push_back({pConstituent, minDistanceToTauAxes});
    }

    //  Sort constituents by increasing distance for each subjet.
    for (unsigned int i = 0; i < numSubjets; ++i) {
	std::sort(sortedConstituents[i].begin(),
                sortedConstituents[i].end(),
                [](const std::pair<TLorentzVector, double> &pair1,
                   const std::pair<TLorentzVector, double> &pair2) {
		      return pair1.second < pair2.second;
		  });
    }
    
    return sortedConstituents;
}

tSub TelescopingJets::telescopeSubjets(unsigned int numSubjets, tSub result, std::vector<double> subjetRadii, std::vector<std::vector<std::pair<TLorentzVector, double>>> constituents) {
    //  Vector of subjet masses.
    std::vector<TLorentzVector> TSubjets(numSubjets);

    //  Vectors for mass and pT raw data for jet and subjets.
    std::vector<std::vector<double>> telescopingMasses(numSubjets + 1, std::vector<double>(subjetRadii.size()));
    std::vector<std::vector<double>> telescopingPTs(numSubjets + 1, std::vector<double>(subjetRadii.size()));

    for (unsigned int i = 0; i < subjetRadii.size(); ++i) {
	for (unsigned int j = 0; j < constituents.size(); ++j) {
	    for (auto it = constituents[j].begin(); it != constituents[j].end(); ++it) {
		//  Add constituent mass to subjet mass if within range.
		if (it->second <= subjetRadii[i]) {
		    TSubjets[j] += it->first;
		    constituents[j].erase(it);
		    --it;
		}
		else break;
	    }
	}

	//  Compute total jet momentum.
	TLorentzVector totalJetMomentum = std::accumulate(TSubjets.begin(), TSubjets.end(), TLorentzVector(0.0, 0.0, 0.0, 0.0));

	//  Store raw data.
	if (totalJetMomentum.M() > 0.01) {
            for (unsigned int k = 0; k < numSubjets; ++k) {
	        //  Store the total jet mass and pT.
                if (k == 0) {
                    telescopingMasses[k][i] = totalJetMomentum.M();
	            telescopingPTs[k][i] = totalJetMomentum.Pt();
                }
                //  Store the subjet mass and pT.
                else {
                    telescopingMasses[k][i]  = TSubjets[k - 1].M();
                    telescopingPTs[k][i] = TSubjets[k - 1].Pt();
                }
            }
        }
    }

    //  Compute and store volatilities.
    if (!telescopingMasses[0].empty()) {
	result.massVolatility = getVolatility(telescopingMasses[0]);
	result.pTVolatility = getVolatility(telescopingPTs[0]);
	result.masses = telescopingMasses;
        result.pTs = telescopingPTs;
    }

    return result;
}

//  Determine if all vectors are empty.
bool TelescopingJets::emptyTelescopingMasses(std::vector<double> totalTelescopingMass, std::vector<std::vector<double>> combinedTelescopingMasses) {
    if (totalTelescopingMass.empty()) return false;

    for (auto const &masses : combinedTelescopingMasses) {
        if (!masses.empty()) return false;
    }
    
    return true;
}

tSub TelescopingJets::telescopeSubjets(unsigned int numSubjets, tSub result, std::vector<double> subjetRadii, std::vector<std::vector<std::pair<TLorentzVector, double>>> constituents, double targetMass) {
    //  Vector of subjet masses.
    std::vector<TLorentzVector> TSubjets(numSubjets);
    
    //  Vectors for mass and pT raw data for jet and subjets.
    std::vector<std::vector<double>> telescopingMasses(numSubjets + 1, std::vector<double>(subjetRadii.size()));
    std::vector<std::vector<double>> telescopingPTs(numSubjets + 1, std::vector<double>(subjetRadii.size()));

    //  Vectors for combinations of sums of masses and pTs.
    unsigned int numCombinations = getChoose(numSubjets, numSubjets - 1) + 1;
    std::vector<std::vector<double>> combinedTelescopingMasses(numCombinations, std::vector<double>(subjetRadii.size()));
    std::vector<std::vector<double>> combinedTelescopingPTs(numCombinations, std::vector<double>(subjetRadii.size()));
    
    for (unsigned int i = 0; i < subjetRadii.size(); ++i) {
        for (unsigned int j = 0; j < constituents.size(); ++j) {
	    for (auto it = constituents[j].begin(); it != constituents[j].end(); ++it) {
		//  Add constituent mass to subjet mass if within range.
		if (it->second <= subjetRadii[i]) {
		    TSubjets[j] += it->first;
		    constituents[j].erase(it);
		    --it;
		}
		else break;
	    }
	}

        //  Compute total jet momentum.
        TLorentzVector totalJetMomentum = std::accumulate(TSubjets.begin(), TSubjets.end(), TLorentzVector(0.0, 0.0, 0.0, 0.0));

        if (totalJetMomentum.M() > 0.01) {
            //  Store raw data.
            for (unsigned int k = 0; k < numSubjets; ++k) {
                //  Store the total jet mass and pT.
                if (k == 0) {
                    telescopingMasses[k][i] = totalJetMomentum.M();
                    telescopingMasses[k][i] = totalJetMomentum.Pt();
                }
                //  Store the subjet mass and pT.
                else {
                    telescopingMasses[k][i] = TSubjets[k - 1].M();
                    telescopingPTs[k][i] = TSubjets[k - 1].Pt();
                }
            }
            
            //  Calculate combinations of masses and pTs.
            //  TODO Add functionality for different number of subjets
            //  instead of just one less jet than the total. Also add
            //  for more than one target mass.
            for (unsigned int k = 0; k < numSubjets; ++ k) {
                TLorentzVector partialJetMomentum = totalJetMomentum - TSubjets[k];
                combinedTelescopingMasses[k][i] = partialJetMomentum.M();
                combinedTelescopingPTs[k][i] = partialJetMomentum.Pt();
            }
        }
    }

    if (!emptyTelescopingMasses(telescopingMasses[0], combinedTelescopingMasses)) {
        result.massVolatility = getVolatility(telescopingMasses[0]);
        result.pTVolatility = getVolatility(telescopingPTs[0]);
        result.masses = telescopingMasses;
        result.pTs = telescopingPTs;

        std::vector<double> residualMass(combinedTelescopingMasses.size());
        std::transform(combinedTelescopingMasses.begin(), combinedTelescopingMasses.end(), std::back_inserter(residualMass), 
                [](std::vector<double> combinedTelescopingMass) -> double { return combinedTelescopingMass.back(); });

        //  Subtract the target mass from all elements in
        //  the masses vector except for the first.
        std::transform(residualMass.begin(), residualMass.end(), residualMass.begin(),
                       bind2nd(std::minus<double>(), targetMass));

        //  Get the distance from the target mass by 
        //  taking the absolute value of the elements
        //  from which we subtracted the target mass.
        std::transform(residualMass.begin(), residualMass.end(), residualMass.begin(),
                       static_cast<double (*)(double)>(&std::abs));
        
        //  Check THIS TODO
        auto targetMassPredictionIt = std::min_element(residualMass.cbegin(), residualMass.cend());
        unsigned int targetMassPredictionIndex = std::distance(residualMass.cbegin(), targetMassPredictionIt);

        result.targetMass = combinedTelescopingMasses[targetMassPredictionIndex].back();
        result.targetMassVolatility = getVolatility(combinedTelescopingMasses[targetMassPredictionIndex]);
        result.targetPTVolatility = getVolatility(combinedTelescopingPTs[targetMassPredictionIndex]);
    }

    return result;
}

tSub TelescopingJets::tNSubjet(unsigned int numSubjets, double minRadius, double maxRadius, unsigned int numRadii, double targetMass = 0.0) {
    tSub result;

    double beta = 1.0;
    std::vector<TLorentzVector> pTauAxes = getTauAxes(numSubjets, beta);
    std::vector<double> anglesBetweenTauAxes = getAnglesBetweenTauAxes(numSubjets, pTauAxes);

    //  Store subjet axes.
    result.tauAxes = pTauAxes;

    //  Collect the minimum angle if there is one.
    if (!anglesBetweenTauAxes.empty()) result.minAngle = *(std::min_element(anglesBetweenTauAxes.cbegin(), anglesBetweenTauAxes.cend()));

    //  Get constituents sorted by subjet and increasing distance from subjet center.
    std::vector<std::vector<std::pair<TLorentzVector, double>>> constituents = sortConstituents(numSubjets, pTauAxes);

    std::vector<double> subjetRadii = getTelescopingParameterSet(minRadius, maxRadius, numRadii);
    
    //  Telescope through the jet radius.
    if (targetMass == 0) result = telescopeSubjets(numSubjets, result, subjetRadii, constituents);
    else result = telescopeSubjets(numSubjets, result, subjetRadii, constituents, targetMass);

    if (result.massVolatility == -1) {
	std::cout << "WARNING zero entries for TNSubjet!   numSubjets: " << numSubjets <<
	    "   minRadius: " << minRadius <<
	    "   maxRadius: " << maxRadius <<
	    "   numRadii: " << numRadii <<
	    "   input mass: " << input.m() << std::endl;
    }       

    return result;
}

///=========================================
/// Telescoping N-subjettiness
///=========================================
double TelescopingJets::tNsubjettiness(int numSubjets, double minBeta, double maxBeta, unsigned int numBetas) {
    std::vector<double> betas = getTelescopingParameterSet(minBeta, maxBeta, numBetas);
    std::vector<double> telescopingTaus(betas.size());

    for (unsigned int i = 0; i < betas.size(); ++i) {
        fastjet::contrib::UnnormalizedMeasure nSubMeasure(betas[i]);
        fastjet::contrib::Nsubjettiness nSub(numSubjets, *axesDefinition, nSubMeasure);
        telescopingTaus[i] = nSub(input);
    }

    if (!telescopingTaus.empty()) return getVolatility(telescopingTaus);
    
    std::cout <<"WARNING zero entries for T_Nsubjettiness! minBeta: " << minBeta 
        << "\tmaxBeta: "<< maxBeta << 
        "\tnumBetas: " << numBetas << std::endl;
    return -1;
}

double TelescopingJets::tNsubjettinessRatio(int nNumerator, int nDemoninator, double minBeta, double maxBeta, unsigned int numBetas) {
    std::vector<double> betas = getTelescopingParameterSet(minBeta, maxBeta, numBetas);
    std::vector<double> telescopingTaus(betas.size());

    for (unsigned int i = 0; i < betas.size(); ++i) {
        fastjet::contrib::UnnormalizedMeasure nsubMeasure(betas[i]);
        fastjet::contrib::Nsubjettiness nSubNumerator(nNumerator, *axesDefinition, nsubMeasure);
        fastjet::contrib::Nsubjettiness nSubDenominator(nDemoninator, *axesDefinition, nsubMeasure);
      
        double numerator = nSubNumerator(input);
        double denominator = nSubDenominator(input);
      
        if (denominator != 0) telescopingTaus[i] = (numerator / denominator);
        else telescopingTaus[i] = -1.0;
    }
    
    if (!telescopingTaus.empty() || std::accumulate(telescopingTaus.begin(), telescopingTaus.end(), 0.0) > 0) return getVolatility(telescopingTaus);
  
    std::cout << "WARNING zero entries for T_NsubjettinessRatio! minBeta: " << minBeta <<
        "\tmaxBeta: " << maxBeta <<
        "\tnumBetas: " << numBetas <<  std::endl;
    return -1;
}


///=========================================
/// Telescoping Energy Correlators
///=========================================
double TelescopingJets::tEnergyCorrelator_C2 (double minBeta, double maxBeta, unsigned int numBetas) {
    std::vector<double> betas = getTelescopingParameterSet(minBeta, maxBeta, numBetas);
    std::vector<double> telescopingEcfs(betas.size());

    for (unsigned int i = 0; i < betas.size(); ++i) {
	fastjet::contrib::EnergyCorrelatorC2 ecf(betas[i]);
	telescopingEcfs[i] = ecf(input);
    }

    if(!telescopingEcfs.empty()) return getVolatility(telescopingEcfs);
    
    std::cout << "WARNING zero entries for T_EnergyCorrelator_C2! minBeta: "<< minBeta <<
        "\tmaxBeta: " << maxBeta <<
        "\tnumBetas: " << numBetas <<  std::endl;
    return -1;   
}

double TelescopingJets::tEnergyCorrelator_D2(double minBeta, double maxBeta, unsigned int numBetas) {
    std::vector<double> betas = getTelescopingParameterSet(minBeta, maxBeta, numBetas);
    std::vector<double> telescopingEcfs(betas.size());

    for (unsigned int i = 0; i < betas.size(); i++) {
	fastjet::contrib::EnergyCorrelatorD2 ecf(betas[i]);
	telescopingEcfs[i] = ecf(input);
    }
    
    if(!telescopingEcfs.empty()) return getVolatility(telescopingEcfs);
    
    std::cout <<  "WARNING zero entries for T_EnergyCorrelator_C2! minBeta: " << minBeta <<
        "\tmaxBeta: " << maxBeta <<
        "\tnumBetas: " << numBetas <<  std::endl;
    return -1;   
}

double TelescopingJets::tEnergyCorrelator_C3(double minBeta, double maxBeta, unsigned int numBetas) {
    std::vector<double> betas = getTelescopingParameterSet(minBeta, maxBeta, numBetas);
    std::vector<double> telescopingEcfs(betas.size());

    for (unsigned int i = 0; i < betas.size(); ++i) {
	fastjet::contrib::EnergyCorrelatorDoubleRatio ecf(3, betas[i]);
	telescopingEcfs[i] = ecf(input);
    }
    
    if (telescopingEcfs.empty()) return getVolatility(telescopingEcfs);
    
    std::cout <<"WARNING zero entries for T_EnergyCorrelator_C3! minBeta: " << minBeta <<
        "\tmaxBeta: " << maxBeta <<
        "\tnumBetas: " << numBetas <<  std::endl;
    return -1;
}
