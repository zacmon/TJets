//  Telescoping Jets Package
//
//  Alexander Emerman, Yang-Ting Chien, Shih-Chieh Hsu, Zachary Montague
//

#include "TelescopingJets.hh"

TelescopingJets::TelescopingJets(const fastjet::PseudoJet& input) {
    this->input = input;
}

std::vector<TLorentzVector> TelescopingJets::convertPseudoJet2TLV(std::vector<fastjet::PseudoJet> pseudoJet) {
    std::vector<TLorentzVector> tLorentzVector(pseudoJet.size());

    for (unsigned i = 0; i < pseudoJet.size(); i++) {
	tLorentzVector[i].SetPxPyPzE(pseudoJet[i].px(), pseudoJet[i].py(), pseudoJet[i].pz(), pseudoJet[i].e());
    }

    return tLorentzVector;
}

std::vector<double> TelescopingJets::getTelescopingParameterSet(double minParameter, double maxParameter, int numParameter, int stepScale = 0) {
    double deltaParameter = -999;
    if (stepScale == 0) deltaParameter = (maxParameter - minParameter) / (numParameter);
    else if (stepScale == 1) deltaParameter = log10(maxParameter / minParameter) / (numParameter);
    
    std::vector<double> parameterSet(numParameter);
    
    for (int i = 0; i <= numParameter; ++i) {
	if (stepScale == 0) parameterSet[i] = minParameter + i * deltaParameter;
	else if (stepScale == 1) parameterSet[i] = minParameter * pow(10, i * deltaParamter);
    }

    return parameterSet;
}

// Useful mathematical functions to find the average, rms and volatility
// of a set of values
double getAverage(const std::vector<double>& values) {
    if (values.empty()) throw std::length_error("Asking for average of empty vector.\n");
    double v = 0;
    for (unsigned int i=0; i < values.size(); i++) {
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
    return getRMS(values)/getAverage(values);
}
    
///=========================================
/// Telescoping Pruning
///=========================================
double tPruning(double minDCut, double maxDCut, int numDCuts) {
    //  Single choice of zcut but can be further telescoped.
    double zCut = 0.1;
    
    std::vector<double> telescopingMasses;

    std::vector<double> dCuts = getTelescopingParameterSet(minDCut, maxDCut, numDCuts);

    //  Obtain a set of jet masses by stepping through the DCut
    //  pruning parameter.
    for (unsigned int i = 0; i <= dCuts.size(); ++i) {
        fastjet::Pruner pruner(fastjet::cambridge_algorithm, zCut, dCuts[i]);
        fastjet::PseudoJet prunedJet = pruner(input);
        double mass = prunedJet.m();
        if (mass > M0) telescopingMasses.push_back(mass);
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
double tTrimming(double minFCut, double maxFCut, int numFCuts) {
    //  The subjet radius used should depend on the input jet's
    //  pT--higher pT should have a smaller subjet radius.
    double Rfilt = 0.1; 
    if (input.pt() > 500) Rfilt = 0.05;
    
    std::vector<double> telescopingMasses;

    std::vector<double> fCuts = getTelescopingParameterSet(minFCut, maxFCut, numFCuts);
    
    //  Obtain a set of jet masses by stepping through the f trimming parameter.
    for (unsigned int i = 0; i <= fCuts.size(); ++i) {
        fastjet::Filter trimmer(Rfilt, fastjet::SelectorPtFractionMin(fCuts[i]));
        fastjet::PseudoJet trimmedJet = trimmer(input);
        double mass = trimmedJet.m();
        if (mass > M0) telescopingMasses.push_back(mass);
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
double tReclustering(int algorithm, double minRadius, double maxRadius, int numRadii) {
    std::vector<double> telescopingMass; 

    std::vector<double> subjetRadii = getTelescopingParameterSet(minRadius, maxRadius, numRadii, stepScale);

    for (unsigned int i = 0; i <= subjetRadius.size(); ++i) {
        fastjet::JetDefinition TjetDefinition(fastjet::JetAlgorithm(algorithm), subjetRadii[i]);
        fastjet::ClusterSequence TClusterSequence(input.constituents(), TjetDefinition);
        std::vector<fastjet::PseudoJet> recoTJets = sorted_by_pt(TClusterSequence.inclusive_jets());
	
        if (recoTJets.empty()) {
            std::cout << "WARNING! Recluster number of subjets is " << recoTJets.size() << std::endl;
            continue;
        }
        if (recoTJets.size() == 1) {
            Tsubjet1.SetPxPyPzE(recoTJets[0].px(), recoTJets[0].py(), recoTJets[0].pz(), recoTJets[0].e());
            double mass = Tsubjet1.M();
            if (mass > M0) telescopingMass.push_back(mass);
        }
        else if (recoTJets.size() >= 2) {
            Tsubjet1.SetPxPyPzE(recoTJets[0].px(), recoTJets[0].py(), recoTJets[0].pz(), recoTJets[0].e());
            Tsubjet2.SetPxPyPzE(recoTJets[1].px(), recoTJets[1].py(), recoTJets[1].pz(), recoTJets[1].e());
            double mass = (Tsubjet1 + Tsubjet2).M();
            if (mass > M0) telescopingMass.push_back(mass);
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

std::vector<TLorentzVector> getTauAxes(unsigned int numSubjets, double beta, int axesType) {   
    //  Create unnormalized measure.
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);

    //  TODO
    //  Better way to differentiate between axes?
    //  Create nSubjettiness axes.
    if (axesType == 0) fastjet::contrib::Nsubjettiness nSubjettiness(numSubjets, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    else if (axesType == 1) fastjet::contrib::Nsubjettiness nSubjettiness(numSubjets, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);
    else std::cout << "Error! Invalid axesType entered. Please enter \"0\" for OnePass_KT_Axes or \"1\" for WTA_KT_Axes." << std::endl;

    //  Get tau axes.
    double tauN = nSubjettiness.result(input);
    std::vector<fastjet::PseudoJet> tauAxes = nSubjettiness.currentAxes();

    std::vector<TLorentzVector> pTauAxes = convertPseudoJet2TLV(tauAxes);
    return pTauAxes;
}

std::vector<double> getAnglesBetweenTauAxes(unsigned int numSubjets, std::vector<TLorentzVector> pTauAxes) {
    std::vector<double> angles;

    //  Calculate deltaR between each pair of tau axes.
    for (unsigned int i = 0; i < numSubjets; ++i) {
	for (unsigned int j = i + 1; j < numSubjets; ++j) {
	    angles.push_back(pTauAxes[i].DeltaR(pTauAxes[j]));
	}
    }

    return angles;
}

std::vector<std::vector<TLorentzVector>> sortConstituents(unsigned int numSubjets, std::vector<TLorentzVector> pTauAxes) {
    //  Get constituents of the input jet.
    std::vector<fastjet::PseudoJet> constituents = input.constituents();

    //  Initialize container that will hold constiuents sorted by
    //  their distance to some tau axis and then sorted by increasing distance.
    std::vector<std::vector<std::pair<TLorentzVector, double>>> tempConstituents(numSubjets);

    //  Vector to be returned without the measure to tau axes.
    std::vector<std::vector<TLorentzVector>> sortedConstituents(numSubjets);
    
    for (auto const &constituent: constituents) {
	//  Convert PseudoJet four-vector to TLorentz four-vector.
	TLorentzVector pConstituent(constituent.px(), constituent.py(), constituent.pz(), constituent.e());

	//  Get distances to different tau axes
	std::vector<double> distanceToTauAxes(numSubjets);
	for (unsigned int i = 0; i < numSubjets; ++i) {
	    distanceToTauAxes[i] = pConstituent.DeltaR(pTauAxes[i]);
	}

	//  Identify tau axes associated with minium distance
	auto minDistanceToTauAxesIt = std::min_element(distanceToTauAxes.cbegin(), distanceToTauAxes.cend());
	double minDistanceToTauAxes = *minDistanceToTauAxesIt;
	unsigned int minDistanceToTauAxesIndex = std::distance(distanceToTauAxes.cbegin(), minDistanceToTauAxesIt);
	
	tempConstituents[minDistanceToTauAxesIndex].push_back({pConstituent, minDistanceToTauAxes});
    }

    //  Sort constituents by increasing distance for each subjet.
    for (int i = 0; i < numSubjets; ++i) {	
	std::sort(tempConstituents[i].begin(), tempConstituents[i].end(), [](const std::pair<TLorentzVector, double> &pair1,
									   const std::pair<TLorentzVector, double> &pair2) {
		      return pair1.second < pair2.second;
		  });

	//  Store only the constituents' four-vectors in the vector
	//  to be returned.
	std::transform(tempConstituents[i].begin(), tempConstituents[i].end(), std::back_inserter(sortedConstituents),
	   (const TLorentzVector& (*)(const std::pair<TLorentzVector, double>&))std::get<0>);
    }
    
    return sortedConstituents;
}



tSub telescopeSubjets(usigned int numSubjets, tSub result, std::vector<double> subjetRadii, std::vector<std::vector<TLorentzVector>> constituents) {
    //  Vector of subjet masses.
    std::vector<TLorentzVector> TSubjets(numSubjets);

    //  Vectors for mass and pT raw data.
    std::vector<double> telescopingMasses;
    std::vector<double> telescopingPT;

    for (int i = 0; i <= subjetRadius.size(); ++i) {
	for (unsigned int j = 0; j < constituents.size(); ++j) {
	    for (auto it = constituents[j].begin(); it != constituents[j].end(); ++it) {
		//  Add constituent mass to subjet mass if within range.
		if (it->second <= subjetRadius[i]) {
		    TSubjets[j] += it->first;
		    constituents[j].erase(it);
		    --it;
		}
		else break;
	    }
	}

	//  Compute total jet mass.
	TLorentzVector sumTSubjet;
	for (auto const &TSubjet : TSubjets) {
	    sumTSubjet += TSubjet;
	}

	//  Store raw data.
	if (sumTSubjet.M() > M0) {
	    telescopingMasses.push_back(sumTSubjet.M());
	    telescopingPT.push_back(sumTSubjet.Pt());
	}
    }

    //  Compute and store volatilities.
    if (!telescopingMasses.empty()) {
	result.massVolatility = getVolatility(telescopingMasses);
	result.pTVolatility = getVolatility(telescopingPT);
	result.masses = telescopingMasses;
    }
    else {
	std::cout << "WARNING zero entries for TNSubjet!   numSubjets: " << numSubjets <<
	    "   minRadius: " << minRadius <<
	    "   maxRadius: " << maxRadius <<
	    "   numRadii: " << numRadii <<
	    "   input mass: " << input.m() << std::endl;
    }
}

tSub telescopeMass(unsigned int numSubjets, tSub result, std::vector<double> subjetRadii, std::vector<std::vector<TLorentzVector>> constituents, double targetMass) {
    //  Vector of subjet masses.
    std::vector<TLorentzVector> TSubjets(numSubjets);

    //  TODO: think of better name
    //  TODO: implement combination instead of +1
    unsigned int A = numSubjets + 1;
    if (numSubjets == 1) A = 1;
    
    //  Vectors for mass and pT raw data.
    std::vector<double> telescopingMasses(A);
    std::vector<double> telescopingPTs(A);

    //  Vectors for combinations of sums of masses and pTs.
    std::vector<double> masses(A);
    std::vector<double> pTs(A);

    for (unsigned int i = 0; i <= subjetRadius.size(); ++i) {
	for (unsigned int j = 0; j < constituents.size(); ++j) {
	    for (auto it = constituents[j].begin(); it != constituents[j].end(); ++it) {
		//  Add constituent mass to subjet mass if within range.
		if (it->second <= subjetRadius[i]) {
		    TSubjets[j] += it->first;
		    constituents[j].erase(it);
		    --it;
		}
		else break;
	    }
	}

	//  Compute total jet mass and pT.
	for (auto const &TSubjet : TSubjets) {
	    sumTSubjet += TSubjet;
	}
	masses[0] = sumTSubjet.M();
	pTs[0] = sumTSubjet.Pt();

	//  TODO: adjust for different amount of partial sums
	//  Compute partial jet mass and PT.
	for (unsigned int k = 0; k < A; ++k) {
	    

	//  Store raw data.
	if (sumTSubjet.M() > M0) {
	    telescopingMasses.push_back(sumTSubjet.M());
	    telescopingPT.push_back(sumTSubjet.Pt());
	}
    }

    //  Compute and store volatilities.
    if (!telescopingMasses.empty()) {
	result.massVolatility = getVolatility(telescopingMasses);
	result.pTVolatility = getVolatility(telescopingPT);
	result.masses = telescopingMasses;
    }
    else {
	std::cout << "WARNING zero entries for TNSubjet!   numSubjets: " << numSubjets <<
	    "   minRadius: " << minRadius <<
	    "   maxRadius: " << maxRadius <<
	    "   numRadii: " << numRadii <<
	    "   input mass: " << input.m() << std::endl;
    }
}

tSub tNSubjet(unsigned int numSubjets, double minRadius, double maxRadius, int numRadii, int stepScale, int axesType = 0, double targetMass = 0.0) {
    tSub result;

    double beta = 1.0;
    std::vector<TLorentzVector> pTauAxes = getTauAxes(numSubjets, beta, axesType);
    std::vector<double> anglesBetweenTauAxes = getAnglesBetweenTauAxes(pTauAxes);

    //  Collect the minimum angle if there is one.
    if (!anglesBetweenTauAxes.empty()) result.minAngle = *(std::min_element(anglesBetweenTauAxes.cbegin(), anglesBetweenTauAxes.cend()));

    //  Get constituents sorted by subjet and increasing distance from subjet center.
    std::vector<std::vector<TLorentzVector>> constituents = sortConstituents(pTauAxes);

    std::vector<double> subjetRadii = getTelescopingParameterSet(minRadius, maxRadius, numRadii, stepScale);
    
    //  Telescope through the jet radius.
    if (targetMass == 0) telescopeSubjets(result, subjetRadii, constituents);
    else telescopeSubjets(result, subjetRadii, constituents, targetMass);
    
    return result;
}


T2Sub T2Subjet(double minRadius, double maxRadius, int numRadii, int stepScale) {
    T2Sub result;
    
    double beta = 1.0;
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
    //fastjet::contrib::Nsubjettiness nSubjettiness(2, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    fastjet::contrib::Nsubjettiness nSubjettiness(2, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);

    double tau2 = nSubjettiness.result(input);
    std::vector<fastjet::PseudoJet> tau2axes = nSubjettiness.currentAxes();

    std::vector<TLorentzVector> pTauAxes(2);
    for (unsigned int i = 0; i < 2; ++i) {
	pTauAxes[i].SetPxPyPzE(tau2axes[i].px(), tau2axes[i].py(), tau2axes[i].pz(), tau2axes[i].e());
    }

    result.minAngle = pTauAxes[0].DeltaR(pTauAxes[1]);

    std::vector<fastjet::PseudoJet> constituents = input.constituents();
    std::vector<std::vector<std::pair<TLorentzVector, double>>> sortedConstituents(2);

    for (auto const &constituent : constituents) {
	TLorentzVector pConstituent(constituent.px(), constituent.py(), constituent.pz(), constituent.e());

	std::vector<double> distanceToTauAxes(2);
	for (unsigned int i = 0; i < 2; ++i) {
	    distanceToTauAxes[i] = pConstituent.DeltaR(pTauAxes[i]);
	}

	auto minDistanceToTauAxesIt = std::min_element(distanceToTauAxes.cbegin(), distanceToTauAxes.cend());
	double minDistanceToTauAxes = *minDistanceToTauAxesIt;
	unsigned int minDistanceToTauAxesIndex = std::distance(distanceToTauAxes.cbegin(), minDistanceToTauAxesIt);

	sortedConstituents[minDistanceToTauAxesIndex].push_back({pConstituent, minDistanceToTauAxes});
    }
    
    for (auto &subjetConstituents : sortedConstituents) {
	std::sort(subjetConstituents.begin(), subjetConstituents.end(), [](const std::pair<TLorentzVector, double> &pair1,
									   const std::pair<TLorentzVector, double> &pair2) {
		      return pair1.second < pair2.second;
		  });
    }
    
    std::vector<TLorentzVector> TSubjets(2);
    std::vector<std::vector<double>> telescopingMasses(3);
    std::vector<std::vector<double>> telescopingPTs(3);
    std::vector<double> masses(3);
    std::vector<double> pTs(3);

    double deltaR = -999;
    if (stepScale == 0) deltaR = (maxRadius - minRadius) / (numRadii);
    else if (stepScale == 1) deltaR = log10(maxRadius / minRadius) / (numRadii);

    for (int i = 0; i <= numRadii; ++i) {
	double r = -999;
	if (stepScale == 0) r = minRadius + i * deltaR;
	else if (stepScale == 1) r = minRadius * pow(10, i * deltaR);

	for (unsigned int j = 0; j < sortedConstituents.size(); ++j) {
	    for (auto it = sortedConstituents[j].begin(); it != sortedConstituents[j].end(); ++it) {
		if (it->second <= r) {
		    TSubjets[j] += it->first;
		    sortedConstituents[j].erase(it);
		    --it;
		}
		else break;
	    }
	}
	masses[0] = (TSubjets[0] + TSubjets[1]).M();
	masses[1] = TSubjets[0].M();
	masses[2] = TSubjets[1].M();

	pTs[0] = (TSubjets[0] + TSubjets[1]).Pt();
	pTs[1] = TSubjets[0].Pt();
	pTs[2] = TSubjets[1].Pt();

	for (unsigned int i = 0; i < masses.size(); ++i) {
	    telescopingMasses[i].push_back(masses[i]);
	    telescopingPTs[i].push_back(pTs[i]);
	}
    }
    
    if (!telescopingMasses[0].empty() && !telescopingMasses[1].empty() && !telescopingMasses[2].empty()) {
	result.massVolatility = getVolatility(telescopingMasses[0]);
	result.pTVolatility = getVolatility(telescopingPTs[0]);
	result.masses = telescopingMasses[0];

	std::vector<double> residualWMass = masses;

	std::transform(residualWMass.begin() + 1, residualWMass.end(), residualWMass.begin() + 1,
		       bind2nd(std::minus<double>(), massW));
	std::transform(residualWMass.begin() + 1, residualWMass.end(), residualWMass.begin() + 1,
		       static_cast<double (*)(double)>(&std::abs));

	auto wMassPredictionIt = std::min_element(residualWMass.cbegin() + 1, residualWMass.cend());
	unsigned int wMassPredictionIndex = std::distance(residualWMass.cbegin(), wMassPredictionIt);

	result.wMass = masses[wMassPredictionIndex];
	result.wMassVolatility = getVolatility(telescopingMasses[wMassPredictionIndex]);
    }

    else {
	std::cout << "WARNING zero entries for T2Subjet!   minRadius: " << minRadius <<
	    "   maxRadius: " << maxRadius <<
	    "   numRadii: " << numRadii <<
	    "   input mass: " << input.m() << std::endl;
    }
    return result;
}

T3Sub T3Subjet(double minRadius, double maxRadius, int numRadii, int stepScale) {
    T3Sub result;

    double beta = 1.0;
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
    //fastjet::contrib::Nsubjettiness nSubjettiness(3, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    fastjet::contrib::Nsubjettiness nSubjettiness(3, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);
    double tau3 = nSubjettiness.result(input);
    std::vector<fastjet::PseudoJet> tau3axes = nSubjettiness.currentAxes();

    std::vector<TLorentzVector> pTauAxes(3);
    for (unsigned int i = 0; i < 3; ++i) {
	pTauAxes[i].SetPxPyPzE(tau3axes[i].px(), tau3axes[i].py(), tau3axes[i].pz(), tau3axes[i].e());
    }
    
    std::vector<double> anglesBetweenTauAxes;
    for (unsigned int i = 0; i < 3; ++i) {
	for (unsigned int j = i + 1; j < 3; ++j) {
	    anglesBetweenTauAxes.push_back(pTauAxes[i].DeltaR(pTauAxes[j]));
	}
    }

    std::vector<double> tempAnglesBetweenTauAxes = anglesBetweenTauAxes;
    std::sort(tempAnglesBetweenTauAxes.begin(), tempAnglesBetweenTauAxes.end());

    result.minAngle = tempAnglesBetweenTauAxes[0];
    result.midAngle = tempAnglesBetweenTauAxes[1];
    result.maxAngle = tempAnglesBetweenTauAxes[2];

    std::vector<fastjet::PseudoJet> constituents = input.constituents();
    std::vector<std::vector<std::pair<TLorentzVector, double>>> sortedConstituents(3);

    for (auto const &constituent: constituents) {
        TLorentzVector pConstituent(constituent.px(), constituent.py(), constituent.pz(), constituent.e());

        std::vector<double> distanceToTauAxes(3);
        for (unsigned int i = 0; i < 3; ++i) {
            distanceToTauAxes[i] = pConstituent.DeltaR(pTauAxes[i]);
        }

	auto minDistanceToTauAxesIt = std::min_element(distanceToTauAxes.cbegin(), distanceToTauAxes.cend());
        double minDistanceToTauAxes = *minDistanceToTauAxesIt;
        unsigned int minDistanceToTauAxesIndex = std::distance(distanceToTauAxes.cbegin(), minDistanceToTauAxesIt);

        sortedConstituents[minDistanceToTauAxesIndex].push_back({pConstituent, minDistanceToTauAxes});
    }

    for (auto &subjetConstituents : sortedConstituents) {
	std::sort(subjetConstituents.begin(), subjetConstituents.end(), [](const std::pair<TLorentzVector, double> &pair1, const std::pair<TLorentzVector, double> &pair2) {
                return pair1.second < pair2.second;
            });
    }

    std::vector<TLorentzVector> TSubjets(3);
    std::vector<std::vector<double>> telescopingMasses(4);
    std::vector<std::vector<double>> telescopingPTs(4);
    std::vector<double> masses(4);
    std::vector<double> pTs(4);

    double deltaR = -999;
    if (stepScale == 0) deltaR = (maxRadius - minRadius) / (numRadii);
    else if (stepScale == 1) deltaR = log10(maxRadius / minRadius) / (numRadii);
    
    for (int i = 0; i <= numRadii; ++i) {
	double r = -999;
	if (stepScale == 0) r = minRadius + i * deltaR;
	else if (stepScale == 1) r = minRadius * pow(10, i * deltaR);

	for (unsigned int j = 0; j < sortedConstituents.size(); ++j) {
	    for (auto it = sortedConstituents[j].begin(); it != sortedConstituents[j].end(); ++it) {
		if (it->second <= r) {
		    TSubjets[j] += it->first;
		    sortedConstituents[j].erase(it);
		    --it;
		}
		else break;
	    }
	}
	
	masses[0] = (TSubjets[0] + TSubjets[1] + TSubjets[2]).M();
	masses[1] = (TSubjets[0] + TSubjets[1]).M();
	masses[2] = (TSubjets[0] + TSubjets[2]).M();
	masses[3] = (TSubjets[1] + TSubjets[2]).M();

	pTs[0] = (TSubjets[0] + TSubjets[1] + TSubjets[2]).Pt();
	pTs[1] = (TSubjets[0] + TSubjets[1]).Pt();
	pTs[2] = (TSubjets[0] + TSubjets[2]).Pt();
	pTs[3] = (TSubjets[1] + TSubjets[2]).Pt();


	for (unsigned int i = 0; i < masses.size(); ++i) {
	    telescopingMasses[i].push_back(masses[i]);
	    telescopingPTs[i].push_back(pTs[i]);
	}
    }
    
    if (!telescopingMasses[0].empty() && !telescopingMasses[1].empty() && !telescopingMasses[2].empty() && !telescopingMasses[3].empty()) {
        result.massVolatility = getVolatility(telescopingMasses[0]);
	result.pTVolatility = getVolatility(telescopingPTs[0]);
	result.masses = telescopingMasses[0];

	//  TODO
	//  Is there a better way to recover W masses?
	std::vector<double> residualWMass = masses;
	
	std::transform(residualWMass.begin() + 1, residualWMass.end(), residualWMass.begin() + 1,
		  bind2nd(std::minus<double>(), massW));
	std::transform(residualWMass.begin() + 1, residualWMass.end(), residualWMass.begin() + 1,
		       static_cast<double (*)(double)>(&std::abs));

	auto wMassPredictionIt = std::min_element(residualWMass.cbegin() + 1, residualWMass.cend());
	unsigned int wMassPredictionIndex = std::distance(residualWMass.cbegin(), wMassPredictionIt);

	result.wMass = masses[wMassPredictionIndex];
	result.wMassVolatility = getVolatility(telescopingMasses[wMassPredictionIndex]);
	result.wPTVolatility = getVolatility(telescopingPTs[wMassPredictionIndex]);
    }

    else {
        std::cout << "WARNING zero entries for T3Subjet!   minRadius: " << minRadius <<
	    "   maxRadius: " << maxRadius <<
	    "   numRadii: " << numRadii <<
	    "   input mass: " << input.m() << std::endl;
    }    
    return result;
}

///=========================================
/// Telescoping N-subjettiness
///=========================================
double T_Nsubjettiness(int numSubjets, double minBeta, double maxBeta, int numBetas) {
  std::vector<double> taus;
  double betaStep = (maxBeta - minBeta) / (numBetas);
  
  for (double beta = minBeta; beta < maxBeta; beta += betaStep) {
      fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
      fastjet::contrib::Nsubjettiness nSub(N, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
      //fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
      taus.push_back(nSub(input));
  }
  if (!taus.empty()) return getVolatility(taus);
  
  std::cout <<"WARNING zero entries for T_Nsubjettiness! minBeta: " << minBeta << "\tmaxBeta: "<< maxBeta << "\tnumBetas: " << numBetas << std::endl;
  return -1;
}

double T_NsubjettinessRatio(int nNumerator, int nDemoninator, fastjet::PseudoJet& input, double minBeta, double maxBeta, int numBetas) {
  std::vector<double> taus;
  double betaStep = (maxBeta - minBeta) / (numBetas);
  
  for (double beta = minBeta; beta < maxBeta + betaStep; beta += betaStep) {
      fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
      
      fastjet::contrib::Nsubjettiness nSubNumerator(nNumerator, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
      fastjet::contrib::Nsubjettiness nSubDenominator(nDemoninator, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
      
      double numerator = nSubNumerator(input);
      double denominator = nSubDenominator(input);
      
      if (denominator != 0) taus.push_back(numerator / denominator);
      else taus.push_back(-1.0);
  }
  
  if (!taus.empty()) return getVolatility(taus);
  
  std::cout <<"WARNING: zero entries for T_NsubjettinessRatio "<< minBeta <<" "<< maxBeta <<" "<< numBetas <<" "<<  std::endl;
  return -1;
}


///=========================================
/// Telescoping Energy Correlators
///=========================================
double T_EnergyCorrelator_C2(fastjet::PseudoJet& input, double minBeta, double maxBeta, int numBetas) {
    std::vector<double> ecfs; ecfs.clear();
    for (int i = 0; i < numBetas; i++) {
	double beta = minBeta + i*(maxBeta - minBeta)/(numBetas-1);
	fastjet::contrib::EnergyCorrelatorC2 ecf(beta);
	ecfs.push_back(ecf(input));
    }
    // getVolatility function provided by TelescopingJets
    if(ecfs.size()>0)
        return getVolatility(ecfs);
    
    std::cout <<"WARNING: zero entries for T_EnergyCorrelator_C2 "<< minBeta <<" "<< maxBeta <<" "<< numBetas <<" "<<  std::endl;
    return -1;   
}

double T_EnergyCorrelator_D2(fastjet::PseudoJet& input, double minBeta, double maxBeta, int numBetas) {
    std::vector<double> ecfs; ecfs.clear();
    for (int i = 0; i < numBetas; i++) {
	double beta = minBeta + i*(maxBeta - minBeta)/(numBetas-1);
	fastjet::contrib::EnergyCorrelatorD2 ecf(beta);
	ecfs.push_back(ecf(input));
    }
    // getVolatility function provided by TelescopingJets
    if(ecfs.size()>0)
        return getVolatility(ecfs);
    
    std::cout <<"WARNING: zero entries for T_EnergyCorrelator_C2 "<< minBeta <<" "<< maxBeta <<" "<< numBetas <<" "<<  std::endl;
    return -1;   
}

double T_EnergyCorrelator_C3(fastjet::PseudoJet& input, double minBeta, double maxBeta, int numBetas) {
    std::vector<double> ecfs; ecfs.clear();
    for (int i = 0; i < numBetas; i++) {
	double beta = minBeta + i*(maxBeta - minBeta)/(numBetas-1);
	fastjet::contrib::EnergyCorrelatorDoubleRatio ecf(3, beta);
	ecfs.push_back(ecf(input));
    }
    // getVolatility function provided by TelescopingJets
    if(ecfs.size()>0)
        return getVolatility(ecfs);
    
    std::cout <<"WARNING: zero entries for T_EnergyCorrelator_C3 "<< minBeta <<" "<< maxBeta <<" "<< numBetas <<" "<<  std::endl;
    return -1;
}


