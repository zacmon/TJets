/*----------------------------------------------------------------------
  
Title : NTupler.cc

DESCRIPTION : Takes input final state particle ntuple from Ana_EventGeneration
and outputs a flat ntuple that contains jet properties created via fastjet
in this code to be used for post-analysis using Ana_MiniNTupleAnalysis. (NOTE!!)
All of the TelescopingJets code from the fastjet contrib should be contained within this code

COMPILE :
$ source compile.sh

RUN :
$ ./NTupler <type> <input> <output>

type   : 0 = dijet , 1 = G*->W+W- , 2 = ttbar
input  : Input file from Ana_MiniNTupleAnalysis
output : Anything you want - but being logical

------------------------------------------------------------------------*/

#include "NTupler.h"

int main (int argc, char* argv[]) {
    //  Exit if you don't pass a run card.
    if (argc < 4) {
	
	std::cout << "You need to specify more arguments" << std::endl;
	std::cout << "Arg1 = process type" << std::endl;
	std::cout << "Arg2 = input file path and name" << std::endl;
	std::cout << "Arg3 = output file path and name" << std::endl;
	std::cout << "Arg4 = scanning regime (\"linear\" or \"log\") (optional)" << std::endl;
	std::cout << "Arg5 = debug flag (optional)" << std::endl;
	return 1;
    }
    
    std::string ProcessType = argv[1];
    std::string InputFile = argv[2];
    std::string OutputFile = argv[3];
    std::string SetScale = "";

    if (argc >= 5) {
	SetScale = argv[4];
	if (SetScale != "log" && SetScale != "linear") {
	    std::cout << "Please use \"log\" or \"linear\" as your fourth input." << std::endl;
	    return 1;
	}
    }
    
    bool debug = false;
    if (argc >= 6) {
	std::string argdebug = argv[5];
	if (argdebug == "debug") debug=true;
	else {
	    std::cout << "Did you mean to enable the debug flag? Type \"debug\" on your next try." << std::endl;
	    return -1;
	}
    }
    
    //print out the input arguments
    std::cout << "ProcessType: " << ProcessType << "\tInputFile: " << InputFile << "\tOutputFile: " << OutputFile << "\tDebug: " << debug << std::endl;
    
    //dR truth matching
    dR_match = 1.0;
    
    //////////////////////////////////////////////
    //INPUT
    //////////////////////////////////////////////
    //get input file and tree
    filein = new TFile( InputFile.c_str() );
    treein = (TTree*)filein->Get( "tree" );
    if(debug) treein->Print();
    
    //set up branch linking to addresses
    treein->SetBranchAddress("fspart_id", &finalStateParticles_id);
    treein->SetBranchAddress("fspart_pt", &finalStateParticles_pt);
    treein->SetBranchAddress("fspart_eta",&finalStateParticles_eta);
    treein->SetBranchAddress("fspart_phi",&finalStateParticles_phi);
    treein->SetBranchAddress("fspart_m",  &finalStateParticles_m);
    
    treein->SetBranchAddress("truth_q1_pt",  &truth_q1_pt);
    treein->SetBranchAddress("truth_q1_eta", &truth_q1_eta);
    treein->SetBranchAddress("truth_q1_phi", &truth_q1_phi);
    treein->SetBranchAddress("truth_q1_m",   &truth_q1_m);
    
    treein->SetBranchAddress("truth_q2_pt",  &truth_q2_pt);
    treein->SetBranchAddress("truth_q2_eta", &truth_q2_eta);
    treein->SetBranchAddress("truth_q2_phi", &truth_q2_phi);
    treein->SetBranchAddress("truth_q2_m",   &truth_q2_m);
    
    treein->SetBranchAddress("truth_t1_pt",  &truth_t1_pt);
    treein->SetBranchAddress("truth_t1_eta", &truth_t1_eta);
    treein->SetBranchAddress("truth_t1_phi", &truth_t1_phi);
    treein->SetBranchAddress("truth_t1_m",   &truth_t1_m);
    
    treein->SetBranchAddress("truth_t2_pt",  &truth_t2_pt);
    treein->SetBranchAddress("truth_t2_eta", &truth_t2_eta);
    treein->SetBranchAddress("truth_t2_phi", &truth_t2_phi);
    treein->SetBranchAddress("truth_t2_m",   &truth_t2_m);
    
    treein->SetBranchAddress("truth_W1_pt",  &truth_W1_pt);
    treein->SetBranchAddress("truth_W1_eta", &truth_W1_eta);
    treein->SetBranchAddress("truth_W1_phi", &truth_W1_phi);
    treein->SetBranchAddress("truth_W1_m",   &truth_W1_m);
    
    treein->SetBranchAddress("truth_W2_pt",  &truth_W2_pt);
    treein->SetBranchAddress("truth_W2_eta", &truth_W2_eta);
    treein->SetBranchAddress("truth_W2_phi", &truth_W2_phi);
    treein->SetBranchAddress("truth_W2_m",   &truth_W2_m);
    
    treein->SetBranchAddress("truth_H_pt",  &truth_H_pt);
    treein->SetBranchAddress("truth_H_eta", &truth_H_eta);
    treein->SetBranchAddress("truth_H_phi", &truth_H_phi);
    treein->SetBranchAddress("truth_H_m",   &truth_H_m);
    
    //////////////////////////////////////////////
    //OUTPUT
    //////////////////////////////////////////////
    fileout = new TFile( OutputFile.c_str() ,"RECREATE");
    
    treeout = new TTree("JetTree","JetTree");
    
    treeout->Branch("NumberOfVertices",    &NumberOfVertices);
    
    treeout->Branch("TruthRaw_flavor",        &TruthRaw_flavor);
    treeout->Branch("TruthRaw_pt",            &TruthRaw_pt);
    treeout->Branch("TruthRaw_eta",           &TruthRaw_eta);
    treeout->Branch("TruthRaw_phi",           &TruthRaw_phi);
    treeout->Branch("TruthRaw_m",             &TruthRaw_m);
    treeout->Branch("TruthRaw_Tpruning",      &TruthRaw_Tpruning);
    treeout->Branch("TruthRaw_Ttrimming",     &TruthRaw_Ttrimming);
    
    treeout->Branch("TruthRawTrim_flavor",        &TruthRawTrim_flavor);
    treeout->Branch("TruthRawTrim_pt",            &TruthRawTrim_pt);
    treeout->Branch("TruthRawTrim_eta",           &TruthRawTrim_eta);
    treeout->Branch("TruthRawTrim_phi",           &TruthRawTrim_phi);
    treeout->Branch("TruthRawTrim_m",             &TruthRawTrim_m);
    treeout->Branch("TruthRawTrim_Tau21",         &TruthRawTrim_Tau21);
    treeout->Branch("TruthRawTrim_Tau32",         &TruthRawTrim_Tau32);
    treeout->Branch("TruthRawTrim_D2",            &TruthRawTrim_D2);
    treeout->Branch("TruthRawTrim_T1jet_angle",   &TruthRawTrim_T1jet_angle);
    treeout->Branch("TruthRawTrim_T1jet",         &TruthRawTrim_T1jet);
    treeout->Branch("TruthRawTrim_T2jet_angle",   &TruthRawTrim_T2jet_angle);
    treeout->Branch("TruthRawTrim_T2jet",         &TruthRawTrim_T2jet);
    treeout->Branch("TruthRawTrim_T3jet_minAngle",   &TruthRawTrim_T3jet_angle);
    treeout->Branch("TruthRawTrim_T3jet",         &TruthRawTrim_T3jet);
    treeout->Branch("TruthRawTrim_T3jet_Wmass",       &TruthRawTrim_T3jet_W);
    treeout->Branch("TruthRawTrim_T3jet_WmassVolatility",      &TruthRawTrim_T3jet_mW);
    treeout->Branch("TruthRawTrim_T4jet_angle",   &TruthRawTrim_T4jet_angle);
    treeout->Branch("TruthRawTrim_T4jet",         &TruthRawTrim_T4jet);
    treeout->Branch("TruthRawTrim_T5jet_angle",   &TruthRawTrim_T5jet_angle);
    treeout->Branch("TruthRawTrim_T5jet",         &TruthRawTrim_T5jet);
    
    treeout->Branch("TruthRawTrim_T1jet_pt", &TruthRawTrim_T1jet_pt);
    treeout->Branch("TruthRawTrim_T2jet_pt", &TruthRawTrim_T2jet_pt);
    treeout->Branch("TruthRawTrim_T3jet_pt", &TruthRawTrim_T3jet_pt);
    treeout->Branch("TruthRawTrim_T4jet_pt", &TruthRawTrim_T4jet_pt);
    treeout->Branch("TruthRawTrim_T5jet_pt", &TruthRawTrim_T5jet_pt);
    
    treeout->Branch("TruthRawTrim_T1masses", &TruthRawTrim_T1masses);
    treeout->Branch("TruthRawTrim_T2masses", &TruthRawTrim_T2masses);
    treeout->Branch("TruthRawTrim_T3masses", &TruthRawTrim_T3masses);
    
    treeout->Branch("CaloRaw_flavor",        &CaloRaw_flavor);
    treeout->Branch("CaloRaw_pt",            &CaloRaw_pt);
    treeout->Branch("CaloRaw_eta",           &CaloRaw_eta);
    treeout->Branch("CaloRaw_phi",           &CaloRaw_phi);
    treeout->Branch("CaloRaw_m",             &CaloRaw_m);
    treeout->Branch("CaloRaw_Tpruning",      &CaloRaw_Tpruning);
    treeout->Branch("CaloRaw_Ttrimming",     &CaloRaw_Ttrimming);
    
    treeout->Branch("CaloTrim_flavor",        &CaloTrim_flavor);
    treeout->Branch("CaloTrim_pt",            &CaloTrim_pt);
    treeout->Branch("CaloTrim_eta",           &CaloTrim_eta);
    treeout->Branch("CaloTrim_phi",           &CaloTrim_phi);
    treeout->Branch("CaloTrim_m",             &CaloTrim_m);
    treeout->Branch("CaloTrim_Tau21",         &CaloTrim_Tau21);
    treeout->Branch("CaloTrim_Tau32",         &CaloTrim_Tau32);
    treeout->Branch("CaloTrim_D2",            &CaloTrim_D2);
    treeout->Branch("CaloTrim_T1jet_angle",   &CaloTrim_T1jet_angle);
    treeout->Branch("CaloTrim_T1jet",         &CaloTrim_T1jet);
    treeout->Branch("CaloTrim_T2jet_angle",   &CaloTrim_T2jet_angle);
    treeout->Branch("CaloTrim_T2jet",         &CaloTrim_T2jet);
    treeout->Branch("CaloTrim_T3jet_minAngle",   &CaloTrim_T3jet_angle);
    treeout->Branch("CaloTrim_T3jet",         &CaloTrim_T3jet);
    treeout->Branch("CaloTrim_T3jet_Wmass",       &CaloTrim_T3jet_W);
    treeout->Branch("CaloTrim_T3jet_WmassVolatility",      &CaloTrim_T3jet_mW);
    treeout->Branch("CaloTrim_T4jet_angle",   &CaloTrim_T4jet_angle);
    treeout->Branch("CaloTrim_T4jet",         &CaloTrim_T4jet);
    treeout->Branch("CaloTrim_T5jet_angle",   &CaloTrim_T5jet_angle);
    treeout->Branch("CaloTrim_T5jet",         &CaloTrim_T5jet);
    
    treeout->Branch("CaloTrim_T1jet_pt", &CaloTrim_T1jet_pt);
    treeout->Branch("CaloTrim_T2jet_pt", &CaloTrim_T2jet_pt);
    treeout->Branch("CaloTrim_T3jet_pt", &CaloTrim_T3jet_pt);
    treeout->Branch("CaloTrim_T4jet_pt", &CaloTrim_T4jet_pt);
    treeout->Branch("CaloTrim_T5jet_pt", &CaloTrim_T5jet_pt);
    
    //  Prepare to analyze events.
    nEvents = treein->GetEntries();
    std::cout << "Number of events: " << nEvents << std::endl;

    //  Establish definition of jet.
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.0);
    
    //  Energy correlation functions
    fastjet::contrib::EnergyCorrelatorC2 ecfC2(1);
    fastjet::contrib::EnergyCorrelatorD2 ecfD2(1);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecfC3(2, 1);
    
    //  Create objects for trimming. (arXiv:0912.1342)
    double Rfilt0 = 0.3;
    double fcut0 = 0.05;
    fastjet::Transformer *trimmer = new fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, Rfilt0), fastjet::SelectorPtFractionMin(fcut0) );
    const fastjet::Transformer &f = *trimmer;
    
    //  Loop over events.
    for (Long64_t jentry = 0; jentry < nEvents; jentry++) {
	if (jentry % 1000 == 0) std::cout << jentry << " processed." << std::endl;
	
	event = jentry;
	filein->cd();
	treein->GetEntry(jentry);
	
	ResetBranches();
	
	std::vector<fastjet::PseudoJet> inputParticles;
	
	int numFinalStateParticles = finalStateParticles_id->size();

	//  Convert all final-state particles TLorentzVectors to
	//  PseudoJet four-vectors.
	for (int i = 0; i < numFinalStateParticles; ++i) {
	    
	    if (debug) {
		std::cout << finalStateParticles_id->at(i) << " "
			  << finalStateParticles_pt->at(i) << " "
			  << finalStateParticles_eta->at(i) << " "
			  << finalStateParticles_phi->at(i) << " "
			  << finalStateParticles_m->at(i) << " " << std::endl;
	    }
	    
	    TLorentzVector finalStateParticleMomentum;
	    finalStateParticleMomentum.SetPtEtaPhiM(finalStateParticles_pt->at(i),
				 finalStateParticles_eta->at(i),
				 finalStateParticles_phi->at(i),
				 finalStateParticles_m->at(i));

	    inputParticles.push_back(fastjet::PseudoJet(finalStateParticleMomentum.Px(),
							finalStateParticleMomentum.Py(),
							finalStateParticleMomentum.Pz(),
							finalStateParticleMomentum.E()));
	}
	

	//  Get the resulting jets ordered in pT.
	fastjet::ClusterSequence clust_seq_TruthRaw(inputParticles, jet_def);
	std::vector<fastjet::PseudoJet> jetsTruthRaw = sorted_by_pt(clust_seq_TruthRaw.inclusive_jets(5.0));
	
	if(debug){
	    //  Label the columns.
	    std::cout << "jet#  pt  eta  phi  mass" << std::endl;
	    std::cout << "Inclusive" << std::endl;
	    //  Print out the details for each jet
	    for (unsigned int i = 0; i < jetsTruthRaw.size(); i++) {
		std::cout<< i << "  " << jetsTruthRaw[i].pt()
			 << "  " << jetsTruthRaw[i].eta()
			 << "  " << jetsTruthRaw[i].phi()
			 << "  " << jetsTruthRaw[i].m() << std::endl;
	    }
	}
	

	//  Get truth objects for truth matching.
	truth_q1.SetPtEtaPhiM(truth_q1_pt,truth_q1_eta,truth_q1_phi,truth_q1_m);
	truth_q2.SetPtEtaPhiM(truth_q2_pt,truth_q2_eta,truth_q2_phi,truth_q2_m);
	truth_t1.SetPtEtaPhiM(truth_t1_pt,truth_t1_eta,truth_t1_phi,truth_t1_m);
	truth_t2.SetPtEtaPhiM(truth_t2_pt,truth_t2_eta,truth_t2_phi,truth_t2_m);
	truth_W1.SetPtEtaPhiM(truth_W1_pt,truth_W1_eta,truth_W1_phi,truth_W1_m);
	truth_W2.SetPtEtaPhiM(truth_W2_pt,truth_W2_eta,truth_W2_phi,truth_W2_m);
	truth_H.SetPtEtaPhiM(truth_H_pt,truth_H_eta,truth_H_phi,truth_H_m);
	
	//  NOTE! All TJet errors ensue from fastjet grooming away the
	//  entire jet. This yields a zero mass or even a negative mass
	//  jet. So when you see a WARNING zero entries . . . error, it is
	//  because fastjet has groomed away the entire jet for whatever
	//  reason. This happens now and again but does not affect the
	//  statistics of the ntupled sample heavily.

	//  Set parameters for telescoping subjet algorithm.
	double minR = 0.1;
	double maxR = 1.0;
	int numRadii = 36;
	int stepScale = 0;

	//  Change appropriate settings for specific scaling regime.
	if (SetScale.compare("linear") == 0) {
	    stepScale = 0;
	    numRadii = 36;
	}
	else if (SetScale.compare("log") == 0) {
	    stepScale = 1;
	    numRadii = 32;
	}

	//  Loop over truth jets in event.
	for (unsigned int iJet = 0; iJet < jetsTruthRaw.size(); iJet++) {
	    //  Convert ungroomed truth PseudoJet momentum back to
	    //  TLorentzVector momentum.
	    TLorentzVector tempJet;
	    tempJet.SetPtEtaPhiM(jetsTruthRaw[iJet].pt(),
				 jetsTruthRaw[iJet].eta(),
				 jetsTruthRaw[iJet].phi(),
				 jetsTruthRaw[iJet].m());
	    
	    jetflavor = GetJetTruthFlavor(tempJet, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
	    if (jetflavor == -1) continue;

	    //  Fill output variables.
	    TruthRaw_flavor.push_back(jetflavor);	    
	    TruthRaw_pt.push_back(tempJet.Pt());
	    TruthRaw_eta.push_back(tempJet.Eta());
	    TruthRaw_phi.push_back(tempJet.Phi());
	    TruthRaw_m.push_back(tempJet.M());

	    //  Trimmed truth jets.
	    fastjet::PseudoJet groomedJet = f(jetsTruthRaw[iJet]);

	    //  Convert trimmed truth PseudoJet momentum back to
	    //  TLorentzVector momentum.
	    tempJet.SetPtEtaPhiM(groomedJet.pt(),
				 groomedJet.eta(),
				 groomedJet.phi(),
				 groomedJet.m());
	    
	    jetflavor = GetJetTruthFlavor(tempJet, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
	    if (debug) std::cout << "FillingJet Trimmed: flav=" << jetflavor << "  pt=" <<tempJet.Pt() << "  m=" << tempJet.M() << std::endl;
	    
	    if (jetflavor == -1) continue;
            
            TelescopingJets* telescopeTruthGroomedJet = new TelescopingJets(groomedJet);

	    //  Run telescoping subjet algorithm with trimmed truth jet.
	    tSub T1SubOutputTrim = telescopeTruthGroomedJet->tNSubjet(1, minR, maxR, numRadii, stepScale, 0, 0.0);
	    tSub T2SubOutputTrim = telescopeTruthGroomedJet->tNSubjet(2, minR, maxR, numRadii, stepScale, 0, 0.0);
	    tSub T3SubOutputTrim = telescopeTruthGroomedJet->tNSubjet(3, minR, maxR, numRadii, stepScale, 0, 80.4);
	    tSub T4SubOutputTrim = telescopeTruthGroomedJet->tNSubjet(4, minR, maxR, numRadii, stepScale, 0, 0.0);
	    tSub T5SubOutputTrim = telescopeTruthGroomedJet->tNSubjet(5, minR, maxR, numRadii, stepScale, 0, 0.0);

	    //  Fill output variables.
	    TruthRawTrim_flavor.push_back(jetflavor);
	    
	    TruthRawTrim_pt.push_back(tempJet.Pt());
	    TruthRawTrim_eta.push_back(tempJet.Eta());
	    TruthRawTrim_phi.push_back(tempJet.Phi());
	    TruthRawTrim_m.push_back(tempJet.M());
	    
	    TruthRawTrim_Tau21.push_back(GetTau21(groomedJet));
	    TruthRawTrim_Tau32.push_back(GetTau32(groomedJet));
	    
	    TruthRawTrim_T1jet_angle.push_back(T1SubOutputTrim.minAngle);
	    TruthRawTrim_T1jet.push_back(T1SubOutputTrim.massVolatility);
	    TruthRawTrim_T1masses.push_back(T1SubOutputTrim.masses);
	    TruthRawTrim_T1jet_pt.push_back(T1SubOutputTrim.pTVolatility);
	    
	    TruthRawTrim_T2jet_angle.push_back(T2SubOutputTrim.minAngle);
	    TruthRawTrim_T2jet.push_back(T2SubOutputTrim.massVolatility);
	    TruthRawTrim_T2masses.push_back(T2SubOutputTrim.masses);
	    TruthRawTrim_T2jet_pt.push_back(T2SubOutputTrim.pTVolatility);
	    
	    TruthRawTrim_T3jet_angle.push_back(T3SubOutputTrim.minAngle);
	    TruthRawTrim_T3jet.push_back(T3SubOutputTrim.massVolatility);
	    TruthRawTrim_T3jet_W.push_back(T3SubOutputTrim.targetMass);
	    TruthRawTrim_T3jet_mW.push_back(T3SubOutputTrim.targetMassVolatility);
	    TruthRawTrim_T3masses.push_back(T3SubOutputTrim.masses);
	    TruthRawTrim_T3jet_pt.push_back(T3SubOutputTrim.pTVolatility);
	    
	    TruthRawTrim_T4jet.push_back(T4SubOutputTrim.massVolatility);
	    TruthRawTrim_T4jet_pt.push_back(T4SubOutputTrim.pTVolatility);
	    TruthRawTrim_T5jet.push_back(T5SubOutputTrim.massVolatility);
	    TruthRawTrim_T5jet_pt.push_back(T5SubOutputTrim.pTVolatility);
	    
	}

	//  Put the final-state particles into pseudo-calorimeter
	//  clusters and retrieve the calorimeter jets.
	std::vector<fastjet::PseudoJet> caloClusters = ToyCalorimeter(inputParticles);
	fastjet::ClusterSequence clusterSequenceCalo(caloClusters, jet_def);
	std::vector<fastjet::PseudoJet> caloJets = sorted_by_pt(clusterSequenceCalo.inclusive_jets(5.0));

	//  Loop over calorimeter jets in event.
	for (unsigned int i = 0; i < caloJets.size(); ++i) {
	    //  Convert ungroomed calorimeter jet PseudoJet
	    //  four-vector to TLorentzVector
	    TLorentzVector tempJet;
	    tempJet.SetPtEtaPhiM(caloJets[i].pt(),
				 caloJets[i].eta(),
				 caloJets[i].phi(),
				 caloJets[i].m());

	    jetflavor = GetJetTruthFlavor(tempJet, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
	    if (jetflavor == -1) continue;

	    //  Fill output variables.
	    CaloRaw_flavor.push_back(jetflavor);
	    CaloRaw_pt.push_back(tempJet.Pt());
	    CaloRaw_eta.push_back(tempJet.Eta());
	    CaloRaw_phi.push_back(tempJet.Phi());
	    CaloRaw_m.push_back(tempJet.M());
	    
	    //  Trim the calorimeter jet.
	    fastjet::PseudoJet groomedCaloJet = f(caloJets[i]);
	    //  Convert trimmed calorimeter jet PseudoJet
            //  four-vector to TLorentzVector
	    tempJet.SetPtEtaPhiM(groomedCaloJet.pt(),
				 groomedCaloJet.eta(),
				 groomedCaloJet.phi(),
				 groomedCaloJet.m());
	    
	    jetflavor = GetJetTruthFlavor(tempJet, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
	    if (jetflavor == -1) continue;

            TelescopingJets* telescopeCaloGroomedJet = new TelescopingJets(groomedCaloJet);

	    //  Run telescoping subjet algorithm with trimmed
	    //  calorimeter jet.
	    tSub T1CaloJetTrim = telescopeCaloGroomedJet->tNSubjet(1, minR, maxR, numRadii, stepScale, 0, 0.0);
	    tSub T2CaloJetTrim = telescopeCaloGroomedJet->tNSubjet(2, minR, maxR, numRadii, stepScale, 0, 0.0);
	    tSub T3CaloJetTrim = telescopeCaloGroomedJet->tNSubjet(3, minR, maxR, numRadii, stepScale, 0, 80.4);
	    tSub T4CaloJetTrim = telescopeCaloGroomedJet->tNSubjet(4, minR, maxR, numRadii, stepScale, 0, 0.0);
	    tSub T5CaloJetTrim = telescopeCaloGroomedJet->tNSubjet(5, minR, maxR, numRadii, stepScale, 0, 0.0);

	    //  Fill output variables.
	    CaloTrim_flavor.push_back(jetflavor);
	    
	    CaloTrim_pt.push_back(tempJet.Pt());
	    CaloTrim_eta.push_back(tempJet.Eta());
	    CaloTrim_phi.push_back(tempJet.Phi());
	    CaloTrim_m.push_back(tempJet.M());
	    
	    CaloTrim_Tau21.push_back(GetTau21(groomedCaloJet));
	    CaloTrim_Tau32.push_back(GetTau32(groomedCaloJet));
	    
	    CaloTrim_T1jet_angle.push_back(T1CaloJetTrim.minAngle);
	    CaloTrim_T1jet.push_back(T1CaloJetTrim.massVolatility);
	    CaloTrim_T1jet.push_back(T1CaloJetTrim.pTVolatility);
	    
	    CaloTrim_T2jet_angle.push_back(T2CaloJetTrim.minAngle);
	    CaloTrim_T2jet.push_back(T2CaloJetTrim.massVolatility);
	    CaloTrim_T2jet_pt.push_back(T2CaloJetTrim.pTVolatility);
	    
	    CaloTrim_T3jet_angle.push_back(T3CaloJetTrim.minAngle);
	    CaloTrim_T3jet.push_back(T3CaloJetTrim.massVolatility);
	    CaloTrim_T3jet_W.push_back(T3CaloJetTrim.targetMass);
	    CaloTrim_T3jet_mW.push_back(T3CaloJetTrim.targetMassVolatility);
	    CaloTrim_T3jet_pt.push_back(T3CaloJetTrim.pTVolatility);
	    
	    CaloTrim_T4jet.push_back(T4CaloJetTrim.massVolatility);
	    CaloTrim_T4jet_pt.push_back(T4CaloJetTrim.pTVolatility);
	    CaloTrim_T5jet.push_back(T5CaloJetTrim.massVolatility);
	    CaloTrim_T5jet_pt.push_back(T5CaloJetTrim.pTVolatility);
	}
	
	if(debug) std::cout<<"Filling Tree"<< std::endl;
	treeout->Fill();
    }

    //  Write output to TTree.
    fileout->cd();
    treeout->Write();
    fileout->Close();    
    return 0;    
}

///=========================================
/// Reset Branches
///=========================================
void ResetBranches(){
    NumberOfVertices = 0;
    
    TruthRaw_flavor.clear();
    TruthRaw_pt.clear();
    TruthRaw_eta.clear();
    TruthRaw_phi.clear();
    TruthRaw_m.clear();
    TruthRaw_Ttrimming.clear();
    
    TruthRawTrim_flavor.clear();
    TruthRawTrim_pt.clear();
    TruthRawTrim_eta.clear();
    TruthRawTrim_phi.clear();
    TruthRawTrim_m.clear();
    TruthRawTrim_Tau21.clear();
    TruthRawTrim_Tau32.clear();
    TruthRawTrim_D2.clear();
    TruthRawTrim_T1jet_angle.clear();
    TruthRawTrim_T1jet.clear();
    TruthRawTrim_T2jet_angle.clear();
    TruthRawTrim_T2jet.clear();
    TruthRawTrim_T3jet_angle.clear();
    TruthRawTrim_T3jet_angle1.clear();
    TruthRawTrim_T3jet_angle2.clear();
    TruthRawTrim_T3jet.clear();
    TruthRawTrim_T3jet_W.clear();
    TruthRawTrim_T3jet_mW.clear();
    TruthRawTrim_T4jet_angle.clear();
    TruthRawTrim_T4jet.clear();
    TruthRawTrim_T5jet_angle.clear();
    TruthRawTrim_T5jet.clear();
    
    TruthRawTrim_T1jet_pt.clear();
    TruthRawTrim_T2jet_pt.clear();
    TruthRawTrim_T3jet_pt.clear();
    TruthRawTrim_T4jet_pt.clear();
    TruthRawTrim_T5jet_pt.clear();
    
    TruthRawTrim_T1masses.clear();
    TruthRawTrim_T2masses.clear();
    TruthRawTrim_T3masses.clear();
    
    TruthRawTrim_T2jet_massW.clear();
    TruthRawTrim_T2jet_volatilityW.clear();
    
    CaloRaw_flavor.clear();
    CaloRaw_pt.clear();
    CaloRaw_eta.clear();
    CaloRaw_phi.clear();
    CaloRaw_m.clear();
    
    CaloTrim_flavor.clear();
    CaloTrim_pt.clear();
    CaloTrim_eta.clear();
    CaloTrim_phi.clear();
    CaloTrim_m.clear();
    CaloTrim_Tau21.clear();
    CaloTrim_Tau32.clear();
    CaloTrim_D2.clear();
    CaloTrim_T1jet_angle.clear();
    CaloTrim_T1jet.clear();
    CaloTrim_T2jet_angle.clear();
    CaloTrim_T2jet.clear();
    CaloTrim_T3jet_angle.clear();
    CaloTrim_T3jet_angle1.clear();
    CaloTrim_T3jet_angle2.clear();
    CaloTrim_T3jet.clear();
    CaloTrim_T3jet_W.clear();
    CaloTrim_T3jet_mW.clear();
    CaloTrim_T4jet_angle.clear();
    CaloTrim_T4jet.clear();
    CaloTrim_T5jet_angle.clear();
    CaloTrim_T5jet.clear();
    CaloTrim_Tpruning.clear();
    CaloTrim_Ttrimming.clear();
    CaloTrim_Taktreclustering.clear();
    CaloTrim_Tktreclustering.clear();
    CaloTrim_TJet_m1.clear();
    CaloTrim_TJet_m2.clear();
    
    CaloTrim_T1jet_pt.clear();
    CaloTrim_T2jet_pt.clear();
    CaloTrim_T3jet_pt.clear();
    CaloTrim_T4jet_pt.clear();
    CaloTrim_T5jet_pt.clear();
}

///=========================================
/// Calorimeter Simulation
///=========================================
std::vector<fastjet::PseudoJet> ToyCalorimeter(std::vector<fastjet::PseudoJet> truthParticles) {
    //  Set resolution of the calorimeter.
    const double pi = 3.14159265359;
    const double etaLim = 5.0;
    const int nEta = 100;
    const int nPhi = 63;
    double dEta = 2 * etaLim / nEta;
    double dPhi = 2 * pi / nPhi;

    //  Create calorimeter.
    double tower[nEta][nPhi];
    for (int i = 0; i < nEta; i++) {
	for (int j = 0; j < nPhi; j++)  {
	    tower[i][j] = -0.001;
	}
    }
    

    //  Add energy of final-state particle to the appropriate tower.
    for (int p = 0; p < (int)truthParticles.size(); ++p) {
	fastjet::PseudoJet particle = truthParticles.at(p);
	
	int etaCell = int((particle.eta() +etaLim) / dEta);
	int phiCell = int(particle.phi() / dPhi);
	if (etaCell >= 0 && etaCell < nEta && phiCell >=0 && phiCell < nPhi){
	    tower[etaCell][phiCell] += particle.e();
	}
    }

    //  Smear four-momentum with local eta, phi, and theta and scale
    //  by the total energy in the tower to generate calorimeter particles.
    std::vector<fastjet::PseudoJet> cellParticles;
    for (int i = 0; i < nEta; ++i) {
	for (int j = 0; j < nPhi; ++j) {
	    if (tower[i][j] > 0) {
		double etaLocal = -etaLim + dEta*(i+0.5);
		double phiLocal = dPhi*(j+0.5);
		double thetaLocal = 2*atan(exp(-etaLocal));
		cellParticles.push_back(tower[i][j] * fastjet::PseudoJet(sin(thetaLocal)*cos(phiLocal), 
									  sin(thetaLocal)*sin(phiLocal),
									  cos(thetaLocal),
									  1));
	    }
	}
    }
    
    return cellParticles;
}


///========================================
int GetJetTruthFlavor(TLorentzVector tempJet,
                      TLorentzVector truth_t1,
                      TLorentzVector truth_t2,
                      TLorentzVector truth_W1,
                      TLorentzVector truth_W2,
                      TLorentzVector truth_H,
                      int debug){
  if(debug){
      std::cout<<"DeltaR:   "<< std::endl
	       <<"dRMatch:  "<<dR_match<< std::endl
	       <<"q1:       "<<tempJet.DeltaR(truth_q1)<< std::endl
	       <<"q2:       "<<tempJet.DeltaR(truth_q2)<< std::endl
	       <<"W1:       "<<tempJet.DeltaR(truth_W1)<< std::endl
	       <<"W2:       "<<tempJet.DeltaR(truth_W2)<< std::endl
	       <<"H:        "<<tempJet.DeltaR(truth_H)<< std::endl
	       <<"t1:       "<<tempJet.DeltaR(truth_t1)<< std::endl
	       <<"t2:       "<<tempJet.DeltaR(truth_t2)<< std::endl;
  }
  int jetflavor = -1;
  if     (tempJet.DeltaR(truth_q1)<dR_match || tempJet.DeltaR(truth_q2)<dR_match){
      jetflavor = 0;
  }
  else if(tempJet.DeltaR(truth_W1)<dR_match || tempJet.DeltaR(truth_W2)<dR_match){
      jetflavor = 1;
  }
  else if(tempJet.DeltaR(truth_t1)<dR_match || tempJet.DeltaR(truth_t2)<dR_match){
      jetflavor = 3;
  }
  else if(tempJet.DeltaR(truth_H)<dR_match){
      jetflavor = 3;
  }
  else{
      jetflavor = -1;
  }
  
  if(debug) std::cout<<"Found jet truth flavor: "<<jetflavor<< std::endl;
  
  return jetflavor;
}


double GetTau21(fastjet::PseudoJet& input){
    float tau21=-1;
    
    //N-subjettiness
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(1.);
    fastjet::contrib::Nsubjettiness nsub1(1, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
    fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
    //  fastjet::contrib::Nsubjettiness nsub1(1, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    //  fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    
    float tau1 = nsub1(input);
    float tau2 = nsub2(input);
    
    if(tau1>0)
	tau21 = tau2/tau1;
    
    return tau21;
    
}

double GetTau32(fastjet::PseudoJet& input){
    float tau32=-1;
    
    //N-subjettiness
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(1.);
    fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
    fastjet::contrib::Nsubjettiness nsub3(3, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
    //  fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    //  fastjet::contrib::Nsubjettiness nsub3(3, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    
    float tau2 = nsub2(input);
    float tau3 = nsub3(input);
    
    if(tau2 > 0) tau32 = tau3 / tau2;
    
    return tau32;    
}

