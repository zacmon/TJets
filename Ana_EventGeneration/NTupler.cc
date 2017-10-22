/*----------------------------------------------------------------------

TITLE : NTupler.cc

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

  //exit if you dont pass a run card
  if (argc < 4) {
    std::cout << "You need to specify more arguments" << std::endl;
    std::cout << "Arg1 = process type" << std::endl;
    std::cout << "Arg2 = input file path and name" << std::endl;
    std::cout << "Arg3 = output file path and name" << std::endl;
    std::cout << "Arg4 = debug flag (optional)" << std::endl;
    return 1;
  }

  std::string ProcessType = argv[1];
  std::string InputFile = argv[2];
  std::string OutputFile = argv[3];

  //debug flag
  bool debug = false;
  if (argc >= 5) {
      std::string argdebug = argv[4];
      if (argdebug == "debug") debug=true;
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
  treein->SetBranchAddress("fspart_id", &fspart_id);
  treein->SetBranchAddress("fspart_pt", &fspart_pt);
  treein->SetBranchAddress("fspart_eta",&fspart_eta);
  treein->SetBranchAddress("fspart_phi",&fspart_phi);
  treein->SetBranchAddress("fspart_m",  &fspart_m);

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
  treeout->Branch("TruthRaw_Tau21",         &TruthRaw_Tau21);
  treeout->Branch("TruthRaw_Tau32",         &TruthRaw_Tau32);
  treeout->Branch("TruthRaw_D2",            &TruthRaw_D2);
  treeout->Branch("TruthRaw_T1jet_angle",   &TruthRaw_T1jet_angle);
  treeout->Branch("TruthRaw_T1jet",         &TruthRaw_T1jet);
  treeout->Branch("TruthRaw_T2jet_angle",   &TruthRaw_T2jet_angle);
  treeout->Branch("TruthRaw_T2jet",         &TruthRaw_T2jet);
  treeout->Branch("TruthRaw_T3jet_minAngle",   &TruthRaw_T3jet_angle);
  treeout->Branch("TruthRaw_T3jet",         &TruthRaw_T3jet);
  treeout->Branch("TruthRaw_T3jet_Wmass",       &TruthRaw_T3jet_W);
  treeout->Branch("TruthRaw_T3jet_WmassVolatility",      &TruthRaw_T3jet_mW);
  treeout->Branch("TruthRaw_T4jet_angle",   &TruthRaw_T4jet_angle);
  treeout->Branch("TruthRaw_T4jet",         &TruthRaw_T4jet);
  treeout->Branch("TruthRaw_T5jet_angle",   &TruthRaw_T5jet_angle);
  treeout->Branch("TruthRaw_T5jet",         &TruthRaw_T5jet);
  treeout->Branch("TruthRaw_Tpruning",      &TruthRaw_Tpruning);
  treeout->Branch("TruthRaw_Ttrimming",     &TruthRaw_Ttrimming);
  treeout->Branch("TruthRaw_Taktreclustering",	&TruthRaw_Taktreclustering);
  treeout->Branch("TruthRaw_Tktreclustering",	&TruthRaw_Tktreclustering);
  treeout->Branch("TruthRaw_TJet_m1",       &TruthRaw_TJet_m1);
  treeout->Branch("TruthRaw_TJet_m2",       &TruthRaw_TJet_m2);

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
  treeout->Branch("TruthRawTrim_Tpruning",      &TruthRawTrim_Tpruning);
  treeout->Branch("TruthRawTrim_Ttrimming",     &TruthRawTrim_Ttrimming);
  treeout->Branch("TruthRawTrim_Taktreclustering",	&TruthRawTrim_Taktreclustering);
  treeout->Branch("TruthRawTrim_Tktreclustering",	&TruthRawTrim_Tktreclustering);
  treeout->Branch("TruthRawTrim_TJet_m1",       &TruthRawTrim_TJet_m1);
  treeout->Branch("TruthRawTrim_TJet_m2",       &TruthRawTrim_TJet_m2);

  treeout->Branch("TruthRawTrim_T1masses", &TruthRawTrim_T1masses);
  treeout->Branch("TruthRawTrim_T2masses", &TruthRawTrim_T2masses);
  treeout->Branch("TruthRawTrim_T3masses", &TruthRawTrim_T3masses);

  treeout->Branch("CaloRaw_flavor",        &CaloRaw_flavor);
  treeout->Branch("CaloRaw_pt",            &CaloRaw_pt);
  treeout->Branch("CaloRaw_eta",           &CaloRaw_eta);
  treeout->Branch("CaloRaw_phi",           &CaloRaw_phi);
  treeout->Branch("CaloRaw_m",             &CaloRaw_m);
  treeout->Branch("CaloRaw_Tau21",         &CaloRaw_Tau21);
  treeout->Branch("CaloRaw_Tau32",         &CaloRaw_Tau32);
  treeout->Branch("CaloRaw_D2",            &CaloRaw_D2);
  treeout->Branch("CaloRaw_T1jet_angle",   &CaloRaw_T1jet_angle);
  treeout->Branch("CaloRaw_T1jet",         &CaloRaw_T1jet);
  treeout->Branch("CaloRaw_T2jet_angle",   &CaloRaw_T2jet_angle);
  treeout->Branch("CaloRaw_T2jet",         &CaloRaw_T2jet);
  treeout->Branch("CaloRaw_T3jet_minAngle",   &CaloRaw_T3jet_angle);
  treeout->Branch("CaloRaw_T3jet",         &CaloRaw_T3jet);
  treeout->Branch("CaloRaw_T3jet_Wmass",       &CaloRaw_T3jet_W);
  treeout->Branch("CaloRaw_T3jet_WmassVolatility",      &CaloRaw_T3jet_mW);
  treeout->Branch("CaloRaw_T4jet_angle",   &CaloRaw_T4jet_angle);
  treeout->Branch("CaloRaw_T4jet",         &CaloRaw_T4jet);
  treeout->Branch("CaloRaw_T5jet_angle",   &CaloRaw_T5jet_angle);
  treeout->Branch("CaloRaw_T5jet",         &CaloRaw_T5jet);
  treeout->Branch("CaloRaw_Tpruning",      &CaloRaw_Tpruning);
  treeout->Branch("CaloRaw_Ttrimming",     &CaloRaw_Ttrimming);
  treeout->Branch("CaloRaw_Taktreclustering",	&CaloRaw_Taktreclustering);
  treeout->Branch("CaloRaw_Tktreclustering",	&CaloRaw_Tktreclustering);
  treeout->Branch("CaloRaw_TJet_m1",       &CaloRaw_TJet_m1);
  treeout->Branch("CaloRaw_TJet_m2",       &CaloRaw_TJet_m2);
  
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
  treeout->Branch("CaloTrim_Tpruning",      &CaloTrim_Tpruning);
  treeout->Branch("CaloTrim_Ttrimming",     &CaloTrim_Ttrimming);
  treeout->Branch("CaloTrim_Taktreclustering",	&CaloTrim_Taktreclustering);
  treeout->Branch("CaloTrim_Tktreclustering",	&CaloTrim_Tktreclustering);
  treeout->Branch("CaloTrim_TJet_m1",       &CaloTrim_TJet_m1);
  treeout->Branch("CaloTrim_TJet_m2",       &CaloTrim_TJet_m2);
  
  //////////////////////////////////////////////
  //random number generator for pileup
  //////////////////////////////////////////////
  TRandom3 *rand_pileup = new TRandom3();

  //////////////////////////////////////////////
  //main event loop
  //////////////////////////////////////////////
  nEvents = treein->GetEntries();
  
  std::cout << "Number of events: " << nEvents << std::endl;

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.0);
  
  //  Energy correlation functions
  fastjet::contrib::EnergyCorrelatorC2 ecfC2(1);
  fastjet::contrib::EnergyCorrelatorD2 ecfD2(1);
  fastjet::contrib::EnergyCorrelatorDoubleRatio ecfC3(2, 1);
  
  //  Filtering with a pt cut as for trimming (arXiv:0912.1342)
  double Rfilt0 = 0.3;
  double fcut0 = 0.05;
  fastjet::Transformer *trimmer = new fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, Rfilt0), fastjet::SelectorPtFractionMin(fcut0) );
  const fastjet::Transformer &f = *trimmer;
  
  for (Long64_t jentry=0; jentry<nEvents; jentry++) {
      
      if (jentry % 1000 == 0) std::cout << jentry << " processed." << std::endl;

      filein->cd();
      treein->GetEntry(jentry);
      
      ResetBranches();

      std::vector<fastjet::PseudoJet> input_particles;
      
      int n_fspart = fspart_id->size();
      for (int i_fspart=0; i_fspart<n_fspart; i_fspart++) {
	  
	  if (debug) {
	      std::cout << fspart_id->at(i_fspart) << " "
			<< fspart_pt->at(i_fspart) << " "
			<< fspart_eta->at(i_fspart) << " "
			<< fspart_phi->at(i_fspart) << " "
			<< fspart_m->at(i_fspart) << " " << std::endl;
	  }
	  
	  TLorentzVector temp_p4;
	  temp_p4.SetPtEtaPhiM(fspart_pt->at(i_fspart),
			       fspart_eta->at(i_fspart),
			       fspart_phi->at(i_fspart),
			       fspart_m->at(i_fspart));
	  	  
	  input_particles.push_back(fastjet::PseudoJet(temp_p4.Px(),temp_p4.Py(),temp_p4.Pz(),temp_p4.E()));
      }

      //////////////////////////////////////////////
      // get the resulting jets ordered in pt
      //////////////////////////////////////////////
      
      fastjet::ClusterSequence clust_seq_TruthRaw(input_particles, jet_def);
      std::vector<fastjet::PseudoJet> jetsTruthRaw = sorted_by_pt(clust_seq_TruthRaw.inclusive_jets(5.0));
      
      if(debug){
	  // label the columns
	  std::cout<<"jet#  pt  eta  phi  mass"<< std::endl;
	  std::cout<<"Inclusive"<< std::endl;
	  // print out the details for each jet
	  for (unsigned int i = 0; i < jetsTruthRaw.size(); i++) {
	      std::cout<<i<<"  "<<jetsTruthRaw[i].pt()
		       <<"  "<<jetsTruthRaw[i].eta()
		       <<"  "<<jetsTruthRaw[i].phi()
		       <<"  "<<jetsTruthRaw[i].m()<< std::endl;
	  }
      }
      
    /////////////////////////////////////////////
    //Get truth objects for truth matching
    /////////////////////////////////////////////
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

    double minR = 0.1;
    double maxR = 1.0;
    int numRadii = 36;
    
    for (unsigned int ijet = 0; ijet < jetsTruthRaw.size(); ijet++) {

	/////////////////////////////
	//TruthRaw
	/////////////////////////////
	TLorentzVector tempJet;
	tempJet.SetPtEtaPhiM(jetsTruthRaw[ijet].pt(),
			     jetsTruthRaw[ijet].eta(),
			     jetsTruthRaw[ijet].phi(),
			     jetsTruthRaw[ijet].m());

	jetflavor = GetJetTruthFlavor(tempJet, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
	if (jetflavor == -1) continue;

	// TSub T1TruthRaw = TNSubjet(jetsTruthRaw[ijet], 1, minR, maxR, numRadii);
	// T2Sub T2TruthRaw = T2Subjet(jetsTruthRaw[ijet], minR, maxR, numRadii);
	// T3Sub T3TruthRaw = T3Subjet(jetsTruthRaw[ijet], minR, maxR, numRadii);

	TruthRaw_flavor.push_back(jetflavor);
	
	TruthRaw_pt.push_back(tempJet.Pt());
	TruthRaw_eta.push_back(tempJet.Eta());
	TruthRaw_phi.push_back(tempJet.Phi());
	TruthRaw_m.push_back(tempJet.M());
	
	// TruthRaw_Tau21.push_back(GetTau21(jetsTruthRaw[ijet]));
	// TruthRaw_Tau32.push_back(GetTau32(jetsTruthRaw[ijet]));
	
	// TruthRaw_T1jet_angle.push_back(T1TruthRaw.minAngle);
	// TruthRaw_T1jet.push_back(T1TruthRaw.volatility);
	
	// TruthRaw_T2jet_angle.push_back(T2TruthRaw.minAngle);
	// TruthRaw_T2jet.push_back(T2TruthRaw.volatility);
	
	// TruthRaw_T3jet_angle.push_back(T3TruthRaw.minAngle);
	// TruthRaw_T3jet.push_back(T3TruthRaw.volatility);
	// TruthRaw_T3jet_W.push_back(T3TruthRaw.massW);
	// TruthRaw_T3jet_mW.push_back(T3TruthRaw.volatilityW);

	
	double TruthRawTPruning = T_Pruning(jetsTruthRaw[ijet], 0.1, 2.0, 20);
	TruthRaw_Tpruning.push_back(TruthRawTPruning);

	double TruthRawTTrimming =	T_Trimming(jetsTruthRaw[ijet], 0.0, 0.1, 20);
	TruthRaw_Ttrimming.push_back(TruthRawTTrimming);
	
	/////////////////////////////
	//TruthRawTrim
	/////////////////////////////
	fastjet::PseudoJet groomedJet = f(jetsTruthRaw[ijet]);
	
	tempJet.SetPtEtaPhiM(groomedJet.pt(),
			     groomedJet.eta(),
			     groomedJet.phi(),
			     groomedJet.m());
	
	/////////////////////////////////
	//Getting truth label for filling into ntuple
	/////////////////////////////////
	jetflavor = GetJetTruthFlavor(tempJet, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
	if (debug) std::cout<<"FillingJet Trimmed: flav="<<jetflavor<<"  pt="<<tempJet.Pt()<<"  m="<<tempJet.M()<< std::endl;
	
	if (jetflavor == -1) continue;
	
	/////////////////////////////////
	//Fill variables that will go into ntuple
	/////////////////////////////////
	TSub  T1SubOutputTrim  = TNSubjet(groomedJet, 1, minR, maxR, numRadii);
	T2Sub  T2SubOutputTrim  = T2Subjet(groomedJet, minR, maxR, numRadii);
	T3Sub  T3SubOutputTrim = T3Subjet(groomedJet, minR, maxR, numRadii);
	
	TruthRawTrim_flavor.push_back(jetflavor);
	
	TruthRawTrim_pt.push_back(tempJet.Pt());
	TruthRawTrim_eta.push_back(tempJet.Eta());
	TruthRawTrim_phi.push_back(tempJet.Phi());
	TruthRawTrim_m.push_back(tempJet.M());
	
	TruthRawTrim_Tau21.push_back(GetTau21(groomedJet));
	TruthRawTrim_Tau32.push_back(GetTau32(groomedJet));
	
	TruthRawTrim_T1jet_angle.push_back(T1SubOutputTrim.minAngle);
	TruthRawTrim_T1jet.push_back(T1SubOutputTrim.volatility);
	TruthRawTrim_T1masses.push_back(T1SubOutputTrim.masses);
	
	TruthRawTrim_T2jet_angle.push_back(T2SubOutputTrim.minAngle);
	TruthRawTrim_T2jet.push_back(T2SubOutputTrim.volatility);
	TruthRawTrim_T2masses.push_back(T2SubOutputTrim.masses);
	
	TruthRawTrim_T3jet_angle.push_back(T3SubOutputTrim.minAngle);
	TruthRawTrim_T3jet.push_back(T3SubOutputTrim.volatility);
	TruthRawTrim_T3jet_W.push_back(T3SubOutputTrim.massW);
	TruthRawTrim_T3jet_mW.push_back(T3SubOutputTrim.volatilityW);
	TruthRawTrim_T3masses.push_back(T3SubOutputTrim.masses);

	double TruthRawTrimTPruning = T_Pruning(groomedJet, 0.1, 2.0, 20);
	TruthRawTrim_Tpruning.push_back(TruthRawTrimTPruning);

	double TruthRawTrimTTrimming = T_Trimming(groomedJet, 0.0, 0.1, 20);
	TruthRawTrim_Ttrimming.push_back(TruthRawTrimTTrimming);
    }

    std::vector<fastjet::PseudoJet> caloClusters = ToyCalorimeter(input_particles);
    fastjet::ClusterSequence clusterSequenceCalo(caloClusters, jet_def);
    std::vector<fastjet::PseudoJet> caloJets = sorted_by_pt(clusterSequenceCalo.inclusive_jets(5.0));
    
    for (unsigned int i = 0; i < caloJets.size(); ++i) {
	//  Ungroomed calo jet.
	TLorentzVector tempJet;
	tempJet.SetPtEtaPhiM(caloJets[i].pt(),
			     caloJets[i].eta(),
			     caloJets[i].phi(),
			     caloJets[i].m());
	
	jetflavor = GetJetTruthFlavor(tempJet, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
	if (jetflavor == -1) continue;

	// TSub T1CaloJetRaw = TNSubjet(caloJets[i], 1, minR, maxR, numRadii);
	// TSub T2CaloJetRaw = TNSubjet(caloJets[i], 2, minR, maxR, numRadii);
	// T3Sub T3CaloJetRaw = T3Subjet(caloJets[i], minR, maxR, numRadii);

	CaloRaw_flavor.push_back(jetflavor);

	CaloRaw_pt.push_back(tempJet.Pt());
	CaloRaw_eta.push_back(tempJet.Eta());
	CaloRaw_phi.push_back(tempJet.Phi());
	CaloRaw_m.push_back(tempJet.M());
	
	// CaloRaw_Tau21.push_back(GetTau21(caloJets[i]));
	// CaloRaw_Tau32.push_back(GetTau32(caloJets[i]));
	
	// CaloRaw_T1jet_angle.push_back(T1CaloJetRaw.minAngle);
	// CaloRaw_T1jet.push_back(T1CaloJetRaw.volatility);
	
	// CaloRaw_T2jet_angle.push_back(T2CaloJetRaw.minAngle);
	// CaloRaw_T2jet.push_back(T2CaloJetRaw.volatility);
	
	// CaloRaw_T3jet_angle.push_back(T3CaloJetRaw.minAngle);
	// CaloRaw_T3jet.push_back(T3CaloJetRaw.volatility);
	// CaloRaw_T3jet_W.push_back(T3CaloJetRaw.massW);
	// CaloRaw_T3jet_mW.push_back(T3CaloJetRaw.volatilityW);
	
	double CaloRawTPruning = T_Pruning(caloJets[i], 0.1, 2.0, 20);
	CaloRaw_Tpruning.push_back(CaloRawTPruning);

	double CaloRawTTrimming =	T_Trimming(caloJets[i], 0.0, 0.1, 20);
	CaloRaw_Ttrimming.push_back(CaloRawTTrimming);
	
	//  Trimmed calo jet.
	fastjet::PseudoJet groomedCaloJet = f(caloJets[i]);
	
	tempJet.SetPtEtaPhiM(groomedCaloJet.pt(),
			     groomedCaloJet.eta(),
			     groomedCaloJet.phi(),
			     groomedCaloJet.m());
	
	jetflavor = GetJetTruthFlavor(tempJet, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
	if (jetflavor == -1) continue;

	TSub T1CaloJetTrim = TNSubjet(groomedCaloJet, 1, minR, maxR, numRadii);
	TSub T2CaloJetTrim = TNSubjet(groomedCaloJet, 2, minR, maxR, numRadii);
	T3Sub T3CaloJetTrim = T3Subjet(groomedCaloJet, minR, maxR, numRadii);

	CaloTrim_flavor.push_back(jetflavor);

	CaloTrim_pt.push_back(tempJet.Pt());
	CaloTrim_eta.push_back(tempJet.Eta());
	CaloTrim_phi.push_back(tempJet.Phi());
	CaloTrim_m.push_back(tempJet.M());
	
	CaloTrim_Tau21.push_back(GetTau21(groomedCaloJet));
	CaloTrim_Tau32.push_back(GetTau32(groomedCaloJet));
	
	CaloTrim_T1jet_angle.push_back(T1CaloJetTrim.minAngle);
	CaloTrim_T1jet.push_back(T1CaloJetTrim.volatility);
	
	CaloTrim_T2jet_angle.push_back(T2CaloJetTrim.minAngle);
	CaloTrim_T2jet.push_back(T2CaloJetTrim.volatility);
	
	CaloTrim_T3jet_angle.push_back(T3CaloJetTrim.minAngle);
	CaloTrim_T3jet.push_back(T3CaloJetTrim.volatility);
	CaloTrim_T3jet_W.push_back(T3CaloJetTrim.massW);
	CaloTrim_T3jet_mW.push_back(T3CaloJetTrim.volatilityW);

	double CaloTrimTPruning = T_Pruning(groomedCaloJet, 0.1, 2.0, 20);
	CaloTrim_Tpruning.push_back(CaloTrimTPruning);
	
	double CaloTrimTTrimming =	T_Trimming(groomedCaloJet, 0.0, 0.1, 20);
	CaloTrim_Ttrimming.push_back(CaloTrimTTrimming);
    }
    
    if(debug) std::cout<<"Filling Tree"<< std::endl;
    treeout->Fill();
  }
  
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
  TruthRaw_Tau21.clear();
  TruthRaw_Tau32.clear();
  TruthRaw_D2.clear();
  TruthRaw_T1jet_angle.clear();
  TruthRaw_T1jet.clear();
  TruthRaw_T2jet_angle.clear();
  TruthRaw_T2jet.clear();
  TruthRaw_T3jet_angle.clear();
  TruthRaw_T3jet_angle1.clear();
  TruthRaw_T3jet_angle2.clear();
  TruthRaw_T3jet.clear();
  TruthRaw_T3jet_W.clear();
  TruthRaw_T3jet_mW.clear();
  TruthRaw_T4jet_angle.clear();
  TruthRaw_T4jet.clear();
  TruthRaw_T5jet_angle.clear();
  TruthRaw_T5jet.clear();
  TruthRaw_Tpruning.clear();
  TruthRaw_Ttrimming.clear();
  TruthRaw_Taktreclustering.clear();
  TruthRaw_Tktreclustering.clear();
  TruthRaw_TJet_m1.clear();
  TruthRaw_TJet_m2.clear();

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
  TruthRawTrim_Tpruning.clear();
  TruthRawTrim_Ttrimming.clear();
  TruthRawTrim_Taktreclustering.clear();
  TruthRawTrim_Tktreclustering.clear();
  TruthRawTrim_TJet_m1.clear();
  TruthRawTrim_TJet_m2.clear();

  TruthRawTrim_T1Volatility.clear();
  TruthRawTrim_T2Volatility.clear();
  TruthRawTrim_T3Volatility.clear();

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
  CaloRaw_Tau21.clear();
  CaloRaw_Tau32.clear();
  CaloRaw_D2.clear();
  CaloRaw_T1jet_angle.clear();
  CaloRaw_T1jet.clear();
  CaloRaw_T2jet_angle.clear();
  CaloRaw_T2jet.clear();
  CaloRaw_T3jet_angle.clear();
  CaloRaw_T3jet_angle1.clear();
  CaloRaw_T3jet_angle2.clear();
  CaloRaw_T3jet.clear();
  CaloRaw_T3jet_W.clear();
  CaloRaw_T3jet_mW.clear();
  CaloRaw_T4jet_angle.clear();
  CaloRaw_T4jet.clear();
  CaloRaw_T5jet_angle.clear();
  CaloRaw_T5jet.clear();
  CaloRaw_Tpruning.clear();
  CaloRaw_Ttrimming.clear();
  CaloRaw_Taktreclustering.clear();
  CaloRaw_Tktreclustering.clear();
  CaloRaw_TJet_m1.clear();
  CaloRaw_TJet_m2.clear();
  
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
}

///=========================================
/// Calorimeter Simulation
///=========================================
std::vector<fastjet::PseudoJet> ToyCalorimeter(std::vector<fastjet::PseudoJet> truth_particles) {
    const double pi = 3.14159265359;
    const double etaLim = 5.0;
    const int nEta = 100;
    const int nPhi = 63;
    double dEta = 2*etaLim/nEta;
    double dPhi = 2*pi/nPhi;
    
    double tower[nEta][nPhi];
    for (int i = 0; i < nEta; i++) {
	for (int j = 0; j < nPhi; j++)  {
	    tower[i][j] = -0.001;
	}
    }
    
    std::vector<fastjet::PseudoJet> cell_particles;
    for (int p = 0; p < (int)truth_particles.size(); p++) {
	fastjet::PseudoJet part = truth_particles.at(p);
	
	int etaCell = int((part.eta()+etaLim)/dEta);
	int phiCell = int(part.phi()/dPhi);
	if (etaCell >= 0 && etaCell < nEta && phiCell >=0 && phiCell < nPhi){
	    tower[etaCell][phiCell] += part.e();
	}
    }
    
    for (int i = 0; i < nEta; i++) {
	for (int j = 0; j < nPhi; j++) {
	    if (tower[i][j] > 0) {
		double etaLocal = -etaLim + dEta*(i+0.5);
		double phiLocal = dPhi*(j+0.5);
		double thetaLocal = 2*atan(exp(-etaLocal));
		cell_particles.push_back(tower[i][j] * fastjet::PseudoJet(sin(thetaLocal)*cos(phiLocal), 
									  sin(thetaLocal)*sin(phiLocal),
									  cos(thetaLocal),
									  1));
	    }
	}
    }
    return cell_particles;
}
    
///=========================================
/// Telescoping Pruning
///=========================================
double T_Pruning(fastjet::PseudoJet& input, double minDCut, double maxDCut, int numDCuts) {
    double zCut = 0.1; // Single choice of zcut but can be further telescoped.
    std::vector<double> telescopingMasses;

    double dCutStep = (maxDCut - minDCut) / (numDCuts);

    for (double dCut = minDCut; dCut < maxDCut + dCutStep; dCut += dCutStep) {
        fastjet::Pruner pruner(fastjet::cambridge_algorithm, zCut, dCut);
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
double T_Trimming(fastjet::PseudoJet& input, double minFCut, double maxFCut, int numFCuts) {
    double Rfilt = 0.2; // single choice of Rfilt. can be further telescoped.
    //  Use Rfilt = 0.1 for higher pT jets, Rfilt = 0.2 for lower pT jets.
    if (input.pt() > 500) Rflit = 0.1;
    
    std::vector<double> telescopingMasses;
    
    double fCutStep = (maxFCut - minFCut) / (numFCuts);
    
    for (double fCut = minFCut; fCut < maxFCut + fCutStep; fCut += fCutStep) {
        fastjet::Filter trimmer(Rfilt, fastjet::SelectorPtFractionMin(fCut));
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
/// for antiki, set algorithm = 2
///=========================================
double T_Reclustering(fastjet::PseudoJet& input, int algorithm, double minRadius, double maxRadius, int numRadii) {
    std::vector<double> telescopingMass; 
    double deltaR = (maxRadius - minRadius) / (numRadii);

    for (double r = minRadius; r < maxRadius + deltaR; r += deltaR) {
        fastjet::JetDefinition TjetDefinition(fastjet::JetAlgorithm(algorithm), r);
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
TSub TNSubjet(fastjet::PseudoJet& input, unsigned int numSubjets, double minRadius, double maxRadius, int numRadii) {
    TSub result;
    
    double beta = 1.0;    
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSubjettiness(numSubjets, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    //fastjet::contrib::Nsubjettiness nSubjettiness(2, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);

    double tauN = nSubjettiness.result(input);
    std::vector<fastjet::PseudoJet> tauAxes = nSubjettiness.currentAxes();
    std::vector<TLorentzVector> pTauAxes(numSubjets);

    for (unsigned int i = 0; i < numSubjets; ++i) {
	pTauAxes[i].SetPxPyPzE(tauAxes[i].px(), tauAxes[i].py(), tauAxes[i].pz(), tauAxes[i].e());
    }

    std::vector<double> anglesBetweenTauAxes;
    for (unsigned int i = 0; i < numSubjets; ++i) {
	for (unsigned int j = i + 1; j < numSubjets; ++j) {
	    anglesBetweenTauAxes.push_back(pTauAxes[i].DeltaR(pTauAxes[j]));
	}
    }

    if (!anglesBetweenTauAxes.empty()) result.minAngle = *(std::min_element(anglesBetweenTauAxes.cbegin(), anglesBetweenTauAxes.cend()));
    
    std::vector<fastjet::PseudoJet> constituents = input.constituents();
    std::vector<std::vector<std::pair<TLorentzVector, double>>> sortedConstituents(numSubjets);
    
    for (auto const &constituent: constituents) {
	TLorentzVector pConstituent(constituent.px(), constituent.py(), constituent.pz(), constituent.e());

	std::vector<double> distanceToTauAxes(numSubjets);
	for (unsigned int i = 0; i < numSubjets; ++i) {
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
    
    std::vector<TLorentzVector> TSubjets(numSubjets);
    std::vector<double> telescopingMasses;

    double deltaR = (maxRadius - minRadius) / (numRadii);
    
    for (double r = minRadius; r < maxRadius + deltaR; r += deltaR) {
	for (unsigned int i = 0; i < sortedConstituents.size(); ++i) {
	    for (auto it = sortedConstituents[i].begin(); it != sortedConstituents[i].end(); ++it) {
		if (it->second <= r) {
		    TSubjets[i] += it->first;
		    sortedConstituents[i].erase(it);
		    --it;
		}
		else break;
	    }
	}
	
	TLorentzVector sumTSubjet;
	for (auto const &TSubjet : TSubjets) {
	    sumTSubjet += TSubjet;
	}
	
	if (sumTSubjet.M() > M0) {
	    telescopingMasses.push_back(sumTSubjet.M());
	    result.volVec.push_back(getVolatility(telescopingMasses));
	}
    }

    if (!telescopingMasses.empty()) {
	result.volatility = getVolatility(telescopingMasses);
	result.masses = telescopingMasses;
    }
    else {
	std::cout << "WARNING zero entries for TNSubjet!   numSubjets: " << numSubjets <<
	    "   minRadius: " << minRadius <<
	    "   maxRadius: " << maxRadius <<
	    "   numRadii: " << numRadii <<
	    "   input mass: " << input.m() << std::endl;
    }
    return result;
}


T2Sub T2Subjet(fastjet::PseudoJet& input, double minRadius, double maxRadius, int numRadii) {
    T2Sub result;
    
    double beta = 1.0;
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSubjettiness(2, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    //fastjet::contrib::Nsubjettiness nSubjettiness(2, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);

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
    std::vector<double> masses(3);

    double deltaR = (maxRadius - minRadius) / (numRadii);
    for (double r = minRadius; r < maxRadius + deltaR; r += deltaR) {
	for (unsigned int i = 0; i < sortedConstituents.size(); ++i) {
	    for (auto it = sortedConstituents[i].begin(); it != sortedConstituents[i].end(); ++it) {
		if (it->second <= r) {
		    TSubjets[i] += it->first;
		    sortedConstituents[i].erase(it);
		    --it;
		}
		else break;
	    }
	}
    
	masses[0] = (TSubjets[0] + TSubjets[1]).M();
	masses[1] = TSubjets[0].M();
	masses[2] = TSubjets[1].M();

	for (unsigned int i = 0; i < masses.size(); ++i) {
	    telescopingMasses[i].push_back(masses[i]);
	    if (i == 0) result.volVec.push_back(getVolatility(telescopingMasses[i]));
	}
    }
    
    if (!telescopingMasses[0].empty() && !telescopingMasses[1].empty() && !telescopingMasses[2].empty()) {
	result.volatility = getVolatility(telescopingMasses[0]);
	result.masses = telescopingMasses[0];

	std::vector<double> residualWMass = masses;

	std::transform(residualWMass.begin() + 1, residualWMass.end(), residualWMass.begin() + 1,
		       bind2nd(std::minus<double>(), massW));
	std::transform(residualWMass.begin() + 1, residualWMass.end(), residualWMass.begin() + 1,
		       static_cast<double (*)(double)>(&std::abs));

	auto wMassPredictionIt = std::min_element(residualWMass.cbegin() + 1, residualWMass.cend());
	unsigned int wMassPredictionIndex = std::distance(residualWMass.cbegin(), wMassPredictionIt);

	result.massW = masses[wMassPredictionIndex];
	result.volatilityW = getVolatility(telescopingMasses[wMassPredictionIndex]);
    }

    else {
	std::cout << "WARNING zero entries for T2Subjet!   minRadius: " << minRadius <<
	    "   maxRadius: " << maxRadius <<
	    "   numRadii: " << numRadii <<
	    "   input mass: " << input.m() << std::endl;
    }
    return result;
}

T3Sub T3Subjet(fastjet::PseudoJet& input, double minRadius, double maxRadius, int numRadii) {
    T3Sub result;

    double beta = 1.0;
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSubjettiness(3, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    //fastjet::contrib::Nsubjettiness nSubjettiness(2, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);
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
    std::vector<double> masses(4);

    double deltaR = (maxRadius - minRadius) / (numRadii);
    for (double r = minRadius; r < maxRadius + deltaR; r += deltaR) {
	for (unsigned int i = 0; i < sortedConstituents.size(); ++i) {
	    for (auto it = sortedConstituents[i].begin(); it != sortedConstituents[i].end(); ++it) {
		if (it->second <= r) {
		    TSubjets[i] += it->first;
		    sortedConstituents[i].erase(it);
		    --it;
		}
		else break;
	    }
	}
	
	masses[0] = (TSubjets[0] + TSubjets[1] + TSubjets[2]).M();
	masses[1] = (TSubjets[0] + TSubjets[1]).M();
	masses[2] = (TSubjets[0] + TSubjets[2]).M();
	masses[3] = (TSubjets[1] + TSubjets[2]).M();

	for (unsigned int i = 0; i < masses.size(); ++i) {
	    telescopingMasses[i].push_back(masses[i]);
	    if (i == 0) result.volVec.push_back(getVolatility(telescopingMasses[i]));
	}
    }
    
    if (!telescopingMasses[0].empty() && !telescopingMasses[1].empty() && !telescopingMasses[2].empty() && !telescopingMasses[3].empty()) {
        result.volatility = getVolatility(telescopingMasses[0]);
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

	result.massW = masses[wMassPredictionIndex];
	result.volatilityW = getVolatility(telescopingMasses[wMassPredictionIndex]);
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
double T_Nsubjettiness(int N, fastjet::PseudoJet& input, double minBeta, double maxBeta, int numBetas) {
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

