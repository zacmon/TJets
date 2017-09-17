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

//----------------------------------------------------------------------*/


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
  treeout->Branch("TruthRaw_T3jet_angle",   &TruthRaw_T3jet_angle);
  treeout->Branch("TruthRaw_T3jet_angle1",  &TruthRaw_T3jet_angle1);
  treeout->Branch("TruthRaw_T3jet_angle2",  &TruthRaw_T3jet_angle2);
  treeout->Branch("TruthRaw_T3jet",         &TruthRaw_T3jet);
  treeout->Branch("TruthRaw_T3jet_W",       &TruthRaw_T3jet_W);
  treeout->Branch("TruthRaw_T3jet_mW",      &TruthRaw_T3jet_mW);
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
  treeout->Branch("TruthRawTrim_T3jet_angle",   &TruthRawTrim_T3jet_angle);
  treeout->Branch("TruthRawTrim_T3jet_angle1",  &TruthRawTrim_T3jet_angle1);
  treeout->Branch("TruthRawTrim_T3jet_angle2",  &TruthRawTrim_T3jet_angle2);
  treeout->Branch("TruthRawTrim_T3jet",         &TruthRawTrim_T3jet);
  treeout->Branch("TruthRawTrim_T3jet_W",       &TruthRawTrim_T3jet_W);
  treeout->Branch("TruthRawTrim_T3jet_mW",      &TruthRawTrim_T3jet_mW);
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

  treeout->Branch("TruthRawTrim_T1Volatility", &TruthRawTrim_T1Volatility);
  treeout->Branch("TruthRawTrim_T2Volatility", &TruthRawTrim_T2Volatility);
  treeout->Branch("TruthRawTrim_T3Volatility", &TruthRawTrim_T3Volatility);

  treeout->Branch("TruthRawTrim_v32", &TruthRawTrim_v32);

  treeout->Branch("TruthRawTrim_T1masses", &TruthRawTrim_T1masses);
  treeout->Branch("TruthRawTrim_T2masses", &TruthRawTrim_T2masses);
  treeout->Branch("TruthRawTrim_T3masses", &TruthRawTrim_T3masses);
  
  //////////////////////////////////////////////
  //random number generator for pileup
  //////////////////////////////////////////////
  TRandom3 *rand_pileup = new TRandom3();

  //////////////////////////////////////////////
  //main event loop
  //////////////////////////////////////////////
  nEvents = treein->GetEntries();
  
  std::cout << "Number of events: " << nEvents << std::endl;

  for (Long64_t jentry=0; jentry<nEvents; jentry++) {
      
      if (jentry % 10 == 0) std::cout << jentry << " processed." << std::endl;

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

      std::cout << "got input particles" << std::endl;
      
      //////////////////////////////////////////////
      // get the resulting jets ordered in pt
      //////////////////////////////////////////////
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.0);
      
      fastjet::ClusterSequence clust_seq_TruthRaw(input_particles, jet_def);
      std::vector<fastjet::PseudoJet> inclusive_jets_TruthRaw = sorted_by_pt(clust_seq_TruthRaw.inclusive_jets(5.0));
      
      if(debug){
	  // label the columns
	  std::cout<<"jet#  pt  eta  phi  mass"<< std::endl;
	  std::cout<<"Inclusive"<< std::endl;
	  // print out the details for each jet
	  for (unsigned int i = 0; i < inclusive_jets_TruthRaw.size(); i++) {
	      std::cout<<i<<"  "<<inclusive_jets_TruthRaw[i].pt()
		       <<"  "<<inclusive_jets_TruthRaw[i].eta()
		       <<"  "<<inclusive_jets_TruthRaw[i].phi()
		       <<"  "<<inclusive_jets_TruthRaw[i].m()<< std::endl;
	  }
      }

    //////////////////////////////////////////////
    //Setup tools for substructure calculation
    //////////////////////////////////////////////

    //Telescoping jets (this looks like the Telescoping reclustering)
    fastjet::contrib::KT_Axes axes_def;
    std::vector<double> r_values;
    int N_r = 20;
    double r_min = 0.1;
    double r_max = 0.6;
    for (int i = 0; i < N_r; i++) {
      r_values.push_back(r_min + i * (r_max - r_min) / (N_r - 1));
    }
    TelescopingJets T_Mass(axes_def,r_values);

    //Energy correlation functions
    fastjet::contrib::EnergyCorrelatorC2 ecfC2(1);
    fastjet::contrib::EnergyCorrelatorD2 ecfD2(1);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecfC3(2, 1);

    // Filtering with a pt cut as for trimming (arXiv:0912.1342)
    double Rfilt0 = 0.3;
    double fcut0 = 0.05;
    fastjet::Transformer *trimmer = new fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, Rfilt0), fastjet::SelectorPtFractionMin(fcut0) );
    const fastjet::Transformer &f = *trimmer;

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

    /////////////////////////////
    //TruthRawTrim
    /////////////////////////////
    for (unsigned int ijet = 0; ijet < inclusive_jets_TruthRaw.size(); ijet++) {
      fastjet::PseudoJet groomed_jet = f(inclusive_jets_TruthRaw[ijet]);

      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(groomed_jet.pt(),
                           groomed_jet.eta(),
                           groomed_jet.phi(),
                           groomed_jet.m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      jetflavor = GetJetTruthFlavor(jettemp, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
      if(debug) std::cout<<"FillingJet Trimmed: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<< std::endl;

      if(jetflavor == -1) continue;

      std::cout << "helllooo" << std::endl;
      
      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      TSub  T1SubOutputTrim  = TNSubjet(groomed_jet, 1, 0.1, 1.0, 40);
      TSub  T2SubOutputTrim  = TNSubjet(groomed_jet, 2, 0.1, 1.0, 40);
      TSub T2Test = T_2Subjet(groomed_jet, 0.1, 1.0, 40);
      T3Sub  T3SubOutputTrim = T_3Subjet(groomed_jet, 0.1, 1.0, 40);

      std::cout << "new min angle: " << T2SubOutputTrim.min_angle << std::endl;
      std::cout << "old min angle: " << T2Test.min_angle << std::endl;
      std::cout << "new volatility: " << T2SubOutputTrim.volatility << std::endl;
      std::cout << "old volatility: " << T2Test.volatility << std::endl;

      
      TruthRawTrim_flavor.push_back(jetflavor);

      TruthRawTrim_pt.push_back(jettemp.Pt());
      TruthRawTrim_eta.push_back(jettemp.Eta());
      TruthRawTrim_phi.push_back(jettemp.Phi());
      TruthRawTrim_m.push_back(jettemp.M());
      
      TruthRawTrim_Tau21.push_back(GetTau21(groomed_jet));
      TruthRawTrim_Tau32.push_back(GetTau32(groomed_jet));
      
      TruthRawTrim_T1jet_angle.push_back(T1SubOutputTrim.min_angle);
      TruthRawTrim_T1jet.push_back(T1SubOutputTrim.volatility);
      TruthRawTrim_T1Volatility.push_back(T1SubOutputTrim.volVec);
      TruthRawTrim_T1masses.push_back(T1SubOutputTrim.masses);
      
      TruthRawTrim_T2jet_angle.push_back(T2SubOutputTrim.min_angle);
      TruthRawTrim_T2jet.push_back(T2SubOutputTrim.volatility);
      TruthRawTrim_T2Volatility.push_back(T2SubOutputTrim.volVec);
      TruthRawTrim_T2masses.push_back(T2SubOutputTrim.masses);
      
      TruthRawTrim_T3jet_angle.push_back(T3SubOutputTrim.min_angle);
      TruthRawTrim_T3jet.push_back(T3SubOutputTrim.volatility);
      TruthRawTrim_T3jet_W.push_back(T3SubOutputTrim.mass_W);
      TruthRawTrim_T3jet_mW.push_back(T3SubOutputTrim.volatility_mass_W);
      TruthRawTrim_T3Volatility.push_back(T3SubOutputTrim.volVec);
      TruthRawTrim_T3masses.push_back(T3SubOutputTrim.volVec);
      
      TruthRawTrim_v32.push_back(T3SubOutputTrim.volatility / T2SubOutputTrim.volatility);
    }


    //////////////////////////////////////
    //Fill event into tree
    //////////////////////////////////////
    if(debug) std::cout<<"Filling Tree"<< std::endl;
    treeout->Fill();
  }


  /////////////////////////////////
  //Write the output TTree to the OutputFile
  /////////////////////////////////
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

  TruthRawTrim_v32.clear();
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
  for (int i = 0; i < nEta; i++)  for (int j = 0; j < nPhi; j++)  tower[i][j] = -0.001;

  std::vector<fastjet::PseudoJet> cell_particles;
  for (int p = 0; p < (int)truth_particles.size(); p++) {
    fastjet::PseudoJet part = truth_particles.at(p);

    int etaCell = int((part.eta()+etaLim)/dEta);
    int phiCell = int(part.phi()/dPhi);
    if (etaCell >= 0 && etaCell < nEta && phiCell >=0 && phiCell < nPhi){
      tower[etaCell][phiCell] += part.e();
    }
  }

  for (int i = 0; i < nEta; i++)  for (int j = 0; j < nPhi; j++) {
    if (tower[i][j] > 0) {
      double etaLocal = -etaLim + dEta*(i+0.5);
      double phiLocal = dPhi*(j+0.5);
      double thetaLocal = 2*atan(exp(-etaLocal));
      cell_particles.push_back(fastjet::PseudoJet(sin(thetaLocal)*cos(phiLocal),sin(thetaLocal)*sin(phiLocal),cos(thetaLocal),1)*tower[i][j]);
    }
  }
  return cell_particles;
}


///=========================================
/// Telescoping Pruning
///=========================================
double T_Pruning(fastjet::PseudoJet& input, double dcut_min, double dcut_max, int N_dcut) {
    double zcut = 0.1; // single choice of zcut. can be further telescoped
    double Tmass = 0;
    std::vector<double> ms; ms.clear();
    for (int i = 0; i < N_dcut; i++){
        double dcut = dcut_min + (dcut_max-dcut_min)*i/(N_dcut-1);
        fastjet::Pruner pruner(fastjet::cambridge_algorithm, zcut, dcut);
        fastjet::PseudoJet pruned_jet = pruner(input);
        Tmass = pruned_jet.m();
        if(Tmass > M0){ms.push_back(Tmass);}
    }
    // getVolatility function provided by TelescopingJets
    if(ms.size()>0)
        return getVolatility(ms);
    
    std::cout <<"WARNING: zero entries for T_Pruning "<< dcut_min <<" "<< dcut_max <<" "<< N_dcut <<" "<<  std::endl;
    return -1;
    
}


///=========================================
/// Telescoping Trimming
///=========================================
double T_Trimming(fastjet::PseudoJet& input, double fcut_min, double fcut_max, int N_fcut) {
    double Rfilt = 0.1; // single choice of Rfilt. can be further telescoped.
    // used Rfilt = 0.1 for higher pT jets and Rfilt = 0.2 for lower pT jets.
    double Tmass = 0;
    std::vector<double> ms; ms.clear();
    for (int i = 0; i < N_fcut; i++){
        double fcut = fcut_min + (fcut_max-fcut_min)*i/(N_fcut-1);
        fastjet::Filter trimmer(Rfilt,fastjet::SelectorPtFractionMin(fcut));
        fastjet::PseudoJet trimmed_jet = trimmer(input);
        Tmass = trimmed_jet.m();
        if(Tmass > M0){ms.push_back(Tmass);}
    }
    // getVolatility function provided by TelescopingJets
    if(ms.size()>0)
        return getVolatility(ms);
    
    std::cout <<"WARNING: zero entries for T_Trimming "<< fcut_min <<" "<< fcut_max <<" "<< N_fcut <<" "<<  std::endl;
    return -1;
    
}


///=========================================
/// Telescoping reclustering
/// for kt, set algorithm = 0
/// for antiki, set algorithm = 2
///=========================================
double T_Reclustering(fastjet::PseudoJet& input, int algorithm, double minRadius, double maxRadius, int numRadii) {
    std::vector<double> telescopingMass; 
    double mass = 0;
    double deltaR = (maxRadius - minRadius) / (numRadii - 1);

    for (double r = minRadius; r <= maxRadius; r += deltaR) {
        fastjet::JetDefinition TjetDefinition(fastjet::JetAlgorithm(algorithm), r);
        fastjet::ClusterSequence TClusterSequence(input.constituents(), TjetDefinition);

        std::vector<fastjet::PseudoJet> recoTJets = sorted_by_pt(TClusterSequence.inclusive_jets());
        if (recoTJets.size() < 1) {
            std::cout << "Warning: recluster number of subjet is " << recoTJets.size() << std::endl;
            continue;
        }
        if (recoTJets.size() == 1) {
            Tsubjet1.SetPxPyPzE(recoTJets[0].px(), recoTJets[0].py(), recoTJets[0].pz(), recoTJets[0].e());
            mass = Tsubjet1.M();
            if (mass > M0) telescopingMass.push_back(mass);
        }
        else if (recoTJets.size() >= 2) {
            Tsubjet1.SetPxPyPzE(recoTJets[0].px(), recoTJets[0].py(), recoTJets[0].pz(), recoTJets[0].e());
            Tsubjet2.SetPxPyPzE(recoTJets[1].px(), recoTJets[1].py(), recoTJets[1].pz(), recoTJets[1].e());
            mass = (Tsubjet1 + Tsubjet2).M();
            if (mass > M0) telescopingMass.push_back(mass);
        }
    }
    
    if (telescopingMass.size() > 0) return getVolatility(telescopingMass);
    
    std::cout << "WARNING: zero entries for T_reclustering." <<
	"\tAlgorithm: " << algorithm << "\tminRadius: " << minRadius << "\tmaxRadius: " << maxRadius << "\tnumRadii: " << numRadii <<  std::endl;
    return -1;
}

///=========================================
/// Telescoping Subjet
///=========================================
TSub TNSubjet(fastjet::PseudoJet& input, unsigned int numSubjets, double minRadius, double maxRadius, int numRadii) {
    TSub result;

    std::vector<double> telescopingMasses;
    double mass = 0;
    double beta = 1.0;

    fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSubjettiness(numSubjets, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);

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

    double minTauAxisAngle = 0;
    if (numSubjets > 1) {
	auto minTauAxisAngleIt = std::min_element(anglesBetweenTauAxes.cbegin(), anglesBetweenTauAxes.cend());
	minTauAxisAngle = *minTauAxisAngleIt;
    }
    
    double deltaR = (maxRadius - minRadius) / (numRadii);

    std::vector<fastjet::PseudoJet> constituents = input.constituents();
    std::vector<TLorentzVector> TSubjets(numSubjets);

    for (double r = minRadius; r <= maxRadius; r += deltaR) {
	for (auto it = constituents.cbegin(); it != constituents.cend(); ++it) {
	    const fastjet::PseudoJet constituent = *it;
	    TLorentzVector pConstituent(constituent.px(), constituent.py(), constituent.pz(), constituent.e());
	    
	    std::vector<double> distanceToTauAxes(numSubjets);
	    for (unsigned int i = 0; i < numSubjets; ++i) {
		distanceToTauAxes[i] = pConstituent.DeltaR(pTauAxes[i]);
	    }
	    
	    auto minDistanceToTauAxesIt = std::min_element(distanceToTauAxes.cbegin(), distanceToTauAxes.cend());
	    double minDistanceToTauAxes = *minDistanceToTauAxesIt;
	    
	    if (minDistanceToTauAxes <= r) {
		unsigned int minDistanceToTauAxesIndex = std::distance(distanceToTauAxes.cbegin(), minDistanceToTauAxesIt);
		TSubjets[minDistanceToTauAxesIndex] += pConstituent;

		constituents.erase(it);
		--it;
	    }
	}
	
	TLorentzVector sumTSubjet;
	for (auto const &TSubjet : TSubjets) {
	    sumTSubjet += TSubjet;
	}
	mass = sumTSubjet.M();

	if (mass > M0) {
	    telescopingMasses.push_back(mass);
	    result.volVec.push_back(getVolatility(telescopingMasses));
	}
    }
    
    result.min_angle = minTauAxisAngle;
    if (telescopingMasses.size() > 0) {
	result.volatility = getVolatility(telescopingMasses);
	result.masses = telescopingMasses;
    }
    else {
	std::cout << "WARNING: zero entries for TNSubjet: \tNumSubjets: " << numSubjets << "\tminRadius: " << minRadius << "\tmaxRadius: " << maxRadius << "\tnumRadii: " << numRadii <<  std::endl;
        result.volatility = -1;
    }
    return result;
}

T3Sub T_3Subjet(fastjet::PseudoJet& input, double R_min, double R_max, int N_R){
    T3Sub result;
    std::vector<double> m123s; m123s.clear();
    std::vector<double> m12s; m12s.clear();
    std::vector<double> m13s; m13s.clear();
    std::vector<double> m23s; m23s.clear();
    // m123 is the invariant mass of the three subjets
    double beta = 1.0;
    double Tmass = 0;
    double Tmass_12 = 0;
    double Tmass_13 = 0;
    double Tmass_23 = 0;
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSub(3, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
    double tau3 = nSub.result(input);
    std::vector<fastjet::PseudoJet> tau3axes = nSub.currentAxes();
    tau_axis1.SetPxPyPzE(tau3axes[0].px(),tau3axes[0].py(),tau3axes[0].pz(),tau3axes[0].e());
    tau_axis2.SetPxPyPzE(tau3axes[1].px(),tau3axes[1].py(),tau3axes[1].pz(),tau3axes[1].e());
    tau_axis3.SetPxPyPzE(tau3axes[2].px(),tau3axes[2].py(),tau3axes[2].pz(),tau3axes[2].e());
    double tau3_D12 = tau_axis1.DeltaR(tau_axis2);
    double tau3_D13 = tau_axis1.DeltaR(tau_axis3);
    double tau3_D23 = tau_axis2.DeltaR(tau_axis3);
    double D_min = tau3_D12;
    double D_mid = tau3_D13;
    double D_max = tau3_D23;
    double D_temp;
    if(D_mid < D_min ) {D_temp = D_min; D_min = D_mid; D_mid = D_temp;}
    if(D_max < D_min ) {D_temp = D_min; D_min = D_max; D_max = D_temp;}
    if(D_max < D_mid ) {D_temp = D_mid; D_mid = D_max; D_max = D_temp;}
    
    double d1, d2, d3;
    double deltaR = (R_max - R_min)/(N_R-1);
    double R      = R_min;
    
    for (int i = 0; i < N_R; i++){
	//    R = R_min + i*deltaR;
	Tsubjet1.SetPxPyPzE(0,0,0,0);
	Tsubjet2.SetPxPyPzE(0,0,0,0);
	Tsubjet3.SetPxPyPzE(0,0,0,0);
	
	for (unsigned int c = 0; c < input.constituents().size(); c++) {
	    particle.SetPxPyPzE(input.constituents()[c].px(),input.constituents()[c].py(),input.constituents()[c].pz(),input.constituents()[c].e());
	    d1 = particle.DeltaR(tau_axis1);
	    d2 = particle.DeltaR(tau_axis2);
	    d3 = particle.DeltaR(tau_axis3);
	    if (d1 <= R && d1 < d2 && d1 < d3){
		Tsubjet1 = Tsubjet1 + particle;
	    }
	    else if (d2 <= R && d2 < d1 && d2 < d3){
		Tsubjet2 = Tsubjet2 + particle;
	    }
	    else if (d3 <= R && d3 < d1 && d3 < d2){
		Tsubjet3 = Tsubjet3 + particle;
	    }
	}
	Tmass = (Tsubjet1 + Tsubjet2 + Tsubjet3).M();
	Tmass_12 = (Tsubjet1 + Tsubjet2).M();
	Tmass_13 = (Tsubjet1 + Tsubjet3).M();
	Tmass_23 = (Tsubjet2 + Tsubjet3).M();
	if(Tmass > M0){
	    m123s.push_back(Tmass);
	    result.volVec.push_back(getVolatility(m123s));
	}
	if(Tmass_12 > M0){m12s.push_back(Tmass_12);}
	if(Tmass_13 > M0){m13s.push_back(Tmass_13);}
	if(Tmass_23 > M0){m23s.push_back(Tmass_23);}
	R += deltaR;
    }
    result.min_angle = D_min;
    result.mid_angle = D_mid;
    result.max_angle = D_max;
    
    if(m123s.size() > 0 && m12s.size() > 0 && m13s.size() > 0 && m23s.size() > 0){
        result.volatility = getVolatility(m123s);
	result.masses = m123s;
        if ( std::abs (Tmass_12 - mW) < std::abs (Tmass_13 - mW) && std::abs (Tmass_12 - mW) < std::abs (Tmass_23 - mW)){
            result.mass_W = Tmass_12;
            result.volatility_mass_W = getVolatility(m12s);
        }
        else if ( std::abs (Tmass_13 - mW) < std::abs (Tmass_12 - mW) && std::abs (Tmass_13 - mW) < std::abs (Tmass_23 - mW)){
            result.mass_W = Tmass_13;
            result.volatility_mass_W = getVolatility(m13s);
        }
        else if ( std::abs (Tmass_23 - mW) < std::abs (Tmass_12 - mW) && std::abs (Tmass_23 - mW) < std::abs (Tmass_13 - mW)){
            result.mass_W = Tmass_23;
            result.volatility_mass_W = getVolatility(m23s);
        }
    }else{
        std::cout <<"WARNING: zero entries for T_3Subjet "<< R_min <<" "<< R_max <<" "<< N_R <<" "<<  std::endl;
	result.mass_W = -1;
	result.volatility = -1;
	result.volatility_mass_W = -1;
    }
    
    return result;
}

///=========================================
/// Telescoping N-subjettiness
///=========================================
double T_Nsubjettiness(int N, fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  std::vector<double> taus; taus.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
    fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
//    fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    taus.push_back(nsub(input));
  }
  // getVolatility function provided by TelescopingJets
    if(taus.size()>0)
        return getVolatility(taus);
    
    std::cout <<"WARNING: zero entries for T_Nsubjettiness "<< beta_min <<" "<< beta_max <<" "<< N_beta <<" "<<  std::endl;
    return -1;
    
}

double T_NsubjettinessRatio(int N_num, int N_den, fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  std::vector<double> taus; taus.clear();
  for (int i = 0; i < N_beta; i++) {

    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);

    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);

    fastjet::contrib::Nsubjettiness nsub_num(N_num, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    fastjet::contrib::Nsubjettiness nsub_den(N_den, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);

    double num=nsub_num(input);
    double den=nsub_den(input);

    if(den!=0)
      taus.push_back(num/den);
    else
      taus.push_back(-1.0);

  }
    if(taus.size()==0){
        std::cout <<"WARNING: zero entries for T_NsubjetinessRatio "<< beta_min <<" "<< beta_max <<" "<< N_beta <<" "<<  std::endl;
        return -1;
    }
    // getVolatility function provided by TelescopingJets
    if(taus.size()>0)
        return getVolatility(taus);
    
    std::cout <<"WARNING: zero entries for T_NsubjettinessRatio "<< beta_min <<" "<< beta_max <<" "<< N_beta <<" "<<  std::endl;
    return -1;
    
}


///=========================================
/// Telescoping Energy Correlators
///=========================================
double T_EnergyCorrelator_C2(fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  std::vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::EnergyCorrelatorC2 ecf(beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
    if(ecfs.size()>0)
        return getVolatility(ecfs);
    
    std::cout <<"WARNING: zero entries for T_EnergyCorrelator_C2 "<< beta_min <<" "<< beta_max <<" "<< N_beta <<" "<<  std::endl;
    return -1;
    
}

double T_EnergyCorrelator_D2(fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  std::vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::EnergyCorrelatorD2 ecf(beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
    if(ecfs.size()>0)
        return getVolatility(ecfs);
    
    std::cout <<"WARNING: zero entries for T_EnergyCorrelator_C2 "<< beta_min <<" "<< beta_max <<" "<< N_beta <<" "<<  std::endl;
    return -1;
    
}

double T_EnergyCorrelator_C3(fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  std::vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecf(3, beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
    if(ecfs.size()>0)
        return getVolatility(ecfs);
    
    std::cout <<"WARNING: zero entries for T_EnergyCorrelator_C3 "<< beta_min <<" "<< beta_max <<" "<< N_beta <<" "<<  std::endl;
    return -1;
    
}

///========================================
int GetJetTruthFlavor(TLorentzVector jettemp,
                      TLorentzVector truth_t1,
                      TLorentzVector truth_t2,
                      TLorentzVector truth_W1,
                      TLorentzVector truth_W2,
                      TLorentzVector truth_H,
                      int debug){
  if(debug){
    std::cout<<"DeltaR:   "<< std::endl
        <<"dRMatch:  "<<dR_match<< std::endl
        <<"q1:       "<<jettemp.DeltaR(truth_q1)<< std::endl
        <<"q2:       "<<jettemp.DeltaR(truth_q2)<< std::endl
        <<"W1:       "<<jettemp.DeltaR(truth_W1)<< std::endl
        <<"W2:       "<<jettemp.DeltaR(truth_W2)<< std::endl
        <<"H:        "<<jettemp.DeltaR(truth_H)<< std::endl
        <<"t1:       "<<jettemp.DeltaR(truth_t1)<< std::endl
        <<"t2:       "<<jettemp.DeltaR(truth_t2)<< std::endl;
  }
  int jetflavor = -1;
  if     (jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){
    jetflavor = 0;
  }
  else if(jettemp.DeltaR(truth_W1)<dR_match || jettemp.DeltaR(truth_W2)<dR_match){
    jetflavor = 1;
  }
  else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){
    jetflavor = 3;
  }
  else if(jettemp.DeltaR(truth_H)<dR_match){
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

  if(tau2>0)
    tau32 = tau3/tau2;

  return tau32;

}

