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


#include "NTuplerTest.h"

int main(int argc, char* argv[]){

  //exit if you dont pass a run card
  if(argc<4){
    cout<<"You need to specify more arguments"<<endl;
    cout<<"Arg1 = process type"<<endl;
    cout<<"Arg2 = input file path and name"<<endl;
    cout<<"Arg3 = output file path and name"<<endl;
    cout<<"Arg4 = debug flag (optional)"<<endl;
    return 1;
  }

  //process
  string ProcessType = argv[1];

  //inputfile
  string InputFile = argv[2];

  //outputfile
  string OutputFile = argv[3];

  //debug flag
  bool debug=false;
  if(argc>=5){
    string argdebug = argv[4];
    if(argdebug=="debug")
      debug=true;
  }

  //print out the input arguments
  cout<<"InputArguments:  ProcessType="<<ProcessType<<"   InputFile="<<InputFile<<"  OutputFile="<<OutputFile<<"  Debug="<<debug<<endl;

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
  treeout->Branch("TruthRaw_T3jet",         &TruthRaw_T3jet);
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
  treeout->Branch("TruthRawTrim_T3jet",         &TruthRawTrim_T3jet);
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

  //////////////////////////////////////////////
  //test NTupler.cc using Alex’s samples
  //////////////////////////////////////////////
  //read N tuple from root files
  //////////////////////////////////////////////

  TChain *sig_t1 = new TChain("Tree");
  TChain *sig_t2 = new TChain("Tree");
  TChain *bkg_t = new TChain("Tree");

//  sig_t1->Add("AlexSample/mc12_8TeV.110903.Pythia8_AU2MSTW2008LO_zprime1000_tt.leading_jets.root");
//  sig_t1->Add("AlexSample/mc12_8TeV.110907.Pythia8_AU2MSTW2008LO_zprime2000_tt.leading_jets.root");

  sig_t1->Add("AlexSample/mc12_8TeV.158864.Pythia8_AU2MSTW2008LO_Wprime_WZ_llqq_m1000.leading_jets.root");
//  sig_t1->Add("AlexSample/mc12_8TeV.158874.Pythia8_AU2MSTW2008LO_Wprime_WZ_llqq_m2000.leading_jets.root");
  bkg_t->Add("AlexSample/mc12_8TeV.14791X.Pythia8_AU2CT10_jetjet_JZXW.leading_jets.root");

    
//  data_file.open("/Users/ytchien/Research/Deep_learning/Discretized_Jets_new/discretized_gluons_200GeV_100k.txt");
//  Tjet_variable_file.open("/Users/ytchien/Research/Deep_learning/Tjet_variables/Tjet_gluon_4J_test.txt");

  Tjet_variable_file.open("/Users/ytchien/Research/Deep_learning/Tjet_variables/Tjet_W_3J.txt");
//  Tjet_variable_file.open("/Users/ytchien/Research/Deep_learning/Tjet_variables/Tjet_QCD_3J.txt");
    
    
  //////////////////////////////////////////////
  //define variables to read the tree information
  //////////////////////////////////////////////

  Int_t sigEventNumber1;
  Int_t ChannelNumber1;
  Double_t sigEventWeight1, sigCrossSection1;//, sigPileupWeight1;
  Bool_t PassEventSelection1;

  int sigTruth_n1;//, sig_dRmatch1;
  float sig_pt1, sig_eta1, sig_phi1, sig_m1;
  vector<float> *sigTruth_pt1 = 0;
  vector<float> *sigTruth_eta1 = 0;
  vector<float> *sigTruth_phi1 = 0;
  vector<float> *sigTruth_m1 = 0;
//  vector<int> *sigTruth_status1 = 0;
//  vector<int> *sigTruth_pdgId1 = 0;
  int sigReco_n1;//, sig_dRmatch1;
  float sigR_pt1, sigR_eta1, sigR_phi1, sigR_m1;
  vector<float> *sigReco_pt1 = 0;
  vector<float> *sigReco_eta1 = 0;
  vector<float> *sigReco_phi1 = 0;
  vector<float> *sigReco_m1 = 0;
//  vector<int> *sigReco_status1 = 0;
//  vector<int> *sigReco_pdgId1 = 0;

  Float_t	Truth_W_pt1;
  Float_t	Truth_W_eta1;
  Float_t	Truth_W_phi1;
  Float_t	Truth_W_m1;
  Int_t 	Truth_W_status1;
  Int_t		Truth_W_pdgId1;	
  Float_t	Truth_T_pt1;
  Float_t	Truth_T_eta1;
  Float_t	Truth_T_phi1;
  Float_t	Truth_T_m1;
  Int_t 	Truth_T_status1;
  Int_t		Truth_T_pdgId1;	

  vector<float>	*Truth_Wdaughter_pt1 = 0;
  vector<float>	*Truth_Wdaughter_eta1 = 0;
  vector<float>	*Truth_Wdaughter_phi1 = 0;
  vector<float>	*Truth_Wdaughter_m1 = 0;
  vector<int>	*Truth_Wdaughter_status1 = 0;
  vector<int>	*Truth_Wdaughter_pdgId1 = 0;

  vector<float>	*Truth_Tdaughter_pt1 = 0;
  vector<float>	*Truth_Tdaughter_eta1 = 0;
  vector<float>	*Truth_Tdaughter_phi1 = 0;
  vector<float>	*Truth_Tdaughter_m1 = 0;
  vector<int>	*Truth_Tdaughter_status1 = 0;
  vector<int>	*Truth_Tdaughter_pdgId1 = 0;

  Int_t sigEventNumber2;
  Int_t ChannelNumber2;
  Double_t sigEventWeight2, sigCrossSection2;//, sigPileupWeight2;
  Bool_t PassEventSelection2;

  int sigTruth_n2;//, sig_dRmatch2;
  float sig_pt2, sig_eta2, sig_phi2, sig_m2;
  vector<float> *sigTruth_pt2 = 0;
  vector<float> *sigTruth_eta2 = 0;
  vector<float> *sigTruth_phi2 = 0;
  vector<float> *sigTruth_m2 = 0;
//  vector<int> *sigTruth_status2 = 0;
//  vector<int> *sigTruth_pdgId2 = 0;
  int sigReco_n2;//, sig_dRmatch2;
  float sigR_pt2, sigR_eta2, sigR_phi2, sigR_m2;
  vector<float> *sigReco_pt2 = 0;
  vector<float> *sigReco_eta2 = 0;
  vector<float> *sigReco_phi2 = 0;
  vector<float> *sigReco_m2 = 0;
//  vector<int> *sigReco_status2 = 0;
//  vector<int> *sigReco_pdgId2 = 0;

  Float_t	Truth_W_pt2;
  Float_t	Truth_W_eta2;
  Float_t	Truth_W_phi2;
  Float_t	Truth_W_m2;
  Int_t 	Truth_W_status2;
  Int_t		Truth_W_pdgId2;	
  Float_t	Truth_T_pt2;
  Float_t	Truth_T_eta2;
  Float_t	Truth_T_phi2;
  Float_t	Truth_T_m2;
  Int_t 	Truth_T_status2;
  Int_t		Truth_T_pdgId2;	

  vector<float>	*Truth_Wdaughter_pt2 = 0;
  vector<float>	*Truth_Wdaughter_eta2 = 0;
  vector<float>	*Truth_Wdaughter_phi2 = 0;
  vector<float>	*Truth_Wdaughter_m2 = 0;
  vector<int>	*Truth_Wdaughter_status2 = 0;
  vector<int>	*Truth_Wdaughter_pdgId2 = 0;

  vector<float>	*Truth_Tdaughter_pt2 = 0;
  vector<float>	*Truth_Tdaughter_eta2 = 0;
  vector<float>	*Truth_Tdaughter_phi2 = 0;
  vector<float>	*Truth_Tdaughter_m2 = 0;
  vector<int>	*Truth_Tdaughter_status2 = 0;
  vector<int>	*Truth_Tdaughter_pdgId2 = 0;

  Int_t bkgEventNumber;
  Double_t bkgEventWeight, bkgCrossSection, bkgPileupWeight;

  int bkgTruth_n;//, bkg_dRmatch;
  float bkg_pt, bkg_eta, bkg_phi, bkg_m;
  vector<float> *bkgTruth_pt = 0;
  vector<float> *bkgTruth_eta = 0;
  vector<float> *bkgTruth_phi = 0;
  vector<float> *bkgTruth_m = 0;
//  vector<int> *bkgTruth_status = 0;
//  vector<int> *bkgTruth_pdgId = 0;
  int bkgReco_n;//, bkg_dRmatch;
  float bkgR_pt, bkgR_eta, bkgR_phi, bkgR_m;
  vector<float> *bkgReco_pt = 0;
  vector<float> *bkgReco_eta = 0;
  vector<float> *bkgReco_phi = 0;
  vector<float> *bkgReco_m = 0;
//  vector<int> *bkgReco_status = 0;
//  vector<int> *bkgReco_pdgId = 0;

  //////////////////////////////////////////////
  //read the branch information
  //////////////////////////////////////////////

  sig_t1->SetBranchAddress("EventNumber",	  	 &sigEventNumber1);
  sig_t1->SetBranchAddress("EventWeight",	  	 &sigEventWeight1);
  sig_t1->SetBranchAddress("CrossSection",	  	 &sigCrossSection1);
//  sig_t1->SetBranchAddress("PileupWeight",	  	 &sigPileupWeight1);
  sig_t1->SetBranchAddress("ChannelNumber",	  	 &ChannelNumber1);
  sig_t1->SetBranchAddress("PassEventSelection",	 &PassEventSelection1);

  sig_t1->SetBranchAddress("AntiKt10Truth_pt",		 &sig_pt1);
  sig_t1->SetBranchAddress("AntiKt10Truth_eta",		 &sig_eta1);
  sig_t1->SetBranchAddress("AntiKt10Truth_phi",		 &sig_phi1);
  sig_t1->SetBranchAddress("AntiKt10Truth_m",		 &sig_m1);
  sig_t1->SetBranchAddress("AntiKt10Truth_constit_n",	 &sigTruth_n1);
  sig_t1->SetBranchAddress("AntiKt10Truth_constit_pt",	 &sigTruth_pt1);
  sig_t1->SetBranchAddress("AntiKt10Truth_constit_eta",	 &sigTruth_eta1);
  sig_t1->SetBranchAddress("AntiKt10Truth_constit_phi",	 &sigTruth_phi1);
  sig_t1->SetBranchAddress("AntiKt10Truth_constit_m",	 &sigTruth_m1);

  sig_t1->SetBranchAddress("AntiKt10Reco_pt",		 &sigR_pt1);
  sig_t1->SetBranchAddress("AntiKt10Reco_eta",		 &sigR_eta1);
  sig_t1->SetBranchAddress("AntiKt10Reco_phi",		 &sigR_phi1);
  sig_t1->SetBranchAddress("AntiKt10Reco_m",		 &sigR_m1);
  sig_t1->SetBranchAddress("AntiKt10Reco_constit_n",	 &sigReco_n1);
  sig_t1->SetBranchAddress("AntiKt10Reco_constit_pt",	 &sigReco_pt1);
  sig_t1->SetBranchAddress("AntiKt10Reco_constit_eta",	 &sigReco_eta1);
  sig_t1->SetBranchAddress("AntiKt10Reco_constit_phi",	 &sigReco_phi1);
  sig_t1->SetBranchAddress("AntiKt10Reco_constit_m",	 &sigReco_m1);

/*sig_t1->SetBranchAddress("AntiKt10LCTopo_pt",		 &sig_pt1);
  sig_t1->SetBranchAddress("AntiKt10LCTopo_eta",	 &sig_eta1);
  sig_t1->SetBranchAddress("AntiKt10LCTopo_phi",	 &sig_phi1);
  sig_t1->SetBranchAddress("AntiKt10LCTopo_m",		 &sig_m1);
  sig_t1->SetBranchAddress("AntiKt10LCTopo_constit_n",	 &sigTruth_n1);
  sig_t1->SetBranchAddress("AntiKt10LCTopo_constit_pt",	 &sigTruth_pt1);
  sig_t1->SetBranchAddress("AntiKt10LCTopo_constit_eta", &sigTruth_eta1);
  sig_t1->SetBranchAddress("AntiKt10LCTopo_constit_phi", &sigTruth_phi1);
  sig_t1->SetBranchAddress("AntiKt10LCTopo_constit_m",	 &sigTruth_m1);*/

  sig_t1->SetBranchAddress("Truth_W_pt",	  	 &Truth_W_pt1);
  sig_t1->SetBranchAddress("Truth_W_eta",	  	 &Truth_W_eta1);
  sig_t1->SetBranchAddress("Truth_W_phi",	  	 &Truth_W_phi1);
  sig_t1->SetBranchAddress("Truth_W_m",		 	 &Truth_W_m1);
  sig_t1->SetBranchAddress("Truth_W_status",		 &Truth_W_status1);
  sig_t1->SetBranchAddress("Truth_W_pdgId",		 &Truth_W_pdgId1);
  sig_t1->SetBranchAddress("Truth_Wdaughter_pt",	 &Truth_Wdaughter_pt1);
  sig_t1->SetBranchAddress("Truth_Wdaughter_eta",	 &Truth_Wdaughter_eta1);
  sig_t1->SetBranchAddress("Truth_Wdaughter_phi",	 &Truth_Wdaughter_phi1);
  sig_t1->SetBranchAddress("Truth_Wdaughter_m",		 &Truth_Wdaughter_m1);
  sig_t1->SetBranchAddress("Truth_Wdaughter_status",     &Truth_Wdaughter_status1);
  sig_t1->SetBranchAddress("Truth_Wdaughter_pdgId",      &Truth_Wdaughter_pdgId1);

/*sig_t1->SetBranchAddress("Truth_T_pt",	 	 &Truth_T_pt1);
  sig_t1->SetBranchAddress("Truth_T_eta",	  	 &Truth_T_eta1);
  sig_t1->SetBranchAddress("Truth_T_phi",	  	 &Truth_T_phi1);
  sig_t1->SetBranchAddress("Truth_T_m",	 		 &Truth_T_m1);
  sig_t1->SetBranchAddress("Truth_T_status",	 	 &Truth_T_status1);
  sig_t1->SetBranchAddress("Truth_T_pdgId",	 	 &Truth_T_pdgId1);
  sig_t1->SetBranchAddress("Truth_Tdaughter_pt",	 &Truth_Tdaughter_pt1);
  sig_t1->SetBranchAddress("Truth_Tdaughter_eta",	 &Truth_Tdaughter_eta1);
  sig_t1->SetBranchAddress("Truth_Tdaughter_phi",	 &Truth_Tdaughter_phi1);
  sig_t1->SetBranchAddress("Truth_Tdaughter_m",		 &Truth_Tdaughter_m1);
  sig_t1->SetBranchAddress("Truth_Tdaughter_status",     &Truth_Tdaughter_status1);
  sig_t1->SetBranchAddress("Truth_Tdaughter_pdgId",      &Truth_Tdaughter_pdgId1);*/

  bkg_t->SetBranchAddress("EventNumber",	  	 &bkgEventNumber);
  bkg_t->SetBranchAddress("EventWeight",	  	 &bkgEventWeight);
  bkg_t->SetBranchAddress("CrossSection",	  	 &bkgCrossSection);
//  bkg_t->SetBranchAddress("PileupWeight",	  	 &bkgPileupWeight);

  bkg_t->SetBranchAddress("AntiKt10Truth_pt",		 &bkg_pt);
  bkg_t->SetBranchAddress("AntiKt10Truth_eta",		 &bkg_eta);
  bkg_t->SetBranchAddress("AntiKt10Truth_phi",		 &bkg_phi);
  bkg_t->SetBranchAddress("AntiKt10Truth_m",		 &bkg_m);
  bkg_t->SetBranchAddress("AntiKt10Truth_constit_n",	 &bkgTruth_n);
  bkg_t->SetBranchAddress("AntiKt10Truth_constit_pt",	 &bkgTruth_pt);
  bkg_t->SetBranchAddress("AntiKt10Truth_constit_eta",	 &bkgTruth_eta);
  bkg_t->SetBranchAddress("AntiKt10Truth_constit_phi",	 &bkgTruth_phi);
  bkg_t->SetBranchAddress("AntiKt10Truth_constit_m",	 &bkgTruth_m);

  bkg_t->SetBranchAddress("AntiKt10Reco_pt",		 &bkgR_pt);
  bkg_t->SetBranchAddress("AntiKt10Reco_eta",		 &bkgR_eta);
  bkg_t->SetBranchAddress("AntiKt10Reco_phi",		 &bkgR_phi);
  bkg_t->SetBranchAddress("AntiKt10Reco_m",		 &bkgR_m);
  bkg_t->SetBranchAddress("AntiKt10Reco_constit_n",	 &bkgReco_n);
  bkg_t->SetBranchAddress("AntiKt10Reco_constit_pt",	 &bkgReco_pt);
  bkg_t->SetBranchAddress("AntiKt10Reco_constit_eta",	 &bkgReco_eta);
  bkg_t->SetBranchAddress("AntiKt10Reco_constit_phi",	 &bkgReco_phi);
  bkg_t->SetBranchAddress("AntiKt10Reco_constit_m",	 &bkgReco_m);

/*bkg_t->SetBranchAddress("AntiKt10LCTopo_pt",		 &bkg_pt);
  bkg_t->SetBranchAddress("AntiKt10LCTopo_eta",		 &bkg_eta);
  bkg_t->SetBranchAddress("AntiKt10LCTopo_phi",		 &bkg_phi);
  bkg_t->SetBranchAddress("AntiKt10LCTopo_m",		 &bkg_m);
  bkg_t->SetBranchAddress("AntiKt10LCTopo_constit_n",	 &bkgTruth_n);
  bkg_t->SetBranchAddress("AntiKt10LCTopo_constit_pt",	 &bkgTruth_pt);
  bkg_t->SetBranchAddress("AntiKt10LCTopo_constit_eta",	 &bkgTruth_eta);
  bkg_t->SetBranchAddress("AntiKt10LCTopo_constit_phi",	 &bkgTruth_phi);
  bkg_t->SetBranchAddress("AntiKt10LCTopo_constit_m",	 &bkgTruth_m);*/

  TH1D *Tprun_volatility_sigPlot = new TH1D("Histogram","T-pruning Volatility", 50, 0, 1.0); 
  TH1D *Tprun_volatility_bkgPlot = new TH1D("Histogram","T-pruning Volatility", 50, 0, 1.0); 
  TH1D *Ttrim_volatility_sigPlot = new TH1D("Histogram","T-trimming Volatility", 50, 0, 1.0); 
  TH1D *Ttrim_volatility_bkgPlot = new TH1D("Histogram","T-trimming Volatility", 50, 0, 1.0); 
  TH1D *TakTrecl_volatility_sigPlot = new TH1D("Histogram","T-antikT reclustering Volatility", 50, 0, 1.0); 
  TH1D *TakTrecl_volatility_bkgPlot = new TH1D("Histogram","T-antikT reclustering Volatility", 50, 0, 1.0);
  TH1D *TkTrecl_volatility_sigPlot = new TH1D("Histogram","T-kT reclustering Volatility", 50, 0, 1.0); 
  TH1D *TkTrecl_volatility_bkgPlot = new TH1D("Histogram","T-kT reclustering Volatility", 50, 0, 1.0);  
  TH1D *Tsubj_volatility_sigPlot = new TH1D("Histogram","T-subjet Volatility", 50, 0, 1.0); 
  TH1D *Tsubj_volatility_bkgPlot = new TH1D("Histogram","T-subjet Volatility", 50, 0, 1.0); 
  TH1D *Tsubj_angle_sigPlot = new TH1D("Histogram","T-subjet Angle", 50, 0, 1.0); 
  TH1D *Tsubj_angle_bkgPlot = new TH1D("Histogram","T-subjet Angle", 50, 0, 1.0); 
  TH1D *Ttau2_volatility_sigPlot = new TH1D("Histogram","T-2subjettiness Volatility", 50, 0, 3.0); 
  TH1D *Ttau2_volatility_bkgPlot = new TH1D("Histogram","T-2subjettiness Volatility", 50, 0, 3.0); 
  TH1D *tau21_sigPlot = new TH1D("Histogram","tau21", 50, 0, 1.0); 
  TH1D *tau21_bkgPlot = new TH1D("Histogram","tau21", 50, 0, 1.0);

/*
    //////////////////////////////////////////////////
    //test telescoping deconstruction in q/g tagging
    //////////////////////////////////////////////////

    int id;
    int N_jet = 0;
    int N_particle = 0;
    double px, py, pz, e;
    
    while(1){
        data_file >> N_particle;
        if (!data_file.good()) break;
        N_jet = N_jet + 1;
        
        vector<fastjet::PseudoJet> input_particles;
        
        for (int n = 0; n < N_particle; n++){
            data_file >> px >> py >> pz >> e;
            //   cout << px << py << pz << e << endl;
            fastjet::PseudoJet part(px,py,pz,e);
            input_particles.push_back(part);
            input_particles.back().set_user_index(id);
        }
        //   cout << endl;
        //   if (N_jet == 100000) break;
        
        double R0 = 1.0;
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R0);
        fastjet::ClusterSequence cs(input_particles, jet_def);
        vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

        TSub  T1SubOutputTrim = T_1Subjet(jets[0], 0.025, 0.3, 12);
        TSub  T2SubOutputTrim = T_2Subjet(jets[0], 0.025, 0.3, 12);
        TSub  T3SubOutputTrim = T_3Subjet(jets[0], 0.025, 0.3, 12);
        TSub  T4SubOutputTrim = T_4Subjet(jets[0], 0.025, 0.3, 12);
        Tjet_variable_file << endl;
    }
*/
    

  //////////////////////////////////////////////
  //signal
  //////////////////////////////////////////////

for (int event = 0; event < sig_t1->GetEntries(); event++) {
    sig_t1->GetEntry(event);

//	Truth_T.SetPtEtaPhiM(Truth_T_pt1,Truth_T_eta1,Truth_T_phi1,Truth_T_m1);
//	Tdaughter1.SetPtEtaPhiM(Truth_Tdaughter_pt1->at(0),Truth_Tdaughter_eta1->at(0),Truth_Tdaughter_phi1->at(0),Truth_Tdaughter_m1->at(0));
//	Tdaughter2.SetPtEtaPhiM(Truth_Tdaughter_pt1->at(1),Truth_Tdaughter_eta1->at(1),Truth_Tdaughter_phi1->at(1),Truth_Tdaughter_m1->at(1));
//	fastjet::PseudoJet truth_t(Truth_T.Px(),Truth_T.Py(),Truth_T.Pz(),Truth_T.E());

	Truth_W.SetPtEtaPhiM(Truth_W_pt1,Truth_W_eta1,Truth_W_phi1,Truth_W_m1);
	Wdaughter1.SetPtEtaPhiM(Truth_Wdaughter_pt1->at(0),Truth_Wdaughter_eta1->at(0),Truth_Wdaughter_phi1->at(0),Truth_Wdaughter_m1->at(0));
	Wdaughter2.SetPtEtaPhiM(Truth_Wdaughter_pt1->at(1),Truth_Wdaughter_eta1->at(1),Truth_Wdaughter_phi1->at(1),Truth_Wdaughter_m1->at(1));
	fastjet::PseudoJet truth_w(Truth_W.Px(),Truth_W.Py(),Truth_W.Pz(),Truth_W.E());

	double weight = sigEventWeight1 * sigCrossSection1;// * sigPileupWeight1;

	total_truth.SetPxPyPzE(0,0,0,0);//to add up truth particles
	total_reco.SetPxPyPzE(0,0,0,0);//to add up reco particles
	total0_truth.SetPtEtaPhiM(sig_pt1,sig_eta1,sig_phi1,sig_m1);//use default truth jet information

	vector<fastjet::PseudoJet> input_particles_truth;//to push back truth particles
	vector<fastjet::PseudoJet> input_particles_reco;//to push back reco particles

    for (int c = 0; c < sigTruth_n1; c++) {
	particle.SetPtEtaPhiM(sigTruth_pt1->at(c),sigTruth_eta1->at(c),sigTruth_phi1->at(c),sigTruth_m1->at(c));
	total_truth = total_truth + particle;
	fastjet::PseudoJet part(particle.Px(),particle.Py(),particle.Pz(),particle.E());
	input_particles_truth.push_back(part);
	input_particles_truth.back().set_user_index(id);
//	if (event == sig_event_select){
//		EventDisplay->Fill(particle.Eta(),particle.Phi(),particle.E());
//	}
    }// read in .root file and push back the truth particles

    for (int c = 0; c < sigReco_n1; c++) {
	particle.SetPtEtaPhiM(sigReco_pt1->at(c),sigReco_eta1->at(c),sigReco_phi1->at(c),sigReco_m1->at(c));
	total_reco = total_reco + particle;
	fastjet::PseudoJet part(particle.Px(),particle.Py(),particle.Pz(),particle.E());
	input_particles_reco.push_back(part);
	input_particles_reco.back().set_user_index(id);
//	if (event == sig_event_select){
//		EventDisplay->Fill(particle.Eta(),particle.Phi(),particle.E());
//	}
    }// read in .root file and push back the reco particles

//	cout << "Channel " << ChannelNumber << " jet " << event << endl;
//	cout << "Truth T " << Truth_T_pt << " " << Truth_T_eta << " " << Truth_T_phi << " " << Truth_T_m << endl;
//	cout << "Truth W " << Truth_W_pt << " " << Truth_W_eta << " " << Truth_W_phi << " " << Truth_W_m << endl;
//	cout << "mass " << total.M() << endl;
//	cout << "number " << sigTruth_n << endl;

      //////////////////////////////////////////////////////////////////////////////////////
      //if(input_particles_truth.size() == 0){
      //	cout << "jet" << event << " has 0 input particles" << endl;
      //	continue;
      //} // tested. All jets have non-zero number of constituents
      //////////////////////////////////////////////////////////////////////////////////////
      //float dR = total0_truth.DeltaR(total_truth);
      //if (dR > 0.05){
      //cout << "jet" << event << " does not match consituents and the sum" << endl;
      //	continue;
      //} // tested. All jets have constituents matching the sum 
      //////////////////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////////////////
      //dRmatch < 0.75 between Truth_T and jet
      //
      //
      //float dR = Truth_T.DeltaR(total_truth);
      //	cout << "jet" << event << " with dR = " << dR << endl;
      //if (dR > 0.75){
      //	cout << "jet" << event << " does not match Truth T, with dR = " << dR << endl;
      //	continue;
      //}
      //bool dRmatch = true;
      //////////////////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////////////////
      //dRmatch between Truth_Tdaughters and Truth_T
      //
      //for (unsigned int d = 0; d < Truth_Tdaughter_pt1->size(); d++){
      //truth.SetPtEtaPhiM(Truth_Tdaughter_pt1->at(d),Truth_Tdaughter_eta1->at(d),Truth_Tdaughter_phi1->at(d),Truth_Tdaughter_m1->at(d));
      // //truth = Tdaughter1 + Tdaughter2;
      //float dR = truth.DeltaR(Truth_T);
      // //float dR = truth.DeltaR(total_truth);
      //if (dR > Rfat){
      //cout << "For jet " << event << ", the distance between T daughter " << d << " and Truth_T is " << dR << endl;   
      // //cout << "For jet " << event << ", the distance between the sum of T daughters and Truth_T is " << dR << endl;   
      //dRmatch = false;
      //} //this for loop goes through Tdaughters
      //cout << "Truth Tdaughter " << Truth_Tdaughter_pt1->at(d) << " " << Truth_Tdaughter_eta1->at(d) << " " << Truth_Tdaughter_phi1->at(d) << " " << Truth_Tdaughter_m1->at(d) << endl;
      //}
      //////////////////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////////////////
      //dRmatch between Truth_Wdaughters and Truth_W
      //
      //for (unsigned int d = 0; d < Truth_Wdaughter_pt1->size(); d++){
      //truth.SetPtEtaPhiM(Truth_Wdaughter_pt1->at(d),Truth_Wdaughter_eta1->at(d),Truth_Wdaughter_phi1->at(d),Truth_Wdaughter_m1->at(d));
      // //truth = Wdaughter1 + Wdaughter2;
      //float dR = truth.DeltaR(Truth_W);
      // //float dR = truth.DeltaR(total_truth);
      //if (dR > Rfat){
      //cout << "For jet " << event << ", the distance between W daughter " << d << " and Truth_W is " << dR << endl;   
      //dRmatch = false;
      //}
      //cout << "Truth Wdaughter " << Truth_Wdaughter_pt1->at(d) << " " << Truth_Wdaughter_eta1->at(d) << " " << Truth_Wdaughter_phi1->at(d) << " " << Truth_Wdaughter_m1->at(d) << endl;
      //}
      //////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////
      //if (dRmatch == false){
      //	cout << "jet" << event << " does not dR match for the Truth T/W and T/W daughters" << endl;
      //	continue;
      //}
      //////////////////////////////////////////////////////////////////////////////////////

	fastjet::JetDefinition jet_def1(fastjet::antikt_algorithm, Rfat);
	fastjet::ClusterSequence cs1(input_particles_truth,jet_def1);
	fastjet::ClusterSequence cs1R(input_particles_reco,jet_def1);
	vector<fastjet::PseudoJet> unpruned_jets = sorted_by_pt(cs1.inclusive_jets());
	vector<fastjet::PseudoJet> unpruned_jetsR = sorted_by_pt(cs1R.inclusive_jets());

      //////////////////////////////////////////////////////////////////////////////////////
      //basic pruning on the leading jet
      //////////////////////////////////////////////////////////////////////////////////////
	fastjet::Pruner prunercut(fastjet::cambridge_algorithm, zcut, dcut0);
	fastjet::PseudoJet pruned_jet_cut = prunercut(unpruned_jets[0]);
//	vector<fastjet::PseudoJet> pruned_jet = pruned_jet_cut.constituents();	
	pruned_mjet = pruned_jet_cut.m();
	pruned_ptjet = pruned_jet_cut.perp();
	pruned_etajet = pruned_jet_cut.eta();

      //////////////////////////////////////////////////////////////////////////////////////
      //basic trimming on the leading jet
      //////////////////////////////////////////////////////////////////////////////////////
	fastjet::JetDefinition trimjet_def(fastjet::kt_algorithm, Rfilt0);
	fastjet::Filter trimmercut(trimjet_def, fastjet::SelectorPtFractionMin(fcut0));
//	fastjet::Filter trimmercut(Rfilt0, fastjet::SelectorPtFractionMin(fcut0));
	fastjet::PseudoJet trimmed_jet_cut = trimmercut(unpruned_jets[0]);
	fastjet::PseudoJet trimmed_jet_cutR = trimmercut(unpruned_jetsR[0]);
//	vector<fastjet::PseudoJet> kept = trimmed_jet_cut.pieces();   //this is to bin on the number of subjets		
//	vector<fastjet::PseudoJet> trimmed_jet = trimmed_jet_cut.constituents();
//	vector<fastjet::PseudoJet> trimmed_jet = trimmed_jet_cutR.constituents();
	trimmed_mjet = trimmed_jet_cut.m();
	trimmed_ptjet = trimmed_jet_cut.perp();
	trimmed_etajet = trimmed_jet_cut.eta();
	Truth_jet.SetPtEtaPhiM(trimmed_jet_cut.pt(),trimmed_jet_cut.eta(),trimmed_jet_cut.phi(),trimmed_jet_cut.m());

	fastjet::PseudoJet ungroomed_jet = unpruned_jets[0];
	fastjet::PseudoJet groomed_jet = trimmed_jet_cut;

//	fastjet::PseudoJet ungroomed_jet = unpruned_jetsR[0];
//	fastjet::PseudoJet groomed_jet = trimmed_jet_cutR;



      //////////////////////////////////////////////////////////////////////////////////////
      //used in ATLAS notes to do dR matching
      //////////////////////////////////////////////////////////////////////////////////////
//	double dR = Truth_T.DeltaR(Truth_jet);// uses ROOT to calculate dR between truth jet and truth W/T 
	double dR = truth_w.delta_R(trimmed_jet_cut);// uses fastjet to calculate dR between truth jet and truth w/t
	double dRJ = trimmed_jet_cutR.delta_R(trimmed_jet_cut);// uses fastjet to calculate dR between truth jet and reco jet

	mjet = trimmed_mjet;
	ptjet = trimmed_ptjet;
	etajet = trimmed_etajet;
	mjetR = trimmed_jet_cutR.m();
	ptjetR = trimmed_jet_cutR.perp();
	etajetR = trimmed_jet_cutR.eta();

      if (ptjet > jet_pt_cut_low && ptjet < jet_pt_cut_up && mjet > jet_mass_cut_low && mjet < jet_mass_cut_up && etajet > jet_eta_cut_low && etajet < jet_eta_cut_up && Truth_W_pt1 > 150 && Truth_W_eta1 > -1.2 && Truth_W_eta1 < 1.2 && dR <= 0.75 && ptjetR > jet_pt_cut_low && ptjetR < jet_pt_cut_up && etajetR > jet_eta_cut_low && etajetR < jet_eta_cut_up && dRJ <= 0.75 ) {

//        TSub  T2SubOutputTrim = T_2Subjet(groomed_jet, 0.1, 0.6, 20);

        TSub  T1SubOutputTrim = T_1Subjet(groomed_jet, 0.05, 0.6, 12);
        TSub  T2SubOutputTrim = T_2Subjet(groomed_jet, 0.05, 0.6, 12);
        TSub  T3SubOutputTrim = T_3Subjet(groomed_jet, 0.05, 0.6, 12);
//        TSub  T4SubOutputTrim = T_4Subjet(groomed_jet, 0.05, 0.6, 12);
        Tjet_variable_file << endl;
          
 	Tprun_volatility = T_Pruning (ungroomed_jet, 0.1, 2.0, 20);
 	Ttrim_volatility = T_Trimming(ungroomed_jet, 0.0, 0.1, 20);
 	TakTrecl_volatility = T_AkTreclustering(groomed_jet, 0.1 ,0.6, 20);
 	TkTrecl_volatility = T_kTreclustering(groomed_jet, 0.1 ,0.6, 20);
 	Tsubj_volatility = T2SubOutputTrim.volatility;
 	Tsubj_angle = T2SubOutputTrim.min_angle;
	Ttau2_volatility = T_Nsubjettiness(2, groomed_jet, 1.0, 3.0, 20);

	Tprun_volatility_sigPlot->Fill(Tprun_volatility,weight);
	Ttrim_volatility_sigPlot->Fill(Ttrim_volatility,weight);
	TakTrecl_volatility_sigPlot->Fill(TakTrecl_volatility,weight);
	TkTrecl_volatility_sigPlot->Fill(TkTrecl_volatility,weight);
	Tsubj_volatility_sigPlot->Fill(Tsubj_volatility,weight);
    Tsubj_angle_sigPlot->Fill(Tsubj_angle,weight);
	Ttau2_volatility_sigPlot->Fill(Ttau2_volatility,weight);
	tau21_sigPlot->Fill(GetTau21(groomed_jet),weight);
        
      }// cut selections
}// loop over events
/*
  TFile *Tprun_plot = new TFile("Plots/W_Tprun_volatility_AKt10_300_500_trim.root", "UPDATE");
  Tprun_volatility_sigPlot->Write();
  Tprun_plot->Close();

  TFile *Ttrim_plot = new TFile("Plots/W_Ttrim_volatility_AKt10_300_500_trim.root", "UPDATE");
  Ttrim_volatility_sigPlot->Write();
  Ttrim_plot->Close();

  TFile *TakTrecl_plot = new TFile("Plots/W_TakTrecl_volatility_AKt10_300_500_trim.root", "UPDATE");
  TakTrecl_volatility_sigPlot->Write();
  TakTrecl_plot->Close();

  TFile *TkTrecl_plot = new TFile("Plots/W_TkTrecl_volatility_AKt10_300_500_trim.root", "UPDATE");
  TkTrecl_volatility_sigPlot->Write();
  TkTrecl_plot->Close();

  TFile *Tsubj_plot1 = new TFile("Plots/W_Tsubj_volatility_AKt10_300_500_trim.root", "UPDATE");
  Tsubj_volatility_sigPlot->Write();
  Tsubj_plot1->Close();

  TFile *Tsubj_plot2 = new TFile("Plots/W_Tsubj_angle_beta1_AKt10_300_500_trim.root", "UPDATE");
  Tsubj_angle_sigPlot->Write();
  Tsubj_plot2->Close();

  TFile *Ttau2_plot = new TFile("Plots/W_Ttau2_volatility_AKt10_300_500_trim.root", "UPDATE");
  Ttau2_volatility_sigPlot->Write();
  Ttau2_plot->Close();

  TFile *tau21_plot = new TFile("Plots/W_tau21_AKt10_300_500_trim.root", "UPDATE");
  tau21_sigPlot->Write();
  tau21_plot->Close();


  TFile *Tprun_plot = new TFile("Plots/W_Tprun_volatility_AKt10_800_1000_trim.root", "UPDATE");
  Tprun_volatility_sigPlot->Write();
  Tprun_plot->Close();

  TFile *Ttrim_plot = new TFile("Plots/W_Ttrim_volatility_AKt10_800_1000_trim.root", "UPDATE");
  Ttrim_volatility_sigPlot->Write();
  Ttrim_plot->Close();

  TFile *TakTrecl_plot = new TFile("Plots/W_TakTrecl_volatility_AKt10_800_1000_trim.root", "UPDATE");
  TakTrecl_volatility_sigPlot->Write();
  TakTrecl_plot->Close();

  TFile *TkTrecl_plot = new TFile("Plots/W_TkTrecl_volatility_AKt10_800_1000_trim.root", "UPDATE");
  TkTrecl_volatility_sigPlot->Write();
  TkTrecl_plot->Close();

  TFile *Tsubj_plot1 = new TFile("Plots/W_Tsubj_volatility_AKt10_800_1000_trim.root", "UPDATE");
  Tsubj_volatility_sigPlot->Write();
  Tsubj_plot1->Close();

  TFile *Tsubj_plot2 = new TFile("Plots/W_Tsubj_angle_beta1_AKt10_800_1000_trim.root", "UPDATE");
  Tsubj_angle_sigPlot->Write();
  Tsubj_plot2->Close();

  TFile *Ttau2_plot = new TFile("Plots/W_Ttau2_volatility_AKt10_800_1000_trim.root", "UPDATE");
  Ttau2_volatility_sigPlot->Write();
  Ttau2_plot->Close();

  TFile *tau21_plot = new TFile("Plots/W_tau21_AKt10_800_1000_trim.root", "UPDATE");
  tau21_sigPlot->Write();
  tau21_plot->Close();
*/
/*
  //////////////////////////////////////////////
  //background
  //////////////////////////////////////////////

for (int event = 0; event < bkg_t->GetEntries(); event++) {
    bkg_t->GetEntry(event);

	double weight = bkgEventWeight * bkgCrossSection;// * bkgPileupWeight;

	total_truth.SetPxPyPzE(0,0,0,0);//to add up truth particles
	total_reco.SetPxPyPzE(0,0,0,0);//to add up reco particles
	total0_truth.SetPtEtaPhiM(bkg_pt,bkg_eta,bkg_phi,bkg_m);

	vector<fastjet::PseudoJet> input_particles_truth;//to push back truth particles
	vector<fastjet::PseudoJet> input_particles_reco;//to push back reco particles

    for (int c = 0; c < bkgTruth_n; c++) {
	particle.SetPtEtaPhiM(bkgTruth_pt->at(c),bkgTruth_eta->at(c),bkgTruth_phi->at(c),bkgTruth_m->at(c));
	total_truth = total_truth + particle;
	fastjet::PseudoJet part(particle.Px(),particle.Py(),particle.Pz(),particle.E());
	input_particles_truth.push_back(part);
	input_particles_truth.back().set_user_index(id);
//	if (event == bkg_event_select){
//		EventDisplay->Fill(particle.Eta(),particle.Phi(),particle.E());
//	}
    }// read in .root file and push back the truth particles

    for (int c = 0; c < bkgReco_n; c++) {
	particle.SetPtEtaPhiM(bkgReco_pt->at(c),bkgReco_eta->at(c),bkgReco_phi->at(c),bkgReco_m->at(c));
	total_reco = total_reco + particle;
	fastjet::PseudoJet part(particle.Px(),particle.Py(),particle.Pz(),particle.E());
	input_particles_reco.push_back(part);
	input_particles_reco.back().set_user_index(id);
//	if (event == bkg_event_select){
//		EventDisplay->Fill(particle.Eta(),particle.Phi(),particle.E());
//	}
    }// read in .root file and push back the reco particles

        if(input_particles_truth.size() == 0)	continue;

	fastjet::JetDefinition jet_def1(fastjet::antikt_algorithm, Rfat);
	fastjet::ClusterSequence cs1(input_particles_truth,jet_def1);
	fastjet::ClusterSequence cs1R(input_particles_reco,jet_def1);
	vector<fastjet::PseudoJet> unpruned_jets = sorted_by_pt(cs1.inclusive_jets());
	vector<fastjet::PseudoJet> unpruned_jetsR = sorted_by_pt(cs1R.inclusive_jets());

      //////////////////////////////////////////////////////////////////////////////////////
      //basic pruning on the leading jet
      //////////////////////////////////////////////////////////////////////////////////////
	fastjet::Pruner prunercut(fastjet::cambridge_algorithm, zcut, dcut0);
	fastjet::PseudoJet pruned_jet_cut = prunercut(unpruned_jets[0]);
//	vector<fastjet::PseudoJet> pruned_jet = pruned_jet_cut.constituents();	
	pruned_mjet = pruned_jet_cut.m();
	pruned_ptjet = pruned_jet_cut.perp();
	pruned_etajet = pruned_jet_cut.eta();

      //////////////////////////////////////////////////////////////////////////////////////
      //basic trimming on the leading jet
      //////////////////////////////////////////////////////////////////////////////////////
	fastjet::JetDefinition trimjet_def(fastjet::kt_algorithm, Rfilt0);
	fastjet::Filter trimmercut(trimjet_def, fastjet::SelectorPtFractionMin(fcut0));
//	fastjet::Filter trimmercut(Rfilt0, fastjet::SelectorPtFractionMin(fcut0));
	fastjet::PseudoJet trimmed_jet_cut = trimmercut(unpruned_jets[0]);
	fastjet::PseudoJet trimmed_jet_cutR = trimmercut(unpruned_jetsR[0]);
//	vector<fastjet::PseudoJet> kept = trimmed_jet_cut.pieces();   //this is to bin on the number of subjets		
//	vector<fastjet::PseudoJet> trimmed_jet = trimmed_jet_cut.constituents();
//	vector<fastjet::PseudoJet> trimmed_jet = trimmed_jet_cutR.constituents();
	trimmed_mjet = trimmed_jet_cut.m();
	trimmed_ptjet = trimmed_jet_cut.perp();
	trimmed_etajet = trimmed_jet_cut.eta();
	Truth_jet.SetPtEtaPhiM(trimmed_jet_cut.pt(),trimmed_jet_cut.eta(),trimmed_jet_cut.phi(),trimmed_jet_cut.m());

	fastjet::PseudoJet ungroomed_jet = unpruned_jets[0];
	fastjet::PseudoJet groomed_jet = trimmed_jet_cut;

//	fastjet::PseudoJet ungroomed_jet = unpruned_jetsR[0];
//	fastjet::PseudoJet groomed_jet = trimmed_jet_cutR;

      //////////////////////////////////////////////////////////////////////////////////////
      //used in ATLAS notes to do dR matching
      //////////////////////////////////////////////////////////////////////////////////////
	double dRJ = trimmed_jet_cutR.delta_R(trimmed_jet_cut);// uses fastjet to calculate dR between truth jet and reco jet

	mjet = trimmed_mjet;
	ptjet = trimmed_ptjet;
	etajet = trimmed_etajet;
	mjetR = trimmed_jet_cutR.m();
	ptjetR = trimmed_jet_cutR.perp();
	etajetR = trimmed_jet_cutR.eta();

      if (ptjet > jet_pt_cut_low && ptjet < jet_pt_cut_up && mjet > jet_mass_cut_low && mjet < jet_mass_cut_up && etajet > jet_eta_cut_low && etajet < jet_eta_cut_up && ptjetR > jet_pt_cut_low && ptjetR < jet_pt_cut_up && etajetR > jet_eta_cut_low && etajetR < jet_eta_cut_up && dRJ <= 0.75 ) {

//        TSub  T2SubOutputTrim = T_2Subjet(groomed_jet, 0.05, 0.6, 20);

          TSub  T1SubOutputTrim = T_1Subjet(groomed_jet, 0.05, 0.6, 12);
          TSub  T2SubOutputTrim = T_2Subjet(groomed_jet, 0.05, 0.6, 12);
          TSub  T3SubOutputTrim = T_3Subjet(groomed_jet, 0.05, 0.6, 12);
//          TSub  T4SubOutputTrim = T_4Subjet(groomed_jet, 0.05, 0.6, 12);
          Tjet_variable_file << endl;
          
 	Tprun_volatility = T_Pruning (ungroomed_jet, 0.1, 2.0, 20);
 	Ttrim_volatility = T_Trimming(ungroomed_jet, 0.0, 0.1, 20);
 	TakTrecl_volatility = T_AkTreclustering(groomed_jet, 0.1, 0.6, 20);
 	TkTrecl_volatility = T_kTreclustering(groomed_jet, 0.1, 0.6, 20);
 	Tsubj_volatility = T2SubOutputTrim.volatility;
 	Tsubj_angle = T2SubOutputTrim.min_angle;
	Ttau2_volatility = T_Nsubjettiness(2, groomed_jet, 1.0, 3.0, 20);

	Tprun_volatility_bkgPlot->Fill(Tprun_volatility,weight);
	Ttrim_volatility_bkgPlot->Fill(Ttrim_volatility,weight);
	TakTrecl_volatility_bkgPlot->Fill(TakTrecl_volatility,weight);
	TkTrecl_volatility_bkgPlot->Fill(TkTrecl_volatility,weight);
	Tsubj_volatility_bkgPlot->Fill(Tsubj_volatility,weight);
        Tsubj_angle_bkgPlot->Fill(Tsubj_angle,weight);
	Ttau2_volatility_bkgPlot->Fill(Ttau2_volatility,weight);
	tau21_bkgPlot->Fill(GetTau21(groomed_jet),weight);

      }// cut selections
}// loop over events

  TFile *Tprun_bkg_plot = new TFile("Plots/W_Tprun_volatility_AKt10_300_500_trim_bkg.root", "UPDATE");
  Tprun_volatility_bkgPlot->Write();
  Tprun_bkg_plot->Close();

  TFile *Ttrim_bkg_plot = new TFile("Plots/W_Ttrim_volatility_AKt10_300_500_trim_bkg.root", "UPDATE");
  Ttrim_volatility_bkgPlot->Write();
  Ttrim_bkg_plot->Close();

  TFile *TakTrecl_bkg_plot = new TFile("Plots/W_TakTrecl_volatility_AKt10_300_500_trim_bkg.root", "UPDATE");
  TakTrecl_volatility_bkgPlot->Write();
  TakTrecl_bkg_plot->Close();

  TFile *TkTrecl_bkg_plot = new TFile("Plots/W_TkTrecl_volatility_AKt10_300_500_trim_bkg.root", "UPDATE");
  TkTrecl_volatility_bkgPlot->Write();
  TkTrecl_bkg_plot->Close();

  TFile *Tsubj_bkg_plot1 = new TFile("Plots/W_Tsubj_volatility_AKt10_300_500_trim_bkg.root", "UPDATE");
  Tsubj_volatility_bkgPlot->Write();
  Tsubj_bkg_plot1->Close();

  TFile *Tsubj_bkg_plot2 = new TFile("Plots/W_Tsubj_angle_beta1_AKt10_300_500_trim_bkg.root", "UPDATE");
  Tsubj_angle_bkgPlot->Write();
  Tsubj_bkg_plot2->Close();

  TFile *Ttau2_bkg_plot = new TFile("Plots/W_Ttau2_volatility_AKt10_300_500_trim_bkg.root", "UPDATE");
  Ttau2_volatility_bkgPlot->Write();
  Ttau2_bkg_plot->Close();

  TFile *tau21_bkg_plot = new TFile("Plots/W_tau21_AKt10_300_500_trim_bkg.root", "UPDATE");
  tau21_bkgPlot->Write();
  tau21_bkg_plot->Close();

  TFile *Tprun_bkg_plot = new TFile("Plots/W_Tprun_volatility_AKt10_800_1000_trim_bkg.root", "UPDATE");
  Tprun_volatility_bkgPlot->Write();
  Tprun_bkg_plot->Close();

  TFile *Ttrim_bkg_plot = new TFile("Plots/W_Ttrim_volatility_AKt10_800_1000_trim_bkg.root", "UPDATE");
  Ttrim_volatility_bkgPlot->Write();
  Ttrim_bkg_plot->Close();

  TFile *TakTrecl_bkg_plot = new TFile("Plots/W_TakTrecl_volatility_AKt10_800_1000_trim_bkg.root", "UPDATE");
  TakTrecl_volatility_bkgPlot->Write();
  TakTrecl_bkg_plot->Close();

  TFile *TkTrecl_bkg_plot = new TFile("Plots/W_TkTrecl_volatility_AKt10_800_1000_trim_bkg.root", "UPDATE");
  TkTrecl_volatility_bkgPlot->Write();
  TkTrecl_bkg_plot->Close();

  TFile *Tsubj_bkg_plot1 = new TFile("Plots/W_Tsubj_volatility_AKt10_800_1000_trim_bkg.root", "UPDATE");
  Tsubj_volatility_bkgPlot->Write();
  Tsubj_bkg_plot1->Close();

  TFile *Tsubj_bkg_plot2 = new TFile("Plots/W_Tsubj_angle_beta1_AKt10_800_1000_trim_bkg.root", "UPDATE");
  Tsubj_angle_bkgPlot->Write();
  Tsubj_bkg_plot2->Close();

  TFile *Ttau2_bkg_plot = new TFile("Plots/W_Ttau2_volatility_AKt10_800_1000_trim_bkg.root", "UPDATE");
  Ttau2_volatility_bkgPlot->Write();
  Ttau2_bkg_plot->Close();

  TFile *tau21_bkg_plot = new TFile("Plots/W_tau21_AKt10_800_1000_trim_bkg.root", "UPDATE");
  tau21_bkgPlot->Write();
  tau21_bkg_plot->Close();
*/


  //////////////////////////////////////////////
  //random number generator for pileup
  //////////////////////////////////////////////
  TRandom3 *rand_pileup = new TRandom3();


  //////////////////////////////////////////////
  //main event loop
  //////////////////////////////////////////////
  nEvents = treein->GetEntries();
  cout<<"Number of events: "<<nEvents<<endl;

  for (Long64_t jentry=0; jentry<nEvents; jentry++) {

    if(jentry%10==0)
      cout<<"NTupler: ProcessType="<<ProcessType<<"  entry="<<jentry<<endl;

    //Get next event from input ntuple
    filein->cd();
    treein->GetEntry(jentry);

    /////////////////////////////
    //Reset branches for next event
    /////////////////////////////
    ResetBranches();

    ///////////////////////////////////////////////////
    //read in all final state particles for jet building from pythia input
    ///////////////////////////////////////////////////
    vector<PseudoJet> input_particles;
    input_particles.clear();

    int n_fspart = fspart_id->size();
    for(int i_fspart=0; i_fspart<n_fspart; i_fspart++){

      if(debug){
        cout<<fspart_id->at(i_fspart)<<"  "
            <<fspart_pt->at(i_fspart)<<"  "
            <<fspart_eta->at(i_fspart)<<"  "
            <<fspart_phi->at(i_fspart)<<"  "
            <<fspart_m->at(i_fspart)<<"  "<<endl;
      }

      TLorentzVector temp_p4;
      temp_p4.SetPtEtaPhiM(fspart_pt ->at(i_fspart),
                           fspart_eta->at(i_fspart),
                           fspart_phi->at(i_fspart),
                           fspart_m  ->at(i_fspart));


      input_particles.push_back(PseudoJet(temp_p4.Px(),temp_p4.Py(),temp_p4.Pz(),temp_p4.E()));
    }

/*
    //////////////////////////////////////////////
    //make new input particles collection with pileup
    //////////////////////////////////////////////
    //this will be using min bias events from simulations

    vector<PseudoJet> input_particles_Pileup;
    input_particles_Pileup.clear();
    for(int ipart=0; ipart<n_fspart; ipart++){
      input_particles_Pileup.push_back(input_particles.at(ipart));
    }

    int n_pileup_vertices      = (int)rand_pileup->Poisson(10);
    int n_particles_per_vertex = 5;
    int n_pileup_particles = n_pileup_vertices*n_particles_per_vertex;

    NumberOfVertices = n_pileup_vertices;

    if(debug) cout<<"Pileup: "<<NumberOfVertices<<"  "<<n_particles_per_vertex<<"  "<<n_pileup_particles<<endl;

    for(int ipart=0; ipart<n_pileup_particles; ipart++){

      double m  = 0.0;
      double px = rand_pileup->Gaus(0,5.0);
      double py = rand_pileup->Gaus(0,5.0);
      double pz = rand_pileup->Gaus(0,5.0);
      double E  = pow( m*m + px*px + py*py + pz*pz , 0.5);

      if(debug) cout<<"Pileup: "<<ipart<<"  "<<px<<"  "<<py<<"  "<<pz<<"  "<<E<<endl;

      input_particles_Pileup.push_back(PseudoJet(px,py,pz,E));

    }


    //////////////////////////////////////////////
    //make pseudocalorimeter cells
    //////////////////////////////////////////////
    vector<PseudoJet> calo_cells        = ToyCalorimeter(input_particles);
    vector<PseudoJet> calo_cells_Pileup = ToyCalorimeter(input_particles_Pileup);
*/

    //////////////////////////////////////////////
    // get the resulting jets ordered in pt
    //////////////////////////////////////////////
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.0);

    fastjet::ClusterSequence clust_seq_TruthRaw(input_particles, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_TruthRaw = sorted_by_pt(clust_seq_TruthRaw.inclusive_jets(5.0));
/*
    fastjet::ClusterSequence clust_seq_TruthPileup(input_particles_Pileup, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_TruthPileup = sorted_by_pt(clust_seq_TruthPileup.inclusive_jets(5.0));

    fastjet::ClusterSequence clust_seq_RecoRaw(calo_cells, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_RecoRaw = sorted_by_pt(clust_seq_RecoRaw.inclusive_jets(5.0));

    fastjet::ClusterSequence clust_seq_RecoPileup(calo_cells_Pileup, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_RecoPileup = sorted_by_pt(clust_seq_RecoPileup.inclusive_jets(5.0));
*/


    if(debug){
      // label the columns
      cout<<"jet#  pt  eta  phi  mass"<<endl;
      cout<<"Inclusive"<<endl;
      // print out the details for each jet
      for (unsigned int i = 0; i < inclusive_jets_TruthRaw.size(); i++) {
        cout<<i<<"  "<<inclusive_jets_TruthRaw[i].pt()
               <<"  "<<inclusive_jets_TruthRaw[i].eta()
               <<"  "<<inclusive_jets_TruthRaw[i].phi()
               <<"  "<<inclusive_jets_TruthRaw[i].m()<<endl;
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
    for(int i=0; i < N_r; i++){
      r_values.push_back( r_min+i*(r_max-r_min)/(N_r-1) );
    }
    TelescopingJets T_Mass(axes_def,r_values);

    //Energy correlation functions
    fastjet::contrib::EnergyCorrelatorC2 ecfC2(1.);
    fastjet::contrib::EnergyCorrelatorD2 ecfD2(1.);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecfC3(2, 1.);

    // Filtering with a pt cut as for trimming (arXiv:0912.1342)
    double Rfilt0 = 0.3;
    double fcut0 = 0.05;
    Transformer *trimmer = new Filter(JetDefinition(kt_algorithm, Rfilt0), SelectorPtFractionMin(fcut0) );
    const Transformer &f = *trimmer;

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
    //TruthRaw
    /////////////////////////////
    if(debug) cout<<"TruthRaw jet"<<endl;
    for(int ijet=0; ijet<inclusive_jets_TruthRaw.size(); ijet++){
      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(inclusive_jets_TruthRaw.at(ijet).pt(),
                           inclusive_jets_TruthRaw.at(ijet).eta(),
                           inclusive_jets_TruthRaw.at(ijet).phi(),
                           inclusive_jets_TruthRaw.at(ijet).m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      jetflavor = GetJetTruthFlavor(jettemp, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
      if(debug) cout<<"FillingJet Raw   : flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau21          = GetTau21(inclusive_jets_TruthRaw[ijet]);
      tempJet_Tau32          = GetTau32(inclusive_jets_TruthRaw[ijet]);
      tempJet_D2             = ecfD2(inclusive_jets_TruthRaw[ijet]);

      TSub  T1SubOutputTrim = T_1Subjet(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_T1jet_angle    = T1SubOutput.min_angle;
      tempJet_T1jet          = T1SubOutput.volatility;
        
      TSub  T2SubOutput     = T_2Subjet(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_T2jet_angle    = T2SubOutput.min_angle;
      tempJet_T2jet          = T2SubOutput.volatility;

      TSub  T3SubOutput     = T_3Subjet(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_T3jet_angle    = T3SubOutput.min_angle;
      tempJet_T3jet          = T3SubOutput.volatility;
        
      TSub  T4SubOutput     = T_4Subjet(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_T4jet_angle    = T4SubOutput.min_angle;
      tempJet_T4jet          = T4SubOutput.volatility;

      TSub  T5SubOutput     = T_5Subjet(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_T5jet_angle    = T5SubOutput.min_angle;
      tempJet_T5jet          = T5SubOutput.volatility;

      tempJet_Tpruning       = T_Pruning (inclusive_jets_TruthRaw[ijet], 0.1, 2.0, 20);
      tempJet_Ttrimming      = T_Trimming(inclusive_jets_TruthRaw[ijet], 0.0, 0.1, 20);
      tempJet_Taktreclustering = T_AkTreclustering(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_Tktreclustering = T_kTreclustering(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_TJet_m1        = T_Mass(1,inclusive_jets_TruthRaw[ijet]);
      tempJet_TJet_m2        = T_Mass(2,inclusive_jets_TruthRaw[ijet]);


      if(tempJet_flavor==-1)
        continue;

      TruthRaw_flavor     .push_back(tempJet_flavor);
      TruthRaw_pt         .push_back(tempJet_pt);
      TruthRaw_eta        .push_back(tempJet_eta);
      TruthRaw_phi        .push_back(tempJet_phi);
      TruthRaw_m          .push_back(tempJet_m);
      TruthRaw_Tau21      .push_back(tempJet_Tau21);
      TruthRaw_Tau32      .push_back(tempJet_Tau32);
      TruthRaw_D2         .push_back(tempJet_D2);
      TruthRaw_T1jet_angle.push_back(tempJet_T1jet_angle);
      TruthRaw_T1jet      .push_back(tempJet_T1jet);
      TruthRaw_T2jet_angle.push_back(tempJet_T2jet_angle);
      TruthRaw_T2jet      .push_back(tempJet_T2jet);
      TruthRaw_T3jet_angle.push_back(tempJet_T3jet_angle);
      TruthRaw_T3jet      .push_back(tempJet_T3jet);
      TruthRaw_T4jet_angle.push_back(tempJet_T4jet_angle);
      TruthRaw_T4jet      .push_back(tempJet_T4jet);
      TruthRaw_T5jet_angle.push_back(tempJet_T5jet_angle);
      TruthRaw_T5jet      .push_back(tempJet_T5jet);
      TruthRaw_Tpruning   .push_back(tempJet_Tpruning);
      TruthRaw_Ttrimming  .push_back(tempJet_Ttrimming);
      TruthRaw_Taktreclustering.push_back(tempJet_Taktreclustering);
      TruthRaw_Tktreclustering.push_back(tempJet_Tktreclustering);
      TruthRaw_TJet_m1    .push_back(tempJet_TJet_m1);
      TruthRaw_TJet_m2    .push_back(tempJet_TJet_m2);
    }


    /////////////////////////////
    //TruthRawTrim
    /////////////////////////////
    for (unsigned int ijet = 0; ijet < inclusive_jets_TruthRaw.size(); ijet++) {
      PseudoJet groomed_jet = f(inclusive_jets_TruthRaw[ijet]);

      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(groomed_jet.pt(),
                           groomed_jet.eta(),
                           groomed_jet.phi(),
                           groomed_jet.m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      jetflavor = GetJetTruthFlavor(jettemp, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
      if(debug) cout<<"FillingJet Trimmed: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      if(tempJet_flavor==-1)
        continue;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau21          = GetTau21(groomed_jet);
      tempJet_Tau32          = GetTau32(groomed_jet);
      tempJet_D2             = ecfD2(groomed_jet);

      TSub  T1SubOutputTrim = T_1Subjet(groomed_jet, 0.05, 0.6, 20);
      tempJet_T1jet_angle    = T1SubOutputTrim.min_angle;
      tempJet_T1jet          = T1SubOutputTrim.volatility;
        
      TSub  T2SubOutputTrim = T_2Subjet(groomed_jet, 0.05, 0.6, 20);
      tempJet_T2jet_angle    = T2SubOutputTrim.min_angle;
      tempJet_T2jet          = T2SubOutputTrim.volatility;

      TSub  T3SubOutputTrim = T_3Subjet(groomed_jet, 0.05, 0.6, 20);
      tempJet_T3jet_angle    = T3SubOutputTrim.min_angle;
      tempJet_T3jet          = T3SubOutputTrim.volatility;

      TSub  T4SubOutputTrim = T_4Subjet(groomed_jet, 0.05, 0.6, 20);
      tempJet_T4jet_angle    = T4SubOutputTrim.min_angle;
      tempJet_T4jet          = T4SubOutputTrim.volatility;

      TSub  T5SubOutputTrim = T_5Subjet(groomed_jet, 0.05, 0.6, 20);
      tempJet_T5jet_angle    = T5SubOutputTrim.min_angle;
      tempJet_T5jet          = T5SubOutputTrim.volatility;

      tempJet_Tpruning       = T_Pruning (groomed_jet, 0.1, 2.0, 20);
      tempJet_Ttrimming      = T_Trimming(groomed_jet, 0.0, 0.1, 20);
      tempJet_Taktreclustering = T_AkTreclustering(groomed_jet, 0.05, 0.6, 20);
      tempJet_Tktreclustering = T_kTreclustering(groomed_jet, 0.05, 0.6, 20);
      tempJet_TJet_m1        = T_Mass(1,groomed_jet);
      tempJet_TJet_m2        = T_Mass(2,groomed_jet);

      TruthRawTrim_flavor     .push_back(tempJet_flavor);
      TruthRawTrim_pt         .push_back(tempJet_pt);
      TruthRawTrim_eta        .push_back(tempJet_eta);
      TruthRawTrim_phi        .push_back(tempJet_phi);
      TruthRawTrim_m          .push_back(tempJet_m);
      TruthRawTrim_Tau21      .push_back(tempJet_Tau21);
      TruthRawTrim_Tau32      .push_back(tempJet_Tau32);
      TruthRawTrim_D2         .push_back(tempJet_D2);
      TruthRawTrim_T1jet_angle.push_back(tempJet_T1jet_angle);
      TruthRawTrim_T1jet      .push_back(tempJet_T1jet);
      TruthRawTrim_T2jet_angle.push_back(tempJet_T2jet_angle);
      TruthRawTrim_T2jet      .push_back(tempJet_T2jet);
      TruthRawTrim_T3jet_angle.push_back(tempJet_T3jet_angle);
      TruthRawTrim_T3jet      .push_back(tempJet_T3jet);
      TruthRawTrim_T4jet_angle.push_back(tempJet_T4jet_angle);
      TruthRawTrim_T4jet      .push_back(tempJet_T4jet);
      TruthRawTrim_T5jet_angle.push_back(tempJet_T5jet_angle);
      TruthRawTrim_T5jet      .push_back(tempJet_T5jet);
      TruthRawTrim_Tpruning   .push_back(tempJet_Tpruning);
      TruthRawTrim_Ttrimming  .push_back(tempJet_Ttrimming);
      TruthRawTrim_Taktreclustering .push_back(tempJet_Taktreclustering);
      TruthRawTrim_Tktreclustering .push_back(tempJet_Tktreclustering);
      TruthRawTrim_TJet_m1    .push_back(tempJet_TJet_m1);
      TruthRawTrim_TJet_m2    .push_back(tempJet_TJet_m2);
    }


    //////////////////////////////////////
    //Fill event into tree
    //////////////////////////////////////
    if(debug) cout<<"Filling Tree"<<endl;
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
  TruthRaw_T3jet.clear();
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
  TruthRawTrim_T3jet.clear();
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

}


///=========================================
/// Calorimeter Simulation
///=========================================
vector<PseudoJet> ToyCalorimeter(vector<PseudoJet> truth_particles) {
  const double pi = 3.14159265359;
  const double etaLim = 5.0;
  const int nEta = 100;
  const int nPhi = 63;
  double dEta = 2*etaLim/nEta;
  double dPhi = 2*pi/nPhi;

  double tower[nEta][nPhi];
  for (int i = 0; i < nEta; i++)  for (int j = 0; j < nPhi; j++)  tower[i][j] = -0.001;

  vector<fastjet::PseudoJet> cell_particles;
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
double T_Pruning(PseudoJet& input, double dcut_min, double dcut_max, int N_dcut) {
  double zcut = 0.1; // single choice of zcut. can be further telescoped
  double Tmass = 0;
  vector<double> ms; ms.clear();
  for (int i = 0; i < N_dcut; i++){
  double dcut = dcut_min + (dcut_max-dcut_min)*i/(N_dcut-1);
  fastjet::Pruner pruner(fastjet::cambridge_algorithm, zcut, dcut);
  fastjet::PseudoJet pruned_jet = pruner(input);
  Tmass = pruned_jet.m();
  if(Tmass > M0){ms.push_back(Tmass);}
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ms);
}


///=========================================
/// Telescoping Trimming
///=========================================
double T_Trimming(PseudoJet& input, double fcut_min, double fcut_max, int N_fcut) {
  double Rfilt = 0.2; // single choice of Rfilt. can be further telescoped.
  // used Rfilt = 0.1 for higher pT jets and Rfilt = 0.2 for lower pT jets. 
  double Tmass = 0;
  vector<double> ms; ms.clear();
  for (int i = 0; i < N_fcut; i++){
  double fcut = fcut_min + (fcut_max-fcut_min)*i/(N_fcut-1);
  fastjet::Filter trimmer(Rfilt,fastjet::SelectorPtFractionMin(fcut));
  fastjet::PseudoJet trimmed_jet = trimmer(input);
  Tmass = trimmed_jet.m();
  if(Tmass > M0){ms.push_back(Tmass);}
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ms);
}


///=========================================
/// Telescoping Anti-kT reclustering
///=========================================
double T_AkTreclustering(PseudoJet& input, double R_min, double R_max, int N_R) {

  vector<double> ms; 

  double Tmass = 0;
  double deltaR = (R_max - R_min)/(N_R-1);
  double R      = R_min;
  // used R_min = 0.05 for higher pT jets and R_min = 0.1 for lower pT jets.
  for(int i = 0; i < N_R; i++){

//    R = R_min + i*deltaR;
    fastjet::JetDefinition Tjet_def(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence Tcs(input.constituents(), Tjet_def);

    vector<fastjet::PseudoJet> recoTjets = sorted_by_pt(Tcs.inclusive_jets());
    if(recoTjets.size() < 1) {
        std::cout <<"Warning: recluster number of subjet is "<< recoTjets.size() << std::endl;
        continue;
    }
    if(recoTjets.size() == 1){
        Tsubjet1.SetPxPyPzE(recoTjets[0].px(),recoTjets[0].py(),recoTjets[0].pz(),recoTjets[0].e());
        Tmass = Tsubjet1.M();
	if(Tmass > M0){ms.push_back(Tmass);}
    }
    else if(recoTjets.size() >= 2){
        Tsubjet1.SetPxPyPzE(recoTjets[0].px(),recoTjets[0].py(),recoTjets[0].pz(),recoTjets[0].e());
        Tsubjet2.SetPxPyPzE(recoTjets[1].px(),recoTjets[1].py(),recoTjets[1].pz(),recoTjets[1].e());
        Tmass = (Tsubjet1+Tsubjet2).M();
	if(Tmass > M0){ms.push_back(Tmass);}
    }
    R += deltaR;
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ms);
}


///=========================================
/// Telescoping kT reclustering
///=========================================
double T_kTreclustering(PseudoJet& input, double R_min, double R_max, int N_R) {

  vector<double> ms; 

  double Tmass = 0;
  double deltaR = (R_max - R_min)/(N_R-1);
  double R      = R_min;
  for(int i = 0; i < N_R; i++){

//    R = R_min + i*deltaR;
    fastjet::JetDefinition Tjet_def(fastjet::kt_algorithm, R);
    fastjet::ClusterSequence Tcs(input.constituents(), Tjet_def);

    vector<fastjet::PseudoJet> recoTjets = sorted_by_pt(Tcs.inclusive_jets());
    if(recoTjets.size() < 1) {
        std::cout <<"Warning: recluster number of subjet is "<< recoTjets.size() << std::endl;
        continue;
    }
    if(recoTjets.size() == 1){
        Tsubjet1.SetPxPyPzE(recoTjets[0].px(),recoTjets[0].py(),recoTjets[0].pz(),recoTjets[0].e());
        Tmass = Tsubjet1.M();
	if(Tmass > M0){ms.push_back(Tmass);}
    }
    else if(recoTjets.size() >= 2){
        Tsubjet1.SetPxPyPzE(recoTjets[0].px(),recoTjets[0].py(),recoTjets[0].pz(),recoTjets[0].e());
        Tsubjet2.SetPxPyPzE(recoTjets[1].px(),recoTjets[1].py(),recoTjets[1].pz(),recoTjets[1].e());
        Tmass = (Tsubjet1+Tsubjet2).M();
	if(Tmass > M0){ms.push_back(Tmass);}
    }
    R += deltaR;     
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ms);
}


///=========================================
/// Telescoping Subjet
///=========================================

TSub T_1Subjet(PseudoJet& input, double R_min, double R_max, int N_R){
    vector<double> ms;
    double beta = 1.0;
    double Tmass = 0;
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSub(1, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
    double tau1 = nSub.result(input);
    std::vector<fastjet::PseudoJet> tau1axes = nSub.currentAxes();
//    Tjet_variable_file << tau1axes[0].eta() << " " << tau1axes[0].phi_std() << endl;
    tau_axis1.SetPxPyPzE(tau1axes[0].px(),tau1axes[0].py(),tau1axes[0].pz(),tau1axes[0].e());
    double D_min = tau_axis1.DeltaR(tau_axis1);
    
    double d1;
    double deltaR = (R_max - R_min)/(N_R-1);
    double R      = R_min;
    
    for (int i = 0; i < N_R; i++){
        //    R = R_min + i*deltaR;
        Tsubjet1.SetPxPyPzE(0,0,0,0);
        
        for (unsigned int c = 0; c < input.constituents().size(); c++) {
            particle.SetPxPyPzE(input.constituents()[c].px(),input.constituents()[c].py(),input.constituents()[c].pz(),input.constituents()[c].e());
            d1 = particle.DeltaR(tau_axis1);
            if (d1 <= R){
                Tsubjet1 = Tsubjet1 + particle;
            }
        }
//        Tjet_variable_file << R << " " << Tsubjet1.Perp() << " " << Tsubjet1.M() << endl;
        Tmass = Tsubjet1.M();
        if(Tmass > M0){ms.push_back(Tmass);}
        R += deltaR;
    }
    TSub result;
    result.min_angle = D_min;
    result.volatility = getVolatility(ms);
    return result;
}


TSub T_2Subjet(PseudoJet& input, double R_min, double R_max, int N_R){
  vector<double> m12s;
  double beta = 1.0;
  double Tmass = 0;
  fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
  fastjet::contrib::Nsubjettiness nSub(2, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
  double tau2 = nSub.result(input);
  std::vector<fastjet::PseudoJet> tau2axes = nSub.currentAxes();
//  Tjet_variable_file << tau2axes[0].eta() << " " << tau2axes[0].phi_std() << endl;
//  Tjet_variable_file << tau2axes[1].eta() << " " << tau2axes[1].phi_std() << endl;
  tau_axis1.SetPxPyPzE(tau2axes[0].px(),tau2axes[0].py(),tau2axes[0].pz(),tau2axes[0].e());
  tau_axis2.SetPxPyPzE(tau2axes[1].px(),tau2axes[1].py(),tau2axes[1].pz(),tau2axes[1].e());
  double D_min = tau_axis1.DeltaR(tau_axis2);

  double d1,d2;
  double deltaR = (R_max - R_min)/(N_R-1);
  double R      = R_min;
 
  for (int i = 0; i < N_R; i++){
//    R = R_min + i*deltaR;
    Tsubjet1.SetPxPyPzE(0,0,0,0);
    Tsubjet2.SetPxPyPzE(0,0,0,0);

    for (unsigned int c = 0; c < input.constituents().size(); c++) {
      particle.SetPxPyPzE(input.constituents()[c].px(),input.constituents()[c].py(),input.constituents()[c].pz(),input.constituents()[c].e());
      d1 = particle.DeltaR(tau_axis1);
      d2 = particle.DeltaR(tau_axis2);
      if (d1 <= R && d1 < d2 ){
        Tsubjet1 = Tsubjet1 + particle;
      }
      else if (d2 <= R && d2 < d1 ){
        Tsubjet2 = Tsubjet2 + particle;
      }
    }
//    Tjet_variable_file << R << " " << Tsubjet1.Perp() << " " << Tsubjet1.M() << " " << Tsubjet2.Perp() << " " << Tsubjet2.M() << endl;
    Tmass = (Tsubjet1 + Tsubjet2).M();
    if(Tmass > M0){m12s.push_back(Tmass);}
    R += deltaR;
  }
  TSub result;
  result.min_angle = D_min;
  result.volatility = getVolatility(m12s);
  return result;
}



T3Sub T_3Subjet(PseudoJet& input, double R_min, double R_max, int N_R){
  vector<double> m123s; m123s.clear();
  vector<double> m12s; m12s.clear();
  vector<double> m13s; m13s.clear();
  vector<double> m23s; m23s.clear();
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
//  Tjet_variable_file << tau3axes[0].eta() << " " << tau3axes[0].phi_std() << endl;
//  Tjet_variable_file << tau3axes[1].eta() << " " << tau3axes[1].phi_std() << endl;
//  Tjet_variable_file << tau3axes[2].eta() << " " << tau3axes[2].phi_std() << endl;
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
//    Tjet_variable_file << R << " " << Tsubjet1.Perp() << " " << Tsubjet1.M() << " " << Tsubjet2.Perp() << " " << Tsubjet2.M() << " " << Tsubjet3.Perp() << " " << Tsubjet3.M() << endl;
    Tmass = (Tsubjet1 + Tsubjet2 + Tsubjet3).M();
    Tmass_12 = (Tsubjet1 + Tsubjet2).M();
    Tmass_13 = (Tsubjet1 + Tsubjet3).M();
    Tmass_23 = (Tsubjet2 + Tsubjet3).M();
    if(Tmass > M0){m123s.push_back(Tmass);}
    if(Tmass_12 > M0){m12s.push_back(Tmass_12);}
    if(Tmass_13 > M0){m13s.push_back(Tmass_13);}
    if(Tmass_23 > M0){m23s.push_back(Tmass_23);}
      R += deltaR;
  }
    
  T3Sub result;
  result.min_angle = D_min;
  result.mid_angle = D_mid;
  result.max_angle = D_max;

  result.volatility = getVolatility(m123s);
    
    if ( abs (Tmass_12 - mW) < abs (Tmass_13 - mW) && abs (Tmass_12 - mW) < abs (Tmass_23 - mW)){
        result.mass_W = Tmass_12;
        result.volatility_mass_W = getVolatility(m12s);
    }
    else if ( abs (Tmass_13 - mW) < abs (Tmass_12 - mW) && abs (Tmass_13 - mW) < abs (Tmass_23 - mW)){
        result.mass_W = Tmass_13;
        result.volatility_mass_W = getVolatility(m13s);
    }
    else if ( abs (Tmass_23 - mW) < abs (Tmass_12 - mW) && abs (Tmass_23 - mW) < abs (Tmass_13 - mW)){
        result.mass_W = Tmass_23;
        result.volatility_mass_W = getVolatility(m23s);
    }
    
  return result;
}


TSub T_4Subjet(PseudoJet& input, double R_min, double R_max, int N_R){
    vector<double> m1234s; m1234s.clear();
    // m1234 is the invariant mass of the four subjets
    double beta = 1.0;
    double Tmass = 0;
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSub(4, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
    double tau4 = nSub.result(input);
    std::vector<fastjet::PseudoJet> tau4axes = nSub.currentAxes();
//    Tjet_variable_file << tau4axes[0].eta() << " " << tau4axes[0].phi_std() << endl;
//    Tjet_variable_file << tau4axes[1].eta() << " " << tau4axes[1].phi_std() << endl;
//    Tjet_variable_file << tau4axes[2].eta() << " " << tau4axes[2].phi_std() << endl;
//    Tjet_variable_file << tau4axes[3].eta() << " " << tau4axes[3].phi_std() << endl;
    tau_axis1.SetPxPyPzE(tau4axes[0].px(),tau4axes[0].py(),tau4axes[0].pz(),tau4axes[0].e());
    tau_axis2.SetPxPyPzE(tau4axes[1].px(),tau4axes[1].py(),tau4axes[1].pz(),tau4axes[1].e());
    tau_axis3.SetPxPyPzE(tau4axes[2].px(),tau4axes[2].py(),tau4axes[2].pz(),tau4axes[2].e());
    tau_axis4.SetPxPyPzE(tau4axes[3].px(),tau4axes[3].py(),tau4axes[3].pz(),tau4axes[3].e());
    double tau4_D12 = tau_axis1.DeltaR(tau_axis2);
    double tau4_D13 = tau_axis1.DeltaR(tau_axis3);
    double tau4_D14 = tau_axis1.DeltaR(tau_axis4);
    double tau4_D23 = tau_axis2.DeltaR(tau_axis3);
    double tau4_D24 = tau_axis2.DeltaR(tau_axis4);
    double tau4_D34 = tau_axis3.DeltaR(tau_axis4);
    double D_min = tau4_D12;
    if(tau4_D13 < D_min ) D_min = tau4_D13;
    if(tau4_D14 < D_min ) D_min = tau4_D14;
    if(tau4_D23 < D_min ) D_min = tau4_D23;
    if(tau4_D24 < D_min ) D_min = tau4_D24;
    if(tau4_D34 < D_min ) D_min = tau4_D34;

    double d1, d2, d3, d4;
    double deltaR = (R_max - R_min)/(N_R-1);
    double R      = R_min;
    
    for (int i = 0; i < N_R; i++){
        //    R = R_min + i*deltaR;
        Tsubjet1.SetPxPyPzE(0,0,0,0);
        Tsubjet2.SetPxPyPzE(0,0,0,0);
        Tsubjet3.SetPxPyPzE(0,0,0,0);
        Tsubjet4.SetPxPyPzE(0,0,0,0);
        
        for (unsigned int c = 0; c < input.constituents().size(); c++) {
            particle.SetPxPyPzE(input.constituents()[c].px(),input.constituents()[c].py(),input.constituents()[c].pz(),input.constituents()[c].e());
            d1 = particle.DeltaR(tau_axis1);
            d2 = particle.DeltaR(tau_axis2);
            d3 = particle.DeltaR(tau_axis3);
            d4 = particle.DeltaR(tau_axis4);
            if (d1 < R && d1 < d2 && d1 < d3 && d1 < d4){
                Tsubjet1 = Tsubjet1 + particle;
            }
            else if (d2 < R && d2 < d1 && d2 < d3 && d2 < d4){
                Tsubjet2 = Tsubjet2 + particle;
            }
            else if (d3 < R && d3 < d1 && d3 < d2 && d3 < d4){
                Tsubjet3 = Tsubjet3 + particle;
            }
            else if (d4 < R && d4 < d1 && d4 < d2 && d4 < d3){
                Tsubjet4 = Tsubjet4 + particle;
            }
        }
//        Tjet_variable_file << R << " " << Tsubjet1.Perp() << " " << Tsubjet1.M() << " " << Tsubjet2.Perp() << " " << Tsubjet2.M() << " " << Tsubjet3.Perp() << " " << Tsubjet3.M() << " " << Tsubjet4.Perp() << " " << Tsubjet4.M() << endl;
        Tmass = (Tsubjet1 + Tsubjet2 + Tsubjet3 + Tsubjet4).M();
        if(Tmass > M0){m1234s.push_back(Tmass);}
        R += deltaR;
    }
    TSub result;
    result.min_angle = D_min;
    result.volatility = getVolatility(m1234s);
    return result;
}


TSub T_5Subjet(PseudoJet& input, double R_min, double R_max, int N_R){
    vector<double> m12345s; m12345s.clear();
    // m12345 is the invariant mass of the five subjets
    double beta = 1.0;
    double Tmass = 0;
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSub(5, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
    double tau5 = nSub.result(input);
    std::vector<fastjet::PseudoJet> tau5axes = nSub.currentAxes();
//    Tjet_variable_file << tau5axes[0].eta() << " " << tau5axes[0].phi_std() << endl;
//    Tjet_variable_file << tau5axes[1].eta() << " " << tau5axes[1].phi_std() << endl;
//    Tjet_variable_file << tau5axes[2].eta() << " " << tau5axes[2].phi_std() << endl;
//    Tjet_variable_file << tau5axes[3].eta() << " " << tau5axes[3].phi_std() << endl;
//    Tjet_variable_file << tau5axes[4].eta() << " " << tau5axes[4].phi_std() << endl;
    tau_axis1.SetPxPyPzE(tau5axes[0].px(),tau5axes[0].py(),tau5axes[0].pz(),tau5axes[0].e());
    tau_axis2.SetPxPyPzE(tau5axes[1].px(),tau5axes[1].py(),tau5axes[1].pz(),tau5axes[1].e());
    tau_axis3.SetPxPyPzE(tau5axes[2].px(),tau5axes[2].py(),tau5axes[2].pz(),tau5axes[2].e());
    tau_axis4.SetPxPyPzE(tau5axes[3].px(),tau5axes[3].py(),tau5axes[3].pz(),tau5axes[3].e());
    tau_axis5.SetPxPyPzE(tau5axes[4].px(),tau5axes[4].py(),tau5axes[4].pz(),tau5axes[4].e());
    double tau5_D12 = tau_axis1.DeltaR(tau_axis2);
    double tau5_D13 = tau_axis1.DeltaR(tau_axis3);
    double tau5_D14 = tau_axis1.DeltaR(tau_axis4);
    double tau5_D15 = tau_axis1.DeltaR(tau_axis5);
    double tau5_D23 = tau_axis2.DeltaR(tau_axis3);
    double tau5_D24 = tau_axis2.DeltaR(tau_axis4);
    double tau5_D25 = tau_axis2.DeltaR(tau_axis5);
    double tau5_D34 = tau_axis3.DeltaR(tau_axis4);
    double tau5_D35 = tau_axis3.DeltaR(tau_axis5);
    double tau5_D45 = tau_axis4.DeltaR(tau_axis5);
    double D_min = tau5_D12;
    if(tau5_D13 < D_min ) D_min = tau5_D13;
    if(tau5_D14 < D_min ) D_min = tau5_D14;
    if(tau5_D15 < D_min ) D_min = tau5_D15;
    if(tau5_D23 < D_min ) D_min = tau5_D23;
    if(tau5_D24 < D_min ) D_min = tau5_D24;
    if(tau5_D25 < D_min ) D_min = tau5_D25;
    if(tau5_D34 < D_min ) D_min = tau5_D34;
    if(tau5_D35 < D_min ) D_min = tau5_D35;
    if(tau5_D45 < D_min ) D_min = tau5_D45;
    
    double d1, d2, d3, d4, d5;
    double deltaR = (R_max - R_min)/(N_R-1);
    double R      = R_min;
    
    for (int i = 0; i < N_R; i++){
        //    R = R_min + i*deltaR;
        Tsubjet1.SetPxPyPzE(0,0,0,0);
        Tsubjet2.SetPxPyPzE(0,0,0,0);
        Tsubjet3.SetPxPyPzE(0,0,0,0);
        Tsubjet4.SetPxPyPzE(0,0,0,0);
        Tsubjet5.SetPxPyPzE(0,0,0,0);
        
        for (unsigned int c = 0; c < input.constituents().size(); c++) {
            particle.SetPxPyPzE(input.constituents()[c].px(),input.constituents()[c].py(),input.constituents()[c].pz(),input.constituents()[c].e());
            d1 = particle.DeltaR(tau_axis1);
            d2 = particle.DeltaR(tau_axis2);
            d3 = particle.DeltaR(tau_axis3);
            d4 = particle.DeltaR(tau_axis4);
            d5 = particle.DeltaR(tau_axis5);
            if (d1 < R && d1 < d2 && d1 < d3 && d1 < d4 && d1 < d5){
                Tsubjet1 = Tsubjet1 + particle;
            }
            else if (d2 < R && d2 < d1 && d2 < d3 && d2 < d4 && d2 < d5){
                Tsubjet2 = Tsubjet2 + particle;
            }
            else if (d3 < R && d3 < d1 && d3 < d2 && d3 < d4 && d3 < d5){
                Tsubjet3 = Tsubjet3 + particle;
            }
            else if (d4 < R && d4 < d1 && d4 < d2 && d4 < d3 && d4 < d5){
                Tsubjet4 = Tsubjet4 + particle;
            }
            else if (d5 < R && d5 < d1 && d5 < d2 && d5 < d3 && d5 < d4){
                Tsubjet5 = Tsubjet5 + particle;
            }
        }
//        Tjet_variable_file << R << " " << Tsubjet1.Perp() << " " << Tsubjet1.M() << " " << Tsubjet2.Perp() << " " << Tsubjet2.M() << " " << Tsubjet3.Perp() << " " << Tsubjet3.M() << " " << Tsubjet4.Perp() << " " << Tsubjet4.M() << " " << Tsubjet5.Perp() << " " << Tsubjet5.M() << endl;
        Tmass = (Tsubjet1 + Tsubjet2 + Tsubjet3 + Tsubjet4 + Tsubjet5).M();
        if(Tmass > M0){m12345s.push_back(Tmass);}
        R += deltaR;
    }
    TSub result;
    result.min_angle = D_min;
    result.volatility = getVolatility(m12345s);
    return result;
}



///=========================================
/// Telescoping N-subjettiness
///=========================================
double T_Nsubjettiness(int N, PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> taus; taus.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
    fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
//    fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    taus.push_back(nsub(input));
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(taus);
}

double T_NsubjettinessRatio(int N_num, int N_den, PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> taus; taus.clear();
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
  // getVolatility function provided by TelescopingJets
  return getVolatility(taus);
}


///=========================================
/// Telescoping Energy Correlators
///=========================================
double T_EnergyCorrelator_C2(PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::EnergyCorrelatorC2 ecf(beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ecfs);
}

double T_EnergyCorrelator_D2(PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::EnergyCorrelatorD2 ecf(beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ecfs);
}

double T_EnergyCorrelator_C3(PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecf(3, beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ecfs);
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
    cout<<"DeltaR:   "<<endl
        <<"dRMatch:  "<<dR_match<<endl
        <<"q1:       "<<jettemp.DeltaR(truth_q1)<<endl
        <<"q2:       "<<jettemp.DeltaR(truth_q2)<<endl
        <<"W1:       "<<jettemp.DeltaR(truth_W1)<<endl
        <<"W2:       "<<jettemp.DeltaR(truth_W2)<<endl
        <<"H:        "<<jettemp.DeltaR(truth_H)<<endl
        <<"t1:       "<<jettemp.DeltaR(truth_t1)<<endl
        <<"t2:       "<<jettemp.DeltaR(truth_t2)<<endl;
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

  if(debug) cout<<"Found jet truth flavor: "<<jetflavor<<endl;

  return jetflavor;
}


double GetTau21(PseudoJet& input){

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

double GetTau32(PseudoJet& input){

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

