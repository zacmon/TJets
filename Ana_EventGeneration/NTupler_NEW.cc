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
    
    //  Exit if you don't pass a run card.
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

    bool debug = false;
    if (argc >= 5) {
        std::string argdebug = argv[4];
        if (argdebug == "debug") debug=true;
    }

    std::cout << "InputArguments. ProcessType: " << ProcessType << "\tInputFile: " << InputFile << "\tOutputFile: " << OutputFile << "\tDebug: " << debug << std::endl;

    dR_match = 1.0;

    //////////////////////////////////////////////
    //INPUT
    //////////////////////////////////////////////
    //  Get input file and tree.
    filein = new TFile(InputFile.c_str());
    treein = (TTree*)filein->Get("tree");
    if (debug) treein->Print();

    //  Set up branch linking to addresses.
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
    fileout = new TFile(OutputFile.c_str(), "RECREATE");
    treeout = new TTree("JetTree", "JetTree");

    treeout->Branch("NumberOfVertices",    &NumberOfVertices);

    treeout->Branch("TruthRaw_flavor", &TruthRaw_flavor);
    treeout->Branch("TruthRaw_pt", &TruthRaw_pt);
    treeout->Branch("TruthRaw_eta", &TruthRaw_eta);
    treeout->Branch("TruthRaw_phi", &TruthRaw_phi);
    treeout->Branch("TruthRaw_m", &TruthRaw_m);
    treeout->Branch("TruthRaw_Tau21", &TruthRaw_Tau21);
    treeout->Branch("TruthRaw_Tau32", &TruthRaw_Tau32);
    treeout->Branch("TruthRaw_D2", &TruthRaw_D2);
    treeout->Branch("TruthRaw_T1jet_angle", &TruthRaw_T1jet_angle);
    treeout->Branch("TruthRaw_T1jet_Volatility", &TruthRaw_T1jet);
    treeout->Branch("TruthRaw_T2jet_angle", &TruthRaw_T2jet_angle);
    treeout->Branch("TruthRaw_T2jet_Volatility", &TruthRaw_T2jet);
    treeout->Branch("TruthRaw_T3jet_angle", &TruthRaw_T3jet_angle);
    treeout->Branch("TruthRaw_T3jet_angle1", &TruthRaw_T3jet_angle1);
    treeout->Branch("TruthRaw_T3jet_angle2", &TruthRaw_T3jet_angle2);
    treeout->Branch("TruthRaw_T3jet_Volatility", &TruthRaw_T3jet);
    treeout->Branch("TruthRaw_T3jet_WMass", &TruthRaw_T3jet_W);
    treeout->Branch("TruthRaw_T3jet_WVolatility", &TruthRaw_T3jet_mW);
    treeout->Branch("TruthRaw_T4jet_angle", &TruthRaw_T4jet_angle);
    treeout->Branch("TruthRaw_T4jet_Volatility", &TruthRaw_T4jet);
    treeout->Branch("TruthRaw_T5jet_angle", &TruthRaw_T5jet_angle);
    treeout->Branch("TruthRaw_T5jet_Volatility", &TruthRaw_T5jet);
    treeout->Branch("TruthRaw_Tpruning", &TruthRaw_Tpruning);
    treeout->Branch("TruthRaw_Ttrimming", &TruthRaw_Ttrimming);
    treeout->Branch("TruthRaw_Taktreclustering", &TruthRaw_Taktreclustering);
    treeout->Branch("TruthRaw_Tktreclustering", &TruthRaw_Tktreclustering);
    treeout->Branch("TruthRaw_TJet_m1", &TruthRaw_TJet_m1);
    treeout->Branch("TruthRaw_TJet_m2", &TruthRaw_TJet_m2);

    treeout->Branch("TruthRawTrim_flavor", &TruthRawTrim_flavor);
    treeout->Branch("TruthRawTrim_pt", &TruthRawTrim_pt);
    treeout->Branch("TruthRawTrim_eta", &TruthRawTrim_eta);
    treeout->Branch("TruthRawTrim_phi", &TruthRawTrim_phi);
    treeout->Branch("TruthRawTrim_m", &TruthRawTrim_m);
    treeout->Branch("TruthRawTrim_Tau21", &TruthRawTrim_Tau21);
    treeout->Branch("TruthRawTrim_Tau32", &TruthRawTrim_Tau32);
    treeout->Branch("TruthRawTrim_D2", &TruthRawTrim_D2);
    treeout->Branch("TruthRawTrim_T1jet_angle", &TruthRawTrim_T1jet_angle);
    treeout->Branch("TruthRawTrim_T1jet", &TruthRawTrim_T1jet);
    treeout->Branch("TruthRawTrim_T2jet_angle", &TruthRawTrim_T2jet_angle);
    treeout->Branch("TruthRawTrim_T2jet", &TruthRawTrim_T2jet);
    treeout->Branch("TruthRawTrim_T3jet_angle", &TruthRawTrim_T3jet_angle);
    treeout->Branch("TruthRawTrim_T3jet_angle1", &TruthRawTrim_T3jet_angle1);
    treeout->Branch("TruthRawTrim_T3jet_angle2", &TruthRawTrim_T3jet_angle2);
    treeout->Branch("TruthRawTrim_T3jet", &TruthRawTrim_T3jet);
    treeout->Branch("TruthRawTrim_T3jet_W", &TruthRawTrim_T3jet_W);
    treeout->Branch("TruthRawTrim_T3jet_mW", &TruthRawTrim_T3jet_mW);
    treeout->Branch("TruthRawTrim_T4jet_angle", &TruthRawTrim_T4jet_angle);
    treeout->Branch("TruthRawTrim_T4jet", &TruthRawTrim_T4jet);
    treeout->Branch("TruthRawTrim_T5jet_angle", &TruthRawTrim_T5jet_angle);
    treeout->Branch("TruthRawTrim_T5jet", &TruthRawTrim_T5jet);
    treeout->Branch("TruthRawTrim_Tpruning", &TruthRawTrim_Tpruning);
    treeout->Branch("TruthRawTrim_Ttrimming", &TruthRawTrim_Ttrimming);
    treeout->Branch("TruthRawTrim_Taktreclustering", &TruthRawTrim_Taktreclustering);
    treeout->Branch("TruthRawTrim_Tktreclustering", &TruthRawTrim_Tktreclustering);
    treeout->Branch("TruthRawTrim_TJet_m1", &TruthRawTrim_TJet_m1);
    treeout->Branch("TruthRawTrim_TJet_m2", &TruthRawTrim_TJet_m2);

    
    //////////////////////////////////////////////
    //random number generator for pileup
    //////////////////////////////////////////////
    TRandom3 *rand_pileup = new TRandom3();

    //////////////////////////////////////////////
    //main event loop
    //////////////////////////////////////////////
    nEvents = treein->GetEntries();
    //  nEvents = 10000;
    std::cout << "Number of events: " << nEvents << std::endl;

    for (Long64_t jentry = 0; jentry < nEvents; jentry++) {
        if (jentry % 10 == 0) std::cout << "NTupler: ProcessType=" << ProcessType << "  entry=" << jentry << std::endl;

        //  Get next event from input ntuple.
        filein->cd();
        treein->GetEntry(jentry);

        /////////////////////////////
        //Reset branches for next event
        /////////////////////////////
        ResetBranches();

        ///////////////////////////////////////////////////
        //read in all final state particles for jet building from pythia input
        ///////////////////////////////////////////////////
	std::vector<fastjet::PseudoJet> inputParticles;
        inputParticles.clear();

        int numFinalStateParticles = fspart_id->size();
        for (int i = 0; i < numFinalStateParticles; i++) {

            if (debug) {
                std::cout << fspart_id->at(i) << "  "
			  << fspart_pt->at(i) << "  "
			  << fspart_eta->at(i) << "  "
			  << fspart_phi->at(i) << "  "
			  << fspart_m->at(i) << "  " << std::endl;
            }

            TLorentzVector tempP;
            tempP.SetPtEtaPhiM(fspart_pt->at(i),
                                 fspart_eta->at(i),
                                 fspart_phi->at(i),
                                 fspart_m->at(i));


            inputParticles.push_back(fastjet::PseudoJet(tempP.Px(),tempP.Py(),tempP.Pz(),tempP.E()));
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

if(debug) std::cout<< "Pileup: " <<NumberOfVertices<< "  " <<n_particles_per_vertex<< "  " <<n_pileup_particles<<std::endl;

for(int ipart=0; ipart<n_pileup_particles; ipart++){

double m  = 0.0;
double px = rand_pileup->Gaus(0,5.0);
double py = rand_pileup->Gaus(0,5.0);
double pz = rand_pileup->Gaus(0,5.0);
double E  = pow( m*m + px*px + py*py + pz*pz , 0.5);

if(debug) std::cout<< "Pileup: " <<ipart<< "  " <<px<< "  " <<py<< "  " <<pz<< "  " <<E<<std::endl;

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
        fastjet::JetDefinition jetDefinition(fastjet::antikt_algorithm, 1.0);

        fastjet::ClusterSequence clusterSequenceTruthRaw(inputParticles, jetDefinition);
        std::vector<fastjet::PseudoJet> inclusiveJetsTruthRaw = sorted_by_pt(clusterSequenceTruthRaw.inclusive_jets(5.0));
        /*
          fastjet::ClusterSequence clust_seq_TruthPileup(inputParticles_Pileup, jet_def);
          vector<fastjet::PseudoJet> inclusive_jets_TruthPileup = sorted_by_pt(clust_seq_TruthPileup.inclusive_jets(5.0));

          fastjet::ClusterSequence clust_seq_RecoRaw(calo_cells, jet_def);
          vector<fastjet::PseudoJet> inclusive_jets_RecoRaw = sorted_by_pt(clust_seq_RecoRaw.inclusive_jets(5.0));

          fastjet::ClusterSequence clust_seq_RecoPileup(calo_cells_Pileup, jet_def);
          vector<fastjet::PseudoJet> inclusive_jets_RecoPileup = sorted_by_pt(clust_seq_RecoPileup.inclusive_jets(5.0));
        */


        if (debug) {
            // label the columns
            std::cout<< "jet#  pt  eta  phi  mass" <<std::endl;
            std::cout<< "Inclusive" <<std::endl;
            // print out the details for each jet
            for (unsigned int i = 0; i < inclusiveJetsTruthRaw.size(); i++) {
                std::cout << i << "  " << inclusiveJetsTruthRaw[i].pt()
                         << "  " << inclusiveJetsTruthRaw[i].eta()
                         << "  " << inclusiveJetsTruthRaw[i].phi()
                         << "  " << inclusiveJetsTruthRaw[i].m() << std::endl;
            }
        }

        //////////////////////////////////////////////
        //Setup tools for substructure calculation
        //////////////////////////////////////////////

        //Telescoping jets (this looks like the Telescoping reclustering)
        fastjet::contrib::KT_Axes axesDefinition;
        std::vector<double> radiusValues;
        int numRadii = 20;
        double minRadius = 0.1;
        double maxRadius = 0.6;
        for (int i = 0; i < numRadii; i++) {
            radiusValues.push_back(minRadius + i * (maxRadius-minRadius) / (numRadii - 1));
        }
        TelescopingJets T_Mass(axesDefinition,radiusValues);

        //Energy correlation functions
        fastjet::contrib::EnergyCorrelatorC2 ecfC2(1.0);
        fastjet::contrib::EnergyCorrelatorD2 ecfD2(1.0);
        fastjet::contrib::EnergyCorrelatorDoubleRatio ecfC3(2, 1.0);

        // Filtering with a pt cut as for trimming (arXiv:0912.1342)
        double Rfilt0 = 0.3;
        double fcut0 = 0.05;
	fastjet::Transformer *trimmer = new fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, Rfilt0), fastjet::SelectorPtFractionMin(fcut0));
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
        /*
/////////////////////////////
//TruthRaw
/////////////////////////////
if(debug) std::cout<< "TruthRaw jet" <<std::endl;
for(int ijet=0; ijet<inclusiveJetsTruthRaw.size(); ijet++){
TLorentzVector jettemp;
jettemp.SetPtEtaPhiM(inclusiveJetsTruthRaw.at(ijet).pt(),
inclusiveJetsTruthRaw.at(ijet).eta(),
inclusiveJetsTruthRaw.at(ijet).phi(),
inclusiveJetsTruthRaw.at(ijet).m());

/////////////////////////////////
//Getting truth label for filling into ntuple
/////////////////////////////////
jetflavor = GetJetTruthFlavor(jettemp, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
if(debug) std::cout<< "FillingJet Raw   : flav=" <<jetflavor<< "  pt=" <<jettemp.Pt()<< "  m=" <<jettemp.M()<<std::endl;

if(jetflavor==-1)
continue;

/////////////////////////////////
//Fill variables that will go into ntuple
/////////////////////////////////
tempJet_flavor         = jetflavor;
tempJet_pt             = jettemp.Pt();
tempJet_eta            = jettemp.Eta();
tempJet_phi            = jettemp.Phi();
tempJet_m              = jettemp.M();
tempJet_Tau21          = GetTau21(inclusiveJetsTruthRaw[ijet]);
tempJet_Tau32          = GetTau32(inclusiveJetsTruthRaw[ijet]);
//      tempJet_D2             = ecfD2(inclusiveJetsTruthRaw[ijet]);

TSub  T2SubOutput     = T_2Subjet(inclusiveJetsTruthRaw[ijet], 0.05, 0.6, 20);
tempJet_T2jet_angle    = T2SubOutput.min_angle;
tempJet_T2jet          = T2SubOutput.volatility;
        
TSub  T3SubOutput     = T_3Subjet(inclusiveJetsTruthRaw[ijet], 0.05, 0.6, 20);
tempJet_T3jet_angle    = T3SubOutput.min_angle;
tempJet_T3jet          = T3SubOutput.volatility;
        
//      TSub  T4SubOutput     = T_4Subjet(inclusiveJetsTruthRaw[ijet], 0.05, 0.6, 20);
//      tempJet_T4jet_angle    = T4SubOutput.min_angle;
//      tempJet_T4jet          = T4SubOutput.volatility;

//      TSub  T5SubOutput     = T_5Subjet(inclusiveJetsTruthRaw[ijet], 0.05, 0.6, 20);
//      tempJet_T5jet_angle    = T5SubOutput.min_angle;
//      tempJet_T5jet          = T5SubOutput.volatility;

tempJet_Tpruning       = T_Pruning (inclusiveJetsTruthRaw[ijet], 0.1, 2.0, 20);
tempJet_Ttrimming      = T_Trimming(inclusiveJetsTruthRaw[ijet], 0.0, 0.1, 20);
tempJet_Taktreclustering = T_AkTreclustering(inclusiveJetsTruthRaw[ijet], 0.05, 0.6, 20);
//      tempJet_Tktreclustering = T_kTreclustering(inclusiveJetsTruthRaw[ijet], 0.05, 0.6, 20);
//      tempJet_TJet_m1        = T_Mass(1,inclusiveJetsTruthRaw[ijet]);
//      tempJet_TJet_m2        = T_Mass(2,inclusiveJetsTruthRaw[ijet]);

//      if(tempJet_flavor==-1)
//        continue;

TruthRaw_flavor     .push_back(tempJet_flavor);
TruthRaw_pt         .push_back(tempJet_pt);
TruthRaw_eta        .push_back(tempJet_eta);
TruthRaw_phi        .push_back(tempJet_phi);
TruthRaw_m          .push_back(tempJet_m);
TruthRaw_Tau21      .push_back(tempJet_Tau21);
TruthRaw_Tau32      .push_back(tempJet_Tau32);
//      TruthRaw_D2         .push_back(tempJet_D2);
TruthRaw_T2jet_angle.push_back(tempJet_T2jet_angle);
TruthRaw_T2jet      .push_back(tempJet_T2jet);
TruthRaw_T3jet_angle.push_back(tempJet_T3jet_angle);
TruthRaw_T3jet      .push_back(tempJet_T3jet);
//      TruthRaw_T4jet_angle.push_back(tempJet_T4jet_angle);
//      TruthRaw_T4jet      .push_back(tempJet_T4jet);
//      TruthRaw_T5jet_angle.push_back(tempJet_T5jet_angle);
//      TruthRaw_T5jet      .push_back(tempJet_T5jet);
TruthRaw_Tpruning   .push_back(tempJet_Tpruning);
TruthRaw_Ttrimming  .push_back(tempJet_Ttrimming);
TruthRaw_Taktreclustering.push_back(tempJet_Taktreclustering);
//      TruthRaw_Tktreclustering.push_back(tempJet_Tktreclustering);
//      TruthRaw_TJet_m1    .push_back(tempJet_TJet_m1);
//      TruthRaw_TJet_m2    .push_back(tempJet_TJet_m2);
}
        */

        /////////////////////////////
        //TruthRawTrim
        /////////////////////////////
        for (unsigned int iJet = 0; iJet < inclusiveJetsTruthRaw.size(); iJet++) {
            fastjet::PseudoJet groomedJet = f(inclusiveJetsTruthRaw[iJet]);

            TLorentzVector pJet;
            pJet.SetPtEtaPhiM(groomedJet.pt(),
                                 groomedJet.eta(),
                                 groomedJet.phi(),
                                 groomedJet.m());

            /////////////////////////////////
            //Getting truth label for filling into ntuple
            /////////////////////////////////
            jetflavor = GetJetTruthFlavor(pJet, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
            if(debug) std::cout << "FillingJet Trimmed: flav=" << jetflavor << "  pt=" << pJet.Pt() << "  m=" << pJet.M() <<std::endl;

             if(jetflavor == -1) continue;
        
            /////////////////////////////////
            //Fill variables that will go into ntuple
            /////////////////////////////////
            tempJet_flavor         = jetflavor;
            tempJet_pt             = pJet.Pt();
            tempJet_eta            = pJet.Eta();
            tempJet_phi            = pJet.Phi();
            tempJet_m              = pJet.M();
            tempJet_Tau21          = GetTau21(groomedJet);
            tempJet_Tau32          = GetTau32(groomedJet);
            //      tempJet_D2             = ecfD2(groomedJet);

            TSub  T1SubOutputTrim  = TNSubjet(groomedJet, 1, 0.1, 0.6, 20);
            tempJet_T1jet_angle    = T1SubOutputTrim.min_angle;
            tempJet_T1jet          = T1SubOutputTrim.volatility;
        
            TSub  T2SubOutputTrim  = TNSubjet(groomedJet, 2, 0.1, 0.6, 20);
            tempJet_T2jet_angle    = T2SubOutputTrim.min_angle;
            tempJet_T2jet          = T2SubOutputTrim.volatility;

            T3Sub  T3SubOutputTrim = T_3Subjet(groomedJet, 0.1, 0.6, 20);
            tempJet_T3jet_angle    = T3SubOutputTrim.min_angle;
            //      tempJet_T3jet_angle1    = T3SubOutputTrim.mid_angle;
            //      tempJet_T3jet_angle2    = T3SubOutputTrim.max_angle;
            tempJet_T3jet          = T3SubOutputTrim.volatility;
            tempJet_T3jet_W        = T3SubOutputTrim.mass_W;
            tempJet_T3jet_mW       = T3SubOutputTrim.volatility_mass_W;
            //      TSub  T4SubOutputTrim = T_4Subjet(groomedJet, 0.05, 0.6, 20);
            //      tempJet_T4jet_angle    = T4SubOutputTrim.min_angle;
            //      tempJet_T4jet          = T4SubOutputTrim.volatility;

            //      TSub  T5SubOutputTrim = T_5Subjet(groomedJet, 0.05, 0.6, 20);
            //      tempJet_T5jet_angle    = T5SubOutputTrim.min_angle;
            //      tempJet_T5jet          = T5SubOutputTrim.volatility;

            //      tempJet_Tpruning       = T_Pruning (groomedJet, 0.1, 2.0, 20);
            //      tempJet_Ttrimming      = T_Trimming(groomedJet, 0.0, 0.1, 20);
            //      tempJet_Taktreclustering = T_AkTreclustering(groomedJet, 0.05, 0.6, 20);
            //      tempJet_Tktreclustering = T_kTreclustering(groomedJet, 0.05, 0.6, 20);
            //      tempJet_TJet_m1        = T_Mass(1,groomedJet);
            //      tempJet_TJet_m2        = T_Mass(2,groomedJet);

            TruthRawTrim_flavor.push_back(tempJet_flavor);
            TruthRawTrim_pt.push_back(tempJet_pt);
            TruthRawTrim_eta.push_back(tempJet_eta);
            TruthRawTrim_phi.push_back(tempJet_phi);
            TruthRawTrim_m.push_back(tempJet_m);
            TruthRawTrim_Tau21.push_back(tempJet_Tau21);
            TruthRawTrim_Tau32.push_back(tempJet_Tau32);
            //      TruthRawTrim_D2         .push_back(tempJet_D2);
            TruthRawTrim_T1jet_angle.push_back(tempJet_T1jet_angle);
            TruthRawTrim_T1jet.push_back(tempJet_T1jet);
            TruthRawTrim_T2jet_angle.push_back(tempJet_T2jet_angle);
            TruthRawTrim_T2jet.push_back(tempJet_T2jet);
            TruthRawTrim_T3jet_angle.push_back(tempJet_T3jet_angle);
            //      TruthRawTrim_T3jet_angle1.push_back(tempJet_T3jet_angle1);
            //      TruthRawTrim_T3jet_angle2.push_back(tempJet_T3jet_angle2);
            TruthRawTrim_T3jet.push_back(tempJet_T3jet);
            TruthRawTrim_T3jet_W.push_back(tempJet_T3jet_W);
            TruthRawTrim_T3jet_mW.push_back(tempJet_T3jet_mW);
            //      TruthRawTrim_T4jet_angle.push_back(tempJet_T4jet_angle);
            //      TruthRawTrim_T4jet      .push_back(tempJet_T4jet);
            //      TruthRawTrim_T5jet_angle.push_back(tempJet_T5jet_angle);
            //      TruthRawTrim_T5jet      .push_back(tempJet_T5jet);
            //      TruthRawTrim_Tpruning   .push_back(tempJet_Tpruning);
            //      TruthRawTrim_Ttrimming  .push_back(tempJet_Ttrimming);
            //      TruthRawTrim_Taktreclustering .push_back(tempJet_Taktreclustering);
            //      TruthRawTrim_Tktreclustering .push_back(tempJet_Tktreclustering);
            //      TruthRawTrim_TJet_m1    .push_back(tempJet_TJet_m1);
            //      TruthRawTrim_TJet_m2    .push_back(tempJet_TJet_m2);
        }

        if (debug) std::cout<< "Filling Tree" <<std::endl;
        treeout->Fill();
    }

    fileout->cd();
    treeout->Write();
    fileout->Close();

    return 0;
}

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
}

///=========================================
/// Calorimeter Simulation
///=========================================
std::vector<fastjet::PseudoJet> ToyCalorimeter(std::vector<fastjet::PseudoJet> truthParticles) {
    const double pi = 3.14159265359;
    const double etaLim = 5.0;
    const int nEta = 100;
    const int nPhi = 63;
    double dEta = 2 * etaLim / nEta;
    double dPhi = 2 * pi / nPhi;

    double tower[nEta][nPhi];
    for (int i = 0; i < nEta; i++) {
	for (int j = 0; j < nPhi; j++) {
	    tower[i][j] = -0.001;
	}
    }

    std::vector<fastjet::PseudoJet> cellParticles;
    for (unsigned int p = 0; p < truthParticles.size(); p++) {
        fastjet::PseudoJet part = truthParticles.at(p);
	
        int etaCell = int((part.eta() + etaLim) / dEta);
        int phiCell = int(part.phi() / dPhi);
        if (etaCell >= 0 && etaCell < nEta && phiCell >=0 && phiCell < nPhi){
            tower[etaCell][phiCell] += part.e();
        }
    }

    for (int i = 0; i < nEta; i++) {
	for (int j = 0; j < nPhi; j++) {
            if (tower[i][j] > 0) {
                double etaLocal = -etaLim + dEta * (i + 0.5);
                double phiLocal = dPhi * (j + 0.5);
                double thetaLocal = 2 * atan(exp(-etaLocal));
                cellParticles.push_back(tower[i][j] * fastjet::PseudoJet(sin(thetaLocal) * cos(phiLocal),
									 sin(thetaLocal) * sin(phiLocal),
									 cos(thetaLocal), 1));
            }
        }
    }
    return cellParticles;
}

///=========================================
/// Telescoping Pruning
///=========================================
double T_Pruning(fastjet::PseudoJet& input, double minDCut, double maxDCut, int numDCuts) {
    double zCut = 0.1; //  Single choice of zCut, but can be further telescoped.
    double mass = 0;
    std::vector<double> telescopingMass;

    for (int i = 0; i < numDCuts; i++) {
        double dCut = minDCut + (maxDCut - minDCut) * i / (numDCuts - 1);
        fastjet::Pruner pruner(fastjet::cambridge_algorithm, zCut, dCut);
        fastjet::PseudoJet prunedJet = pruner(input);
        mass = prunedJet.m();
        if (mass > M0) telescopingMass.push_back(mass);
    }

    if (telescopingMass.size() > 0) return getVolatility(telescopingMass);
    
    std::cout << "WARNING: zero entries for T_Pruning " << minDCut << " " << maxDCut << " " << numDCuts << " " <<  std::endl;
    return -1;   
}


///=========================================
/// Telescoping Trimming
///=========================================
double T_Trimming(fastjet::PseudoJet& input, double minFCut, double maxFCut, int numFCuts) {
    double Rfilt = 0.1; //  Single choice of Rfilt but can be further telescoped.
    //  Used Rfilt = 0.1 for higher pT jets and Rfilt = 0.2 for lower pT jets.
    double mass = 0;
    std::vector<double> telescopingMass;
    
    for (int i = 0; i < numFCuts; i++) {
        double fCut = minFCut + (maxFCut - minFCut) * i / (numFCuts - 1);
        fastjet::Filter trimmer(Rfilt,fastjet::SelectorPtFractionMin(fCut));
        fastjet::PseudoJet trimmedJet = trimmer(input);
        mass = trimmedJet.m();
        if (mass > M0) telescopingMass.push_back(mass); 
    }
    
    if (telescopingMass.size() > 0) return getVolatility(telescopingMass);
    
    std::cout << "WARNING: zero entries for T_Trimming " << minFCut << " " << maxFCut << " " << numFCuts << " " << std::endl;
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
    if (numSubjets > 0) {
	auto minTauAxisAngleIt = std::min_element(anglesBetweenTauAxes.cbegin(), anglesBetweenTauAxes.cend());
	minTauAxisAngle = *minTauAxisAngleIt;
    }
    
    double deltaR = (maxRadius - minRadius) / (numRadii - 1);

    //  TODO
    //  Exploit that constituents at smaller R are contained in the
    //  larger R subjets. No need to loop and add and calculate all
    //  the necessary items for those constituents. Sort vector by
    //  what metric? Store in another vector?
    for (double r = minRadius; r <= maxRadius; r += deltaR) {
	std::vector<TLorentzVector> TSubjets(numSubjets);
	for (auto const &constituent : input.constituents()) {
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
	    }
	}

	TLorentzVector sumTSubjet;
	for (auto const &TSubjet : TSubjets) {
	    sumTSubjet += TSubjet;
	}
	mass = sumTSubjet.M();

	if (mass > M0) telescopingMasses.push_back(mass);
    }

    TSub result;
    result.min_angle = minTauAxisAngle;
    if (telescopingMasses.size() > 0) result.volatility = getVolatility(telescopingMasses);
    else {
	std::cout << "WARNING: zero entries for TNSubjet: \tNumSubjets: " << numSubjets << "\tminRadius: " << minRadius << "\tmaxRadius: " << maxRadius << "\tnumRadii: " << numRadii <<  std::endl;
        result.volatility = -1;
    }
    return result;
}	

T3Sub T_3Subjet(fastjet::PseudoJet& input, double minRadius, double maxRadius, int numRadii) {
    std::vector<double> m123s; m123s.clear();
    std::vector<double> m12s; m12s.clear();
    std::vector<double> m13s; m13s.clear();
    std::vector<double> m23s; m23s.clear();
    // m123 is the invariant mass of the three subjets
    double beta = 1.0;
    double mass = 0;
    double mass_12 = 0;
    double mass_13 = 0;
    double mass_23 = 0;
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
    fastjet::contrib::Nsubjettiness nSubjettiness(3, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    double tau3 = nSubjettiness.result(input);
    std::vector<fastjet::PseudoJet> tau3axes = nSubjettiness.currentAxes();
    tau_axis1.SetPxPyPzE(tau3axes[0].px(),tau3axes[0].py(),tau3axes[0].pz(),tau3axes[0].e());
    tau_axis2.SetPxPyPzE(tau3axes[1].px(),tau3axes[1].py(),tau3axes[1].pz(),tau3axes[1].e());
    tau_axis3.SetPxPyPzE(tau3axes[2].px(),tau3axes[2].py(),tau3axes[2].pz(),tau3axes[2].e());
    double tau3DistanceBetweenTauAxis12 = tau_axis1.DeltaR(tau_axis2);
    double tau3DistanceBetweenTauAxis13 = tau_axis1.DeltaR(tau_axis3);
    double tau3DistanceBetweenTauAxis23 = tau_axis2.DeltaR(tau_axis3);
    double minTauAxisAngle = tau3DistanceBetweenTauAxis12;
    double D_mid = tau3DistanceBetweenTauAxis13;
    double D_max = tau3DistanceBetweenTauAxis23;
    double D_temp;
    if (D_mid < minTauAxisAngle ) {D_temp = minTauAxisAngle; minTauAxisAngle = D_mid; D_mid = D_temp;}
    if (D_max < minTauAxisAngle ) {D_temp = minTauAxisAngle; minTauAxisAngle = D_max; D_max = D_temp;}
    if (D_max < D_mid ) {D_temp = D_mid; D_mid = D_max; D_max = D_temp;}
    
    double d1, d2, d3;
    double deltaR = (maxRadius - minRadius)/(numRadii-1);
    double R      = minRadius;

    for (double r = minRadius; r <= maxRadius; r += deltaR) {
        //    R = minRadius + i*deltaR;
        Tsubjet1.SetPxPyPzE(0,0,0,0);
        Tsubjet2.SetPxPyPzE(0,0,0,0);
        Tsubjet3.SetPxPyPzE(0,0,0,0);

        for (unsigned int c = 0; c < input.constituents().size(); c++) {
            particle.SetPxPyPzE(input.constituents()[c].px(),input.constituents()[c].py(),input.constituents()[c].pz(),input.constituents()[c].e());
            d1 = particle.DeltaR(tau_axis1);
            d2 = particle.DeltaR(tau_axis2);
            d3 = particle.DeltaR(tau_axis3);
            if (d1 <= R && d1 < d2 && d1 < d3) {
                Tsubjet1 = Tsubjet1 + particle;
            }
            else if (d2 <= R && d2 < d1 && d2 < d3) {
                Tsubjet2 = Tsubjet2 + particle;
            }
            else if (d3 <= R && d3 < d1 && d3 < d2) {
                Tsubjet3 = Tsubjet3 + particle;
            }
        }
        mass = (Tsubjet1 + Tsubjet2 + Tsubjet3).M();
        mass_12 = (Tsubjet1 + Tsubjet2).M();
        mass_13 = (Tsubjet1 + Tsubjet3).M();
        mass_23 = (Tsubjet2 + Tsubjet3).M();
        if (mass > M0) {m123s.push_back(mass);}
        if (mass_12 > M0) {m12s.push_back(mass_12);}
        if (mass_13 > M0) {m13s.push_back(mass_13);}
        if (mass_23 > M0) {m23s.push_back(mass_23);}
        R += deltaR;
    }
    T3Sub result;
    result.min_angle = minTauAxisAngle;
    result.mid_angle = D_mid;
    result.max_angle = D_max;
    
    if (m123s.size() > 0 && m12s.size() > 0 && m13s.size() > 0 && m23s.size() > 0) {
        result.volatility = getVolatility(m123s);
        if ( std::abs (mass_12 - mW) < std::abs (mass_13 - mW) && std::abs (mass_12 - mW) < std::abs (mass_23 - mW)) {
            result.mass_W = mass_12;
            result.volatility_mass_W = getVolatility(m12s);
        }
        else if ( std::abs (mass_13 - mW) < std::abs (mass_12 - mW) && std::abs (mass_13 - mW) < std::abs (mass_23 - mW)) {
            result.mass_W = mass_13;
            result.volatility_mass_W = getVolatility(m13s);
        }
        else if ( std::abs (mass_23 - mW) < std::abs (mass_12 - mW) && std::abs (mass_23 - mW) < std::abs (mass_13 - mW)) {
            result.mass_W = mass_23;
            result.volatility_mass_W = getVolatility(m23s);
        }
    }else{
        std::cout << "WARNING: zero entries for T_3Subjet " << minRadius << " " << maxRadius << " " << numRadii << " " <<  std::endl;
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
        fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);
        fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
        //    fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);
        taus.push_back(nsub(input));
    }
    // getVolatility function provided by TelescopingJets
    if (taus.size()>0)
        return getVolatility(taus);
    
    std::cout << "WARNING: zero entries for T_Nsubjettiness " << beta_min << " " << beta_max << " " << N_beta << " " <<  std::endl;
    return -1;
    
}

double T_NsubjettinessRatio(int N_num, int N_den, fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta) {
    std::vector<double> taus; taus.clear();
    for (int i = 0; i < N_beta; i++) {

        double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);

        fastjet::contrib::UnnormalizedMeasure nSubMeasure(beta);

        fastjet::contrib::Nsubjettiness nsub_num(N_num, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);
        fastjet::contrib::Nsubjettiness nsub_den(N_den, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);

        double num=nsub_num(input);
        double den=nsub_den(input);

        if (den!=0)
            taus.push_back(num/den);
        else
            taus.push_back(-1.0);

    }
    if (taus.size()==0) {
        std::cout << "WARNING: zero entries for T_NsubjetinessRatio " << beta_min << " " << beta_max << " " << N_beta << " " <<  std::endl;
        return -1;
    }
    // getVolatility function provided by TelescopingJets
    if (taus.size()>0)
        return getVolatility(taus);
    
    std::cout << "WARNING: zero entries for T_NsubjettinessRatio " << beta_min << " " << beta_max << " " << N_beta << " " <<  std::endl;
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
    if (ecfs.size()>0)
        return getVolatility(ecfs);
    
    std::cout << "WARNING: zero entries for T_EnergyCorrelator_C2 " << beta_min << " " << beta_max << " " << N_beta << " " <<  std::endl;
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
    if (ecfs.size()>0)
        return getVolatility(ecfs);
    
    std::cout << "WARNING: zero entries for T_EnergyCorrelator_C2 " << beta_min << " " << beta_max << " " << N_beta << " " <<  std::endl;
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
    if (ecfs.size()>0)
        return getVolatility(ecfs);
    
    std::cout << "WARNING: zero entries for T_EnergyCorrelator_C3 " << beta_min << " " << beta_max << " " << N_beta << " " <<  std::endl;
    return -1;
    
}

///========================================
int GetJetTruthFlavor(TLorentzVector pJet,
                      TLorentzVector truth_t1,
                      TLorentzVector truth_t2,
                      TLorentzVector truth_W1,
                      TLorentzVector truth_W2,
                      TLorentzVector truth_H,
                      int debug) {
    if (debug) {
        std::cout<< "DeltaR:   " <<std::endl
                 << "dRMatch:  " <<dR_match<<std::endl
                 << "q1:       " <<pJet.DeltaR(truth_q1)<<std::endl
                 << "q2:       " <<pJet.DeltaR(truth_q2)<<std::endl
                 << "W1:       " <<pJet.DeltaR(truth_W1)<<std::endl
                 << "W2:       " <<pJet.DeltaR(truth_W2)<<std::endl
                 << "H:        " <<pJet.DeltaR(truth_H)<<std::endl
                 << "t1:       " <<pJet.DeltaR(truth_t1)<<std::endl
                 << "t2:       " <<pJet.DeltaR(truth_t2)<<std::endl;
    }
    int jetflavor = -1;
    if     (pJet.DeltaR(truth_q1)<dR_match || pJet.DeltaR(truth_q2)<dR_match) {
        jetflavor = 0;
    }
    else if (pJet.DeltaR(truth_W1)<dR_match || pJet.DeltaR(truth_W2)<dR_match) {
        jetflavor = 1;
    }
    else if (pJet.DeltaR(truth_t1)<dR_match || pJet.DeltaR(truth_t2)<dR_match) {
        jetflavor = 3;
    }
    else if (pJet.DeltaR(truth_H)<dR_match) {
        jetflavor = 3;
    }
    else{
        jetflavor = -1;
    }

    if (debug) std::cout<< "Found jet truth flavor: " <<jetflavor<<std::endl;

    return jetflavor;
}


double GetTau21(fastjet::PseudoJet& input) {

    float tau21=-1;

    //N-subjettiness
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(1.);
    fastjet::contrib::Nsubjettiness nsub1(1, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    //  fastjet::contrib::Nsubjettiness nsub1(1, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);
    //  fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);

    float tau1 = nsub1(input);
    float tau2 = nsub2(input);

    if (tau1>0)
        tau21 = tau2/tau1;

    return tau21;

}

double GetTau32(fastjet::PseudoJet& input) {

    float tau32=-1;

    //N-subjettiness
    fastjet::contrib::UnnormalizedMeasure nSubMeasure(1.);
    fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    fastjet::contrib::Nsubjettiness nsub3(3, fastjet::contrib::OnePass_KT_Axes(), nSubMeasure);
    //  fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);
    //  fastjet::contrib::Nsubjettiness nsub3(3, fastjet::contrib::WTA_KT_Axes(), nSubMeasure);

    float tau2 = nsub2(input);
    float tau3 = nsub3(input);

    if (tau2>0)
        tau32 = tau3/tau2;

    return tau32;

}

