#ifndef NTUPLER_H
#define NTUPLER_H

#define CERRD cout<<"Problem on "<<__FILE__<<"  "<<__LINE__<<endl;

#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TChain.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstdio>

#include "TelescopingJets.hh"

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>

#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

using namespace fastjet;
using namespace std;


///////////////////////////
//input file and tree
///////////////////////////
TTree *treein;
TFile *filein;

int nEvents;

// delta R for truth labelling
double dR_match;

///////////////////////////
//input tree branches
///////////////////////////
double truth_q1_pt;
double truth_q1_eta;
double truth_q1_phi;
double truth_q1_m;

double truth_q2_pt;
double truth_q2_eta;
double truth_q2_phi;
double truth_q2_m;

double truth_t1_pt;
double truth_t1_eta;
double truth_t1_phi;
double truth_t1_m;

double truth_t2_pt;
double truth_t2_eta;
double truth_t2_phi;
double truth_t2_m;

double truth_W1_pt;
double truth_W1_eta;
double truth_W1_phi;
double truth_W1_m;

double truth_W2_pt;
double truth_W2_eta;
double truth_W2_phi;
double truth_W2_m;

double truth_H_pt;
double truth_H_eta;
double truth_H_phi;
double truth_H_m;

vector<int>*    fspart_id;
vector<double>* fspart_pt;
vector<double>* fspart_eta;
vector<double>* fspart_phi;
vector<double>* fspart_m;



///////////////////////////
//output file and tree
///////////////////////////
TTree *treeout;
TFile *fileout;

///////////////////////////
//for temporary storage
///////////////////////////
TLorentzVector truth_q1;
TLorentzVector truth_q2;
TLorentzVector truth_t1;
TLorentzVector truth_t2;
TLorentzVector truth_W1;
TLorentzVector truth_W2;
TLorentzVector truth_H;

TLorentzVector tau_axis1;
TLorentzVector tau_axis2;
TLorentzVector tau_axis3;
TLorentzVector Tsubjet1;
TLorentzVector Tsubjet2;
TLorentzVector Tsubjet3;
TLorentzVector particle;
TLorentzVector truth;

///////////////////////////////////////////////////////////////
TLorentzVector Truth_W;
TLorentzVector Truth_T;
TLorentzVector Truth_jet;
TLorentzVector Wdaughter1;
TLorentzVector Wdaughter2;
TLorentzVector Tdaughter1;
TLorentzVector Tdaughter2;

TLorentzVector total_truth;
TLorentzVector total_reco;
TLorentzVector total0_truth;
TLorentzVector total0_reco;

int id;
double pruned_mjet, pruned_ptjet, pruned_etajet;
double trimmed_mjet, trimmed_ptjet, trimmed_etajet;
double mjet, ptjet, etajet;
double mjetR, ptjetR, etajetR;
double Tprun_volatility, Ttrim_volatility, TakTrecl_volatility, TkTrecl_volatility, Tsubj_volatility, Ttau2_volatility;
double Tsubj_angle;

const double M0 = 0.1;
const double Rfat = 1.0;
const double zcut = 0.1, dcut0 = 0.5;
const double Rfilt0 = 0.3, fcut0 = 0.05;
const double jet_pt_cut_low = 800, jet_pt_cut_up = 1000;
const double jet_eta_cut_low = -1.2, jet_eta_cut_up = 1.2;
const double jet_mass_cut_low = 70, jet_mass_cut_up = 90;
///////////////////////////////////////////////////////////////


int jetflavor;

int    tempJet_flavor;
double tempJet_pt;
double tempJet_eta;
double tempJet_phi;
double tempJet_m;
double tempJet_Tau21;
double tempJet_Tau32;
double tempJet_D2;
double tempJet_TJet_m1;
double tempJet_TJet_m2;

double tempJet_T2jet_angle;
double tempJet_T2jet;
double tempJet_T3jet_angle;
double tempJet_T3jet;
double tempJet_Tpruning;
double tempJet_Ttrimming;
double tempJet_Taktreclustering;
double tempJet_Tktreclustering;

///////////////////////////
//output tree branches
///////////////////////////
int    NumberOfVertices;

vector<int>    TruthRaw_flavor;
vector<double> TruthRaw_pt;
vector<double> TruthRaw_eta;
vector<double> TruthRaw_phi;
vector<double> TruthRaw_m;
vector<double> TruthRaw_Tau21;
vector<double> TruthRaw_Tau32;
vector<double> TruthRaw_D2;
vector<double> TruthRaw_TJet_m1;
vector<double> TruthRaw_TJet_m2;
vector<double> TruthRaw_T2jet_angle;
vector<double> TruthRaw_T2jet;
vector<double> TruthRaw_T3jet_angle;
vector<double> TruthRaw_T3jet;
vector<double> TruthRaw_Tpruning;
vector<double> TruthRaw_Ttrimming;
vector<double> TruthRaw_Taktreclustering;
vector<double> TruthRaw_Tktreclustering;

vector<int>    TruthRawTrim_flavor;
vector<double> TruthRawTrim_pt;
vector<double> TruthRawTrim_eta;
vector<double> TruthRawTrim_phi;
vector<double> TruthRawTrim_m;
vector<double> TruthRawTrim_Tau21;
vector<double> TruthRawTrim_Tau32;
vector<double> TruthRawTrim_D2;
vector<double> TruthRawTrim_TJet_m1;
vector<double> TruthRawTrim_TJet_m2;
vector<double> TruthRawTrim_T2jet_angle;
vector<double> TruthRawTrim_T2jet;
vector<double> TruthRawTrim_T3jet_angle;
vector<double> TruthRawTrim_T3jet;
vector<double> TruthRawTrim_Tpruning;
vector<double> TruthRawTrim_Ttrimming;
vector<double> TruthRawTrim_Taktreclustering;
vector<double> TruthRawTrim_Tktreclustering;

///////////////////////////
//extra functions
///////////////////////////
void ResetBranches();

int GetJetTruthFlavor(TLorentzVector jettemp,
                      TLorentzVector truth_t1,
                      TLorentzVector truth_t2,
                      TLorentzVector truth_W,
                      TLorentzVector truth_Z,
                      TLorentzVector truth_H,
                      int debug);

vector<PseudoJet> ToyCalorimeter(vector<PseudoJet> truth_particles);

double GetTau21(PseudoJet& input);
double GetTau32(PseudoJet& input);

double T_Nsubjettiness(int N, PseudoJet& input, double beta_min, double beta_max, int N_beta);
double T_NsubjettinessRatio(int N_num, int N_den, PseudoJet& input, double beta_min, double beta_max, int N_beta);

double T_EnergyCorrelator_C2(PseudoJet& input, double beta_min, double beta_max, int N_beta);
double T_EnergyCorrelator_D2(PseudoJet& input, double beta_min, double beta_max, int N_beta);
double T_EnergyCorrelator_C3(PseudoJet& input, double beta_min, double beta_max, int N_beta);

///=========================================
/// Telescoping Subjet
///=========================================

struct TSub{
  double min_angle;
  double volatility;
};

TSub T_2Subjet(PseudoJet& input, double R_min, double R_max, int N_R);
TSub T_3Subjet(PseudoJet& input, double R_min, double R_max, int N_R);

double T_Pruning(PseudoJet& input, double dcut_min, double dcut_max, int N_dcut);
double T_Trimming(PseudoJet& input, double fcut_min, double fcut_max, int N_fcut);
double T_AkTreclustering(PseudoJet& input, double R_min, double R_max, int N_R);
double T_kTreclustering(PseudoJet& input, double R_min, double R_max, int N_R);

#endif
