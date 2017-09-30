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
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdio>
#include <algorithm>
#include <iterator>

#include "TelescopingJets.hh"

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>

#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"


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

std::vector<int>*    fspart_id;
std::vector<double>* fspart_pt;
std::vector<double>* fspart_eta;
std::vector<double>* fspart_phi;
std::vector<double>* fspart_m;



///////////////////////////
//output file and tree
///////////////////////////
TTree *treeout;
TFile *fileout;

std::ifstream data_file;
std::ofstream Tjet_variable_file;

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
TLorentzVector tau_axis4;
TLorentzVector tau_axis5;
TLorentzVector Tsubjet1;
TLorentzVector Tsubjet2;
TLorentzVector Tsubjet3;
TLorentzVector Tsubjet4;
TLorentzVector Tsubjet5;

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
const double massW = 80.4;
const double Rfat = 1.0;
const double zcut = 0.1, dcut0 = 0.5;
const double Rfilt0 = 0.3, fcut0 = 0.05;
const double jet_pt_cut_low = 800, jet_pt_cut_up = 1000;
const double jet_eta_cut_low = -1.2, jet_eta_cut_up = 1.2;
const double jet_mass_cut_low = 150, jet_mass_cut_up = 190;
///////////////////////////////////////////////////////////////


int jetflavor;

///////////////////////////
//output tree branches
///////////////////////////
int    NumberOfVertices;

std::vector<int>    TruthRaw_flavor;
std::vector<double> TruthRaw_pt;
std::vector<double> TruthRaw_eta;
std::vector<double> TruthRaw_phi;
std::vector<double> TruthRaw_m;
std::vector<double> TruthRaw_Tau21;
std::vector<double> TruthRaw_Tau32;
std::vector<double> TruthRaw_D2;
std::vector<double> TruthRaw_TJet_m1;
std::vector<double> TruthRaw_TJet_m2;
std::vector<double> TruthRaw_T1jet_angle;
std::vector<double> TruthRaw_T1jet;
std::vector<double> TruthRaw_T2jet_angle;
std::vector<double> TruthRaw_T2jet;
std::vector<double> TruthRaw_T3jet_angle;
std::vector<double> TruthRaw_T3jet_angle1;
std::vector<double> TruthRaw_T3jet_angle2;
std::vector<double> TruthRaw_T3jet;
std::vector<double> TruthRaw_T3jet_W;
std::vector<double> TruthRaw_T3jet_mW;
std::vector<double> TruthRaw_T4jet_angle;
std::vector<double> TruthRaw_T4jet;
std::vector<double> TruthRaw_T5jet_angle;
std::vector<double> TruthRaw_T5jet;
std::vector<double> TruthRaw_Tpruning;
std::vector<double> TruthRaw_Ttrimming;
std::vector<double> TruthRaw_Taktreclustering;
std::vector<double> TruthRaw_Tktreclustering;

std::vector<int>    TruthRawTrim_flavor;
std::vector<double> TruthRawTrim_pt;
std::vector<double> TruthRawTrim_eta;
std::vector<double> TruthRawTrim_phi;
std::vector<double> TruthRawTrim_m;
std::vector<double> TruthRawTrim_Tau21;
std::vector<double> TruthRawTrim_Tau32;
std::vector<double> TruthRawTrim_D2;
std::vector<double> TruthRawTrim_TJet_m1;
std::vector<double> TruthRawTrim_TJet_m2;
std::vector<double> TruthRawTrim_T1jet_angle;
std::vector<double> TruthRawTrim_T1jet;
std::vector<double> TruthRawTrim_T2jet_angle;
std::vector<double> TruthRawTrim_T2jet;
std::vector<double> TruthRawTrim_T3jet_angle;
std::vector<double> TruthRawTrim_T3jet_angle1;
std::vector<double> TruthRawTrim_T3jet_angle2;
std::vector<double> TruthRawTrim_T3jet;
std::vector<double> TruthRawTrim_T3jet_W;
std::vector<double> TruthRawTrim_T3jet_mW;
std::vector<double> TruthRawTrim_T4jet_angle;
std::vector<double> TruthRawTrim_T4jet;
std::vector<double> TruthRawTrim_T5jet_angle;
std::vector<double> TruthRawTrim_T5jet;
std::vector<double> TruthRawTrim_Tpruning;
std::vector<double> TruthRawTrim_Ttrimming;
std::vector<double> TruthRawTrim_Taktreclustering;
std::vector<double> TruthRawTrim_Tktreclustering;

std::vector<std::vector<double>> TruthRawTrim_T1Volatility;
std::vector<std::vector<double>> TruthRawTrim_T2Volatility;
std::vector<std::vector<double>> TruthRawTrim_T3Volatility;

std::vector<double> TruthRawTrim_v32;
std::vector<std::vector<double>> TruthRawTrim_T1masses;
std::vector<std::vector<double>> TruthRawTrim_T2masses;
std::vector<std::vector<double>> TruthRawTrim_T3masses;

std::vector<double> TruthRawTrim_T2jet_massW;
std::vector<double> TruthRawTrim_T2jet_volatilityW;
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

std::vector<fastjet::PseudoJet> ToyCalorimeter(std::vector<fastjet::PseudoJet> truth_particles);

double GetTau21(fastjet::PseudoJet& input);
double GetTau32(fastjet::PseudoJet& input);

double T_Nsubjettiness(int N, fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta);
double T_NsubjettinessRatio(int N_num, int N_den, fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta);

double T_EnergyCorrelator_C2(fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta);
double T_EnergyCorrelator_D2(fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta);
double T_EnergyCorrelator_C3(fastjet::PseudoJet& input, double beta_min, double beta_max, int N_beta);

///=========================================
/// Telescoping Subjet
///=========================================

struct TSub {
    double minAngle = -1;
    double volatility = -1;
    std::vector<double> volVec;
    std::vector<double> masses;
};

struct T2Sub {
    double minAngle = -1;
    double volatility = -1;
    std::vector<double> volVec;
    std::vector<double> masses;

    double massW = -1;
    double volatilityW = -1;
};
    
struct T3Sub {
    double minAngle = -1;
    double volatility = -1;
    std::vector<double> volVec;
    std::vector<double> masses;

    double midAngle = -1;
    double maxAngle = -1;

    double massW = -1;
    double volatilityW = -1;
};

TSub TNSubjet(fastjet::PseudoJet& input, unsigned int numSubjets, double minRadius, double maxRadius, int numRadii);
T2Sub T2Subjet(fastjet::PseudoJet& input, double minRadius, double maxRadius, int numRadii);
T3Sub T3Subjet(fastjet::PseudoJet& input, double minRadius, double maxRadius, int numRadii);
double T_Pruning(fastjet::PseudoJet& input, double dcut_min, double dcut_max, int N_dcut);
double T_Trimming(fastjet::PseudoJet& input, double fcut_min, double fcut_max, int N_fcut);
double T_Reclustering(fastjet::PseudoJet& input, int algorithm, double R_min, double R_max, int N_R);

#endif
