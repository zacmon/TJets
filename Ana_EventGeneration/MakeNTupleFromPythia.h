#ifndef MAKENTUPLEFROMPYTHIA_H
#define MAKENTUPLEFROMPYTHIA_H

#define CERRD cout<<"Problem on "<<__FILE__<<"  "<<__LINE__<<endl;

#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <vector>
#include <string>

#include "Pythia8/Pythia.h"

//used for the output particle ntuple
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

std::vector<int>    fspart_id;
std::vector<double> fspart_pt;
std::vector<double> fspart_eta;
std::vector<double> fspart_phi;
std::vector<double> fspart_m;

#endif