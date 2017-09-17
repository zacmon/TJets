#! /usr/bin/env bash
source compile_pythia.sh
./MakeNTupleFromPythia.exe dijet pythia_dijet_test$1.root $1
source compile_ntupler.sh
./NTupler.exe dijet pythia_dijet_test$1.root ntuple_dijet_test$1.root
