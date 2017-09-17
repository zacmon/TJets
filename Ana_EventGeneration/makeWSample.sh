#! /usr/bin/env bash
source compile_pythia.sh
./MakeNTupleFromPythia.exe ww pythia_tt_test$1.root $1
source compile_ntupler.sh
./NTupler.exe ww pythia_tt_test$1.root ntuple_tt_test$1.root
