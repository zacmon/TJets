#! /usr/bin/env bash
source compile_pythia.sh
./MakeNTupleFromPythia.exe tt pythia_tt_$1.root $1
source compile_ntupler.sh
./NTupler.exe tt pythia_tt_$1.root ntuple_tt_$1.root
