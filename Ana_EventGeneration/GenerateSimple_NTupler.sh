#compile
source compile.sh

#generate dijet
./NTupler.exe 0 ../Ana_EventGeneration/pythia_dijet.root ntuple_dijet.root

#generate W'-->WZ
./NTupler.exe 1 ../Ana_EventGeneration/pythia_ww.root    ntuple_ww.root

#generate ttbar
./NTupler.exe 2 ../Ana_EventGeneration/pythia_ttbar.root ntuple_ttbar.root