#compile
source compile.sh

#generate dijet
./MakeNTupleFromPythia.exe jj pythia_dijet.root     1000

#generate W'-->WZ
./MakeNTupleFromPythia.exe ww pythia_ww.root        1000

#generate ttbar
./MakeNTupleFromPythia.exe tt pythia_ttbar.root     1000