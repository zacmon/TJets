export FJINSTALL=/Users/ytchien/Research/fastjet-install

echo
echo "Compiling with : "
echo "$ROOTSYS    : "${ROOTSYS}
echo "$FJINSTALL : "${FJINSTALL}
echo

g++ -o NTuplerTest.exe NTuplerTest.cc TelescopingJets.cc -I$ROOTSYS/include -I$FJINSTALL/include/fastjet  -rpath $ROOTSYS/lib `$ROOTSYS/bin/root-config --glibs` `$FJINSTALL/bin/fastjet-config --cxxflags --libs --plugins` -lNsubjettiness -lEnergyCorrelator   -std=c++11
