export PYTHIAPATH=/Users/zacmon/Research/pythia8226

echo
echo "Compiling with : "
echo "$ROOTSYS    : "${ROOTSYS}
echo "$PYTHIAPATH : "${PYTHIAPATH}
echo

g++ MakeNTupleFromPythia.cc $PYTHIAPATH/lib/libpythia8.a -o MakeNTupleFromPythia.exe -I$ROOTSYS/include/root  -I$PYTHIAPATH/include  -rpath $ROOTSYS/lib `$ROOTSYS/bin/root-config --glibs` -std=c++11
