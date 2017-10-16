export PYTHIAPATH=/phys/users/zacmon/Packages/pythia8226

echo
echo "Compiling with : "
echo "$ROOTSYS    : "${ROOTSYS}
echo "$PYTHIAPATH : "${PYTHIAPATH}
echo

g++ MakeNTupleFromPythia.cc $PYTHIAPATH/lib/libpythia8.a -o MakeNTupleFromPythia.exe -I$ROOTSYS/include  -I$PYTHIAPATH/include  -Wl,-rpath=$ROOTSYS/lib `$ROOTSYS/bin/root-config --glibs` -std=c++1y
