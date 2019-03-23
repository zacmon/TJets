# Makefile for the DataAnalysis Example
# Eckhard von Toerne, Bonn U., 2009
#
# the makefiles contains all rules to create the target
# usage: make name_of_target
# the default target is "all"
#
# calling make from the command line starts the target creation pocess 
# defined by the makefile.
# The makefile consists of variable definitions and target declarations.
# The target has dependencies which have to be created before the target 
# the dependencies must also have one or several target rules.
# Syntax: 
# target_name: list of dependencies
# -TAB-target rule
# -TAB-target rule
# -TAB-target rule
#
# the make utility proceeds recursively from the main target trough the 
# chain of dependencies and their target rules to the targets that need 
# to be created first - and then creates all targets in the correct order.
#---------------
# $() constructs
# variables: defined by:  VAR = xxx , addressed via $(VAR)
# shell commands:         $(shell command)
# string substitution:    $(subst string,replacement,targetstring)
#----------------
# placeholders
# %.exe : any file that matches *.exe
# $@ Name of target
# $< the first dependence 
# $+ List of all dependencies
# $^ list of all dependencies with double listings removed
#----------------
# Important:
# always start a target rule with a TAB
#---------GENERAL SETTINGS--------------
#PROGRAMS  = $(subst .cxx,.exe,$(shell ls run*.cxx))
PROGRAMS  = MakeNTupleFromPythia.exe NTupler.exe

ROOTSYS   = $(shell echo $ROOTSYS)

# PYTHIAPATH and  FASTJETPATH are set at variables.sh
# -fPIC option needed only for compilation into a shared library
CXXFLAGS = -fPIC #-Wl,--verbose
ROOTFLAGS = $(shell root-config --cflags)
FASTJETFLAGS = $(shell ${FASTJETPATH}/bin/fastjet-config --cxxflags)
PYTHIAFLAGS = $(shell ${PYTHIAPATH}/bin/pythia8-config --cxxflags)

ROOTLIBS = $(shell root-config --glibs) -lMinuit
FASTJETLIBS = $(shell ${FASTJETPATH}/bin/fastjet-config --libs --plugins) -lNsubjettiness -lEnergyCorrelator
PYTHIALIBS = $(shell ${PYTHIAPATH}/bin/pythia8-config --libs)

# all libs
# LIBS     = $(ROOTLIBS) -lMinuit $(PYTHIALIBS) #$(FASTJETLIBS)
OPT      = -ggdb -Wshadow # for debug
#OPT      = -ggdb -Wold-style-cast -Wshadow # for debug
#OPT      = -O3   # optimized
GCC       = g++
SOFLAGS   = -shared 
SHAREDLIB = $(shell echo ${PWD})/ana_event_gen.so
#---------DEFAULT ACTION-------------------------

all: $(PROGRAMS)

# compiling the executable 
MakeNTupleFromPythia.exe:  MakeNTupleFromPythia.cc MakeNTupleFromPythia.h 
	echo creating the executable $@
	$(GCC) $(OPT) $(CXXFLAGS) $(ROOTFLAGS) $(PYTHIAFLAGS) $<  -o $@ $(ROOTLIBS) $(PYTHIALIBS)

NTupler.exe: NTupler.cc NTupler.h TelescopingJets.cc
	echo creating the executable $@
	$(GCC) $(OPT) $(CXXFLAGS) $(ROOTFLAGS) $(FASTJETFLAGS) $+ -o $@ $(ROOTLIBS) $(FASTJETLIBS)

#----------CLEANUP------------------
clean:
	rm -f *.o *.so *.exe *~ $(PROGRAMS)


#--COMPILING if the shared library does not work for your operating system---
# static linking of all object files 
#%.exe: %.o *.h *.cxx $(MYLIBS) DataItemCreator.o
#	$(GCC) $(OPT) $(CXXFLAGS) -L. -o $@ $< $(MYLIBS) $(LIBS) 

