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
PROGRAMS  = MakeNTupleFromPythia.exe 
MYLIBCXX  = $(subst .h,.cc,$(shell ls *.h)) 
MYLIBS    = $(subst .cc,.o,$(shell ls $(MYLIBCXX)))
ROOTSYS   = $(shell echo $ROOTSYS)
# add pythia dependency; the binary and header files -- wesley
PYTHIAPATH = /home/wuy45/Temp/pythia8240
PYTHIALIBS = -L$(PYTHIAPATH)/lib/libpythia8.a
PYTHIAFLAGS = -I$(PYTHIAPATH)/include

# add fastjet dependency;
FASTJETPATH = /home/wuy45/Temp/fastjet-install
FASTJETFLAGS = $(shell $(FASTJETPATH)/bin/fastjet-config --cxxflags) -lNsubjettiness -lEnergyCorrelator
FASTJETFLIBS = $(shell $(FASTJETPATH)/bin/fastjet-config --libs) -lNsubjettiness -lEnergyCorrelator

# -fPIC option needed only for compilation into a shared library
CXXFLAGS = $(shell root-config --cflags) $(PYTHIAFLAGS) $(FASTJETFLAGS) -fPIC 
ROOTLIBS = $(shell root-config --libs)

# all libs
LIBS     = $(ROOTLIBS) -lMinuit $(PYTHIALIBS) $(FASTJETLIBS)
OPT      = -ggdb -Wshadow # for debug
#OPT      = -ggdb -Wold-style-cast -Wshadow # for debug
#OPT      = -O3   # optimized
GCC       = g++
SOFLAGS   = -shared 
SHAREDLIB = $(shell echo ${PWD})/ana_event_gen.so
#---------DEFAULT ACTION-------------------------
all: $(PROGRAMS)

# compiling the executable 
%.exe: %.o *.h *.cc $(SHAREDLIB)
	echo creating the executable $@
	$(GCC) $(OPT) $(CXXFLAGS) -L. -o $@ $< $(SHAREDLIB) $(LIBS) 
# creating the shared library
$(SHAREDLIB): $(MYLIBS) 
	echo "Building shared library ... "
	rm -f $(SHAREDLIB)
	$(GCC) -L$(shell root-config --libdir) $(SOFLAGS) $(MYLIBS) -o $(SHAREDLIB)
	echo "Done"
#
%.o: %.cc *.h 
	$(GCC) $(OPT) $(CXXFLAGS) -o $@ -c $<
#

#----------CLEANUP------------------
clean:
	rm -f *.o *.so *.exe *~ $(PROGRAMS) $(SHAREDLIB)


#--COMPILING if the shared library does not work for your operating system---
# static linking of all object files 
#%.exe: %.o *.h *.cxx $(MYLIBS) DataItemCreator.o
#	$(GCC) $(OPT) $(CXXFLAGS) -L. -o $@ $< $(MYLIBS) $(LIBS) 

