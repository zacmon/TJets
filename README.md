# TelescopingJets
Analysis tools for the Telescoping Jets project.  These tools describe how to set up and generate sets of events using Pythia and fastjet and then analyze these with Root.

# First Time Setup
In addition to this code, you will need to compile and install the following additional packages, whose libraries are used within the generation and analysis framework.

## ROOT: https://root.cern.ch/building-root   

## Pythia8: http://home.thep.lu.se/~torbjorn/Pythia.html 
- This must be compiled/installed in a location that is then accessible when you compile the EventGenerator.  You will need to point the compilation to the Pythia directory.
- If you are not familiar with Pythia then it is worthwhile working through the worksheet outlined (in this tutorial)[http://home.thep.lu.se/~torbjorn/pythia8/worksheet8183.pdf] since the EventGeneration code is built starting from this.

## fastjet : http://fastjet.fr/
- This was developed with fastjet-3.2.1
- - This must be compiled/installed in a location that is then accessible when you compile the next package, fastjet contrib.
- This must be compiled/installed in a location that is then accessible when you compile the EventAnalyzer.  You will need to point the compilation to the fastjet directory so that it can successfully find the configurations

## fastjet contrib : http://fastjet.hepforge.org/contrib/ \\
- It is necessary to point the code to the fastjet install directory as shown in the install directions.  This is done via the `--fastjet-config=/path/to/fastjet-config` argument.  
- If you installed fastjet in the directory above the fjcontrib directory, this means pointing it to the fastjet-config at ../fastjet-3.1.1/fastjet-config.  If you do this correctly, then you should see the /contrib directory housed in `/usr/local/include/fastjet` after running `make install`.


# Running EventGeneration : `Ana_EventGeneration`
This is the first state of analysis in which the Pythia8 event generator is used to produce a set of events for any process you wish and perform the subsequent shower and hadronization to produce a set of final state particles that can then be used for analysis via jet building.

Start by compiling the code in the `Ana_EventGeneration` directory which will generate Pythia events

```
cd Ana_EventGeneration
source compile_pythia.sh
```

Next, execute an example run command for this executable

```
./MakeNTupleFromPythia.exe ww pythia_ww.root 200
```

This will generate 200 events and save them to the pythia_ww.root output file.  The data which is all of the final state particles (those that would be detected) along with the information of the initial particles coming from the collision.  In this case, the process generates two high momentum W bosons, and therefore the four momenta of the two W bosons are saved for future truth matching.

Next, to process these events you just generated and build jets, you must first compile the _NTupler_ as
```
source compile_ntupler.sh
```

after which you will run the executable as
```
./NTupler.exe ww pythia_ww.root ntuple_ww.root
```

which will run jet clustering and labelling (the matching of W particles to the constructed jets) on the `pythia_ww.root` file and save the output in the `ntuple_ww.root` file.

All of this can be done in one command using the `GenerateFull.py` script which itself takes arguments about the process to generate, along with those described by looking at the help menu :
```
meehan:Ana_EventGeneration > python GenerateFull.py --help
usage: GenerateFull.py [-h] [--type TYPE] [--nevents NEVENTS]
                       [--jobstart JOBSTART] [--jobend JOBEND] [--recompile]
                       [--pythia] [--ntuple] [--debug]

optional arguments:
  -h, --help           show this help message and exit
  --type TYPE          Type of the physical process you wish to generate.
                       Possible values are ["jj"=parton jets, "ww"=boosted W
                       bosons], "tt"=boosted top quarks] (default: ww)
  --nevents NEVENTS    The number of events to create in a single job
                       (default: 1000)
  --jobstart JOBSTART  When you want to generate multiple jobs, this is the
                       lower job number in the list of [jobStart,jobEnd] jobs.
                       (default: 0)
  --jobend JOBEND      When you want to generate multiple jobs, this is the
                       upper job number in the list of [jobStart,jobEnd] jobs.
                       (default: 1)
  --recompile          Flag that will recompile the program before running
                       (default: False)
  --pythia             Flag that will run the generation of Pythia events
                       (default: False)
  --ntuple             Flag that will run the NTupler over pre-made Pythia
                       particle ntuples (default: False)
  --debug              Flag that will turn on the debugging printout (default:
                       False)
meehan:Ana_EventGeneration >
```
To compile the _NTuplerTest_, 
```
source compile_ntupler_test.sh
```
after which run the executable as
```
./NTuplerTest.exe ww pythia_ww.root ntuple_ww.root
```
_NTuplerTest_ does exactly the same thing as _NTupler_, with the addition of the part from line 157 to line 906. It reads in Alex's samples and calculates various Telescoping jet volatilities and other jet observables. It also outputs histograms as .root files. Actually, we don't need this part because the MiniNTupleAnalysis will do. That's only for Yang-Ting to cross-check plots in his BOOST 2015 talk.

Alex's samples are available at the UW TeV cluster (tev01.phys.washington.edu) in the folder /ytchien/AlexSample. There are five of them:
```
mc12_8TeV.14791X.Pythia8_AU2CT10_jetjet_JZXW.leading_jets.root
mc12_8TeV.158864.Pythia8_AU2MSTW2008LO_Wprime_WZ_llqq_m1000.leading_jets.root
mc12_8TeV.158874.Pythia8_AU2MSTW2008LO_Wprime_WZ_llqq_m2000.leading_jets.root
mc12_8TeV.110903.Pythia8_AU2MSTW2008LO_zprime1000_tt.leading_jets.root
mc12_8TeV.110907.Pythia8_AU2MSTW2008LO_zprime2000_tt.leading_jets.root
```
for QCD (1st), W (2nd and 3rd), top (4th and 5th) samples. Yang-Ting used the first three to reproduce the BOOST 2015 results. Two pT bins are considered: 1) 300 GeV to 500 Gev, and 2) 800 Gev to 1000 GeV.

For the (300, 500) bin, the sample
```
mc12_8TeV.158864.Pythia8_AU2MSTW2008LO_Wprime_WZ_llqq_m1000.leading_jets.root
```
is used. Also, `Rfilt = 0.2` in T-trimming and `R_min = 0.1` in T-reclustering and T-subjet. For the (800, 1000) bin, the following one
```
mc12_8TeV.158874.Pythia8_AU2MSTW2008LO_Wprime_WZ_llqq_m2000.leading_jets.root
```
is used, and `Rfilt = 0.1` in T-trimming and `R_min = 0.05` in T-reclustering and T-subjet.

# Running MiniNTupleAnalysis

This is the final stage of analysis where the analyzer uses any method they wish to analyze the flat ntuple of jets from the previous stage and produce final results.  The preferred method is to use the (TTree::Draw())[https://root.cern.ch/root/html/TTree.html#TTree:Draw@2] method to quickly make histograms that can be used for further calculation.  This means that in principle there should be *no event loop* in this part of the code.  This is made more amenable by implementing all analysis tools with PyROOT.

All of the existing scripts are stored in the `Ana_MiniNTupleAnalysis`.  You should explore them on your own to understand what they do.

If all you want to do is create a comparison of signal and background, you can execute
```
python SimplePlot.py
```
noting that you will need to point to the proper output file directory after having executed, in the `Ana_EventGeneration` stage, the following two commands
```
python GenerateFull.py --type dijet --nevents 5000 --pythia --ntuple --recompile
python GenerateFull.py --type ww --nevents 5000 --pythia --ntuple --recompile
```
