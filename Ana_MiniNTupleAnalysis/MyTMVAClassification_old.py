#!/usr/bin/env python


#example: python MyTMVAClassification.py   truth  akt10truth_trim_pt,akt10truth_trim_mass "pt>0,pt<1000,mass>0,mass<200,pass_selection==1"  "truth_tau2_WTA,truth_tau1_WTA"  KNN

# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser

# Import ROOT classes
from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut

# check ROOT version, give alarm if 5.18 
if gROOT.GetVersionCode() >= 332288 and gROOT.GetVersionCode() < 332544:
    print "*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA"
    print "*** does not run properly (function calls with enums in the argument are ignored)."
    print "*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),"
    print "*** or use another ROOT version (e.g., ROOT 5.19)."
    sys.exit(1)

# Import TMVA classes from ROOT
from ROOT import TMVA



# Print usage help
def usage():
    print " "
    print "Usage: python %s [options]"

# Main routine
def main():

    NTupName   = "JetTree"
    alg        = "TruthRawTrim"
    cuts       =  ["eta>-1.2","eta<1.2","pt>300","pt<5000","m>61","m<85"]
    vars       = ["Tau21","T2jet"]

    #NTupName   = "varTree"    
    verbose    = True
    
    #alg        = "AK10LCTRIMF5R20"
    spectators = ["m"]
    #cuts       = ["eta>-1.2","eta<1.2","pt>200","pt<350","m>61","m<85","TruthRecoMatch==1"]
    #vars       = ["TauWTA2TauWTA1","ZCUT12","Dip23","TJetVol","ActiveArea","PullC10","Angularity"]
    methods    = "Likelihood"

    print "Starting and getting arguments:"
    allargs = sys.argv[1:]    
    if len(allargs)<5:
        print "You input these args"
        print allargs
        print "Not enough args, please try again"
        return 1
    else:
        alg        = allargs[0]
        spectators = allargs[1].split(",")
        cuts       = allargs[2].split(",")
        vars       = allargs[3].split(",")
        methods    = allargs[4]
    
    print "Running with args:"
    print "  alg        = ",alg        
    print "  spectators = ",spectators 
    print "  cuts       = ",cuts       
    print "  vars       = ",vars       
    print "  methods    = ",methods    
    

    # Print methods
    mlist = methods.replace(' ',',').split(',')
    print "=== TMVAClassification: use method(s)..."
    for m in mlist:
        if m.strip() != '':
            print "=== - <%s>" % m.strip()



    #===============================
    #Read training and test data
    #===============================
    #InputDir = "../gen_20170529/"
    InputDir="../Ana_EventGeneration/"
    #InputDir="~/Downloads/"
    print "Getting inputs from: ",InputDir
    #s1 = TFile(InputDir+"ntuple_ttbar_2000.root");
    #b1 = TFile(InputDir+"ntuple_dijet_800_1400.root");
    s1 = TFile(InputDir+"ntuple_tt_test10000.root");
    b1 = TFile(InputDir+"ntuple_dijet_test10000.root");

    # Output file
    OutFileName="testout.root"
    outputFile = TFile( OutFileName, 'RECREATE' )
    
    # Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for more factory options)
    # All TMVA output can be suppressed by removing the "!" (not) in 
    # front of the "Silent" argument in the option string
    factory = TMVA.Factory( "TMVAClassification", outputFile, 
                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

    # Set verbosity
    factory.SetVerbose( verbose )

    # data loader
    dataloader=TMVA.DataLoader("dataset")

#     weight=""
#     weight+="pass_selection*EventWeight*CrossSection*("
#     weight+=alg+"_pt>"+pt1+" && "
#     weight+=alg+"_pt<"+pt2
#     if m1!="0":
#         weight+=" && "+alg+"_m>"+m1+" && "
#         weight+=alg+"_m<"+m2
#     weight+=")"
#     
#     #Get signal and background histograms
#     if variable=="mass":
#         histname = "akt10"+alg+"_trim_"+variable
#     else:
#         histname = alg+"_"+variable
    
    #======================================
    #Predefined cuts - for isntance on M(j1)
    #======================================
    mycuts = "1.0"
    mycutb = "1.0"

    for cut in cuts:
        placecut=cut
        if cut[:2]=="pt" or cut[:3]=="eta" or cut[:4]=="m":
            placecut = "* ("+alg+"_"+cut+")"
        else:
            placecut="* ("+cut+") "
        mycuts += placecut
        mycutb += placecut

    
    print "MyCutsSig: ",mycuts
    print "MyCutsBkg: ",mycutb

    #===================================
    #Spectator variables from tree
    #=====================================
    for spec in spectators:
        dataloader.AddSpectator( spec, 'F' )
        
    #===================================
    #MVA variables from tree
    #=====================================
    for var in vars:
        dataloader.AddVariable( var , 'F' )

    #===============================
    #Read training and test data
    #===============================
    print "Getting trees ... ",NTupName
    st1 = s1.Get(NTupName)
    bt1 = b1.Get(NTupName)

    #=========================================
    # global event weights per tree (see below for setting event-wise weights)
    #=========================================
    ws1 = 1.0
    wb1 = 1.0

    #=========================================
    # You can add an arbitrary number of signal or background trees
    #=========================================
    dataloader.AddSignalTree    ( st1, ws1 );
    #SCHSU
    #dataloader.SetSignalWeightExpression("EventWeight*CrossSection");

    dataloader.AddBackgroundTree( bt1, wb1 );
    #dataloader.SetBackgroundWeightExpression("EventWeight*CrossSection");
    
    # Here, the relevant variables are copied over in new, slim trees that are
    # used for TMVA training and testing
    # "SplitMode=Random" means that the input events are randomly shuffled before
    # splitting them into training and test samples
    mycutSig = TCut(mycuts)
    mycutBkg = TCut(mycutb)
    
    dataloader.PrepareTrainingAndTestTree( mycutSig, mycutBkg,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" )

    # --------------------------------------------------------------------------------------------------

    # ---- Book MVA methods
    #
    # please lookup the various method configuration options in the corresponding cxx files, eg:
    # src/MethoCuts.cxx, etc, or here: http:#tmva.sourceforge.net/optionRef.html
    # it is possible to preset ranges in the option string in which the cut optimisation should be done:
    # "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

    # Cut optimisation
    if "Cuts" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "Cuts",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" )

    if "CutsD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "CutsD",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" )

    if "CutsPCA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "CutsPCA",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" )

    if "CutsGA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "CutsGA",
                            "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" )

    if "CutsSA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "CutsSA",
                            "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" )

    # Likelihood ("naive Bayes estimator")
    if "Likelihood" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "Likelihood",
                            "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )

    # Decorrelated likelihood
    if "LikelihoodD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodD",
                            "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" )

    # PCA-transformed likelihood
    if "LikelihoodPCA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodPCA",
                            "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ) 

    # Use a kernel density estimator to approximate the PDFs
    if "LikelihoodKDE" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodKDE",
                            "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ) 

    # Use a variable-dependent mix of splines and kernel density estimator
    if "LikelihoodMIX" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodMIX",
                            "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ) 

    # Test the multi-dimensional probability density estimator
    # here are the options strings for the MinMax and RMS methods, respectively:
    #      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
    #      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
    if "PDERS" in mlist:
# DEFAULT
#         factory.BookMethod( dataloader, TMVA.Types.kPDERS, "PDERS",
#                             "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" )
# CHOOSE RIGID VOLUME SO IT DOESNT TAKE SO LONG
        factory.BookMethod( dataloader, TMVA.Types.kPDERS, "PDERS",
                            "!H:!V:NormTree=T:VolumeRangeMode=Unscaled:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" )


    if "PDERSD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kPDERS, "PDERSD",
                            "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" )

    if "PDERSPCA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kPDERS, "PDERSPCA",
                             "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" )

   # Multi-dimensional likelihood estimator using self-adapting phase-space binning
    if "PDEFoam" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kPDEFoam, "PDEFoam",
                            "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" )

    if "PDEFoamBoost" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kPDEFoam, "PDEFoamBoost",
                            "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" )

    # K-Nearest Neighbour classifier (KNN)
    if "KNN" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kKNN, "KNN",
                            "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" )

    # H-Matrix (chi2-squared) method
    if "HMatrix" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kHMatrix, "HMatrix", "!H:!V" )

    # Linear discriminant (same as Fisher discriminant)
    if "LD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

    # Fisher discriminant (same as LD)
    if "Fisher" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

    # Fisher with Gauss-transformed input variables
    if "FisherG" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFisher, "FisherG", "H:!V:VarTransform=Gauss" )

    # Composite classifier: ensemble (tree) of boosted Fisher classifiers
    if "BoostedFisher" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFisher, "BoostedFisher", 
                            "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" )

    # Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    if "FDA_MC" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_MC",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

    if "FDA_GA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_GA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

    if "FDA_SA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_SA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

    if "FDA_MT" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_MT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

    if "FDA_GAMT" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_GAMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

    if "FDA_MCMT" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_MCMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

    # TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
    if "MLP" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" )

    if "MLPBFGS" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" )

    if "MLPBNN" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ) # BFGS training with bayesian regulators

    # CF(Clermont-Ferrand)ANN
    if "CFMlpANN" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ) # n_cycles:#nodes:#nodes:...  

    # Tmlp(Root)ANN
    if "TMlpANN" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ) # n_cycles:#nodes:#nodes:...

    # Support Vector Machine
    if "SVM" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" )

    # Boosted Decision Trees
    if "BDTG" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDTG",
                            "!H:!V:NTrees=1000:MinNodeSize=1.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" )                        

    if "BDT" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" )

    if "BDTB" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" )

    if "BDTD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" )

    # RuleFit -- TMVA implementation of Friedman's method
    if "RuleFit" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kRuleFit, "RuleFit",
                            "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" )

    # --------------------------------------------------------------------------------------------------
            
    # ---- Now you can tell the factory to train, test, and evaluate the MVAs. 

    # Train MVAs
    factory.TrainAllMethods()
    
    # Test MVAs
    factory.TestAllMethods()
    
    # Evaluate MVAs
    factory.EvaluateAllMethods()    
    
    # Save the output.
    outputFile.Close()
    
    print "=== wrote root file %s\n" % OutFileName
    print "=== TMVAClassification is done!\n"
    
    # open the GUI for the result macros    
    #gROOT.ProcessLine( "TMVAGui(\"%s\")" % OutFileName )
    
    # keep the ROOT thread running
    #gApplication.Run() 

# ----------------------------------------------------------

if __name__ == "__main__":
    main()
