import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)

def InterpolateGraph(gr, targetseff):
    #print "Interpolating graph ...",targetseff
    
    n = gr.GetN()
    
    tempsigeff = Double()
    tempbkgrej = Double()
    tempsigefferr = Double()
    tempbkgrejerr = Double()
    sigeff = 1.0
    bkgeff = 1.0
    mindistance=10.0
    foundsigeff=1.0
    foundsigefferr=1.0
    foundbkgeff=1.0
    foundbkgefferr=1.0
    itarget=1

    for ipoint in range(n):

        gr.GetPoint(ipoint,tempsigeff,tempbkgrej)
        tempsigefferr = gr.GetErrorX(ipoint)
        tempbkgrejerr = gr.GetErrorY(ipoint)
        
        #print "Int: ",ipoint,tempsigeff,tempbkgrej

        distance = abs(targetseff-tempsigeff)
        
        #print distance,mindistance

        if distance<mindistance:
        #    print "switching"
            itarget         = ipoint
            mindistance     = distance
            foundsigeff     = float(tempsigeff)
            foundsigefferr  = float(tempsigefferr)
            foundbkgrej     = float(tempbkgrej)
            foundbkgrejerr  = float(tempbkgrejerr)
            
           
#    print "DONE: ",foundsigeff,foundsigefferr,foundbkgrej,foundbkgrejerr
    return foundsigeff,foundsigefferr,foundbkgrej,foundbkgrejerr
    
def ReadMassEffs(alg, pt1, pt2, inputsinglevariable):

    f=TFile(inputsinglevariable+"MassEff_"+alg+"_pt"+pt1+pt2+".root")
    h=f.Get("masseff")
    
    effsig=h.GetBinContent(1)
    effsig_err=h.GetBinContent(2)
    effbkg=h.GetBinContent(3)
    effbkg_err=h.GetBinContent(4)
    
    return effsig,effsig_err,effbkg,effbkg_err
    
    

############################################################A
############################################################A
############################################################A
############################################################A
############################################################A

#==========================
#Set output directory name
#==========================
InputDir="Data/20150709/"

outputdir1 = "OutputSingleMatrix/"
outputdir1 = MakeNewDir(outputdir1)

inputsinglevariable="OutputSingleVariable/20150710/"


#ALGORITHMS
algs=[]
algs.append("truth")
algs.append("reco")
algs.append("fullsim")

# VARIABLES AND RANGES
VarsAndRanges={}
VarsAndRanges["tau21_WTA"]=                        [0, "100,0,1", "100,0,1" ,"L"]
VarsAndRanges["tnsub_beta10to20_tau21_vol"]=       [0, "100,0,1", "100,0,1" ,"R"]
VarsAndRanges["D2"]=                               [0, "100,0,5", "100,0,5" ,"L"]
VarsAndRanges["C2"]=                               [0, "100,0,1", "100,0,1" ,"L"]


ForCutTable = 0



    
wp_number=0.50
wp_label="50"


CutRegions=[]
CutRegions.append("1")
#CutRegions.append("2")

ROCTypes=[]
#ROCTypes.append("SBOrdered")
ROCTypes.append("LeftRight")

dictCuts={}

for CutRegion in CutRegions:

    dictCuts.clear()
    for ROCType in ROCTypes:

        DictMatrix={}
        DictMatrixErr={}
        DictMatrixCutVal={}
        DictMatrixCutDir={}
        DictMatrixMassLow={}
        DictMatrixMassHigh={}
        DictMaxLoc={}
            
        for alg in algs:
    
            if CutRegion=="1" or CutRegion=="4": 
                pt1="350"; pt2="500";  m1="60"; m2="100"; 
            if CutRegion=="2" or CutRegion=="5": 
                pt1="800"; pt2="1000";  m1="60"; m2="100"; 
            
            masssigeff,masssigefferr,massbkgeff,massbkgefferr = ReadMassEffs(alg, pt1, pt2, inputsinglevariable)
            print "MassEff: ",masssigeff,masssigefferr,massbkgeff,massbkgefferr


            
            ###################################
            #Rank Vars
            vars=[]
            bgrejs=[]
            bgrejserr=[]
            bgcutval=[]
            bgcutdir=[]
            masslow=[]
            masshigh=[]
            
            sortedVars=sorted(VarsAndRanges.keys())
            for var in sortedVars:
                #print alg,var
            
                checkpath=inputsinglevariable+"ROC_"+alg+"_"+var+"_pt"+pt1+pt2+".root"
                if not os.path.isfile(checkpath):
                    print "File not existent"
        
                f = TFile(checkpath)
                f.ls()
                #print "LeftRight roc"
                gr = f.Get("ROC_"+VarsAndRanges[var][3])

                #grfolded = FoldMassEffGraphs( gr, masssigeff, massbkgeff)
                grfolded = FoldMassEffGraphs_WithUncer( gr, masssigeff, masssigefferr, massbkgeff, massbkgefferr)
                
                #Get the error by finding the closest graph point close to 50%
                sigInt,sigIntErr,bgInt,bgIntErr = InterpolateGraph(grfolded, wp_number)
                
                #print "Precise: ",sigInt,sigIntErr,bgInt,bgIntErr
                
                #Get the value using the interpolation
                bgeffvar = grfolded.Eval(wp_number)
                bgeffvarerr = bgIntErr
                
                #print "Interpolated: ",bgeffvar,bgeffvarerr
                    
                #print alg,var,bgeffvar,bgeffvarerr
                vars     .append(var)
                bgrejs   .append(1.0/bgeffvar)
                bgrejserr.append(bgeffvarerr/(bgeffvar)**2) ##HACK 
        
            print "\nPRINTING REJECTION:"
            for var,bgrej,bgrejerr in zip(vars,bgrejs,bgrejserr):
                print var,bgrej,bgrejerr
    

            AllVarsSorted = []
            DictAllVarsSorted = {}
            DictAllVarsSortedErr = {}
            for x in sorted(zip(bgrejs,bgrejserr,vars)):
                AllVarsSorted.append([x[0],x[1],x[2]])
                DictAllVarsSorted[x[2]]   =x[0]   
                DictAllVarsSortedErr[x[2]]=x[1]   

            tempalg={}
            tempalgerr={}
            for key in vars:
                #print key,DictAllVarsSorted[key],DictAllVarsSortedErr[key]
                tempalg[key]=DictAllVarsSorted[key]
                tempalgerr[key]=DictAllVarsSortedErr[key]

            DictMatrix[alg]=tempalg
            DictMatrixErr[alg]=tempalgerr

        #print "\n\nMAKING MATRIX: ",cutslabel
        nalgs=len(algs)
        nvars=len(vars)
    
        hmatrix = TH2F("hmatrix","hmatrix",nvars,0,nvars,nalgs,0,nalgs)
        algcount=0
        
        DictOptList={}
        
        for alg in algs:
            algcount+=1
            #print "\n\nALGORITHM: ",alg
        
            algVars = sorted(DictMatrix[alg].keys())
            
            varcount=0
            for key in algVars:
                varcount+=1
                hmatrix.SetBinContent(varcount,algcount,round(DictMatrix[alg][key],1))
                hmatrix.SetBinError(varcount,algcount,round(DictMatrixErr[alg][key],2))
            
                hmatrix.GetYaxis().SetBinLabel(algcount,TranslateAlg(alg))
                hmatrix.GetXaxis().SetBinLabel(varcount,TranslateVar(key))
            
            
        c = TCanvas("c","c",2000,1500)
        c.SetRightMargin(0.14)
        c.SetBottomMargin(0.15)
        c.SetLeftMargin(0.25)
        c.SetTopMargin(0.18)

        hmatrix.GetZaxis().SetTitle("QCD Rejection @ #epsilon_{ W}^{ G & T} = 50%")
        hmatrix.GetZaxis().SetTitleSize(0.04)
        hmatrix.GetZaxis().SetTitleOffset(1.0)
        hmatrix.GetZaxis().SetLabelSize(0.03)

#        hmatrix.GetXaxis().LabelsOption("v")
        hmatrix.GetXaxis().SetLabelSize(0.04)
    
        hmatrix.GetYaxis().SetLabelSize(0.03)
    
        hmatrix.SetMarkerSize(1.5)

        hmatrix.Draw("colz text E")


        ATLASLabel(   0.25,0.95,1,0.07,0.05,"#sqrt{s}=13 TeV")
        myText(       0.50,0.95,1,0.05, "Summary(W-jet vs. QCD-jet)")
        myText(       0.25,0.87,1,0.05, TranslateRegion(pt1,pt2,m1,m2))
        c.SaveAs(outputdir1+"RankMatrix_pt"+pt1+pt2+".eps")
