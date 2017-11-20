import os
import sys
from ROOT import *
from AtlasStyle import *
from numpy import log,arange
SetAtlasStyle();
gStyle.SetPalette(1)
TH1.SetDefaultSumw2()


def MakeNewDir(outputdir):
    time = TDatime();
    date = time.GetDate();
    outputdir+=str(date);
    outputdir+="/";
    os.system("mkdir -p "+outputdir);
    
    return outputdir

def GetHist1D(file, ntupname, var, rangex, weight):
    print "Getting 1D Hist: ",file, ntupname, var, rangex, weight
    
    ftemp    = TFile(file)
    ntuptemp = ftemp.Get(ntupname)
    ntuptemp.Draw(var+">>htemp("+rangex+")",weight,"e")
    htemp = gDirectory.Get("htemp")
    htemp.SetDirectory(0)

    ftemp.Close()
    return htemp

def GetHist2D(file, ntupname, varx, vary, rangex, rangey, weight):
    print "Getting 2D Hist: ",file, ntupname, varx, vary, rangex, rangey, weight
    ftemp    = TFile(file)
    ntuptemp = ftemp.Get(ntupname)

    varstring = vary+":"+varx+">>htemp("+rangex+","+rangey+")"

    ntuptemp.Draw(varstring,weight,"e")
    htemp = gDirectory.Get("htemp")

    htemp.SetDirectory(0)

    ftemp.Close()
    return htemp

def NormalizeHist( hist, norm=1.0):

    tot = hist.Integral()

    if tot<=0:
        print "Something is wrong with histogram"
        return hist
    else:
        hist.Scale(norm/tot)
        return hist

def GetMaxVal(hists):
    max = 0
    
    for hist in hists:
        if hist.GetBinContent(hist.GetMaximumBin())>max:
            max=hist.GetBinContent(hist.GetMaximumBin())

    print "MaxVal = ",max

    return max

def TranslateRegion(pt1,pt2,m1,m2):
    print "TranslatingRegion:",pt1,pt2,m1,m2
    region  = "|#eta|<1.2 , p_{T}=["+pt1+","+pt2+"] GeV , m=["+m1+","+m2+"] GeV"
    return region  
    
def TranslateAlg(alg):
    print "Translating alg: ",alg
    output = ""
    if   alg=="TruthRawTrim":    output = "Trimmed truth anti-k_{t} R = 1.0"    
    elif alg=="CaloTrim":     output = "Reco anti-k_{t} R = 1.0"   
    return output

def TranslateVar(var):
    #print "Translating var: ",var

    AllVars={}
    AllVars["m"]                         = "M [GeV]"
    AllVars["Tau21"]                     = "#tau_{21}"
    AllVars["Tau32"]                     = "#tau_{32}"
    AllVars["T2jet"]                     = "v_{2}"
    AllVars["T3jet"]                     = "v_{3}"
    AllVars["T2jet_angle"]               = "#theta_{2}"
    AllVars["T3jet_angle"]               = "#theta_{min}"
    AllVars["tau21_WTA"]                 = "#tau_{21}"
    AllVars["tnsub_beta10to20_tau21_vol"]= "T(#tau_{2}^{wta}, #beta=[1.0,2.0])"
    AllVars["C2"]= "C_{2}"
    AllVars["D2"]= "D_{2}"
    AllVars["v32"] = "v_{32}"
    AllVars["v42"] = "v_{42}"
    AllVars["v52"] = "v_{52}"
    AllVars["v21"] = "v_{21}"
    AllVars["Wmass"] = "M_{W}"
    AllVars["WmassVolatility"] = "v_{W}"
    AllVars["Ttrimming"]                 = "v_{trim}"
    AllVars["Tpruning"]                  = "v_{prune}"
    
    if var in AllVars.keys():
        varout=AllVars[var]
    else:
        print "This var is not in your list: ",var
        varout=var

    return varout

def GetRange(var):
    if var == "T2jet": return 0.5
    elif var == "T3jet": return 0.2
    elif var == "v32": return 1 
    else: return 0

def MakeReferenceGraph( roctype ):
    '''Make reference graph for plotting'''

    graph = TGraph(1)

    graph.SetTitle("")
    graph.GetXaxis().SetTitleOffset(1.2)
    graph.GetYaxis().SetTitleOffset(1.3)

    graph.GetXaxis().SetRangeUser(0.0,1.0)

    if roctype==0:
        graph.GetXaxis().SetTitle("#epsilon ^{FullTag}_{W Jets}")
        graph.GetYaxis().SetTitle("1- #epsilon ^{FullTag}_{QCD Jets}")
        graph.SetMinimum(0.0)
        graph.SetMaximum(1.0)
    elif roctype==1:
        graph.GetXaxis().SetTitle(Form("#font[132]{Signal Efficiency}"))
        graph.GetYaxis().SetTitle(Form("#font[132]{Background Efficiency}"))
        graph.SetMinimum(0.0)
        graph.SetMaximum(1.0)
    elif roctype==2:
        graph.GetXaxis().SetTitle("#epsilon ^{FullTag}_{W Jets}")
        graph.GetYaxis().SetTitle("1 / #epsilon ^{FullTag}_{QCD Jets}")
        graph.SetMinimum(0.0)
        graph.SetMaximum(1000.0)

    return graph

def ConvertROC_BGRej_BGPow_WithUncer( gr ):

    grTemp = gr
    n = grTemp.GetN()
    
    print "Type: ",type(gr)
    
    print "Converting NPoints ",n

    xinit = Double()
    yinit = Double()
    mindistance=10.0
    foundseff=1.0
    foundbgeff=1.0
    itarget=1

    for i in range(n):
    
        sigeff_init = Double()
        bkgrej_init = Double()
        grTemp.GetPoint(i,sigeff_init,bkgrej_init)
        
        sigefferr_init = Double()
        bkgrejerr_init = Double()
        sigefferr_init = gr.GetErrorX(i)
        bkgrejerr_init = gr.GetErrorY(i)


        #no conversion for signal efficiency
        sigeff_final    = sigeff_init
        sigefferr_final = sigefferr_init

        #convert bg rejection
        if bkgrej_init >0.999:
            bkgrej_final  = 1000.0
            bkgrejerr_final = 100.0
        else:
            h_one    = TH1D("h_one","h_one",1,0,1)
            h_one.SetBinContent(1,1.0)
            h_one.SetBinError(1,0.0)
    
            h_bkgrej = TH1D("h_bkgrej","h_bkgrej",1,0,1)
            h_bkgrej.SetBinContent(1,bkgrej_init)
            h_bkgrej.SetBinError(1,bkgrejerr_init)            
        
            h_bkgrej.Add(h_bkgrej,h_one,-1.0,1.0)
            h_one.Divide(h_one,h_bkgrej)
        
            bkgrej_final    = h_one.GetBinContent(1)
            bkgrejerr_final = h_one.GetBinError(1)

        print "Power converted: ",i,sigeff_final,sigefferr_final,bkgrej_final,bkgrejerr_final
        gr.SetPoint(i,sigeff_final,bkgrej_final)
        gr.SetPointError(i,sigefferr_final,bkgrejerr_final)

    return gr

def FoldMassEffGraphs_WithUncer( gr, masssigeff, masssigefferr, massbkgeff, massbkgefferr):
    #print "Fold mass effs",masssigeff,massbkgeff
    n = gr.GetN();
    grOut = TGraphErrors(n)

    #for calculations    
    h_one    = TH1D("h_one","h_one",1,0,1)
    h_one.SetBinContent(1,1.0)
    h_one.SetBinError(1,0.0)
    
    h_masssigeff = TH1D("h_masssigeff","h_masssigeff",1,0,1)
    h_masssigeff.SetBinContent(1,masssigeff)
    h_masssigeff.SetBinError(1,masssigefferr)
    
    h_massbkgeff = TH1D("h_massbkgeff","h_massbkgeff",1,0,1)
    h_massbkgeff.SetBinContent(1,massbkgeff)
    h_massbkgeff.SetBinError(1,massbkgefferr)
    
    h_sig = TH1D("h_sig","h_sig",1,0,1)
    h_bkg = TH1D("h_bkg","h_bkg",1,0,1)
    
    for i in range(n):

        sigeffinit = Double()
        bkgeffinit = Double()
        gr.GetPoint(i,sigeffinit,bkgeffinit)
        
        sigeffiniterr = Double()
        bkgeffiniterr = Double()
        sigeffiniterr = gr.GetErrorX(i)
        bkgeffiniterr = gr.GetErrorY(i)

        #print "Init: ",i,"  ",sigeffinit,"  ",sigeffiniterr,"  ",bkgeffinit,"  ",bkgeffiniterr

        h_sig.SetBinContent(1,sigeffinit)
        h_sig.SetBinError(1,sigeffiniterr)
        h_sig.Multiply(h_masssigeff)
        sigefffinal    = h_sig.GetBinContent(1)
        sigefffinalerr = h_sig.GetBinError(1)
        #print "FromHist: ",sigefffinal,sigefffinalerr


        h_bkg.SetBinContent(1,bkgeffinit)
        h_bkg.SetBinError(1,bkgeffiniterr)
        h_bkg.Multiply(h_massbkgeff)
        bkgefffinal    = h_bkg.GetBinContent(1)
        bkgefffinalerr = h_bkg.GetBinError(1)
        #print "FromHist: ",bkgefffinal,bkgefffinalerr
        
        

        #print "Final: ",i,"  ",sigefffinal,"  ",bkgrejfinal

        grOut.SetPoint(i,sigefffinal,bkgefffinal)
        grOut.SetPointError(i,sigefffinalerr,bkgefffinalerr)

#         grOut.SetMinimum(0.0)
#         grOut.SetMaximum(0.0)
#         grOut.GetXaxis().SetRangeUser(0.0,1.0)
#         grOut.SetLineColor(2)
#         grOut.SetFillColor(2) 
#         grOut.Draw("CE3Same")

    return grOut

def RocCurve_SoverBOrdered_WithUncer(sig, bg ):
    print 'Make ROC curve using S over B ordering'
    
    binnum=[]
    s=[]
    serr=[]
    b=[]
    berr=[]
    r=[]

    for i in range(1,sig.GetNbinsX()+1):
        #print "bin ",i

        binnumtemp = i
        stemp = sig.GetBinContent(i)
        serrtemp = sig.GetBinError(i)
        btemp = bg.GetBinContent(i)
        berrtemp = bg.GetBinError(i)

        if btemp==0:
            rtemp=0.0
        else:
            rtemp = stemp/btemp

        binnum.append(binnumtemp)
        s.append(stemp)
        serr.append(serrtemp)
        b.append(btemp)
        berr.append(berrtemp)
        r.append(rtemp)

    for i in range(len(s)):
        ifix = len(s)-i
        #print ifix
        for j in range(0,ifix-1):
            if r[j]<r[j+1]:

                binnum[j],binnum[j+1]=binnum[j+1],binnum[j]
                b[j],b[j+1]=b[j+1],b[j]
                berr[j],berr[j+1]=berr[j+1],berr[j]
                s[j],s[j+1]=s[j+1],s[j]
                serr[j],serr[j+1]=serr[j+1],serr[j]
                r[j],r[j+1]=r[j+1],r[j]


    #make reordered histograms
    n = len(s)
    print n
    
    hsigout = TH1F("hsigout","hsigout",n,0,1)
    hsigout.SetDirectory(0)
    hbkgout = TH1F("hbkgout","hbkgout",n,0,1)
    hbkgout.SetDirectory(0)
    for i in range(0,n):
        hsigout.SetBinContent(i+1,s[i])
        hsigout.SetBinError(i+1,serr[i])
        hbkgout.SetBinContent(i+1,b[i])
        hbkgout.SetBinError(i+1,berr[i])


    
    siglow  = sig.GetXaxis().GetXmin()
    sighigh = sig.GetXaxis().GetXmax()
    h1 = TH1F("h1","h1",n,siglow,sighigh)
    h1.SetDirectory(0)
    
    print "Npoints: ",n
    gr = TGraphErrors(n)
    for i in range(1,n):

        #integrate from 0 to i
        myBerr=Double()
        mySerr=Double()
        myB = hbkgout.IntegralAndError(0,i,myBerr)
        myS = hsigout.IntegralAndError(0,i,mySerr)
        print i,"  myS=",myS,mySerr,"  myB=",myB,myBerr
        gr.SetPoint(i, myS, myB)
        gr.SetPointError(i, mySerr, myBerr)

        #get histograms that are colored for signal efficiency at 50%
        if myS<=0.73:
            print binnum[i]
            print s[i]
            h1.SetBinContent(binnum[i], s[i])
        
    return gr,hsigout,hbkgout,h1
  
def RocCurve_SingleSided_WithUncer(sig, bkg, rightleft):
    print "\n\nMake ROC curve using right/left cut",rightleft

    n = bkg.GetNbinsX()
    print "NBins",n

    totalBerr=Double()
    totalSerr=Double()
    totalB = bkg.IntegralAndError(0,n,totalBerr)
    totalS = sig.IntegralAndError(0,n,totalSerr)
    
    siglow  = sig.GetXaxis().GetXmin()
    sighigh = sig.GetXaxis().GetXmax()
    hsigreg50 = TH1F("hsigreg50","hsigreg50",n,siglow,sighigh)
    hsigreg50.SetDirectory(0)
    hcutval50 = TH1F("hcutval50","hcutval50",5,0,5)
    hcutval50.SetDirectory(0)
    hcutval50.GetXaxis().SetBinLabel(1,"Left(0) , Right(1)")
    hcutval50.GetXaxis().SetBinLabel(2,"LowerCut")
    hcutval50.GetXaxis().SetBinLabel(3,"UpperCut")
    hsigreg25 = TH1F("hsigreg25","hsigreg25",n,siglow,sighigh)
    hsigreg25.SetDirectory(0)
    hcutval25 = TH1F("hcutval25","hcutval25",5,0,5)
    hcutval25.SetDirectory(0)
    hcutval25.GetXaxis().SetBinLabel(1,"Left(0) , Right(1)")
    hcutval25.GetXaxis().SetBinLabel(2,"LowerCut")
    hcutval25.GetXaxis().SetBinLabel(3,"UpperCut")
    if rightleft=="R":
        hcutval50.SetBinContent(1,1)
        hcutval50.SetBinContent(3,sig.GetXaxis().GetBinLowEdge(n)+sig.GetXaxis().GetBinWidth(n))
        extrema50 = 100000
        hcutval25.SetBinContent(1,1)
        hcutval25.SetBinContent(3,sig.GetXaxis().GetBinLowEdge(n)+sig.GetXaxis().GetBinWidth(n))
        extrema25 = 100000
    elif rightleft=="L":
        hcutval50.SetBinContent(1,0)
        hcutval50.SetBinContent(2,sig.GetXaxis().GetBinLowEdge(1))
        extrema50 = -100000
        hcutval25.SetBinContent(1,0)
        hcutval25.SetBinContent(2,sig.GetXaxis().GetBinLowEdge(1))
        extrema25 = -100000

    gr = TGraphErrors(n-1)
    for i in range(1,n):
        myS = 0.
        myB = 0.

        if rightleft=="R":
            #loop grom i to end
            myBerr=Double()
            mySerr=Double()
            myB = bkg.IntegralAndError(i,n,myBerr)
            myS = sig.IntegralAndError(i,n,mySerr)
            print i,"  myS=",myS,"  myB=",myB
            gr.SetPoint(i, myS, myB)
            gr.SetPointError(i, mySerr, myBerr)
            if myS<=0.73:
                hsigreg50.SetBinContent(i, sig.GetBinContent(i))
                tempex=sig.GetXaxis().GetBinLowEdge(i)
                print tempex,extrema50
                if tempex<extrema50:
                    extrema50 = tempex
                    print "found extrema R: ",extrema50
            if myS<=0.36:
                hsigreg25.SetBinContent(i, sig.GetBinContent(i))
                tempex=sig.GetXaxis().GetBinLowEdge(i)
                print tempex,extrema25
                if tempex<extrema25:
                    extrema25 = tempex
                    print "found extrema R: ",extrema50

            #artificially set the first point to (1,1) to avoid overflow issues    
            gr.SetPoint(0, 1.0, 1.0)
            gr.SetPointError(0, 0.00001, 0.00001)            

        elif rightleft=="L":
            #loop grom 0 to i
            myBerr=Double()
            mySerr=Double()
            myB = bkg.IntegralAndError(1,i,myBerr)
            myS = sig.IntegralAndError(1,i,mySerr)
            print i,"  myS=",myS,"  myB=",myB
            gr.SetPoint(i, myS, myB)
            gr.SetPointError(i, mySerr, myBerr)
            if myS<=0.73:
                hsigreg50.SetBinContent(i, sig.GetBinContent(i))
                tempex=sig.GetXaxis().GetBinLowEdge(i)+sig.GetXaxis().GetBinWidth(i)
                print tempex,extrema50
                if tempex>extrema50:
                    extrema50 = tempex
                    print "found extrema L: ",extrema50
            if myS<=0.36:
                hsigreg25.SetBinContent(i, sig.GetBinContent(i))
                tempex=sig.GetXaxis().GetBinLowEdge(i)+sig.GetXaxis().GetBinWidth(i)
                print tempex,extrema25
                if tempex>extrema25:
                    extrema25 = tempex
                    print "found extrema L: ",extrema50
            
            #artificially set the first point to (1,1) to avoid overflow issues    
            gr.SetPoint(0, 0.0, 0.0)
            gr.SetPointError(0, 0.00001, 0.00001)            
        else:
            print "You did not choose a left or right handed cut - EXITTING ..."
            sys.exit()
            

            
    ctest = TCanvas("ctest","ctest",400,400)
    gr.SetMinimum(-5.0)
    gr.SetMaximum(5.0)
    gr.GetXaxis().SetRangeUser(-5.0,5.0)
    gr.Draw("AE3")
    
    if rightleft=="R":
        hcutval50.SetBinContent(2,extrema50)
        hcutval25.SetBinContent(2,extrema25)
    elif rightleft=="L":
        hcutval50.SetBinContent(3,extrema50)
        hcutval25.SetBinContent(3,extrema25)


    print "RETURNING"
    return gr,hsigreg50,hcutval50,hsigreg25,hcutval25

def GetBGRej( gr, targetseff ):

    n = gr.GetN();

    xinit       = Double()
    yinit       = Double()
    mindistance = 10.0
    foundseff   = 1.0
    foundbgeff=1.0
    itarget=1

    for i in range(n):
        gr.GetPoint(i,xinit,yinit)
        distance = abs(targetseff-xinit)
        #print "FindBR: ",i," ",xinit," ",yinit," ",distance,"   -   ",foundseff,"  ",foundbgeff


        #print distance,"   -  ",mindistance
        if distance<mindistance:
            #print "SWITCH"
            mindistance = distance
            #Tricky - need to cast them as floats or pyroot treats them as pointer references
            foundseff   = float(xinit)
            foundbgeff  = float(yinit)
            itarget     = i



    #print itarget,"  ",mindistance,"  ",foundseff,"  ",foundbgeff

    return foundseff,foundbgeff



