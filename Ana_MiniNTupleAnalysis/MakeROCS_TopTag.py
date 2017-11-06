import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)

sigFile="ntuple_ttLo.root"
bkgFile="ntuple_dijetLo.root"
#sigFile="ntuple_ttHi.v2.root"
#bkgFile="ntuple_dijetHi.v2.root"

def SignalBGCompare1D(InputDir, alg, variable, range, logy, pt1, pt2, m1, m2, outputdir):
    '''Implementation of simple signal and background comparison'''
    print "Making 1D comparison: ",alg,variable
    
    c = TCanvas("c","c",300,300)
    c.SetBottomMargin(0.15);
    c.SetLeftMargin(0.15);
    c.SetTopMargin(0.05);
    c.SetRightMargin(0.05);

    dry=0.045

    weight=""
    weight+="("
    weight+=alg+"_pt>"+pt1+" && "
    weight+=alg+"_pt<"+pt2
    if m1!="0":
        weight+=" && "+alg+"_m>"+m1+" && "
        weight+=alg+"_m<"+m2
    weight+=")"
    
    #Get signal and background histograms
    histname = alg+"_"+variable
    if variable == "v32":
        histname = alg+"_T3jet / "+alg+"_T2jet"
    if variable == "v21":
        histname = alg+"_T2jet / "+alg+"_T1jet"

    hsig = GetHist1D(InputDir+sigFile, "JetTree", histname, range, weight)#+"*("+alg+"_flavor==3)")
    hbkg = GetHist1D(InputDir+bkgFile, "JetTree", histname, range, weight)#+"*("+alg+"_flavor==0)")

    hsig.Rebin()
    hbkg.Rebin()

    #Normalize them to unity
    hsig = NormalizeHist(hsig)
    hbkg = NormalizeHist(hbkg)
    
    #==============================
    #Make all rocs
    rocL,hsigregL,hcutvalL,hsigregL25,hcutvalL25 = RocCurve_SingleSided_WithUncer(hsig, hbkg, "L")
    rocR,hsigregR,hcutvalR,hsigregR25,hcutvalR25 = RocCurve_SingleSided_WithUncer(hsig, hbkg, "R")
    rocSB,sigordered,bkgordered,h1 = RocCurve_SoverBOrdered_WithUncer(hsig, hbkg)

    print alg+"_"+variable
    f = TFile(outputdir+"Top_ROC_"+alg+"_"+variable.replace("\\","")+"_pt"+pt1+pt2+".root","RECREATE")
    rocL.Write("ROC_L")
    rocR.Write("ROC_R")
    rocSB.Write("ROC_SoverB")
    f.Close()
    
    c.cd()
    rocSB.Draw("AC*")
    c.SaveAs(outputdir+"Top_"+"ROCDraw_"+alg+"_"+variable.replace("\\","")+"_pt"+pt1+pt2+".eps")
    
    #################################

    #copy for line drawing
    #hbkg.Rebin(4)
    #sh.Rebin(4)
    hbkgline = hbkg.Clone("hbkgline")
    hsigline = hsig.Clone("hsigline")

    hbkgline.SetFillStyle(0)
    hbkgline.SetLineColor(2)
    hbkgline.SetLineStyle(1)
    hbkgline.SetLineWidth(3)
    hsigline.SetFillStyle(0)
    hsigline.SetLineColor(4)
    hsigline.SetLineStyle(1)
    hsigline.SetLineWidth(3)

    hbkg.SetLineColor(2)
    hbkg.SetLineWidth(4)
    hbkg.SetLineStyle(1)
    hbkg.SetFillColor(2)
    hbkg.SetFillStyle(3002)
    hbkg.SetMarkerSize(0)

    hsig.SetLineColor(4)
    hsig.SetLineWidth(4)
    hsig.SetLineStyle(1)
    hsig.SetFillColor(4)
    hsig.SetFillStyle(3002)
    hsig.SetMarkerSize(0)

    # DRAW
    print variable
    maxval = GetMaxVal([hbkg, hsig])
    hbkg.SetMaximum(maxval*1.5)
    hbkg.SetMinimum(0.001)

    hbkg.Draw("E2")

    hbkg.GetXaxis().SetTitle(TranslateVar(variable))
    hbkg.GetXaxis().SetTitleOffset(1.2)
    hbkg.GetYaxis().SetTitleOffset(1.6)
    hbkg.GetYaxis().SetTitle("Normalised Entries")
    #hbkg.GetYaxis().SetTitle("")
    hbkg.GetXaxis().SetLabelSize(0.04)
    hbkg.GetYaxis().SetLabelSize(0.04)
    axrange = GetRange(variable)
    if axrange != 0:
        hbkg.GetXaxis().SetRangeUser(0,axrange)

    hbkgline.Draw("samehist")
    hsig.Draw("E2same")
    hsigline.Draw("samehist")

    ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
    myText(       0.20,0.83,1,0.03, TranslateAlg(alg))
    myText(       0.20,0.78,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
    rocbox3=myLineBoxText(0.70, 0.88, 4, 1, 1, 0, 0.1, 0.08, "Top Jets")
    rocbox4=myLineBoxText(0.70, 0.83, 2, 1, 1, 0, 0.1, 0.08, "QCD Jets")
    c.SetLogy(logy)
    c.SaveAs(outputdir+"Top_"+"SignalBGCompare_"+alg+"_"+variable+"_pt"+pt1+pt2+".eps")

def GetMassEffs(InputDir, alg, m1, m2, outputdir):

    weight=""
    weight+="("
    weight+=alg+"_pt>"+pt1+" && "
    weight+=alg+"_pt<"+pt2+")"
    
    range="100,0,200"
    
    #Get signal and background histograms
    histname=alg+"_m"
    print histname
    hsig = GetHist1D(InputDir+sigFile,    "JetTree", histname, range, weight+"*("+alg+"_flavor==3)")
    hbkg = GetHist1D(InputDir+bkgFile, "JetTree", histname, range, weight+"*("+alg+"_flavor==0)")

    #Normalize them to unity
    hsig = NormalizeHist(hsig)
    hbkg = NormalizeHist(hbkg)
    
    #Integrate in window to get efficiencies
    bin1=hsig.FindBin(float(m1))
    bin2=hsig.FindBin(float(m2))

    effsig_err = Double()
    effsig = hsig.IntegralAndError(bin1,bin2,effsig_err)
    effbkg_err = Double()
    effbkg = hbkg.IntegralAndError(bin1,bin2,effbkg_err)
    
    h=TH1D()
    h.Fill("effsig",effsig)
    h.Fill("effsig_err",effsig_err)
    h.Fill("effbkg",effbkg)
    h.Fill("effbkg_err",effbkg_err)
    fout = TFile(outputdir+"Top_MassEff_"+alg+"_pt"+pt1+pt2+".root","RECREATE")
    h.Write("masseff")
    fout.Close()

    return effsig,effsig_err,effbkg,effbkg_err

def Make2DROC(alg, varX, varXcutdir, varY, varYcutdir, sig, bkg, rankMetric, statErrorOnRatioThreshold, flagFilter, outputdir):

    print alg, varX, varXcutdir, varY, varYcutdir, sig, bkg, rankMetric, statErrorOnRatioThreshold, flagFilter, outputdir


    #output label flag that labels if stats filter is applied
    flagFilterLabel="noFilter"
    if flagFilter:
        flagFilterLabel="Filter"  

    #Rebin if need be
    rebinx=2
    rebiny=2
    sig.RebinX(rebinx)
    sig.RebinY(rebiny)
    bkg.RebinX(rebinx)
    bkg.RebinY(rebiny)
    
    corrsig = sig.GetCorrelationFactor()
    corrbg  = bkg.GetCorrelationFactor()
    
    #1D Projections and cration of ROC curves for comparison at end
    sig1DvarX = sig.ProjectionX("sig1DvarX",0,-1,"e")
    bkg1DvarX = bkg.ProjectionX("bkg1DvarX",0,-1,"e")
    sig1DvarY = sig.ProjectionY("sig1DvarY",0,-1,"e")
    bkg1DvarY = bkg.ProjectionY("bkg1DvarY",0,-1,"e")
    
    sig1DvarX = NormalizeHist(sig1DvarX)
    bkg1DvarX = NormalizeHist(bkg1DvarX)
    sig1DvarY = NormalizeHist(sig1DvarY)
    bkg1DvarY = NormalizeHist(bkg1DvarY)
    
    c1d = TCanvas("c1d","c1d",400,400)   
    sig1DvarX.SetLineColor(4)
    sig1DvarX.SetLineWidth(2)
    bkg1DvarX.SetLineColor(2)
    bkg1DvarX.SetLineWidth(2)
    sig1DvarX.GetXaxis().SetTitle(TranslateVar(varX))
    sig1DvarX.GetYaxis().SetTitle("Normalized Units")
    sig1DvarX.Draw("hist")
    bkg1DvarX.Draw("histsame")
    ATLASLabel(   0.70,0.90,1,0.1,0.03,"#sqrt{s}=13 TeV")
    bx1=myLineBoxText(0.70, 0.80, 2, 1, 2, 3004, 0.1, 0.1, "Signal")
    bx2=myLineBoxText(0.70, 0.75, 4, 1, 4, 3005, 0.1, 0.1, "Background")
    #c1d.SaveAs(outputdir+"Make2DROC_"+varX+"_pt"+pt1+pt2+"_"+flagFilterLabel+".eps")

    sig1DvarY.SetLineColor(4)
    sig1DvarY.SetLineWidth(2)
    bkg1DvarY.SetLineColor(2)
    bkg1DvarY.SetLineWidth(2)
    sig1DvarY.GetXaxis().SetTitle(TranslateVar(varY))
    sig1DvarY.GetYaxis().SetTitle("Normalized Units")
    sig1DvarY.Draw("hist")
    bkg1DvarY.Draw("histsame")
    ATLASLabel(   0.70,0.90,1,0.1,0.03,"#sqrt{s}=13 TeV")
    bx3=myLineBoxText(0.70, 0.80, 2, 1, 2, 3004, 0.1, 0.1, "Signal")
    bx4=myLineBoxText(0.70, 0.75, 4, 1, 4, 3005, 0.1, 0.1, "Background")
    #c1d.SaveAs(outputdir+"Make2DROC_"+varY+"_pt"+pt1+pt2+"_"+flagFilterLabel+".eps")
    
    print "Single sided rocs"
    rocvarX,hsigregX,hcutvalX,hsigregX25,hcutvalX25 = RocCurve_SingleSided_WithUncer(sig1DvarX, bkg1DvarX, varXcutdir)
    rocvarX.SetLineColor(2)
    rocvarY,hsigregY,hcutvalY,hsigregX25,hcutvalX25 = RocCurve_SingleSided_WithUncer(sig1DvarY, bkg1DvarY, varYcutdir)
    rocvarY.SetLineColor(4)        
    
    file1droc = TFile(outputdir+"Top_ROC2DOutput_"+varX+"_"+varY+"_pt"+pt1+pt2+"_"+flagFilterLabel+".root","RECREATE")
    rocvarX.Write("roc_"+varX)
    rocvarY.Write("roc_"+varY)
    file1droc.Close()

    #Make ratio of S/B for making likelihood ordering
    c0 = TCanvas("c0","c0",300,300);
    
    if rankMetric=="SoverSplusB":
        metriclabel="L = S / (S+B)"
        rat = sig.Clone("rat")
        rat.Sumw2()
        ratdiv = bkg.Clone("ratdiv")
        ratdiv.Add(sig)
        rat.Divide(ratdiv)        
    elif rankMetric=="SoverB":
        metriclabel="L = S / B"
        rat = sig.Clone("rat")
        rat.Sumw2()
        ratdiv = bkg.Clone("ratdiv")
        rat.Divide(ratdiv)
    else:
        print "The ranking metric is not correct - EXITTING ..."
        sys.exit()


    #convert hist to 1D array and store all the bin information and error information separately
    binnumX=[]
    binnumY=[]
    s=[]
    serr=[]
    b=[]
    berr=[]
    r=[]
    rerr=[]
    rerrfrac=[]    
    bzero=[]
    bbigerr=[]
    szero=[]
    sbigerr=[]

    for i in range(1,sig.GetNbinsX()+1):
        for j in range(1,sig.GetNbinsY()+1):
            #print "bin ",i,j
            binnumtempX = i
            binnumtempY = j
            stemp    = sig.GetBinContent(i,j)
            serrtemp = sig.GetBinError(i,j)
            btemp    = bkg.GetBinContent(i,j)
            berrtemp = bkg.GetBinError(i,j)
            rtemp    = rat.GetBinContent(i,j)
            rerrtemp = rat.GetBinError(i,j)

            binnumX.append(binnumtempX)
            binnumY.append(binnumtempY)
            s.append(stemp)
            serr.append(serrtemp)
            b.append(btemp)
            berr.append(berrtemp)
            r.append(rtemp)
            rerr.append(rerrtemp)
            if rtemp==0:
                rerrfrac.append(1.0)
            else:
                rerrfrac.append(rerrtemp/rtemp)
            
            bzero.append(0)
            bbigerr.append(0)
            if b[-1]==0:
                bzero[-1]=1
            elif berr[-1]/b[-1]>0.7:
                bbigerr[-1]=1
        
            szero.append(0)
            sbigerr.append(0)
            if s[-1]==0:
                szero[-1]=1
            elif serr[-1]/s[-1]>0.7:
                sbigerr[-1]=1


    #Initial histograms
    cstep = TCanvas("cstep","cstep",1000,1000)   
    cstep.Divide(1,5)
    hab = TH1F("hab","hab",len(s),0,len(s))
    hab.SetMaximum( 1.4*max([max(b),max(s)]) )
    hab.SetLineColor(4)
    hab.GetXaxis().SetTitle("Bin Number")
    hab.GetXaxis().SetTitleSize(0.1)
    hab.GetXaxis().SetTitleOffset(0.7)
    hab.GetXaxis().SetLabelSize(0.1)
    hab.GetYaxis().SetLabelSize(0.1)
    hab.GetYaxis().SetTitle("Fraction Events")
    hab.GetYaxis().SetTitleOffset(0.5)
    hab.GetYaxis().SetTitleSize(0.1)
    hab.GetYaxis().CenterTitle()
    has = TH1F("has","has",len(s),0,len(s)) 
    has.SetLineColor(2)
    har = TH1F("har","har",len(s),0,len(s))
    har.SetLineColor(1)
    har.SetLineWidth(3)
    har.SetMaximum(2.0)
    har.GetXaxis().SetTitle("Bin Number")
    har.GetXaxis().SetTitleSize(0.1)
    har.GetXaxis().SetTitleOffset(0.7)
    har.GetXaxis().SetLabelSize(0.1)
    har.GetYaxis().SetLabelSize(0.1)
    har.GetYaxis().SetTitle("Ratio and Error")
    har.GetYaxis().SetTitleOffset(0.5)
    har.GetYaxis().SetTitleSize(0.1)
    har.GetYaxis().CenterTitle()
    hare = TH1F("hare","hare",len(s),0,len(s))
    hare.SetLineColor(3)
    hare.SetMaximum(3.0)
    haref = TH1F("haref","haref",len(s),0,len(s))
    haref.SetLineColor(7)
    haref.SetMaximum(2.0)
    haref.GetXaxis().SetTitle("Bin Number")
    haref.GetXaxis().SetTitleSize(0.1)
    haref.GetXaxis().SetTitleOffset(0.7)
    haref.GetXaxis().SetLabelSize(0.1)
    haref.GetYaxis().SetLabelSize(0.1)
    haref.GetYaxis().SetTitle("Frac. Error")
    haref.GetYaxis().SetTitleOffset(0.5)
    haref.GetYaxis().SetTitleSize(0.1)
    haref.GetYaxis().CenterTitle()
    hbwithzero = TH1F("hbwithzero","hbwithzero",len(s),0,len(s))
    hbwithzero.SetLineColor(8)
    hbwithzero.SetMaximum(2.0)
    hbwithzero.GetXaxis().SetTitle("Bin Number")
    hbwithzero.GetXaxis().SetTitleSize(0.1)
    hbwithzero.GetXaxis().SetTitleOffset(0.7)
    hbwithzero.GetXaxis().SetLabelSize(0.1)
    hbwithzero.GetYaxis().SetLabelSize(0.1)
    hbwithzero.GetYaxis().SetTitle("Yes(1) or No(0)")
    hbwithzero.GetYaxis().SetTitleOffset(0.5)
    hbwithzero.GetYaxis().SetTitleSize(0.1)
    hbwithzero.GetYaxis().CenterTitle()
    hbwithone = TH1F("hbwithone","hbwithone",len(s),0,len(s))
    hbwithone.SetLineColor(9)
    hswithzero = TH1F("hswithzero","hswithzero",len(s),0,len(s))
    hswithzero.SetLineColor(5)
    hswithzero.SetMaximum(2.0)
    hswithzero.GetXaxis().SetTitle("Bin Number")
    hswithzero.GetXaxis().SetTitleSize(0.1)
    hswithzero.GetXaxis().SetTitleOffset(0.7)
    hswithzero.GetXaxis().SetLabelSize(0.1)
    hswithzero.GetYaxis().SetLabelSize(0.1)
    hswithzero.GetYaxis().SetTitle("Yes(1) or No(0)")
    hswithzero.GetYaxis().SetTitleOffset(0.5)
    hswithzero.GetYaxis().SetTitleSize(0.1)
    hswithzero.GetYaxis().CenterTitle()
    hswithone = TH1F("hswithone","hswithone",len(s),0,len(s))
    hswithone.SetLineColor(6)



    for i in range(len(s)):
        #print s[i],b[i],r[i]
        has.SetBinContent(i+1,s[i])
        has.SetBinError(i+1,serr[i])
        hab.SetBinContent(i+1,b[i])
        hab.SetBinError(i+1,berr[i])
        har.SetBinContent(i+1,r[i])
        hare.SetBinContent(i+1,rerr[i])
        haref.SetBinContent(i+1,rerrfrac[i])
        hbwithzero.SetBinContent(i+1,bzero[i])
        hbwithone.SetBinContent(i+1 ,bbigerr[i])
        hswithzero.SetBinContent(i+1,szero[i])
        hswithone.SetBinContent(i+1 ,sbigerr[i])
    cstep.cd(1)
    hab.Draw("hist")
    has.Draw("histsame")
    myText(        0.70,0.85,1,0.1,varX+" & "+varY)
    myLineBoxText(0.30, 0.85, 2, 1, 0, 0, 0.1, 0.2, "Signal")
    myLineBoxText(0.50, 0.85, 4, 1, 0, 0, 0.1, 0.2, "Background")
    cstep.cd(2)
    har.Draw("hist")
    hare.Draw("histsame")
    myLineBoxText(0.30, 0.85, 1, 1, 0, 0, 0.1, 0.2, "Ratio")
    myLineBoxText(0.50, 0.85, 3, 1, 0, 0, 0.1, 0.2, "Ratio Error")
    cstep.cd(3)
    haref.Draw("hist")
    myLineBoxText(0.30, 0.85, 7, 1, 0, 0, 0.1, 0.2, "Fractional Ratio Error")
    cstep.cd(4)
    hbwithzero.Draw("hist")
    hbwithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 8, 1, 0, 0, 0.1, 0.2, "N(b)==0")
    myLineBoxText(0.50, 0.85, 9, 1, 0, 0, 0.1, 0.2, "frac_err(b)>0.7")
    cstep.cd(5)
    hswithzero.Draw("hist")
    hswithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 5, 1, 0, 0, 0.1, 0.2, "N(s)==0")
    myLineBoxText(0.50, 0.85, 6, 1, 0, 0, 0.1, 0.2, "frac_err(s)>0.7")
    #cstep.SaveAs(outputdir+"BinReodering_"+alg+"_"+varX+"_"+varY+"_stepA_pt"+pt1+pt2+"_"+flagFilterLabel+".eps")

    #if background have large statistical error then pop them to the back
    if flagFilter:
        i=0
        count=0
        while count<len(s):
            if bzero[i]==1 or bbigerr[i]==1:
                binnumX.append(binnumX.pop(i))
                binnumY.append(binnumY.pop(i))

                s.append(s.pop(i))
                serr.append(serr.pop(i))
                
                b.append(b.pop(i))
                berr.append(berr.pop(i))
                
                r.append(r.pop(i))
                rerr.append(rerr.pop(i))
                rerrfrac.append(rerrfrac.pop(i))

                bzero.append(bzero.pop(i))
                bbigerr.append(bbigerr.pop(i))
                szero.append(szero.pop(i))
                sbigerr.append(sbigerr.pop(i))
            else:
                i+=1
                
            count+=1



    for i in range(len(s)):
        #print s[i],b[i],r[i]
        has.SetBinContent(i+1,s[i])
        has.SetBinError(i+1,serr[i])
        hab.SetBinContent(i+1,b[i])
        hab.SetBinError(i+1,berr[i])
        har.SetBinContent(i+1,r[i])
        hare.SetBinContent(i+1,rerr[i])
        haref.SetBinContent(i+1,rerrfrac[i])
        hbwithzero.SetBinContent(i+1,bzero[i])
        hbwithone.SetBinContent(i+1 ,bbigerr[i])
        hswithzero.SetBinContent(i+1,szero[i])
        hswithone.SetBinContent(i+1 ,sbigerr[i])
    cstep.cd(1)
    hab.Draw("hist")
    has.Draw("histsame")
    myText(        0.70,0.85,1,0.1,varX+" & "+varY)
    myText(        0.70,0.75,1,0.1,"Move Large BG Err.")
    myLineBoxText(0.30, 0.85, 2, 1, 0, 0, 0.1, 0.3, "Signal")
    myLineBoxText(0.50, 0.85, 4, 1, 0, 0, 0.1, 0.3, "Background")
    cstep.cd(2)
    har.Draw("hist")
    hare.Draw("histsame")
    myLineBoxText(0.30, 0.85, 1, 1, 0, 0, 0.1, 0.3, "Ratio")
    myLineBoxText(0.50, 0.85, 3, 1, 0, 0, 0.1, 0.3, "Ratio Error")
    cstep.cd(3)
    haref.Draw("hist")
    myLineBoxText(0.30, 0.85, 7, 1, 0, 0, 0.1, 0.3, "Fractional Ratio Error")
    cstep.cd(4)
    hbwithzero.Draw("hist")
    hbwithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 8, 1, 0, 0, 0.1, 0.3, "N(b)==0")
    myLineBoxText(0.50, 0.85, 9, 1, 0, 0, 0.1, 0.3, "frac_err(b)>0.7")
    cstep.cd(5)
    hswithzero.Draw("hist")
    hswithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 5, 1, 0, 0, 0.1, 0.3, "N(s)==0")
    myLineBoxText(0.50, 0.85, 6, 1, 0, 0, 0.1, 0.3, "frac_err(s)>0.7")
    #cstep.SaveAs(outputdir+"BinReordering_"+alg+"_"+varX+"_"+varY+"_stepB1_pt"+pt1+pt2+"_"+flagFilterLabel+".eps")



    #if background have large statistical error then pop them to the back
    if flagFilter:
        i=0
        count=0
        while count<len(s):
            if szero[i]==1 or sbigerr[i]==1:
                binnumX.append(binnumX.pop(i))
                binnumY.append(binnumY.pop(i))

                s.append(s.pop(i))
                serr.append(serr.pop(i))
                
                b.append(b.pop(i))
                berr.append(berr.pop(i))
                
                r.append(r.pop(i))
                rerr.append(rerr.pop(i))
                rerrfrac.append(rerrfrac.pop(i))

                bzero.append(bzero.pop(i))
                bbigerr.append(bbigerr.pop(i))
                szero.append(szero.pop(i))
                sbigerr.append(sbigerr.pop(i))
            else:
                i+=1
                
            count+=1

    for i in range(len(s)):
        #print s[i],b[i],r[i]
        has.SetBinContent(i+1,s[i])
        has.SetBinError(i+1,serr[i])
        hab.SetBinContent(i+1,b[i])
        hab.SetBinError(i+1,berr[i])
        har.SetBinContent(i+1,r[i])
        hare.SetBinContent(i+1,rerr[i])
        haref.SetBinContent(i+1,rerrfrac[i])
        hbwithzero.SetBinContent(i+1,bzero[i])
        hbwithone.SetBinContent(i+1 ,bbigerr[i])
        hswithzero.SetBinContent(i+1,szero[i])
        hswithone.SetBinContent(i+1 ,sbigerr[i])
    cstep.cd(1)
    hab.Draw("hist")
    has.Draw("histsame")
    myText(        0.70,0.85,1,0.1,varX+" & "+varY)
    myText(        0.70,0.75,1,0.1,"Move Large BG Err.")
    myText(        0.70,0.65,1,0.1,"Move Large Sig Err.")
    myLineBoxText(0.30, 0.85, 2, 1, 0, 0, 0.1, 0.3, "Signal")
    myLineBoxText(0.50, 0.85, 4, 1, 0, 0, 0.1, 0.3, "Background")
    cstep.cd(2)
    har.Draw("hist")
    hare.Draw("histsame")
    myLineBoxText(0.30, 0.85, 1, 1, 0, 0, 0.1, 0.3, "Ratio")
    myLineBoxText(0.50, 0.85, 3, 1, 0, 0, 0.1, 0.3, "Ratio Error")
    cstep.cd(3)
    haref.Draw("hist")
    myLineBoxText(0.30, 0.85, 7, 1, 0, 0, 0.1, 0.3, "Fractional Ratio Error")
    cstep.cd(4)
    hbwithzero.Draw("hist")
    hbwithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 8, 1, 0, 0, 0.1, 0.3, "N(b)==0")
    myLineBoxText(0.50, 0.85, 9, 1, 0, 0, 0.1, 0.3, "frac_err(b)>0.7")
    cstep.cd(5)
    hswithzero.Draw("hist")
    hswithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 5, 1, 0, 0, 0.1, 0.3, "N(s)==0")
    myLineBoxText(0.50, 0.85, 6, 1, 0, 0, 0.1, 0.3, "frac_err(s)>0.7")
    #cstep.SaveAs(outputdir+"BinReordering_"+alg+"_"+varX+"_"+varY+"_stepB2_pt"+pt1+pt2+"_"+flagFilterLabel+".eps")
    

    #if bin has statistical uncertainty greater than 80% then remove it from the possible signal region
    if flagFilter:
        i=0
        for i in range(len(s)):
            if r[i]!=0:
                if rerrfrac[i]>statErrorOnRatioThreshold:
                    r[i]=0
                    #rerrfrac[i]=0

    for i in range(len(s)):
        #print s[i],b[i],r[i]
        has.SetBinContent(i+1,s[i])
        has.SetBinError(i+1,serr[i])
        hab.SetBinContent(i+1,b[i])
        hab.SetBinError(i+1,berr[i])
        har.SetBinContent(i+1,r[i])
        hare.SetBinContent(i+1,rerr[i])
        haref.SetBinContent(i+1,rerrfrac[i])
        hbwithzero.SetBinContent(i+1,bzero[i])
        hbwithone.SetBinContent(i+1 ,bbigerr[i])
        hswithzero.SetBinContent(i+1,szero[i])
        hswithone.SetBinContent(i+1 ,sbigerr[i])
    cstep.cd(1)
    hab.Draw("hist")
    has.Draw("histsame")
    myText(        0.70,0.85,1,0.1,varX+" & "+varY)
    myText(        0.70,0.75,1,0.1,"Move Large BG Err.")
    myText(        0.70,0.65,1,0.1,"Move Large Sig Err.")
    myText(        0.70,0.55,1,0.1,"Move Large Ratio Err.")
    myLineBoxText(0.30, 0.85, 2, 1, 0, 0, 0.1, 0.3, "Signal")
    myLineBoxText(0.50, 0.85, 4, 1, 0, 0, 0.1, 0.3, "Background")
    cstep.cd(2)
    har.Draw("hist")
    hare.Draw("histsame")
    myLineBoxText(0.30, 0.85, 1, 1, 0, 0, 0.1, 0.3, "Ratio")
    myLineBoxText(0.50, 0.85, 3, 1, 0, 0, 0.1, 0.3, "Ratio Error")
    cstep.cd(3)
    haref.Draw("hist")
    myLineBoxText(0.30, 0.85, 7, 1, 0, 0, 0.1, 0.3, "Fractional Ratio Error")
    cstep.cd(4)
    hbwithzero.Draw("hist")
    hbwithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 8, 1, 0, 0, 0.1, 0.3, "N(b)==0")
    myLineBoxText(0.50, 0.85, 9, 1, 0, 0, 0.1, 0.3, "frac_err(b)>0.7")
    cstep.cd(5)
    hswithzero.Draw("hist")
    hswithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 5, 1, 0, 0, 0.1, 0.3, "N(s)==0")
    myLineBoxText(0.50, 0.85, 6, 1, 0, 0, 0.1, 0.3, "frac_err(s)>0.7")
    #cstep.SaveAs(outputdir+"BinReordering_"+alg+"_"+varX+"_"+varY+"_stepB3_pt"+pt1+pt2+"_"+flagFilterLabel+".eps")

    #order bins based on the ratio histogram content
    for i in range(len(s)):
        ifix = len(s)-i
        for j in range(1,ifix-1):
            if r[j]<r[j+1]:
                binnumX[j],binnumX[j+1]=binnumX[j+1],binnumX[j]
                binnumY[j],binnumY[j+1]=binnumY[j+1],binnumY[j]

                s[j],s[j+1] = s[j+1],s[j]
                serr[j],serr[j+1] = serr[j+1],serr[j]
                
                b[j],b[j+1] = b[j+1],b[j]
                berr[j],berr[j+1] = berr[j+1],berr[j]
                
                r[j],r[j+1] = r[j+1],r[j]
                rerr[j],rerr[j+1] = rerr[j+1],rerr[j]
                rerrfrac[j],rerrfrac[j+1] = rerrfrac[j+1],rerrfrac[j]

                bzero[j],bzero[j+1] = bzero[j+1],bzero[j]
                bbigerr[j],bbigerr[j+1]   = bbigerr[j+1],bbigerr[j]
                szero[j],szero[j+1] = szero[j+1],szero[j]
                sbigerr[j],sbigerr[j+1]   = sbigerr[j+1],sbigerr[j]

    for i in range(len(s)):
        #print s[i],b[i],r[i]
        has.SetBinContent(i+1,s[i])
        has.SetBinError(i+1,serr[i])
        hab.SetBinContent(i+1,b[i])
        hab.SetBinError(i+1,berr[i])
        har.SetBinContent(i+1,r[i])
        hare.SetBinContent(i+1,rerr[i])
        haref.SetBinContent(i+1,rerrfrac[i])
        hbwithzero.SetBinContent(i+1,bzero[i])
        hbwithone.SetBinContent(i+1 ,bbigerr[i])
        hswithzero.SetBinContent(i+1,szero[i])
        hswithone.SetBinContent(i+1 ,sbigerr[i])
    cstep.cd(1)
    hab.Draw("hist")
    has.Draw("histsame")
    myText(        0.70,0.85,1,0.1,varX+" & "+varY)
    myText(        0.70,0.75,1,0.1,"Move Large BG Err.")
    myText(        0.70,0.65,1,0.1,"Move Large Sig Err.")
    myText(        0.70,0.55,1,0.1,"Move Large Ratio Err.")
    myText(        0.70,0.45,1,0.1,"Reorder by Ratio")
    myLineBoxText(0.30, 0.85, 2, 1, 0, 0, 0.1, 0.3, "Signal")
    myLineBoxText(0.50, 0.85, 4, 1, 0, 0, 0.1, 0.3, "Background")
    cstep.cd(2)
    har.Draw("hist")
    hare.Draw("histsame")
    myLineBoxText(0.30, 0.85, 1, 1, 0, 0, 0.1, 0.3, "Ratio")
    myLineBoxText(0.50, 0.85, 3, 1, 0, 0, 0.1, 0.3, "Ratio Error")
    cstep.cd(3)
    haref.Draw("hist")
    myLineBoxText(0.30, 0.85, 7, 1, 0, 0, 0.1, 0.3, "Fractional Ratio Error")
    cstep.cd(4)
    hbwithzero.Draw("hist")
    hbwithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 8, 1, 0, 0, 0.1, 0.3, "N(b)==0")
    myLineBoxText(0.50, 0.85, 9, 1, 0, 0, 0.1, 0.3, "frac_err(b)>0.7")
    cstep.cd(5)
    hswithzero.Draw("hist")
    hswithone.Draw("histsame")
    myLineBoxText(0.30, 0.85, 5, 1, 0, 0, 0.1, 0.3, "N(s)==0")
    myLineBoxText(0.50, 0.85, 6, 1, 0, 0, 0.1, 0.3, "frac_err(s)>0.7")
    #cstep.SaveAs(outputdir+"BinReordering_"+alg+"_"+varX+"_"+varY+"_stepC_pt"+pt1+pt2+"_"+flagFilterLabel+".eps")

    #find value of first bin with zero ratio value
    zerobin=0
    for i in range(len(s)):
        zerobin=i
        if r[i]==0:
            break
    print "Found first zero at: ",zerobin


    totalB = sum(b)
    totalS = sum(s)
    n = len(s)

    hsigout = TH1F("hsigout","hsigout",n,0,1)
    hsigout.SetDirectory(0)
    hbkgout = TH1F("hbkgout","hbkgout",n,0,1)
    hbkgout.SetDirectory(0)

    siglowX  = sig.GetXaxis().GetXmin()
    sighighX = sig.GetXaxis().GetXmax()
    nbinsX   = sig.GetXaxis().GetNbins()
    siglowY  = sig.GetYaxis().GetXmin()
    sighighY = sig.GetYaxis().GetXmax()
    nbinsY   = sig.GetYaxis().GetNbins()

    h1 = TH2F("h1","h1",nbinsX,siglowX,sighighX,nbinsY,siglowY,sighighY)
    h1.SetDirectory(0)

    gr = TGraphErrors(n)
    for i in range(1,n+1):
        myS = 0.
        myB = 0.

        myBerr=Double()
        mySerr=Double()
        myB = hab.IntegralAndError(1,i,myBerr)
        myS = has.IntegralAndError(1,i,mySerr)
        
        print i,"  myS=",myS,"  myB=",myB,"  Errors: ",mySerr, myBerr
        
        gr.SetPoint(i, myS, myB)
        gr.SetPointError(i, mySerr, myBerr)
        
        if myS<=0.73:
            h1.SetBinContent(binnumX[i], binnumY[i], 1)

    #get histograms that are colored for signal efficiency at 50%
    cc = TCanvas("cc","cc",500,500);
    cc.Divide(2,2);
    cb1 = cc.cd(1);
    cb2 = cc.cd(2);
    cb3 = cc.cd(3);
    cb4 = cc.cd(4);
    cb1.cd()
    bkg.GetXaxis().SetTitle(TranslateVar(varX))
    bkg.GetYaxis().SetTitle(TranslateVar(varY))
    bkg.Draw("colz")
    ATLASLabel(   0.20,0.90,1,0.12,0.04,"#sqrt{s}=13 TeV")
    myText(       0.20,0.80,1,0.04, "QCD Jets")
    myText(       0.50,0.90,1,0.04, TranslateAlg(alg))
    myText(       0.20,0.70,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
    myText(       0.50,0.80,1,0.04, "corr = "+str(round(corrbg,4)))
    cb2.cd()
    sig.GetXaxis().SetTitle(TranslateVar(varX))
    sig.GetYaxis().SetTitle(TranslateVar(varY))
    sig.Draw("colz")
    ATLASLabel(   0.20,0.90,1,0.12,0.04,"#sqrt{s}=13 TeV")
    myText(       0.20,0.80,1,0.04, "Top Jets")
    myText(       0.50,0.90,1,0.04, TranslateAlg(alg))
    myText(       0.20,0.70,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
    myText(       0.50,0.80,1,0.04, "corr = "+str(round(corrsig,4)))
    cb3.cd()
    cb3.SetLogz()
    rat.GetXaxis().SetTitle(TranslateVar(varX))
    rat.GetYaxis().SetTitle(TranslateVar(varY))
    rat.Draw("colz")
    ATLASLabel(   0.20,0.90,1,0.14,0.04,"#sqrt{s}=13 TeV")
    myText(       0.20,0.85,1,0.04,metriclabel)
    cb4.cd()
    h1.GetXaxis().SetTitle(TranslateVar(varX))
    h1.GetYaxis().SetTitle(TranslateVar(varY))
    h1.Draw("colz")
    ATLASLabel(   0.20,0.90,1,0.14,0.04,"#sqrt{s}=13 TeV")  
    myText(       0.20,0.85,1,0.04,"Accepted Region")
    myText(       0.20,0.80,1,0.04,"at 50% Signal Eff.")
    cc.SaveAs(outputdir+"Final2DSpectrum_"+alg+"_"+varX+"_"+varY+"_pt"+pt1+pt2+"_"+flagFilterLabel+".eps")
    
    
    
    #PRINT OUT SEPARATELY
    if False:
        #get histograms that are colored for signal efficiency at 50%
        csep = TCanvas("csep","csep",444,400);
        csep.SetRightMargin(0.15)
        bkg.GetXaxis().SetTitle(TranslateVar(varX))
        bkg.GetYaxis().SetTitle(TranslateVar(varY))
        bkg.Draw("colz")
        ATLASLabel(   0.18,0.90,1,0.09,0.030,"#sqrt{s}=13 TeV")
        myText(       0.18,0.86,1,0.030,"#sqrt{s}=13 TeV")
        myText(       0.20,0.70,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
        myText(       0.46,0.90,1,0.030, "QCD Jets")
        myText(       0.46,0.86,1,0.030, TranslateAlg(alg))
        csep.SetLogz(1)
        csep.SaveAs(outputdir+"2DSpectrum_"+alg+"_"+varX+"_"+varY+"_pt"+pt1+pt2+"_bkg_log.eps")

        sig.GetXaxis().SetTitle(TranslateVar(varX))
        sig.GetYaxis().SetTitle(TranslateVar(varY))
        sig.Draw("colz")
        ATLASLabel(   0.18,0.90,1,0.09,0.030,"#sqrt{s}=13 TeV")
        myText(       0.20,0.70,1,0.04, TranslateRegion(pt1,pt2,m1,m2))        
        myText(       0.46,0.90,1,0.030, "Top Jets")
        myText(       0.46,0.86,1,0.030, TranslateAlg(alg))
        csep.SetLogz(1)
        csep.SaveAs(outputdir+"2DSpectrum_"+alg+"_"+varX+"_"+varY+"_pt"+pt1+pt2+"_sig_log.eps")

        rat.GetXaxis().SetTitle(TranslateVar(varX))
        rat.GetYaxis().SetTitle(TranslateVar(varY))
        rat.Draw("colz")
        ATLASLabel(   0.20,0.90,1,0.14,0.04,"#sqrt{s}=13 TeV")
        myText(       0.48,0.90,1,0.04, metriclabel)
        myText(       0.48,0.84,1,0.04, TranslateAlg(alg))
        myText(       0.20,0.70,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
        csep.SetLogz(0)
        csep.SaveAs(outputdir+"2DSpectrum_"+alg+"_"+varX+"_"+varY+"_pt"+pt1+pt2+"_likelihood.eps")
    
        h1.GetXaxis().SetTitle(TranslateVar(varX))
        h1.GetYaxis().SetTitle(TranslateVar(varY))
        h1.Draw("colz")
        ATLASLabel(   0.20,0.90,1,0.14,0.04,"#sqrt{s}=13 TeV")  
        myText(       0.48,0.90,1,0.04, "Signal Region")
        myText(       0.48,0.84,1,0.04, TranslateAlg(alg))
        myText(       0.20,0.70,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
        csep.SetLogz(0)
        csep.SaveAs(outputdir+"2DSpectrum_"+alg+"_"+varX+"_"+varY+"_pt"+pt1+pt2+"_signalregion.eps")
    
    cgr = TCanvas("cgr","cgr",500,500);
    gr.GetXaxis().SetTitle("Signal Efficiency")
    gr.GetYaxis().SetTitle("Background Efficiency")
    gr.SetMarkerStyle(8)
    gr.SetMarkerSize(0.3)
    gr.SetMinimum(0.0)
    gr.SetMaximum(1.0)
    gr.Draw("ACE3")
    
    rocvarX.SetFillColor(2)
    rocvarX.SetFillStyle(3001)
    rocvarX.Draw("CE3same")
    rocvarY.SetFillColor(4)
    rocvarY.SetFillStyle(3001)
    rocvarY.Draw("CE3same")
    
    line = TLine(0.73,0.0,0.73,1.0)
    line.SetLineColor(1)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw("same")
    
    ATLASLabel(   0.20,0.90,1,0.1,0.03,"#sqrt{s}=13 TeV")
    myText(       0.20,0.85,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
    rocbox0=myLineBoxText(0.26, 0.80, 1, 1, 1, 0, 0.1, 0.08, TranslateVar(varX)+" & "+TranslateVar(varY))
    rocbox1=myLineBoxText(0.26, 0.75, 2, 2, 2, 3001, 0.1, 0.08, TranslateVar(varX))
    rocbox2=myLineBoxText(0.26, 0.70, 4, 1, 4, 3001, 0.1, 0.08, TranslateVar(varY))
    myText(       0.20,0.65,1,0.03,"Stat Err Threshold = "+str(statErrorOnRatioThreshold))
    cgr.SaveAs(outputdir+"Top_"+"ROCComparison_"+alg+"_"+varX+"_"+varY+"_pt"+pt1+pt2+"_"+flagFilterLabel+".eps")

    return gr,has,hab,h1,zerobin

def SignalBGCompare2D(InputDir, alg, var1, var1cutdir, var2, var2cutdir, range1, range2, pt1, pt2, m1, m2, outputdir):
    '''Implementation of simple signal and background comparison'''
    print "getting correlation: ",alg,var1,var2
    
    c = TCanvas("c","c",300,300)

    dry=0.045

    weight=""
    weight+="("
    weight+=alg+"_pt>"+pt1+" && "
    weight+=alg+"_pt<"+pt2
    if m1!="0":
        weight+=" && "+alg+"_m>"+m1+" && "
        weight+=alg+"_m<"+m2
    weight+=")"
    
    #Get signal and background histograms
    hsig = GetHist2D(InputDir+sigFile, "JetTree", alg+"_"+var1, alg+"_"+var2, range1, range2, weight+"*("+alg+"_flavor==3)")
    hbkg = GetHist2D(InputDir+bkgFile, "JetTree", alg+"_"+var1, alg+"_"+var2, range1, range2, weight+"*("+alg+"_flavor==0)")

    # rebin and normalize
    rebinX=2
    rebinY=2
    hsig.RebinX(rebinX)
    hsig.RebinY(rebinY)
    hbkg.RebinX(rebinX)
    hbkg.RebinY(rebinY)
    
    hsig = NormalizeHist(hsig)
    hbkg = NormalizeHist(hbkg)
    
    corrsig = hsig.GetCorrelationFactor()
    corrbkg = hbkg.GetCorrelationFactor()

    c = TCanvas("c","c",600,300);
    c.Divide(2,1);
    c1 = c.cd(1);
    c2 = c.cd(2);
    c1.cd()
    hbkg.GetXaxis().SetTitle(var1)
    hbkg.GetYaxis().SetTitle(var2)
    hbkg.Draw("colz")
    ATLASLabel(   0.20,0.90,1,0.12,0.04,"#sqrt{s}=13 TeV")
    myText(       0.20,0.80,1,0.04, "QCD Jets")
    myText(       0.20,0.75,1,0.04, TranslateAlg(alg))
    myText(       0.20,0.80,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
    myText(       0.20,0.65,1,0.04, "corr = "+str(corrbkg))
    c2.cd()
    hsig.Draw("colz")
    hsig.GetXaxis().SetTitle(var1)
    hsig.GetYaxis().SetTitle(var2)
    ATLASLabel(   0.20,0.90,1,0.12,0.04,"#sqrt{s}=13 TeV")
    myText(       0.20,0.80,1,0.04, "Top Jets")
    myText(       0.20,0.75,1,0.04, TranslateAlg(alg))
    myText(       0.20,0.80,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
    myText(       0.20,0.65,1,0.04, "corr = "+str(corrsig))
    c.SaveAs(outputdir+"Top_"+"Correlation_"+alg+"_"+var1+"_"+var2+"_pt"+pt1+pt2+".eps")
    

    #FILTER
    flagFilter=1
    flagFilterLabel="Filter"
    roc,sigordered,bkgordered,h50selected,zerobin = Make2DROC(alg, var1, var1cutdir, var2, var2cutdir, hsig, hbkg, "SoverSplusB", 0.7, flagFilter, outputdir)
    f = TFile(outputdir+"Top_ROC2D_"+alg+"_"+var1+"_"+var2+"_pt"+pt1+pt2+"_"+flagFilterLabel+".root","RECREATE")
    roc.Write("ROC")
    f.Close()
    
    f = TFile(outputdir+"Top_SignalRegion2D_"+alg+"_"+var1+"_"+var2+"_pt"+pt1+pt2+"_"+flagFilterLabel+".root","RECREATE")
    h50selected.Write("SignalRegion")
    f.Close()

    #Likelihood for cutting on - reordered hists
    cc = TCanvas("cc","cc",600,600)
    sigordered.SetLineColor(4)
    sigordered.SetLineWidth(3)
    bkgordered.SetLineColor(2)
    bkgordered.SetLineWidth(3)

    sigordered.Rebin(20)
    bkgordered.Rebin(20)

    sigordered = NormalizeHist(sigordered)
    bkgordered = NormalizeHist(bkgordered)

    maxval = GetMaxVal([bkgordered, sigordered])
    bkgordered.SetMaximum(1.5*maxval)

    bkgordered.GetXaxis().SetTitle("Likelihood")
    bkgordered.GetXaxis().SetTitleSize(0.05)
    bkgordered.GetXaxis().SetTitleOffset(1.2)
    bkgordered.GetXaxis().SetLabelSize(0.05)
    bkgordered.GetYaxis().SetTitle("Normalized Events")
    bkgordered.GetYaxis().SetTitleSize(0.05)
    bkgordered.GetYaxis().SetTitleOffset(1.5)
    bkgordered.GetYaxis().SetLabelSize(0.05)
    bkgordered.Draw("hist")
    sigordered.Draw("histsame")

    line = TLine(zerobin,0.0,zerobin,1.5*maxval)
    line.SetLineColor(1)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw("same")

    ATLASLabel(    0.20,0.90,1,0.1,0.03,"#sqrt{s}=13 TeV")
    myText(       0.20,0.80,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
    box0=myLineBoxText(0.67, 0.80, 2, 2, 2, 3004, 0.1, 0.1, "QCD Jets")
    box1=myLineBoxText(0.67, 0.75, 4, 1, 4, 3005, 0.1, 0.1, "Top Jets")
    myText(        0.60,0.70,1,0.03, "2D Likelihood")
    myText(        0.60,0.65,1,0.03, "0 = Signal-like")
    myText(        0.60,0.60,1,0.03, "1 = Background-like")

def GetCorrelationAndSeparationSigBG(InputDir, alg, var1, var2, range1, range2, pt1, pt2, m1, m2, outputdir):
    '''Implementation of simple signal and background comparison'''
    print "getting correlation: ",alg,var1,var2
    
    debug=0
    
    outfname1=outputdir+"Top_CorrelationCalc_"+alg+"_"+var1+"_"+var2+"_pt"+pt1+pt2
    outfname2=outputdir+"Top_CorrelationCalc_"+alg+"_"+var2+"_"+var1+"_pt"+pt1+pt2
   
    print outfname1
    print outfname2
   
    #first check if the output file exists
    #if it does then take the correlation from that
    if os.path.isfile(outfname1+".root"):
        print "FOUND CORRELATION FILE 1"
        ftemp = TFile(outfname1+".root")
        htemp = ftemp.Get("corr")
        corrsig = htemp.GetBinContent(1)
        corrbg = htemp.GetBinContent(2)
        return corrsig,corrbg
    
    if os.path.isfile(outfname2+".root"):
        print "FOUND CORRELATION FILE 1"
        ftemp = TFile(outfname2+".root")
        htemp = ftemp.Get("corr")
        corrsig = htemp.GetBinContent(1)
        corrbg = htemp.GetBinContent(2)
        return corrsig,corrbg
    
    #if neither exists you must go and recreate them
    c = TCanvas("c","c",300,300)

    dry=0.045

    weight=""
    weight+="("
    weight+=alg+"_pt>"+pt1+" && "
    weight+=alg+"_pt<"+pt2+")"
    
    #Get signal and background histograms
    hsig = GetHist2D(InputDir+sigFile, "JetTree", alg+"_"+var1, alg+"_"+var2, range1, range2, weight+"*("+alg+"_flavor==3)")
    hbkg = GetHist2D(InputDir+bkgFile, "JetTree", alg+"_"+var1, alg+"_"+var2, range1, range2, weight+"*("+alg+"_flavor==0)")

    
    # NORMALIZE
    hsig = NormalizeHist(hsig)
    hbkg = NormalizeHist(hbkg)    
    
    # CORRELATION
    corrsig = hsig.GetCorrelationFactor()
    corrbkg = hbkg.GetCorrelationFactor()

    c = TCanvas("csig","csig",600,300);
    c.Divide(2,1);
    c1 = c.cd(1);
    c2 = c.cd(2);
    c1.cd()
    hbkg.GetXaxis().SetTitle(var1)
    hbkg.GetYaxis().SetTitle(var2)
    hbkg.Draw("colz")
    ATLASLabel(   0.20,0.90,1,0.12,0.04,"#sqrt{s}=13 TeV")
    myText(       0.20,0.85,1,0.04,"#sqrt{s}=13 TeV")
    myText(       0.20,0.80,1,0.04, "QCD Jets")
    myText(       0.20,0.75,1,0.04, TranslateAlg(alg))
    myText(       0.20,0.80,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
    myText(       0.20,0.65,1,0.04, "corr = "+str(corrbkg))
    c2.cd()
    hsig.Draw("colz")
    hsig.GetXaxis().SetTitle(var1)
    hsig.GetYaxis().SetTitle(var2)
    ATLASLabel(   0.20,0.90,1,0.12,0.04,"#sqrt{s}=13 TeV")
    myText(       0.20,0.80,1,0.04, "Top Jets")
    myText(       0.20,0.75,1,0.04, TranslateAlg(alg))
    myText(       0.20,0.80,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
    myText(       0.20,0.65,1,0.04, "corr = "+str(corrsig))
    #c.SaveAs(outfname1+".eps")
    
    # FISHER DISCRIMINANT
    proj1_sig = hsig.ProjectionX("proj1_sig")
    proj1_bkg = hbkg.ProjectionX("proj1_bkg")
    proj2_sig = hsig.ProjectionY("proj2_sig")
    proj2_bkg = hbkg.ProjectionY("proj2_bkg")
    
    c1.cd()
    proj1_bkg.SetLineColor(4)
    proj1_sig.SetLineColor(2)
    proj1_bkg.GetXaxis().SetTitle(var1)
    proj1_bkg.Draw("hist")
    proj1_sig.Draw("histsame")
    ATLASLabel(   0.20,0.90,1,0.12,0.04,"#sqrt{s}=13 TeV")
    bx1=myLineBoxText(0.56, 0.80, 2, 2, 0, 0, 0.08, 0.08, "Top Jets")
    bx2=myLineBoxText(0.75, 0.75, 4, 1, 0, 0, 0.08, 0.08, "QCD Jets")
    myText(       0.20,0.75,1,0.04, TranslateAlg(alg))
    myText(       0.20,0.65,1,0.04, "corr = "+str(corrbkg))
    c2.cd()
    proj2_bkg.SetLineColor(4)
    proj2_sig.SetLineColor(2)
    proj2_bkg.GetXaxis().SetTitle(var2)
    proj2_bkg.Draw("hist")
    proj2_sig.Draw("histsame")
    ATLASLabel(   0.20,0.90,1,0.12,0.04,"#sqrt{s}=13 TeV")
    bx1=myLineBoxText(0.50, 0.80, 2, 2, 0, 0, 0.08, 0.08, "Top Jets")
    bx2=myLineBoxText(0.50, 0.75, 4, 1, 0, 0, 0.08, 0.08, "QCD Jets")
    myText(       0.20,0.75,1,0.04, TranslateAlg(alg))
    myText(       0.20,0.80,1,0.04, TranslateRegion(pt1,pt2,m1,m2))
    #c.SaveAs(outfname1+"_fishercheck.eps")    

    mean1_sig = proj1_sig.GetMean()    
    mean1_bkg = proj1_bkg.GetMean()    
    rms1_sig = proj1_sig.GetRMS()    
    rms1_bkg = proj1_bkg.GetRMS()    
    print "GettingFisher 1: ",mean1_sig,mean1_bkg,rms1_sig,rms1_bkg 
    fisher1=10000.0
    if rms1_sig>0 and rms1_bkg>0:
        fisher1 = abs(mean1_sig-mean1_bkg)/(rms1_sig+rms1_bkg)
    print "Fisher 1: ",fisher1
    
    mean2_sig = proj2_sig.GetMean()    
    mean2_bkg = proj2_bkg.GetMean()    
    rms2_sig = proj2_sig.GetRMS()    
    rms2_bkg = proj2_bkg.GetRMS()    
    print "GettingFisher 2: ",mean2_sig,mean2_bkg,rms2_sig,rms2_bkg 
    fisher2=10000.0
    if rms2_sig>0 and rms2_bkg>0:
        fisher2 = abs(mean2_sig-mean2_bkg)/(rms2_sig+rms2_bkg)
    print "Fisher 2: ",fisher2
    
    #output file of correlation and fisher
    h1 = TH1F()
    h1.Fill("corrsig",corrsig)
    h1.Fill("corrbg",corrbkg)

    h1.Fill("mean1_sig",mean1_sig)
    h1.Fill("mean1_bkg",mean1_bkg)
    h1.Fill("rms1_sig",rms1_sig)
    h1.Fill("rms1_bkg",rms1_bkg)
    h1.Fill("fisher1",fisher1)

    h1.Fill("mean2_sig",mean2_sig)
    h1.Fill("mean2_bkg",mean2_bkg)
    h1.Fill("rms2_sig",rms2_sig)
    h1.Fill("rms2_bkg",rms2_bkg)
    h1.Fill("fisher2",fisher2)

    fout1 = TFile(outfname1+".root","RECREATE")
    h1.Write("corr")
    fout1.Close()

    fout2 = TFile(outfname2+".root","RECREATE")
    h1.Write("corr")
    fout2.Close()

    return corrsig,corrbkg
    
    
def MakeMVAROCS(outfilename,outputdir,mvatypes):
 
    print outfilename
    fin = TFile(outfilename)
    fin.ls()
    PlotTree = fin.Get("dataset/TestTree")
    PlotTree.Print()

    c = TCanvas("c","c",300,300)
    
    for mvatype in mvatypes.split(","):
        print "MVAType: ",mvatype

        #signal
        PlotTree.Draw(mvatype+">>hb1(200,-2,2)","(classID==1)*weight")
        #background
        PlotTree.Draw(mvatype+">>hb2(200,-2,2)","(classID==0)*weight")

        hsig = gDirectory.Get("hb1")
        hbkg = gDirectory.Get("hb2")
    
        hsig = NormalizeHist(hsig)
        hbkg = NormalizeHist(hbkg)

        hsig.Rebin(2)
        hbkg.Rebin(2)
        hsig.SetLineColor(4)
        hbkg.SetLineColor(2)
    
        max=GetMaxVal([hsig,hbkg])

        hsig.SetMaximum(1.8*max)
        hsig.GetXaxis().SetTitle(mvatype+" Output")
        hsig.GetYaxis().SetTitle("Arbitrary Units")
        hsig.Draw("hist")
        hbkg.Draw("histsame")

        ATLASLabel(   0.20,0.85,1,0.12,0.04,"#sqrt{s}=13 TeV")
        myText(       0.20,0.80,1,0.03, TranslateAlg(alg))
        myText(       0.20,0.75,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
        bx1=myLineBoxText(0.70, 0.85, 2, 1, 0, 0, 0.08, 0.1, "Top Jets")
        bx2=myLineBoxText(0.70, 0.80, 4, 1, 0, 0, 0.08, 0.1, "QCD Jets")

        c.SaveAs(outputdir+"Top_"+"SigBG_pt"+pt1+pt2+"_"+mvatype+".eps")

        rocL,hsigregL,hcutvalL,hsigregL25,hcutvalL25 = RocCurve_SingleSided_WithUncer(hsig, hbkg, "L")
        rocR,hsigregR,hcutvalR,hsigregR25,hcutvalR25 = RocCurve_SingleSided_WithUncer(hsig, hbkg, "R")
        rocSB,sigordered,bkgordered,h1 = RocCurve_SoverBOrdered_WithUncer(hsig, hbkg)

        fout = TFile(outfilename.replace(".root","_ROCS"+mvatype+".root"),"RECREATE")
        rocL.Write("ROC_L")
        rocR.Write("ROC_R")
        rocSB.Write("ROC_SoverB")
        fout.Close()
        
        
def OverlayROCS(outputdir1,outputdir2,outputdir3,outputdir4,alg,var0,var1,pt1,pt2,m1,m2,mvatypes,VarsAndRanges):
    print outputdir1,outputdir2,outputdir3,outputdir4,alg,var0,var1,pt1,pt2,m1,m2,mvatypes
    
    path = outputdir1+"Top_ROC_"+alg+"_"+var0+"_pt"+pt1+pt2+".root"
    f1   = TFile(path)
    roc1 = f1.Get("ROC_"+VarsAndRanges[var0][3])
    roc1.SetFillColor(2)
    roc1.SetLineColor(2)
    roc1.SetFillStyle(3001)
    
    path = outputdir1+"Top_ROC_"+alg+"_"+var1+"_pt"+pt1+pt2+".root"
    f2   = TFile(path)
    roc2 = f2.Get("ROC_"+VarsAndRanges[var1][3])
    roc2.SetFillColor(4)
    roc2.SetLineColor(4)
    roc2.SetFillStyle(3001)
    
    path = outputdir2+"Top_ROC2D_"+alg+"_"+var0+"_"+var1+"_pt"+pt1+pt2+"_Filter.root"
    f3   = TFile(path)
    roc3 = f3.Get("ROC")
    roc3.SetFillColor(3)
    roc3.SetLineColor(3)
    roc3.SetFillStyle(3001)
    
    path = outputdir3+"Top_TMVAOutput__"+alg+"__"+var0+"."+var1+"__pt"+pt1+pt2+"_ROCSLikelihood.root"
    f4   = TFile(path)
    roc4 = f4.Get("ROC_L")
    roc4.SetFillColor(1)
    roc4.SetLineColor(1)
    roc4.SetFillStyle(3001)
    
    path = outputdir3+"Top_TMVAOutput__"+alg+"__"+var0+"."+var1+"__pt"+pt1+pt2+"_ROCSMLP.root"
    f5   = TFile(path)
    roc5 = f5.Get("ROC_L")
    roc5.SetFillColor(6)
    roc5.SetLineColor(6)
    roc5.SetFillStyle(3001)
    
    path = outputdir3+"Top_TMVAOutput__"+alg+"__"+var0+"."+var1+"__pt"+pt1+pt2+"_ROCSKNN.root"
    f6   = TFile(path)
    roc6 = f6.Get("ROC_L")
    roc6.SetFillColor(9)
    roc6.SetLineColor(9)
    roc6.SetFillStyle(3001)
    
    path = outputdir3+"Top_TMVAOutput__"+alg+"__"+var0+"."+var1+"__pt"+pt1+pt2+"_ROCSBDT.root"
    f7   = TFile(path)
    roc7 = f7.Get("ROC_L")
    roc7.SetFillColor(95)
    roc7.SetLineColor(95)
    roc7.SetFillStyle(3001)
    
    cgr = TCanvas("cgr","cgr",500,500);
    gr=MakeReferenceGraph(1)
    gr.Draw("ACE3")
    

    roc1.Draw("CE3same")
    roc2.Draw("CE3same")
    roc3.Draw("CE3same")
    roc4.Draw("CE3same")
    roc5.Draw("CE3same")
    roc6.Draw("CE3same")
    roc7.Draw("CE3same")
    
    ATLASLabel(   0.20,0.90,1,0.1,0.03,"#sqrt{s}=13 TeV")
    myText(       0.20,0.85,1,0.03, TranslateAlg(alg))
    myText(       0.20,0.80,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
    rocbox1=myLineBoxText(0.26, 0.75, 2, 1, 2, 0, 0.1, 0.08, TranslateVar(var0))
    rocbox2=myLineBoxText(0.26, 0.70, 4, 1, 4, 0, 0.1, 0.08, TranslateVar(var1))
    rocbox3=myLineBoxText(0.26, 0.65, 3, 1, 1, 0, 0.1, 0.08, "Simple Likelihood (by hand)")
    rocbox4=myLineBoxText(0.26, 0.60, 1, 1, 1, 0, 0.1, 0.08, "Bayes Likelihood")
    rocbox4=myLineBoxText(0.26, 0.55, 6, 1, 1, 0, 0.1, 0.08, "Multi Layer Perceptron")
    rocbox4=myLineBoxText(0.26, 0.50, 9, 1, 1, 0, 0.1, 0.08, "K Nearest Neighbors")
    rocbox4=myLineBoxText(0.26, 0.45, 95, 1, 1, 0, 0.1, 0.08, "Boosted Decision Tree")
    cgr.SaveAs(outputdir4+"Top_"+"FullROCComparison_"+alg+"_"+var0+"_"+var1+"_pt"+pt1+pt2+".eps")


def OverlayTJetROCS(outputdir1,outputdir2,outputdir3,outputdir4,pt1,pt2,m1,m2,mvatypes,VarsAndRanges):

    algs=[]
    algs.append("TruthRawTrim")
    algs.append("CaloTrim")

    path = outputdir1+"Top_ROC_"+algs[0]+"_Tau32_pt"+pt1+pt2+".root"
    f1   = TFile(path)
    roc1 = f1.Get("ROC_L")
    roc1.SetLineColor(1)
    roc1.SetLineStyle(1)

    path = outputdir1+"Top_ROC_"+algs[1]+"_Tau32_pt"+pt1+pt2+".root"
    f2   = TFile(path)
    roc2 = f2.Get("ROC_L")
    roc2.SetLineColor(1)
    roc2.SetLineStyle(2)

    path = outputdir1+"Top_ROC_"+algs[0]+"_T3jet_minAngle_pt"+pt1+pt2+".root"
    f3   = TFile(path)
    roc3 = f3.Get("ROC_SoverB")
    roc3.SetLineColor(8)
    roc3.SetLineStyle(1)

    path = outputdir1+"Top_ROC_"+algs[1]+"_T3jet_minAngle_pt"+pt1+pt2+".root"
    f4   = TFile(path)
    roc4 = f4.Get("ROC_SoverB")
    roc4.SetLineColor(8)
    roc4.SetLineStyle(2)

    path = outputdir1 + "Top_ROC_" + algs[0] + "_v32_pt" + pt1 + pt2 + ".root"
    f5   = TFile(path)
    roc5 = f5.Get("ROC_L")
    roc5.SetLineColor(4)
    roc5.SetLineStyle(1)

    path = outputdir1 + "Top_ROC_" + algs[1] + "_v32_pt" + pt1 + pt2 + ".root"
    f6   = TFile(path)
    roc6 = f6.Get("ROC_L")
    roc6.SetLineColor(4)
    roc6.SetLineStyle(2)

    path = outputdir3+"Top_TMVAOutput__"+algs[0]+"__BDTAllTjet__pt"+pt1+pt2+"_ROCSBDT.root"
    f7   = TFile(path)
    roc7 = f7.Get("ROC_L")
    roc7.SetLineColor(95)
    roc7.SetLineStyle(1)

    path = outputdir3+"Top_TMVAOutput__"+algs[1]+"__BDTAllTjet__pt"+pt1+pt2+"_ROCSBDT.root"
    f8   = TFile(path)
    roc8 = f8.Get("ROC_L")
    roc8.SetLineColor(95)
    roc8.SetLineStyle(2)

    path = outputdir3+"Top_TMVAOutput__"+algs[0]+"__BDTMinimal__pt"+pt1+pt2+"_ROCSBDT.root"
    f9   = TFile(path)
    roc9 = f9.Get("ROC_L")
    roc9.SetLineColor(2)
    roc9.SetLineStyle(1)

    path = outputdir3+"Top_TMVAOutput__"+algs[1]+"__BDTMinimal__pt"+pt1+pt2+"_ROCSBDT.root"
    f10   = TFile(path)
    roc10 = f10.Get("ROC_L")
    roc10.SetLineColor(2)
    roc10.SetLineStyle(2)
    
    cgr = TCanvas("cgr","cgr",500,500);
    #cgr.SetLogy()
    #cgr.SetLogx()
    gr=MakeReferenceGraph(1)
    gr.Draw("ACX") #ACE3

    roc1.Draw("CXsame")
    #roc2.Draw("CXsame")
    roc3.Draw("CXsame")
    #roc4.Draw("CXsame")
    roc5.Draw("CXsame")
    #roc6.Draw("CXsame")
    roc7.Draw("CXsame")
    #roc8.Draw("CXsame")
    roc9.Draw("CXsame")
    #roc10.Draw("CXsame")
    
    ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
    #myText(       0.20,0.85,1,0.03, alg)
    myText(       0.20,0.85,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
    rocbox1=myLineBoxText(0.26, 0.80, 1, 1, 2, 0, 0.1, 0.08, TranslateVar("Tau32"))
    rocbox3=myLineBoxText(0.26, 0.76, 8, 1, 1, 0, 0.1, 0.08, TranslateVar("T3jet_angle"))
    rocbox5=myLineBoxText(0.26, 0.72, 4, 1, 1, 0, 0.1, 0.08, TranslateVar("v32"))
    rocbox7=myLineBoxText(0.26, 0.68, 95, 1, 1, 0, 0.1, 0.08, "Full Boosted Decision Tree")
    rocbox9=myLineBoxText(0.26, 0.64, 2, 1, 1, 0, 0.1, 0.08, "Minimal Boosted Decision Tree")
    cgr.SaveAs(outputdir4+"Top_"+"truthROC_tau32_v32_theta3_BDTboth_pt"+pt1+pt2+".eps")


############################
#
#
# MAIN CODE STARTS HERE
#
#
############################

flag_singlevariable  = False 
flag_2variable_hand  = False
flag_AllTjet_tmva    = True
flag_rocoverlay      = False

#==========================
#Set output directory name
#==========================
#InputDir="~/Downloads/
#InputDir="../Ana_EventGeneration/"
InputDir="/afs/cern.ch/work/a/aemerman/TJets2017/gen_20171024/"
#InputDir="/afs/cern.ch/work/a/aemerman/TJets2017/gen_20171011/"

outputdir1 = "OutputSingleVariable/"
outputdir2 = "OutputTwoVariableByHand/"
outputdir3 = "OutputTwoVariableTMVA/"
outputdir4 = "OutputROCOverlay/"

#subdir = "20171023/"
#outputdir1 += subdir 
#outputdir2 += subdir 
#outputdir3 += subdir 
#outputdir4 += subdir 

outputdir1 = MakeNewDir(outputdir1)
outputdir2 = MakeNewDir(outputdir2)
outputdir3 = MakeNewDir(outputdir3)
outputdir4 = MakeNewDir(outputdir4)

#ALGORITHMS
algs=[]
algs.append("TruthRawTrim")
#algs.append("CaloTrim")
#algs.append("TruthRaw")
#algs.append("CaloRaw")

# VARIABLES AND RANGES
VarsAndRanges={}
VarsAndRanges["Tau21"] = [0, "100, 0, 0.8", "100, 0, 0.8" ,"R"]
VarsAndRanges["Tau32"] = [0, "100, 0.1, 1", "100, 0.1, 1.0" ,"L"]
VarsAndRanges["T1jet"] = [0, "100, 0.3, 1.4", "100, 0.2, 1.0","R"]
VarsAndRanges["T2jet"] = [0, "100, 0, 0.45", "100, 0, 0.5","R"]
VarsAndRanges["T2jet_angle"]  = [0, "100, 0.275 , 1.15", "100, 0.075 , 1.0","L"]
VarsAndRanges["T3jet"] = [0, "100, 0, 0.25", "100, 0, 0.22","L"]
VarsAndRanges["T3jet_Wmass"] = [0, "100,40,120", "100,40,120","L"]
VarsAndRanges["T3jet_WmassVolatility"] = [0, "100, 0, 0.35", "100, 0, 0.25","L"]
VarsAndRanges["T3jet_minAngle"]  = [0, "100, 0, 0.9", "100, 0, 0.45","R"]
VarsAndRanges["v32"] = [0, "100, 0, 1.05", "100, 0, 1.05","L"]
VarsAndRanges["Ttrimming"] = [0, "100,0,1.5", "100,0,1.5", "L"]
VarsAndRanges["Tpruning"] = [0, "100,0,1.5", "100,0,1.5", "R"]
#VarsAndRanges["play"] = [0, "100,0,5", "100,0,5", "100,0,1.2","L"]
#VarsAndRanges["T3jet_angle1"]  = [0, "100,0,0.5", "100,0,0.5" ,"L"]
#VarsAndRanges["T3jet_angle2"]  = [0, "100,0,0.5", "100,0,0.5" ,"L"]
# VarsAndRanges["D2"]         = [0, "100,0,5", "100,0,5" ,"L"]
# VarsAndRanges["C2"]         = [0, "100,0,1", "100,0,1" ,"L"]
# VarsAndRanges["TJet_Tau21"]      = [0, "100,0,1", "100,0,1" ,"L"]
# VarsAndRanges["TJet_D2"]         = [0, "100,0,5", "100,0,5" ,"L"]
# VarsAndRanges["TJet_C2"]         = [0, "100,0,1", "100,0,1" ,"L"]

#################################
# Loop over algorithms
#################################
for alg in algs:
    print "RUNNING: ",alg

    print "\n\nGetting mass optimization"
    CutRegions=[]
    CutRegions.append("1")
    #CutRegions.append("2")

    for CutRegion in CutRegions:

        if CutRegion=="1": pt1="350"; pt2="500";   m1="160"; m2="190";
        if CutRegion=="2": pt1="800"; pt2="1000";  m1="160"; m2="190";
        ##################################################################
        # Get mass window efficiency
        ##################################################################
        SignalBGCompare1D(      InputDir, alg, "m", "100,0,200", 0, pt1, pt2, "0", "200", outputdir1)        
        effsig,effsig_err,effbkg,effbkg_err = GetMassEffs(InputDir, alg, m1, m2, outputdir1)
        print "Optimal Mass Cuts Sam:   ",m1,m2,effsig,effsig_err,effbkg,effbkg_err
    

        ##################################################################
        # ROC curves for single variables with and without mass window cut
        ##################################################################
        if flag_singlevariable:
            for var in VarsAndRanges.keys():
                print "LOOK AT ME", VarsAndRanges[var][int(CutRegion)]
                SignalBGCompare1D(InputDir, alg, var, VarsAndRanges[var][int(CutRegion)], VarsAndRanges[var][0], pt1, pt2, m1, m2, outputdir1)    

        ##################################################################
        # ROC curves for two variable combinations defined below
        ##################################################################
        TwoVarCombos = []
        #TwoVarCombos.append(["Tau21","T2jet_angle"])
        TwoVarCombos.append(["Tau32","T3jet"])

        for combo in TwoVarCombos:
        #if 1:
            var0 = combo[0]
            var1 = combo[1]

            print var0,var1

            if flag_2variable_hand:
                #######################
                #Get basic 2D correlation plots
                GetCorrelationAndSeparationSigBG(InputDir, alg, var0, var1, VarsAndRanges[var0][int(CutRegion)], VarsAndRanges[var1][int(CutRegion)], pt1, pt2, m1, m2, outputdir2)    
            
                #######################
                #Do simple likelihood combination
                SignalBGCompare2D(InputDir, alg, var0, VarsAndRanges[var0][3], var1, VarsAndRanges[var1][3], VarsAndRanges[var0][int(CutRegion)], VarsAndRanges[var1][int(CutRegion)], pt1, pt2, m1, m2, outputdir2)

            if flag_AllTjet_tmva:
                #######################
                #TMVA combination
                #mvatypes="Likelihood,MLP,BDT,KNN"
                mvatypes="BDT"
 
                tmvacommand =  " python MyTMVAClassification.py "
                tmvacommand += " "+alg+" "
                tmvacommand += alg+"_pt,"+alg+"_m"
                tmvacommand += " \"pt>"+str(pt1)+",pt<"+str(pt2)+","+alg+"_m>"+str(m1)+","+alg+"_m<"+str(m2)+"\" "
                #tmvacommand +=  " \"" + alg + "_T2jet, " + alg + "_T3jet, " + alg + "_T3jet_Wmass\" "
                tmvacommand +=  " \"" + alg + "_T3jet_WmassVolatility, " + alg + "_T2jet, " + alg + "_T3jet, " + alg + "_T2jet_angle, " + alg + "_T3jet_minAngle, " + alg + "_T3jet_Wmass\" "
                tmvacommand += " "+mvatypes+" "
                tmvacommand += " " + InputDir + sigFile + " "
                tmvacommand += " " + InputDir + bkgFile + " " 

                #outfilename=outputdir3+"Top_TMVAOutput__"+alg+"__BDTMinimal__pt"+pt1+pt2+".root"
                outfilename=outputdir3+"Top_TMVAOutput__"+alg+"__BDTAllTjet__pt"+pt1+pt2+".root"
                print "Running TMVA: ",tmvacommand
                os.system(tmvacommand)
 
                copycommand = "cp testout.root "+outfilename
                print "Copying: ",copycommand
                os.system(copycommand)
                os.system("rm testout.root")
            
                #open output root file and draw ROC curves as previously
                MakeMVAROCS(outfilename,outputdir3,mvatypes)
        
            if flag_rocoverlay:
                mvatypes="BDT"
                #Overlay all ROC curves relevant here
                #OverlayROCS(outputdir1,outputdir2,outputdir3,outputdir4,alg,var0,var1,pt1,pt2,m1,m2,mvatypes,VarsAndRanges)
                OverlayTJetROCS(outputdir1,outputdir2,outputdir3,outputdir4,pt1,pt2,m1,m2,mvatypes,VarsAndRanges) 
