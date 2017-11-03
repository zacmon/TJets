import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)

sigFileHi="ntuple_ttHi.v2.root"
bkgFileHi="ntuple_dijetHi.v2.root"
sigFileLo="ntuple_ttLo.root"
bkgFileLo="ntuple_dijetLo.root"
InputDirHi="/afs/cern.ch/work/a/aemerman/TJets2017/gen_20171011/"
InputDirLo="/afs/cern.ch/work/a/aemerman/TJets2017/gen_20171024/"

def SignalBGCompare1D(alg, cutregion, outputdir):
    '''Implementation of simple signal and background comparison'''
    VarsAndRanges={}
    VarsAndRanges["T2jet"] = [0, "100, 0, 0.45", "100, 0, 0.5","R"]
    VarsAndRanges["T3jet"] = [0, "100, 0, 0.25", "100, 0, 0.22","L"]
    VarsAndRanges["v32"] = [0, "100, 0, 1.05", "100, 0, 1.05","L"]

    if cutregion==1: pt1="350"; pt2="500";   m1="160"; m2="190"; sigFile = sigFileLo; bkgFile = bkgFileLo; InputDir = InputDirLo;
    if cutregion==2: pt1="800"; pt2="1000";  m1="160"; m2="190"; sigFile = sigFileHi; bkgFile = bkgFileHi; InputDir = InputDirHi;
    
    c = TCanvas("c","c",1500,500)
    c.SetBottomMargin(0.15);
    c.SetLeftMargin(0.15);
    c.SetTopMargin(0.05);
    c.SetRightMargin(0.05);
    p1 = TPad("p1","",0,0,0.34,1)
    p1.SetBottomMargin(0.15)
    p1.SetTopMargin(0.05)
    p1.SetLeftMargin(0.15)
    p1.SetRightMargin(0.05)
    p1.Draw()
    p2 = TPad("p2","",0.34,0,0.67,1)
    p2.SetBottomMargin(0.15)
    p2.SetTopMargin(0.05)
    p2.SetLeftMargin(0.15)
    p2.SetRightMargin(0.05)
    p2.Draw()
    p3 = TPad("p3","",0.67,0,1,1)
    p3.SetBottomMargin(0.15)
    p3.SetTopMargin(0.05)
    p3.SetLeftMargin(0.15)
    p3.SetRightMargin(0.05)
    p3.Draw()
    p1.cd()

    for variable in VarsAndRanges.keys():
        range = VarsAndRanges[variable][cutregion]
        logy = VarsAndRanges[variable][0]
        if variable == "T2jet": p1.cd(); p1.SetLogy(logy);
        elif variable == "T3jet": p2.cd(); p2.SetLogy(logy);
        elif variable == "v32": p3.cd(); p3.SetLogy(logy);

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

        hsig = GetHist1D(InputDir+sigFile, "JetTree", histname, range, weight)#+"*("+alg+"_flavor==3)")
        hbkg = GetHist1D(InputDir+bkgFile, "JetTree", histname, range, weight)#+"*("+alg+"_flavor==0)")

        hsig.Rebin()
        hbkg.Rebin()

        #Normalize them to unity
        hsig = NormalizeHist(hsig)
        hbkg = NormalizeHist(hbkg)

        #copy for line drawing
        #hbkgline = hbkg.Clone("hbkgline")
        #hsigline = hsig.Clone("hsigline")

        hbkg.SetFillStyle(0)
        hbkg.SetLineColor(2)
        hbkg.SetLineStyle(1)
        hbkg.SetLineWidth(2)
        hsig.SetFillStyle(0)
        hsig.SetLineColor(4)
        hsig.SetLineStyle(1)
        hsig.SetLineWidth(2)

        #hbkg.SetLineColor(2)
        #hbkg.SetLineWidth(4)
        #hbkg.SetLineStyle(1)
        #hbkg.SetFillColor(2)
        #hbkg.SetFillStyle(3002)
        #hbkg.SetMarkerSize(0)

        #hsig.SetLineColor(4)
        #hsig.SetLineWidth(4)
        #hsig.SetLineStyle(1)
        #hsig.SetFillColor(4)
        #hsig.SetFillStyle(3002)
        #hsig.SetMarkerSize(0)

        # DRAW
        print variable
        maxval = GetMaxVal([hbkg, hsig])
        hbkg.SetMaximum(maxval*1.5)
        hbkg.SetMinimum(0.001)

        hbkg.Draw("hist")

        #hbkg.GetYaxis().SetTitle("")
        #if variable == "T2jet": 
        hbkg.GetYaxis().SetTitle(Form("#font[42]{Normalised Entries}"))
        hbkg.GetYaxis().SetTitleOffset(1.5)
        hbkg.GetYaxis().SetLabelSize(0.04)
        hbkg.GetXaxis().SetTitle(Form("#font[42]{"+TranslateVar(variable)+"}"))
        hbkg.GetXaxis().SetTitleOffset(1.2)
        hbkg.GetXaxis().SetLabelSize(0.04)
        if variable == "T3jet": 
            hbkg.GetXaxis().SetRangeUser(0,0.205)

        #hsig.Draw("E2same")
        #hbkgline.Draw("samehist")
        hsig.Draw("samehist")

        #ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
        myText(       0.20,0.88,1,0.04, "#font[22]{Top Tagging (Pythia8)}")
        myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0, |#eta| < 1.2}")
        #myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0}")
        #myText(       0.20,0.78,1,0.04, "#font[132]{|#eta| < 1.2}")
        myText(       0.20,0.78,1,0.04, "#font[132]{"+pt1+" GeV < p_{T} < "+pt2+" GeV}")
        myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" GeV < m < "+m2+" GeV}")
        #myText(       0.20,0.88,1,0.04, "#font[22]{Pythia8,}")
        #myText(       0.35,0.88,1,0.04, "#font[132]{ #sqrt{s} = 13 TeV}")
        #myText(       0.20,0.83,1,0.04, "#font[132]{"+TranslateAlg(alg)+"}")
        #myText(       0.20,0.78,1,0.04, "#font[132]{|#eta|<1.2, "+pt1+" < p_{T} < "+pt2+" GeV}")
        #myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" < m < "+m2+" GeV}")
    #region  = "|#eta|<1.2 , p_{T}=["+pt1+","+pt2+"] GeV , m=["+m1+","+m2+"] GeV"
        rocbox3=myLineBoxText(0.75, 0.88, 4, 1, 1, 0, 0.1, 0.08, "#font[132]{Top Jets}")
        rocbox4=myLineBoxText(0.75, 0.83, 2, 1, 1, 0, 0.1, 0.08, "#font[132]{QCD Jets}")

    c.SaveAs(outputdir+"Top_SignalBGCompare_"+alg+"_pt"+pt1+pt2+".eps")

def OverlayTJetROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges):

    algs=[]
    algs.append("TruthRawTrim")
    algs.append("CaloTrim")
    CutRegions=[]
    CutRegions.append("1")
    CutRegions.append("2")

    cgr = TCanvas("cgr","cgr",1000,500)
    cgr.SetBottomMargin(0.15)
    cgr.SetTopMargin(0.05)
    cgr.SetLeftMargin(0.15)
    cgr.SetRightMargin(0.05)
    p1 = TPad("p1","",0,0,0.5,1)
    p1.SetBottomMargin(0.15)
    p1.SetTopMargin(0.05)
    p1.SetLeftMargin(0.15)
    p1.SetRightMargin(0.05)
    p1.Draw()
    p2 = TPad("p2","",0.5,0,1,1)
    p2.SetBottomMargin(0.15)
    p2.SetTopMargin(0.05)
    p2.SetLeftMargin(0.15)
    p2.SetRightMargin(0.05)
    p2.Draw()
    p1.cd()
    gr=MakeReferenceGraph(1)

    for CutRegion in CutRegions:
        if CutRegion=="1": pt1="350"; pt2="500";   m1="160"; m2="190"; p1.cd();
        if CutRegion=="2": pt1="800"; pt2="1000";  m1="160"; m2="190"; p2.cd();

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
        
        #ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
        #myText(       0.20,0.85,1,0.03, alg)
        #myText(       0.20,0.85,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
        myText(       0.20,0.88,1,0.04, "#font[22]{Top tagging (Pythia8)}")
        myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0}")
        myText(       0.20,0.78,1,0.04, "#font[132]{#sqrt{s} = 13 TeV}, |#eta| < 1.2")
        myText(       0.20,0.73,1,0.04, "#font[132]{"+pt1+" GeV < p_{T} < "+pt2+" GeV}")
        myText(       0.20,0.68,1,0.04, "#font[132]{"+m1+" GeV < m < "+m2+" GeV}")
        rocbox1=myLineBoxText(0.26, 0.63, 1, 1, 1, 0, 0.1, 0.08, TranslateVar("Tau32"))
        rocbox3=myLineBoxText(0.26, 0.43, 8, 1, 1, 0, 0.1, 0.08, TranslateVar("T3jet_angle"))
        rocbox5=myLineBoxText(0.26, 0.58, 4, 1, 1, 0, 0.1, 0.08, TranslateVar("v32"))
        rocbox7=myLineBoxText(0.26, 0.53, 95, 1, 1, 0, 0.1, 0.08, "#font[132]{BDT(6)}")
        rocbox9=myLineBoxText(0.26, 0.48, 2, 1, 1, 0, 0.1, 0.08, "#font[132]{BDT(3)}")

    #cgr.SaveAs(outputdir4+"Top_"+"ROC_tau32_v32_BDTboth.eps")
    cgr.SaveAs(outputdir4+"Top_"+"truthROC_tau32_v32_theta3_BDTboth.eps")

def OverlayThreePanelROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges):

    algs=[]
    algs.append("TruthRawTrim")
    algs.append("CaloTrim")
    CutRegions=[]
    CutRegions.append("1")
    CutRegions.append("2")

    cgr = TCanvas("cgr","cgr",1500,500)
    cgr.SetBottomMargin(0.15)
    cgr.SetTopMargin(0.05)
    cgr.SetLeftMargin(0.15)
    cgr.SetRightMargin(0.05)
    p1 = TPad("p1","",0,0,0.33,1)
    p1.SetBottomMargin(0.15)
    p1.SetTopMargin(0.05)
    p1.SetLeftMargin(0.15)
    p1.SetRightMargin(0.05)
    p1.Draw()
    p2 = TPad("p2","",0.33,0,0.66,1)
    p2.SetBottomMargin(0.15)
    p2.SetTopMargin(0.05)
    p2.SetLeftMargin(0.15)
    p2.SetRightMargin(0.05)
    p2.Draw()
    p3 = TPad("p3","",0.66,0,1,1)
    p3.SetBottomMargin(0.15)
    p3.SetTopMargin(0.05)
    p3.SetLeftMargin(0.15)
    p3.SetRightMargin(0.05)
    p3.Draw()
    p1.cd()
    gr=MakeReferenceGraph(1)

    for CutRegion in CutRegions:
        if CutRegion=="1": pt1="350"; pt2="500";   m1="160"; m2="190"; p1.cd();
        if CutRegion=="2": pt1="800"; pt2="1000";  m1="160"; m2="190"; p2.cd();

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
        roc3.SetLineColor(95)
        roc3.SetLineStyle(1)

        path = outputdir1+"Top_ROC_"+algs[1]+"_T3jet_minAngle_pt"+pt1+pt2+".root"
        f4   = TFile(path)
        roc4 = f4.Get("ROC_SoverB")
        roc4.SetLineColor(95)
        roc4.SetLineStyle(2)

        path = outputdir1 + "Top_ROC_" + algs[0] + "_v32_pt" + pt1 + pt2 + ".root"
        f5   = TFile(path)
        roc5 = f5.Get("ROC_L")
        roc5.SetLineColor(2)
        roc5.SetLineStyle(1)

        path = outputdir1 + "Top_ROC_" + algs[1] + "_v32_pt" + pt1 + pt2 + ".root"
        f6   = TFile(path)
        roc6 = f6.Get("ROC_L")
        roc6.SetLineColor(2)
        roc6.SetLineStyle(2)

        path = outputdir3+"Top_TMVAOutput__"+algs[0]+"__BDTAllTjet__pt"+pt1+pt2+"_ROCSBDT.root"
        f7   = TFile(path)
        roc7 = f7.Get("ROC_L")
        roc7.SetLineColor(8)
        roc7.SetLineStyle(1)

        path = outputdir3+"Top_TMVAOutput__"+algs[1]+"__BDTAllTjet__pt"+pt1+pt2+"_ROCSBDT.root"
        f8   = TFile(path)
        roc8 = f8.Get("ROC_L")
        roc8.SetLineColor(8)
        roc8.SetLineStyle(2)

        path = outputdir3+"Top_TMVAOutput__"+algs[0]+"__BDTMinimal__pt"+pt1+pt2+"_ROCSBDT.root"
        f9   = TFile(path)
        roc9 = f9.Get("ROC_L")
        roc9.SetLineColor(4)
        roc9.SetLineStyle(1)

        path = outputdir3+"Top_TMVAOutput__"+algs[1]+"__BDTMinimal__pt"+pt1+pt2+"_ROCSBDT.root"
        f10   = TFile(path)
        roc10 = f10.Get("ROC_L")
        roc10.SetLineColor(4)
        roc10.SetLineStyle(2)

        gr.Draw("ACX") #ACE3

        roc7.Draw("CXsame")
        #roc8.Draw("CXsame")
        roc9.Draw("CXsame")
        #roc10.Draw("CXsame")
        roc3.Draw("CXsame")
        #roc4.Draw("CXsame")
        roc5.Draw("CXsame")
        #roc6.Draw("CXsame")
        roc1.Draw("CXsame")
        #roc2.Draw("CXsame")
        
        #ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
        #myText(       0.20,0.85,1,0.03, alg)
        #myText(       0.20,0.85,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
        myText(       0.20,0.88,1,0.04, "#font[22]{Top tagging (Pythia8)}")
        myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0, |#eta| < 1.2}")
        #myText(       0.20,0.78,1,0.04, "#font[132]{|#eta| < 1.2}")
        myText(       0.20,0.78,1,0.04, "#font[132]{"+pt1+" GeV < p_{T} < "+pt2+" GeV}")
        myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" GeV < m < "+m2+" GeV}")
        #myText(       0.20,0.83,1,0.04, "#font[132]{#sqrt{s} = 13 TeV}, |#eta|<1.2")
        #myText(       0.20,0.78,1,0.04, "#font[132]{"+pt1+" < p_{T} < "+pt2+" GeV}")
        #myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" < m < "+m2+" GeV}")
        rocbox1=myLineBoxText(0.26, 0.68, 1, 1, 1, 0, 0.1, 0.08, TranslateVar("Tau32"))
        rocbox3=myLineBoxText(0.26, 0.58, 95, 1, 1, 0, 0.1, 0.08, TranslateVar("T3jet_angle"))
        rocbox5=myLineBoxText(0.26, 0.63, 2, 1, 1, 0, 0.1, 0.08, TranslateVar("v32"))
        rocbox7=myLineBoxText(0.26, 0.48, 8, 1, 1, 0, 0.1, 0.08, "#font[132]{BDT(6)}")
        rocbox9=myLineBoxText(0.26, 0.53, 4, 1, 1, 0, 0.1, 0.08, "#font[132]{BDT(3)}")

    pt1="800"; pt2="1000";  m1="160"; m2="190"; p3.cd();

    path = outputdir1+"Top_ROC_"+algs[0]+"_Tau32_pt"+pt1+pt2+".root"
    f1   = TFile(path)
    roc1 = f1.Get("ROC_L")
    roc1.SetLineColor(1)
    roc1.SetLineStyle(1)

    path = outputdir1+"Top_ROC_"+algs[1]+"_Tau32_pt"+pt1+pt2+".root"
    f2   = TFile(path)
    roc2 = f2.Get("ROC_L")
    roc2.SetLineColor(923)
    roc2.SetLineStyle(2)

    path = outputdir1 + "Top_ROC_" + algs[0] + "_v32_pt" + pt1 + pt2 + ".root"
    f5   = TFile(path)
    roc5 = f5.Get("ROC_L")
    roc5.SetLineColor(2)
    roc5.SetLineStyle(1)

    path = outputdir1 + "Top_ROC_" + algs[1] + "_v32_pt" + pt1 + pt2 + ".root"
    f6   = TFile(path)
    roc6 = f6.Get("ROC_L")
    roc6.SetLineColor(95)
    roc6.SetLineStyle(2)

    gr.Draw("ACX") #ACE3

    roc1.Draw("CXsame")
    roc5.Draw("CXsame")
    roc2.Draw("CXsame")
    roc6.Draw("CXsame")
    
    #ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
    #myText(       0.20,0.85,1,0.03, alg)
    #myText(       0.20,0.85,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
    myText(       0.20,0.88,1,0.04, "#font[22]{Top tagging (Pythia8)}")
    #myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0}")
    myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0, |#eta| < 1.2}")
    #myText(       0.20,0.78,1,0.04, "#font[132]{|#eta| < 1.2}")
    myText(       0.20,0.78,1,0.04, "#font[132]{"+pt1+" GeV < p_{T} < "+pt2+" GeV}")
    myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" GeV < m < "+m2+" GeV}")
    #myText(       0.20,0.83,1,0.04, "#font[132]{#sqrt{s} = 13 TeV}, |#eta|<1.2")
    rocbox5=myLineBoxText(0.26, 0.68, 1, 1, 1, 0, 0.1, 0.08,   "#font[132]{"+TranslateVar("Tau32")+" (truth)}")
    rocbox6=myLineBoxText(0.26, 0.63, 923, 2, 1, 0, 0.1, 0.08, "#font[132]{"+TranslateVar("Tau32")+" (cell)}")
    rocbox3=myLineBoxText(0.26, 0.58, 2, 1, 1, 0, 0.1, 0.08,  "#font[132]{"+TranslateVar("v32")+" (truth)}")
    rocbox4=myLineBoxText(0.26, 0.53, 95, 2, 1, 0, 0.1, 0.08,"#font[132]{"+TranslateVar("v32")+" (cell)}")

    #cgr.SaveAs(outputdir4+"Top_"+"ROC_tau32_v32_BDTboth.eps")
    cgr.SaveAs(outputdir4+"Top_"+"tripartROC_tau32_v32_theta3_BDTboth.eps")

############################
#
#
# MAIN CODE STARTS HERE
#
#
############################

flag_singlevariable  = True 
flag_2variable_hand  = False
flag_AllTjet_tmva    = False
flag_rocoverlay      = True

#==========================
#Set output directory name
#==========================
#InputDir="~/Downloads/
#InputDir="../Ana_EventGeneration/"

outputdir1 = "OutputSingleVariable/"
outputdir2 = "OutputTwoVariableByHand/"
outputdir3 = "OutputTwoVariableTMVA/"
outputdir4 = "OutputROCOverlay/"

subdir = "20171026/"
outputdir1 += subdir 
outputdir2 += subdir 
outputdir3 += subdir 
outputdir4 += subdir 

#outputdir1 = MakeNewDir(outputdir1)
#outputdir2 = MakeNewDir(outputdir2)
#outputdir3 = MakeNewDir(outputdir3)
#outputdir4 = MakeNewDir(outputdir4)

#ALGORITHMS
algs=[]
algs.append("TruthRawTrim")
algs.append("CaloTrim")
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
if flag_rocoverlay:
    mvatypes="BDT"
    #Overlay all ROC curves relevant here
    #OverlayTJetROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges) 
    OverlayThreePanelROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges) 
for alg in algs:
    print "RUNNING: ",alg

    if flag_singlevariable:
        SignalBGCompare1D(alg, 1, outputdir4)        
        SignalBGCompare1D(alg, 2, outputdir4)        
