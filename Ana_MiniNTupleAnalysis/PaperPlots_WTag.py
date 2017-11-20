import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)

#sigFile="GenNTuple/20171015/ntuple_wwLowPt_0.root"
#bkgFile="GenNTuple/20171015/ntuple_dijetLowPt_0.root"
#sigFile="GenNTuple/20171016/ntuple_ww_0.root"
#bkgFile="GenNTuple/20171016/ntuple_dijet_0.root"
sigFile="ntuple_wwLo.v2.root"
bkgFile="ntuple_dijetLo.v2.root"

def OverlayCaloROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges):

    #ALGORITHMS
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
        if CutRegion=="1": pt1="350"; pt2="500";   m1="70"; m2="90"; p1.cd();
        if CutRegion=="2": pt1="800"; pt2="1000";  m1="70"; m2="90"; p2.cd();

        path = outputdir1+"W_ROC_"+algs[0]+"_T3jet_pt"+pt1+pt2+".root"
        f3   = TFile(path)
        roc3 = f3.Get("ROC_L")
        roc3.SetLineColor(4)
        roc3.SetLineStyle(1)

        path = outputdir1+"W_ROC_"+algs[1]+"_T3jet_pt"+pt1+pt2+".root"
        f4   = TFile(path)
        roc4 = f4.Get("ROC_L")
        roc4.SetLineColor(4)
        roc4.SetLineStyle(2)

        path = outputdir1+"W_ROC_"+algs[0]+"_Tau21_pt"+pt1+pt2+".root"
        f5   = TFile(path)
        roc5 = f5.Get("ROC_L")
        roc5.SetLineColor(1)
        roc5.SetLineStyle(1)

        path = outputdir1+"W_ROC_"+algs[1]+"_Tau21_pt"+pt1+pt2+".root"
        f6   = TFile(path)
        roc6 = f6.Get("ROC_L")
        roc6.SetLineColor(1)
        roc6.SetLineStyle(2)

        gr.Draw("ACX") #ACE3

        roc3.Draw("CXsame")
        roc4.Draw("CXsame")
        roc5.Draw("CXsame")
        roc6.Draw("CXsame")
        
        #ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
        #myText(       0.20,0.85,1,0.03, alg)
        #myText(       0.20,0.85,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
        myText(       0.20,0.88,1,0.04, "#font[22]{W tagging (Pythia8)}")
        myText(       0.20,0.83,1,0.04, "#font[132]{|#eta|<1.2, "+pt1+" < p_{T} < "+pt2+" GeV}")
        myText(       0.20,0.78,1,0.04, "#font[132]{"+m1+" < m < "+m2+" GeV}")
        #myText(       0.20,0.83,1,0.04, "#font[132]{#sqrt{s} = 13 TeV}, |#eta|<1.2")
        #myText(       0.20,0.78,1,0.04, "#font[132]{"+pt1+" < p_{T} < "+pt2+" GeV}")
        #myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" < m < "+m2+" GeV}")
        rocbox5=myLineBoxText(0.26, 0.73, 1, 1, 1, 0, 0.1, 0.08, "#font[132]{"+TranslateVar("Tau21")+"}")
        rocbox3=myLineBoxText(0.26, 0.68, 4, 1, 1, 0, 0.1, 0.08, "#font[132]{"+TranslateVar("T3jet")+"}")

    cgr.SaveAs(outputdir4+"W_"+"ROC_tau21_v3.eps")
    #cgr.SaveAs(outputdir4+"W_"+"truthROC_tau21_v3_Ttrim_Tprun_BDT.eps")

def OverlayTruthROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges):

    #ALGORITHMS
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
        if CutRegion=="1": pt1="350"; pt2="500";   m1="70"; m2="90"; p1.cd();
        if CutRegion=="2": pt1="800"; pt2="1000";  m1="70"; m2="90"; p2.cd();

        path = outputdir1+"W_ROC_"+algs[0]+"_T3jet_pt"+pt1+pt2+".root"
        f3   = TFile(path)
        roc3 = f3.Get("ROC_L")
        roc3.SetLineColor(4)
        roc3.SetLineStyle(1)

        path = outputdir1+"W_ROC_"+algs[1]+"_T3jet_pt"+pt1+pt2+".root"
        f4   = TFile(path)
        roc4 = f4.Get("ROC_L")
        roc4.SetLineColor(4)
        roc4.SetLineStyle(2)

        path = outputdir1+"W_ROC_"+algs[0]+"_Tau21_pt"+pt1+pt2+".root"
        f5   = TFile(path)
        roc5 = f5.Get("ROC_L")
        roc5.SetLineColor(1)
        roc5.SetLineStyle(1)

        path = outputdir1+"W_ROC_"+algs[1]+"_Tau21_pt"+pt1+pt2+".root"
        f6   = TFile(path)
        roc6 = f6.Get("ROC_L")
        roc6.SetLineColor(1)
        roc6.SetLineStyle(2)

        path = outputdir3+"W_TMVAOutput__"+algs[0]+"__BDTAllTjet__pt"+pt1+pt2+"_ROCSBDT.root"
        f7   = TFile(path)
        roc7 = f7.Get("ROC_L")
        roc7.SetLineColor(8) #95
        roc7.SetLineStyle(1)

        path = outputdir3+"W_TMVAOutput__"+algs[1]+"__BDTAllTjet__pt"+pt1+pt2+"_ROCSBDT.root"
        f8   = TFile(path)
        roc8 = f8.Get("ROC_L")
        roc8.SetLineColor(8) #95
        roc8.SetLineStyle(2)
        
        path = outputdir1+"W_ROC_TruthRaw_Tpruning_pt"+pt1+pt2+".root"
        f9   = TFile(path)
        roc9 = f9.Get("ROC_L")
        roc9.SetLineColor(2)
        roc9.SetLineStyle(1)
        
        path = outputdir1+"W_ROC_CaloRaw_Tpruning_pt"+pt1+pt2+".root"
        f10   = TFile(path)
        roc10 = f10.Get("ROC_L")
        roc10.SetLineColor(2)
        roc10.SetLineStyle(2)

        path = outputdir1+"W_ROC_TruthRaw_Ttrimming_pt"+pt1+pt2+".root"
        f11   = TFile(path)
        roc11 = f11.Get("ROC_L")
        roc11.SetLineColor(95)
        roc11.SetLineStyle(1)
        
        path = outputdir1+"W_ROC_CaloRaw_Ttrimming_pt"+pt1+pt2+".root"
        f12   = TFile(path)
        roc12 = f12.Get("ROC_L")
        roc12.SetLineColor(95)
        roc12.SetLineStyle(2)

        path = outputdir1+"W_ROC_"+algs[0]+"_v21_pt"+pt1+pt2+".root"
        f13   = TFile(path)
        roc13 = f13.Get("ROC_L")
        roc13.SetLineColor(9)
        roc13.SetLineStyle(1)

        path = outputdir1+"W_ROC_"+algs[1]+"_v21_pt"+pt1+pt2+".root"
        f14   = TFile(path)
        roc14 = f14.Get("ROC_L")
        roc14.SetLineColor(9)
        roc14.SetLineStyle(2)


        gr.Draw("ACX") #ACE3

        roc3.Draw("CXsame")
        #roc4.Draw("CXsame")
        roc5.Draw("CXsame")
        #roc6.Draw("CXsame")
        roc7.Draw("CXsame")
        #roc8.Draw("CXsame")
        roc9.Draw("CXsame")
        #roc10.Draw("CXsame")
        roc11.Draw("CXsame")
        #roc12.Draw("CXsame")
        roc13.Draw("CXsame")
        #roc14.Draw("CXsame")
        
        #ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
        #myText(       0.20,0.85,1,0.03, alg)
        #myText(       0.20,0.85,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
        myText(       0.20,0.88,1,0.04, "#font[22]{W tagging (Pythia8)}")
        myText(       0.20,0.83,1,0.04, "#font[132]{|#eta|<1.2, "+pt1+" < p_{T} < "+pt2+" GeV}")
        myText(       0.20,0.78,1,0.04, "#font[132]{"+m1+" < m < "+m2+" GeV}")
        #myText(       0.20,0.83,1,0.04, "#font[132]{#sqrt{s} = 13 TeV}, |#eta|<1.2")
        #myText(       0.20,0.78,1,0.04, "#font[132]{"+pt1+" < p_{T} < "+pt2+" GeV}")
        #myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" < m < "+m2+" GeV}")
        rocbox5=myLineBoxText(0.26,  0.73, 1, 1, 1, 0, 0.1, 0.08, "#font[132]{"+TranslateVar("Tau21")+"}")
        rocbox3=myLineBoxText(0.26,  0.68, 4, 1, 1, 0, 0.1, 0.08, "#font[132]{"+TranslateVar("T3jet")+"}")
        rocbox13=myLineBoxText(0.26, 0.63, 9, 1, 1, 0, 0.1, 0.08, "#font[132]{"+TranslateVar("v21")+"}")
        rocbox7=myLineBoxText(0.26,  0.58, 8, 1, 1, 0, 0.1, 0.08, "#font[132]{BDT(3)}")
        rocbox9=myLineBoxText(0.26,  0.53, 2, 1, 1, 0, 0.1, 0.08, "#font[132]{"+TranslateVar("Tpruning")+"}")
        rocbox11=myLineBoxText(0.26, 0.48, 95, 1, 1, 0, 0.1, 0.08, "#font[132]{"+TranslateVar("Ttrimming")+"}")

    #cgr.SaveAs(outputdir4+"W_"+"ROC_tau21_v3_Ttrim_Tprun_BDT.eps")
    cgr.SaveAs(outputdir4+"W_"+"truthROC_tau21_v3_v21_Ttrim_Tprun_BDT.eps")

def OverlayGroomingROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges):

    #ALGORITHMS
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
        if CutRegion=="1": pt1="350"; pt2="500";   m1="70"; m2="90"; p1.cd();
        if CutRegion=="2": pt1="800"; pt2="1000";  m1="70"; m2="90"; p2.cd();

        path = outputdir1+"W_ROC_"+algs[0]+"_T3jet_pt"+pt1+pt2+".root"
        f3   = TFile(path)
        roc3 = f3.Get("ROC_L")
        roc3.SetLineColor(1)
        roc3.SetLineStyle(1)
        
        path = outputdir1+"W_ROC_TruthRaw_Tpruning_pt"+pt1+pt2+".root"
        f9   = TFile(path)
        roc9 = f9.Get("ROC_L")
        roc9.SetLineColor(4)
        roc9.SetLineStyle(2)

        path = outputdir1+"W_ROC_TruthRaw_Ttrimming_pt"+pt1+pt2+".root"
        f11   = TFile(path)
        roc11 = f11.Get("ROC_L")
        roc11.SetLineColor(2)
        roc11.SetLineStyle(4)

        gr.Draw("ACX") #ACE3

        roc9.Draw("CXsame")
        roc11.Draw("CXsame")
        roc3.Draw("CXsame")
        
        #ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
        #myText(       0.20,0.85,1,0.03, alg)
        #myText(       0.20,0.85,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
        myText(       0.20,0.88,1,0.04, "#font[22]{W tagging (Pythia8)}")
        myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0, |#eta| < 1.2}")
        #myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0}")
        #myText(       0.20,0.78,1,0.04, "#font[132]{|#eta| < 1.2}")
        myText(       0.20,0.78,1,0.04, "#font[132]{"+pt1+" GeV < p_{T} < "+pt2+" GeV}")
        myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" GeV < m < "+m2+" GeV}")
        rocbox3=myLineBoxText(0.30,  0.68, 1, 1, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("T3jet")+"}")
        rocbox9=myLineBoxText(0.30,  0.63, 4, 2, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("Tpruning")+"}")
        rocbox11=myLineBoxText(0.30, 0.58, 2, 4, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("Ttrimming")+"}")

    cgr.SaveAs(outputdir4+"W_"+"truthROC_v3_Ttrim_Tprun.eps")

def OverlayThreePanelROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges):

    #ALGORITHMS
    algs=[]
    algs.append("TruthRawTrim")
    algs.append("CaloTrim")
    CutRegions=[]
    CutRegions.append("1")
    CutRegions.append("2")

    cgr = TCanvas("cgr","cgr",1500,500)
    #cgr.SetGrayscale()
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
        if CutRegion=="1": pt1="350"; pt2="500";   m1="70"; m2="90"; p1.cd();
        if CutRegion=="2": pt1="800"; pt2="1000";  m1="70"; m2="90"; p2.cd();

        path = outputdir1+"W_ROC_"+algs[0]+"_T3jet_pt"+pt1+pt2+".root"
        f3   = TFile(path)
        roc3 = f3.Get("ROC_L")
        roc3.SetLineColor(2)
        roc3.SetLineStyle(1)

        path = outputdir1+"W_ROC_"+algs[0]+"_Tau21_pt"+pt1+pt2+".root"
        f5   = TFile(path)
        roc5 = f5.Get("ROC_L")
        roc5.SetLineColor(1)
        roc5.SetLineStyle(1)

        path = outputdir3+"W_TMVAOutput__"+algs[0]+"__BDTAllTjet__pt"+pt1+pt2+"_ROCSBDT.root"
        f7   = TFile(path)
        roc7 = f7.Get("ROC_L")
        roc7.SetLineColor(4)
        roc7.SetLineStyle(4)
        
        path = outputdir1+"W_ROC_TruthRaw_Tpruning_pt"+pt1+pt2+".root"
        f9   = TFile(path)
        roc9 = f9.Get("ROC_L")
        roc9.SetLineColor(95)
        roc9.SetLineStyle(1)
        
        path = outputdir1+"W_ROC_TruthRaw_Ttrimming_pt"+pt1+pt2+".root"
        f11   = TFile(path)
        roc11 = f11.Get("ROC_L")
        roc11.SetLineColor(417)
        roc11.SetLineStyle(3)
        
        path = outputdir1+"W_ROC_"+algs[0]+"_T2jet_pt"+pt1+pt2+".root"
        #path = outputdir1+"W_ROC_"+algs[0]+"_v21_pt"+pt1+pt2+".root"
        f13   = TFile(path)
        roc13 = f13.Get("ROC_L")
        roc13.SetLineColor(617)
        roc13.SetLineStyle(2)


        gr.Draw("ACX") #ACE3

        roc9.Draw("CXsame")
        roc11.Draw("CXsame")
        roc7.Draw("CXsame")
        roc13.Draw("CXsame")
        roc3.Draw("CXsame")
        roc5.Draw("CXsame")
        
        #ATLASLabel(   0.20,0.88,1,0.1,0.03,"#sqrt{s}=13 TeV")
        #myText(       0.20,0.85,1,0.03, alg)
        #myText(       0.20,0.85,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
        myText(       0.20,0.88,1,0.04, "#font[22]{W tagging (Pythia8)}")
        myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0, |#eta| < 1.2}")
        #myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0}")
        #myText(       0.20,0.78,1,0.04, "#font[132]{|#eta| < 1.2}")
        myText(       0.20,0.78,1,0.04, "#font[132]{"+pt1+" GeV < p_{T} < "+pt2+" GeV}")
        myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" GeV < m < "+m2+" GeV}")
        #myText(       0.20,0.83,1,0.04, "#font[132]{#sqrt{s} = 13 TeV}, |#eta|<1.2")
        rocbox5= myLineBoxText(0.30, 0.68, 1, 1, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("Tau21")+"}")
        rocbox3= myLineBoxText(0.30, 0.63, 2, 1, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("T3jet")+"}")
        rocbox13=myLineBoxText(0.30, 0.58, 617, 2, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("T2jet")+"}")
        #rocbox13=myLineBoxText(0.30, 0.58, 95, 2, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("v21")+"}")
        rocbox7= myLineBoxText(0.30, 0.53, 4, 4, 1, 0, 0.15, 0.08, "#font[132]{BDT(3)}")
        rocbox9= myLineBoxText(0.30, 0.48, 95, 1, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("Tpruning")+"}")
        rocbox11=myLineBoxText(0.30, 0.43, 417, 3, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("Ttrimming")+"}")

    pt1="800"; pt2="1000";  m1="70"; m2="90"; p3.cd();

    path = outputdir1+"W_ROC_"+algs[0]+"_T3jet_pt"+pt1+pt2+".root"
    f3   = TFile(path)
    roc3 = f3.Get("ROC_L")
    roc3.SetLineColor(2)
    roc3.SetLineStyle(1)

    path = outputdir1+"W_ROC_"+algs[1]+"_T3jet_pt"+pt1+pt2+".root"
    f4   = TFile(path)
    roc4 = f4.Get("ROC_L")
    roc4.SetLineColor(2)#95
    roc4.SetLineStyle(2)

    path = outputdir1+"W_ROC_"+algs[0]+"_Tau21_pt"+pt1+pt2+".root"
    f5   = TFile(path)
    roc5 = f5.Get("ROC_L")
    roc5.SetLineColor(1)
    roc5.SetLineStyle(1)

    path = outputdir1+"W_ROC_"+algs[1]+"_Tau21_pt"+pt1+pt2+".root"
    f6   = TFile(path)
    roc6 = f6.Get("ROC_L")
    roc6.SetLineColor(1)#923
    roc6.SetLineStyle(2)

    gr.Draw("ACX") #ACE3

    roc3.Draw("CXsame")
    roc5.Draw("CXsame")
    roc4.Draw("CXsame")
    roc6.Draw("CXsame")

    myText(       0.20,0.88,1,0.04, "#font[22]{W tagging (Pythia8)}")
    myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0, |#eta| < 1.2}")
    #myText(       0.20,0.83,1,0.04, "#font[132]{Anti-k_{T} R = 1.0}")
    #myText(       0.20,0.78,1,0.04, "#font[132]{|#eta| < 1.2}")
    myText(       0.20,0.78,1,0.04, "#font[132]{"+pt1+" GeV < p_{T} < "+pt2+" GeV}")
    myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" GeV < m < "+m2+" GeV}")
    #myText(       0.20,0.83,1,0.04, "#font[132]{#sqrt{s} = 13 TeV}, |#eta|<1.2")
    #myText(       0.20,0.73,1,0.04, "#font[132]{"+m1+" < m < "+m2+" GeV}")
    rocbox5=myLineBoxText(0.30, 0.68, 1, 1, 1, 0, 0.15, 0.08,   "#font[132]{"+TranslateVar("Tau21")+" (truth)}")
    rocbox6=myLineBoxText(0.30, 0.63, 1, 2, 1, 0, 0.15, 0.08, "#font[132]{"+TranslateVar("Tau21")+" (cell)}")
    rocbox3=myLineBoxText(0.30, 0.58, 2, 1, 1, 0, 0.15, 0.08,  "#font[132]{"+TranslateVar("T3jet")+" (truth)}")
    rocbox4=myLineBoxText(0.30, 0.53, 2, 2, 1, 0, 0.15, 0.08,"#font[132]{"+TranslateVar("T3jet")+" (cell)}")

    #cgr.SaveAs(outputdir4+"W_"+"tripartROC_tau21_v3_v21_Ttrim.eps")
    #cgr.SaveAs(outputdir4+"W_"+"tripartROC_tau21_v3_v2_Ttrim_BDT.eps")
    #cgr.SaveAs(outputdir4+"W_"+"tripartROC_tau21_v3_v2_BDT.eps")
    #cgr.SaveAs(outputdir4+"W_"+"tripartROC_tau21_v3_v2_Ttrim_BDT_grayscale.eps")
    cgr.SaveAs(outputdir4+"W_"+"tripartROC_tau21_v3_v2_Ttrim_Tprun_BDT.eps")

############################
#
#
# MAIN CODE STARTS HERE
#
#
############################

flag_singlevariable  = False     
flag_2variable_hand  = False
flag_AllTjet_tmva    = False
flag_rocoverlay      = True

#==========================
#Set output directory name
#==========================
#InputDir="~/Downloads/
#InputDir="../Ana_EventGeneration/"
InputDir="/afs/cern.ch/work/a/aemerman/TJets2017/gen_20171011/"

outputdir1 = "OutputSingleVariable/"
outputdir2 = "OutputTwoVariableByHand/"
outputdir3 = "OutputTwoVariableTMVA/"
outputdir4 = "OutputROCOverlay/"

subdir = "20171102/"
outputdir1 += subdir 
outputdir2 += subdir 
outputdir3 += subdir 
outputdir4 += subdir

#ALGORITHMS
algs=[]
algs.append("TruthRawTrim")
algs.append("CaloTrim")
#algs.append("TruthRaw")
#algs.append("CaloRaw")

CutRegions=[]
CutRegions.append("1")
CutRegions.append("2")

# VARIABLES AND RANGES
VarsAndRanges={}
VarsAndRanges["Tau21"] = [0, "100, 0, 0.9", "100, 0, 1" ,"R"]
VarsAndRanges["Tau32"] = [0, "100, 0.2, 1", "100, 0, 1" ,"L"]
VarsAndRanges["T1jet"] = [0, "100, 0.2, 0.9", "100, 0, 0.5","R"]
VarsAndRanges["T2jet"] = [0, "100, 0, 0.45", "100, 0, 0.3","R"]
VarsAndRanges["T2jet_angle"] = [0, "100, 0, 1.0", "100, 0, 0.45","L"]
VarsAndRanges["T3jet"] = [0, "100, 0, 0.21", "100, 0, 0.25","L"]
VarsAndRanges["T3jet_Wmass"] = [0, "100,45,90", "100,40,90","L"]
VarsAndRanges["T3jet_WmassVolatility"] = [0, "100, 0, 0.25", "100, 0, 0.25","L"]
VarsAndRanges["T3jet_minAngle"] = [0, "100, 0, 0.4", "100, 0, 0.3","R"]
VarsAndRanges["v32"] = [0, "100, 0, 1.05", "100, 0, 1.05","L"]
VarsAndRanges["Ttrimming"] = [0, "100,0,2", "300,0,1", "L"]
VarsAndRanges["Tpruning"] = [0, "100,0,2", "300,0,1.1", "L"]
VarsAndRanges["v21"] = [0, "100,0,5", "100,0,5", "100,0,1.2","L"]
#VarsAndRanges["T3jet_angle1"]  = [0, "100,0,0.5", "100,0,0.5" ,"L"]
#VarsAndRanges["T3jet_angle2"]  = [0, "100,0,0.5", "100,0,0.5" ,"L"]
# VarsAndRanges["D2"]         = [0, "100,0,5", "100,0,5" ,"L"]
# VarsAndRanges["C2"]         = [0, "100,0,1", "100,0,1" ,"L"]
# VarsAndRanges["TJet_Tau21"]      = [0, "100,0,1", "100,0,1" ,"L"]
# VarsAndRanges["TJet_D2"]         = [0, "100,0,5", "100,0,5" ,"L"]
# VarsAndRanges["TJet_C2"]         = [0, "100,0,1", "100,0,1" ,"L"]

mvatypes="BDT"
#Overlay all ROC curves relevant here
#OverlayTruthROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges) 
#OverlayCaloROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges) 
OverlayThreePanelROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges) 
#OverlayGroomingROCS(outputdir1,outputdir2,outputdir3,outputdir4,mvatypes,VarsAndRanges) 
