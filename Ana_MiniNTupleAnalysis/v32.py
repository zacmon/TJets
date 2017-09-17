import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)

InputDir = "../Ana_EventGeneration/"
sigFile="ntuple_tt_50000v32.root"
bkgFile="ntuple_dijet_50000v32.root"

variable = "v32"
alg = "TruthRawTrim"
pt1 = "800"
pt2 = "1000"
m1 = "160"
m2 = "190"

range = ""
logy=0
outputdir="~/"

'''Implementation of simple signal and background comparison'''
print "Making 1D comparison: ",alg,variable

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
histname = alg+"_"+variable
if variable == "v32":
    histname = variable
    
print histname
hsig = GetHist1D(InputDir+sigFile, "JetTree", histname, range, weight+"*("+alg+"_flavor==3)")
hbkg = GetHist1D(InputDir+bkgFile, "JetTree", histname, range, weight+"*("+alg+"_flavor==0)")

#Normalize them to unity
hsig = NormalizeHist(hsig)
hbkg = NormalizeHist(hbkg)

rocL,hsigregL,hcutvalL,hsigregL25,hcutvalL25 = RocCurve_SingleSided_WithUncer(hsig, hbkg, "L")
rocR,hsigregR,hcutvalR,hsigregR25,hcutvalR25 = RocCurve_SingleSided_WithUncer(hsig, hbkg, "R")
rocSB,sigordered,bkgordered,h1 = RocCurve_SoverBOrdered_WithUncer(hsig, hbkg)

print alg+"_"+variable
f = TFile(outputdir+"ROC_"+alg+"_"+variable.replace("\\","")+"_pt"+pt1+pt2+".root","RECREATE")
rocL.Write("ROC_L")
rocR.Write("ROC_R")
rocSB.Write("ROC_SoverB")
f.Close()

c.cd()
rocSB.Draw("AC*")
c.SaveAs(outputdir+"ROCDraw_"+alg+"_"+variable.replace("\\","")+"_pt"+pt1+pt2+".eps")

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

print variable
hbkg.GetXaxis().SetTitle(variable)
hbkg.GetXaxis().SetTitleOffset(1.2)
hbkg.GetYaxis().SetTitleOffset(1.7)
hbkg.GetYaxis().SetTitle("Normalised Entries")
hbkg.SetLineColor(2)
hbkg.SetLineWidth(4)
hbkg.SetLineStyle(2)
hbkg.GetXaxis().SetTitle(variable)
hbkg.GetYaxis().SetTitleOffset(1.6)
hbkg.SetFillColor(2)
hbkg.SetLineColor(2)
hbkg.SetLineStyle(1)
hbkg.SetFillStyle(3001)
hbkg.SetMarkerSize(0)

hsig.SetLineColor(4)
hsig.SetLineWidth(4)
hsig.GetXaxis().SetTitle(variable)
hsig.GetYaxis().SetTitleOffset(1.6)
hsig.SetFillColor(4)
hsig.SetLineColor(4)
hsig.SetLineStyle(1)
hsig.SetFillStyle(3002)
hsig.SetMarkerSize(0)

maxval = GetMaxVal([hbkg, hsig])
hbkg.SetMaximum(maxval*2.0)
hbkg.SetMinimum(0.001)

hbkg.Draw("E2")
hbkgline.Draw("samehist")
hsig.Draw("E2same")
hsigline.Draw("samehist")

ATLASLabel(   0.20,0.85,1,0.1,0.03,"#sqrt{s}=13 TeV")
myText(       0.20,0.80,1,0.03, TranslateAlg(alg))
myText(       0.20,0.75,1,0.03, TranslateRegion(pt1,pt2,m1,m2))
rocbox3=myLineBoxText(0.70, 0.85, 4, 1, 1, 0, 0.1, 0.08, "Top Jets")
rocbox4=myLineBoxText(0.70, 0.80, 2, 1, 1, 0, 0.1, 0.08, "QCD Jets")
c.SetLogy(logy)
c.SaveAs(outputdir+"SignalBGCompare_"+alg+"_"+variable+"_pt"+pt1+pt2+".eps")
