import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)

outputDir = "OutputSingleVariable/20170925/"

file1 = TFile(outputDir + "ROC_TruthRawTrim_v32_pt8001000.root")
roc1 = file1.Get("ROC_L")
roc1.SetFillColor(7)
roc1.SetLineColor(7)
roc1.SetLineWidth(5)
roc1.SetFillStyle(3001)

file2 = TFile(outputDir + "ROC_TruthRawTrim_Tau32_pt8001000.root")
roc2 = file2.Get("ROC_L")
roc2.SetFillColor(2)
roc2.SetLineColor(2)
roc2.SetLineWidth(5)
roc2.SetFillStyle(3001)

canvas = TCanvas("c", "c", 500, 500)
graph = MakeReferenceGraph(1)
graph.Draw("ACE3")

roc1.Draw("CE3same")
roc2.Draw("CE3same")

ATLASLabel(0.20,0.90,1,0.1,0.03,"#sqrt{s}=13 TeV")

rocbox1=myLineBoxText(0.26, 0.75, 2, 1, 2, 0, 0.1, 0.08, "#tau_{32}")
rocbox2=myLineBoxText(0.26, 0.70, 7, 1, 1, 0, 0.1, 0.08, "T_{32}")
canvas.SaveAs("~/BDTOVERLAY1.eps")
