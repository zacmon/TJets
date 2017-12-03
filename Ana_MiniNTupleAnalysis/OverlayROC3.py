import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)


scale="Linear"
outputDir1 = "OutputSingleVariable/20171120/"
outputDir2 = "OutputSingleVariable/20171128/LinearScale/"
alg = "CaloTrim"

file1 = TFile(outputDir1 + "ROC_"+alg+"_tau32_pt8001000.root")
roc1 = file1.Get("ROC_L")
roc1.SetFillColor(1)
roc1.SetLineColor(1)
roc1.SetLineWidth(3)
roc1.SetFillStyle(3001)

file2 = TFile(outputDir1 + "ROC_"+alg+"_v32_pt8001000.root")
roc2 = file2.Get("ROC_L")
roc2.SetFillColor(7)
roc2.SetLineColor(7)
roc2.SetLineWidth(3)
roc2.SetFillStyle(3001)

file3 = TFile(outputDir1 + "ROC_"+alg+"_v42_pt8001000.root")
roc3 = file3.Get("ROC_L")
roc3.SetFillColor(3)
roc3.SetLineColor(3)
roc3.SetLineWidth(3)
roc3.SetFillStyle(3001)

file4 = TFile(outputDir2 + "ROC_"+alg+"_tau32_pt8001000.root")
roc4 = file4.Get("ROC_L")
roc4.SetFillColor(4)
roc4.SetLineColor(4)
roc4.SetLineWidth(1)
roc4.SetFillStyle(3001)

file5 = TFile(outputDir2 + "ROC_"+alg+"_v32_pt8001000.root")
roc5 = file5.Get("ROC_L")
roc5.SetFillColor(95)
roc5.SetLineColor(95)
roc5.SetLineWidth(3)
roc5.SetFillStyle(3001)

file6 = TFile(outputDir2 + "ROC_"+alg+"_v42_pt8001000.root")
roc6= file6.Get("ROC_L")
roc6.SetFillColor(28)
roc6.SetLineColor(28)
roc6.SetLineWidth(3)
roc6.SetFillStyle(3001)

canvas = TCanvas("c", "c", 500, 500)
graph = MakeReferenceGraph(1)
graph.Draw("ACE3")

roc1.Draw("Csame")
roc2.Draw("Csame")
roc3.Draw("Csame")
roc4.Draw("Csame")
roc5.Draw("Csame")
roc6.Draw("Csame")

ATLASLabel(   0.70,0.90,1,0.1,0.03, alg)
rocbox1=myLineBoxText(0.26, 0.75, 1, 1, 2, 0, 0.1, 0.08, "1-pass #tau_{32}")
rocbox2=myLineBoxText(0.26, 0.70, 7, 1, 1, 0, 0.1, 0.08, "1-pass v_{32}")
rocbox3=myLineBoxText(0.26, 0.65, 3, 1, 2, 0, 0.1, 0.08, "1-pass v_{42}")
rocbox4=myLineBoxText(0.26, 0.60, 4, 1, 1, 0, 0.1, 0.08, "WTA #tau_{32}")
rocbox5=myLineBoxText(0.26, 0.55, 95, 1, 1, 0, 0.1, 0.08, "WTA pt_{32}")
rocbox6=myLineBoxText(0.26, 0.50, 28, 1, 2, 0, 0.1, 0.08, "WTA pt_{42}")
canvas.SaveAs("~/BDTOVERLAY1.eps")
