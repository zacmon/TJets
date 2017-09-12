import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)


outputDir = "OutputTwoVariableTMVA/20170912/"
preFileName = "TMVAOutput__TruthRawTrim__BDTAll"
postFileName = "Tjet__pt8001000_ROCSBDT.root"
rocs = []

for x in range(7):
    file = TFile(outputDir + preFileName + str(x) + postFileName)
    roc = file.Get("ROC_SoverB")
    roc.SetFillColor(x + 1)
    roc.SetLineColor(x + 1)
    roc.SetLineWidth(1)
    roc.SetFillStyle(3001)
    rocs.append(roc)

canvas = TCanvas("c", "c", 500, 500)
graph = MakeReferenceGraph(1)
graph.Draw("ACE3")

for roc in rocs:
    roc.Draw("CE3same")

ATLASLabel(0.20,0.90,1,0.1,0.03,"#sqrt{s}=13 TeV")

rocbox1=myLineBoxText(0.26, 0.75, 1, 1, 2, 0, 0.1, 0.08, "All attributes")
rocbox2=myLineBoxText(0.26, 0.70, 2, 1, 1, 0, 0.1, 0.08, "No T3jet_mW")
rocbox3=myLineBoxText(0.26, 0.65, 3, 1, 4, 0, 0.1, 0.08, "No T2jet")
rocbox4=myLineBoxText(0.26, 0.60, 4, 1, 1, 0, 0.1, 0.08, "No T2jet_angle")
rocbox5=myLineBoxText(0.26, 0.55, 5, 1, 1, 0, 0.1, 0.08, "No T3jet")
rocbox6=myLineBoxText(0.26, 0.50, 6, 1, 1, 0, 0.1, 0.08, "No T3jet_angle")
rocbox7=myLineBoxText(0.26, 0.45, 7, 1, 1, 0, 0.1, 0.08, "No T3jet_W")
canvas.SaveAs("~/BDTOVERLAY.eps")
