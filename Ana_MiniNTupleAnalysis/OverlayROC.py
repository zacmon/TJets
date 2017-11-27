import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)


outputDir = "OutputTwoVariableTMVA/20170922/"
preFileName = "TMVAOutput__TruthRawTrim__BDTAll"
postFileName = "Tjet__pt8001000_ROCSBDT.root"
rocs = []
 
# for x in range(7):
#     file = TFile(outputDir + preFileName + str(x) + postFileName)
#     roc = file.Get("ROC_SoverB")
#     roc.SetFillColor(x + 1)
#     roc.SetLineColor(x + 1)
#     roc.SetLineWidth(1)
#     roc.SetFillStyle(3001)
#     rocs.append(roc)

file = TFile("OutputTwoVariableTMVA/20171120/TMVAOutput__TruthRawTrim__BDTAllTjet__pt8001000_7inputs_ROCSBDT.root")
roc1 = file.Get("ROC_L")
roc1.SetFillColor(1)
roc1.SetLineColor(1)
roc1.SetLineWidth(1)
roc1.SetFillStyle(3001)

file = TFile("OutputTwoVariableTMVA/20171120/TMVAOutput__TruthRawTrim__BDTAllTjet__pt8001000_6inputs_ROCSBDT.root")
roc2 = file.Get("ROC_L")
roc2.SetFillColor(2)
roc2.SetLineColor(2)
roc2.SetLineWidth(1)
roc2.SetFillStyle(3001)

file = TFile("OutputTwoVariableTMVA/20171120/TMVAOutput__TruthRawTrim__BDTAllTjet__pt8001000_v2-v3-mW2_ROCSBDT.root")
roc3 = file.Get("ROC_L")
roc3.SetFillColor(3)
roc3.SetLineColor(3)
roc3.SetLineWidth(1)
roc3.SetFillStyle(3001)

file = TFile("OutputTwoVariableTMVA/20171120/TMVAOutput__TruthRawTrim__BDTAllTjet__pt8001000_v2-v4_ROCSBDT.root")
roc4 = file.Get("ROC_L")
roc4.SetFillColor(4)
roc4.SetLineColor(4)
roc4.SetLineWidth(1)
roc4.SetFillStyle(3001)

file = TFile("OutputTwoVariableTMVA/20171120/TMVAOutput__TruthRawTrim__BDTAllTjet__pt8001000_v2-v3_ROCSBDT.root")
roc5 = file.Get("ROC_L")
roc5.SetFillColor(28)
roc5.SetLineColor(28)
roc5.SetLineWidth(1)
roc5.SetFillStyle(3001)

file = TFile("OutputSingleVariable/20171120/ROC_TruthRawTrim_v32_pt8001000.root")
roc6 = file.Get("ROC_L")
roc6.SetFillColor(6)
roc6.SetLineColor(6)
roc6.SetLineWidth(1)
roc6.SetFillStyle(3001)

file = TFile("OutputSingleVariable/20171120/ROC_TruthRawTrim_v42_pt8001000.root")
roc7 = file.Get("ROC_L")
roc7.SetFillColor(7)
roc7.SetLineColor(7)
roc7.SetLineWidth(1)
roc7.SetFillStyle(3001)

file = TFile("OutputSingleVariable/20171120/ROC_TruthRawTrim_Tau32_pt8001000.root")
roc8 = file.Get("ROC_L")
roc8.SetFillColor(8)
roc8.SetLineColor(8)
roc8.SetLineWidth(1)
roc8.SetFillStyle(3001)

file = TFile("OutputTwoVariableTMVA/20171120/TMVAOutput__TruthRawTrim__BDTAllTjet__pt8001000_v2-v3-v4_ROCSBDT.root")
roc9 = file.Get("ROC_L")
roc9.SetFillColor(95)
roc9.SetLineColor(95)
roc9.SetLineWidth(1)
roc9.SetFillStyle(3001)

canvas = TCanvas("c", "c", 500, 500)
graph = MakeReferenceGraph(1)
graph.Draw("ACE3")

# for roc in rocs:
#     roc.Draw("CE3same")
roc1.Draw("Csame")
roc2.Draw("Csame")
roc3.Draw("Csame")
roc4.Draw("Csame")
roc5.Draw("Csame")
roc6.Draw("Csame")
roc7.Draw("Csame")
roc8.Draw("Csame")
roc9.Draw("Csame")

rocbox1=myLineBoxText(0.26, 0.75, 1, 1, 2, 0, 0.1, 0.08, "{mW2, v_{mW2}, #theta_{2} #theta_{min}, v_{2}, v_{3}, v_{4}}")
rocbox2=myLineBoxText(0.26, 0.70, 2, 1, 1, 0, 0.1, 0.08, "{mW2, v_{mW2}, #theta_{2} #theta_{min}, v_{2}, v_{3}}")
rocbox3=myLineBoxText(0.26, 0.65, 3, 1, 4, 0, 0.1, 0.08, "{v2, v3, mW2}")
rocbox4=myLineBoxText(0.26, 0.60, 4, 1, 4, 0, 0.1, 0.08, "{v2, v4}")
rocbox5=myLineBoxText(0.26, 0.55, 28, 1, 4, 0, 0.1, 0.08, "{v2, v3}")
rocbox6=myLineBoxText(0.26, 0.50, 6, 1, 1, 0, 0.1, 0.08, "v_{32}")
rocbox7=myLineBoxText(0.26, 0.45, 7, 1, 1, 0, 0.1, 0.08, "v_{42}")
rocbox8=myLineBoxText(0.26, 0.40, 8, 1, 1, 0, 0.1, 0.08, "#tau_{32}")
rocbox9=myLineBoxText(0.26, 0.35, 95, 1, 1, 0, 0.1, 0.08, "{v2, v3, v4}")
canvas.SaveAs("~/BDTOVERLAY.eps")
