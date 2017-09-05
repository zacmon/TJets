import os
import sys
import itertools
from ROOT import *
from MyLocalFunctions import *
from AtlasStyle import *
SetAtlasStyle();
gStyle.SetPalette(1)


#where are the data stored
InputDir="../Ana_EventGeneration/GenNTuple/20170428/"

#where will you put the plots
outputdir = "myPlots/"
outputdir = MakeNewDir(outputdir)

c = TCanvas("c","c",200,200)

variable = "TruthRaw_TJet_m1"
range    = "100,0,4"

h1 = GetHist1D(InputDir+"ntuple_ww_0.root",    "JetTree", variable, range, "TruthRawTrim_m>60")
h2 = GetHist1D(InputDir+"ntuple_dijet_0.root", "JetTree", variable, range, "TruthRawTrim_m>60")

#Normalize them to unity
h1 = NormalizeHist(h1)
h2 = NormalizeHist(h2)

h2line = h2.Clone("h2line")
h1line = h1.Clone("h1line")

h2line.SetFillStyle(0)
h2line.SetLineColor(2)
h2line.SetLineStyle(1)
h2line.SetLineWidth(3)
h1line.SetFillStyle(0)
h1line.SetLineColor(4)
h1line.SetLineStyle(1)
h1line.SetLineWidth(3)


# DRAW
h2.GetXaxis().SetTitle("Mass [GeV]")
h2.GetXaxis().SetTitleOffset(1.2)
h2.GetYaxis().SetTitleOffset(1.7)
h2.GetYaxis().SetTitle("Normalised Entries")
h2.SetLineColor(2)
h2.SetLineWidth(4)
h2.SetLineStyle(2)
h2.GetXaxis().SetTitle("Mass [GeV]")
h2.GetYaxis().SetTitleOffset(1.6)
h2.SetFillColor(2)
h2.SetLineColor(2)
h2.SetLineStyle(1)
h2.SetFillStyle(3001)
h2.SetMarkerSize(0)

h1.SetLineColor(4)
h1.SetLineWidth(4)
h1.GetXaxis().SetTitle("Mass [GeV]")
h1.GetYaxis().SetTitleOffset(1.6)
h1.SetFillColor(4)
h1.SetLineColor(4)
h1.SetLineStyle(1)
h1.SetFillStyle(3002)
h1.SetMarkerSize(0)


maxval = GetMaxVal([h2, h1])
h2.SetMaximum(maxval*2.0)
h2.SetMinimum(0.001)

h2.Draw("E2")
h2line.Draw("samehist")
h1.Draw("E2same")
h1line.Draw("samehist")

rocbox3=myLineBoxText(0.70, 0.85, 4, 1, 1, 0, 0.1, 0.08, "W jets")
rocbox4=myLineBoxText(0.70, 0.80, 2, 1, 1, 0, 0.1, 0.08, "QCD jets")
c.SaveAs(outputdir+"validation.eps")
