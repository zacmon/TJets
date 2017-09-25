#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TTree.h>
#include <THStack.h>

bool existsFile(const char* pathname) {
    struct stat info;
    return stat(pathname, &info) == 0 && S_ISREG(info.st_mode);
}

int markerstyle[] = {20, 21, 22, 34, 33, 39, 41};
int colorarray[] = {1, 2, 4, 3, 6, 7, 9};
int lineStyleArray[] = {1, 2, 5, 7, 3, 8, 10};

void overlayHistograms(std::vector<TH1D*> histograms, TString fileName, TString title, int log=0) {
    THStack* thstack = new THStack("thstack", "");
    TLegend* legendOverlay = new TLegend(0.99, 0.98, 0.80, 0.85);
    legendOverlay->SetTextSize(0.03);
	
    for (unsigned int i(0); i < histograms.size(); ++i) {
	std::cout << "Hist name: " << histograms[i]->GetName() << "\tEntries: " << histograms[i]->GetEntries() << std::endl;

	histograms[i]->SetTitle(title);
	histograms[i]->SetMarkerStyle(markerstyle[i]);
	histograms[i]->SetMarkerColor(colorarray[i]);
	histograms[i]->SetLineColor(colorarray[i]);
	histograms[i]->SetLineStyle(1);//lineStyleArray[i]);
	histograms[i]->SetMarkerSize(2);
        histograms[i]->SetLineWidth(3);
	histograms[i]->SetStats(0);

	thstack->Add(histograms[i]);
	
        legendOverlay -> AddEntry(histograms[i], histograms[i]->GetName(), "l");
    }

    TCanvas* canvasOverlay = new TCanvas("canvasOverlay", "canvasOverlay", 1500, 1500);
    thstack->Draw("nostack");
    legendOverlay->Draw();

    thstack->GetXaxis()->SetTitle(title);
    thstack->GetYaxis()->SetTitle("# of events");
    
    if (log == 1) gPad->SetLogy(1);

    canvasOverlay->SaveAs(fileName);
    std::cout << "Successfully saved!\n" << std::endl;
}

void ComparePythiaWCommands() {
    TString dir = "../Ana_EventGeneration/";
    
    TFile* topWithWCmnd = new TFile(dir + TString("ntuple_tt_test10000.root"), "READ");
    TTree* topWithWCmndTree = (TTree*)topWithWCmnd->Get("JetTree");
    
    TFile* topNoWCmnd = new TFile(dir + TString("ntuple_tt_test10000_NO_W_CMND.root"), "READ");
    TTree* topNoWCmndTree = (TTree*)topNoWCmnd->Get("JetTree");
    
    TFile* dijet = new TFile(dir + TString("ntuple_dijet_test10000.root"), "READ");
    TTree* dijetTree = (TTree*)dijet->Get("JetTree");

    
    TString algorithm = "TruthRawTrim_";
    std::vector<TString> histNames = {"eta", "phi", "m",
				      "Tau21", "Tau32",
				      "T1jet",
				      "T2jet_angle", "T2jet",
				      "T3jet_angle", "T3jet", "T3jet_W", "T3jet_mW"};
    
    for (auto const &histName : histNames) {
	TString name = algorithm + histName;
	std::cout << name << std::endl;
	std::vector<TH1D*> histograms;
	TString saveName = name + ".png";
	
	topWithWCmndTree->Draw(name + ">>h1", "", "goff");
	TH1D *h1 = (TH1D*)gDirectory->Get("h1");
	h1->SetName("topWithWCmnd");
	topNoWCmndTree->Draw(name + ">>h2", "", "goff");
	TH1D *h2 = (TH1D*)gDirectory->Get("h2");
	h2->SetName("topNoWCmnd");
	dijetTree->Draw(name + ">>h3", "", "goff");
	TH1D *h3 = (TH1D*)gDirectory->Get("h3");
	h3->SetName("dijet");
	overlayHistograms({h1, h2, h3}, name + TString(".png"), name);
    }   
}
