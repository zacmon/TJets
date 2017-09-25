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
int colorarray[] = {1, 2, 4, 8, 95, 6, 7, 9};
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
        histograms[i]->SetLineWidth(5);
	histograms[i]->SetStats(0);

	thstack->Add(histograms[i]);
	
        legendOverlay -> AddEntry(histograms[i], histograms[i]->GetName(), "l");
    }

    TCanvas* canvasOverlay = new TCanvas("canvasOverlay", "canvasOverlay", 1500, 1500);
    thstack->Draw("nostack");
    legendOverlay->Draw();

    thstack->GetXaxis()->SetTitle("m_{jj}");
    thstack->GetYaxis()->SetTitle("# of events / 2 GeV");
    
    if (log == 1) gPad->SetLogy(1);

    canvasOverlay->SaveAs(fileName);
    std::cout << "Successfully saved!\n" << std::endl;
}

TH1D* getHistogramShape(TH1D* hist) {
    double normalTerm = hist->Integral();
    hist->Scale(1/normalTerm);
    return hist;
}

void CompareVolatility() {
    TString dir = "../Ana_EventGeneration/";
    
    TFile* file = new TFile(dir + TString("ntuple_tt_1.root"), "READ");
    TTree* tree = (TTree*)file->Get("JetTree");
        
    TString algorithm = "TruthRawTrim_";
    std::vector<TString> histNames = {"T2masses[0]"};

    TString cuts = "(TruthRawTrim_m>160&&TruthRawTrim_m<190&&TruthRawTrim_pt>800&&TruthRawTrim_m<1000)";
    
    for (auto const &histName : histNames) {
	TString name = algorithm + histName;
	std::cout << name << std::endl;
	std::vector<TH1D*> histograms;
	TString saveName = name + ".png";
	
	tree->Draw(name + "[0]>>h1(100,50,250)", cuts);
	TH1D *h1 = (TH1D*)gDirectory->Get("h1");
	std::cout << h1->GetName() << std::endl;
	h1->SetName("R=0.1");
	tree->Draw(name + "[5]>>h2(100,50,250)", cuts);
	TH1D *h2 = (TH1D*)gDirectory->Get("h2");
	h2->SetName("R=0.25");

	tree->Draw(name + "[11]>>h3(100,50,250)", cuts);
	TH1D *h3 = (TH1D*)gDirectory->Get("h3");
	h3->SetName("R=0.4");

	tree->Draw(name + "[17]>>h4(100,50,250)", cuts);
	TH1D *h4 = (TH1D*)gDirectory->Get("h4");
	h4->SetName("R=0.55");

	tree->Draw(name + "[35]>>h5(100,50,250)", cuts);
	TH1D *h5 = (TH1D*)gDirectory->Get("h5");
	h5->SetName("R=1.0");

	overlayHistograms({h1, h2, h3, h4, h5}, name + TString(".eps"), name);
    }   
}
