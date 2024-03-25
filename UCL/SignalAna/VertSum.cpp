#include <iostream>
#include <vector>
#include <TH1.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>
#include <THStack.h>
using namespace std;

void getTotalEntriesExceptFirstBin(const std::vector<TH1*>& histograms) {
    for (TH1* histogram : histograms) {
        int totalEntries = 0;
        // Loop over all bins except the first one
        for (int bin = 0; bin <= histogram->GetNbinsX(); ++bin) {
            totalEntries += histogram->GetBinContent(bin);
        }
        std::cout << "Total entries in all bins except the first one for "
                  << histogram->GetName() << ": " << totalEntries << std::endl;
    }
}

void example() {
    // Create multiple histograms
    TFile *file = new TFile("merged_signal_histograms_nvert.root", "READ");
    TH1 *h1 = (TH1D*)file->Get("nvert");
    TH1 *h2 = (TH1D*)file->Get("nvert_chi2_cut");
    TH1 *h3 = (TH1D*)file->Get("nvert_chi2_mass");
    TH1 *h4 = (TH1D*)file->Get("nvert_chi2_mass_mom");
    TH1 *h5 = (TH1D*)file->Get("nvert_chi2_mass_mom_mep1");
    TH1 *h6 = (TH1D*)file->Get("nvert_TDR");

    std::vector<TH1*> histograms = {h1,h2,h3,h4,h5,h6};
    // Call the function to get the total entries in all bins except the first one for each histogram
    getTotalEntriesExceptFirstBin(histograms);
}

int main() {
    example();
    return 0;
}