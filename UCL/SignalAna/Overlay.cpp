#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>
#include <THStack.h>
using namespace std;

int main(){
    TFile *file = new TFile("merged_signal_histograms2.root", "READ");
    TH1 *h1 = (TH1D*)file->Get("mmass_4hit_1");
    TH1 *h2 = (TH1D*)file->Get("mmass_4hit_2");
    TH1 *h3 = (TH1D*)file->Get("mmass_4hit_3");


h1->SetLineColor(kRed);
h2->SetLineColor(kBlue);
h3->SetLineColor(kGreen);
double h1mean = h1->GetMean(1);
double h1std = h1->GetStdDev(1);
double h2mean = h2->GetMean(1);
double h2std = h2->GetStdDev(1);
double h3mean = h3->GetMean(1);
double h3std = h3->GetStdDev(1);

    TCanvas *c1 = new TCanvas();
c1->cd();
h1->DrawCopy();
h2->DrawCopy("SAME");
h3->DrawCopy("SAME");


    TLegend *legend = new TLegend(0.5,0.7,1.0,1.0);
    legend->SetHeader("","C"); // option "C" allows to center the header
    legend->AddEntry(h1,Form("Reconstructed mass with 4 hits allowed for particle 1, #mu = %.2f MeV, #sigma = %.2f MeV", h1mean, h1std),"l");
    legend->AddEntry(h2,Form("Reconstructed mass with 4 hits allowed for particle 2, #mu = %.2f MeV, #sigma = %.2f MeV", h2mean, h2std),"l");
    legend->AddEntry(h3,Form("Reconstructed mass with 4 hits allowed for particle 3, #mu = %.2f MeV, #sigma = %.2f MeV",h3mean,h3std),"l");
    legend->Draw();


c1->Print("Mass_with_4hits_allowed_on_1_track.pdf");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TH1 *h4 = (TH1D*)file->Get("mmass_4hit_12");
    TH1 *h5 = (TH1D*)file->Get("mmass_4hit_13");
    TH1 *h6 = (TH1D*)file->Get("mmass_4hit_23");

//h4->Fit("gaus","WW");
//h5->Fit("gaus","WW");
//h6->Fit("gaus","WW");
//h4->GetFunction("gaus")->SetLineColor(1);
//h5->GetFunction("gaus")->SetLineColor(2);
//h6->GetFunction("gaus")->SetLineColor(3);
h4->SetLineColor(kOrange);
h5->SetLineColor(kViolet);
h6->SetLineColor(kBlack);
double h4mean = h4->GetMean(1);
double h4std = h4->GetStdDev(1);
double h5mean = h5->GetMean(1);
double h5std = h5->GetStdDev(1);
double h6mean = h6->GetMean(1);
double h6std = h6->GetStdDev(1);

    TCanvas *c2 = new TCanvas();
c2->cd();
h4->DrawCopy();
h5->DrawCopy("SAME");
h6->DrawCopy("SAME");


    TLegend *legend2 = new TLegend(0.5,0.7,1.0,1.0);
    legend2->SetHeader("","C"); // option "C" allows to center the header
    legend2->AddEntry(h4,Form("Reconstructed mass with 4 hits allowed for particles 1 and 2, #mu = %.2f MeV, #sigma = %.2f MeV", h4mean, h4std),"l");
    legend2->AddEntry(h5,Form("Reconstructed mass with 4 hits allowed for particles 1 and 3, #mu = %.2f, #sigma = %.2f MeV", h5mean, h5std),"l");
    legend2->AddEntry(h6,Form("Reconstructed mass with 4 hits allowed for particles 2 and 3, #mu = %.2f, #sigma = %.2f MeV", h6mean, h6std),"l");
    legend2->Draw();

c2->Print("Mass_with_4hits_allowed_on_2_tracks.pdf");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TH1 *h7 = (TH1D*)file->Get("mmass_4hit_123");

h7->SetLineColor(kCyan);
double h7mean = h7->GetMean(1);
double h7std = h7->GetStdDev(1);

    TCanvas *c3 = new TCanvas();
c3->cd();
h7->DrawCopy();
h1->DrawCopy("SAME");
h2->DrawCopy("SAME");
h3->DrawCopy("SAME");
h4->DrawCopy("SAME");
h5->DrawCopy("SAME");
h6->DrawCopy("SAME");


    TLegend *legend3 = new TLegend(0.5,0.7,1.0,1.0);
    legend3->SetHeader("","C"); // option "C" allows to center the header
    legend3->AddEntry(h1,Form("4 hits allowed for track 1, #mu = %.2f MeV, #sigma = %.2f MeV", h1mean, h1std),"l");
    legend3->AddEntry(h2,Form("4 hits allowed for track 2, #mu = %.2f MeV, #sigma = %.2f MeV", h2mean, h2std),"l");
    legend3->AddEntry(h3,Form("4 hits allowed for track 3, #mu = %.2f MeV, #sigma = %.2f MeV", h3mean, h3std),"l");
    legend3->AddEntry(h4,Form("4 hits allowed for tracks 1 and 2, #mu = %.2f MeV, #sigma = %.2f MeV", h4mean, h4std),"l");
    legend3->AddEntry(h5,Form("4 hits allowed for tracks 1 and 3, #mu = %.2f, #sigma = %.2f MeV", h5mean, h5std),"l");
    legend3->AddEntry(h6,Form("4 hits allowed for tracks 2 and 3, #mu = %.2f, #sigma = %.2f MeV", h6mean, h6std),"l");
    legend3->AddEntry(h7,Form("4 hits allowed for tracks 1, 2 and 3, #mu = %.2f, #sigma = %.2f MeV", h7mean, h7std),"l");
    legend3->Draw();

c3->Print("Overlaid_mass_with_4_hits.pdf");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

THStack pm("pm","");

    TH1 *h8 = (TH1D*)file->Get("pm_4hit_1");
    TH1 *h9 = (TH1D*)file->Get("pm_4hit_2");
    TH1 *h10 = (TH1D*)file->Get("pm_4hit_3");
    TH1 *h11 = (TH1D*)file->Get("pm_4hit_12");
    TH1 *h12 = (TH1D*)file->Get("pm_4hit_13");
    TH1 *h13 = (TH1D*)file->Get("pm_4hit_23");
    TH1 *h14 = (TH1D*)file->Get("pm_4hit_123");

h8->SetLineColor(kRed);
h9->SetLineColor(kBlue);
h10->SetLineColor(kGreen);
h11->SetLineColor(kOrange);
h12->SetLineColor(kViolet);
h13->SetLineColor(kBlack);
h14->SetLineColor(kCyan);
double h8mean = h8->GetMean(1);
double h8std = h8->GetStdDev(1);
double h9mean = h9->GetMean(1);
double h9std = h9->GetStdDev(1);
double h10mean = h10->GetMean(1);
double h10std = h10->GetStdDev(1);
double h11mean = h11->GetMean(1);
double h11std = h11->GetStdDev(1);
double h12mean = h12->GetMean(1);
double h12std = h12->GetStdDev(1);
double h13mean = h13->GetMean(1);
double h13std = h13->GetStdDev(1);
double h14mean = h14->GetMean(1);
double h14std = h14->GetStdDev(1);

pm.Add(h8);
pm.Add(h9);
pm.Add(h10);
pm.Add(h11);
pm.Add(h12);
pm.Add(h13);
pm.Add(h14);


    TCanvas *c4 = new TCanvas();

c4->cd();
pm.Draw("nostack");
//h8->DrawCopy();
//h9->DrawCopy("SAME");
//h10->DrawCopy("SAME");
//h11->DrawCopy("SAME");
//h12->DrawCopy("SAME");
//h13->DrawCopy("SAME");
//h14->DrawCopy("SAME");


    TLegend *legend4 = new TLegend(0.5,0.7,1.0,1.0);
    legend4->SetHeader("","C"); // option "C" allows to center the header
    legend4->AddEntry(h8,Form("4 hits allowed for track 1, #mu = %.2f MeV, #sigma = %.2f MeV", h8mean, h8std),"l");
    legend4->AddEntry(h9,Form("4 hits allowed for track 2, #mu = %.2f MeV, #sigma = %.2f MeV", h9mean, h9std),"l");
    legend4->AddEntry(h10,Form("4 hits allowed for track 3, #mu = %.2f MeV, #sigma = %.2f MeV", h10mean, h10std),"l");
    legend4->AddEntry(h11,Form("4 hits allowed for tracks 1 and 2, #mu = %.2f MeV, #sigma = %.2f MeV", h11mean, h11std),"l");
    legend4->AddEntry(h12,Form("4 hits allowed for tracks 1 and 3, #mu = %.2f, #sigma = %.2f MeV", h12mean, h12std),"l");
    legend4->AddEntry(h13,Form("4 hits allowed for tracks 2 and 3, #mu = %.2f, #sigma = %.2f MeV", h13mean, h13std),"l");
    legend4->AddEntry(h14,Form("4 hits allowed for tracks 1, 2 and 3, #mu = %.2f, #sigma = %.2f MeV", h14mean, h14std),"l");
    legend4->Draw();


c4->Print("Overlaid_momentum_with_4_hits.pdf");


cout<<"Overlaid plots have been written"<<endl;

}
