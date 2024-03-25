#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>
#include <TStyle.h>
#include <string>


using namespace std;

int main(int argc, char *argv[]){

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <signal_file> <IC_file> <bhabha_file> " << std::endl;
        return 1;
    }

    TFile *signal = new TFile(argv[1],"READ");
    TFile *IC = new TFile(argv[2],"READ");
    TFile *bhabha = new TFile(argv[3],"READ");

    cout<<"Files Read"<<endl;

    gStyle->SetOptStat(0);
    cout<<"Style set"<<endl;

    TH1 *s1 = (TH1D*)signal->Get("mmass_TDRcuts");
    TH1 *ic1 = (TH1D*)IC->Get("mmass_TDRcuts");
    TH1 *b1 = (TH1D*)bhabha->Get("mmass_TDRcuts");
    cout<<"Mass plots read"<<endl;

    ic1->SetLineColor(kBlue);
    
    s1->SetLineColor(kRed);
    
    b1->SetLineColor(kGreen);
    

    double s1mean = s1->GetMean(1);
    double s1std = s1->GetStdDev(1);
    double ic1mean = ic1->GetMean(1);
    double ic1std = ic1->GetStdDev(1);
    double b1mean = b1->GetMean(1);
    double b1std = b1->GetStdDev(1);
    cout<<"Mass statistics calculated"<<endl;

        TCanvas *c_mmass = new TCanvas();
    c_mmass->cd();
    s1->Scale(0.01);

    s1->DrawCopy();

    ic1->DrawCopy("SAME");

    b1->DrawCopy("SAME");


    TLegend *legend_mmass = new TLegend(0.5,0.7,1.0,1.0);
    legend_mmass->SetHeader("","C"); // option "C" allows to center the header
    legend_mmass->AddEntry(s1,Form("m_{rec, signal} (scaled by 1/100) #mu = %.2f MeV, #sigma = %.2f MeV", s1mean, s1std),"l");
    legend_mmass->AddEntry(ic1,Form("m_{rec, IC} #mu = %.2f MeV, #sigma = %.2f MeV", ic1mean, ic1std),"l");
    legend_mmass->AddEntry(b1,Form("m_{rec, Bhabha} #mu = %.2f MeV, #sigma = %.2f MeV",b1mean,b1std),"l");
    legend_mmass->Draw();


    c_mmass->Print("mmass_TDRcuts.pdf");


    TH1 *s2 = (TH1D*)signal->Get("pm_TDRcuts");
    TH1 *ic2 = (TH1D*)IC->Get("pm_TDRcuts");
    TH1 *b2 = (TH1D*)bhabha->Get("pm_TDRcuts");
    cout<<"Momentum plots read"<<endl;

    ic2->SetLineColor(kBlue);
    
    s2->SetLineColor(kRed);
    
    b2->SetLineColor(kGreen);
    

    double s2mean = s2->GetMean(1);
    double s2std = s2->GetStdDev(1);
    double ic2mean = ic2->GetMean(1);
    double ic2std = ic2->GetStdDev(1);
    double b2mean = b2->GetMean(1);
    double b2std = b2->GetStdDev(1);
    cout<<"Momentum statistics calculated"<<endl;

        TCanvas *c_pm = new TCanvas();
    c_pm->cd();
    s2->Scale(0.01);
    
    s2->DrawCopy();

    ic2->DrawCopy("SAME");

    b2->DrawCopy("SAME");


    TLegend *legend_pm = new TLegend(0.5,0.7,1.0,1.0);
    legend_pm->SetHeader("","C"); // option "C" allows to center the header
    legend_pm->AddEntry(s2,Form("pm_{rec, signal} (scaled by 1/100) #mu = %.2f MeV, #sigma = %.2f MeV", s2mean, s2std),"l");
    legend_pm->AddEntry(ic2,Form("pm_{rec, IC} #mu = %.2f MeV, #sigma = %.2f MeV", ic2mean, ic2std),"l");
    legend_pm->AddEntry(b2,Form("pm_{rec, Bhabha} #mu = %.2f MeV, #sigma = %.2f MeV",b2mean,b2std),"l");
    legend_pm->Draw();


    c_pm->Print("pm_TDRcuts.pdf");



}