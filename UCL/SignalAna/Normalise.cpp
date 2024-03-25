#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>
#include <TStyle.h>
#include <string>
#include <TMath.h>
#include <THStack.h>


int main(){
// IC Normalisation
float GenEff=2.486041e-06;
float BrIC=3.4e-5;
float NmuStops=2.5e+15;
//double ICsumOfWeights = 2.1378e+12;
double ICsumOfWeights = 44737066;
float IC_event_norm = (NmuStops*BrIC*GenEff)/ICsumOfWeights;

using namespace std;

cout<<"The normalisation scaling for internal conversion is "<<IC_event_norm<<endl;

// Signal Normalisation

double N_signal = 8.76223e+07;
float Br1 = 1e-12;
float Br2 = 1e-13;
float Br3 = 1e-14;
float Br4 = 1e-15;
float signal_event_norm1=(NmuStops*Br1)/N_signal;
float signal_event_norm2=(NmuStops*Br2)/N_signal;
float signal_event_norm3=(NmuStops*Br3)/N_signal;
float signal_event_norm4=(NmuStops*Br4)/N_signal;

cout<<"The normalisation scaling for the signal with B = 1e-12 is "<<signal_event_norm1<<endl;
cout<<"The normalisation scaling for the signal with B = 1e-13 is "<<signal_event_norm2<<endl;
cout<<"The normalisation scaling for the signal with B = 1e-14 is "<<signal_event_norm3<<endl;
cout<<"The normalisation scaling for the signal with B = 1e-15 is "<<signal_event_norm4<<endl;

//Bhabha Normalisation

float bhabha_frac = 7.7e-5;
float timing_suppression = 70;
float num_total_frames=174062382;
float avg_bhabha_weight_per_frame=4.1224e-06;
float event_norm_BB = (NmuStops*bhabha_frac)/(num_total_frames*avg_bhabha_weight_per_frame*timing_suppression);

cout<<"The normalisation scaling for Bhabha is "<<event_norm_BB<<endl;


TFile *signal = new TFile("merged_signal_histograms_nvert.root", "READ");
TFile *IC = new TFile("merged_IC_histograms2.root", "READ");
TFile *bhabha = new TFile("merged_bhabha_histograms2.root", "READ");

gStyle->SetOptStat(0);

TH1 *mmass_signal1 = (TH1D*)signal->Get("mmass_TDR_to_normalise");
mmass_signal1->SetName("signal_mass");
TH1 *mmass_signal2 = (TH1D*)signal->Get("mmass_TDR_to_normalise");
mmass_signal2->SetName("signal_mass");
TH1 *mmass_signal3 = (TH1D*)signal->Get("mmass_TDR_to_normalise");
mmass_signal3->SetName("signal_mass");
TH1 *mmass_signal4 = (TH1D*)signal->Get("mmass_TDR_to_normalise");
mmass_signal4->SetName("signal_mass");


TH1 *mmass_IC = (TH1D*)IC->Get("mmass_TDR_without_weights");
mmass_IC->SetName("IC_mass");


TH1 *mmass_bhabha = (TH1D*)bhabha->Get("mmass_TDR_to_normalise");
mmass_bhabha->SetName("bhabha_mass");

mmass_signal1->SetLineColor(kRed);
mmass_signal2->SetLineColor(kRed);
mmass_signal3->SetLineColor(kRed);
mmass_signal4->SetLineColor(kRed);
mmass_IC->SetLineColor(kBlue);
mmass_bhabha->SetLineColor(kGreen);

double mmass_signal1_mean = mmass_signal1->GetMean(1);
double mmass_signal1_std = mmass_signal1->GetStdDev(1);
double mmass_signal2_mean = mmass_signal2->GetMean(1);
double mmass_signal2_std = mmass_signal2->GetStdDev(1);
double mmass_signal3_mean = mmass_signal3->GetMean(1);
double mmass_signal3_std = mmass_signal3->GetStdDev(1);
double mmass_signal4_mean = mmass_signal4->GetMean(1);
double mmass_signal4_std = mmass_signal4->GetStdDev(1);
double mmass_IC_mean = mmass_IC->GetMean(1);
double mmass_IC_std= mmass_IC->GetStdDev(1);
double mmass_bhabha_mean = mmass_bhabha->GetMean(1);
double mmass_bhabha_std= mmass_bhabha->GetStdDev(1);


TCanvas *c_mmass1 = new TCanvas();
    c_mmass1->cd();
    c_mmass1->SetLogy();
    mmass_signal1->Scale(signal_event_norm1);
    mmass_signal2->Scale(signal_event_norm2);
    mmass_signal3->Scale(signal_event_norm3);
    mmass_signal4->Scale(signal_event_norm4);
    mmass_IC->Scale(IC_event_norm);
    mmass_bhabha->Scale(event_norm_BB);
    mmass_signal1->DrawCopy();
    mmass_signal2->DrawCopy("SAME");
    mmass_signal3->DrawCopy("SAME");
    mmass_signal4->DrawCopy("SAME");
    mmass_IC->DrawCopy("SAME");
    mmass_bhabha->DrawCopy("SAME");


TLegend *legend_mmass = new TLegend(0.5,0.7,1.0,1.0);
    legend_mmass->SetHeader("","C"); // option "C" allows to center the header
    legend_mmass->AddEntry(mmass_signal1,Form("m_{rec, signal} (Br = #10^{-12}) #mu = %.2f MeV, #sigma = %.2f MeV", mmass_signal1_mean, mmass_signal1_std),"l");
    legend_mmass->AddEntry(mmass_signal2,Form("m_{rec, signal} (Br = #10^{-13}) #mu = %.2f MeV, #sigma = %.2f MeV", mmass_signal2_mean, mmass_signal2_std),"l");
    legend_mmass->AddEntry(mmass_signal3,Form("m_{rec, signal} (Br = #10^{-14}) #mu = %.2f MeV, #sigma = %.2f MeV", mmass_signal3_mean, mmass_signal3_std),"l");
    legend_mmass->AddEntry(mmass_signal4,Form("m_{rec, signal} (Br = #10^{-15}) #mu = %.2f MeV, #sigma = %.2f MeV", mmass_signal4_mean, mmass_signal4_std),"l");
    legend_mmass->AddEntry(mmass_IC,Form("m_{rec, IC} #mu = %.2f MeV, #sigma = %.2f MeV", mmass_IC_mean, mmass_IC_std),"l");
    legend_mmass->AddEntry(mmass_bhabha,Form("m_{rec, Bhabha} #mu = %.2f MeV, #sigma = %.2f MeV", mmass_bhabha_mean, mmass_bhabha_std),"l");
    //legend_mmass->Draw();


    c_mmass1->Print("mmass_TDRcuts_normalised.pdf");
    

THStack pm("pm","stacked pm histograms");


TH1 *pm_signal1 = (TH1D*)signal->Get("pm_TDR_to_normalise1");

TH1 *pm_signal2 = (TH1D*)signal->Get("pm_TDR_to_normalise2");

TH1 *pm_signal3 = (TH1D*)signal->Get("pm_TDR_to_normalise3");

TH1 *pm_signal4 = (TH1D*)signal->Get("pm_TDR_to_normalise4");

TH1 *pm_IC = (TH1D*)IC->Get("pm_TDR_without_weights");

TH1 *pm_bhabha = (TH1D*)bhabha->Get("pm_TDR_to_normalise");


pm_signal1->SetLineColor(kRed);
pm_signal2->SetLineColor(kRed);
pm_signal3->SetLineColor(kRed);
pm_signal4->SetLineColor(kRed);
pm_IC->SetLineColor(kBlue);
pm_bhabha->SetLineColor(kGreen);

//TCanvas *c_pm = new TCanvas();
    //c_pm->cd();
    //c_pm->SetLogy();
    pm_signal1->Scale(signal_event_norm1);
    pm_signal2->Scale(signal_event_norm2);
    pm_signal3->Scale(signal_event_norm3);
    pm_signal4->Scale(signal_event_norm4);
    pm_IC->Scale(IC_event_norm);
    pm_bhabha->Scale(event_norm_BB);

    pm.Add(pm_signal1);
    pm.Add(pm_signal2);
    pm.Add(pm_signal3);
    pm.Add(pm_signal4);
    pm.Add(pm_IC);
    pm.Add(pm_bhabha);
    pm.SetMaximum(100.);

    TCanvas *c_pm = new TCanvas();
    
    c_pm->cd();
    c_pm->SetLogy();
    pm.Draw("nostack");
    //pm_signal1->DrawCopy();
    //pm_signal2->DrawCopy();
    //pm_signal3->DrawCopy("SAME");
    //pm_signal4->DrawCopy("SAME");
    //pm_IC->DrawCopy("SAME");
    //pm_bhabha->DrawCopy("SAME");

    
c_pm->Print("pm_TDRcuts_normalised.pdf");



}
