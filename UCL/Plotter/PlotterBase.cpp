#define PlotterBase_cxx

#include "PlotterBase.hpp"
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <TChain.h>

using std::cout;
using std::endl;
using std::vector;
using std::to_string;


//=========================================================
//tree setup etc
//=========================================================

PlotterBase::PlotterBase(TString output_file) : output_file_name(output_file) {
  plots1D.clear();
  plots2D.clear();
}

PlotterBase::PlotterBase(TChain* chain, TString output_file) : output_file_name(output_file) {
  plots1D.clear();
  plots2D.clear();
}

PlotterBase::~PlotterBase()
{

  std::cout<<std::endl<<"Finished!"<<std::endl;
// Change directory depending on sample type
  mkdir("/unix/muons/mu3e/Alex/Run/plot_files/parallelrun/bhabha_weightcut_nvert/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
 
  TFile * output_file = new TFile("/unix/muons/mu3e/Alex/Run/plot_files/parallelrun/bhabha_weightcut_nvert/"+output_file_name+".root", "RECREATE");
  output_file->cd();
  
  
  for(int i=0; i<plots1D.size(); i++) {
    if( ((TString)plots1D[i]->GetName()).Contains("CutFlow")) PrintCutFlow(plots1D[i]);
   std::cout <<"1D plot "<<plots1D[i]->GetName();  
    if(plots1D[i]->GetEntries()==0) {
      std::cout <<" not used."<<std::endl;
      continue;
    }
    plots1D[i]->Write();
    std::cout<<" written."<<std::endl;
  }
  
  for(int i=0; i<plots2D.size(); i++) {
    if(plots2D[i]->GetEntries()==0) std::cout <<"2D Plot not used: "<<plots2D[i]->GetName()<<std::endl;
    else{
    plots2D[i]->Write();
    std::cout<<"2D plot "<<plots2D[i]->GetName()<<" written."<<std::endl;
   }
  }


 for(int i=0; i<profiles.size(); i++) {
    if(profiles[i]->GetEntries()==0) std::cout <<"Profile not used: "<<profiles[i]->GetName()<<std::endl;
    else profiles[i]->Write();
  }


  output_file->Close();
  std::cout<<"Plots written to "<<output_file->GetName()<<std::endl;
  delete output_file;


  for(int i=0; i<plots1D.size(); i++) delete plots1D[i];
  for(int i=0; i<plots2D.size(); i++) delete plots2D[i];
  for(int i=0; i<profiles.size(); i++) delete profiles[i];
  plots1D.clear();
  plots2D.clear();
  profiles.clear();
  
}


//=========================================================
//defining and filling plots:
//=========================================================

void PlotterBase::plot1D(TString name, TString title, int nbins, float xlow, float xhigh){
  plots1D.push_back(plot1DMaker(name, title, nbins, xlow, xhigh));
  cout<<name<<" "<<plots1D.size()<<endl;
}
void PlotterBase::plot1DHitCuts(TString name, TString title, int nbins, float xlow, float xhigh){
  plots1D.push_back(plot1DMaker(name, title, nbins, xlow, xhigh));
  plots1D.push_back(plot1DMaker(name+"_nh4", title, nbins, xlow, xhigh));
  plots1D.push_back(plot1DMaker(name+"_nh6", title, nbins, xlow, xhigh));
  plots1D.push_back(plot1DMaker(name+"_nh8", title, nbins, xlow, xhigh));
}
TH1 * PlotterBase::plot1DMaker(TString name, TString title, int nbins, float xlow, float xhigh){
  TH1 * p = new TH1D(name, title, nbins, xlow, xhigh);
  p->Sumw2();
  p->SetDirectory(0);
  p->GetXaxis()->SetLabelSize(0.04);
  p->GetYaxis()->SetLabelSize(0.04);
  p->GetXaxis()->SetTitleOffset(1.2);
  p->GetYaxis()->SetTitleOffset(1.2);
  p->GetXaxis()->SetLabelFont(42);
  p->GetYaxis()->SetLabelFont(42);
  p->GetXaxis()->SetTitleFont(42);
  p->GetYaxis()->SetTitleFont(42);
  //  p->SetTitle(0);
  p->SetLineWidth(2);
  p->SetMarkerSize(1.5);
  //int binmax = p->FindLastBinAbove(0,1,1,-1);
  //cout<<"Maximum bin for "+name+" is "<< binmax <<endl;
  return p;
}

void PlotterBase::plot2D(TString name, TString title, int xnbins, float xlow, float xhigh, int ynbins, float ylow, float yhigh){
  plots2D.push_back(plot2DMaker(name, title, xnbins, xlow, xhigh, ynbins, ylow, yhigh));
}
TH2 * PlotterBase::plot2DMaker(TString name, TString title, int nbinsx, float xlow, float xhigh, int nbinsy, float ylow, float yhigh){
  TH2 * p = new TH2D(name, title, nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
  p->Sumw2();
  p->SetDirectory(0);
  p->GetXaxis()->SetLabelSize(0.04);
  p->GetYaxis()->SetLabelSize(0.04);
  p->GetXaxis()->SetTitleOffset(1.2);
  p->GetYaxis()->SetTitleOffset(1.2);
  p->GetXaxis()->SetLabelFont(42);
  p->GetYaxis()->SetLabelFont(42);
  p->GetXaxis()->SetTitleFont(42);
  p->GetYaxis()->SetTitleFont(42);
  // p->SetTitle(0);
  p->SetLineWidth(2);
  p->SetMarkerSize(1.5);
  p->SetContour(99);
  return p;
}


void PlotterBase::profile(TString name, TString title, int nbins, float xlow, float xhigh){
  profiles.push_back(profileMaker(name, title, nbins, xlow, xhigh));
}
TProfile * PlotterBase::profileMaker(TString name, TString title, int nbins, float xlow, float xhigh){
  TProfile * p = new TProfile(name, title, nbins, xlow, xhigh);
  p->Sumw2();
  p->SetDirectory(0);
  p->GetXaxis()->SetLabelSize(0.04);
  p->GetYaxis()->SetLabelSize(0.04);
  p->GetXaxis()->SetTitleOffset(1.2);
  p->GetYaxis()->SetTitleOffset(1.2);
  p->GetXaxis()->SetLabelFont(42);
  p->GetYaxis()->SetLabelFont(42);
  p->GetXaxis()->SetTitleFont(42);
  p->GetYaxis()->SetTitleFont(42);
  //  p->SetTitle(0);
  p->SetLineWidth(2);
  p->SetMarkerSize(1.5);
  return p;
}


//=============================================================

void PlotterBase::Fill1D(TString name, double val, double weight) {
  for(int i=0; i<plots1D.size(); i++){
    if(name.CompareTo(plots1D[i]->GetName())) continue;
    plots1D[i]->Fill(val, weight);

    //int binmax = plots1D[-1]->FindLastBinAbove(0,1,1,-1);
    //cout<<"Maximum bin for "+name+" is "<< binmax <<endl;
    return;
  }
  
  std::cout<<"Plot1D "<<name<<" not found!"<<std::endl;
}

void PlotterBase::Fill1DHitCuts(TString name, int nhits, double val, double weight) {
  TString stub = "";
  if(nhits==4) stub="_nh4";
  else if(nhits==6) stub="_nh6";
  else if(nhits==8) stub="_nh8";

  for(int i=0; i<plots1D.size(); i++){
    if( ! name.CompareTo(plots1D[i]->GetName())) plots1D[i]->Fill(val, weight);
    else if( ! (name+stub).CompareTo(plots1D[i]->GetName())) plots1D[i]->Fill(val, weight);
  }
  // std::cout<<"Plot1D "<<name<<" not found!"<<std::endl;
}


void PlotterBase::Fill2D(TString name, double xval, double yval, double weight) {
  for(int i=0; i<plots2D.size(); i++){
    if(name.CompareTo(plots2D[i]->GetName())) continue;
    plots2D[i]->Fill(xval, yval, weight);
    return;
  }
  std::cout<<"Plot2D "<<name<<" not found!"<<std::endl;
}


void PlotterBase::FillProfile(TString name, double val, double weight) {
  for(int i=0; i<profiles.size(); i++){
    if(name.CompareTo(profiles[i]->GetName())) continue;
    profiles[i]->Fill(val, weight);
    return;
  }
  std::cout<<"Profile "<<name<<" not found!"<<std::endl;
}



//=============================================================

 
template<typename CUTST> void PlotterBase::SetCutBit(int nCut, CUTST &cutStatus, bool passCut){
  //void PlotterBase::SetCutBit(int nCut, int &cutStatus, bool failCut){
  if(passCut){
    cutStatus |= (1ULL << nCut);
  } else {
    cutStatus &= ~(1ULL << nCut);
  }
}
template void PlotterBase::SetCutBit<int>(int nCut, int &cutStatus, bool failCut);
//template void PlotterBase::SetCutBit<long>(int nCut, long &cutStatus, bool failCut);


template<class CUTST> bool PlotterBase::didCutPass(int nCut, CUTST cutStatus){
  return (cutStatus >> nCut) & 1U;
}
template bool PlotterBase::didCutFail<int>(int nCut, int cutStatus);

template<class CUTST> bool PlotterBase::didCutFail(int nCut, CUTST cutStatus){
  return !didCutPass(nCut, cutStatus);
}
template bool PlotterBase::didCutPass<int>(int nCut, int cutStatus);


//this will return -1 if 0 or >1 cuts fail
int PlotterBase::nCutsFailed(int cutStatus){
  int nfail=0;
  for(int i=0; i<16; i++) {
    if(didCutFail(i, cutStatus) ) {
      nfail++;
    }
  }
  return nfail;
}

TString PlotterBase::BitPrint(int cutStatus){
  TString cs="";
  for(int i=0; i<16;i++) cs += didCutPass(i, cutStatus);
  return cs;
}


void PlotterBase::PrintCutFlow(TH1 *plot) {

  float all = plot->GetBinContent(1);
  if(all<1) return;
  cout<<endl<<endl<<"Printing cut-flow ("<<plot->GetName()<<")"<<endl;

  float lastcut = all;
  TAxis *xax = plot->GetXaxis();

  cout<<std::setw(30) <<"Cut Name"<<"|"<<std::setw(10) <<"Number"<<"|" <<std::setw(10) <<"Abs eff"<<"|"<<"Rel eff"<<endl;

  for(int ibin=1; ibin<plot->GetNbinsX()+1; ibin++){
    float entries = plot->GetBinContent(ibin);
    if(entries==0) {
      ibin=plot->GetNbinsX()+2;
      continue;
    }
    cout<<std::setw(30) <<xax->GetBinLabel(ibin)<<"|"<<std::setw(10) <<entries<<"|"<<std::setw(10) <<entries / all<<"|"<<entries / lastcut<<endl;
    lastcut = entries;
  }
  
}


/*
void PlotterBase::printCutStatus(long long qualityStatus){
for(int i = 0; i < 16; i++){
if(cutDescriptions_[i] != ""){
if(std::find(ignoredCuts_.begin(), ignoredCuts_.end(), i) != ignoredCuts_.end()) {
std::cout << Form("%02d",i) << " : " << cutDescriptions_[i] << " : Ignored" << std::endl;
} else {
std::cout << Form("%02d",i) << " : " << cutDescriptions_[i] << " : " << (didCutPass(i,qualityStatus) ? "Pass" : "Fail") << std::endl;
}
}
}
}
*/
