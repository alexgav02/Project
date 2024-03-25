//////////////////////////////////////////////////////////
// This class implements the basic THn defining, filling and writing
// The specific tree setup is handeled in derived classes
// gavin.hesketh@ucl.ac.uk, 09/2019
//////////////////////////////////////////////////////////
#ifndef PLOTTERBASE_H_INCLUDED
#define PLOTTERBASE_H_INCLUDED

#include "TString.h"
#include <vector>
#include <TChain.h>

class TString;
class TH1;
class TH2;
class TProfile;

class PlotterBase {
 public :
  
  PlotterBase(TString output_file_name);
  PlotterBase(TChain* chain, TString output_file_name);
  ~PlotterBase();
  
  TProfile* profileMaker(TString name, TString title, int nbins, float xlow, float xhigh);  
  TH1* plot1DMaker(TString name, TString title, int nbins, float xlow, float xhigh);  
  TH2* plot2DMaker(TString name, TString title, int nbinsx, float xlow, float xhigh, int nbinsy, float ylow, float yhigh);

  void plot1D(TString name, TString title, int nbins, float xlow, float xhigh);  
  void plot1DHitCuts(TString name, TString title, int nbins, float xlow, float xhigh);  
  void Fill1D(TString name, double val, double weight=1);
  void Fill1DHitCuts(TString name, int nhits, double val, double weight=1);
  inline void Fill(TString name, double val, double weight=1) {
    Fill1D(name, val, weight);
  }

  
  void plot2D(TString name, TString title, int xnbins, float xlow, float xhigh, int ynbins, float ylow, float yhigh);
  void Fill2D(TString name, double xval, double yval, double weight=1);

  void profile(TString name, TString title, int nbins, float xlow, float xhigh);  
  void FillProfile(TString name, double val, double weight=1);


  template<typename CUTST> void SetCutBit(int nCut, CUTST &cutStatus, bool failCut);
  //  void SetCutBit(int nCut, int &cutStatus, bool failCut);
  template<class CUTST> bool didCutFail(int nCut, CUTST qualityStatus);
  template<class CUTST> bool didCutPass(int nCut, CUTST qualityStatus);
  int nCutsFailed(int cutStatus);
  void PrintCutFlow(TH1 *plot);
  TString BitPrint(int cutStatus);

  std::vector<TH1*> plots1D;
  std::vector<TProfile*> profiles;
  std::vector<TH2*> plots2D;
  
  TString output_file_name;

};

#endif
