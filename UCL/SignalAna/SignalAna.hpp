#pragma once

#include <vector>
#include <iostream>
#include "PlotterBase.hpp"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include <cmath>
#include <TSystem.h>
#include <string>
#include <TSystemDirectory.h>

using std::vector;


class SignalAna : public PlotterBase  {
public :

  SignalAna(TChain* chain, TString output_file_name, TString treename) : PlotterBase(chain,output_file_name), fChain(chain), nentries_(0), treeName_(treename) {
    std::cout<<"Reading "<<chain<<std::endl;
    if (!fChain) return;
    InitHistos();
    Init();
  };

  SignalAna(TString filename, TString output_file_name, TString treename) : PlotterBase(output_file_name), fChain(0), nentries_(0), treeName_(treename) {
    std::cout<<"Reading "<<filename<< " and tree "<<treename<<std::endl;
    TFile *f = new TFile(filename);
    f->GetObject(treename, fChain);
    if (!fChain) return;
    InitHistos();
    Init();
  };
  ~SignalAna(){
    EndJob();
    if (!fChain) return;
    delete fChain->GetCurrentFile();
  };

 // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   UInt_t          ncomb;
   UInt_t          nrightchargecomb;
   UInt_t          nvert;
   UInt_t          ntwoone;
   Double_t        weight;
   vector<double>  *vx;
   vector<double>  *vy;
   vector<double>  *vz;
   vector<double>  *chi2;
   vector<double>  *true_vx;
   vector<double>  *true_vy;
   vector<double>  *true_vz;
   vector<double>  *true_mass;
   vector<double>  *true_pxm;
   vector<double>  *true_pym;
   vector<double>  *true_pzm;
   vector<double>  *guess_vx;
   vector<double>  *guess_vy;
   vector<double>  *guess_vz;
   vector<double>  *px1;
   vector<double>  *py1;
   vector<double>  *pz1;
   vector<double>  *true_px1;
   vector<double>  *true_py1;
   vector<double>  *true_pz1;
   vector<double>  *r1;
   vector<double>  *rerr21;
   vector<int>     *id1;
   vector<int>     *pid1;
   vector<int>     *type1;
   vector<int>     *mid1;
   vector<int>     *hid1;
   vector<int>     *prime1;
   vector<double>  *px2;
   vector<double>  *py2;
   vector<double>  *pz2;
   vector<double>  *true_px2;
   vector<double>  *true_py2;
   vector<double>  *true_pz2;
   vector<double>  *r2;
   vector<double>  *rerr22;
   vector<int>     *id2;
   vector<int>     *pid2;
   vector<int>     *type2;
   vector<int>     *mid2;
   vector<int>     *hid2;
   vector<int>     *prime2;
   vector<double>  *px3;
   vector<double>  *py3;
   vector<double>  *pz3;
   vector<double>  *true_px3;
   vector<double>  *true_py3;
   vector<double>  *true_pz3;
   vector<double>  *r3;
   vector<double>  *rerr23;
   vector<int>     *id3;
   vector<int>     *pid3;
   vector<int>     *type3;
   vector<int>     *mid3;
   vector<int>     *hid3;
   vector<int>     *prime3;
   vector<double>  *pxm;
   vector<double>  *pym;
   vector<double>  *pzm;
   vector<double>  *mmass;
   vector<double>  *pcoplanar;
   vector<double>  *pacoplanar;
   vector<double>  *pm;
   vector<double>  *targetdist;
   vector<double>  *mep1;
   vector<double>  *mep2;
   vector<int>     *nhit1;
   vector<int>     *nhit2;
   vector<int>     *nhit3;
   vector<int>     *nsharedhit1;
   vector<int>     *nsharedhit2;
   vector<int>     *nsharedhit3;
   vector<int>     *nsharedseg1;
   vector<int>     *nsharedseg2;
   vector<int>     *nsharedseg3;
   vector<double>  *guessx;
   vector<double>  *guessy;
   vector<double>  *guessz;
   vector<double>  *guessdist;
   vector<double>  *true_targetdist;
   vector<int>     *mc_signal;
   vector<double>  *fit_tolerance;
   vector<int>     *fit_steps;
   vector<double>  *cosalpha;
   vector<double>  *theta_kuno;
   vector<double>  *phi_kuno;


  void     InitHistos();
  void     Init();
  void     Loop();
  void     FillPlots();
  void     ReadEvent(int jentry);
  void     EndJob() {};

  template<class VAR> void Enable(TString name, VAR &var);
  
  TTree          *fChain;
  int nentries_;
  TString treeName_;

  //cut flags:
  int SetRecoCuts(int i);
  void RecoNMinusOnePlots(int ivertex, int cutStatus);
  void CutCorrelation(int cutStatus, TString plotname);
  void RecoCutFlow(int cutStatus, TString type);
  int SetFrameCuts(int cutStatus);

  int SetTruthCuts(int i);
  void TruthCutFlow(int cutStatus);
  void AKCutFlow(int cutStatus, int cutStatusTrue);

  float true_vr;

  //can have 16 cuts in an int cutstring!
  enum RCUTS {CHI2, MASSLOOSE, MASS, MOM,
    NHIT1, NHIT2, NHIT3,
    NSHAHIT1, NSHAHIT2, NSHAHIT3,
    NSHASEG1, NSHASEG2, NSHASEG3,
    MEP1, P1, P2};
  
  enum TCUTS { TVTX, TVR, TVZ,
    THID1, THID2, THID3, TMID,
    TPRIME1, TPRIME2, TPRIME3,
    TPID1, TPID2, TPID3, TTID, TTYPE};

  
};

template<class VAR>
void SignalAna::Enable(TString name, VAR &var) { 
    std::cout<<"Activating branch "<<name<<std::endl;
    var=0;
    fChain->SetBranchStatus(name,1);
    fChain->SetBranchAddress(name, &var);
}


void SignalAna::Init() {
   fChain->SetMakeClass(1);
   fChain->SetBranchStatus("*",0);

   //enable branches you want to use:
  
   Enable("chi2", chi2);
   Enable("pm", pm);
   Enable("mmass", mmass);
   Enable("nhit1", nhit1);
   Enable("nhit2", nhit2);
   Enable("nhit3", nhit3);
   Enable("mep1", mep1);
   Enable("mep2", mep2);
   Enable("px1", px1);
   Enable("py1", py1);
   Enable("pz1", pz1);
   Enable("px2", px2);
   Enable("py2", py2);
   Enable("pz2", pz2);
   Enable("px3", px3);
   Enable("py3", py3);
   Enable("pz3", pz3);
   Enable("vx",vx);
   Enable("vy",vy);
   Enable("vz",vz);
   Enable("targetdist",targetdist);
   Enable("pxm",pxm);
   Enable("pym",pym);
   Enable("pzm",pzm);
   //Enable("p1",p1);
   //Enable("p2",p2);
   Enable("pcoplanar",pcoplanar);
   Enable("pacoplanar",pacoplanar);
   Enable("theta_kuno", theta_kuno);
   Enable("phi_kuno", phi_kuno);
   Enable("cosalpha", cosalpha);
   Enable("weight", weight);
   

   
   Enable("nsharedhit1", nsharedhit1);
   Enable("nsharedhit2", nsharedhit2);
   Enable("nsharedhit3", nsharedhit3);
   Enable("nsharedseg1", nsharedseg1);
   Enable("nsharedseg2", nsharedseg2);
   Enable("nsharedseg3", nsharedseg3);
   
   //truth cuts:
   Enable("true_vx", true_vx);
   Enable("true_vy", true_vy);
   Enable("true_vz", true_vz);
   Enable("true_pxm",true_pxm);
   Enable("true_pym",true_pym);
   Enable("true_pzm",true_pzm);
   Enable("true_targetdist", true_targetdist);
   Enable("true_mass", true_mass);
   Enable("true_px1", true_px1);
   Enable("true_py1", true_py1);
   Enable("true_pz1", true_pz1);
   Enable("true_px2", true_px2);
   Enable("true_py2", true_py2);
   Enable("true_pz2", true_pz2);
   Enable("true_px3", true_px3);
   Enable("true_py3", true_py3);
   Enable("true_pz3", true_pz3);
   Enable("pid1", pid1);
   Enable("pid2", pid2);
   Enable("pid3", pid3);

   Enable("hid1", hid1);
   Enable("hid2", hid2);
   Enable("hid3", hid3);

   Enable("mid1", mid1);
   Enable("mid2", mid2);
   Enable("mid3", mid3);

   Enable("id1", id1);
   Enable("id2", id2);
   Enable("id3", id3);

   Enable("prime1", prime1);
   Enable("prime2", prime2);
   Enable("prime3", prime3);

   Enable("type1", type1);
   Enable("type2", type2);
   Enable("type3", type3);

   nentries_ = fChain->GetEntriesFast();

   std::cout<<"Tree contains "<<nentries_<<" entries."<<std::endl;

}

//============================================


void SignalAna::ReadEvent(int jentry) {
  const int frac = 100*(jentry+1) / nentries_;
  if(frac%5 == 0) {
    //std::cout<<"\r"<<"Processed "<<frac<<"% of "<<nentries_<<" "<<std::flush;
  }
  fChain->LoadTree(jentry);
  fChain->GetEntry(jentry);

  true_vr=-1;

  
}




//============================================================================

int SignalAna::SetRecoCuts(int ivertex) {

  int cutStatus=0;

  //apply cuts: set bit to 1 if passed
  SetCutBit(CHI2, cutStatus,  chi2->at(ivertex) < 15 );
  SetCutBit(MASSLOOSE, cutStatus,  mmass->at(ivertex) > 80  && mmass->at(ivertex) < 115 );
  SetCutBit(MASS, cutStatus,  mmass->at(ivertex) > 103  && mmass->at(ivertex) < 110 );
  SetCutBit(MOM, cutStatus,  pm->at(ivertex) <4 );
  
  SetCutBit(NHIT1, cutStatus,  (nhit1->at(ivertex) > 4) );
  SetCutBit(NHIT2, cutStatus,  (nhit2->at(ivertex) > 4) );
  SetCutBit(NHIT3, cutStatus,  (nhit3->at(ivertex) > 4) );

  SetCutBit(NSHAHIT1, cutStatus, nsharedhit1->at(ivertex) ==0 );
  SetCutBit(NSHAHIT2, cutStatus, nsharedhit2->at(ivertex) ==0 );
  SetCutBit(NSHAHIT3, cutStatus, nsharedhit3->at(ivertex) ==0 );
  SetCutBit(NSHASEG1, cutStatus, nsharedseg1->at(ivertex) ==0 );
  SetCutBit(NSHASEG2, cutStatus, nsharedseg2->at(ivertex) ==0 );
  SetCutBit(NSHASEG3, cutStatus, nsharedseg3->at(ivertex) ==0 );

  SetCutBit(MEP1, cutStatus, mep1->at(ivertex)<5 || mep1->at(ivertex)>10 );

  //if( p1<0 ) p1 = sqrt(px1->at(ivertex)*px1->at(ivertex) + py1->at(ivertex)*py1->at(ivertex) + pz1->at(ivertex)*pz1->at(ivertex));
  //if( p2<0 ) p2 = sqrt(px2->at(ivertex)*px2->at(ivertex) + py2->at(ivertex)*py2->at(ivertex) + pz2->at(ivertex)*pz2->at(ivertex));

  //SetCutBit(P1, cutStatus, p1->at(ivertex)<50 );
  //SetCutBit(P2, cutStatus, p2->at(ivertex)<50 );
  SetCutBit(P1, cutStatus, sqrt(px1->at(ivertex)*px1->at(ivertex) + py1->at(ivertex)*py1->at(ivertex) + pz1->at(ivertex)*pz1->at(ivertex))<50 );
  SetCutBit(P2, cutStatus, sqrt(px2->at(ivertex)*px2->at(ivertex) + py2->at(ivertex)*py2->at(ivertex) + pz2->at(ivertex)*pz2->at(ivertex))<50 );

  //a cutStatus of zero means all cuts passed!  
  return cutStatus;
  
}


//============================================

int SignalAna::SetFrameCuts(int cutStatus) {

       int tmp_cutStatus =0;

       bool rec1 = didCutPass(NHIT1, cutStatus);
       bool rec2 = didCutPass(NHIT2, cutStatus);
       bool rec3 = didCutPass(NHIT3, cutStatus);
      
       SetCutBit(NHIT1, tmp_cutStatus, rec1 );
       SetCutBit(NHIT2, tmp_cutStatus, rec2 );
       SetCutBit(NHIT3, tmp_cutStatus, rec3 );

       if(!rec1 || !rec2 || !rec3) return tmp_cutStatus;

       if( !didCutPass(CHI2, cutStatus) ) return tmp_cutStatus;
       SetCutBit(CHI2, tmp_cutStatus,  true );

       if( !didCutPass(MOM, cutStatus) ) return tmp_cutStatus;
       SetCutBit(MOM, tmp_cutStatus,  true );

       if( !didCutPass(MEP1, cutStatus) ) return tmp_cutStatus;
       SetCutBit(MEP1, tmp_cutStatus,  true );

       if( !didCutPass(MASS, cutStatus) ) return tmp_cutStatus;
       SetCutBit(MASS, tmp_cutStatus,  true );

       return tmp_cutStatus;
}

//============================================

void SignalAna::RecoCutFlow(int cutStatus, TString type) {
  //the ordering used in the TDR:
  // long tracks
  // recurler rejection / fit convergence ??
  // vertex chi^2 < 15
  // p < 4
  // mep1 <5 || >10
  // vertex mass 103 - 110

  int icut=0;
  
  Fill("RecoCutFlow_"+type,  icut++);

  int nrecurl= didCutPass(NHIT1, cutStatus);
  nrecurl+= didCutPass(NHIT2, cutStatus);
  nrecurl+= didCutPass(NHIT3, cutStatus);

  if(nrecurl>0)   Fill("RecoCutFlow_"+type,  icut++);
  else return;
  if(nrecurl>1)   Fill("RecoCutFlow_"+type,  icut++);
  else return;
  if(nrecurl>2)   Fill("RecoCutFlow_"+type,  icut++);
  else return;
  if( didCutPass(CHI2, cutStatus) )   Fill("RecoCutFlow_"+type,  icut++);
  else return;
  if( didCutPass(MOM, cutStatus) )   Fill("RecoCutFlow_"+type,  icut++);
  else return;
  if( didCutPass(MEP1, cutStatus) )   Fill("RecoCutFlow_"+type,  icut++);
  else return;
  if( didCutPass(MASS, cutStatus) )   Fill("RecoCutFlow_"+type,  icut++);

  return;
}


//======================================================================


int SignalAna::SetTruthCuts(int ivertex) {

  int cutStatus=0;  
  //apply cuts:

  //  BaseCut & truevxcut_r = true_vr < 19 && abs(true_vz) < 50 && pid1 == -11 && pid2 == -11 && pid3 == 11 && prime1 > 0 && prime2 > 0 && prime3 > 0 && threedifferent;


  const float vx = true_vx->at(ivertex);
  const float vy = true_vy->at(ivertex);
  const float vz = true_vz->at(ivertex);

  SetCutBit(TVTX, cutStatus, (vx>-888 && vx< 9000) &&  (vy>-888 && vy< 9000) &&  (vz>-888 && vz< 9000) );

  if(true_vr<0) true_vr = sqrt(vx*vx + vy*vy);
  
  SetCutBit(TVR, cutStatus,  true_vr < 19 );
  SetCutBit(TVZ, cutStatus,  fabs( vz ) < 50 );
  
  //is this not the same as mcprime?
  SetCutBit(THID1, cutStatus,  hid1->at(ivertex) < 3 && hid1->at(ivertex) > -3 );
  SetCutBit(THID2, cutStatus,  hid2->at(ivertex) < 3 && hid2->at(ivertex) > -3 );
  SetCutBit(THID3, cutStatus,  hid3->at(ivertex) < 3 && hid3->at(ivertex) > -3 );
  
  SetCutBit(TMID, cutStatus,  (mid1->at(ivertex) == mid2->at(ivertex) )  &&
	    (mid1->at(ivertex) == mid3->at(ivertex) ) );
  
  SetCutBit(TPRIME1, cutStatus,  prime1->at(ivertex) == 1 );
  SetCutBit(TPRIME2, cutStatus,  prime2->at(ivertex) == 1 );
  SetCutBit(TPRIME3, cutStatus,  prime3->at(ivertex) == 1 );
  
  SetCutBit(TPID1, cutStatus,  pid1->at(ivertex) ==-11 );
  SetCutBit(TPID2, cutStatus,  pid2->at(ivertex) ==-11 );
  SetCutBit(TPID3, cutStatus,  pid3->at(ivertex) ==11 );

  const int cid1 = id1->at(ivertex);
  const int cid2 = id2->at(ivertex);
  const int cid3 = id3->at(ivertex);
 
  SetCutBit(TTID, cutStatus, (cid1!= cid2) && (cid1 != cid3) && (cid2 != cid3) );

  SetCutBit(TTYPE, cutStatus, type1->at(ivertex)==91 || type1->at(ivertex)==92 ||
	    type2->at(ivertex)==91 || type2->at(ivertex)==92 ||
	    type3->at(ivertex)==91 || type3->at(ivertex)==92 );

  
  return cutStatus;
  }
//============================================

void SignalAna::TruthCutFlow(int cutStatus) {

  int icut=0;

  Fill("TruthCutFlow_vertices", icut++);

  if(didCutPass(TVTX, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;

  if(didCutPass(TVR, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;

  if(didCutPass(TVZ, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;

  if(didCutPass(TPRIME1, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;
  if(didCutPass(TPRIME2, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;
  if(didCutPass(TPRIME3, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;

  if(didCutPass(THID1, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;
  if(didCutPass(THID2, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;
  if(didCutPass(THID3, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;
  
  if(didCutPass(TMID, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;
  
  if(didCutPass(TPID1, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;
  if(didCutPass(TPID2, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;
  if(didCutPass(TPID3, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;

  if(didCutPass(TTID, cutStatus) )  Fill("TruthCutFlow_vertices", icut++);
  else return;
  
}
//============================================

void SignalAna::AKCutFlow(int cutStatus, int cutStatusTrue) {

  int icut=0;

  Fill("AKCutFlow", icut++);

  if(didCutPass(TTYPE, cutStatusTrue) )  Fill("AKCutFlow", icut++);
  else return;

  if(didCutPass(TVR, cutStatusTrue) && didCutPass(TVZ, cutStatusTrue) )  Fill("AKCutFlow", icut++);
  else return;

  if(didCutPass(TTID, cutStatusTrue) )  Fill("AKCutFlow", icut++);
  else return;

  // this gives nGenerated
  
  if(didCutPass(TPID1, cutStatusTrue) && didCutPass(TPID2, cutStatusTrue) && didCutPass(TPID3, cutStatusTrue) )  Fill("AKCutFlow", icut++);
  else return;
  
  if(didCutPass(TMID, cutStatusTrue) )  Fill("AKCutFlow", icut++);
  else return;

  if(didCutPass(TPRIME1, cutStatusTrue) && didCutPass(TPRIME2, cutStatusTrue)  && didCutPass(TPRIME3, cutStatusTrue) )  Fill("AKCutFlow", icut++);
  else return;

  if(didCutPass(THID1, cutStatusTrue) && didCutPass(THID2, cutStatusTrue) && didCutPass(THID3, cutStatusTrue) )  Fill("AKCutFlow", icut++);
  else return;

  int nrecurl= didCutPass(NHIT1, cutStatus);
  nrecurl+= didCutPass(NHIT2, cutStatus);
  nrecurl+= didCutPass(NHIT3, cutStatus);
  if(nrecurl>2)   Fill("AKCutFlow", icut++);
  else return;

  if( didCutPass(CHI2, cutStatus) )   Fill("AKCutFlow",  icut++);
  else return;
  if( didCutPass(MOM, cutStatus) )   Fill("AKCutFlow",  icut++);
  else return;
  if( didCutPass(MASS, cutStatus) )   Fill("AKCutFlow",  icut++);
  else return;
  if( didCutPass(MEP1, cutStatus) )   Fill("AKCutFlow",  icut++);
  else return;


  
}

//======================================================================

//fill a 2D plot showing which cuts fail:
void SignalAna::CutCorrelation(int cutStatus, TString plotname) {
  int ncuts = 0;
  for (int i=0; i<plots2D.size() && ncuts<1; i++) {
    if(plotname.CompareTo(plots2D[i]->GetName())) continue;
    ncuts = plots2D[i]->GetNbinsX();
  }
  for(int i=0; i<ncuts; i++) {
    if(didCutPass(i, cutStatus)) continue;
    for(int j=i; j<ncuts; j++) {
      if(didCutPass(j, cutStatus)) continue;
      Fill2D(plotname, i,j);
    }
  }
}



//======================================================================


void SignalAna::RecoNMinusOnePlots(int ivertex, int cutStatus) {
  
  int failCut = nCutsFailed(cutStatus); 
  //reject if failing more than one:
  if(failCut>1) return; 
  
  
  if(failCut==0 || didCutFail(0, cutStatus) )  Fill("nMinusOne_chi2", chi2->at(ivertex));
  if(failCut==0 || didCutFail(2, cutStatus) )  Fill("nMinusOne_mass", mmass->at(ivertex));
  if(failCut==0 || didCutFail(3, cutStatus) )  Fill("nMinusOne_p", pm->at(ivertex));
  if(failCut==0 || didCutFail(4, cutStatus) )  Fill("nMinusOne_nHit1", nhit1->at(ivertex));
  if(failCut==0 || didCutFail(5, cutStatus) )  Fill("nMinusOne_nHit2", nhit2->at(ivertex));
  if(failCut==0 || didCutFail(6, cutStatus) )  Fill("nMinusOne_nHit3", nhit3->at(ivertex));
  if(failCut==0 || didCutFail(7, cutStatus) )  Fill("nMinusOne_nSharedHit1", nsharedhit1->at(ivertex));
  if(failCut==0 || didCutFail(8, cutStatus) )  Fill("nMinusOne_nSharedHit2", nsharedhit2->at(ivertex));
  if(failCut==0 || didCutFail(9, cutStatus) )  Fill("nMinusOne_nSharedHit3", nsharedhit3->at(ivertex));
  if(failCut==0 || didCutFail(10, cutStatus) )  Fill("nMinusOne_nSharedSeg1", nsharedseg1->at(ivertex));
  if(failCut==0 || didCutFail(11, cutStatus) )  Fill("nMinusOne_nSharedSeg2", nsharedseg2->at(ivertex));
  if(failCut==0 || didCutFail(12, cutStatus) )  Fill("nMinusOne_nSharedSeg3", nsharedseg3->at(ivertex));
  if(failCut==0 || didCutFail(13, cutStatus) )  Fill("nMinusOne_mep", mep1->at(ivertex));
  if(failCut==0 || didCutFail(14, cutStatus) )  Fill("nMinusOne_p1", mep1->at(ivertex));
  if(failCut==0 || didCutFail(15, cutStatus) )  Fill("nMinusOne_p2", mep1->at(ivertex));

}


void addSkimFilesToChain(const std::string& folderPath, TChain& chain, int startNum, int numFiles) {
    std::vector<std::string> excludedFiles = {
"skim_100000.root",
"skim_100007.root",
"skim_100017.root",
"skim_100023.root",
"skim_100028.root",
"skim_100031.root",
"skim_100035.root",
"skim_100043.root",
"skim_100052.root",
"skim_100060.root",
"skim_100069.root",
"skim_100085.root",
"skim_100100.root",
"skim_100101.root",
"skim_100104.root",
"skim_100106.root",
"skim_100118.root",
"skim_100122.root",
"skim_100135.root",
"skim_100150.root",
"skim_100172.root",
"skim_100195.root",
"skim_100205.root",
"skim_100206.root",
"skim_100212.root",
"skim_100223.root",
"skim_100234.root",
"skim_100270.root",
"skim_100271.root",
"skim_100274.root",
"skim_100275.root",
"skim_100282.root",
"skim_100288.root",
"skim_100289.root",
"skim_100297.root",
"skim_100305.root",
"skim_100306.root",
"skim_100333.root",
"skim_100335.root",
"skim_100340.root",
"skim_100381.root",
"skim_100408.root",
"skim_100423.root",
"skim_100431.root",
"skim_100446.root",
"skim_100464.root",
"skim_100468.root",
"skim_100476.root",
"skim_100490.root",
"skim_100491.root",
"skim_100495.root",
"skim_100498.root",
"skim_100504.root",
"skim_100505.root",
"skim_100518.root",
"skim_100519.root",
"skim_100524.root",
"skim_100535.root",
"skim_100551.root",
"skim_100575.root",
"skim_100579.root",
"skim_100580.root",
"skim_100594.root",
"skim_100605.root",
"skim_100609.root",
"skim_100644.root",
"skim_100656.root",
"skim_100673.root",
"skim_100675.root",
"skim_100681.root",
"skim_100701.root",
"skim_100703.root",
"skim_100706.root",
"skim_100710.root",
"skim_100715.root",
"skim_100728.root",
"skim_100755.root",
"skim_100779.root",
"skim_100791.root",
"skim_100793.root",
"skim_100798.root",
"skim_100820.root",
"skim_100831.root",
"skim_100846.root",
"skim_100858.root",
"skim_100867.root",
"skim_100869.root",
"skim_100877.root",
"skim_100878.root",
"skim_100887.root",
"skim_100896.root",
"skim_100899.root",
"skim_100911.root",
"skim_100914.root",
"skim_100959.root",
"skim_100963.root",
"skim_101012.root",
"skim_101017.root",
"skim_101033.root",
"skim_101047.root",
"skim_101074.root",
"skim_101076.root",
"skim_101080.root",
"skim_101088.root",
"skim_101094.root",
"skim_101099.root",
"skim_101121.root",
"skim_101128.root",
"skim_101138.root",
"skim_101167.root",
"skim_101179.root",
"skim_101183.root",
"skim_101184.root",
"skim_101187.root",
"skim_101195.root",
"skim_101198.root",
"skim_101199.root",
"skim_101201.root",
"skim_101216.root",
"skim_101248.root",
"skim_101258.root",
"skim_101261.root",
"skim_101289.root",
"skim_101303.root",
"skim_101305.root",
"skim_101333.root",
"skim_101345.root",
"skim_101355.root",
"skim_101363.root",
"skim_101374.root",
"skim_101388.root",
"skim_101404.root",
"skim_101407.root",
"skim_101419.root",
"skim_101438.root",
"skim_101457.root",
"skim_101491.root",
"skim_101496.root",
"skim_101504.root",
"skim_101515.root",
"skim_101519.root",
"skim_101533.root",
"skim_101535.root",
"skim_101554.root",
"skim_101560.root",
"skim_101568.root",
"skim_101581.root",
"skim_101584.root",
"skim_101585.root",
"skim_101592.root",
"skim_101606.root",
"skim_101623.root",
"skim_101624.root",
"skim_101632.root",
"skim_101646.root",
"skim_101650.root",
"skim_101666.root",
"skim_101683.root",
"skim_101685.root",
"skim_101692.root",
"skim_101695.root",
"skim_101697.root",
"skim_101708.root",
"skim_101719.root",
"skim_101732.root",
"skim_101749.root",
"skim_101750.root",
"skim_101772.root",
"skim_101792.root",
"skim_101794.root",
"skim_101798.root",
"skim_101814.root",
"skim_101816.root",
"skim_101829.root",
"skim_101834.root",
"skim_101836.root",
"skim_101838.root",
"skim_101848.root",
"skim_101852.root",
"skim_101864.root",
"skim_101865.root",
"skim_101872.root",
"skim_101873.root",
"skim_101874.root",
"skim_101880.root",
"skim_101888.root",
"skim_101894.root",
"skim_101921.root",
"skim_101925.root",
"skim_101930.root",
"skim_101937.root",
"skim_101976.root",
"skim_101977.root",
"skim_101983.root",
"skim_102008.root",
"skim_102012.root",
"skim_102026.root",
"skim_102028.root",
"skim_102033.root",
"skim_102049.root",
"skim_102057.root",
"skim_102065.root",
"skim_102077.root",
"skim_102095.root",
"skim_102098.root",
"skim_102105.root",
"skim_102109.root",
"skim_102115.root",
"skim_102120.root",
"skim_102191.root",
"skim_102199.root",
"skim_102223.root",
"skim_102245.root",
"skim_102280.root",
"skim_102296.root",
"skim_102310.root",
"skim_102314.root",
"skim_102326.root",
"skim_102338.root",
"skim_102340.root",
"skim_102342.root",
"skim_102346.root",
"skim_102350.root",
"skim_102359.root",
"skim_102371.root",
"skim_102376.root",
"skim_102384.root",
"skim_102386.root",
"skim_102399.root",
"skim_102409.root",
"skim_102414.root",
"skim_102431.root",
"skim_102432.root",
"skim_102457.root",
"skim_102459.root",
"skim_102471.root",
"skim_102476.root",
"skim_102519.root",
"skim_102520.root",
"skim_102529.root",
"skim_102538.root",
"skim_102539.root",
"skim_102545.root",
"skim_102577.root",
"skim_102586.root",
"skim_102601.root",
"skim_102618.root",
"skim_102629.root",
"skim_102646.root",
"skim_102650.root",
"skim_102652.root",
"skim_102653.root",
"skim_102657.root",
"skim_102658.root",
"skim_102664.root",
"skim_102667.root",
"skim_102694.root",
"skim_102696.root",
"skim_102701.root",
"skim_102737.root",
"skim_102785.root",
"skim_102791.root",
"skim_102801.root",
"skim_102802.root",
"skim_102820.root",
"skim_102823.root",
"skim_102846.root",
"skim_102860.root",
"skim_102868.root",
"skim_102870.root",
"skim_102872.root",
"skim_102875.root",
"skim_102880.root",
"skim_102881.root",
"skim_102903.root",
"skim_102910.root",
"skim_102915.root",
"skim_102924.root",
"skim_102937.root",
"skim_102938.root",
"skim_102940.root",
"skim_102957.root",
"skim_103011.root",
"skim_103142.root",
"skim_103144.root",
"skim_103155.root",
"skim_103171.root",
"skim_103191.root",
"skim_103192.root",
"skim_103205.root",
"skim_103215.root",
"skim_103245.root",
"skim_103248.root",
"skim_103250.root",
"skim_103252.root",
"skim_103266.root",
"skim_103269.root",
"skim_103271.root",
"skim_103276.root",
"skim_103283.root",
"skim_103296.root",
"skim_103357.root",
"skim_103359.root",
"skim_103407.root",
"skim_103429.root",
"skim_103446.root",
"skim_103453.root",
"skim_103463.root",
"skim_103478.root",
"skim_103484.root",
"skim_103508.root",
"skim_103519.root",
"skim_103524.root",
"skim_103526.root",
"skim_103529.root",
"skim_103540.root",
"skim_103546.root",
"skim_103558.root",
"skim_103563.root",
"skim_103568.root",
"skim_103569.root",
"skim_103573.root",
"skim_103580.root",
"skim_103587.root",
"skim_103589.root",
"skim_103598.root",
"skim_103630.root",
"skim_103640.root",
"skim_103675.root",
"skim_103679.root",
"skim_103688.root",
"skim_103691.root",
"skim_103692.root",
"skim_103695.root",
"skim_103698.root",
"skim_103713.root",
"skim_103727.root",
"skim_103728.root",
"skim_103730.root",
"skim_103740.root",
"skim_103765.root",
"skim_103777.root",
"skim_103791.root",
"skim_103807.root",
"skim_103821.root",
"skim_103825.root",
"skim_103835.root",
"skim_103839.root",
"skim_103840.root",
"skim_103864.root",
"skim_103880.root",
"skim_103893.root",
"skim_103899.root",
"skim_103905.root",
"skim_103906.root",
"skim_103911.root",
"skim_103918.root",
"skim_103920.root",
"skim_103927.root",
"skim_103936.root",
"skim_103937.root",
"skim_103959.root",
"skim_103966.root",
"skim_103972.root",
"skim_103977.root",
"skim_103988.root",
"skim_104009.root",
"skim_104010.root",
"skim_104041.root",
"skim_104044.root",
"skim_104061.root",
"skim_104081.root",
"skim_104100.root",
"skim_104124.root",
"skim_104143.root",
"skim_104144.root",
"skim_104150.root",
"skim_104178.root",
"skim_104180.root",
"skim_104186.root",
"skim_104196.root",
"skim_104215.root",
"skim_104232.root",
"skim_104247.root",
"skim_104278.root",
"skim_104292.root",
"skim_104295.root",
"skim_104302.root",
"skim_104318.root",
"skim_104325.root",
"skim_104332.root",
"skim_104369.root",
"skim_104386.root",
"skim_104403.root",
"skim_104406.root",
"skim_104408.root",
"skim_104446.root",
"skim_104453.root",
"skim_104458.root",
"skim_104508.root",
"skim_104524.root",
"skim_104542.root",
"skim_104589.root",
"skim_104607.root",
"skim_104623.root",
"skim_104630.root",
"skim_104639.root",
"skim_104646.root",
"skim_104651.root",
"skim_104654.root",
"skim_104664.root",
"skim_104678.root",
"skim_104683.root",
"skim_104696.root",
"skim_104700.root",
"skim_104701.root",
"skim_104702.root",
"skim_104742.root",
"skim_104752.root",
"skim_104753.root",
"skim_104760.root",
"skim_104765.root",
"skim_104773.root",
"skim_104784.root",
"skim_104785.root",
"skim_104799.root",
"skim_104801.root",
"skim_104802.root",
"skim_104812.root",
"skim_104824.root",
"skim_104829.root",
"skim_104836.root",
"skim_104844.root",
"skim_104863.root",
"skim_104872.root",
"skim_104876.root",
"skim_104891.root",
"skim_104903.root",
"skim_104912.root",
"skim_104968.root",
"skim_104971.root",
"skim_104976.root",
"skim_104989.root",
"skim_105004.root",
"skim_105007.root",
"skim_105009.root",
"skim_105011.root",
"skim_105023.root",
"skim_105042.root",
"skim_105061.root",
"skim_105068.root",
"skim_105087.root",
"skim_105127.root",
"skim_105135.root",
"skim_105147.root",
"skim_105187.root",
"skim_105197.root",
"skim_105217.root",
"skim_105227.root",
"skim_105260.root",
"skim_105269.root",
"skim_105275.root",
"skim_105276.root",
"skim_105287.root",
"skim_105290.root",
"skim_105294.root",
"skim_105296.root",
"skim_105308.root",
"skim_105318.root",
"skim_105321.root",
"skim_105334.root",
"skim_105345.root",
"skim_105348.root",
"skim_105352.root",
"skim_105358.root",
"skim_105370.root",
"skim_105403.root",
"skim_105408.root",
"skim_105424.root",
"skim_105438.root",
"skim_105445.root",
"skim_105461.root",
"skim_105464.root",
"skim_105493.root",
"skim_105496.root",
"skim_105502.root",
"skim_105504.root",
"skim_105508.root",
"skim_105520.root",
"skim_105521.root",
"skim_105525.root",
"skim_105544.root",
"skim_105547.root",
"skim_105554.root",
"skim_105561.root",
"skim_105592.root",
"skim_105605.root",
"skim_105606.root",
"skim_105624.root",
"skim_105627.root",
"skim_105642.root",
"skim_105649.root",
"skim_105656.root",
"skim_105660.root",
"skim_105686.root",
"skim_105690.root",
"skim_105691.root",
"skim_105714.root",
"skim_105719.root",
"skim_105720.root",
"skim_105726.root",
"skim_105735.root",
"skim_105736.root",
"skim_105742.root",
"skim_105743.root",
"skim_105762.root",
"skim_105769.root",
"skim_105770.root",
"skim_105779.root",
"skim_105782.root",
"skim_105795.root",
"skim_105830.root",
"skim_105841.root",
"skim_105847.root",
"skim_105854.root",
"skim_105860.root",
"skim_105867.root",
"skim_105887.root",
"skim_105902.root",
"skim_105903.root",
"skim_105918.root",
"skim_105920.root",
"skim_105926.root",
"skim_105934.root",
"skim_105935.root",
"skim_105939.root",
"skim_105955.root",
"skim_105966.root",
"skim_105968.root",
"skim_105991.root",
"skim_106012.root",
"skim_106016.root",
"skim_106036.root",
"skim_106042.root",
"skim_106050.root",
"skim_106060.root",
"skim_106062.root",
"skim_106064.root",
"skim_106075.root",
"skim_106077.root",
"skim_106090.root",
"skim_106118.root",
"skim_106127.root",
"skim_106131.root",
"skim_106132.root",
"skim_106175.root",
"skim_106191.root",
"skim_106192.root",
"skim_106194.root",
"skim_106205.root",
"skim_106207.root",
"skim_106209.root",
"skim_106227.root",
"skim_106230.root",
"skim_106231.root",
"skim_106233.root",
"skim_106240.root",
"skim_106253.root",
"skim_106265.root",
"skim_106275.root",
"skim_106284.root",
"skim_106287.root",
"skim_106321.root",
"skim_106338.root",
"skim_106346.root",
"skim_106364.root",
"skim_106370.root",
"skim_106374.root",
"skim_106412.root",
"skim_106429.root",
"skim_106438.root",
"skim_106439.root",
"skim_106444.root",
"skim_106473.root",
"skim_106497.root",
"skim_106509.root",
"skim_106519.root",
"skim_106544.root",
"skim_106550.root",
"skim_106552.root",
"skim_106557.root",
"skim_106558.root",
"skim_106561.root",
"skim_106563.root",
"skim_106587.root",
"skim_106593.root",
"skim_106624.root",
"skim_106630.root",
"skim_106634.root",
"skim_106652.root",
"skim_106659.root",
"skim_106677.root",
"skim_106680.root",
"skim_106681.root",
"skim_106682.root",
"skim_106686.root",
"skim_106693.root",
"skim_106710.root",
"skim_106718.root",
"skim_106730.root",
"skim_106745.root",
"skim_106756.root",
"skim_106758.root",
"skim_106780.root",
"skim_106781.root",
"skim_106786.root",
"skim_106808.root",
"skim_106818.root",
"skim_106823.root",
"skim_106843.root",
"skim_106856.root",
"skim_106864.root",
"skim_106867.root",
"skim_106872.root",
"skim_106880.root",
"skim_106885.root",
"skim_106895.root",
"skim_106900.root",
"skim_106924.root",
"skim_106935.root",
"skim_106943.root",
"skim_106950.root",
"skim_106978.root",
"skim_106980.root",
"skim_107010.root",
"skim_107015.root",
"skim_107022.root",
"skim_107028.root",
"skim_107029.root",
"skim_107034.root",
"skim_107043.root",
"skim_107055.root",
"skim_107059.root",
"skim_107065.root",
"skim_107071.root",
"skim_107073.root",
"skim_107082.root",
"skim_107084.root",
"skim_107086.root",
"skim_107092.root",
"skim_107097.root",
"skim_107123.root",
"skim_107158.root",
"skim_107163.root",
"skim_107167.root",
"skim_107187.root",
"skim_107189.root",
"skim_107202.root",
"skim_107213.root",
"skim_107218.root",
"skim_107222.root",
"skim_107251.root",
"skim_107253.root",
"skim_107257.root",
"skim_107265.root",
"skim_107274.root",
"skim_107278.root",
"skim_107284.root",
"skim_107286.root",
"skim_107290.root",
"skim_107322.root",
"skim_107325.root",
"skim_107330.root",
"skim_107345.root",
"skim_107346.root",
"skim_107353.root",
"skim_107360.root",
"skim_107365.root",
"skim_107390.root",
"skim_107405.root",
"skim_107420.root",
"skim_107439.root",
"skim_107451.root",
"skim_107452.root",
"skim_107476.root",
"skim_107480.root",
"skim_107509.root",
"skim_107528.root",
"skim_107553.root",
"skim_107581.root",
"skim_107591.root",
"skim_107614.root",
"skim_107628.root",
"skim_107641.root",
"skim_107648.root",
"skim_107670.root",
"skim_107677.root",
"skim_107714.root",
"skim_107738.root",
"skim_107747.root",
"skim_107765.root",
"skim_107769.root",
"skim_107770.root",
"skim_107782.root",
"skim_107804.root",
"skim_107833.root",
"skim_107836.root",
"skim_107844.root",
"skim_107846.root",
"skim_107847.root",
"skim_107859.root",
"skim_107863.root",
"skim_107879.root",
"skim_107880.root",
"skim_107890.root",
"skim_107894.root",
"skim_107911.root",
"skim_107915.root",
"skim_107923.root",
"skim_107925.root",
"skim_107926.root",
"skim_107930.root",
"skim_107988.root",
"skim_107991.root",
"skim_107998.root",
"skim_108014.root",
"skim_108023.root",
"skim_108026.root",
"skim_108047.root",
"skim_108051.root",
"skim_108063.root",
"skim_108080.root",
"skim_108092.root",
"skim_108137.root",
"skim_108147.root",
"skim_108190.root",
"skim_108205.root",
"skim_108216.root",
"skim_108254.root",
"skim_108295.root",
"skim_108311.root",
"skim_108312.root",
"skim_108353.root",
"skim_108367.root",
"skim_108370.root",
"skim_108383.root",
"skim_108387.root",
"skim_108418.root",
"skim_108419.root",
"skim_108427.root",
"skim_108428.root",
"skim_108449.root",
"skim_108454.root",
"skim_108460.root",
"skim_108514.root",
"skim_108528.root",
"skim_108533.root",
"skim_108539.root",
"skim_108557.root",
"skim_108563.root",
"skim_108579.root",
"skim_108580.root",
"skim_108583.root",
"skim_108596.root",
"skim_108601.root",
"skim_108613.root",
"skim_108623.root",
"skim_108642.root",
"skim_108650.root",
"skim_108659.root",
"skim_108663.root",
"skim_108675.root",
"skim_108678.root",
"skim_108697.root",
"skim_108707.root",
"skim_108738.root",
"skim_108741.root",
"skim_108748.root",
"skim_108768.root",
"skim_108783.root",
"skim_108785.root",
"skim_108789.root",
"skim_108795.root",
"skim_108807.root",
"skim_108815.root",
"skim_108817.root",
"skim_108841.root",
"skim_108842.root",
"skim_108844.root",
"skim_108860.root",
"skim_108873.root",
"skim_108875.root",
"skim_108880.root",
"skim_108886.root",
"skim_108909.root",
"skim_108914.root",
"skim_108922.root",
"skim_108933.root",
"skim_108950.root",
"skim_108987.root",
"skim_109010.root",
"skim_109015.root",
"skim_109031.root",
"skim_109035.root",
"skim_109040.root",
"skim_109043.root",
"skim_109044.root",
"skim_109051.root",
"skim_109062.root",
"skim_109063.root",
"skim_109096.root",
"skim_109097.root",
"skim_109103.root",
"skim_109113.root",
"skim_109126.root",
"skim_109128.root",
"skim_109135.root",
"skim_109155.root",
"skim_109165.root",
"skim_109170.root",
"skim_109174.root",
"skim_109186.root",
"skim_109188.root",
"skim_109197.root",
"skim_109199.root",
"skim_109209.root",
"skim_109212.root",
"skim_109216.root",
"skim_109232.root",
"skim_109233.root",
"skim_109239.root",
"skim_109254.root",
"skim_109257.root",
"skim_109261.root",
"skim_109270.root",
"skim_109302.root",
"skim_109331.root",
"skim_109351.root",
"skim_109365.root",
"skim_109369.root",
"skim_109383.root",
"skim_109389.root",
"skim_109401.root",
"skim_109404.root",
"skim_109411.root",
"skim_109450.root",
"skim_109451.root",
"skim_109467.root",
"skim_109480.root",
"skim_109481.root",
"skim_109501.root",
"skim_109505.root",
"skim_109514.root",
"skim_109517.root",
"skim_109527.root",
"skim_109537.root",
"skim_109544.root",
"skim_109558.root",
"skim_109562.root",
"skim_109582.root",
"skim_109583.root",
"skim_109611.root",
"skim_109615.root",
"skim_109617.root",
"skim_109624.root",
"skim_109628.root",
"skim_109629.root",
"skim_109630.root",
"skim_109631.root",
"skim_109643.root",
"skim_109645.root",
"skim_109648.root",
"skim_109650.root",
"skim_109671.root",
"skim_109676.root",
"skim_109683.root",
"skim_109702.root",
"skim_109705.root",
"skim_109722.root",
"skim_109735.root",
"skim_109749.root",
"skim_109752.root",
"skim_109753.root",
"skim_109765.root",
"skim_109772.root",
"skim_109775.root",
"skim_109797.root",
"skim_109798.root",
"skim_109804.root",
"skim_109806.root",
"skim_109808.root",
"skim_109813.root",
"skim_109816.root",
"skim_109827.root",
"skim_109828.root",
"skim_109835.root",
"skim_109836.root",
"skim_109838.root",
"skim_109854.root",
"skim_109855.root",
"skim_109861.root",
"skim_109862.root",
"skim_109880.root",
"skim_109892.root",
"skim_109896.root",
"skim_109902.root",
"skim_109905.root",
"skim_109917.root",
"skim_109928.root",
"skim_109938.root",
"skim_109939.root",
"skim_109951.root",
"skim_109961.root",
"skim_109977.root",
"skim_109988.root",
"skim_109991.root",
"skim_110007.root",
"skim_110014.root",
"skim_110021.root",
"skim_110023.root",
"skim_110025.root",
"skim_110038.root",
"skim_110052.root",
"skim_110059.root",
"skim_110060.root",
"skim_110073.root",
"skim_110097.root",
"skim_110109.root",
"skim_110110.root",
"skim_110116.root",
"skim_110131.root",
"skim_110133.root",
"skim_110138.root",
"skim_110148.root",
"skim_110161.root",
"skim_110163.root",
"skim_110167.root",
"skim_110171.root",
"skim_110185.root",
"skim_110187.root",
"skim_110189.root",
"skim_110199.root",
"skim_110202.root",
"skim_110211.root",
"skim_110225.root",
"skim_110248.root",
"skim_110265.root",
"skim_110294.root",
"skim_110305.root",
"skim_110307.root",
"skim_110313.root",
"skim_110322.root",
"skim_110332.root",
"skim_110339.root",
"skim_110346.root",
"skim_110374.root",
"skim_110375.root",
"skim_110379.root",
"skim_110382.root",
"skim_110385.root",
"skim_110390.root",
"skim_110391.root",
"skim_110399.root",
"skim_110409.root",
"skim_110412.root",
"skim_110420.root",
"skim_110423.root",
"skim_110438.root",
"skim_110440.root",
"skim_110441.root",
"skim_110466.root",
"skim_110474.root",
"skim_110505.root",
"skim_110524.root",
"skim_110530.root",
"skim_110542.root",
"skim_110546.root",
"skim_110562.root",
"skim_110577.root",
"skim_110583.root",
"skim_110588.root",
"skim_110599.root",
"skim_110606.root",
"skim_110632.root",
"skim_110647.root",
"skim_110670.root",
"skim_110675.root",
"skim_110677.root",
"skim_110707.root",
"skim_110712.root",
"skim_110714.root",
"skim_110729.root",
"skim_110757.root",
"skim_110767.root",
"skim_110785.root",
"skim_110787.root",
"skim_110794.root",
"skim_110832.root",
"skim_110833.root",
"skim_110845.root",
"skim_110846.root",
"skim_110886.root",
"skim_110889.root",
"skim_110890.root",
"skim_110918.root",
"skim_110925.root",
"skim_110927.root",
"skim_110939.root",
"skim_110958.root",
"skim_110973.root",
"skim_110974.root",
"skim_110994.root",
"skim_110997.root",
"skim_110999.root",
"skim_111003.root",
"skim_111004.root",
"skim_111006.root",
"skim_111009.root",
"skim_111018.root",
"skim_111019.root",
"skim_111022.root",
"skim_111025.root",
"skim_111050.root",
"skim_111055.root",
"skim_111062.root",
"skim_111064.root",
"skim_111066.root",
"skim_111079.root",
"skim_111094.root",
"skim_111099.root",
"skim_111111.root",
"skim_111116.root",
"skim_111118.root",
"skim_111121.root",
"skim_111126.root",
"skim_111129.root",
"skim_111151.root",
"skim_111161.root",
"skim_111164.root",
"skim_111194.root",
"skim_111217.root",
"skim_111236.root",
"skim_111238.root",
"skim_111265.root",
"skim_111274.root",
"skim_111284.root",
"skim_111286.root",
"skim_111303.root",
"skim_111313.root",
"skim_111314.root",
"skim_111315.root",
"skim_111321.root",
"skim_111377.root",
"skim_111379.root",
"skim_111395.root",
"skim_111418.root",
"skim_111428.root",
"skim_111440.root",
"skim_111444.root",
"skim_111447.root",
"skim_111452.root",
"skim_111455.root",
"skim_111489.root",
"skim_111491.root",
"skim_111518.root",
"skim_111521.root",
"skim_111523.root",
"skim_111533.root",
"skim_111551.root",
"skim_111563.root",
"skim_111573.root",
"skim_111576.root",
"skim_111590.root",
"skim_111616.root",
"skim_111635.root",
"skim_111638.root",
"skim_111671.root",
"skim_111683.root",
"skim_111688.root",
"skim_111705.root",
"skim_111719.root",
"skim_111731.root",
"skim_111762.root",
"skim_111767.root",
"skim_111784.root",
"skim_111797.root",
"skim_111814.root",
"skim_111817.root",
"skim_111819.root",
"skim_111835.root",
"skim_111857.root",
"skim_111859.root",
"skim_111867.root",
"skim_111869.root",
"skim_111877.root",
"skim_111908.root",
"skim_111918.root",
"skim_111940.root",
"skim_111951.root",
"skim_111952.root",
"skim_111954.root",
"skim_111957.root",
"skim_111959.root",
"skim_111968.root",
"skim_111974.root",
"skim_111976.root",
"skim_111978.root",
"skim_111982.root",
"skim_111989.root",
"skim_111994.root",
"skim_112005.root",
"skim_112012.root",
"skim_112013.root",
"skim_112035.root",
"skim_112037.root",
"skim_112039.root",
"skim_112054.root",
"skim_112055.root",
"skim_112060.root",
"skim_112067.root",
"skim_112070.root",
"skim_112078.root",
"skim_112106.root",
"skim_112108.root",
"skim_112123.root",
"skim_112146.root",
"skim_112160.root",
"skim_112164.root",
"skim_112173.root",
"skim_112175.root",
"skim_112189.root",
"skim_112203.root",
"skim_112204.root",
"skim_112205.root",
"skim_112214.root",
"skim_112219.root",
"skim_112220.root",
"skim_112233.root",
"skim_112235.root",
"skim_112257.root",
"skim_112274.root",
"skim_112288.root",
"skim_112300.root",
"skim_112301.root",
"skim_112317.root",
"skim_112322.root",
"skim_112327.root",
"skim_112329.root",
"skim_112334.root",
"skim_112338.root",
"skim_112347.root",
"skim_112362.root",
"skim_112372.root",
"skim_112395.root",
"skim_112396.root",
"skim_112400.root",
"skim_112401.root",
"skim_112405.root",
"skim_112407.root",
"skim_112425.root",
"skim_112426.root",
"skim_112435.root",
"skim_112440.root",
"skim_112441.root",
"skim_112444.root",
"skim_112452.root",
"skim_112472.root",
"skim_112480.root",
"skim_112508.root",
"skim_112511.root",
"skim_112518.root",
"skim_112538.root",
"skim_112546.root",
"skim_112551.root",
"skim_112571.root",
"skim_112572.root",
"skim_112585.root",
"skim_112586.root",
"skim_112588.root",
"skim_112593.root",
"skim_112625.root",
"skim_112629.root",
"skim_112630.root",
"skim_112651.root",
"skim_112653.root",
"skim_112662.root",
"skim_112670.root",
"skim_112683.root",
"skim_112696.root",
"skim_112700.root",
"skim_112715.root",
"skim_112717.root",
"skim_112720.root",
"skim_112723.root",
"skim_112729.root",
"skim_112739.root",
"skim_112752.root",
"skim_112762.root",
"skim_112770.root",
"skim_112785.root",
"skim_112813.root",
"skim_112818.root",
"skim_112831.root",
"skim_112846.root",
"skim_112884.root",
"skim_112896.root",
"skim_112903.root",
"skim_112905.root",
"skim_112912.root",
"skim_112915.root",
"skim_112929.root",
"skim_112930.root",
"skim_112931.root",
"skim_112933.root",
"skim_112955.root",
"skim_112961.root",
"skim_112965.root",
"skim_113001.root",
"skim_113012.root",
"skim_113024.root",
"skim_113053.root",
"skim_113088.root",
"skim_113096.root",
"skim_113101.root",
"skim_113113.root",
"skim_113114.root",
"skim_113133.root",
"skim_113140.root",
"skim_113143.root",
"skim_113174.root",
"skim_113180.root",
"skim_113191.root",
"skim_113195.root",
"skim_113201.root",
"skim_113203.root",
"skim_113205.root",
"skim_113213.root",
"skim_113217.root",
"skim_113238.root",
"skim_113264.root",
"skim_113276.root",
"skim_113279.root",
"skim_113280.root",
"skim_113288.root",
"skim_113290.root",
"skim_113303.root",
"skim_113308.root",
"skim_113312.root",
"skim_113321.root",
"skim_113326.root",
"skim_113327.root",
"skim_113330.root",
"skim_113333.root",
"skim_113363.root",
"skim_113364.root",
"skim_113366.root",
"skim_113368.root",
"skim_113376.root",
"skim_113389.root",
"skim_113397.root",
"skim_113407.root",
"skim_113428.root",
"skim_113435.root",
"skim_113438.root",
"skim_113442.root",
"skim_113444.root",
"skim_113452.root",
"skim_113457.root",
"skim_113460.root",
"skim_113466.root",
"skim_113473.root",
"skim_113484.root",
"skim_113503.root",
"skim_113505.root",
"skim_113512.root",
"skim_113524.root",
"skim_113536.root",
"skim_113540.root",
"skim_113559.root",
"skim_113562.root",
"skim_113576.root",
"skim_113581.root",
"skim_113596.root",
"skim_113597.root",
"skim_113612.root",
"skim_113623.root",
"skim_113654.root",
"skim_113658.root",
"skim_113679.root",
"skim_113697.root",
"skim_113708.root",
"skim_113717.root",
"skim_113729.root",
"skim_113749.root",
"skim_113750.root",
"skim_113788.root",
"skim_113807.root",
"skim_113817.root",
"skim_113819.root",
"skim_113830.root",
"skim_113849.root",
"skim_113868.root",
"skim_113879.root",
"skim_113892.root",
"skim_113900.root",
"skim_113910.root",
"skim_113917.root",
"skim_113921.root",
"skim_113944.root",
"skim_113954.root",
"skim_113972.root",
"skim_113982.root",
"skim_113985.root",
"skim_113990.root",
"skim_114000.root",
"skim_114005.root",
"skim_114021.root",
"skim_114026.root",
"skim_114029.root",
"skim_114030.root",
"skim_114039.root",
"skim_114050.root",
"skim_114061.root",
"skim_114062.root",
"skim_114074.root",
"skim_114078.root",
"skim_114087.root",
"skim_114091.root",
"skim_114110.root",
"skim_114121.root",
"skim_114125.root",
"skim_114182.root",
"skim_114206.root",
"skim_114236.root",
"skim_114237.root",
"skim_114240.root",
"skim_114244.root",
"skim_114260.root",
"skim_114265.root",
"skim_114266.root",
"skim_114267.root",
"skim_114288.root",
"skim_114301.root",
"skim_114307.root",
"skim_114324.root",
"skim_114339.root",
"skim_114344.root",
"skim_114346.root",
"skim_114358.root",
"skim_114360.root",
"skim_114361.root",
"skim_114391.root",
"skim_114408.root",
"skim_114445.root",
"skim_114462.root",
"skim_114476.root",
"skim_114495.root",
"skim_114497.root",
"skim_114519.root",
"skim_114539.root",
"skim_114542.root",
"skim_114565.root",
"skim_114583.root",
"skim_114585.root",
"skim_114587.root",
"skim_114595.root",
"skim_114600.root",
"skim_114601.root",
"skim_114609.root",
"skim_114638.root",
"skim_114642.root",
"skim_114646.root",
"skim_114657.root",
"skim_114660.root",
"skim_114661.root",
"skim_114670.root",
"skim_114671.root",
"skim_114713.root",
"skim_114727.root",
"skim_114739.root",
"run000347-ana_vertexfit.root",
"run000365-ana_vertexfit.root",
"run000402-ana_vertexfit.root",
"run000435-ana_vertexfit.root",
"run000495-ana_vertexfit.root",
"run000820-ana_vertexfit.root",
"run000822-ana_vertexfit.root",
"run000840-ana_vertexfit.root",
"run000848-ana_vertexfit.root",
"run000873-ana_vertexfit.root",
"run000906-ana_vertexfit.root",
"run000916-ana_vertexfit.root",
"run000917-ana_vertexfit.root",
"run000926-ana_vertexfit.root",
"run000930-ana_vertexfit.root",
"run000933-ana_vertexfit.root",
"run000941-ana_vertexfit.root",
"run000943-ana_vertexfit.root",
"run000950-ana_vertexfit.root",
"run000953-ana_vertexfit.root",
"run000957-ana_vertexfit.root",
"run000960-ana_vertexfit.root",
"run000967-ana_vertexfit.root",
"run000968-ana_vertexfit.root",
"run000971-ana_vertexfit.root",
"run000984-ana_vertexfit.root",
"run000991-ana_vertexfit.root",
"run000993-ana_vertexfit.root",
"run000995-ana_vertexfit.root",
"run000996-ana_vertexfit.root",
    };

    for (int i = 50000+startNum; i < 50000+startNum + numFiles; ++i) {
      // Change filename to names in specific directory (i.e. Signal, Bhabha, IC)
        std::string fileName = "run0" + std::to_string(i) + "-ana_vertexfit.root";
        std::string filePath = folderPath + "/" + fileName;

        // Check if the file should be excluded
        bool excludeFile = false;
        for (const std::string& excludedFile : excludedFiles) {
            if (fileName == excludedFile) {
                excludeFile = true;
                break;
            }
        }

        if (excludeFile) {
            std::cout << "File " << fileName << " excluded from chain." << std::endl;
            continue; // Skip this file
        }

        chain.Add(filePath.c_str(), 0);
        std::cout << "File " << fileName << " added to chain." << std::endl;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <input_folder> <output_file_name> <start_num> <num_files>" << std::endl;
        return 1;
    }

    std::string folderPath = argv[1];
    std::string outputFileName = argv[2];
    int startNum = std::atoi(argv[3]);
    int numFiles = std::atoi(argv[4]);

    TChain chain("vertex");
    addSkimFilesToChain(folderPath, chain, startNum, numFiles);

    SignalAna * a = new SignalAna(&chain, outputFileName, "vertex");
    a->Loop();
    std::cout << "Loop completed" << std::endl;
    delete a;

    return 0;
}


