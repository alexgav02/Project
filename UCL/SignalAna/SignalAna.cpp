#include "SignalAna.hpp"
#include <iomanip>
#include "TAxis.h"
#include "TH2.h"
#include <cmath>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>

using std::cout;
using std::endl;


//this function defines the histograms we wiil fill:
void SignalAna::InitHistos() {

  plot1D("nvert", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("pm", "; vertex momentum [MeV] ; entries", 100, 0,100);
  plot1D("mmass", "; vertex mass [MeV] ; entries", 150, 0,150);
  plot1D("chi2", "; #chi^2 ; entries", 100, 0,100);
  plot1D("nrecurl", "; # recurl tracks ; entries", 3, -0.5,2.5);


  ////////////////////////////
  //Applying all cuts in TDR//
  ////////////////////////////

  plot1D("mmass_TDRcuts", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("pm_TDRcuts", "; p_{rec} [MeV] ; entries", 40,0,4);

  plot1D("p1_TDR", "; p_{1} [MeV/c] ; entries", 100, 0, 100);
  plot1D("p2_TDR", "; p_{2} [MeV/c] ; entries", 100, 0, 100);
  plot1D("p3_TDR", "; p_{3} [MeV/c] ; entries", 100, 0, 100);

  plot1D("chi2_TDRcuts", "; #chi^2 ; entries", 100, 0,100);
  plot1D("nhit1_TDRcuts", "; # hits for particle 1 ; entries", 13, -0.5, 12.5);
  plot1D("nhit2_TDRcuts", "; # hits for particle 1 ; entries", 13, -0.5, 12.5);
  plot1D("nhit3_TDRcuts", "; # hits for particle 1 ; entries", 13, -0.5, 12.5);



  plot1D("mmass_TDR_to_normalise","; m_{rec} [MeV] ; entries", 150, 95, 110);
  plot1D("pm_TDR_to_normalise1", "; p_{rec} [MeV] ; entries", 200,0,20);
  plot1D("pm_TDR_to_normalise2", "; p_{rec} [MeV] ; entries", 200,0,20);
  plot1D("pm_TDR_to_normalise3", "; p_{rec} [MeV] ; entries", 200,0,20);
  plot1D("pm_TDR_to_normalise4", "; p_{rec} [MeV] ; entries", 200,0,20);


  plot1D("mmass_manualTDR", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("pm_manualTDR", "; p_{rec} [MeV] ; entries", 40,0,4);

  plot1D("mmass_TDR_without_weights","; m_{rec} [MeV] ; entries", 150, 95, 110);
  plot1D("pm_TDR_without_weights","; p_{rec} [MeV] ; entries", 200,0,20);


  plot1D("nsharedhit1_TDRcuts", "; # shared hits for particle 1 ; entries", 9, -0.5, 8.5);
  plot1D("nsharedhit2_TDRcuts", "; # shared hits for particle 3 ; entries", 9, -0.5, 8.5);
  plot1D("nsharedhit3_TDRcuts", "; # shared hits for particle 3 ; entries", 9, -0.5, 8.5);


  plot1D("lambda1_TDRcuts", "; Angle #lambda for particle 1 ; entries",100,-M_PI/2,M_PI/2);
  plot1D("lambda2_TDRcuts", "; Angle #lambda for particle 2 ; entries",100,-M_PI/2,M_PI/2);
  plot1D("lambda3_TDRcuts", "; Angle #lambda for particle 3 ; entries",100,-M_PI/2,M_PI/2);

  plot1D("phi1_TDRcuts", "; Angle |#phi| for particle 1  ; entries",100,-M_PI,M_PI);
  plot1D("phi2_TDRcuts", "; Angle |#phi| for particle 2 ; entries",100,-M_PI,M_PI);
  plot1D("phi3_TDRcuts", "; Angle |#phi| for particle 3 ; entries",100,-M_PI,M_PI);

  plot1D("pT1_TDRcuts", "; p_{T,1} [MeV/c] ; entries", 60, 0, 60);
  plot1D("pT2_TDRcuts", "; p_{T,2} [MeV/c] ; entries", 60, 0, 60);
  plot1D("pT3_TDRcuts", "; p_{T,3} [MeV/c] ; entries", 60, 0, 60);

  plot2D("mmass_vs_lambda1_TDRcuts", "; m_{rec} [MeV] ; Angle #lambda between particle 1 track and transverse plane ; entries",70, 103, 110, 180,-90,90);
  plot2D("mmass_vs_lambda2_TDRcuts", "; m_{rec} [MeV] ; Angle #lambda between particle 2 track and transverse plane ; entries",70, 103, 110, 180,-90,90);
  plot2D("mmass_vs_lambda3_TDRcuts", "; m_{rec} [MeV] ; Angle #lambda between particle 3 track and transverse plane ; entries",70, 103, 110, 180,-90,90);

  plot2D("mmass_vs_momentum_TDRcuts", "; m_{rec} [MeV] ; p_{rec} ; entries", 70, 103, 110, 40, 0, 4);

  plot2D("chi2_vs_mass_TDRcuts","; mmass [MeV] ; #chi2 ; entries", 70, 103, 110, 15, 0, 15);
  plot2D("chi2_vs_momentum_TDRcuts", "; pm_{rec} [MeV] ; #chi2 ; entries", 40, 0, 4, 15, 0, 15);

  plot1D("weight_TDRcuts", "; weight per vertex ; entries", 1000, 0, 1000);


  plot1D("vx_res_TDRcuts", "; x_{rec} - x_{true} [mm] ; entries",100, -5,5);
  plot1D("vy_res_TDRcuts", "; y_{rec} - y_{true} [mm] ; entries",100, -5,5);
  plot1D("vz_res_TDRcuts", "; z_{rec} - z_{true} [mm] ; entries",100, -5,5);
  plot1D("pm_res_TDRcuts", "; |pm_{rec} - pm_{true}| [MeV]; entries", 250, 0, 50);
  plot1D("p1_res_TDRcuts", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p2_res_TDRcuts", ";|pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p3_res_TDRcuts", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("target_dist_res_TDRcuts", "; d_{target, rec} - d_{target, true} [mm] ; entries", 100, -5, 5);
  plot1D("mmass_res_TDRcuts", "; |m_{rec} - m_{true}| [MeV] ; entries", 250, 0, 50);

  plot1D("true_p1_TDRcuts", "; p_{rec} [MeV] ; entries", 100,0,100);
  plot1D("true_p2_TDRcuts", "; p_{rec} [MeV] ; entries", 100,0,100);
  plot1D("true_p3_TDRcuts", "; p_{rec} [MeV] ; entries", 100,0,100);


  ////////////////////////////////////////
  //Investgating second peak in momentum//
  ////////////////////////////////////////


  plot1D("pm_secondpeak","; vertex momentum [MeV] ; entries", 90,10,100);
  plot1D("mmass_secondpeak","; vertex mass [MeV] ; entries", 150, 0,150);
  plot1D("chi_secondpeak", "; #chi^2 ; entries", 100, 0, 100);
  plot1D("nrecurl_secondpeak","; # recurl tracks ; entries",3,-0.5,2.5);
  plot1D("pxm_secondpeak", "; p_{x, cms} [MeV] ; entries", 100, -10, 10);
  plot1D("pym_secondpeak", "; p_{y, cms} [MeV] ; entries", 100, -10, 10);
  plot1D("pzm_secondpeak", "; p_{z, cms} [MeV] ; entries", 100, -10, 10);


  ////////////////////////////////////////////////
  //Plotting reconstructed parameter resolutions//
  ////////////////////////////////////////////////
  

  plot1D("vx_res", "; x_{rec} - x_{true} [mm] ; entries",100, -5,5);
  plot1D("vy_res", "; y_{rec} - y_{true} [mm] ; entries",100, -5,5);
  plot1D("vz_res", "; z_{rec} - z_{true} [mm] ; entries",100, -5,5);
  plot1D("pm_res", "; |pm_{rec} - pm_{true}| [MeV]; entries", 250, 0, 50);
  plot1D("target_dist_res", "; d_{target, rec} - d_{target, true} [mm] ; entries", 100, -5, 5);
  plot1D("mmass_res", "; |m_{rec} - m_{true}| [MeV] ; entries", 250, 0, 50);

  //////////////////////////////////////////////////
  //Plotting muon momentum components in COM frame//
  //////////////////////////////////////////////////

  plot1D("px_cms", "; p_{x,cms} [MeV/c] ; entries", 100, -10, 10);
  plot1D("py_cms", "; p_{y,cms} [MeV/c] ; entries", 100, -10, 10);
  plot1D("pz_cms", "; p_{z,cms} [MeV/c] ; entries", 100, -10, 10);

  /////////////////////////////////////////////////
  ////Plotting z components of particle momenta////
  /////////////////////////////////////////////////

  plot1D("pz1", "; p_{z,1} [MeV/c] ; entries", 120, -60, 60);
  plot1D("pz2", "; p_{z,2} [MeV/c] ; entries", 120, -60, 60);
  plot1D("pz3", "; p_{z,3} [MeV/c] ; entries", 120, -60, 60);


  ///////////////////////////////////
  //Plotting particle total momenta//
  ///////////////////////////////////

  plot1D("p1", "; p_{1} [MeV/c] ; entries", 100, 0, 100);
  plot1D("p2", "; p_{2} [MeV/c] ; entries", 100, 0, 100);
  plot1D("p3", "; p_{3}, [MeV/c] ; entries", 100, 0, 100);


  /////////////////////////////////////
  //Plotting pcoplanar and pacoplanar//
  /////////////////////////////////////

  plot1D("pcoplanar", "; p_{coplanar} [MeV] ; entries", 100, 0, 100);
  plot1D("pacoplanar", "; p_{acoplanar} [MeV] ; entries", 140, -70, 70);


  ////////////////////////////////////////////////////////////////////////////
  //Plotting reconstructed parameters allowing different tracks with 4 hits///
  ////////////////////////////////////////////////////////////////////////////
  
  plot1D("mmass_4hit_1", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("mmass_4hit_2", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("mmass_4hit_3", "; m_{rec} [MeV] ; entries", 70, 103, 110);


  plot1D("pm_4hit_1", "; p_{rec} [MeV] ; entries", 40, 0, 4);
  plot1D("pm_4hit_2", "; p_{rec} [MeV] ; entries", 40, 0, 4);
  plot1D("pm_4hit_3", "; p_{rec} [MeV] ; entries", 40, 0, 4);

  plot1D("pm_res_4hit_1", "; |pm_{rec} - pm_{true}| [MeV]; entries", 250, 0, 50);
  plot1D("pm_res_4hit_2", "; |pm_{rec} - pm_{true}| [MeV]; entries", 250, 0, 50);
  plot1D("pm_res_4hit_3", "; |pm_{rec} - pm_{true}| [MeV]; entries", 250, 0, 50);

  plot1D("p1_res_4hit_1", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p2_res_4hit_1", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p3_res_4hit_1", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p1_res_4hit_2", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p2_res_4hit_2", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p3_res_4hit_2", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p1_res_4hit_3", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p2_res_4hit_3", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p3_res_4hit_3", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("mmass_res_4hit_1", "; |m_{rec} - m_{true}| [MeV] ; entries", 250, 0, 50);
  plot1D("mmass_res_4hit_2", "; |m_{rec} - m_{true}| [MeV] ; entries", 250, 0, 50);
  plot1D("mmass_res_4hit_3", "; |m_{rec} - m_{true}| [MeV] ; entries", 250, 0, 50);


  plot1D("nvert_4hit_1", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_4hit_2", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_4hit_3", "; # vertices ; entries", 100, -0.5,99.5);



  plot1D("mmass_4hit_12", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("mmass_4hit_13", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("mmass_4hit_23", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("nvert_4hit_12", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_4hit_13", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_4hit_23", "; # vertices ; entries", 100, -0.5,99.5);


  plot1D("pm_4hit_12", "; p_{rec} [MeV] ; entries", 40, 0, 4);
  plot1D("pm_4hit_13", "; p_{rec} [MeV] ; entries", 40, 0, 4);
  plot1D("pm_4hit_23", "; p_{rec} [MeV] ; entries", 40, 0, 4);

  plot1D("pm_res_4hit_12", ";|pm_{rec} - pm_{true}| [MeV]; entries", 250, 0, 50);
  plot1D("pm_res_4hit_13", "; |pm_{rec} - pm_{true}| [MeV]; entries", 250, 0, 50);
  plot1D("pm_res_4hit_23", "; |pm_{rec} - pm_{true}| [MeV]; entries", 250, 0, 50);

  plot1D("p1_res_4hit_12", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p2_res_4hit_12", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p3_res_4hit_12", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p1_res_4hit_13", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p2_res_4hit_13", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p3_res_4hit_13", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p1_res_4hit_23", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p2_res_4hit_23", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p3_res_4hit_23", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);

  plot1D("mmass_res_4hit_12", "; |m_{rec} - m_{true}| [MeV] ; entries", 250, 0, 50);
  plot1D("mmass_res_4hit_13", "; |m_{rec} - m_{true}| [MeV] ; entries", 250, 0, 50);
  plot1D("mmass_res_4hit_23", "; |m_{rec} - m_{true}| [MeV] ; entries", 250, 0, 50);




  plot1D("mmass_4hit_123", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("pm_4hit_123", "; p_{rec} [MeV] ; entries", 40, 0, 4);

  plot1D("pm_res_4hit_123", "; |pm_{rec} - pm_{true}| [MeV]; entries", 250, 0, 50);
  plot1D("p1_res_4hit_123", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p2_res_4hit_123", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("p3_res_4hit_123", "; |pm_{rec} - pm_{true}| [MeV]; entries", 100, 0, 20);
  plot1D("mmass_res_4hit_123", "; |m_{rec} - m_{true}| [MeV] ; entries", 250, 0, 50);
  plot1D("nvert_4hit_123", "; # vertices ; entries", 100, -0.5,99.5);

  /////////////////////////////////////////////
  //Plotting number of hits for each particle//
  /////////////////////////////////////////////

  plot1D("nhits1", "; # hits for particle 1 ; entries", 5, 3.5, 8.5);
  plot1D("nhits2", "; # hits for particle 3 ; entries", 5, 3.5, 8.5);
  plot1D("nhits3", "; # hits for particle 3 ; entries", 5, 3.5, 8.5);


 /////////////////////////////////////////////////////////////////////////////
 //Plotting transverse momenta against chi^2 and reconstructed muon momentum// 
 /////////////////////////////////////////////////////////////////////////////
 
  plot2D("pT1_vs_chi2", "; p_{T,1} [MeV] ; #chi2 ; entries",100, 0, 60, 100, 0, 100);
  plot2D("pT2_vs_chi2", "; p_{T,2} [MeV] ; #chi2 ; entries",100, 0, 60, 100, 0, 100);
  plot2D("pT3_vs_chi2", "; p_{T,3} [MeV] ; #chi2 ; entries",100, 0, 60, 100, 0, 100); 

  plot2D("pT1_vs_pm", "; p_{T,1} [MeV] ; pm ; entries",60, 0, 60, 30, 0, 30);
  plot2D("pT2_vs_pm", "; p_{T,2} [MeV] ; pm ; entries",60, 0, 60, 30, 0, 30);
  plot2D("pT3_vs_pm", "; p_{T,3} [MeV] ; pm ; entries",60, 0, 60, 30, 0, 30);

  ///////////////////////////////////////////////////////////////////////
  //Plotting chi^2 for tracks with transverse momentum less than 10 MeV//
  ///////////////////////////////////////////////////////////////////////

  plot1D("chi2_pT1<10", "; #chi^2 ; entries", 100, 0,10000);
  plot1D("chi2_pT2<10", "; #chi^2 ; entries", 100, 0,10000);
  plot1D("chi2_pT3<10", "; #chi^2 ; entries", 100, 0,10000);

  plot1D("pT1_nocuts", "; p_{T,1} [MeV] ; entries", 100, 0, 100);
  plot1D("pT2_nocuts", "; p_{T,1} [MeV] ; entries", 100, 0, 100);
  plot1D("pT3_nocuts", "; p_{T,1} [MeV] ; entries", 100, 0, 100);


  plot1D("pT1_chi100", "; p_{T,1} [MeV] ; entries", 100, 0, 100);
  plot1D("pT2_chi100", "; p_{T,2} [MeV] ; entries", 100, 0, 100);
  plot1D("pT3_chi100", "; p_{T,3} [MeV] ; entries", 100, 0, 100);


  /////////////////////////////////////////////////////////////////////////////
  //Plotting nsharedhits for tracks with transverse momentum less than 10 MeV//
  /////////////////////////////////////////////////////////////////////////////

  plot1D("nsharedhits_particle1_pT1<10", "; # shared hits for particle 1 ; entries", 20, -0.5, 19.5);
  plot1D("nsharedhits_particle2_pT2<10", "; # shared hits for particle 2 ; entries", 20, -0.5, 19.5);
  plot1D("nsharedhits_particle3_pT3<10", "; # shared hits for particle 3 ; entries", 20, -0.5, 19.5);

  plot1D("nsharedseg_particle1_pT1<10", "; # shared segments for particle 1 ; entries", 9, -0.5, 8.5);
  plot1D("nsharedseg_particle2_pT1<10", "; # shared segments for particle 2 ; entries", 9, -0.5, 8.5);
  plot1D("nsharedseg_particle3_pT1<10", "; # shared segments for particle 3 ; entries", 9, -0.5, 8.5);
 
  /////////////////////////////////////////
  //Investigating TDR cuts on backgrounds//
  /////////////////////////////////////////

  plot1D("mmass_TDR_loose_mass", " ; m_{rec} [MeV] ; entries", 35, 80, 115);
  plot1D("pm_TDR_loose_mass", " ; p_{rec} [MeV] ; entries", 40, 0, 4);

  plot1D("mmass_TDR_without_mass_cuts", " ; m_{rec} [MeV] ; entries", 150, 0, 150);
  plot1D("pm_TDR_without_mass_cuts", " ; p_{rec} [MeV] ; entries", 40, 0, 4);

  plot1D("mmass_TDR_without_mass_momentum_mep1_cuts", " ; m_{rec} [MeV] ; entries", 150, 0, 150);
  plot1D("pm_TDR_without_mass_momentum_mep1_cuts", " ; p_{rec} [MeV] ; entries", 100, 0, 100);

  plot1D("mmass_TDR_cuts_chi2<100", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("pm_TDR_cuts_chi2<100", "; p_{rec} [MeV] ; entries", 40, 0, 4);
  
  plot1D("mmass_TDR_without_mass_momentum_cuts", "; m_{rec} [MeV] ; entries", 150, 0, 150 );
  plot1D("pm_TDR_without_mass_momentum_cuts", "; p_{rec} [MeV] ; entries", 100, 0, 100);
  plot1D("lambda1_TDR_without_mass_momentum_cuts", "; Angle #lambda for particle 1 ; entries", 100,-M_PI/2,M_PI/2);
  plot1D("lambda2_TDR_without_mass_momentum_cuts", "; Angle #lambda for particle 2 ; entries", 100,-M_PI/2,M_PI/2);
  plot1D("lambda3_TDR_without_mass_momentum_cuts", "; Angle #lambda for particle 3 ; entries", 100,-M_PI/2,M_PI/2);
  plot1D("phi1_TDR_without_mass_momentum_cuts", "; Angle |#phi| for particle 1  ; entries",100,-M_PI,M_PI);
  plot1D("phi2_TDR_without_mass_momentum_cuts", "; Angle |#phi| for particle 2  ; entries",100,-M_PI,M_PI);
  plot1D("phi3_TDR_without_mass_momentum_cuts", "; Angle |#phi| for particle 3  ; entries",100,-M_PI,M_PI);

  
  plot1D("mmass_only_chi2_cut", "; m_{rec} [MeV] ; entries", 150, 0, 150);
  plot1D("pm_only_chi2_cut", "; p_{rec} [MeV] ; entries", 100, 0, 100);

  plot1D("mmass_TDR_without_momentum_cuts", "; m_{rec} [MeV] ; entries", 70, 103, 110);
  plot1D("pm_TDR_without_momentum_cuts", "; p_{rec} [MeV] ; entries", 100, 0, 100);

  plot1D("mmass_TDR_mmass>80||<110", "; m_{rec} [MeV] ; entries", 30, 80, 110);
  plot1D("pm_TDR_mmass>80||<110", "; m_{rec} [MeV] ; entries", 40, 0, 4);
  plot1D("weight_TDR_mmass>80||<110", "; weight per vertex ; entries", 1000, 0, 1000);
  
// Note: Plot for all cuts except momentum (IC typically has higher momentum distribution)

////////////////////
//Plotting weights//
////////////////////

  plot1D("weights", "; weight per vertex ; entries", 1000,0,1000);

  plot2D("mmass_vs_weight", "; m_{rec} [MeV] ; weight ; entries",30, 80, 110, 1000,0,1000);

  plot2D("mmass_vs_mep1_cuts", "; m_{e^{+}e^{-}e^{+}} [MeV] ; m_{e^{+}e^{-}, small} [MeV]", 150, 95, 110, 300, 0, 30);
  plot2D("mmass_vs_mep1", "; m_{e^{+}e^{-}e^{+}} [MeV] ; m_{e^{+}e^{-}, small} [MeV]", 150, 95, 110, 300, 0, 30);


  plot1D("nvert_chi2_cut", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_mass_cut", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_mep1_cut", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_mom_cut", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_nhit_cut", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_massloose_cut", "; # vertices ; entries", 100, -0.5,99.5);


  plot1D("nvert_faketracks", "; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_realtracks", "; # vertices ; entries", 100, -0.5,99.5);


  plot1D("nvert_chi2_mass","; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_chi2_mass_mom","; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_chi2_mass_mom_mep1","; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_TDR","; # vertices ; entries", 100, -0.5,99.5);

  plot1D("nvert_TDR_4hit_1","; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_TDR_4hit_2","; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_TDR_4hit_3","; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_TDR_4hit_12","; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_TDR_4hit_13","; # vertices ; entries", 100, -0.5,99.5);
  plot1D("nvert_TDR_4hit_23","; # vertices ; entries", 100, -0.5,99.5);

  


}




//======================================================================
//======================================================================
//======================================================================

//this function then loops over the entries in the input file:
void SignalAna::Loop()
{
   if (fChain == 0) return;

  //here we loop over the tree, loading each frame in turn:
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries_;jentry++) {
     //for (Long64_t jentry=0; jentry<1000;jentry++) {
     ReadEvent(jentry);

     
     // work out how many vertices there are in this frame:
     int nvert = (*chi2).size();
     Fill("nvert", nvert);

     //now loop over the vertices in this frame:
     for(int i=0; i<nvert; i++) {

      //access some basic variables:
       const double ipm = (*pm)[i];
       const double ichi2 = (*chi2)[i];
       const double immass = (*mmass)[i];
       const int inrecurl = ((*nhit1)[i]>4) + ((*nhit2)[i]>4) + ((*nhit3)[i]>4);
       const double ipxm = (*pxm)[i];
       const double ipym = (*pym)[i];
       const double ipzm = (*pzm)[i];
       const double ipx1 = (*px1)[i];
       const double ipy1 = (*py1)[i];
       const double ipz1 = (*pz1)[i];
       const double ipx2 = (*px2)[i];
       const double ipy2 = (*py2)[i];
       const double ipz2 = (*pz2)[i];
       const double ipx3 = (*px3)[i];
       const double ipy3 = (*py3)[i];
       const double ipz3 = (*pz3)[i];
       const int insharedhit1 = (*nsharedhit1)[i];
       const int insharedhit2 = (*nsharedhit2)[i];
       const int insharedhit3 = (*nsharedhit3)[i];
       const int insharedhit_total = (insharedhit1 + insharedhit2 + insharedhit3);
       const double ivx_res = ((*vx)[i]-(*true_vx)[i]);
       const double ivy_res = ((*vy)[i]-(*true_vy)[i]);
       const double ivz_res = ((*vz)[i]-(*true_vz)[i]);
       const double itrue_pxm = (*true_pxm)[i];
       const double itrue_pym = (*true_pym)[i];
       const double itrue_pzm = (*true_pzm)[i];
       const double itrue_pm = sqrt(itrue_pxm*itrue_pxm + itrue_pym*itrue_pym + itrue_pzm*itrue_pzm);
       const double ipm_res = abs(((*pm)[i]-itrue_pm));
       const double itrue_px1 = (*true_px1)[i];
       const double itrue_px2 = (*true_py2)[i];
       const double itrue_px3 = (*true_pz3)[i];
       const double itrue_py1 = (*true_px1)[i];
       const double itrue_py2 = (*true_py2)[i];
       const double itrue_py3 = (*true_pz3)[i];
       const double itrue_pz1 = (*true_px1)[i];
       const double itrue_pz2 = (*true_py2)[i];
       const double itrue_pz3 = (*true_pz3)[i];
       const double itrue_p1 = sqrt(itrue_px1*itrue_px1 + itrue_py1*itrue_py1 + itrue_pz1*itrue_pz1);
       const double itrue_p2 = sqrt(itrue_px2*itrue_px2 + itrue_py2*itrue_py2 + itrue_pz2*itrue_pz2);
       const double itrue_p3 = sqrt(itrue_px3*itrue_px3 + itrue_py3*itrue_py3 + itrue_pz3*itrue_pz3);
       const double ip1 = sqrt(ipx1*ipx1 + ipy1*ipy1 + ipz1*ipz1);
       const double ip2 = sqrt(ipx2*ipx2 + ipy2*ipy2 + ipz2*ipz2);
       const double ip3 = sqrt(ipx3*ipx3 + ipy3*ipy3 + ipz3*ipz3);   
       const double ip1_res = (itrue_p1-ip1);
       const double ip2_res = (itrue_p2-ip2);
       const double ip3_res = (itrue_p3-ip3);
       const double itarget_dist_res = ((*targetdist)[i] - (*true_targetdist)[i]);
       const double immass_res = abs(((*mmass)[i] - (*true_mass)[i]));
       const double imep1 = (*mep1)[i];
       const double imep2 = (*mep2)[i];
       const double ipT1 = sqrt(ipx1*ipx1 + ipy1*ipy1);
       const double ipT2 = sqrt(ipx2*ipx2 + ipy2*ipy2);
       const double ipT3 = sqrt(ipx3*ipx3 + ipy3*ipy3);
       const double ilambda1 = atan(ipz1/ipT1);
       const double ilambda2 = atan(ipz2/ipT2);
       const double ilambda3 = atan(ipz3/ipT3);
       const double iphi1 = atan2(ipx1,ipy1);
       const double iphi2 = atan2(ipx2,ipy2);
       const double iphi3 = atan2(ipx3,ipy3);
       const double ipcoplanar = (*pcoplanar)[i];
       const double ipacoplanar = (*pacoplanar)[i];
       const int cutStatus = SetRecoCuts(i);
       const int inhit1 = (*nhit1)[i];
       const int inhit2 = (*nhit2)[i];
       const int inhit3 = (*nhit3)[i];
       const int iid1 = (*id1)[i];
       const int iid2 = (*id2)[i];
       const int iid3 = (*id3)[i];
       const int itype1 = (*type1)[i];
       const int itype2 = (*type2)[i];
       const int itype3 = (*type3)[i];
       const int insharedseg1 = (*nsharedseg1)[i];
       const int insharedseg2 = (*nsharedseg2)[i];
       const int insharedseg3 = (*nsharedseg3)[i];
       
      if (weight < 1e+06){
       ////////////////////////////////////////////////////////////////////////////////////
       //Plotting reconstructed parameters while allowing different tracks to have 4 hits//
       ////////////////////////////////////////////////////////////////////////////////////

       // Plotting for 1 track with 4 hits allowed

        if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) &&  didCutPass(NHIT2, cutStatus) &&  didCutPass(NHIT3, cutStatus) ) {
          Fill("mmass_4hit_1", immass, weight);
	        Fill("pm_4hit_1", ipm, weight);
          Fill("pm_res_4hit_1", ipm_res);
          Fill("p1_res_4hit_1", ip1_res);
          Fill("p2_res_4hit_1", ip2_res);
          Fill("p3_res_4hit_1", ip3_res);
          Fill("mmass_res_4hit_1", immass_res);
          Fill("nvert_4hit_1",nvert);
        }

        if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) &&  didCutPass(NHIT1, cutStatus) &&  didCutPass(NHIT3, cutStatus) ) {
          Fill("mmass_4hit_2", immass, weight);
	        Fill("pm_4hit_2", ipm, weight);
          Fill("pm_res_4hit_2", ipm_res);
          Fill("p1_res_4hit_2", ip1_res);
          Fill("p2_res_4hit_2", ip2_res);
          Fill("p3_res_4hit_2", ip3_res);
          Fill("mmass_res_4hit_2", immass_res);
          Fill("nvert_4hit_2",nvert);
        }
       
        if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) &&  didCutPass(NHIT1, cutStatus) &&  didCutPass(NHIT2, cutStatus) ) {
          Fill("mmass_4hit_3", immass, weight);
	        Fill("pm_4hit_3", ipm, weight);
          Fill("pm_res_4hit_3", ipm_res);
          Fill("p1_res_4hit_3", ip1_res);
          Fill("p2_res_4hit_3", ip2_res);
          Fill("p3_res_4hit_3", ip3_res);
          Fill("mmass_res_4hit_3", immass_res);
          Fill("nvert_4hit_3",nvert);
        }

        // Plotting for 2 tracks with 4 hits allowed

        if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(NHIT3, cutStatus) ) {
          Fill("mmass_4hit_12", immass, weight);
	        Fill("pm_4hit_12", ipm, weight);
          Fill("pm_res_4hit_12", ipm_res);
          Fill("p1_res_4hit_12", ip1_res);
          Fill("p2_res_4hit_12", ip2_res);
          Fill("p3_res_4hit_12", ip3_res);
          Fill("mmass_res_4hit_12", immass_res);
          Fill("nvert_4hit_12",nvert);
        }

        if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(NHIT2, cutStatus) ) {
          Fill("mmass_4hit_13", immass, weight);
          Fill("pm_4hit_13", ipm, weight);
          Fill("pm_res_4hit_13", ipm_res);
          Fill("p1_res_4hit_13", ip1_res);
          Fill("p2_res_4hit_13", ip2_res);
          Fill("p3_res_4hit_13", ip3_res);
          Fill("mmass_res_4hit_13", immass_res);
          Fill("nvert_4hit_13",nvert);
	      }

        if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(NHIT1, cutStatus) ) {
          Fill("mmass_4hit_23", immass, weight);
	        Fill("pm_4hit_23", ipm, weight);
          Fill("pm_res_4hit_23", ipm_res);
          Fill("p1_res_4hit_23", ip1_res);
          Fill("p2_res_4hit_23", ip2_res);
          Fill("p3_res_4hit_23", ip3_res);
          Fill("mmass_res_4hit_23", immass_res);
          Fill("nvert_4hit_23",nvert);
        }

        // Plotting for 3 tracks with 4 hits allowed

        if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) ) {
          Fill("mmass_4hit_123", immass, weight);
	        Fill("pm_4hit_123", ipm, weight);
          Fill("pm_res_4hit_123", ipm_res);
          Fill("p1_res_4hit_123", ip1_res);
          Fill("p2_res_4hit_123", ip2_res);
          Fill("p3_res_4hit_123", ip3_res);
          Fill("mmass_res_4hit_123", immass_res);
          Fill("nvert_4hit_123",nvert);
        }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


       ////////////////////////////////////  
       //Plotting parameters with no cuts//
       ////////////////////////////////////

       Fill("pm", ipm, weight);
       Fill("chi2", ichi2, weight);
       Fill("mmass", immass, weight);
       Fill("nrecurl", inrecurl, weight);

       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
       //////////////////////////////////////////////  
       //Plotting parameter resolution with no cuts//
       //////////////////////////////////////////////

       Fill("vx_res", ivx_res, weight);
       Fill("vy_res", ivy_res, weight);
       Fill("vz_res", ivz_res, weight);
       Fill("pm_res", ipm_res, weight);
       Fill("target_dist_res", itarget_dist_res, weight);
       Fill("mmass_res", immass_res, weight);

       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

       /////////////////////////////////////
       //Plotting muon momentum components//
       /////////////////////////////////////

       Fill("px_cms", ipxm, weight);
       Fill("py_cms", ipym, weight);
       Fill("pz_cms", ipzm, weight);
 
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
       ////////////////////////////////////////////
       //Plotting coplanar and acoplanar momentum//
       ////////////////////////////////////////////

       Fill("pcoplanar", ipcoplanar, weight);
       Fill("pacoplanar", ipacoplanar, weight);
   
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

       ///////////////////////////
       //Plotting using TDR cuts//
       ///////////////////////////
  
       if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ) {
       Fill("chi2_TDRcuts",ichi2);
       Fill("nhit1_TDRcuts",inhit1);
       Fill("nhit2_TDRcuts",inhit2);
       Fill("nhit3_TDRcuts",inhit3);

       
       Fill2D("chi2_vs_mass_TDRcuts", immass, ichi2, weight);
       Fill2D("chi2_vs_momentum_TDRcuts", ipm, ichi2, weight);

       Fill2D("mmass_vs_momentum_TDRcuts", immass, ipm, weight);
       Fill("mmass_TDRcuts", immass, weight);
       Fill("pm_TDRcuts", ipm, weight);


       Fill("p1_TDR", ip1, weight);
       Fill("p2_TDR", ip2, weight);
       Fill("p2_TDR", ip3, weight);

       Fill("nsharedhit1_TDRcuts", insharedhit1, weight);
       Fill("nsharedhit2_TDRcuts", insharedhit2, weight);
       Fill("nsharedhit3_TDRcuts", insharedhit3, weight);


       Fill2D("mmass_vs_lambda1_TDRcuts", immass, ilambda1, weight);
       Fill2D("mmass_vs_lambda2_TDRcuts", immass, ilambda2, weight);
       Fill2D("mmass_vs_lambda3_TDRcuts", immass, ilambda3, weight);

       Fill("pT1_TDRcuts", ipT1, weight);
       Fill("pT2_TDRcuts", ipT2, weight);
       Fill("pT3_TDRcuts", ipT3, weight);

       Fill("lambda1_TDRcuts", ilambda1);
       Fill("lambda2_TDRcuts", ilambda2);
       Fill("lambda3_TDRcuts", ilambda3);

       Fill("phi1_TDRcuts", iphi1, weight);
       Fill("phi2_TDRcuts", iphi2, weight);
       Fill("phi3_TDRcuts", iphi3, weight);
      
       Fill("weight_TDRcuts", weight);

       Fill("vx_res_TDRcuts", ivx_res, weight);
       Fill("vy_res_TDRcuts", ivy_res, weight);
       Fill("vz_res_TDRcuts", ivz_res, weight);
       Fill("pm_res_TDRcuts", ipm_res, weight);
       Fill("p1_res_TDRcuts", ip1_res);
       Fill("p2_res_TDRcuts", ip1_res);
       Fill("p3_res_TDRcuts", ip1_res);
       Fill("target_dist_res_TDRcuts", itarget_dist_res, weight);
       Fill("mmass_res_TDRcuts", immass_res, weight);


       Fill("true_p1_TDRcuts",itrue_p1);
       Fill("true_p2_TDRcuts",itrue_p2);
       Fill("true_p3_TDRcuts",itrue_p3);
       }

       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
       /////////////////////////////////////
       //Investigating transverse momentum//
       /////////////////////////////////////

       Fill2D("pT1_vs_chi2", ipT1, ichi2, weight);
       Fill2D("pT2_vs_chi2", ipT2, ichi2, weight);
       Fill2D("pT3_vs_chi2", ipT3, ichi2, weight);

       Fill2D("pT1_vs_pm", ipT1, ipm, weight);
       Fill2D("pT2_vs_pm", ipT2, ipm, weight);
       Fill2D("pT3_vs_pm", ipT3, ipm, weight);


       if ( (ipT1 < 10) ){
        Fill("chi2_pT1<10", ichi2, weight);
        Fill("nsharedhits_particle1_pT1<10",insharedhit1, weight);
        Fill("nsharedseg_particle1_pT1<10", insharedseg1, weight);
        }

        if ( (ipT2 < 10) ){
        Fill("chi2_pT2<10", ichi2, weight);
        Fill("nsharedhits_particle2_pT2<10",insharedhit2, weight);
        Fill("nsharedseg_particle2_pT1<10", insharedseg2, weight);
        }

        if ( (ipT3 < 10) ){
        Fill("chi2_pT3<10", ichi2, weight);
        Fill("nsharedhits_particle3_pT3<10",insharedhit3, weight);
        Fill("nsharedseg_particle3_pT1<10", insharedseg3, weight);
        }


        if ( ichi2 < 100){
        Fill("pT1_chi100",ipT1, weight);
        Fill("pT2_chi100",ipT2, weight);
        Fill("pT3_chi100",ipT3, weight);
        }               

       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
       ///////////////////////////////////////
       //Plotting momentum for each particle//
       ///////////////////////////////////////

       Fill("pz1", ipz1, weight);
       Fill("pz2", ipz2, weight);
       Fill("pz3", ipz3, weight);


       Fill("p1", ip1, weight);
       Fill("p2", ip2, weight);
       Fill("p3", ip3, weight);

       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

       ////////////////////////////////////////
       //Plotting no. hits for each particle //
       ////////////////////////////////////////       

       Fill("nhits1",inhit1, weight);
       Fill("nhits2",inhit2, weight);
       Fill("nhits3",inhit3, weight);


       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       /////////////////////////////////////////
       //Investigating second peak in momentum//
       /////////////////////////////////////////

       if (ipm>=10 && ipm<=90) {
        Fill("pm_secondpeak",ipm, weight);
        Fill("mmass_secondpeak",immass, weight);
        Fill("chi_secondpeak",ichi2, weight);
        Fill("nrecurl_secondpeak",inrecurl, weight);
        Fill("pxm_secondpeak",ipxm, weight);
        Fill("pym_secondpeak",ipym, weight);
        Fill("pzm_secondpeak",ipzm, weight);
       }

       //////////////////////////
       //Cut Flow Investigation//
       //////////////////////////

        Fill("pT1_nocuts", ipT1, weight);
        Fill("pT2_nocuts", ipT2, weight);
        Fill("pT3_nocuts", ipT3, weight);

        if (didCutPass(CHI2, cutStatus)){
          Fill("nvert_chi2_cut", nvert);
        }
        if (didCutPass(MASS, cutStatus)){
          Fill("nvert_mass_cut", nvert);
        }
        if (didCutPass(MEP1, cutStatus)){
          Fill("nvert_mep1_cut", nvert);
        }
        if (didCutPass(MOM, cutStatus)){
          Fill("nvert_mom_cut", nvert);
        }
        if (didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus)){
          Fill("nvert_nhit_cut", nvert);
        }
        if (didCutPass(MASSLOOSE, cutStatus)){
          Fill("nvert_massloose_cut", nvert);
        }
        


        if (didCutPass(CHI2, cutStatus) && didCutPass(MASSLOOSE, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ){
        Fill("mmass_TDR_loose_mass", immass, weight);
        Fill("pm_TDR_loose_mass", ipm, weight);
        }

        if (didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ){
        Fill("mmass_TDR_without_mass_cuts", immass, weight);
        Fill("pm_TDR_without_mass_cuts", ipm, weight);
        }

        if ( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ){
        Fill("mmass_TDR_without_mass_momentum_cuts", immass, weight);
        Fill("pm_TDR_without_mass_momentum_cuts", ipm, weight);
        Fill("lambda1_TDR_without_mass_momentum_cuts", ilambda1, weight);
        Fill("lambda2_TDR_without_mass_momentum_cuts", ilambda2, weight);
        Fill("lambda3_TDR_without_mass_momentum_cuts", ilambda3, weight);
        Fill("phi1_TDR_without_mass_momentum_cuts", iphi1, weight);
        Fill("phi2_TDR_without_mass_momentum_cuts", iphi2, weight);
        Fill("phi3_TDR_without_mass_momentum_cuts", iphi3, weight);
        }

        if (didCutPass(CHI2, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ){
        Fill("mmass_TDR_without_mass_momentum_mep1_cuts", immass, weight);
        Fill("pm_TDR_without_mass_momentum_mep1_cuts", ipm, weight);
        }

        if ( ichi2<100 && didCutPass(MASS, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ){
        Fill("mmass_TDR_cuts_chi2<100", immass, weight);
        Fill("pm_TDR_cuts_chi2<100", ipm, weight);
        }
        
        if ( didCutPass(CHI2, cutStatus) ){
        Fill("mmass_only_chi2_cut", immass, weight);
        Fill("pm_only_chi2_cut", ipm, weight);
        }

        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ){
        Fill("mmass_TDR_without_momentum_cuts", immass, weight);
        Fill("pm_TDR_without_momentum_cuts", ipm, weight);
        }

        if (didCutPass(CHI2, cutStatus) && (immass > 80 || immass < 110) && didCutPass(MEP1, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ){
        Fill("mmass_TDR_mmass>80||<110", immass, weight);
        Fill("pm_TDR_mmass>80||<110", ipm, weight);
        Fill2D("mmass_vs_weight", immass, weight);
        Fill("weight_TDR_mmass>80||<110", weight);
        }

        ////////////////////
        //Plotting weights//
        ////////////////////

        Fill("weights", weight);
        
        ///////////////
        //Manual cuts//
        ///////////////

        if ( ichi2<15 && (immass<110 || immass>103) && ipm<4 && (imep1<5 || imep1>10) && inhit1>4 && inhit2>4 && inhit3>4 ){
        Fill("mmass_manualTDR", immass, weight);
        Fill("pm_manualTDR", immass, weight);
        }

        
        if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && (immass<110 || immass>95) && didCutPass(MOM, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ) {
        Fill("mmass_TDR_to_normalise", immass, weight);
        Fill("mmass_TDR_without_weights", immass);
        
        }

        if( didCutPass(CHI2, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(MASS, cutStatus) && (ipm<20) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ) {
        Fill("pm_TDR_to_normalise1", ipm, weight);
        Fill("pm_TDR_to_normalise2", ipm, weight);
        Fill("pm_TDR_to_normalise3", ipm, weight);
        Fill("pm_TDR_to_normalise4", ipm, weight);
        Fill("pm_TDR_without_weights", ipm);
        }


        if( didCutPass(CHI2, cutStatus) && (immass>95 || immass<110) && (ipm<30) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus) ) {
        Fill2D("mmass_vs_mep1_cuts", immass, imep1);
        }

        Fill2D("mmass_vs_mep1", immass, imep1);

        if( (immass>85 || immass<100)){
          Fill("nvert_faketracks", nvert);
        }

        if( (immass>100 || immass<110)){
          Fill("nvert_realtracks", nvert);
        }

        ///////////////////////
        //nvert investigation//
        ///////////////////////

        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus)){
          Fill("nvert_chi2_mass", nvert);
        }
        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus)){
          Fill("nvert_chi2_mass_mom", nvert);
        }
        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(MEP1, cutStatus)){
          Fill("nvert_chi2_mass_mom_mep1", nvert);
        }
        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus)){
          Fill("nvert_TDR", nvert);
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(NHIT2, cutStatus) && didCutPass(NHIT3, cutStatus)){
          Fill("nvert_TDR_4hit_1", nvert);
        }
        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT3, cutStatus)){
          Fill("nvert_TDR_4hit_2", nvert);
        }
        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(NHIT1, cutStatus) && didCutPass(NHIT2, cutStatus)){
          Fill("nvert_TDR_4hit_3", nvert);
        }
        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(NHIT3, cutStatus)){
          Fill("nvert_TDR_4hit_12", nvert);
        }
        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(NHIT2, cutStatus)){
          Fill("nvert_TDR_4hit_13", nvert);
        }
        if (didCutPass(CHI2, cutStatus) && didCutPass(MASS, cutStatus) && didCutPass(MOM, cutStatus) && didCutPass(MEP1, cutStatus) && didCutPass(NHIT1, cutStatus)){
          Fill("nvert_TDR_4hit_23", nvert);
        }
        


    

      } //end weight if statement
        

     } //end of loop over vertices in frame



   } //end of loop over frames

}
