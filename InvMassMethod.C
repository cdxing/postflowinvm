//gROOT->Reset();
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include "math.h"
#include "string.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
//#include <fstream.h>
//For Residual background >>  fit = 0 : linear, 1 : quadratic  
//Int_t energy = 200;
Int_t energy = 19;
Int_t    fit      = 0;
Int_t    rebin    = 1;// Don't Change  
Int_t opt=0;


//For Background v2 function
//TString  Bgv2     = "Poly3rd"; //"Poly3rd" "Poly2nd" "Poly1st"
TString  Bgv2     = "Poly2nd"; //"Poly3rd" "Poly2nd" "Poly1st"
//TString  Bgv2     = "Poly1st"; //"Poly3rd" "Poly2nd" "Poly1st"
Bool_t   bincount = kFALSE;
//Bool_t   bincount = kTRUE;
Double_t ratio_sig_total[120];
Double_t ratio_bg_total[120];
Double_t ratio_counting[120];
Double_t mass_lowfunc;

//new add
TH1D *hf =new TH1D("hf","",120, 0.98, 1.10);
TH1D *hfb =new TH1D("hfb","",120, 0.98, 1.10);


void InvMassMethod()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetStatColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFillColor(0);
  gStyle->SetLineColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleColor(1); 
  //---------------------------------------------
  //Int_t    centrality[10]       = {0,0,0,0,0,0,0,0,1,1};//0-10
  //Int_t    centrality[10]       = {0,0,0,0,0,1,1,1,0,0};//10-40%
  Int_t    centrality[10]       = {0,0,0,0,1,1,1,0,0,0};//10-40%
  //Int_t    centrality[10]       = {0,1,1,1,1,0,0,0,0,0};//40-80
  //Int_t    centrality[10]       = {0,0,1,1,1,1,1,1,0,0};//10-70%
  //Int_t    centrality[10]       = {0,1,1,1,1,1,1,1,1,1};//0-80%

  //Int_t    centrality[10]       = {0,0,0,0,0,0,1,1,1,1};//0-30
  //Int_t    centrality[10]       = {0,1,1,1,1,1,0,0,0,0};//30-80
  //Int_t    centrality[10]       = {0,1,0,0,0,0,0,0,0,0};//60-80%
  
  Double_t reso = 1.0;//(= 2 for SP Meth else =1)
  Bool_t      CosTerm           = kTRUE;// kTRUE for SP Meth, kFALSE InvMass Meth
  
  Int_t       philow            = 1;
  Int_t       phihi             = 24;
  Bool_t      projection        = kTRUE;
  Int_t       signalturnon      = 0;  
  
  
  //for 7 GeV
  /*  const Int_t nRapBin = 7;
  const Double_t RapBin[nRapBin] = {-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9};
  Double_t    LowNor[nRapBin]    = {1.035, 1.035, 1.035, 1.035, 1.035, 1.035, 1.035};//, 1.035, 1.035};//, 1.035, 1.03, 1.03, 1.035, 1.035, 1.035};
  Double_t    HighNor[nRapBin]   = {1.04, 1.04, 1.05, 1.05, 1.05, 1.05, 1.05};//, 1.05, 1.05};//, 1.05, 1.04, 1.04, 1.05, 1.05, 1.05};
  Int_t    signalbinlow[nRapBin]  = {22, 22, 22, 22, 22, 22, 22};//, 22, 22};//, 22, 22, 22, 22, 22, 22};
  Int_t    signalbinhigh[nRapBin] = {54, 54, 54, 54, 54, 54, 54};//, 54, 54};//, 54, 54, 54, 54, 54, 54};
  Int_t    startfitbin[nRapBin]   = {15, 15, 15, 15, 15, 15, 15};//, 15, 15};//, 15, 15, 15, 15, 15, 15};
  Int_t    fitendbin[nRapBin]     = {109, 109, 109, 109, 109, 109, 109};//, 109, 109};//, 109, 109, 109, 109, 109, 109};
  Double_t lowrange[nRapBin]      = {0.998, 0.998, 0.998, 0.998, 0.998, 0.998, 0.998};//, 0.998, 0.998};//, 0.998, 0.998, 0.998, 0.998, 0.998, 0.998};
  Double_t highrange[nRapBin]     = {1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08};//, 1.08, 1.08};//, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08}; 
 */

    
  //binning Feb 2014
  const Int_t nRapBin = 11;
  //const Double_t RapBin[nRapBin] = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5, 0.7, 0.9};
  const Double_t RapBin[nRapBin] = {0.4, 0.8, 1.2, 1.4, 1.6, 2.0, 2.4, 3.4, 3.6, 3.8, 4.0};
  //const Double_t RapBin[nRapBin] = {-0.7, -0.5, -0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5, 0.7, 0.9};
  Double_t    LowNor[nRapBin]    = {1.035, 1.035, 1.035, 1.035, 1.035, 1.035, 1.035, 1.035, 1.035, 1.035, 1.03};//, 1.03, 1.035, 1.035, 1.035};
  Double_t    HighNor[nRapBin]   = {1.04, 1.04, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.04};//, 1.04, 1.05, 1.05, 1.05};
  Int_t    signalbinlow[nRapBin]  = {22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22};//, 22, 22, 22, 22};
  Int_t    signalbinhigh[nRapBin] = {54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54};//, 54, 54, 54, 54};
  Int_t    startfitbin[nRapBin]   = {15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15};//, 15, 15, 15, 15};
  Int_t    fitendbin[nRapBin]     = {109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109};//, 109, 109, 109, 109};
  Double_t lowrange[nRapBin]      = {0.998, 0.998, 0.998, 0.998, 0.998, 0.998, 0.998, 0.998, 0.998, 0.998, 0.998};//, 0.998, 0.998, 0.998, 0.998};
  Double_t highrange[nRapBin]     = {1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08};//, 1.08, 1.08, 1.08, 1.08}; 
  



  Int_t    xbin[nRapBin]          = {0};
  Double_t NorRangeForCheck[2] = {1.04, 1.06};
  Double_t SigCountRange[2]    = {1.0,1.04};
  //Guess Yield for Fitting
  //Double_t Para1[pTBins]     = {2.4e+05, 5.1e+05, 5.3e+05, 1.5e+05, 1.3e+04, 1.38446e+04};
  Double_t Para1[nRapBin]     = {2.4e+05, 5.1e+05, 5.3e+05, 1.5e+05, 1.3e+04, 1.38446e+04,};

  //---------------------------------------------
  //Int_t CentH[10] = {100,80,70,60,50,40,30,20,10,5};
  //Int_t CentL[10] = {80,70,60,50,40,30,20,10,5,0};
  Int_t CentH[10] = {80,70,60,50,40,30,20,10,5,80};
  Int_t CentL[10] = {70,60,50,40,30,20,10,5,0,0};
  Int_t CentHigh = 5,CentLow = 80;
  for(Int_t i = 0; i < 10; i++)
    {
      if(centrality[i])
	{
	  if(CentLow > CentL[i])
	    CentLow = CentL[i]; 
	  if(CentHigh < CentH[i])
	    CentHigh = CentH[i]; 
	}     
    }


  
  Char_t BgType[56];
  if(fit == 0) sprintf(BgType,"BgPolyOne");
  else if (fit == 1)sprintf(BgType,"BgPolyTwo");
  cout<<"==================================================================="<<endl;
  cout<<"========================== Input parameters ======================="<<endl;
  cout<<"Current Centrality Selection  :  "<<CentLow<<"  "<<CentHigh<<endl;
  cout<<"Background fitting function   :  "<<BgType<<endl;
  cout<<"Pt Rebin                      :  "<<rebin<<endl;
  cout<<"Bin counting                  :  "<<bincount<<endl;
  cout<<"Background V2 fitting function:  "<<Bgv2.Data()<<endl;
  cout<<"==================================================================="<<endl;
  
  TProfile *ProMassV2P[nRapBin-1];
  TH2D *histinv[nRapBin-1];
  Int_t    MassBin[nRapBin-1];
  Int_t    ndf[nRapBin-1];
  Double_t MassLow[nRapBin-1];
  Double_t MassHi[nRapBin-1];
  Double_t NorFact[nRapBin-1];
  Double_t v2obs[nRapBin-1];
  Double_t v2obserr[nRapBin-1];
  Double_t chisquare[nRapBin-1];
  Double_t chisquarendf[nRapBin-1];
  Double_t RawYield[nRapBin-1];
  Double_t ErRawYield[nRapBin-1];
  Double_t Mass[nRapBin-1];
  Double_t ErMass[nRapBin-1];
  Double_t Width[nRapBin-1];
  Double_t ErWidth[nRapBin-1];
  Double_t pT[nRapBin-1];
  Double_t ErpT[nRapBin-1];
  Double_t dNdpTdy[nRapBin-1];
  Double_t ErdNdpTdy[nRapBin-1];
  Double_t BinCount[nRapBin-1];

  Double_t rapi[nRapBin-1];
  Double_t err_rapi[nRapBin-1];
  TH2D *dummy_prof[nRapBin-1];

  Double_t significance[nRapBin-1];
  Double_t significance_err[nRapBin-1]={0,};
  Double_t SigByBg[nRapBin-1];
  Double_t SigByBg_err[nRapBin-1]={0,};
  
  //TFile *f1 = new TFile(Form("data/output_%dGeV_run11_mb6.root", energy));
  TFile *f1 = new TFile(Form("./data/Yields_SE_%dGeV.root", energy));
  TFile *f2 = new TFile(Form("./data/Yields_ME_%dGeV_53.root", energy));

  TH3F *hSigCBkg0  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_0_Phi_SE");    
  TH3F *hCBkg0     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_0_Phi_ME"); 
  TH3F *hSigCBkg1  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_1_Phi_SE");     
  TH3F *hCBkg1     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_1_Phi_ME"); 
  TH3F *hSigCBkg2  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_2_Phi_SE");
  TH3F *hCBkg2     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_2_Phi_ME"); 
  TH3F *hSigCBkg3  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_3_Phi_SE");
  TH3F *hCBkg3     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_3_Phi_ME");
  TH3F *hSigCBkg4  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_4_Phi_SE");
  TH3F *hCBkg4     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_4_Phi_ME"); 
  TH3F *hSigCBkg5  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_5_Phi_SE");
  TH3F *hCBkg5     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_5_Phi_ME");
  TH3F *hSigCBkg6  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_6_Phi_SE");
  TH3F *hCBkg6     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_6_Phi_ME");
  TH3F *hSigCBkg7  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_7_Phi_SE");
  TH3F *hCBkg7     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_7_Phi_ME");
  TH3F *hSigCBkg8  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_8_Phi_SE");
  TH3F *hCBkg8     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_8_Phi_ME");
  //TH3F *hSigCBkg9  = (TH3F *) f1->Get("NumInvMassvsPtPhi_Cen_9_Phi_SE");
  //TH3F *hCBkg9     = (TH3F *) f2->Get("NumInvMassvsPtPhi_Cen_9_Phi_ME"); 
 
  
  TH3F *hNum = new TH3F("hNum","hNum",
			hSigCBkg8->GetNbinsX(),hSigCBkg8->GetXaxis()->GetXmin(),hSigCBkg8->GetXaxis()->GetXmax(),
			hSigCBkg8->GetNbinsY(),hSigCBkg8->GetYaxis()->GetXmin(),hSigCBkg8->GetYaxis()->GetXmax(),
			hSigCBkg8->GetNbinsZ(),hSigCBkg8->GetZaxis()->GetXmin(),hSigCBkg8->GetZaxis()->GetXmax());
  TH3F *hDen = new TH3F("hDen","hDen",
			hCBkg8->GetNbinsX(),hCBkg8->GetXaxis()->GetXmin(),hCBkg8->GetXaxis()->GetXmax(),
			hCBkg8->GetNbinsY(),hCBkg8->GetYaxis()->GetXmin(),hCBkg8->GetYaxis()->GetXmax(),
			hCBkg8->GetNbinsZ(),hCBkg8->GetZaxis()->GetXmin(),hCBkg8->GetZaxis()->GetXmax());
			//hSigCBkg9->GetNbinsX(),hSigCBkg9->GetXaxis()->GetXmin(),hSigCBkg9->GetXaxis()->GetXmax(),
			//hSigCBkg9->GetNbinsY(),hSigCBkg9->GetYaxis()->GetXmin(),hSigCBkg9->GetYaxis()->GetXmax(),
			//hSigCBkg9->GetNbinsZ(),hSigCBkg9->GetZaxis()->GetXmin(),hSigCBkg9->GetZaxis()->GetXmax());
  //TH3F *hDen = new TH3F("hDen","hDen",
			//hCBkg9->GetNbinsX(),hCBkg9->GetXaxis()->GetXmin(),hCBkg9->GetXaxis()->GetXmax(),
			//hCBkg9->GetNbinsY(),hCBkg9->GetYaxis()->GetXmin(),hCBkg9->GetYaxis()->GetXmax(),
			//hCBkg9->GetNbinsZ(),hCBkg9->GetZaxis()->GetXmin(),hCBkg9->GetZaxis()->GetXmax());
  /*
    cout<<hSigCBkg8->GetNbinsX()<<"  "<<
    hSigCBkg8->GetXaxis()->GetXmin()<<"  "<<
    hSigCBkg8->GetXaxis()->GetXmax()<<"  "<<
    hSigCBkg8->GetNbinsY()<<"  "<<
    hSigCBkg8->GetYaxis()->GetXmin()<<"  "<<
    hSigCBkg8->GetYaxis()->GetXmax()<<"  "<<
    hSigCBkg8->GetNbinsZ()<<"  "<<
    hSigCBkg8->GetZaxis()->GetXmin()<<"  "<<
    hSigCBkg8->GetZaxis()->GetXmax()<<"  "<<endl;
    hSigCBkg9->Draw();
  */
  if(centrality[0]){
    hNum->Add(hSigCBkg0); 
    hDen->Add(hCBkg0);
    //cout<<"Centrality  : 80 - 100"<<endl;
    cout<<"Centrality  : 70 - 80"<<endl;
  } 
  if(centrality[1]){
    hNum->Add(hSigCBkg1); 
    hDen->Add(hCBkg1); 
    //cout<<"Centrality  : 70 - 80"<<endl;
    cout<<"Centrality  : 60 - 70"<<endl;
  }
  if(centrality[2]){
    hNum->Add(hSigCBkg2); 
    hDen->Add(hCBkg2); 
    //cout<<"Centrality  : 60 - 70"<<endl;
    cout<<"Centrality  : 50 - 60"<<endl;
  }
  if(centrality[3]){
    hNum->Add(hSigCBkg3); 
    hDen->Add(hCBkg3);
    //cout<<"Centrality  : 50 - 60"<<endl;
    cout<<"Centrality  : 40 - 50"<<endl;
  }
  if(centrality[4]){
    hNum->Add(hSigCBkg4); 
    hDen->Add(hCBkg4); 
    //cout<<"Centrality  : 40 - 50"<<endl;
    cout<<"Centrality  : 30 - 40"<<endl;
  }
  if(centrality[5]){
    hNum->Add(hSigCBkg5); 
    hDen->Add(hCBkg5); 
    //cout<<"Centrality  : 30 - 40"<<endl;
    cout<<"Centrality  : 20 - 30"<<endl;
  }
  if(centrality[6]){
    hNum->Add(hSigCBkg6); 
    hDen->Add(hCBkg6); 
    //cout<<"Centrality  : 20 - 30"<<endl;
    cout<<"Centrality  : 10 - 20"<<endl;
  }
  if(centrality[7]){
    hNum->Add(hSigCBkg7); 
    hDen->Add(hCBkg7);
    //cout<<"Centrality  : 10 - 20"<<endl;
    cout<<"Centrality  : 5 - 10"<<endl;
  }
  if(centrality[8]){
    hNum->Add(hSigCBkg8); 
    hDen->Add(hCBkg8); 
    //cout<<"Centrality  : 5 - 10"<<endl;
    cout<<"Centrality  : 0 - 5"<<endl;
  }
  /*if(centrality[9]){
    hNum->Add(hSigCBkg9); 
    hDen->Add(hCBkg9); 
    cout<<"Centrality  : 0 - 5"<<endl;
  }*/


  TFile *tfout = new TFile(Form("out_phiv1_Bg_%dGeV_cen%d%d.root", energy, CentLow, CentHigh),"recreate");
  tfout->cd();



  Double_t    MeanLow           = -1.0;
  Double_t    MeanHi            = 1.0;
    for(Int_t ip = 0; ip < nRapBin-1; ip++)
    {
      BinCount[ip] = 0;
      xbin[ip] = fitendbin[ip] - startfitbin[ip] + 1;
      MassBin[ip] = xbin[ip];
      MassLow[ip] = 0.98 + (startfitbin[ip]-1)*0.001;;
      MassHi[ip]  = 0.98 + (startfitbin[ip]-1.0 + xbin[ip])*0.001;

      cout<<xbin[ip]<<"  "<<MassBin[ip]<<"  "<<MassLow[ip]<<"  "<<MassHi[ip]<<endl;
      //ProMassV2P[ip]= new TProfile(Form("MassV2P%d",ip+1),Form("MassV2P%d",ip+1),MassBin[ip],MassLow[ip],MassHi[ip],MeanLow,MeanHi);
      //ProMassV2P[ip]->Sumw2();
      
    }
  
  if(projection)
    {
      //TH3F *massv2 = new TH3F("massv2","massv2",20, -1.0, 1.0, 200, -1.0, 1.0, 120, 0.98, 1.10);
      TH3F *massv2 = new TH3F("massv2","massv2",20, -1.0, 1.0, 1200, -12.0, 12.0, 120, 0.98, 1.10);
      if(energy==62 || energy==200) massv2 = new TH3F("massv2","massv2",20, -1.0, 1.0, 1600, -16.0, 16.0, 120, 0.98, 1.10);
      
      if(CosTerm){
	TH3F *m0 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_0_Phi_SE");
	TH3F *m1 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_1_Phi_SE");
	TH3F *m2 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_2_Phi_SE");
	TH3F *m3 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_3_Phi_SE");
	TH3F *m4 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_4_Phi_SE");
	TH3F *m5 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_5_Phi_SE");
	TH3F *m6 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_6_Phi_SE");
	TH3F *m7 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_7_Phi_SE");
	TH3F *m8 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_8_Phi_SE");
	TH3F *m9 = (TH3F *)f1->Get("ebyeInvMassv2_Cen_9_Phi_SE");
      }
      else{
	//with out resolution correction
	TH3F *m0 = (TH3F *)f1->Get("InvMassv1_Cen_0_Phi_SE");
	TH3F *m1 = (TH3F *)f1->Get("InvMassv1_Cen_1_Phi_SE");
	TH3F *m2 = (TH3F *)f1->Get("InvMassv1_Cen_2_Phi_SE");
	TH3F *m3 = (TH3F *)f1->Get("InvMassv1_Cen_3_Phi_SE");
	TH3F *m4 = (TH3F *)f1->Get("InvMassv1_Cen_4_Phi_SE");
	TH3F *m5 = (TH3F *)f1->Get("InvMassv1_Cen_5_Phi_SE");
	TH3F *m6 = (TH3F *)f1->Get("InvMassv1_Cen_6_Phi_SE");
	TH3F *m7 = (TH3F *)f1->Get("InvMassv1_Cen_7_Phi_SE");
	TH3F *m8 = (TH3F *)f1->Get("InvMassv1_Cen_8_Phi_SE");
	TH3F *m9 = (TH3F *)f1->Get("InvMassv1_Cen_9_Phi_SE");
      }
      
      if(centrality[0]){
	massv2->Add(m0);
	//cout<<"Centrality  : 80 - 100"<<endl;
	cout<<"Centrality  : 70 - 80"<<endl;
      } 
      if(centrality[1]){
	massv2->Add(m1);
	//cout<<"Centrality  : 70 - 80"<<endl;
	cout<<"Centrality  : 60 - 70"<<endl;
      } 
      if(centrality[2]){
	massv2->Add(m2); 
	//cout<<"Centrality  : 60 - 70"<<endl;
	cout<<"Centrality  : 50 - 60"<<endl;
      } 
      if(centrality[3]){
	massv2->Add(m3);
	//cout<<"Centrality  : 50 - 60"<<endl;
	cout<<"Centrality  : 40 - 50"<<endl;
      } 
      if(centrality[4]){
	massv2->Add(m4);
	//cout<<"Centrality  : 40 - 50"<<endl;
	cout<<"Centrality  : 30 - 40"<<endl;
      } 
      if(centrality[5]){
	massv2->Add(m5);
	//cout<<"Centrality  : 30 - 40"<<endl;
	cout<<"Centrality  : 20 - 30"<<endl;
      } 
      if(centrality[6]){
	massv2->Add(m6);
	//cout<<"Centrality  : 20 - 30"<<endl;
	cout<<"Centrality  : 10 - 20"<<endl;
      } 
      if(centrality[7]){
	massv2->Add(m7);
	//cout<<"Centrality  : 10 - 20"<<endl;
	cout<<"Centrality  : 5 - 10"<<endl;
      } 
      if(centrality[8]){
	massv2->Add(m8);
	//cout<<"Centrality  : 5 - 10"<<endl;
	cout<<"Centrality  : 0 - 5"<<endl;
      } 
      /*if(centrality[9]){
	massv2->Add(m9);
	cout<<"Centrality  : 0 - 5"<<endl;
      } */
      
      //massv2->Draw();return;
      m0->Reset();
      m1->Reset();
      m2->Reset();
      m3->Reset();
      m4->Reset();
      m5->Reset();
      m6->Reset();
      m7->Reset();
      m8->Reset();
      //m9->Reset();

      /*         
      for(Int_t ip = 0; ip < nRapBin-1; ip++)// pT bins
	{ 
	  //for(Int_t i = 1; i<= 200; i++)//Cos[2(phi-Psi)] Bins
	  for(Int_t i = 1; i<= 1000; i++)//Cos[2(phi-Psi)] Bins
	    { 
	      for(Int_t k = 0; k < xbin[ip]; k++)//invariant mass
		{
		  Int_t    InvMBin  = startfitbin[ip] + k;
		  Double_t MassCent = InvMBin*0.001 + 0.9795;
		  Double_t CosCent  = i*0.01 + (-1.005);// Ajay
		  Int_t RapLow = hNum->GetXaxis()->FindBin(RapBin[ip]);
		  Int_t RapMax = hNum->GetXaxis()->FindBin(RapBin[ip+1]);
		  
		  for(Int_t j = RapLow; j <= RapMax-1; j++) //different Pt Bins
		    {
		      Double_t weight = massv2->GetBinContent(j,i,InvMBin);
		      ProMassV2P[ip]->Fill(MassCent, CosCent, weight); 
		      //cout<<MassCent<<"  "<<CosCent<<"  "<<weight<<endl;
		      //cout<< RapLow << "  "<< RapMax-1 << " "<< MassCent <<"  "<< CosCent<<"  "<<weight<<endl;
		    }//j loop	
		}// k loop
	    }//i loop
	  ProMassV2P[ip]->Sumw2(); 
	  //TCanvas *ajay = new TCanvas(Form("ajay%d",ip),"",10,10,700,700);
	  //ajay->cd(ip);
	  //ProMassV2P[ip]->Draw();
	}//ipt loop
      TFile f2("PbPbMVsV2.root","RECREATE");
      for(Int_t ip = 0; ip < nRapBin-1; ip++)
	{
	  ProMassV2P[ip]->Write();
	}
      f2.Close();

*/

    }//if projection
 
  //return; //subhash checked


  
  //*****************start pt loop*********************//
  TCanvas *cinv[nRapBin-1];
  TCanvas *canfit2[nRapBin-1];
  TH1D *hq[nRapBin-1];
  TH1D *hbq[nRapBin-1];
  TH1D *hDiff[nRapBin-1];
  TH1D    *hRat[nRapBin-1];
  TFitResultPtr r;
  TF1 *fBW[nRapBin-1];
  TF1 *fBg[nRapBin-1];
  
  for(Int_t ip = 0; ip < nRapBin-1; ip++)
  //for(Int_t ip = 0; ip < 1; ip++)
    {
      Int_t raplow = hNum->GetXaxis()->FindBin(RapBin[ip]);
      Int_t raphi = hNum->GetXaxis()->FindBin(RapBin[ip+1]);
      //Int_t    ptlow     = LowPtBin[ip];
      //Int_t    pthi      = HighPtBin[ip];
      Double_t fitlow    = lowrange[ip];
      Double_t fithi     = highrange[ip];
      Double_t lownorm   = LowNor[ip];
      Double_t hinorm    = HighNor[ip];
      Double_t unitmass  = 0.001;
      Double_t scalenorm = 1.0;
      Double_t peakposlow = 1.01;
      Double_t peakposhig = 1.025;
      Double_t ScaleBG    = 0;
      Int_t    Iter       = 5;  
      cout<<"Projecting X and Y bins  : "<<raplow<<"  "<<raphi<<"  "<<philow<<"  "<<phihi<<endl;
      TString hNameS = "hsig";
      hNameS += ip;
      TString hNameB = "hbkg";
      hNameB += ip;
      hq[ip]  = hNum->ProjectionZ(hNameS.Data(), raplow, raphi, philow, phihi,"e"); 
      hq[ip]->SetXTitle("m_{inv} (GeV/c^{2})");
      //hq[ip]->SetTitle("SIG and scaled BG counts vs. m_{inv}");
      hbq[ip] = hDen->ProjectionZ(hNameB.Data(), raplow, raphi, philow, phihi,"e");
      hbq[ip]->SetXTitle("m_{inv} (GeV/c^{2})");
      hbq[ip]->SetLineColor(2);

      hRat[ip] = new TH1D(Form("hNorBg*d",ip),Form("hNorBg*d",ip),120,0.98,1.1);
      hRat[ip] = (TH1D *)hbq[ip]->Clone(hNameB.Data());
      Int_t nbinsx = hq[ip]->GetNbinsX();
      cout << "Nbinsx = " << nbinsx << endl;
      hq[ip]->Rebin(rebin);
      hbq[ip]->Rebin(rebin);
      hRat[ip]->Rebin(rebin);
      cinv[ip] = new TCanvas(Form("cinv%d",ip),Form("cinv%d",ip),10,10,1200,600); 
      cinv[ip]->Divide(2,1);
      cinv[ip]->cd(2);
      if(Iter)
	{
	  hDiff[ip] = MakeDiff(Iter,hq[ip],hbq[ip],peakposlow,peakposhig,ScaleBG,r,fitlow,fithi,Para1[ip]);
	  //hbg[ip]->Scale(ScaleBG);
	  hRat[ip]->Scale(ScaleBG);
	}
      //hDiff[ip]->SetTitle(Form("%3.2f - rap_cent - %3.2f",hNum->GetXaxis()->GetBinCenter(raplow),hNum->GetXaxis()->GetBinCenter(raphi)));
      hDiff[ip]->SetTitle(Form("%3.2f < y < %3.2f", RapBin[ip], RapBin[ip+1]));
      Double_t sigInt1 = hq[ip]->Integral(hq[ip]->GetXaxis()->FindBin(NorRangeForCheck[0]), hq[ip]->GetXaxis()->FindBin(NorRangeForCheck[1]));
      Double_t bgInt1  = hbq[ip]->Integral(hbq[ip]->GetXaxis()->FindBin(NorRangeForCheck[0]), hbq[ip]->GetXaxis()->FindBin(NorRangeForCheck[1]));
      NorFact[ip] = ScaleBG;
      //checks the normalization constant;
      cout<< "SIGNAL ====>  " << sigInt1 <<"   Background =====>  " << bgInt1 << endl;
      cout<< "Signal/Signal+Bkg                   ======= > " << (sigInt1/bgInt1) << endl;
      cout<< "Normalization from iteration method ======= > "  << ScaleBG <<endl;

      
      //fetching the fit parameters
      TF1      *func  = hDiff[ip]->GetFunction(Form("hsig%dfunc%d",ip,Iter-1));
      Double_t *par   = func->GetParameters();
      Double_t *Erpar = func->GetParErrors();
      //Breit-Wigner function alone
      fBW[ip] = new TF1(Form("BW%d",ip),Form("(1.0*0.001*%d/(2*3.14159))*(([0]*[1])/((x-[2])**2 +([1]/2)**2))",rebin), lowrange[ip], highrange[ip]);
      fBW[ip]->SetParameters(&par[0]);
      fBW[ip]->SetParErrors(&Erpar[0]);
      //Background function alone
      if(fit == 0)
	{
	  fBg[ip]   = new TF1(Form("Bg%d",ip),"[0]+[1]*x",lowrange[ip], highrange[ip]);
	}
      else if (fit == 1)
	{
	  fBg[ip]   = new TF1(Form("Bg%d",ip),"[0]+[1]*x+[2]*x*x",lowrange[ip], highrange[ip]);
	}
      fBg[ip]->SetParameters(&par[3]);
      fBg[ip]->SetParErrors(&Erpar[3]);
      
      hDiff[ip]->Draw();
      hDiff[ip]->SetMarkerStyle(20);
      hDiff[ip]->SetMarkerColor(4);
      hDiff[ip]->SetMarkerSize(1.0);
      
      fBg[ip]->SetLineColor(6);
      fBg[ip]->SetLineStyle(2);
      fBg[ip]->Draw("lsame");
      
      fBW[ip]->SetLineColor(2);
      fBW[ip]->SetLineStyle(4);
      fBW[ip]->Draw("lsame");
      TLegend * lsig1 = DrawLegend();
      lsig1->AddEntry(func,"Total Function","l");
      lsig1->AddEntry(fBW[ip],"BW Function only","l");
      //lsig1->AddEntry(fBg[ip],"Bg Function only","l");
      lsig1->Draw();
      
      cinv[ip]->cd(1);
      hq[ip]->SetTitle(Form("%3.2f < y < %3.2f", RapBin[ip], RapBin[ip+1]));
      hq[ip]->DrawCopy();
      hRat[ip]->DrawCopy("same");
      TLegend * lsig = DrawLegend();
      lsig->AddEntry(hq[ip],"Total","l");
      lsig->AddEntry(hRat[ip],"Scaled Bg","l");
      lsig->Draw();

      cinv[ip]->Modified();
      cinv[ip]->Update();
      //cinv[ip]->SaveAs(Form("Plots/InvMassDistributionMethodOne%3.2fPtCent%3.2f_%dYbin%dCent_%d_%d%s.gif",ip,ip+1,philow,phihi,CentLow,CentHigh,BgType));
      cinv[ip]->SaveAs(Form("Plots/SignalBkg_%dRap%d_%dDphi%d_Cent_%d%d_%s.gif",ip,ip+1,philow,phihi,CentLow,CentHigh,BgType));

      cinv[ip]->Write();

      

      //return;//subhash
      
      TMatrixDSym mat = r->GetCovarianceMatrix();
      //mat.Print();
      TMatrixDSym mat1;
      mat.GetSub(0,2,0,2,mat1);
      //mat1.Print();
      Double_t * b    = mat1.GetMatrixArray();
      if(fit == 0){
        TMatrixDSym mat2;
        mat.GetSub(3,4,3,4,mat2);
        //mat2.Print();
        Double_t * a    = mat2.GetMatrixArray();
      }
      else if(fit == 1){
	TMatrixDSym mat2;
	mat.GetSub(3,5,3,5,mat2);
	//mat2.Print();
	Double_t * a    = mat2.GetMatrixArray();
      }


      //signal to background lowrange[ip], highrange[ip]
      Double_t sigInt2 = hq[ip]->Integral(hq[ip]->GetXaxis()->FindBin(lowrange[ip]), hq[ip]->GetXaxis()->FindBin(highrange[ip]));
      Double_t bgInt2  = hbq[ip]->Integral(hbq[ip]->GetXaxis()->FindBin(lowrange[ip]), hbq[ip]->GetXaxis()->FindBin(highrange[ip]));
      SigByBg[ip] = sigInt2/bgInt2;//(S/S+B)
      significance[ip] = sigInt2/sqrt(bgInt2);

      /*    
      //signal to background and significance
      Double_t SG = fBW[ip]->Integral(SigCountRange[0],SigCountRange[1])/(0.001*rebin);
      Double_t BG = hq[ip]->Integral(hq[ip]->GetXaxis()->FindBin(SigCountRange[0]), hq[ip]->GetXaxis()->FindBin(SigCountRange[1]));;
      cout << " SG = " << SG << " BG =   "<<  BG << "  significance = "<< SG/TMath::Sqrt(BG) << "  SG/BG  = " << SG/BG<< endl;
      SigByBg[ip] = SG/BG;
      significance[ip] = SG/sqrt(BG);
      */
      
      rapi[ip]   = (RapBin[ip] + RapBin[ip+1]) / 2.;
      err_rapi[ip]   = TMath::Abs((RapBin[ip+1] - RapBin[ip])/2.);
      //Double_t dpT   = (0.1*((pthi - ptlow) + 1));
      //Double_t dy    = 2.0;
      //pT[ip]         = 0.5*(pthi + (ptlow-1))*0.1;
      //ErpT[ip]       = dpT*0.5;
      RawYield[ip]   = fBW[ip]->Integral(SigCountRange[0],SigCountRange[1])/(0.001*rebin);
      ErRawYield[ip] = fBW[ip]->IntegralError(SigCountRange[0],SigCountRange[1],&par[0],b)/(0.001*rebin);
      Width[ip]      = func->GetParameter(1);
      ErWidth[ip]    = func->GetParError(1);
      Mass[ip]       = func->GetParameter(2);
      ErMass[ip]     = func->GetParError(2);
      //dNdpTdy[ip]    = RawYield[ip]/(Float_t)(TMath::TwoPi()*pT[ip]*dpT*dy);
      //ErdNdpTdy[ip]  = ErRawYield[ip]/(Float_t)(TMath::TwoPi()*pT[ip]*dpT*dy);    
      cout<<par[0]<<"  "<<par[1]<<"  "<<par[2]<<"  "<<par[3]<<"  "<<par[4]<<endl;	 

      //Bin Counting method
      Double_t sigcount = 0;
      Double_t bkgcount = 0;
      for(Int_t bin = hDiff[ip]->GetXaxis()->FindBin(SigCountRange[0]); bin <= hDiff[ip]->GetXaxis()->FindBin(SigCountRange[1]); bin++){
	sigcount = hDiff[ip]->GetBinContent(bin);
	if(fit == 0)
	  {
	    bkgcount = par[3]+par[4]*hDiff[ip]->GetBinCenter(bin);
	  }
	else if(fit == 1)
	  {
	    bkgcount = par[3]+par[4]*hDiff[ip]->GetBinCenter(bin)+par[5]*hDiff[ip]->GetBinCenter(bin)*hDiff[ip]->GetBinCenter(bin);
	  }
	if((sigcount - bkgcount) > 0)
	  BinCount[ip] += TMath::Abs(sigcount - bkgcount);
      }


      Double_t *signal_fit = new Double_t[xbin[ip]];
      Double_t *total      = new Double_t[xbin[ip]];
      Double_t xx; 
      Double_t signal_count   = 0;
      Double_t signalrebin    = 0;


      for(Int_t t = 0; t < xbin[ip]; t++)
	{
	  total[t]      = 0;
	  signal_fit[t] = 0; 
	  Double_t bkgc = 0;
	  Int_t massbin = t +  startfitbin[ip];
	  xx = 0.9795+0.001*massbin;//xx=1.0+0.001+0.002*t;
	  //cout<<massbin<<"   masscenter  "<<xx<<endl;
	  if(signalturnon == 1 && massbin >= signalbinlow[ip] && massbin <= signalbinhigh[ip])
	    {
	      signalrebin  = 1.0*0.001*par[0]*par[1]/(2*3.1415926)/(pow(xx-par[2],2)+par[1]*par[1]/4.);
	      if(fit == 0)
		signal_count = hDiff[ip]->GetBinContent(massbin) - (par[3]+par[4]*xx);
	      if(fit == 1)
		signal_count = hDiff[ip]->GetBinContent(massbin) - (par[3]+par[4]*xx+par[5]*xx*xx);

	      total[t]     = hq[ip]->GetBinContent(massbin);  // from 0.995 to 1.06 bin-width 1MeV	
	      bkgc         = hbq[ip]->GetBinContent(massbin);
	    }
	  
	  if(signalturnon != 1 )
	    {
	      signalrebin  = 1.0*0.001*par[0]*par[1]/(2*3.1415926)/(pow(xx-par[2],2)+par[1]*par[1]/4.);
	      if(fit == 0)
	      signal_count = hDiff[ip]->GetBinContent(massbin)- (par[3]+par[4]*xx);
	      if(fit == 1)
		signal_count = hDiff[ip]->GetBinContent(massbin) - (par[3]+par[4]*xx+par[5]*xx*xx);
	      total[t]     = hq[ip]->GetBinContent(massbin); 
	      bkgc         = hbq[ip]->GetBinContent(massbin);
	    }

	  if(bincount)
	    signal_fit[t]   = signal_count;
	  else
	    signal_fit[t]   = signalrebin;
	  
	  	  	  
	  if(fabs(total[t]) < 1.0e-6) ratio_sig_total[t] = 1.0e-6;
	  else  ratio_sig_total[t] = signal_fit[t]/total[t];
	  
	  if(fabs(total[t]) < 1.0e-6) ratio_bg_total[t] = 1.0e-6;
	  else  ratio_bg_total[t] = bkgc/total[t];	   
	  
	  if(fabs(total[t])<1.0e-6)ratio_counting[t]=1.0e-6;
	  else ratio_counting[t] = fabs(signal_count/total[t]); 
	}
	  
      
      

      /*  
	  TCanvas *cant= new TCanvas(Form("cant%d",ip),"cant",10,10,1000,800);
	  cant->cd();
	  hq->Draw();
	  hbqclone->Draw("SAME");
	  cant->Update();
      */
      delete signal_fit;
      delete total;


      //new add
      hf = hDiff[ip];
      hfb = hDiff[ip];
      hf->Divide(hRat[ip]);
      //===================== Fit inv mass v2 ===============================
      //TCanvas *canfit1= new TCanvas(Form("canfit%d",ip),"canfit1",10,10,800,700);
      //canfit1->cd();
      //canfit1->cd(ip+1);
      canfit2[ip]= new TCanvas(Form("canfit%d",ip),Form("canfit%d",ip),10,10,800,700);
      //canfit2[ip]= new TCanvas(Form("canfit%d",ip),Form("canfit%d",ip),10,10,1200,600);
      //canfit2[ip]->Divide(2,1);
      //canfit2[ip]->cd(1);
      Double_t fitv2low = 0.98 + (startfitbin[ip]-1)*0.001;
      Double_t fitv2hi  = 0.98 + (startfitbin[ip]-1.0 + xbin[ip])*0.001;  
      
      mass_lowfunc = MassLow[ip];
      cout<<"fitv2hi   "<<fitv2hi<<endl;

      dummy_prof[ip] = new TH2D(Form("dummy%d",ip+1),"", 1200, fitv2low, fitv2hi, 200, -0.1, 0.1);
      if(energy==62 || energy==200) dummy_prof[ip] = new TH2D(Form("dummy%d",ip+1),"", 1600, fitv2low, fitv2hi, 200, -0.1, 0.1);
      //dummy_prof[ip] = new TH2D(Form("dummy%d",ip+1),"", 200, fitv2low, fitv2hi, 200, -0.1, 0.1);
      //if(energy==62) dummy_prof[ip] = new TH2D(Form("dummy%d",ip+1),"", 1600, fitv2low, fitv2hi, 200, -0.1, 0.1);

      TF1 *fitfun;
      /*    if(Bgv2.CompareTo("Poly3rd") == 0)
	fitfun = new TF1("fitfun",fitfun_raw,fitv2low,fitv2hi,5);
      if(Bgv2.CompareTo("Poly2nd") == 0)
	fitfun = new TF1("fitfun",fitfun_raw,fitv2low,fitv2hi,4);
      if(Bgv2.CompareTo("Poly1st") == 0)
	fitfun = new TF1("fitfun",fitfun_raw,fitv2low,fitv2hi,3);
     */
      //new add
      if(Bgv2.CompareTo("Poly3rd") == 0)
	fitfun = new TF1("fitfun", Invm_func_tot,fitv2low,fitv2hi,5);
      if(Bgv2.CompareTo("Poly2nd") == 0)
	fitfun = new TF1("fitfun", Invm_func_tot,fitv2low,fitv2hi,4);
      if(Bgv2.CompareTo("Poly1st") == 0)
	fitfun = new TF1("fitfun", Invm_func_tot,fitv2low,fitv2hi,3);
      

      //added subhash
      massv2->GetXaxis()->SetRange(raplow,raphi);
      histinv[ip] = (TH2D) massv2->Project3D("yz");
      ProMassV2P[ip] = (TProfile*) histinv[ip]->ProfileX();

 
      ProMassV2P[ip]->SetMarkerStyle(20);
      ProMassV2P[ip]->SetMarkerColor(2);

      fitfun->SetParName(0,"v_{1}^{obs}"); 
      fitfun->SetParName(1,"p0");
      fitfun->SetParName(2,"p1");
      fitfun->SetParName(3,"p2");
      fitfun->SetParName(4,"p3");
      
      fitfun->SetParameter(0,0.005);
      fitfun->SetParameter(1,1.39);
      fitfun->SetParameter(2,-3.8);
      fitfun->SetParameter(3,3.6);
      fitfun->SetParameter(4,-1.199);
      fitfun->SetNpx(10000);
      fitfun->SetLineColor(1);
      dummy_prof[ip]->SetTitle(Form("%3.2f < y < %3.2f", RapBin[ip], RapBin[ip+1]));
      if(opt==1) 
	ProMassV2P[ip]->Fit("fitfun","ERMI");
      else ProMassV2P[ip]->Fit("fitfun","QER");
      dummy_prof[ip]->Draw("");
      ProMassV2P[ip]->DrawCopy("same");
      //fitfun->DrawCopy("same");
      dummy_prof[ip]->GetYaxis()->SetRangeUser(-0.07, 0.07);
      //To Draw sig and bkg components
      TF1 *fitfunsig;
      fitfunsig = new TF1("fitfunsig",fitfun_sig, fitv2low, fitv2hi,1);
      fitfunsig->SetParameter(0, fitfun->GetParameter(0));
      fitfunsig->SetLineColor(4);
      fitfunsig->SetLineStyle(2);
      fitfunsig->DrawCopy("same");
      
      TF1 *fitfunbkg;
      if(Bgv2.CompareTo("Poly3rd") == 0)
	{fitfunbkg = new TF1("fitfunbkg",fitfun_bkg, fitv2low, fitv2hi,4);
	  fitfunbkg->SetParameters(fitfun->GetParameter(1), fitfun->GetParameter(2), fitfun->GetParameter(3), fitfun->GetParameter(4));}
      if(Bgv2.CompareTo("Poly2nd") == 0)
	{fitfunbkg = new TF1("fitfunbkg",fitfun_bkg, fitv2low, fitv2hi,3);
	  fitfunbkg->SetParameters(fitfun->GetParameter(1), fitfun->GetParameter(2), fitfun->GetParameter(3));}
      if(Bgv2.CompareTo("Poly1st") == 0)
	{fitfunbkg = new TF1("fitfunbkg",fitfun_bkg, fitv2low, fitv2hi,2);
	  fitfunbkg->SetParameters(fitfun->GetParameter(1), fitfun->GetParameter(2));}
      
      
      //fitfunbkg->SetLineColor(kGreen+2);
      fitfunbkg->SetLineColor(6);
      fitfunbkg->SetLineStyle(4);
      fitfunbkg->DrawCopy("same");
   

      //canfit2[ip]->cd(2);
      //dummy_prof[ip]->Draw("");
      //dummy_prof[ip]->GetYaxis()->SetRangeUser(-0.04, 0.04);
      //fitfun->DrawCopy("same l");
      //fitfunsig->DrawCopy("same l");
      //fitfunbkg->DrawCopy("same l");

      TLegend * lsig2 = DrawLegend();
      lsig2->AddEntry(fitfun,"v^{Sig+Bg}_{obs}","l");
      lsig2->AddEntry(fitfunsig,"#frac{Sig}{Sig+Bg} v^{Sig}_{obs}","l");
      lsig2->AddEntry(fitfunbkg,"#frac{Bg}{Sig+Bg} v^{Bg}_{obs}","l");
      lsig2->Draw();
   
      //return;
      Char_t outfig[128];
      //sprintf(outfig,"Plots/InvMass_v2_MethodOneP%dCentrality%d_%d%s_BgV2%s.gif",ip+1,CentLow,CentHigh,BgType,Bgv2.Data());
      //sprintf(outfig,"Plots/InvMass_v1_%Rap%d_Cent%d%s_BgV2%s.gif",ip, ip+1,CentLow,CentHigh,BgType,Bgv2.Data());
      //canfit1->SaveAs(outfig);
      canfit2[ip]->SaveAs(Form("Plots/InvMassv1_%dRap%d_%dDphi%d_Cent_%d%d_%s.gif",ip,ip+1,philow,phihi,CentLow,CentHigh,BgType));
      canfit2[ip]->Write();

      v2obs[ip]=fitfun->GetParameter(0)/reso;
      v2obserr[ip]=fitfun->GetParError(0)/reso;
      chisquare[ip]=fitfun->GetChisquare();
      ndf[ip]=fitfun->GetNDF();
      chisquarendf[ip]=chisquare[ip]/(Double_t)ndf[ip];
      canfit2[ip]->Update();
    }

  //resolution
  /* 
  if(reso == 1.0){
  if(energy == 7 && CentLow==0 && CentHigh==10) reso = 0.225812;
  if(energy == 7 && CentLow==10 && CentHigh==40) reso = 0.534333;
  if(energy == 7 && CentLow==40 && CentHigh==80) reso = 0.586273;
  if(energy == 11 && CentLow==0 && CentHigh==10) reso = 0.22533;
  if(energy == 11 && CentLow==10 && CentHigh==40) reso = 0.518218;
  if(energy == 11 && CentLow==40 && CentHigh==80) reso = 0.531816;
  if(energy == 19 && CentLow==0 && CentHigh==10) reso = 0.192782;;
  if(energy == 19 && CentLow==10 && CentHigh==40) reso = 0.401668;
  if(energy == 19 && CentLow==40 && CentHigh==80) reso = 0.359168;
  if(energy == 27 && CentLow==0 && CentHigh==10) reso = 0.160824;;
  if(energy == 27 && CentLow==10 && CentHigh==40) reso = 0.312083;
  if(energy == 27 && CentLow==40 && CentHigh==80) reso = 0.24177;
  if(energy == 39 && CentLow==0 && CentHigh==10) reso = 0.124849;
  if(energy == 39 && CentLow==10 && CentHigh==40) reso = 0.202612;
  if(energy == 39 && CentLow==40 && CentHigh==80) reso = 0.102032;
  if(energy == 62 && CentLow==0 && CentHigh==10) reso = 0.0;
  if(energy == 62 && CentLow==10 && CentHigh==40) reso = 0.0982692;
  if(energy == 62 && CentLow==40 && CentHigh==80) reso = 0.159668;
  if(energy == 62 && CentLow==10 && CentHigh==70) reso = 0.110072;
  }
*/
  //save v1-data in txt
  ofstream outf(Form("v1_rap_cent%d%d_%d.dat",CentLow, CentHigh, energy));
  



  for(Int_t ip = 0; ip < nRapBin-1; ip++){
    //cout<<    RawYield[ip] <<"  "<<  ErRawYield[ip]<<endl;
    //cout<<pT[ip]<<"  "<< v2obs[ip] <<"  "<<  v2obserr[ip]<<"  "<<NorFact[ip]<<endl;
    //v2obs[ip] = v2obs[ip]/0.893451;
    //v2obserr[ip]= v2obserr[ip]/0.893451;
    //v2obs[ip] = v2obs[ip]/0.87445;
    //v2obserr[ip]= v2obserr[ip]/0.87445;
    //v2obs[ip] = v2obs[ip]/0.820889;
    //v2obserr[ip]= v2obserr[ip]/0.820889;
    //cout<<rapi[ip]<<"  "<< v2obs[ip]/reso <<"  "<< err_rapi[ip]<<" " << v2obserr[ip]/reso<<"  "<<endl;
    //outf<< rapi[ip]<<"  "<< v2obs[ip]/reso <<"  "<< err_rapi[ip]<<" "<< v2obserr[ip]/reso <<"  "<<endl;

    cout<<rapi[ip]<<"  "<< v2obs[ip] <<"  "<< err_rapi[ip]<<" " << v2obserr[ip]<<"  "<<endl;
    outf<< rapi[ip]<<"  "<< v2obs[ip] <<"  "<< err_rapi[ip]<<" "<< v2obserr[ip] <<"  "<<endl;
  }
  //return;//subhash checked

  
  TCanvas *c = new TCanvas("c","c",10,10,600,600);
  c->cd();
  c->SetLeftMargin(0.15);
  
  TGraphErrors *gr = PlotGraph(nRapBin-1, 1, 20, rapi, err_rapi, v2obs, v2obserr);
  gr->GetYaxis()->SetTitle("v_{1}^{obs}");
  gr->GetXaxis()->SetTitle("y");
  gr->SetName("gr_v1_y");
  //gr->GetYaxis()->SetRangeUser(-0.1, 0.1);
  //gr->GetXaxis()->SetRangeUser(-0.7, 0.7);
  gr->GetYaxis()->SetRangeUser(-0.05, 0.2);
  gr->GetXaxis()->SetRangeUser(-0.1, 3.1);
  gr->Draw("AP");
  c->SaveAs(Form("Plots/Obsv1_Cent_%d%d_%d.gif",CentLow,CentHigh,energy));
  TF1 *fslope = new TF1("fslope", "[0]*x", -0.5, 0.5); 
  fslope->SetParameter(0, 0.001);
  //TF1 *fslope = new TF1("fslope", "[0]*x+[1]*x*x*x", -0.6, 0.6); 
  gr->Fit(fslope,"ERIM");
  gr->Write();
  //return;

  TGraphErrors *gryld_y = PlotGraph(nRapBin-1, 1, 24, rapi, err_rapi, RawYield, ErRawYield);
  TCanvas *c11 = new TCanvas("c11","c11",10,10,600,600);
  c11->SetBottomMargin(0.15);
  c11->SetLeftMargin(0.2);
  c11->cd();
  gryld_y->GetYaxis()->SetTitle("raw yield");
  gryld_y->GetXaxis()->SetTitle("y");
  gryld_y->SetName("gr_yld_y");
  gryld_y->Draw("ape");
  gryld_y->Write();

  TGraphErrors *gr_significance = PlotGraph(nRapBin-1, 1, 24, rapi, err_rapi, significance, significance_err);
  TCanvas *c12 = new TCanvas("c12", "c12", 10, 10, 600, 600);
  c12->SetBottomMargin(0.15);
  c12->SetLeftMargin(0.2);
  c12->cd();
  gr_significance->GetYaxis()->SetTitle("Significance");
  gr_significance->GetXaxis()->SetTitle("y");
  gr_significance->SetName("gr_significance");
  gr_significance->Draw("ape");
  gr_significance->Write();

  TGraphErrors *gr_SigByBg = PlotGraph(nRapBin-1, 1, 24, rapi, err_rapi, SigByBg, SigByBg_err);
  TCanvas *c13 = new TCanvas("c13", "c13", 10, 10, 600, 600);
  c13->SetBottomMargin(0.15);
  c13->SetLeftMargin(0.2);
  c13->cd();
  gr_SigByBg->GetYaxis()->SetTitle("Signal/Background");
  gr_SigByBg->GetXaxis()->SetTitle("y");
  gr_SigByBg->SetName("gr_SigByBg");
  gr_SigByBg->Draw("ape");
  gr_SigByBg->Write();
  return;
  /*
  //TGraphErrors *gr1 = PlotGraph(pTBins,2,20,pT,ErpT,dNdpTdy,ErdNdpTdy);
  //gr1->GetXaxis()->SetTitle("p_{T} ( GeV/c)");
  //gr1->GetYaxis()->SetTitle("#frac{1}{2#pi p_{T}}#frac{dN^{2}}{dp_{T}dy}");
  //gr1->Draw("AP");
  //c1->SetLogy();
  
  //
  TCanvas *c2 = new TCanvas("c2","c",10,10,600,600);
  c2->SetBottomMargin(0.15);
  c2->SetLeftMargin(0.2);
  c2->cd();

  TGraphErrors *gr2 = PlotGraph(nRapBin-1,4,20,rapi,err_rapi,Width,ErWidth);
  gr2->SetMaximum(0.01);
  gr2->SetMinimum(0.0);
  gr2->GetXaxis()->SetTitle("y");
  gr2->GetYaxis()->SetTitle("Width (GeV/c^{2})");
  gr2->Draw("AP");
  //TLine * line1 = DrawLine(0.,0.00426,rapi[pTBins-1]+err_rapi[pTBins-1],0.00426,2,2,4);
  //line1->Draw();
  //TLatex* tex1  = DrawText(5.0, 0.00436, 2, "PDG Value");
  //tex1->Draw();
  //
  TCanvas *c3 = new TCanvas("c3","c",10,10,600,600);
  c3->SetBottomMargin(0.15);
  c3->SetLeftMargin(0.2);
  c3->cd();
  TGraphErrors *gr3 = PlotGraph(nRapBin-1,6,20,rapi,err_rapi,Mass,ErMass);
  gr3->SetMaximum(1.025);
  gr3->SetMinimum(1.014);
  gr3->GetXaxis()->SetTitle("p_{T} ( GeV/c)");
  gr3->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
  gr3->Draw("AP");
  //TLine * line2 = DrawLine(0.0,1.019455,pT[pTBins-1]+ErpT[pTBins-1],1.019455,2,2,4);
  //line2->Draw();
  //TLatex* tex2  = DrawText(1.0, 1.02, 2, "PDG Value");
  //tex2->Draw();
  //
  TCanvas *c4 = new TCanvas("c4","c",10,10,600,600);
  c4->SetBottomMargin(0.15);
  c4->SetLeftMargin(0.2);
  c4->cd();
  TGraphErrors *gr4 = PlotGraph(nRapBin-1,6,20,rapi,err_rapi,NorFact,0);
  gr4->SetMaximum(0.5);
  gr4->SetMinimum(0.01);
  gr4->GetXaxis()->SetTitle("p_{T} ( GeV/c)");
  gr4->GetYaxis()->SetTitle("Normalization Factor");
  gr4->Draw("AP");

*/
  /*
  Bool_t SaveDataInFile = kFALSE;
  cout<<"Please Type : 1 (To save data in file) : 0 (Not to save data in file)"<<endl;
  //cin>>SaveDataInFile;
  if(SaveDataInFile)
    {
      ofstream myfile;
      char file[128];
      sprintf(file,"Plots/%d_%d/InvMassMeth/DataInvariantMassMethodOneCentrality%d_%d_Ybin%d_%d%s.txt",CentLow,CentHigh,CentLow,CentHigh,philow,phihi,BgType);
      myfile.open(file);
      for(Int_t i = 0; i < pTBins; i++){
	myfile<<pT[i]<<"  "<<ErpT[i]<<"  "<<dNdpTdy[i]<<"  "<<ErdNdpTdy[i]<<"  "<<Width[i]<<"  "<<ErWidth[i]<<"  "<<Mass[i]<<"  "<<ErMass[i]<<"  "<<RawYield[i]<<"  "<<ErRawYield[i]<<"  "<<NorFact[i]<<"  "<<BinCount[i]<<endl;
	cout<<pT[i]<<"  "<<RawYield[i]<<"  "<<ErRawYield[i]<<"  "<<NorFact[i]<<"  "<<BinCount[i]<<endl;
      }
      myfile.close();

      FILE *file1;
      if(bincount)
	file1 = fopen(Form("InvariantMassV2AtDiffPtBinsMethodOneBinCount_BgV2%s.txt",Bgv2.Data()),"a");
      else
	file1 = fopen(Form("InvariantMassV2AtDiffPtBinsMethodOne%s_BgV2%s.txt",BgType,Bgv2.Data()),"a");
      fprintf(file1,"%6d %6d %6d %6d\n",CentLow, CentHigh,philow,phihi);
      for(Int_t ip = 0; ip < pTBins; ip++){
	fprintf(file1,"%10.5f %10.5f %10.5f %10.5f%10.5f \n",pT[ip],ErpT[ip],v2obs[ip],v2obserr[ip],chisquarendf[ip]);
      }
      fclose(file1);
      c1->cd();
      c1->SaveAs(Form("Plots/%d_%d/InvMassMeth/PtSpectraCentralityBin%d_%d_%dYbin%dInvMassMethOne%s.gif",CentLow,CentHigh,CentLow,CentHigh,philow,phihi,BgType));
      c2->cd();
      c2->SaveAs(Form("Plots/%d_%d/InvMassMeth/PtVsMassCentralityBin%d_%d_%dYbin%dInvMassMethOne%s.gif",CentLow,CentHigh,CentLow,CentHigh,philow,phihi,BgType));
      c3->cd();
      c3->SaveAs(Form("Plots/%d_%d/InvMassMeth/PtVsWidthCentralityBin%d_%d_%dYbin%dInvMassMethOne%s.gif",CentLow,CentHigh,CentLow,CentHigh,philow,phihi,BgType));
      c4->cd();
      c4->SaveAs(Form("Plots/%d_%d/InvMassMeth/PtVsNorFactorCentralityBin%d_%d_%dYbin%dInvMassMethOne%s.gif",CentLow,CentHigh,CentLow,CentHigh,philow,phihi,BgType));
    }
*/
}

Double_t fitfun_raw(Double_t *x, Double_t *para)
{
  Double_t fitval;
  Int_t bin = (Int_t)((x[0] - mass_lowfunc)/(0.001));

  if(Bgv2.CompareTo("Poly3rd") == 0)
    fitval = para[0]*ratio_sig_total[bin]+(para[1]+para[2]*x[0]+para[3]*x[0]*x[0]+para[4]*x[0]*x[0]*x[0])*(1-ratio_sig_total[bin]);
  if(Bgv2.CompareTo("Poly2nd") == 0)
    fitval = para[0]*ratio_sig_total[bin]+(para[1]+para[2]*x[0]+para[3]*x[0]*x[0])*(1-ratio_sig_total[bin]);
  if(Bgv2.CompareTo("Poly1st") == 0)
    fitval = para[0]*ratio_sig_total[bin]+ (para[1]+para[2]*x[0])*(1-ratio_sig_total[bin]);
  return fitval;
}
//---------------------------------------------------------------

Double_t fitfun_sig(Double_t *x, Double_t *para)
{
  Double_t fitval;
  Int_t bin = (Int_t)((x[0] - mass_lowfunc)/(0.001));
  fitval = para[0]*ratio_sig_total[bin];
  return fitval;
}

Double_t fitfun_bkg(Double_t *x, Double_t *para)
{
  Double_t fitval;
  Int_t bin = (Int_t)((x[0] - mass_lowfunc)/(0.001));

  if(Bgv2.CompareTo("Poly3rd") == 0)
    fitval = (para[0]+para[1]*x[0]+para[2]*x[0]*x[0]+para[3]*x[0]*x[0]*x[0])*(1-ratio_sig_total[bin]);
  if(Bgv2.CompareTo("Poly2nd") == 0)
    fitval = (para[0]+para[1]*x[0]+para[2]*x[0]*x[0])*(1-ratio_sig_total[bin]);
  if(Bgv2.CompareTo("Poly1st") == 0)
    fitval = (para[0]+para[1]*x[0])*(1-ratio_sig_total[bin]);
  return fitval;
}
//---------------------------------------------------------------
//=================================================================================================
TH1D * MakeDiff(Int_t reps, TH1D *hsig, TH1D *hbg, Double_t peakL, Double_t peakU, Double_t &scale,TFitResultPtr &r,Double_t lowrange,Double_t highrange, Double_t Para1)
{
  TH1::AddDirectory(0);
  const Int_t iters = reps;
  //Declare arrays of histos to be used in iterations
  TH1D *hdiff[iters]; //difference histo
  TH1D *hbgs[iters];  //scaled bg histo
  
  //Fit functions for difference histos
  TF1 *fBWBg[iters];//Total function (fBW + fBg)
  TF1 *fBW[iters]; //Breit-Wigner function alone
  
  TString diffname = hsig->GetName();
  diffname += "diff";
  TString bgname = hbg->GetName();
  TString formname = hsig->GetName();
  formname += "func";
  TString breitname = hsig->GetName();
  breitname += "BW";
  
  //Arrays to contain:
  Double_t IntS[iters];  //signal integral
  Double_t scale[iters]; //scaling factor
  
  Double_t IntBG = hbg->Integral(hbg->GetXaxis()->FindBin(lowrange), hbg->GetXaxis()->FindBin(highrange), "width");
  cout<<"1) BG integral                  = "<<IntBG<<endl;
  cout<<"   Integrating Background from bin = "<<hbg->GetXaxis()->FindBin(lowrange)<< " to bin = "<<hbg->GetXaxis()->FindBin(highrange)<<endl;

  //Begin iterations
  for(Int_t i = 0; i < iters; i++)
    {
      //Calculate the scaling factor
      if(i == 0)
	{
	  //first time, signal integral = everything under peak
	  IntS[i] =  hsig->Integral(hsig->GetXaxis()->FindBin(lowrange), hsig->GetXaxis()->FindBin(highrange), "width");
	}
      else
	{
	  //every other time, signal integral = total - Breit-Wig fit from time before
	  IntS[i] =  IntS[0] - fBW[i-1]->Integral(lowrange, highrange);
	}
      
      scale[i] = IntS[i]/IntBG;
      cout<<i<<" Signal integral = "<<IntS[i]<<endl;
      cout<<i<<" BG integral     = "<<IntBG<<endl;
      cout<<i<<" Scale factor    = "<<scale[i]<<endl;
      
      //Scale background histo
      TString bname = bgname;
      bname += i;
      hbgs[i] = (TH1D *)hbg->Clone(bname.Data());
      hbgs[i]->Scale(scale[i]);
      
      //Make difference histo
      TString dname = diffname;
      dname += i;
      
      hdiff[i] = (TH1D *)hsig->Clone(dname.Data());
      hdiff[i]->Add(hbgs[i], -1.0);
      hdiff[i]->GetXaxis()->SetRangeUser(0.98, 1.1);
      hdiff[i]->SetLineColor(4);
      hdiff[i]->SetTitle("Signal-BG_{scaled}");
      //hdiff[i]->Sumw2();

      //Fit difference histo with Breit-Wigner + straight line
      //WITH 2PI Factor
      TString fname = formname;
      fname += i;
      if(fit==0){
	fBWBg[i] = new TF1(fname.Data(),Form("1.0*0.001*%d*[0]*[1]/(2*3.1415926)/(pow(x - [2],2) + [1]*[1]/4) + [3]+[4]*x",rebin),lowrange, highrange);   //rebin*unitmass
      }
      else if(fit == 1){
	fBWBg[i] = new TF1(fname.Data(),Form("1.0*0.001*%d*[0]*[1]/(2*3.1415926)/(pow(x - [2],2) + [1]*[1]/4) + [3]+[4]*x+[5]*x*x",rebin),lowrange, highrange[1]);  //rebin*unitmass
      }
      fBWBg[i]->SetParName(0, "BW Area");
      fBWBg[i]->SetParName(1, "#Gamma");
      fBWBg[i]->SetParName(2, "Mass");
      fBWBg[i]->SetParName(3, "pol0");
      fBWBg[i]->SetParName(4, "pol1");
      if(fit == 0){
	//fBWBg[i]->SetParameters(Para1, 0.0055, 1.019, 4.52218e+02, -1.41651e+03);
	fBWBg[i]->SetParameters(IntS[i], 0.0055, 1.019, 4.52218e+02, -1.41651e+03);
      }
      if(fit == 1){
	fBWBg[i]->SetParName(5, "plo2");
	fBWBg[i]->SetParameters(Para1, 0.0055, 1.019, 4.52218e+05, -1.41651e+05, 9.53836e+05);
	//fBWBg[i]->SetParameters(IntS[i], 0.0055, 1.019, -4.52218e+05, 1.41651e+05, -9.53836e+05);
      }
      
      
      cout<<"FUNCTION NAME ================================> "<<fname.Data()<<"  "<<IntS[i]<<endl;
      //fBWBg[i]->FixParameter(1,0.0055);  
      fBWBg[i]->SetParLimits(1, 0.003, 0.01);
      fBWBg[i]->SetParLimits(2, peakL, peakU);
      //fBWBg[i]->SetParLimits(3, 2.82285e+02, 2.82285e+05);
      //fBWBg[i]->SetParLimits(4, -4.60518e+05, 4.60518e+05);

      if(opt==1)
	r = hdiff[i]->Fit(fBWBg[i], "ERS");
      //r = hdiff[i]->Fit(fBWBg[i], "ERMIS");
      else r = hdiff[i]->Fit(fBWBg[i], "ERS");
      fBWBg[i]->SetLineColor(1);
      fBWBg[i]->Draw("LSAME");
      Double_t *PARER = fBWBg[i]->GetParErrors();
      //Make Breit-Wigner function alone
      TString bwname = breitname;
      bwname += i;
      fBW[i] = new TF1(bwname.Data(),Form("(1.0*0.001*%d/(2*3.14159))*(([0]*[1])/((x-[2])**2 +([1]/2)**2))",rebin), lowrange, highrange);
      fBW[i]->SetParameters(fBWBg[i]->GetParameter(0), fBWBg[i]->GetParameter(1), fBWBg[i]->GetParameter(2));
      fBW[i]->SetParErrors(&PARER[0]);
      
    }
  /*
    cout<<"Integrated AREA          ------------>>>>>   "<<
    (fBWBg[reps-1]->Integral(lowrange, highrange)/(0.001*rebin))<<endl;
    cout<<"Integrated AREA  Error   ------------>>>>>   "<<
    (fBWBg[reps-1]->IntegralError(lowrange, highrange)/(0.001*rebin))<<endl;
  */
  scale = scale[reps-1];
  return hdiff[reps-1];
}
//====================================================================================



TGraphErrors * PlotGraph(Int_t NdataPoint, Int_t MarkerColor = 1, Int_t MarkerStyle = 20, Double_t *X, Double_t *ErX,Double_t *Y, Double_t *ErY)
{
  TGraphErrors *gr = new TGraphErrors(NdataPoint, X ,Y , ErX, ErY);
  gr->SetTitle("");
  gr->SetMarkerStyle(MarkerStyle);
  gr->SetMarkerColor(MarkerColor);
  gr->SetMarkerSize(1.5);
  //gr->GetXaxis()->SetTitle("p_{T} ( GeV/c)");
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetLabelFont(42);
  //gr->GetYaxis()->SetTitle("#Phi Mass ( GeV/c^{2})");
  gr->GetXaxis()->CenterTitle(true);
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->CenterTitle(true);
  gr->GetYaxis()->SetTitleOffset(1.7);
  return gr;
}

//==========================================
TLine *DrawLine(Double_t xmin = 0,Double_t ymin = 1, Double_t xmax = 1, 
		Double_t ymax = 1,Int_t lStyle = 1, Int_t lColor = 1, Int_t lWidth = 1)
{
  TLine *line = new TLine(xmin,ymin,xmax,ymax);
  line->SetLineStyle(lStyle);
  line->SetLineColor(lColor);
  line->SetLineWidth(lWidth);
  //line->Draw();
  return line;
}
//==========================================
TLatex *DrawText(Double_t x = 0, Double_t y = 0,Int_t tColor = 2,TString name)
{
  TLatex* tex = new TLatex(x,y,name.Data());
  tex->SetTextSize(0.04);
  tex->SetTextColor(tColor);
  tex->SetTextFont(42);
  //tex->Draw();
  return tex;
}
//==================================
TLegend *DrawLegend()
{
  TLegend *legend = new TLegend(0.6,0.65,0.88,0.85);
  legend->SetTextFont(42);
  // legend->SetTextSize(0.04);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  //legend->AddEntry(gr1,"(0 - 100) %","p");
  //legend->AddEntry(func1,"p_{0}[ 1 + 2 v_{2}^{obs} cos(2(#Phi - #Psi))]","l");
  return legend;
}

//new add
Double_t Invm_func_tot(Double_t *x, Double_t *par)
{
  Int_t bin =hf->FindBin(x[0]);
  Double_t y =hf->GetBinContent(bin);
  //Int_t bin1 =hfb->FindBin(x[0]);
  //Double_t y1 =hfb->GetBinContent(bin1);
  Double_t vnTot;
  //Double_t vnBG = ;
  
  //if (x[0] > fLeftOfSignal(currentMeanPt) && x[0] < fRightOfSignal(currentMeanPt)) 
  if (x[0] > 0.92 && x[0] < 1.12) 
  //if (x[0] > 0.99 && x[0] < 1.044) 
     {
      
       if(Bgv2.CompareTo("Poly1st") == 0)
	 vnTot =  par[0]*y + (par[1] + par[2]*x[0])*(1-y);
       if(Bgv2.CompareTo("Poly2nd") == 0)
	 vnTot =  par[0]*y + (par[1] + par[2]*x[0] + par[3]*x[0]*x[0])*(1-y);
       if(Bgv2.CompareTo("Poly3rd") == 0)
	 vnTot =  par[0]*y + (par[1] + par[2]*x[0] + par[3]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0])*(1-y);
     } 
   else
     {
       if(Bgv2.CompareTo("Poly1st") == 0)
	 vnTot =  par[0]*y + (par[1] + par[2]*x[0]);
       if(Bgv2.CompareTo("Poly2nd") == 0)
	 vnTot =  par[0]*y + (par[1] + par[2]*x[0] + par[3]*x[0]*x[0]);
       if(Bgv2.CompareTo("Poly3rd") == 0)
	 vnTot =  par[0]*y + (par[1] + par[2]*x[0] + par[3]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0]);
     }
   return vnTot;
}

Double_t Invm_func_bkg(Double_t *x, Double_t *par)
{


   Int_t bin =hf->FindBin(x[0]);
   Double_t y =hf->GetBinContent(bin);
   //Int_t bin1 =hfb->FindBin(x[0]);
   //Double_t y1 =hfb->GetBinContent(bin1);
   Double_t vnBG;// = (par[0] + par[1]*x[0] + par[2]*x[0]*x[0])*(1-y-y1);
   
   if(Bgv2.CompareTo("Poly3rd") == 0)
     vnBG = (par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0])*(1-y);
   if(Bgv2.CompareTo("Poly2nd") == 0)
     vnBG = (par[0] + par[1]*x[0] + par[2]*x[0]*x[0])*(1-y);
   if(Bgv2.CompareTo("Poly1st") == 0)
     vnBG = (par[0] + par[1]*x[0])*(1-y);
   return vnBG;

}

Double_t Invm_func_sig(Double_t *x, Double_t *par)
{


   Int_t bin =hf->FindBin(x[0]);
   Double_t y =hf->GetBinContent(bin);
   //Double_t vnBG = (p0_p + p1_p*x[0] + p2_p*x[0]*x[0])(1-y);
  Double_t vnSig = par[0]*y;
  return vnSig;

}
