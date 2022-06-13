// C++ headers
#include <iostream>
#include <fstream>
#include "Riostream.h"
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <map>

//need this stuff to compile in linux:
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStreamerElement.h>
#include <TStyle.h>
#include "TSystemDirectory.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TMath.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TGaxis.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TPaveText.h"
#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "Math/MinimizerOptions.h"

using namespace std;
// -------------------------- set Some fitting prerequsites --------------------

const Double_t _sigmaRange = 5.; // Sigma of the Fitting range
const Double_t _y_CM = -2.03;
const Double_t _n_jkk = 1;
Double_t dParBg[3]; // Bkg fitting parameters
Double_t dParSig[4]; // Sig + Bkg fitting parameters
Double_t proportion(Double_t *x, Double_t *p);
Double_t BackgroundFitting(Double_t *x, Double_t *p);
Double_t TotalFitting(Double_t *x, Double_t *p);

// ======================== (1) Analysis Start =================================
void FlowExtractor( /*TString invMFileName = "./res_sys/result_sys_invM/merged_merged_sys_primary_var0_iter1_.root",*/
                   // TString FlowFileName =  "./res_sys/result_sys_flow/hadd_PhiMesonAna_OUTPUT_sys_primary_var0_iter3_.root" ,
                   TString FlowFileName =  "/mnt/c/Users/pjska/ana/7p2gev_Phi_v2/FlowExtractor/input/merged_merged_PhiMesonAna_OUTPUT_sys_primary_var0_iter3_BFD568E307BF80B18B7C0AD173647968_.root" ,
                    // double inputParameter1 = 0.
				   TString out_file = "minus1_v2_plus1",
                    Int_t   inputp2 = 0, // sysErr cut Indexes 0-15
                    Int_t   inputp3 = 0, // sysErr cut variations, each systematic check has 2 or 3 vertions
                    Int_t   inputp4 = 3 // Iteration of the analysis is. In this analysis, 2 iterations is enough
){
  Int_t sys_cutN = inputp2; // sysErr cut Indexes 0-15
  Int_t sys_varN = inputp3; // sysErr cut variations, each systematic check has 2 or 3 vertions
  Int_t sys_iterN = inputp4; // Iteration of the analysis is. In this analysis, 2 iterations is enough
  string sys_object[17]  = {"primary", "etaGap", "etaRange",
                            "vz", "vr", "dedx", "dca",
                            "nHitsFit", "ratio", "nSigK", "mass2",
                            "pT", "dipAngle", "vtxDiff", "mthdDiff",
                            "binning",
                            "TPCpid"};
  std::cout << "sys_cutN == "<< sys_cutN <<": "<< sys_object[sys_cutN] << std::endl;
  // TString outTxt = "./out_sys/out_sys_Crosscheck/newProd_3ybin_3Sig_phi_v1_y_sys_";
  TString outTxt = "./out/phi_v2_pt_sys_";
  TString outHead = outTxt;
  outTxt.Append(sys_object[sys_cutN]);
  outTxt.Append(Form("_var%d_iter%d_", sys_varN, sys_iterN));
  outTxt.Append(out_file);
  outTxt.Append(".txt");
  std::ofstream flowFile(outTxt,ofstream::out);
  // ---------------------- Analysis Setup -------------------------------------
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  // ----- InvMass plots in different centrality and pT or y bins --------------
  Double_t ptSetA[3]  = {0.6, 1.2, 2.4};
  Double_t ptSetB[5]  = {0.4, 0.7, 1.0, 1.4, 2.0};
  Double_t ptSetC[11] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0};

  Double_t rapSetA[5]  = {-2.0, -1.5, -1.0, -0.5, 0};

  Double_t centSetA[5]  = {0, 10, 40, 60, 80}; // %
  Double_t centSetB[10]  = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80}; // %
  // directed and elliptic flow. Indexes: 0: v1, 1: v2; raw, reso; pT/y SetA; cent SetA; jkk
  Double_t d_FLow_ptSetA_centSetA[11][2][2][2][6]; // pt SetA, cent SetA
  Double_t d_FLow_ptSetA_centSetB[2][2][2][9]; // pt SetA, cent SetB
  Double_t d_FLow_ptSetB_centSetA[2][2][4][6]; // pt SetA, cent SetA
  Double_t d_FLow_ptSetB_centSetB[2][2][4][6]; // pt SetB, cent SetB
  // Double_t d_FLow_ptSetC_centAll[2][2][10][2]; // pt SetC, cent 0-60%, 0-80%
  // Double_t d_FLow_rapSetA_centSetA[2][2][4][6]; // pt SetA, cent 0-60%, 0-80%
  // Double_t d_FLow_rapSetA_centSetA_pTRange[2][2][4][6][3]; // pt SetC, cent 0-60%, 0-80%, 3 pT range, [0.1,1], [1,2], [0.1,2]
  // Double_t d_FLow_rapSetA_centSetB[2][2][4][9]; // pt SetC, cent 0-60%, 0-80%


  Double_t d_Flow_err_ptSetA_centSetA[11][2][2][2][6]; // pt SetA, cent SetA
  Double_t d_Flow_err_ptSetA_centSetB[2][2][2][9]; // pt SetA, cent SetB
  Double_t d_Flow_err_ptSetB_centSetA[2][2][4][6]; // pt SetA, cent SetA
  Double_t d_Flow_err_ptSetB_centSetB[2][2][4][6]; // pt SetB, cent SetB
  // Double_t d_Flow_err_ptSetC_centAll[2][2][10][2]; // pt SetC, cent 0-60%, 0-80%
  // Double_t d_Flow_err_rapSetA_centSetA[2][2][4][6]; // pt SetC, cent 0-60%, 0-80%
  // Double_t d_Flow_err_rapSetA_centSetA_pTRange[2][2][4][6][3]; // pt SetC, cent 0-60%, 0-80%, 3 pT range, [0.1,1], [1,2], [0.1,2]
  // Double_t d_Flow_err_rapSetA_centSetB[2][2][4][9]; // pt SetC, cent 0-60%, 0-80%
  // ---------------------- Input files and plots ------------------------------
  // SE/ME invM input
  // TFile * file_KK_InvM_Input = new TFile(invMFileName,"READ");
  // if( !file_KK_InvM_Input->IsOpen() ) std::cout<<"No SE/ME input!"<<std::endl;
  // if(  file_KK_InvM_Input->IsOpen() ) {
  //     std::cout<<"#phi InvM loaded successfully!"<<std::endl;
  // }
  // flow VS Invariant Mass input
  // TFile * file_flow_invM_Input = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/res/merged_merged_PhiMesonAna_OUTPUT_F7793427B87FC5429328F2DB142A9B34_.root","READ");
  // Default phi-flow
  TFile * file_flow_invM_Input = new TFile(FlowFileName,"READ");
  if( !file_flow_invM_Input->IsOpen() ) std::cout<<"No flow input!"<<std::endl;
  if(  file_flow_invM_Input->IsOpen() ) {
      std::cout<<"flow file loaded successfully!"<<std::endl;
  }
  // pt SetA, cent SetA
  TH1D *mHist_Input_SE_InvM_ptSetA_centSetA[11][2][6];
  TH1D *mHist_Input_ME_InvM_ptSetA_centSetA[11][2][6];
  TProfile *mProfile_Input_v1_raw_ptSetA_centSetA[11][2][6];
  TProfile *mProfile_Input_v1_reso_ptSetA_centSetA[11][2][6];
  TProfile *mProfile_Input_v2_raw_ptSetA_centSetA[11][2][6];
  TProfile *mProfile_Input_v2_reso_ptSetA_centSetA[11][2][6];
  for(int jkk=0; jkk<_n_jkk; jkk++){
    for(int pt=0; pt<2; pt++)
    {
      for(int cent=0; cent<6;cent++){
        mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent] = (TH1D*) file_flow_invM_Input/*file_KK_InvM_Input*/->Get(Form("Hist_SE_InvM_ptSetA%d_centSetA%d",pt,cent));
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] = (TH1D*) file_flow_invM_Input/*file_KK_InvM_Input*/->Get(Form("Hist_rotation_InvM_ptSetA%d_centSetA%d",pt,cent));
        mProfile_Input_v1_raw_ptSetA_centSetA[jkk][pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v1_raw_ptSetA%d_centSetA%d_pfx",pt,cent));
        mProfile_Input_v1_reso_ptSetA_centSetA[jkk][pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v1_reso_ptSetA%d_centSetA%d_pfx",pt,cent));
        mProfile_Input_v2_raw_ptSetA_centSetA[jkk][pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v2_raw_ptSetA%d_centSetA%d_pfx",pt,cent));
        mProfile_Input_v2_reso_ptSetA_centSetA[jkk][pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v2_reso_ptSetA%d_centSetA%d_pfx",pt,cent));
      }
    }
  }
  // pt SetA, cent SetB
  TH1D *mHist_Input_SE_InvM_ptSetA_centSetB[2][9];
  TH1D *mHist_Input_ME_InvM_ptSetA_centSetB[2][9];
  TProfile *mProfile_Input_v1_raw_ptSetA_centSetB[2][9];
  TProfile *mProfile_Input_v1_reso_ptSetA_centSetB[2][9];
  TProfile *mProfile_Input_v2_raw_ptSetA_centSetB[2][9];
  TProfile *mProfile_Input_v2_reso_ptSetA_centSetB[2][9];
  for(int pt=0; pt<2; pt++)
  {
    for(int cent=0; cent<9;cent++){
      mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent] = (TH1D*) file_flow_invM_Input/*file_KK_InvM_Input*/->Get(Form("Hist_SE_InvM_ptSetA%d_centSetB%d",pt,cent));
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] = (TH1D*) file_flow_invM_Input/*file_KK_InvM_Input*/->Get(Form("Hist_rotation_InvM_ptSetA%d_centSetB%d",pt,cent));
      mProfile_Input_v1_raw_ptSetA_centSetB[pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v1_raw_ptSetA%d_centSetB%d_pfx",pt,cent));
      mProfile_Input_v1_reso_ptSetA_centSetB[pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v1_reso_ptSetA%d_centSetB%d_pfx",pt,cent));
      mProfile_Input_v2_raw_ptSetA_centSetB[pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v2_raw_ptSetA%d_centSetB%d_pfx",pt,cent));
      mProfile_Input_v2_reso_ptSetA_centSetB[pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v2_reso_ptSetA%d_centSetB%d_pfx",pt,cent));
    }
  }
  // pt SetB, cent SetA
  TH1D *mHist_Input_SE_InvM_ptSetB_centSetA[4][6];
  TH1D *mHist_Input_ME_InvM_ptSetB_centSetA[4][6];
  TProfile *mProfile_Input_v1_raw_ptSetB_centSetA[4][6];
  TProfile *mProfile_Input_v1_reso_ptSetB_centSetA[4][6];
  TProfile *mProfile_Input_v2_raw_ptSetB_centSetA[4][6];
  TProfile *mProfile_Input_v2_reso_ptSetB_centSetA[4][6];
  for(int pt=0; pt<4; pt++)
  {
    for(int cent=0; cent<6;cent++){
      mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent] = (TH1D*) file_flow_invM_Input/*file_KK_InvM_Input*/->Get(Form("Hist_SE_InvM_ptSetB%d_centSetA%d",pt,cent));
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] = (TH1D*) file_flow_invM_Input/*file_KK_InvM_Input*/->Get(Form("Hist_rotation_InvM_ptSetB%d_centSetA%d",pt,cent));
      mProfile_Input_v1_raw_ptSetB_centSetA[pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v1_raw_ptSetB%d_centSetA%d_pfx",pt,cent));
      mProfile_Input_v1_reso_ptSetB_centSetA[pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v1_reso_ptSetB%d_centSetA%d_pfx",pt,cent));
      mProfile_Input_v2_raw_ptSetB_centSetA[pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v2_raw_ptSetB%d_centSetA%d_pfx",pt,cent));
      mProfile_Input_v2_reso_ptSetB_centSetA[pt][cent] = (TProfile*) file_flow_invM_Input->Get(Form("Hist_v2_reso_ptSetB%d_centSetA%d_pfx",pt,cent));
    }
  }
  // pt SetB, cent SetB
  // TH1D *mHist_Input_SE_InvM_ptSetB_centSetB[4][9];
  // TH1D *mHist_Input_ME_InvM_ptSetB_centSetB[4][9];
  // TProfile *mProfile_Input_v1_raw_ptSetB_centSetB[4][9];
  // TProfile *mProfile_Input_v1_reso_ptSetB_centSetB[4][9];
  // TProfile *mProfile_Input_v2_raw_ptSetB_centSetB[4][9];
  // TProfile *mProfile_Input_v2_reso_ptSetB_centSetB[4][9];

  // pt SetC, cent 0-60%, 0-80%
  // TH1D *mHist_Input_SE_InvM_ptSetC_centAll[10][2];
  // TH1D *mHist_Input_ME_InvM_ptSetC_centAll[10][2];
  // TProfile *mProfile_Input_v1_raw_ptSetC_centAll[10][2];
  // TProfile *mProfile_Input_v1_reso_ptSetC_centAll[10][2];
  // TProfile *mProfile_Input_v2_raw_ptSetC_centAll[10][2];
  // TProfile *mProfile_Input_v2_reso_ptSetC_centAll[10][2];
  // rap SetA, cent SetA
  // pT range cut [0.1,1.0], [1.0, 2.0], [0.1, 2.0]
  // TH1D *mHist_Input_SE_InvM_rapSetA_centSetA_pTRange[4][6][3];
  // TH1D *mHist_Input_ME_InvM_rapSetA_centSetA_pTRange[4][6][3];
  // TH2D *mHist_v1_raw_rapSetA_centSetA_pTRange[4][6][3];
  // TH2D *mHist_v1_reso_rapSetA_centSetA_pTRange[4][6][3];
  // TH2D *mHist_v2_raw_rapSetA_centSetA_pTRange[4][6][3];
  // TH2D *mHist_v2_reso_rapSetA_centSetA_pTRange[4][6][3];
  // TProfile *mProfile_Input_v1_raw_rapSetA_centSetA_pTRange[4][6][3];
  // TProfile *mProfile_Input_v1_reso_rapSetA_centSetA_pTRange[4][6][3];
  // TProfile *mProfile_Input_v2_raw_rapSetA_centSetA_pTRange[4][6][3];
  // TProfile *mProfile_Input_v2_reso_rapSetA_centSetA_pTRange[4][6][3];
  // rap SetA, cent SetB
  // ---------------------- Output files and plots -----------------------------
  out_file.Append(".phiflow.result.root");
  out_file.Prepend(Form("_var%d_iter%d_", sys_varN, sys_iterN));
  out_file.Prepend(sys_object[sys_cutN]);
  out_file.Prepend(outHead);
  // out_file.Prepend("./out_sys/out_sys_Crosscheck/");
  TFile *outputFile = new TFile(out_file,"recreate");
  // pt SetA, cent SetA
  TCanvas *canvas_InvM_ptSetA_centSetA = new TCanvas("canvas_InvM_ptSetA_centSetA","canvas_InvM_ptSetA_centSetA",1920,1080);
  TCanvas *canvas_v1_raw_ptSetA_centSetA = new TCanvas("canvas_v1_raw_ptSetA_centSetA","canvas_v1_raw_ptSetA_centSetA",1920,1080);
  TCanvas *canvas_v1_reso_ptSetA_centSetA = new TCanvas("canvas_v1_reso_ptSetA_centSetA","canvas_v1_reso_ptSetA_centSetA",1920,1080);
  TCanvas *canvas_v2_raw_ptSetA_centSetA = new TCanvas("canvas_v2_raw_ptSetA_centSetA","canvas_v2_raw_ptSetA_centSetA",1920,1080);
  TCanvas *canvas_v2_reso_ptSetA_centSetA = new TCanvas("canvas_v2_reso_ptSetA_centSetA","canvas_v2_reso_ptSetA_centSetA",1920,1080);
  canvas_InvM_ptSetA_centSetA->Divide(6,2);
  canvas_v1_raw_ptSetA_centSetA->Divide(6,2);
  canvas_v1_reso_ptSetA_centSetA->Divide(6,2);
  canvas_v2_raw_ptSetA_centSetA->Divide(6,2);
  canvas_v2_reso_ptSetA_centSetA->Divide(6,2);
  TCanvas *canvas_v1_raw_vs_pT_ptSetA_centSetA = new TCanvas("canvas_v1_raw_vs_pT_ptSetA_centSetA","canvas_v1_raw_vs_pT_ptSetA_centSetA",1920,1080);
  TCanvas *canvas_v1_reso_vs_pT_ptSetA_centSetA = new TCanvas("canvas_v1_reso_vs_pT_ptSetA_centSetA","canvas_v1_reso_vs_pT_ptSetA_centSetA",1920,1080);
  TCanvas *canvas_v2_raw_vs_pT_ptSetA_centSetA = new TCanvas("canvas_v2_raw_vs_pT_ptSetA_centSetA","canvas_v2_raw_vs_pT_ptSetA_centSetA",1920,1080);
  TCanvas *canvas_v2_reso_vs_pT_ptSetA_centSetA = new TCanvas("canvas_v2_reso_vs_pT_ptSetA_centSetA","canvas_v2_reso_vs_pT_ptSetA_centSetA",1920,1080);
  canvas_v1_raw_vs_pT_ptSetA_centSetA->Divide(6);
  canvas_v1_reso_vs_pT_ptSetA_centSetA->Divide(6);
  canvas_v2_raw_vs_pT_ptSetA_centSetA->Divide(6);
  canvas_v2_reso_vs_pT_ptSetA_centSetA->Divide(6);
  TGraphErrors *mTGE_v1_raw_vs_pT_ptSetA_centSetA[6];
  TGraphErrors *mTGE_v1_reso_vs_pT_ptSetA_centSetA[6];
  TGraphErrors *mTGE_v2_raw_vs_pT_ptSetA_centSetA[6];
  TGraphErrors *mTGE_v2_reso_vs_pT_ptSetA_centSetA[6];
  Double_t d_v1_raw_vs_pT_ptSetA_centSetA[2];
  Double_t d_v1_reso_vs_pT_ptSetA_centSetA[2];
  Double_t d_v2_raw_vs_pT_ptSetA_centSetA[2];
  Double_t d_v2_reso_vs_pT_ptSetA_centSetA[2];
  // pt SetA, cent SetB
  TCanvas *canvas_InvM_ptSetA_centSetB = new TCanvas("canvas_InvM_ptSetA_centSetB","canvas_InvM_ptSetA_centSetB",1920,1080);
  TCanvas *canvas_v1_raw_ptSetA_centSetB = new TCanvas("canvas_v1_raw_ptSetA_centSetB","canvas_v1_raw_ptSetA_centSetB",1920,1080);
  TCanvas *canvas_v1_reso_ptSetA_centSetB = new TCanvas("canvas_v1_reso_ptSetA_centSetB","canvas_v1_reso_ptSetA_centSetB",1920,1080);
  TCanvas *canvas_v2_raw_ptSetA_centSetB = new TCanvas("canvas_v2_raw_ptSetA_centSetB","canvas_v2_raw_ptSetA_centSetB",1920,1080);
  TCanvas *canvas_v2_reso_ptSetA_centSetB = new TCanvas("canvas_v2_reso_ptSetA_centSetB","canvas_v2_reso_ptSetA_centSetB",1920,1080);
  canvas_InvM_ptSetA_centSetB->Divide(9,2);
  canvas_v1_raw_ptSetA_centSetB->Divide(9,2);
  canvas_v1_reso_ptSetA_centSetB->Divide(9,2);
  canvas_v2_raw_ptSetA_centSetB->Divide(9,2);
  canvas_v2_reso_ptSetA_centSetB->Divide(9,2);
  TCanvas *canvas_v1_raw_vs_pT_ptSetA_centSetB = new TCanvas("canvas_v1_raw_vs_pT_ptSetA_centSetB","canvas_v1_raw_vs_pT_ptSetA_centSetB",1920,1080);
  TCanvas *canvas_v1_reso_vs_pT_ptSetA_centSetB = new TCanvas("canvas_v1_reso_vs_pT_ptSetA_centSetB","canvas_v1_reso_vs_pT_ptSetA_centSetB",1920,1080);
  TCanvas *canvas_v2_raw_vs_pT_ptSetA_centSetB = new TCanvas("canvas_v2_raw_vs_pT_ptSetA_centSetB","canvas_v2_raw_vs_pT_ptSetA_centSetB",1920,1080);
  TCanvas *canvas_v2_reso_vs_pT_ptSetA_centSetB = new TCanvas("canvas_v2_reso_vs_pT_ptSetA_centSetB","canvas_v2_reso_vs_pT_ptSetA_centSetB",1920,1080);
  canvas_v1_raw_vs_pT_ptSetA_centSetB->Divide(9);
  canvas_v1_reso_vs_pT_ptSetA_centSetB->Divide(9);
  canvas_v2_raw_vs_pT_ptSetA_centSetB->Divide(9);
  canvas_v2_reso_vs_pT_ptSetA_centSetB->Divide(9);
  TGraphErrors *mTGE_v1_raw_vs_pT_ptSetA_centSetB[9];
  TGraphErrors *mTGE_v1_reso_vs_pT_ptSetA_centSetB[9];
  TGraphErrors *mTGE_v2_raw_vs_pT_ptSetA_centSetB[9];
  TGraphErrors *mTGE_v2_reso_vs_pT_ptSetA_centSetB[9];
  // pt SetB, cent SetA
  TCanvas *canvas_InvM_ptSetB_centSetA = new TCanvas("canvas_InvM_ptSetB_centSetA","canvas_InvM_ptSetB_centSetA",1920,1080);
  TCanvas *canvas_v1_raw_ptSetB_centSetA = new TCanvas("canvas_v1_raw_ptSetB_centSetA","canvas_v1_raw_ptSetB_centSetA",1920,1080);
  TCanvas *canvas_v1_reso_ptSetB_centSetA = new TCanvas("canvas_v1_reso_ptSetB_centSetA","canvas_v1_reso_ptSetB_centSetA",1920,1080);
  TCanvas *canvas_v2_raw_ptSetB_centSetA = new TCanvas("canvas_v2_raw_ptSetB_centSetA","canvas_v2_raw_ptSetB_centSetA",1920,1080);
  TCanvas *canvas_v2_reso_ptSetB_centSetA = new TCanvas("canvas_v2_reso_ptSetB_centSetA","canvas_v2_reso_ptSetB_centSetA",1920,1080);
  canvas_InvM_ptSetB_centSetA->Divide(6,4);
  canvas_v1_raw_ptSetB_centSetA->Divide(6,4);
  canvas_v1_reso_ptSetB_centSetA->Divide(6,4);
  canvas_v2_raw_ptSetB_centSetA->Divide(6,4);
  canvas_v2_reso_ptSetB_centSetA->Divide(6,4);
  TCanvas *canvas_v1_raw_vs_pT_ptSetB_centSetA = new TCanvas("canvas_v1_raw_vs_pT_ptSetB_centSetA","canvas_v1_raw_vs_pT_ptSetB_centSetA",1920,1080);
  TCanvas *canvas_v1_reso_vs_pT_ptSetB_centSetA = new TCanvas("canvas_v1_reso_vs_pT_ptSetB_centSetA","canvas_v1_reso_vs_pT_ptSetB_centSetA",1920,1080);
  TCanvas *canvas_v2_raw_vs_pT_ptSetB_centSetA = new TCanvas("canvas_v2_raw_vs_pT_ptSetB_centSetA","canvas_v2_raw_vs_pT_ptSetB_centSetA",1920,1080);
  TCanvas *canvas_v2_reso_vs_pT_ptSetB_centSetA = new TCanvas("canvas_v2_reso_vs_pT_ptSetB_centSetA","canvas_v2_reso_vs_pT_ptSetB_centSetA",1920,1080);
  canvas_v1_raw_vs_pT_ptSetB_centSetA->Divide(6);
  canvas_v1_reso_vs_pT_ptSetB_centSetA->Divide(6);
  canvas_v2_raw_vs_pT_ptSetB_centSetA->Divide(6);
  canvas_v2_reso_vs_pT_ptSetB_centSetA->Divide(6);
  TGraphErrors *mTGE_v1_raw_vs_pT_ptSetB_centSetA[6];
  TGraphErrors *mTGE_v1_reso_vs_pT_ptSetB_centSetA[6];
  TGraphErrors *mTGE_v2_raw_vs_pT_ptSetB_centSetA[6];
  TGraphErrors *mTGE_v2_reso_vs_pT_ptSetB_centSetA[6];
  // pt SetB, cent SetB
  /*
  TCanvas *canvas_InvM_ptSetB_centSetB = new TCanvas("canvas_InvM_ptSetB_centSetB","canvas_InvM_ptSetB_centSetB",1920,1080);
  TCanvas *canvas_v1_raw_ptSetB_centSetB = new TCanvas("canvas_v1_raw_ptSetB_centSetB","canvas_v1_raw_ptSetB_centSetB",1920,1080);
  TCanvas *canvas_v1_reso_ptSetB_centSetB = new TCanvas("canvas_v1_reso_ptSetB_centSetB","canvas_v1_reso_ptSetB_centSetB",1920,1080);
  TCanvas *canvas_v2_raw_ptSetB_centSetB = new TCanvas("canvas_v2_raw_ptSetB_centSetB","canvas_v2_raw_ptSetB_centSetB",1920,1080);
  TCanvas *canvas_v2_reso_ptSetB_centSetB = new TCanvas("canvas_v2_reso_ptSetB_centSetB","canvas_v2_reso_ptSetB_centSetB",1920,1080);
  canvas_InvM_ptSetB_centSetB->Divide(9,4);
  canvas_v1_raw_ptSetB_centSetB->Divide(9,4);
  canvas_v1_reso_ptSetB_centSetB->Divide(9,4);
  canvas_v2_raw_ptSetB_centSetB->Divide(9,4);
  canvas_v2_reso_ptSetB_centSetB->Divide(9,4);
  TCanvas *canvas_v1_raw_vs_pT_ptSetB_centSetB = new TCanvas("canvas_v1_raw_vs_pT_ptSetB_centSetB","canvas_v1_raw_vs_pT_ptSetB_centSetB",1920,1080);
  TCanvas *canvas_v1_reso_vs_pT_ptSetB_centSetB = new TCanvas("canvas_v1_reso_vs_pT_ptSetB_centSetB","canvas_v1_reso_vs_pT_ptSetB_centSetB",1920,1080);
  TCanvas *canvas_v2_raw_vs_pT_ptSetB_centSetB = new TCanvas("canvas_v2_raw_vs_pT_ptSetB_centSetB","canvas_v2_raw_vs_pT_ptSetB_centSetB",1920,1080);
  TCanvas *canvas_v2_reso_vs_pT_ptSetB_centSetB = new TCanvas("canvas_v2_reso_vs_pT_ptSetB_centSetB","canvas_v2_reso_vs_pT_ptSetB_centSetB",1920,1080);
  canvas_v1_raw_vs_pT_ptSetB_centSetB->Divide(9);
  canvas_v1_reso_vs_pT_ptSetB_centSetB->Divide(9);
  canvas_v2_raw_vs_pT_ptSetB_centSetB->Divide(9);
  canvas_v2_reso_vs_pT_ptSetB_centSetB->Divide(9);
  TGraphErrors *mTGE_v1_raw_vs_pT_ptSetB_centSetB[9];
  TGraphErrors *mTGE_v1_reso_vs_pT_ptSetB_centSetB[9];
  TGraphErrors *mTGE_v2_raw_vs_pT_ptSetB_centSetB[9];
  TGraphErrors *mTGE_v2_reso_vs_pT_ptSetB_centSetB[9];
  */
  // rap SetA, cent SetA

  // rap SetA, cent SetA, pTRange [0.1,1.0], [1,2], [0.1,2.0]
  // rap SetA, cent SetB
  // ======================== (2) Fit SE and ME InvM plots =====================
  // ----------------------- Normalization range -------------------------------
  Double_t a_d_int_range[4] ={
    0.99,
    1.014,
    1.026,
    1.09
  };
  // pt SetA, cent SetA
  for(int jkk=0; jkk<_n_jkk; jkk++){
    for(int pt=0; pt<2; pt++)
    {
      for(int cent=0; cent<6;cent++){
        canvas_InvM_ptSetA_centSetA->cd((cent+1)+6*pt);
        mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent]->GetYaxis()->SetRangeUser(-0.1*(Double_t)mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent]->GetMaximum(),1.1*(Double_t)mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent]->GetMaximum());
        mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
        // Get the bin of the Normalization range
        int a_iBin_range[4];
        for(int ijk = 0; ijk < 4; ijk++) a_iBin_range[ijk] =  mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent] -> FindFixBin(a_d_int_range[ijk]);
        //right bg Normalization
        Double_t d_r_area     = a_d_int_range[3]-a_d_int_range[2];
        Double_t d_r_same_int = mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent] -> Integral(a_iBin_range[2],a_iBin_range[3]);
        Double_t d_r_mx_int   = mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> Integral(a_iBin_range[2],a_iBin_range[3]);
        Double_t d_r_norm     = (d_r_mx_int != 0.0) ? d_r_same_int/d_r_mx_int : 1.0;
        cout<<" R: "<<d_r_area<<" : "<<d_r_same_int<<" : "<<d_r_mx_int<<endl;
        //left bg Normalization
        Double_t d_l_area =  a_d_int_range[1]-a_d_int_range[0];
        Double_t d_l_same_int = mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent] -> Integral(a_iBin_range[0],a_iBin_range[1]);
        Double_t d_l_mx_int   = mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> Integral(a_iBin_range[0],a_iBin_range[1]);
        Double_t d_l_norm     = (d_l_mx_int != 0.0) ? d_l_same_int/d_l_mx_int : 1.0;
        cout<<" L: "<<d_l_area<<" : "<<d_l_same_int<<" : "<<d_l_mx_int<<endl;
        cout<<" l: "<<d_l_norm<<" r: "<<d_r_norm<<endl;
        Double_t d_norm = ((d_r_norm/d_r_area) + (d_l_norm/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));
        cout<<" d_norm = "<<d_norm<<endl;

        // Normalize mixed event invariant mass
        mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent] -> SetMarkerStyle(1);
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> Sumw2();
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> Scale(d_norm);
        mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent]->Draw();
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> SetLineColor(kRed);
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> SetFillColor(kRed);
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> SetFillStyle(3002);
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]->Draw("HISTsames");
        // Substract normalized ME from SE to get Signal
        TH1D * HistSignal = (TH1D*) mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent] -> Clone("HistSignal");
        // HistSignal->SetLineColor(kRed);
        HistSignal -> Reset();
        HistSignal -> Sumw2();

        for(int ijk = 1; ijk < HistSignal->GetNbinsX()+1; ijk++)
        {
          Double_t d_center   = mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent] -> GetBinCenter(ijk);
          Double_t d_same     = mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent] -> GetBinContent(ijk);
          Double_t d_same_err = mHist_Input_SE_InvM_ptSetA_centSetA[jkk][pt][cent] -> GetBinError(ijk);
          Double_t d_mx       = mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> GetBinContent(ijk);
          Double_t d_mx_err   = mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> GetBinError(ijk);
          Double_t d_sig      = d_same - d_mx;
          Double_t d_sig_err  = sqrt(d_same_err*d_same_err+d_mx_err*d_mx_err);

          HistSignal -> SetBinContent(ijk,d_sig);
          HistSignal -> SetBinError(ijk,d_sig_err);
        }
        HistSignal->SetMarkerStyle(2);
        HistSignal->SetMarkerColor(kBlue);
        HistSignal->SetLineColor(kBlue);
        HistSignal -> Draw("HISTsamesE");
        gStyle->SetOptFit(1111);

        //Fit function
        //fit Signal with Gauss plus constant
        TFormula * GausPlus = new TFormula("GausPlus","gaus(0)+[3]");
        TF1 * tf1_Signal = new TF1("polygauss_single",GausPlus->GetExpFormula(),0.99,1.09);
        //fit to a simple gauss first to get seed
        TF1 * tf1_gauss = new TF1("tf1_gauss","gaus",0.9,1.1);
        HistSignal -> Fit(tf1_gauss,"0","R",1.01,1.03);
        // seeds
        Double_t d_seeds_p0    = tf1_gauss -> GetParameter(0);
        Double_t d_seeds_mean  = tf1_gauss -> GetParameter(1);
        Double_t d_seeds_sigma = tf1_gauss -> GetParameter(2);

        tf1_Signal -> SetParameter(1,d_seeds_mean);
        tf1_Signal -> SetParameter(2,d_seeds_sigma);
        tf1_Signal -> SetParLimits(1,d_seeds_mean-d_seeds_sigma,d_seeds_mean+d_seeds_sigma);
        tf1_Signal -> SetParLimits(2,0.66*d_seeds_sigma,1.5*d_seeds_sigma);
        tf1_Signal -> SetLineColor(kBlue);

        int FitStatus = HistSignal   -> Fit(tf1_Signal,"E+","R",0.99,1.09);
        tf1_Signal->Draw("same");
        canvas_InvM_ptSetA_centSetA->cd((cent+1)+6*pt)->Update();
        cout << "FitStatus= " << FitStatus << endl;
        // To count how many #phi mesons HistSignal has
        dParSig[0]    = tf1_Signal -> GetParameter(0);
        dParSig[1]    = tf1_Signal -> GetParameter(1);
        dParSig[2]    = tf1_Signal -> GetParameter(2);
        dParSig[3]    = tf1_Signal -> GetParameter(3);
        int iBin_3sigint_low = HistSignal -> FindFixBin(dParSig[1] - (3*dParSig[2]));
        int iBin_3sigint_hi  = HistSignal -> FindFixBin(dParSig[1] + (3*dParSig[2]));
        Double_t d_3sig_integral_error;
        Double_t d_3sig_integral = HistSignal -> IntegralAndError(iBin_3sigint_low,iBin_3sigint_hi,d_3sig_integral_error,"");
        TPaveText * ptext = new TPaveText(0.1,0.65,0.30,0.9,"NDCARC");
        ptext -> AddText(Form("Mean: %.4f",dParSig[1]));
        ptext -> AddText(Form("Sigma: %.4f",dParSig[2]));
        ptext -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral));
        ptext -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error));
        ptext -> Draw("same");
        // Fitting the background
        // TH1D * mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] = (TH1D*) mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> Clone("mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]");
        // mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]   -> SetTitle(Form("%s Background fit with Pol2",mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]->GetTitle() ));
        TF1 * tf1_Background = new TF1("tf1_pol","[0] + [1]*x + [2]*x**2",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> SetFillColor(kRed);
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent] -> SetFillStyle(3002);
        // mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]->SetMaximum(mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]->GetBinContent(mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]->GetMaximumBin())*1.4);
        // mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]->Draw();
        mHist_Input_ME_InvM_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_Background,"E+","R",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
        // sameEventInvM->Draw("same");
        tf1_Background->SetLineColor(kRed);
        tf1_Background->Draw("same");
        dParBg[0]=tf1_Background->GetParameter(0);
        dParBg[1]=tf1_Background->GetParameter(1);
        dParBg[2]=tf1_Background->GetParameter(2);
        TF1 * tf1_backgroundFlow = new TF1("tf1_backgroundFlow",BackgroundFitting,0.99,1.09,/*1*//*2*/3/*4*/);
        TF1 * tf1_totalFlow = new TF1("tf1_totalFlow",TotalFitting,0.99,1.09,/*2*//*3*/4/*5*/);

        canvas_v1_raw_ptSetA_centSetA->cd((cent+1)+6*pt);
        mProfile_Input_v1_raw_ptSetA_centSetA[jkk][pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
        mProfile_Input_v1_raw_ptSetA_centSetA[jkk][pt][cent]->Draw();
        mProfile_Input_v1_raw_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
        Double_t d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
        Double_t d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
        Double_t d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
        tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
        tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
        tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
        tf1_totalFlow -> SetLineColor(kBlue);
        mProfile_Input_v1_raw_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
        d_FLow_ptSetA_centSetA[jkk][0][0][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
        d_Flow_err_ptSetA_centSetA[jkk][0][0][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
        TPaveText * ptextFlow_v1_raw_ptSetA_centSetA = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
        ptextFlow_v1_raw_ptSetA_centSetA -> AddText(Form("v_{1}^{sig}: %.4f %c %.4f",d_FLow_ptSetA_centSetA[jkk][0][0][pt][cent],177,d_Flow_err_ptSetA_centSetA[jkk][0][0][pt][cent]));
        ptextFlow_v1_raw_ptSetA_centSetA -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
        ptextFlow_v1_raw_ptSetA_centSetA->Draw("same");

        canvas_v1_reso_ptSetA_centSetA->cd((cent+1)+6*pt);
        mProfile_Input_v1_reso_ptSetA_centSetA[jkk][pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
        mProfile_Input_v1_reso_ptSetA_centSetA[jkk][pt][cent]->Draw();
        mProfile_Input_v1_reso_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
        d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
        d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
        d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
        tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
        tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
        tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
        tf1_totalFlow -> SetLineColor(kBlue);
        mProfile_Input_v1_reso_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
        d_FLow_ptSetA_centSetA[jkk][0][1][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
        d_Flow_err_ptSetA_centSetA[jkk][0][1][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
        TPaveText * ptextFlow_v1_reso_ptSetA_centSetA = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
        ptextFlow_v1_reso_ptSetA_centSetA -> AddText(Form("v_{1}^{sig}: %.4f %c %.4f",d_FLow_ptSetA_centSetA[jkk][0][1][pt][cent],177,d_Flow_err_ptSetA_centSetA[jkk][0][1][pt][cent]));
        ptextFlow_v1_reso_ptSetA_centSetA -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
        ptextFlow_v1_reso_ptSetA_centSetA->Draw("same");

        canvas_v2_raw_ptSetA_centSetA->cd((cent+1)+6*pt);
        mProfile_Input_v2_raw_ptSetA_centSetA[jkk][pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
        mProfile_Input_v2_raw_ptSetA_centSetA[jkk][pt][cent]->Draw();
        mProfile_Input_v2_raw_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
        d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
        d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
        d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
        tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
        tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
        tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
        tf1_totalFlow -> SetLineColor(kBlue);
        mProfile_Input_v2_raw_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
        d_FLow_ptSetA_centSetA[jkk][1][0][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
        d_Flow_err_ptSetA_centSetA[jkk][1][0][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
        TPaveText * ptextFlow_v2_raw_ptSetA_centSetA = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
        ptextFlow_v2_raw_ptSetA_centSetA -> AddText(Form("v_{2}^{sig}: %.4f %c %.4f",d_FLow_ptSetA_centSetA[jkk][1][0][pt][cent],177,d_Flow_err_ptSetA_centSetA[jkk][1][0][pt][cent]));
        ptextFlow_v2_raw_ptSetA_centSetA -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
        ptextFlow_v2_raw_ptSetA_centSetA->Draw("same");

        canvas_v2_reso_ptSetA_centSetA->cd((cent+1)+6*pt);
        mProfile_Input_v2_reso_ptSetA_centSetA[jkk][pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
        mProfile_Input_v2_reso_ptSetA_centSetA[jkk][pt][cent]->Draw();
        mProfile_Input_v2_reso_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
        TFitResultPtr  bkg_fit_result = mProfile_Input_v2_reso_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_backgroundFlow,"S","R",0.99,1.09);
        d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
        d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
        d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
        Double_t integral_bkg = tf1_backgroundFlow->Integral(1.04,1.09);
        integral_bkg /= 0.05;
        Double_t integralErr_bkg = tf1_backgroundFlow->IntegralError(1.04,1.09,bkg_fit_result->GetParams(), bkg_fit_result->GetCovarianceMatrix().GetMatrixArray() );
        integralErr_bkg /= 0.05;
        tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
        tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
        tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
        tf1_totalFlow -> SetLineColor(kBlue);
        mProfile_Input_v2_reso_ptSetA_centSetA[jkk][pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
        d_FLow_ptSetA_centSetA[jkk][1][1][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
        d_Flow_err_ptSetA_centSetA[jkk][1][1][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
        TPaveText * ptextFlow_v2_reso_ptSetA_centSetA = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
        ptextFlow_v2_reso_ptSetA_centSetA -> AddText(Form("v_{2}^{sig}: %.4f %c %.4f",d_FLow_ptSetA_centSetA[jkk][1][1][pt][cent],177,d_Flow_err_ptSetA_centSetA[jkk][1][1][pt][cent]));
        ptextFlow_v2_reso_ptSetA_centSetA -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
        ptextFlow_v2_reso_ptSetA_centSetA -> AddText(Form("v_{2}^{bkg} : %.3f M_{inv}^{2} + %.3f M_{inv} + %.3f",d_V2_bg_p2,d_V2_bg_p1,d_V2_bg_p0));
        ptextFlow_v2_reso_ptSetA_centSetA -> AddText(Form("Integrated v_{2}^{bkg} (M_{inv}: [1.04, 1.09]) : %.4f %c %.4f ",integral_bkg,177,integralErr_bkg) );
        ptextFlow_v2_reso_ptSetA_centSetA->Draw("same");
      }
    }
  }
  const int n_ptSetA_centSetA = 2;
  for(int jkk=0; jkk<_n_jkk; jkk++){
    for(int cent=0; cent<6;cent++){
      TLine *l1_ptSetA_centSetA = new TLine(0.2,0,2.2,0);
      l1_ptSetA_centSetA->SetLineStyle(2);
      Double_t x[n_ptSetA_centSetA] = {0.9, 1.8};
      Double_t ex[n_ptSetA_centSetA] = {0.3, 0.6};
      Double_t y_v1_raw[n_ptSetA_centSetA] = {d_FLow_ptSetA_centSetA[0][0][0][0][cent], d_FLow_ptSetA_centSetA[0][0][0][1][cent]};
      Double_t ey_v1_raw[n_ptSetA_centSetA] = {d_Flow_err_ptSetA_centSetA[0][0][0][0][cent], d_Flow_err_ptSetA_centSetA[0][0][0][1][cent]};
      canvas_v1_raw_vs_pT_ptSetA_centSetA->cd(cent+1);
      mTGE_v1_raw_vs_pT_ptSetA_centSetA[cent] = new TGraphErrors(n_ptSetA_centSetA,x,y_v1_raw,ex,ey_v1_raw);
      mTGE_v1_raw_vs_pT_ptSetA_centSetA[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
      mTGE_v1_raw_vs_pT_ptSetA_centSetA[cent]->GetYaxis()->SetTitle("v_{1}");
      mTGE_v1_raw_vs_pT_ptSetA_centSetA[cent]->SetMarkerColor(4);
      mTGE_v1_raw_vs_pT_ptSetA_centSetA[cent]->SetMarkerStyle(24);
      mTGE_v1_raw_vs_pT_ptSetA_centSetA[cent]->Draw("AP");
      l1_ptSetA_centSetA->Draw("same");

      Double_t y_v1_reso[n_ptSetA_centSetA] = {d_FLow_ptSetA_centSetA[0][0][1][0][cent], d_FLow_ptSetA_centSetA[0][0][1][1][cent]};
      Double_t ey_v1_reso[n_ptSetA_centSetA] = {d_Flow_err_ptSetA_centSetA[0][0][1][0][cent], d_Flow_err_ptSetA_centSetA[0][0][1][1][cent]};
      canvas_v1_reso_vs_pT_ptSetA_centSetA->cd(cent+1);
      mTGE_v1_reso_vs_pT_ptSetA_centSetA[cent] = new TGraphErrors(n_ptSetA_centSetA,x,y_v1_reso,ex,ey_v1_reso);
      mTGE_v1_reso_vs_pT_ptSetA_centSetA[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
      mTGE_v1_reso_vs_pT_ptSetA_centSetA[cent]->GetYaxis()->SetTitle("v_{1}");
      mTGE_v1_reso_vs_pT_ptSetA_centSetA[cent]->SetMarkerColor(4);
      mTGE_v1_reso_vs_pT_ptSetA_centSetA[cent]->SetMarkerStyle(24);
      mTGE_v1_reso_vs_pT_ptSetA_centSetA[cent]->Draw("AP");
      l1_ptSetA_centSetA->Draw("same");

      Double_t y_v2_raw[n_ptSetA_centSetA] = {d_FLow_ptSetA_centSetA[0][1][0][0][cent], d_FLow_ptSetA_centSetA[0][1][0][1][cent]};
      Double_t ey_v2_raw[n_ptSetA_centSetA] = {d_Flow_err_ptSetA_centSetA[0][1][0][0][cent], d_Flow_err_ptSetA_centSetA[0][1][0][1][cent]};
      canvas_v2_raw_vs_pT_ptSetA_centSetA->cd(cent+1);
      mTGE_v2_raw_vs_pT_ptSetA_centSetA[cent] = new TGraphErrors(n_ptSetA_centSetA,x,y_v2_raw,ex,ey_v2_raw);
      mTGE_v2_raw_vs_pT_ptSetA_centSetA[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
      mTGE_v2_raw_vs_pT_ptSetA_centSetA[cent]->GetYaxis()->SetTitle("v_{2}^{raw}");
      mTGE_v2_raw_vs_pT_ptSetA_centSetA[cent]->SetMarkerColor(4);
      mTGE_v2_raw_vs_pT_ptSetA_centSetA[cent]->SetMarkerStyle(24);
      mTGE_v2_raw_vs_pT_ptSetA_centSetA[cent]->Draw("AP");
      l1_ptSetA_centSetA->Draw("same");

      Double_t y_v2_reso[n_ptSetA_centSetA] = {d_FLow_ptSetA_centSetA[0][1][1][0][cent], d_FLow_ptSetA_centSetA[0][1][1][1][cent]};
      Double_t ey_v2_reso[n_ptSetA_centSetA] = {d_Flow_err_ptSetA_centSetA[0][1][1][0][cent], d_Flow_err_ptSetA_centSetA[0][1][1][1][cent]};
      canvas_v2_reso_vs_pT_ptSetA_centSetA->cd(cent+1);
      mTGE_v2_reso_vs_pT_ptSetA_centSetA[cent] = new TGraphErrors(n_ptSetA_centSetA,x,y_v2_reso,ex,ey_v2_reso);
      mTGE_v2_reso_vs_pT_ptSetA_centSetA[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
      mTGE_v2_reso_vs_pT_ptSetA_centSetA[cent]->GetYaxis()->SetTitle("v_{2}");
      mTGE_v2_reso_vs_pT_ptSetA_centSetA[cent]->SetMarkerColor(4);
      mTGE_v2_reso_vs_pT_ptSetA_centSetA[cent]->SetMarkerStyle(24);
      mTGE_v2_reso_vs_pT_ptSetA_centSetA[cent]->Draw("AP");
      l1_ptSetA_centSetA->Draw("same");
    }
  }
  mTGE_v1_raw_vs_pT_ptSetA_centSetA[0]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[1]));
  mTGE_v1_raw_vs_pT_ptSetA_centSetA[1]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[1],centSetA[2]));
  mTGE_v1_raw_vs_pT_ptSetA_centSetA[2]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[2],centSetA[3]));
  mTGE_v1_raw_vs_pT_ptSetA_centSetA[3]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[2],centSetA[4]));
  mTGE_v1_raw_vs_pT_ptSetA_centSetA[4]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[3]));
  mTGE_v1_raw_vs_pT_ptSetA_centSetA[5]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[4]));

  mTGE_v1_reso_vs_pT_ptSetA_centSetA[0]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[0],centSetA[1]));
  mTGE_v1_reso_vs_pT_ptSetA_centSetA[1]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[1],centSetA[2]));
  mTGE_v1_reso_vs_pT_ptSetA_centSetA[2]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[2],centSetA[3]));
  mTGE_v1_reso_vs_pT_ptSetA_centSetA[3]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[2],centSetA[4]));
  mTGE_v1_reso_vs_pT_ptSetA_centSetA[4]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[0],centSetA[3]));
  mTGE_v1_reso_vs_pT_ptSetA_centSetA[5]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[0],centSetA[4]));

  mTGE_v2_raw_vs_pT_ptSetA_centSetA[0]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[1]));
  mTGE_v2_raw_vs_pT_ptSetA_centSetA[1]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[1],centSetA[2]));
  mTGE_v2_raw_vs_pT_ptSetA_centSetA[2]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[2],centSetA[3]));
  mTGE_v2_raw_vs_pT_ptSetA_centSetA[3]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[2],centSetA[4]));
  mTGE_v2_raw_vs_pT_ptSetA_centSetA[4]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[3]));
  mTGE_v2_raw_vs_pT_ptSetA_centSetA[5]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[4]));

  mTGE_v2_reso_vs_pT_ptSetA_centSetA[0]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[0],centSetA[1]));
  mTGE_v2_reso_vs_pT_ptSetA_centSetA[1]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[1],centSetA[2]));
  mTGE_v2_reso_vs_pT_ptSetA_centSetA[2]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[2],centSetA[3]));
  mTGE_v2_reso_vs_pT_ptSetA_centSetA[3]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[2],centSetA[4]));
  mTGE_v2_reso_vs_pT_ptSetA_centSetA[4]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[0],centSetA[3]));
  mTGE_v2_reso_vs_pT_ptSetA_centSetA[5]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[0],centSetA[4]));
  // pt SetB, cent SetA
  for(int pt=0; pt<4; pt++)
  {
    for(int cent=0; cent<6;cent++){
      canvas_InvM_ptSetB_centSetA->cd((cent+1)+6*pt);
      mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent]->GetYaxis()->SetRangeUser(-0.1*(Double_t)mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent]->GetMaximum(),1.1*(Double_t)mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent]->GetMaximum());
      mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      // Get the bin of the Normalization range
      int a_iBin_range[4];
      for(int ijk = 0; ijk < 4; ijk++) a_iBin_range[ijk] =  mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent] -> FindFixBin(a_d_int_range[ijk]);
      //right bg Normalization
      Double_t d_r_area     = a_d_int_range[3]-a_d_int_range[2];
      Double_t d_r_same_int = mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent] -> Integral(a_iBin_range[2],a_iBin_range[3]);
      Double_t d_r_mx_int   = mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> Integral(a_iBin_range[2],a_iBin_range[3]);
      Double_t d_r_norm     = (d_r_mx_int != 0.0) ? d_r_same_int/d_r_mx_int : 1.0;
      cout<<" R: "<<d_r_area<<" : "<<d_r_same_int<<" : "<<d_r_mx_int<<endl;
      //left bg Normalization
      Double_t d_l_area =  a_d_int_range[1]-a_d_int_range[0];
      Double_t d_l_same_int = mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent] -> Integral(a_iBin_range[0],a_iBin_range[1]);
      Double_t d_l_mx_int   = mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> Integral(a_iBin_range[0],a_iBin_range[1]);
      Double_t d_l_norm     = (d_l_mx_int != 0.0) ? d_l_same_int/d_l_mx_int : 1.0;
      cout<<" L: "<<d_l_area<<" : "<<d_l_same_int<<" : "<<d_l_mx_int<<endl;
      cout<<" l: "<<d_l_norm<<" r: "<<d_r_norm<<endl;
      Double_t d_norm = ((d_r_norm/d_r_area) + (d_l_norm/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));
      cout<<" d_norm = "<<d_norm<<endl;

      // Normalize mixed event invariant mass
      mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent] -> SetMarkerStyle(1);
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> Sumw2();
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> Scale(d_norm);
      mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent]->Draw();
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> SetLineColor(kRed);
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> SetFillColor(kRed);
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> SetFillStyle(3002);
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]->Draw("HISTsames");
      // Substract normalized ME from SE to get Signal
      TH1D * HistSignal = (TH1D*) mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent] -> Clone("HistSignal");
      // HistSignal->SetLineColor(kRed);
      HistSignal -> Reset();
      HistSignal -> Sumw2();

      for(int ijk = 1; ijk < HistSignal->GetNbinsX()+1; ijk++)
      {
        Double_t d_center   = mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent] -> GetBinCenter(ijk);
        Double_t d_same     = mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent] -> GetBinContent(ijk);
        Double_t d_same_err = mHist_Input_SE_InvM_ptSetB_centSetA[pt][cent] -> GetBinError(ijk);
        Double_t d_mx       = mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> GetBinContent(ijk);
        Double_t d_mx_err   = mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> GetBinError(ijk);
        Double_t d_sig      = d_same - d_mx;
        Double_t d_sig_err  = sqrt(d_same_err*d_same_err+d_mx_err*d_mx_err);

        HistSignal -> SetBinContent(ijk,d_sig);
        HistSignal -> SetBinError(ijk,d_sig_err);
      }
      HistSignal->SetMarkerStyle(2);
      HistSignal->SetMarkerColor(kBlue);
      HistSignal->SetLineColor(kBlue);
      HistSignal -> Draw("HISTsamesE");
      gStyle->SetOptFit(1111);

      //Fit function
      //fit Signal with Gauss plus constant
      TFormula * GausPlus = new TFormula("GausPlus","gaus(0)+[3]");
      TF1 * tf1_Signal = new TF1("polygauss_single",GausPlus->GetExpFormula(),0.99,1.09);
      //fit to a simple gauss first to get seed
      TF1 * tf1_gauss = new TF1("tf1_gauss","gaus",0.9,1.1);
      HistSignal -> Fit(tf1_gauss,"0","R",1.01,1.03);
      // seeds
      Double_t d_seeds_p0    = tf1_gauss -> GetParameter(0);
      Double_t d_seeds_mean  = tf1_gauss -> GetParameter(1);
      Double_t d_seeds_sigma = tf1_gauss -> GetParameter(2);

      tf1_Signal -> SetParameter(1,d_seeds_mean);
      tf1_Signal -> SetParameter(2,d_seeds_sigma);
      tf1_Signal -> SetParLimits(1,d_seeds_mean-d_seeds_sigma,d_seeds_mean+d_seeds_sigma);
      tf1_Signal -> SetParLimits(2,0.66*d_seeds_sigma,1.5*d_seeds_sigma);
      tf1_Signal -> SetLineColor(kBlue);

      int FitStatus = HistSignal   -> Fit(tf1_Signal,"E+","R",0.99,1.09);
      tf1_Signal->Draw("same");
      canvas_InvM_ptSetB_centSetA->cd((cent+1)+6*pt)->Update();
      cout << "FitStatus= " << FitStatus << endl;
      // To count how many #phi mesons HistSignal has
      dParSig[0]    = tf1_Signal -> GetParameter(0);
      dParSig[1]    = tf1_Signal -> GetParameter(1);
      dParSig[2]    = tf1_Signal -> GetParameter(2);
      dParSig[3]    = tf1_Signal -> GetParameter(3);
      int iBin_3sigint_low = HistSignal -> FindFixBin(dParSig[1] - (3*dParSig[2]));
      int iBin_3sigint_hi  = HistSignal -> FindFixBin(dParSig[1] + (3*dParSig[2]));
      Double_t d_3sig_integral_error;
      Double_t d_3sig_integral = HistSignal -> IntegralAndError(iBin_3sigint_low,iBin_3sigint_hi,d_3sig_integral_error,"");
      TPaveText * ptext = new TPaveText(0.1,0.65,0.30,0.9,"NDCARC");
      ptext -> AddText(Form("Mean: %.4f",dParSig[1]));
      ptext -> AddText(Form("Sigma: %.4f",dParSig[2]));
      ptext -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral));
      ptext -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error));
      ptext -> Draw("same");
      // Fitting the background
      // TH1D * mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] = (TH1D*) mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> Clone("mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]");
      // mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]   -> SetTitle(Form("%s Background fit with Pol2",mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]->GetTitle() ));
      TF1 * tf1_Background = new TF1("tf1_pol","[0] + [1]*x + [2]*x**2",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> SetFillColor(kRed);
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent] -> SetFillStyle(3002);
      // mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]->SetMaximum(mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]->GetBinContent(mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]->GetMaximumBin())*1.4);
      // mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]->Draw();
      mHist_Input_ME_InvM_ptSetB_centSetA[pt][cent]->Fit(tf1_Background,"E+","R",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
      // sameEventInvM->Draw("same");
      tf1_Background->SetLineColor(kRed);
      tf1_Background->Draw("same");
      dParBg[0]=tf1_Background->GetParameter(0);
      dParBg[1]=tf1_Background->GetParameter(1);
      dParBg[2]=tf1_Background->GetParameter(2);
      TF1 * tf1_backgroundFlow = new TF1("tf1_backgroundFlow",BackgroundFitting,0.99,1.09,/*1*//*2*/3/*4*/);
      TF1 * tf1_totalFlow = new TF1("tf1_totalFlow",TotalFitting,0.99,1.09,/*2*//*3*/4/*5*/);

      canvas_v1_raw_ptSetB_centSetA->cd((cent+1)+6*pt);
      mProfile_Input_v1_raw_ptSetB_centSetA[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mProfile_Input_v1_raw_ptSetB_centSetA[pt][cent]->Draw();
      mProfile_Input_v1_raw_ptSetB_centSetA[pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
      Double_t d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
      Double_t d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
      Double_t d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
      tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
      tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
      tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
      tf1_totalFlow -> SetLineColor(kBlue);
      mProfile_Input_v1_raw_ptSetB_centSetA[pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
      d_FLow_ptSetB_centSetA[0][0][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
      d_Flow_err_ptSetB_centSetA[0][0][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
      TPaveText * ptextFlow_v1_raw_ptSetB_centSetA = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
      ptextFlow_v1_raw_ptSetB_centSetA -> AddText(Form("v_{1}^{sig}: %.4f %c %.4f",d_FLow_ptSetB_centSetA[0][0][pt][cent],177,d_Flow_err_ptSetB_centSetA[0][0][pt][cent]));
      ptextFlow_v1_raw_ptSetB_centSetA -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
      ptextFlow_v1_raw_ptSetB_centSetA->Draw("same");

      canvas_v1_reso_ptSetB_centSetA->cd((cent+1)+6*pt);
      mProfile_Input_v1_reso_ptSetB_centSetA[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mProfile_Input_v1_reso_ptSetB_centSetA[pt][cent]->Draw();
      mProfile_Input_v1_reso_ptSetB_centSetA[pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
      d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
      d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
      d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
      tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
      tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
      tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
      tf1_totalFlow -> SetLineColor(kBlue);
      mProfile_Input_v1_reso_ptSetB_centSetA[pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
      d_FLow_ptSetB_centSetA[0][1][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
      d_Flow_err_ptSetB_centSetA[0][1][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
      TPaveText * ptextFlow_v1_reso_ptSetB_centSetA = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
      ptextFlow_v1_reso_ptSetB_centSetA -> AddText(Form("v_{1}^{sig}: %.4f %c %.4f",d_FLow_ptSetB_centSetA[0][1][pt][cent],177,d_Flow_err_ptSetB_centSetA[0][1][pt][cent]));
      ptextFlow_v1_reso_ptSetB_centSetA -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
      ptextFlow_v1_reso_ptSetB_centSetA->Draw("same");

      canvas_v2_raw_ptSetB_centSetA->cd((cent+1)+6*pt);
      mProfile_Input_v2_raw_ptSetB_centSetA[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mProfile_Input_v2_raw_ptSetB_centSetA[pt][cent]->Draw();
      mProfile_Input_v2_raw_ptSetB_centSetA[pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
      d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
      d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
      d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
      tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
      tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
      tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
      tf1_totalFlow -> SetLineColor(kBlue);
      mProfile_Input_v2_raw_ptSetB_centSetA[pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
      d_FLow_ptSetB_centSetA[1][0][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
      d_Flow_err_ptSetB_centSetA[1][0][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
      TPaveText * ptextFlow_v2_raw_ptSetB_centSetA = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
      ptextFlow_v2_raw_ptSetB_centSetA -> AddText(Form("v_{2}^{sig}: %.4f %c %.4f",d_FLow_ptSetB_centSetA[1][0][pt][cent],177,d_Flow_err_ptSetB_centSetA[1][0][pt][cent]));
      ptextFlow_v2_raw_ptSetB_centSetA -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
      ptextFlow_v2_raw_ptSetB_centSetA->Draw("same");

      canvas_v2_reso_ptSetB_centSetA->cd((cent+1)+6*pt);
      mProfile_Input_v2_reso_ptSetB_centSetA[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mProfile_Input_v2_reso_ptSetB_centSetA[pt][cent]->Draw();
      mProfile_Input_v2_reso_ptSetB_centSetA[pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
      d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
      d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
      d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
      tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
      tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
      tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
      tf1_totalFlow -> SetLineColor(kBlue);
      mProfile_Input_v2_reso_ptSetB_centSetA[pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
      d_FLow_ptSetB_centSetA[1][1][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
      d_Flow_err_ptSetB_centSetA[1][1][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
      TPaveText * ptextFlow_v2_reso_ptSetB_centSetA = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
      ptextFlow_v2_reso_ptSetB_centSetA -> AddText(Form("v_{2}^{sig}: %.4f %c %.4f",d_FLow_ptSetB_centSetA[1][1][pt][cent],177,d_Flow_err_ptSetB_centSetA[1][1][pt][cent]));
      ptextFlow_v2_reso_ptSetB_centSetA -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
      ptextFlow_v2_reso_ptSetB_centSetA->Draw("same");
    }
  }
  const int n_ptSetB_centSetA = 4;
  for(int cent=0; cent<6;cent++){
    TLine *l1_ptSetB_centSetA = new TLine(0.2,0,2.2,0);
    l1_ptSetB_centSetA->SetLineStyle(2);
    Double_t x[n_ptSetB_centSetA] = {0.55, 0.85, 1.2, 1.7};
    Double_t ex[n_ptSetB_centSetA] = {0.15, 0.15, 0.2, 0.3};
    Double_t y_v1_raw[n_ptSetB_centSetA] = {d_FLow_ptSetB_centSetA[0][0][0][cent], d_FLow_ptSetB_centSetA[0][0][1][cent],
      d_FLow_ptSetB_centSetA[0][0][2][cent], d_FLow_ptSetB_centSetA[0][0][3][cent]};
    Double_t ey_v1_raw[n_ptSetB_centSetA] = {d_Flow_err_ptSetB_centSetA[0][0][0][cent], d_Flow_err_ptSetB_centSetA[0][0][1][cent],
      d_Flow_err_ptSetB_centSetA[0][0][2][cent], d_Flow_err_ptSetB_centSetA[0][0][3][cent]};
    canvas_v1_raw_vs_pT_ptSetB_centSetA->cd(cent+1);
    mTGE_v1_raw_vs_pT_ptSetB_centSetA[cent] = new TGraphErrors(n_ptSetB_centSetA,x,y_v1_raw,ex,ey_v1_raw);
    mTGE_v1_raw_vs_pT_ptSetB_centSetA[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
    mTGE_v1_raw_vs_pT_ptSetB_centSetA[cent]->GetYaxis()->SetTitle("v_{1}");
    mTGE_v1_raw_vs_pT_ptSetB_centSetA[cent]->SetMarkerColor(4);
    mTGE_v1_raw_vs_pT_ptSetB_centSetA[cent]->SetMarkerStyle(24);
    mTGE_v1_raw_vs_pT_ptSetB_centSetA[cent]->Draw("AP");
    l1_ptSetB_centSetA->Draw("same");

    Double_t y_v1_reso[n_ptSetB_centSetA] = {d_FLow_ptSetB_centSetA[0][1][0][cent], d_FLow_ptSetB_centSetA[0][1][1][cent],
      d_FLow_ptSetB_centSetA[0][1][2][cent], d_FLow_ptSetB_centSetA[0][1][3][cent]};
    Double_t ey_v1_reso[n_ptSetB_centSetA] = {d_Flow_err_ptSetB_centSetA[0][1][0][cent], d_Flow_err_ptSetB_centSetA[0][1][1][cent],
      d_Flow_err_ptSetB_centSetA[0][1][2][cent], d_Flow_err_ptSetB_centSetA[0][1][3][cent]};
    canvas_v1_reso_vs_pT_ptSetB_centSetA->cd(cent+1);
    mTGE_v1_reso_vs_pT_ptSetB_centSetA[cent] = new TGraphErrors(n_ptSetB_centSetA,x,y_v1_reso,ex,ey_v1_reso);
    mTGE_v1_reso_vs_pT_ptSetB_centSetA[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
    mTGE_v1_reso_vs_pT_ptSetB_centSetA[cent]->GetYaxis()->SetTitle("v_{1}");
    mTGE_v1_reso_vs_pT_ptSetB_centSetA[cent]->SetMarkerColor(4);
    mTGE_v1_reso_vs_pT_ptSetB_centSetA[cent]->SetMarkerStyle(24);
    mTGE_v1_reso_vs_pT_ptSetB_centSetA[cent]->Draw("AP");
    l1_ptSetB_centSetA->Draw("same");

    Double_t y_v2_raw[n_ptSetB_centSetA] = {d_FLow_ptSetB_centSetA[1][0][0][cent], d_FLow_ptSetB_centSetA[1][0][1][cent],
      d_FLow_ptSetB_centSetA[1][0][2][cent], d_FLow_ptSetB_centSetA[1][0][3][cent]};
    Double_t ey_v2_raw[n_ptSetB_centSetA] = {d_Flow_err_ptSetB_centSetA[1][0][0][cent], d_Flow_err_ptSetB_centSetA[1][0][1][cent],
      d_Flow_err_ptSetB_centSetA[1][0][2][cent], d_Flow_err_ptSetB_centSetA[1][0][3][cent]};
    canvas_v2_raw_vs_pT_ptSetB_centSetA->cd(cent+1);
    mTGE_v2_raw_vs_pT_ptSetB_centSetA[cent] = new TGraphErrors(n_ptSetB_centSetA,x,y_v2_raw,ex,ey_v2_raw);
    mTGE_v2_raw_vs_pT_ptSetB_centSetA[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
    mTGE_v2_raw_vs_pT_ptSetB_centSetA[cent]->GetYaxis()->SetTitle("v_{2}^{raw}");
    mTGE_v2_raw_vs_pT_ptSetB_centSetA[cent]->SetMarkerColor(4);
    mTGE_v2_raw_vs_pT_ptSetB_centSetA[cent]->SetMarkerStyle(24);
    mTGE_v2_raw_vs_pT_ptSetB_centSetA[cent]->Draw("AP");
    l1_ptSetB_centSetA->Draw("same");

    Double_t y_v2_reso[n_ptSetB_centSetA] = {d_FLow_ptSetB_centSetA[1][1][0][cent], d_FLow_ptSetB_centSetA[1][1][1][cent],
      d_FLow_ptSetB_centSetA[1][1][2][cent], d_FLow_ptSetB_centSetA[1][1][3][cent]};
    Double_t ey_v2_reso[n_ptSetB_centSetA] = {d_Flow_err_ptSetB_centSetA[1][1][0][cent], d_Flow_err_ptSetB_centSetA[1][1][1][cent],
      d_Flow_err_ptSetB_centSetA[1][1][2][cent], d_Flow_err_ptSetB_centSetA[1][1][3][cent]};
    canvas_v2_reso_vs_pT_ptSetB_centSetA->cd(cent+1);
    mTGE_v2_reso_vs_pT_ptSetB_centSetA[cent] = new TGraphErrors(n_ptSetB_centSetA,x,y_v2_reso,ex,ey_v2_reso);
    mTGE_v2_reso_vs_pT_ptSetB_centSetA[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
    mTGE_v2_reso_vs_pT_ptSetB_centSetA[cent]->GetYaxis()->SetTitle("v_{2}");
    mTGE_v2_reso_vs_pT_ptSetB_centSetA[cent]->SetMarkerColor(4);
    mTGE_v2_reso_vs_pT_ptSetB_centSetA[cent]->SetMarkerStyle(24);
    mTGE_v2_reso_vs_pT_ptSetB_centSetA[cent]->Draw("AP");
    l1_ptSetB_centSetA->Draw("same");
  }
  mTGE_v1_raw_vs_pT_ptSetB_centSetA[0]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[1]));
  mTGE_v1_raw_vs_pT_ptSetB_centSetA[1]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[1],centSetA[2]));
  mTGE_v1_raw_vs_pT_ptSetB_centSetA[2]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[2],centSetA[3]));
  mTGE_v1_raw_vs_pT_ptSetB_centSetA[3]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[2],centSetA[4]));
  mTGE_v1_raw_vs_pT_ptSetB_centSetA[4]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[3]));
  mTGE_v1_raw_vs_pT_ptSetB_centSetA[5]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[4]));

  mTGE_v1_reso_vs_pT_ptSetB_centSetA[0]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[0],centSetA[1]));
  mTGE_v1_reso_vs_pT_ptSetB_centSetA[1]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[1],centSetA[2]));
  mTGE_v1_reso_vs_pT_ptSetB_centSetA[2]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[2],centSetA[3]));
  mTGE_v1_reso_vs_pT_ptSetB_centSetA[3]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[2],centSetA[4]));
  mTGE_v1_reso_vs_pT_ptSetB_centSetA[4]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[0],centSetA[3]));
  mTGE_v1_reso_vs_pT_ptSetB_centSetA[5]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetA[0],centSetA[4]));

  mTGE_v2_raw_vs_pT_ptSetB_centSetA[0]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[1]));
  mTGE_v2_raw_vs_pT_ptSetB_centSetA[1]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[1],centSetA[2]));
  mTGE_v2_raw_vs_pT_ptSetB_centSetA[2]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[2],centSetA[3]));
  mTGE_v2_raw_vs_pT_ptSetB_centSetA[3]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[2],centSetA[4]));
  mTGE_v2_raw_vs_pT_ptSetB_centSetA[4]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[3]));
  mTGE_v2_raw_vs_pT_ptSetB_centSetA[5]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetA[0],centSetA[4]));

  mTGE_v2_reso_vs_pT_ptSetB_centSetA[0]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[0],centSetA[1]));
  mTGE_v2_reso_vs_pT_ptSetB_centSetA[1]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[1],centSetA[2]));
  mTGE_v2_reso_vs_pT_ptSetB_centSetA[2]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[2],centSetA[3]));
  mTGE_v2_reso_vs_pT_ptSetB_centSetA[3]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[2],centSetA[4]));
  mTGE_v2_reso_vs_pT_ptSetB_centSetA[4]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[0],centSetA[3]));
  mTGE_v2_reso_vs_pT_ptSetB_centSetA[5]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetA[0],centSetA[4]));
  // pt SetA, cent SetB
  for(int pt=0; pt<2; pt++)
  {
    for(int cent=0; cent<9;cent++){
      canvas_InvM_ptSetA_centSetB->cd((cent+1)+9*pt);
      mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent]->GetYaxis()->SetRangeUser(-0.1*(Double_t)mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent]->GetMaximum(),1.1*(Double_t)mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent]->GetMaximum());
      mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      // Get the bin of the Normalization range
      int a_iBin_range[4];
      for(int ijk = 0; ijk < 4; ijk++) a_iBin_range[ijk] =  mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent] -> FindFixBin(a_d_int_range[ijk]);
      //right bg Normalization
      Double_t d_r_area     = a_d_int_range[3]-a_d_int_range[2];
      Double_t d_r_same_int = mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent] -> Integral(a_iBin_range[2],a_iBin_range[3]);
      Double_t d_r_mx_int   = mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> Integral(a_iBin_range[2],a_iBin_range[3]);
      Double_t d_r_norm     = (d_r_mx_int != 0.0) ? d_r_same_int/d_r_mx_int : 1.0;
      cout<<" R: "<<d_r_area<<" : "<<d_r_same_int<<" : "<<d_r_mx_int<<endl;
      //left bg Normalization
      Double_t d_l_area =  a_d_int_range[1]-a_d_int_range[0];
      Double_t d_l_same_int = mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent] -> Integral(a_iBin_range[0],a_iBin_range[1]);
      Double_t d_l_mx_int   = mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> Integral(a_iBin_range[0],a_iBin_range[1]);
      Double_t d_l_norm     = (d_l_mx_int != 0.0) ? d_l_same_int/d_l_mx_int : 1.0;
      cout<<" L: "<<d_l_area<<" : "<<d_l_same_int<<" : "<<d_l_mx_int<<endl;
      cout<<" l: "<<d_l_norm<<" r: "<<d_r_norm<<endl;
      Double_t d_norm = ((d_r_norm/d_r_area) + (d_l_norm/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));
      cout<<" d_norm = "<<d_norm<<endl;

      // Normalize mixed event invariant mass
      mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent] -> SetMarkerStyle(1);
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> Sumw2();
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> Scale(d_norm);
      mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent]->Draw();
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> SetLineColor(kRed);
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> SetFillColor(kRed);
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> SetFillStyle(3002);
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]->Draw("HISTsames");
      // Substract normalized ME from SE to get Signal
      TH1D * HistSignal = (TH1D*) mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent] -> Clone("HistSignal");
      // HistSignal->SetLineColor(kRed);
      HistSignal -> Reset();
      HistSignal -> Sumw2();

      for(int ijk = 1; ijk < HistSignal->GetNbinsX()+1; ijk++)
      {
        Double_t d_center   = mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent] -> GetBinCenter(ijk);
        Double_t d_same     = mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent] -> GetBinContent(ijk);
        Double_t d_same_err = mHist_Input_SE_InvM_ptSetA_centSetB[pt][cent] -> GetBinError(ijk);
        Double_t d_mx       = mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> GetBinContent(ijk);
        Double_t d_mx_err   = mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> GetBinError(ijk);
        Double_t d_sig      = d_same - d_mx;
        Double_t d_sig_err  = sqrt(d_same_err*d_same_err+d_mx_err*d_mx_err);

        HistSignal -> SetBinContent(ijk,d_sig);
        HistSignal -> SetBinError(ijk,d_sig_err);
      }
      HistSignal->SetMarkerStyle(2);
      HistSignal->SetMarkerColor(kBlue);
      HistSignal->SetLineColor(kBlue);
      HistSignal -> Draw("HISTsamesE");
      gStyle->SetOptFit(1111);

      //Fit function
      //fit Signal with Gauss plus constant
      TFormula * GausPlus = new TFormula("GausPlus","gaus(0)+[3]");
      TF1 * tf1_Signal = new TF1("polygauss_single",GausPlus->GetExpFormula(),0.99,1.09);
      //fit to a simple gauss first to get seed
      TF1 * tf1_gauss = new TF1("tf1_gauss","gaus",0.9,1.1);
      HistSignal -> Fit(tf1_gauss,"0","R",1.01,1.03);
      // seeds
      Double_t d_seeds_p0    = tf1_gauss -> GetParameter(0);
      Double_t d_seeds_mean  = tf1_gauss -> GetParameter(1);
      Double_t d_seeds_sigma = tf1_gauss -> GetParameter(2);

      tf1_Signal -> SetParameter(1,d_seeds_mean);
      tf1_Signal -> SetParameter(2,d_seeds_sigma);
      tf1_Signal -> SetParLimits(1,d_seeds_mean-d_seeds_sigma,d_seeds_mean+d_seeds_sigma);
      tf1_Signal -> SetParLimits(2,0.66*d_seeds_sigma,1.5*d_seeds_sigma);
      tf1_Signal -> SetLineColor(kBlue);

      int FitStatus = HistSignal   -> Fit(tf1_Signal,"E+","R",0.99,1.09);
      tf1_Signal->Draw("same");
      canvas_InvM_ptSetA_centSetB->cd((cent+1)+9*pt)->Update();
      cout << "FitStatus= " << FitStatus << endl;
      // To count how many #phi mesons HistSignal has
      dParSig[0]    = tf1_Signal -> GetParameter(0);
      dParSig[1]    = tf1_Signal -> GetParameter(1);
      dParSig[2]    = tf1_Signal -> GetParameter(2);
      dParSig[3]    = tf1_Signal -> GetParameter(3);
      int iBin_3sigint_low = HistSignal -> FindFixBin(dParSig[1] - (3*dParSig[2]));
      int iBin_3sigint_hi  = HistSignal -> FindFixBin(dParSig[1] + (3*dParSig[2]));
      Double_t d_3sig_integral_error;
      Double_t d_3sig_integral = HistSignal -> IntegralAndError(iBin_3sigint_low,iBin_3sigint_hi,d_3sig_integral_error,"");
      TPaveText * ptext = new TPaveText(0.1,0.65,0.30,0.9,"NDCARC");
      ptext -> AddText(Form("Mean: %.4f",dParSig[1]));
      ptext -> AddText(Form("Sigma: %.4f",dParSig[2]));
      ptext -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral));
      ptext -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error));
      ptext -> Draw("same");
      // Fitting the background
      // TH1D * mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] = (TH1D*) mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> Clone("mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]");
      // mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]   -> SetTitle(Form("%s Background fit with Pol2",mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]->GetTitle() ));
      TF1 * tf1_Background = new TF1("tf1_pol","[0] + [1]*x + [2]*x**2",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> SetFillColor(kRed);
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent] -> SetFillStyle(3002);
      // mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]->SetMaximum(mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]->GetBinContent(mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]->GetMaximumBin())*1.4);
      // mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]->Draw();
      mHist_Input_ME_InvM_ptSetA_centSetB[pt][cent]->Fit(tf1_Background,"E+","R",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
      // sameEventInvM->Draw("same");
      tf1_Background->SetLineColor(kRed);
      tf1_Background->Draw("same");
      dParBg[0]=tf1_Background->GetParameter(0);
      dParBg[1]=tf1_Background->GetParameter(1);
      dParBg[2]=tf1_Background->GetParameter(2);
      TF1 * tf1_backgroundFlow = new TF1("tf1_backgroundFlow",BackgroundFitting,0.99,1.09,/*1*//*2*/3/*4*/);
      TF1 * tf1_totalFlow = new TF1("tf1_totalFlow",TotalFitting,0.99,1.09,/*2*//*3*/4/*5*/);

      canvas_v1_raw_ptSetA_centSetB->cd((cent+1)+9*pt);
      mProfile_Input_v1_raw_ptSetA_centSetB[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mProfile_Input_v1_raw_ptSetA_centSetB[pt][cent]->Draw();
      mProfile_Input_v1_raw_ptSetA_centSetB[pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
      Double_t d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
      Double_t d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
      Double_t d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
      tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
      tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
      tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
      tf1_totalFlow -> SetLineColor(kBlue);
      mProfile_Input_v1_raw_ptSetA_centSetB[pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
      d_FLow_ptSetA_centSetB[0][0][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
      d_Flow_err_ptSetA_centSetB[0][0][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
      TPaveText * ptextFlow_v1_raw_ptSetA_centSetB = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
      ptextFlow_v1_raw_ptSetA_centSetB -> AddText(Form("v_{1}^{sig}: %.4f %c %.4f",d_FLow_ptSetA_centSetB[0][0][pt][cent],177,d_Flow_err_ptSetA_centSetB[0][0][pt][cent]));
      ptextFlow_v1_raw_ptSetA_centSetB -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
      ptextFlow_v1_raw_ptSetA_centSetB->Draw("same");

      canvas_v1_reso_ptSetA_centSetB->cd((cent+1)+9*pt);
      mProfile_Input_v1_reso_ptSetA_centSetB[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mProfile_Input_v1_reso_ptSetA_centSetB[pt][cent]->Draw();
      mProfile_Input_v1_reso_ptSetA_centSetB[pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
      d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
      d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
      d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
      tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
      tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
      tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
      tf1_totalFlow -> SetLineColor(kBlue);
      mProfile_Input_v1_reso_ptSetA_centSetB[pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
      d_FLow_ptSetA_centSetB[0][1][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
      d_Flow_err_ptSetA_centSetB[0][1][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
      TPaveText * ptextFlow_v1_reso_ptSetA_centSetB = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
      ptextFlow_v1_reso_ptSetA_centSetB -> AddText(Form("v_{1}^{sig}: %.4f %c %.4f",d_FLow_ptSetA_centSetB[0][1][pt][cent],177,d_Flow_err_ptSetA_centSetB[0][1][pt][cent]));
      ptextFlow_v1_reso_ptSetA_centSetB -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
      ptextFlow_v1_reso_ptSetA_centSetB->Draw("same");

      canvas_v2_raw_ptSetA_centSetB->cd((cent+1)+9*pt);
      mProfile_Input_v2_raw_ptSetA_centSetB[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mProfile_Input_v2_raw_ptSetA_centSetB[pt][cent]->Draw();
      mProfile_Input_v2_raw_ptSetA_centSetB[pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
      d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
      d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
      d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
      tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
      tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
      tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
      tf1_totalFlow -> SetLineColor(kBlue);
      mProfile_Input_v2_raw_ptSetA_centSetB[pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
      d_FLow_ptSetA_centSetB[1][0][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
      d_Flow_err_ptSetA_centSetB[1][0][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
      TPaveText * ptextFlow_v2_raw_ptSetA_centSetB = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
      ptextFlow_v2_raw_ptSetA_centSetB -> AddText(Form("v_{2}^{sig}: %.4f %c %.4f",d_FLow_ptSetA_centSetB[1][0][pt][cent],177,d_Flow_err_ptSetA_centSetB[1][0][pt][cent]));
      ptextFlow_v2_raw_ptSetA_centSetB -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
      ptextFlow_v2_raw_ptSetA_centSetB->Draw("same");

      canvas_v2_reso_ptSetA_centSetB->cd((cent+1)+9*pt);
      mProfile_Input_v2_reso_ptSetA_centSetB[pt][cent]->GetXaxis()->SetRangeUser(0.99,1.09);
      mProfile_Input_v2_reso_ptSetA_centSetB[pt][cent]->Draw();
      mProfile_Input_v2_reso_ptSetA_centSetB[pt][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
      d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
      d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
      d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
      tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
      tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
      tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
      tf1_totalFlow -> SetLineColor(kBlue);
      mProfile_Input_v2_reso_ptSetA_centSetB[pt][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
      d_FLow_ptSetA_centSetB[1][1][pt][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
      d_Flow_err_ptSetA_centSetB[1][1][pt][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
      TPaveText * ptextFlow_v2_reso_ptSetA_centSetB = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
      ptextFlow_v2_reso_ptSetA_centSetB -> AddText(Form("v_{2}^{sig}: %.4f %c %.4f",d_FLow_ptSetA_centSetB[1][1][pt][cent],177,d_Flow_err_ptSetA_centSetB[1][1][pt][cent]));
      ptextFlow_v2_reso_ptSetA_centSetB -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
      ptextFlow_v2_reso_ptSetA_centSetB->Draw("same");
    }
  }
  const int n_ptSetA_centSetB = 2;
  for(int cent=0; cent<9;cent++){
    TLine *l1_ptSetA_centSetB = new TLine(0.2,0,2.2,0);
    l1_ptSetA_centSetB->SetLineStyle(2);
    Double_t x[n_ptSetA_centSetB] = {0.8, 1.6};
    Double_t ex[n_ptSetA_centSetB] = {0.4, 0.4};
    Double_t y_v1_raw[n_ptSetA_centSetB] = {d_FLow_ptSetA_centSetB[0][0][0][cent], d_FLow_ptSetA_centSetB[0][0][1][cent]};
    Double_t ey_v1_raw[n_ptSetA_centSetB] = {d_Flow_err_ptSetA_centSetB[0][0][0][cent], d_Flow_err_ptSetA_centSetB[0][0][1][cent]};
    canvas_v1_raw_vs_pT_ptSetA_centSetB->cd(cent+1);
    mTGE_v1_raw_vs_pT_ptSetA_centSetB[cent] = new TGraphErrors(n_ptSetA_centSetB,x,y_v1_raw,ex,ey_v1_raw);
    mTGE_v1_raw_vs_pT_ptSetA_centSetB[cent]->SetTitle(Form("v_{1}^{raw}, %3.f -%3.f%%",centSetB[cent],centSetB[cent+1]));
    mTGE_v1_raw_vs_pT_ptSetA_centSetB[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
    mTGE_v1_raw_vs_pT_ptSetA_centSetB[cent]->GetYaxis()->SetTitle("v_{1}");
    mTGE_v1_raw_vs_pT_ptSetA_centSetB[cent]->SetMarkerColor(4);
    mTGE_v1_raw_vs_pT_ptSetA_centSetB[cent]->SetMarkerStyle(24);
    mTGE_v1_raw_vs_pT_ptSetA_centSetB[cent]->Draw("AP");
    l1_ptSetA_centSetB->Draw("same");

    Double_t y_v1_reso[n_ptSetA_centSetB] = {d_FLow_ptSetA_centSetB[0][1][0][cent], d_FLow_ptSetA_centSetB[0][1][1][cent]};
    Double_t ey_v1_reso[n_ptSetA_centSetB] = {d_Flow_err_ptSetA_centSetB[0][1][0][cent], d_Flow_err_ptSetA_centSetB[0][1][1][cent]};
    canvas_v1_reso_vs_pT_ptSetA_centSetB->cd(cent+1);
    mTGE_v1_reso_vs_pT_ptSetA_centSetB[cent] = new TGraphErrors(n_ptSetA_centSetB,x,y_v1_reso,ex,ey_v1_reso);
    mTGE_v1_reso_vs_pT_ptSetA_centSetB[cent]->SetTitle(Form("v_{1}, %3.f -%3.f%%",centSetB[cent],centSetB[cent+1]));
    mTGE_v1_reso_vs_pT_ptSetA_centSetB[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
    mTGE_v1_reso_vs_pT_ptSetA_centSetB[cent]->GetYaxis()->SetTitle("v_{1}");
    mTGE_v1_reso_vs_pT_ptSetA_centSetB[cent]->SetMarkerColor(4);
    mTGE_v1_reso_vs_pT_ptSetA_centSetB[cent]->SetMarkerStyle(24);
    mTGE_v1_reso_vs_pT_ptSetA_centSetB[cent]->Draw("AP");
    l1_ptSetA_centSetB->Draw("same");

    Double_t y_v2_raw[n_ptSetA_centSetB] = {d_FLow_ptSetA_centSetB[1][0][0][cent], d_FLow_ptSetA_centSetB[1][0][1][cent]};
    Double_t ey_v2_raw[n_ptSetA_centSetB] = {d_Flow_err_ptSetA_centSetB[1][0][0][cent], d_Flow_err_ptSetA_centSetB[1][0][1][cent]};
    canvas_v2_raw_vs_pT_ptSetA_centSetB->cd(cent+1);
    mTGE_v2_raw_vs_pT_ptSetA_centSetB[cent] = new TGraphErrors(n_ptSetA_centSetB,x,y_v2_raw,ex,ey_v2_raw);
    mTGE_v2_raw_vs_pT_ptSetA_centSetB[cent]->SetTitle(Form("v_{2}^{raw}, %3.f -%3.f%%",centSetB[cent],centSetB[cent+1]));
    mTGE_v2_raw_vs_pT_ptSetA_centSetB[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
    mTGE_v2_raw_vs_pT_ptSetA_centSetB[cent]->GetYaxis()->SetTitle("v_{2}^{raw}");
    mTGE_v2_raw_vs_pT_ptSetA_centSetB[cent]->SetMarkerColor(4);
    mTGE_v2_raw_vs_pT_ptSetA_centSetB[cent]->SetMarkerStyle(24);
    mTGE_v2_raw_vs_pT_ptSetA_centSetB[cent]->Draw("AP");
    l1_ptSetA_centSetB->Draw("same");

    Double_t y_v2_reso[n_ptSetA_centSetB] = {d_FLow_ptSetA_centSetB[1][1][0][cent], d_FLow_ptSetA_centSetB[1][1][1][cent]};
    Double_t ey_v2_reso[n_ptSetA_centSetB] = {d_Flow_err_ptSetA_centSetB[1][1][0][cent], d_Flow_err_ptSetA_centSetB[1][1][1][cent]};
    canvas_v2_reso_vs_pT_ptSetA_centSetB->cd(cent+1);
    mTGE_v2_reso_vs_pT_ptSetA_centSetB[cent] = new TGraphErrors(n_ptSetA_centSetB,x,y_v2_reso,ex,ey_v2_reso);
    mTGE_v2_reso_vs_pT_ptSetA_centSetB[cent]->SetTitle(Form("v_{2}, %3.f -%3.f%%",centSetB[cent],centSetB[cent+1]));
    mTGE_v2_reso_vs_pT_ptSetA_centSetB[cent]->GetXaxis()->SetTitle("pT [GeV/c^{2}]");
    mTGE_v2_reso_vs_pT_ptSetA_centSetB[cent]->GetYaxis()->SetTitle("v_{2}");
    mTGE_v2_reso_vs_pT_ptSetA_centSetB[cent]->SetMarkerColor(4);
    mTGE_v2_reso_vs_pT_ptSetA_centSetB[cent]->SetMarkerStyle(24);
    mTGE_v2_reso_vs_pT_ptSetA_centSetB[cent]->Draw("AP");
    l1_ptSetA_centSetB->Draw("same");
  }

  // rap SetA, cent SetA
  // rap SetA, cent SetB
  TF1 * tf1_dv1dy = new TF1("tf1_dv1dy",proportion,0.,2.,1);
  for(int jkk=0;jkk<11;jkk++){
    cout<<"v1 10-40% ptSetA_centSetA_jkk_: " << jkk << " " << d_FLow_ptSetA_centSetA[jkk][0][1][0][1] <<", "<<d_FLow_ptSetA_centSetA[jkk][0][1][1][1] <<endl;
    cout<<"v1 Err 10-40% ptSetA_centSetA_jkk_: "<< jkk << " " << d_Flow_err_ptSetA_centSetA[jkk][0][1][0][1] <<", "<<d_Flow_err_ptSetA_centSetA[jkk][0][1][1][1] <<endl;
    cout<<endl;
  }

  for(int pt=0;pt<3;pt++){
  }
  for(int jkk=0;jkk<11;jkk++){
    cout<<"v2 10-40% ptSetA_centSetA_jkk_" <<jkk <<", " << d_FLow_ptSetA_centSetA[jkk][1][1][0][1] /*<<", "<<d_FLow_ptSetA_centSetA[jkk][1][1][1][1]*/ <<endl;
    // flowFile<<"v2 10-40% ptSetA_centSetA_jkk_" <<jkk <<", " << d_FLow_ptSetA_centSetA[jkk][1][1][0][1] /*<<", "<<d_FLow_ptSetA_centSetA[jkk][1][1][1][1]*/ <<endl;
    cout<<"v2 Err 10-40% ptSetA_centSetA_jkk_" <<jkk <<", " << d_Flow_err_ptSetA_centSetA[jkk][1][1][0][1] /*<<", "<<d_Flow_err_ptSetA_centSetA[jkk][1][1][1][1]*/ <<endl;
    // flowFile<<"v2 Err 10-40% ptSetA_centSetA_jkk_" <<jkk <<", " << d_Flow_err_ptSetA_centSetA[jkk][1][1][0][1] /*<<", "<<d_Flow_err_ptSetA_centSetA[jkk][1][1][1][1]*/ <<endl;
    // cout<<"v2 40-60% ptSetA_centSetA: " << d_FLow_ptSetA_centSetA[jkk][1][1][0][4] <<", "<<d_FLow_ptSetA_centSetA[jkk][1][1][1][4] <<endl;
    // cout<<"v2 Err 40-60% ptSetA_centSetA: " << d_Flow_err_ptSetA_centSetA[jkk][1][1][0][4] <<", "<<d_Flow_err_ptSetA_centSetA[jkk][1][1][1][4] <<endl;
  }
  cout<<endl;
  for(int jkk=0;jkk<11;jkk++){
    cout<<"v2 10-40% ptSetA_centSetA_jkk_" <<jkk <<", " /*<< d_FLow_ptSetA_centSetA[jkk][1][1][0][1] <<", "*/<<d_FLow_ptSetA_centSetA[jkk][1][1][1][1] <<endl;
    // flowFile<<"v2 10-40% ptSetA_centSetA_jkk_" <<jkk <<", " /*<< d_FLow_ptSetA_centSetA[jkk][1][1][0][1] <<", "*/<<d_FLow_ptSetA_centSetA[jkk][1][1][1][1] <<endl;
    cout<<"v2 Err 10-40% ptSetA_centSetA_jkk_" <<jkk <<", " /*<< d_Flow_err_ptSetA_centSetA[jkk][1][1][0][1] <<", "*/<<d_Flow_err_ptSetA_centSetA[jkk][1][1][1][1] <<endl;
    // flowFile<<"v2 Err 10-40% ptSetA_centSetA_jkk_" <<jkk <<", " /*<< d_Flow_err_ptSetA_centSetA[jkk][1][1][0][1] <<", "*/<<d_Flow_err_ptSetA_centSetA[jkk][1][1][1][1] <<endl;
    // cout<<"v2 40-60% ptSetA_centSetA: " << d_FLow_ptSetA_centSetA[jkk][1][1][0][4] <<", "<<d_FLow_ptSetA_centSetA[jkk][1][1][1][4] <<endl;
    // cout<<"v2 Err 40-60% ptSetA_centSetA: " << d_Flow_err_ptSetA_centSetA[jkk][1][1][0][4] <<", "<<d_Flow_err_ptSetA_centSetA[jkk][1][1][1][4] <<endl;
  }
  cout << endl;
  cout<<"v2 10-40% ptSetA_centSetA: " << d_FLow_ptSetA_centSetA[0][1][1][0][1] <<", "<<d_FLow_ptSetA_centSetA[0][1][1][1][1]   <<endl;
  flowFile << d_FLow_ptSetA_centSetA[0][1][1][0][1] <<" "<<d_FLow_ptSetA_centSetA[0][1][1][1][1] <<endl;
  cout<<"v2 Err 10-40% ptSetA_centSetA: " << d_Flow_err_ptSetA_centSetA[0][1][1][0][1] <<", "<<d_Flow_err_ptSetA_centSetA[0][1][1][1][1]   <<endl;
  flowFile << d_Flow_err_ptSetA_centSetA[0][1][1][0][1] <<" "<<d_Flow_err_ptSetA_centSetA[0][1][1][1][1] <<endl;
  outputFile->cd();
  canvas_InvM_ptSetA_centSetA->Write();
  canvas_v1_raw_ptSetA_centSetA->Write();
  canvas_v1_reso_ptSetA_centSetA->Write();
  canvas_v2_raw_ptSetA_centSetA->Write();
  canvas_v2_reso_ptSetA_centSetA->Write();
  canvas_v1_raw_vs_pT_ptSetA_centSetA->Write();
  canvas_v1_reso_vs_pT_ptSetA_centSetA->Write();
  canvas_v2_raw_vs_pT_ptSetA_centSetA->Write();
  canvas_v2_reso_vs_pT_ptSetA_centSetA->Write();
  canvas_InvM_ptSetB_centSetA->Write();
  canvas_v1_raw_ptSetB_centSetA->Write();
  canvas_v1_reso_ptSetB_centSetA->Write();
  canvas_v2_raw_ptSetB_centSetA->Write();
  canvas_v2_reso_ptSetB_centSetA->Write();
  canvas_v1_raw_vs_pT_ptSetB_centSetA->Write();
  canvas_v1_reso_vs_pT_ptSetB_centSetA->Write();
  canvas_v2_raw_vs_pT_ptSetB_centSetA->Write();
  canvas_v2_reso_vs_pT_ptSetB_centSetA->Write();
  canvas_InvM_ptSetA_centSetB->Write();
  canvas_v1_raw_ptSetA_centSetB->Write();
  canvas_v1_reso_ptSetA_centSetB->Write();
  canvas_v2_raw_ptSetA_centSetB->Write();
  canvas_v2_reso_ptSetA_centSetB->Write();
  canvas_v1_raw_vs_pT_ptSetA_centSetB->Write();
  canvas_v1_reso_vs_pT_ptSetA_centSetB->Write();
  canvas_v2_raw_vs_pT_ptSetA_centSetB->Write();
  canvas_v2_reso_vs_pT_ptSetA_centSetB->Write();
  // outputFile->Write();
  flowFile.close();
}

Double_t proportion(Double_t *x, Double_t *p)
{
  return ( p[0] * x[0]);
}

// Fitting functions for Flow VS Invariant Mass
Double_t BackgroundFitting(Double_t *x, Double_t *p)
{
  if( x[0] >= (dParSig[1] - (_sigmaRange*dParSig[2])) && x[0] <= (dParSig[1] + (_sigmaRange*dParSig[2])))
   {
     return (((dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))
     /((dParSig[0]*exp(-0.5*pow(((x[0]-dParSig[1])/dParSig[2]),2))+dParSig[3])+(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3)*/);
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3)*/;
   }
}
Double_t TotalFitting(Double_t *x, Double_t *p)
{
  if( x[0] >= (dParSig[1] - (_sigmaRange*dParSig[2])) && x[0] <= (dParSig[1] + (_sigmaRange*dParSig[2])))
   {
     return (
       (p[/*1*//*2*/3/*4*/]*((dParSig[0]*exp(-0.5*pow(((x[0]-dParSig[1])/dParSig[2]),2))+dParSig[3])
     /((dParSig[0]*exp(-0.5*pow(((x[0]-dParSig[1])/dParSig[2]),2))+dParSig[3])+(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2)))))
     +
       (((dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))
     /((dParSig[0]*exp(-0.5*pow(((x[0]-dParSig[1])/dParSig[2]),2))+dParSig[3])+(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3))*/));
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3))*/;
   }
}
