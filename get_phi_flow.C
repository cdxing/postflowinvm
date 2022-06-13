#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGaxis.h"
#include <fstream>      // std::filebuf


Double_t BackgroundFitting(Double_t *x, Double_t *p);
Double_t TotalFitting(Double_t *x, Double_t *p);
Double_t SignalFitting(Double_t *x, Double_t *p);
Double_t PolyBreitWigner(Double_t *x_val, Double_t *par)
{
    Double_t x = x_val[0];
    Double_t m0 = par[0];
    Double_t Gamma = par[1];
    Double_t Norm = par[2];

    Double_t denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
    Double_t BW = Norm*Gamma/denom;

    Double_t Poly = par[3] + par[4]*x;

    Double_t y = BW + Poly;

    return y;
}

Double_t Poly(Double_t *x_val, Double_t *par)
{
    Double_t x = x_val[0];
    Double_t y = par[0] + par[1]*x;

    return y;
}

Double_t Poly2(Double_t *x_val, Double_t *par)
{
    Double_t x = x_val[0];
    Double_t y = par[0] + par[1]*x + par[2]*x*x;

    return y;
}

Double_t BreitWigner(Double_t *x_val, Double_t *par)
{
    Double_t x = x_val[0];
    Double_t m0 = par[0];
    Double_t Gamma = par[1];
    Double_t Norm = par[2];

    Double_t denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
    Double_t BW = Norm*Gamma/denom;

    return BW;
}

Double_t flow_2(Double_t *x_val, Double_t *par)
{
    Double_t x, y;
    Double_t Ampl, v2;
    x = x_val[0];
    Ampl = par[0];
    v2 = par[1];

    y = Ampl*(1.0 + 2.0*v2*TMath::Cos(2.0*x));

    return y;
}

Double_t flow_3(Double_t *x_val, Double_t *par)
{
    Double_t x, y;
    Double_t Ampl, v3;
    x = x_val[0];
    Ampl = par[0];
    v3 = par[1];

    y = Ampl*(1.0 + 2.0*v3*TMath::Cos(3.0*x));

    return y;
}

Double_t Gaussion(Double_t *x_val, Double_t *par)
{
    Double_t x, y;
    x = x_val[0];

    Double_t mu = par[0];
    Double_t sigma = par[1];
    Double_t norm = par[2];

    y = norm*TMath::Exp(-1.0*(x-mu)*(x-mu)/(2.0*sigma*sigma))/(sigma*TMath::Sqrt(2*TMath::Pi()));

    return y;
}

TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1)
{
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->Draw();
    return text;
}

void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y1_val);
    Zero_line -> SetY2(y2_val);
    Zero_line -> SetLineWidth(LineWidth);
    Zero_line -> SetLineStyle(LineStyle);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}
void defineStyle()
{
    gStyle->SetLabelColor(1,"X");
    gStyle->SetLabelColor(1,"Y");
    gStyle->SetTitleTextColor(1);
    gStyle->SetNdivisions(6, "X");
    gStyle->SetNdivisions(4, "Y");
    gStyle->SetMarkerStyle(8);
    gStyle->SetMarkerSize(1.7);
    gStyle->SetHistLineColor(2);
    gStyle->SetHistLineWidth(1);
    gStyle->SetFuncWidth(3);
    gStyle->SetPadBorderMode(0);
    gStyle->SetHistLineStyle(0);
    gStyle->SetTitleH(0.058);
    gStyle->SetTitleW(0.98);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(0);
    gStyle->SetOptDate(0);
    gStyle->SetTitleFont(42);
    gStyle->SetLabelFont(42);
    TGaxis::SetMaxDigits(3);
}
//const Double_t _sigmaRange = 5.; // Sigma of the Fitting range // use the nSigmaPhi
static const TString Energy[4] = {"200GeV","20GeV","15GeV","19GeV",};
static const Int_t pt_total_phi = 23;
static const Int_t Centrality_total = 4;
static const Int_t Centrality_start = 0;
static const Int_t Centrality_stop = 4;
static const Int_t Mode_total = 2;
static const Int_t EtaGap_total = 1;
static const Int_t EtaGap_start = 0;
static const Int_t EtaGap_stop = 1;
static const Int_t Phi_Psi_total = 7;
static const Float_t nSigmaPhi = 3;
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
static const Int_t cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-30%, 2 = 30-80%
static const Int_t cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-30%, 2 = 30-80%
//static const Int_t cent_low[4] = {1,1,3,6}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
//static const Int_t cent_up[4]  = {9,2,5,9}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
static const Float_t BW_Start = 0.994;
static const Float_t BW_Stop  = 1.050;

//static const Int_t pt_total_New_phi = 11;
//static Float_t pt_low_phi[pt_total_New_phi] =     {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8};
//static Float_t pt_up_phi[pt_total_New_phi]  =     {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.4};
//static Int_t pt_new_bin_start[pt_total_New_phi] = {1  ,2  ,3  ,4  ,5  ,6  ,7  ,8  ,9  ,11 ,13};
//static Int_t pt_new_bin_stop[pt_total_New_phi]  = {1  ,2  ,3  ,4  ,5  ,6  ,7  ,8  ,10 ,12 ,14};

static const Int_t pt_total_New_phi = 8;
static Float_t pt_low_phi[pt_total_New_phi] =     {0.4,0.8,1.0,1.2,1.4,1.6,2.0,2.4};
static Float_t pt_up_phi[pt_total_New_phi]  =     {0.8,1.0,1.2,1.4,1.6,2.0,2.4,3.4};
static Int_t pt_new_bin_start[pt_total_New_phi] = {1  ,3  ,4  ,5  ,6  ,7  ,9  ,11};
static Int_t pt_new_bin_stop[pt_total_New_phi]  = {2  ,3  ,4  ,5  ,6  ,8  ,10 ,14};


Double_t dParBg[3]; // Bkg fitting parameters
Double_t dParSig[5]; // Sig + Bkg fitting parameters

//static const Int_t pt_total_New_phi = 1;
//static Float_t pt_low_phi[pt_total_New_phi] =     {0.6};
//static Float_t pt_up_phi[pt_total_New_phi]  =     {2.0};
//static Int_t pt_new_bin_start[pt_total_New_phi] = {2 };
//static Int_t pt_new_bin_stop[pt_total_New_phi]  = {8 };


void get_phi_flow(const Int_t energy = 3)
{
    defineStyle();
    // open input file for same event and mixed event
    //TString inputfile_SE = Form("data/out_PhiMesonAna_%s_test1.root",Energy[energy].Data());
    TString inputfile_SE = Form("data/Yields_SE_%s.root",Energy[energy].Data());
    TFile *File_Ana_SE = TFile::Open(inputfile_SE.Data());
    //TString inputfile_SE = Form("../../processed_v2_phi_eventplane_%s_FXT/cooked_phi_flow.root",Energy[energy].Data());
    //TFile *File_Ana = TFile::Open(inputfile_SE.Data());
    TString inputfile_ME = Form("data/Yields_ME_%s.root",Energy[energy].Data());
    TFile *File_Ana_ME = TFile::Open(inputfile_ME.Data());
    // TFile *File_Ana = TFile::Open(inputfile_ME.Data());
    cout << " file read successfully! " <<  endl;
    // flowinvm
    //  0: pt; 1: centrality; 2: Eta_gap
    TH1F *h_mMass_Spec_SE[23][4][4];
    TH1F *h_mMass_Spec_ME[23][4][4];
    TProfile *p_mMass2_invMfit[23][4][4];


    //TH1F *h_mMass_Phi2_Ana[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    //TH1F *h_mMass_Phi2_Ana[pt_total_phi][Centrality_total][Mode_total][Phi_Psi_total];

    // read histogram for same event and mixed event
    // TH1F *h_mMass_Phi2_SE[pt_total_phi][Centrality_total][Phi_Psi_total];
    TH1F *h_mMass_Phi2_SE[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  // TH1F *h_mMass_Phi3_SE[pt_total_phi][Centrality_total][Phi_Psi_total];

    TH1F *h_mMass_Phi2_ME[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    // TH1F *h_mMass_Phi3_ME[pt_total_phi][Centrality_total][Phi_Psi_total];
    cout << "test 0 successfully! " <<  endl;
    std::ofstream txtFile("./flow/postflowep.txt",ofstream::out);
    std::ofstream txtFile_invm("./flow/postflowinvm.txt",ofstream::out);
    // flowinv
    // input
    // raw pt spectra; invM Fit flow plots | TODO: use finer pt_bin
    for(Int_t i = 0; i < pt_total_phi; i++) // pt bin
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {   
                TString Mode[2] = {"SE","ME"};
                TString HistName; 
		// SE invarian mass spectrum
		Int_t X_flag = 0;
                HistName = Form("Spec_pt_%d_Centrality_%d_EtaGap_%d_Phi_%s",i,j,l,Mode[X_flag].Data());
                h_mMass_Spec_SE[i][j][l] = (TH1F*)File_Ana_SE->Get(HistName.Data())->Clone();
                
		// invMfit flow
                HistName = Form("InvMfit_pt_%d_Centrality_%d_EtaGap_%d_Phi_%s",i,j,l,Mode[X_flag].Data());
                p_mMass2_invMfit[i][j][l] = (TProfile*)File_Ana_SE->Get(HistName.Data())->Clone();

		// ME invarian mass spectrum
		X_flag = 1;
                HistName = Form("Spec_pt_%d_Centrality_%d_EtaGap_%d_Phi_%s",i,j,l,Mode[X_flag].Data());
                h_mMass_Spec_ME[i][j][l] = (TH1F*)File_Ana_ME->Get(HistName.Data())->Clone();

            }   
        }
    }   


    for(Int_t i = 0; i < pt_total_phi; i++) // pt bin
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
          //for(Int_t X_flag = 0; X_flag < 2; X_flag++) // x_flag
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                  //TString Mode[2] = {"SE","ME"};
                  TString HistName;
                  //HistName = Form("pt_%d_Centrality_%d_phi_Psi_%d_2nd_%s_%s",i,j,m,"Phi",Mode[X_flag].Data());
                  //h_mMass_Phi2_Ana[i][j][X_flag][m] = (TH1F*)File_Ana->Get(HistName.Data())->Clone();

                  //HistName = Form("pt_%d_Centrality_%d_phi_Psi_%d_2nd_Phi_SE",i,j,m);
                  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_Phi_SE",i,j,l,m);
                  h_mMass_Phi2_SE[i][j][l][m] = (TH1F*)File_Ana_SE->Get(HistName.Data())->Clone();
                  //HistName = Form("pt_%d_Centrality_%d_phi_Psi_%d_3rd_Phi_SE",i,j,m);
                  //HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_Phi_SE",i,j,l,m);
                  //h_mMass_Phi3_SE[i][j][l][m] = (TH1F*)File_Ana_SE->Get(HistName.Data())->Clone();

                  //HistName = Form("pt_%d_Centrality_%d_phi_Psi_%d_2nd_Phi_ME",i,j,m);
                  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_Phi_ME",i,j,l,m);
                  h_mMass_Phi2_ME[i][j][l][m] = (TH1F*)File_Ana_ME->Get(HistName.Data())->Clone();
                  //HistName = Form("pt_%d_Centrality_%d_phi_Psi_%d_3rd_Phi_ME",i,j,m);
                  //HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_Phi_ME",i,j,l,m);
                  //h_mMass_Phi3_ME[i][j][l][m] = (TH1F*)File_Ana_ME->Get(HistName.Data())->Clone();
                }
            }
        }
    }

    // calculate scaling factor and subtract mixed event from same event
    TH1F *h_mMass_Phi2_SM[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    // TH1F *h_mMass_Phi3_SM[pt_total_phi][Centrality_total][Phi_Psi_total];
    for(Int_t i = 0; i < pt_total_phi; i++) // pt bin
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
             for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
             {
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    //Int_t bin2_SE_start = h_mMass_Phi2_Ana[i][j][0][m]->FindBin(1.04);
                    //Int_t bin2_SE_stop  = h_mMass_Phi2_Ana[i][j][0][m]->FindBin(1.05);
                    //Float_t Inte2_SE = h_mMass_Phi2_Ana[i][j][0][m]->Integral(bin2_SE_start,bin2_SE_stop);
                    //h_mMass_Phi2_SM[i][j][m] = (TH1F*)h_mMass_Phi2_Ana[i][j][0][m]->Clone();

                    Int_t bin2_SE_start = h_mMass_Phi2_SE[i][j][l][m]->FindBin(1.04);
                    Int_t bin2_SE_stop  = h_mMass_Phi2_SE[i][j][l][m]->FindBin(1.05);
                    Float_t Inte2_SE = h_mMass_Phi2_SE[i][j][l][m]->Integral(bin2_SE_start,bin2_SE_stop);
                    h_mMass_Phi2_SM[i][j][l][m] = (TH1F*)h_mMass_Phi2_SE[i][j][l][m]->Clone();
                    
		    // Int_t bin3_SE_start = h_mMass_Phi3_SE[i][j][m]->FindBin(1.04);
                    // Int_t bin3_SE_stop  = h_mMass_Phi3_SE[i][j][m]->FindBin(1.05);
                    // Float_t Inte3_SE = h_mMass_Phi3_SE[i][j][m]->Integral(bin3_SE_start,bin3_SE_stop);
                    // h_mMass_Phi3_SM[i][j][m] = (TH1F*)h_mMass_Phi3_SE[i][j][m]->Clone();

                    //Int_t bin2_ME_start = h_mMass_Phi2_Ana[i][j][1][m]->FindBin(1.04);
                    //Int_t bin2_ME_stop  = h_mMass_Phi2_Ana[i][j][1][m]->FindBin(1.05);
                    //Float_t Inte2_ME = h_mMass_Phi2_Ana[i][j][1][m]->Integral(bin2_ME_start,bin2_ME_stop);
                    //h_mMass_Phi2_Ana[i][j][1][m]->Scale(Inte2_SE/Inte2_ME);
                    //h_mMass_Phi2_SM[i][j][m]->Add(h_mMass_Phi2_Ana[i][j][1][m],-1.0);

                    Int_t bin2_ME_start = h_mMass_Phi2_ME[i][j][l][m]->FindBin(1.04);
                    Int_t bin2_ME_stop  = h_mMass_Phi2_ME[i][j][l][m]->FindBin(1.05);
                    Float_t Inte2_ME = h_mMass_Phi2_ME[i][j][l][m]->Integral(bin2_ME_start,bin2_ME_stop);
                    h_mMass_Phi2_ME[i][j][l][m]->Scale(Inte2_SE/Inte2_ME);
                    h_mMass_Phi2_SM[i][j][l][m]->Add(h_mMass_Phi2_ME[i][j][l][m],-1.0);
                    
		    // Int_t bin3_ME_start = h_mMass_Phi3_ME[i][j][m]->FindBin(1.04);
                    // Int_t bin3_ME_stop  = h_mMass_Phi3_ME[i][j][m]->FindBin(1.05);
                    // Float_t Inte3_ME = h_mMass_Phi3_ME[i][j][m]->Integral(bin3_ME_start,bin3_ME_stop);
                    // h_mMass_Phi3_ME[i][j][m]->Scale(Inte3_SE/Inte3_ME);
                    // h_mMass_Phi3_SM[i][j][m]->Add(h_mMass_Phi3_ME[i][j][m],-1.0);
                }
           }
        }
    }
    cout << "test 1 ! " <<  endl;

    TH1::SetDefaultSumw2();

    // same event and mixed event distribution for QA usage
    //TH1F *h_mMass_Phi2_Ana_New_Bin[pt_total_New_phi][Centrality_total][Mode_total][Phi_Psi_total];
    //TH1F *h_mMass_Phi2_SE_New_Bin[pt_total_New_phi][Centrality_total][Phi_Psi_total];
    //TH1F *h_mMass_Phi2_ME_New_Bin[pt_total_New_phi][Centrality_total][Phi_Psi_total];
    //TH1F *h_mMass_Phi2_SM_new_bin[pt_total_New_phi][Centrality_total][Phi_Psi_total];
    TH1F *h_mMass_Phi2_SE_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    TH1F *h_mMass_Phi2_ME_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    TH1F *h_mMass_Phi2_SM_new_bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    TH1F *h_mMass_Phi2_SM_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    TF1 *f_gauss_Spec_2[pt_total_New_phi][Centrality_total][EtaGap_total];
    TF1 *f_PolyBW_Spec_total[pt_total_New_phi][Centrality_total][EtaGap_total];
    TF1 * f_Background[pt_total_New_phi][Centrality_total][EtaGap_total];
    TF1 * f_flow_2_Sig[pt_total_New_phi][Centrality_total][EtaGap_total];
    TF1 * f_flow_2_Bg[pt_total_New_phi][Centrality_total][EtaGap_total];
    TF1 * f_global[pt_total_New_phi][Centrality_total][EtaGap_total];

    Float_t ParFit2_Spec_total[pt_total_New_phi][Centrality_total][EtaGap_total][5];
    // flowinvm 
    //  0: pt; 1: centrality; 2: Eta_gap
    TH1F *h_mMass_Spec_SE_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total];
    TH1F *h_mMass_Spec_ME_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total];
    TH1F *h_mMass_Spec_SM_new_bin[pt_total_New_phi][Centrality_total][EtaGap_total];
    TH1F *h_mMass_Spec_SM_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total];
    TProfile *p_mMass2_invMfit_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total];

    TCanvas *c1_invmass_Spec[pt_total_New_phi][Centrality_total][EtaGap_total];
    TCanvas *c1_invmass_Spec_SM[pt_total_New_phi][Centrality_total][EtaGap_total];
    TCanvas *c1_invmass_invMfit[pt_total_New_phi][Centrality_total][EtaGap_total];
    TCanvas *c1_invmass[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    char name[80];
    Float_t gaussSig = 1.5; // dchen
    txtFile_invm << Form("##columns: pT ,  v2 , stat.error, syst.error") << endl;
    txtFile_invm << Form("#---------- phi-meson v2 ------------") << endl;
    txtFile_invm << Form("# invM fitting method ===> ") << endl;
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
    	txtFile_invm << Form("# phi-meson v2 (centrality bin %d) in 19.6 GeV ===> ", j) << endl;
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            //for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
            //{
                for(Int_t i = 0; i < pt_total_New_phi; i++)
                {
                    TString HistName;
                    TString FuncName_g;
                    TString FuncName_bw;
                    TString FuncName_sig;
                    TString FuncName_bg;
                    TString FuncName_flow_sig;
                    TString FuncName_flow_bg;
                    TString FuncName_global;

                    HistName = Form("Sig_Spec_pt_%d_Centality_%d_EtaGap_%d_2nd",i,j,l);
                    FuncName_g = Form("f_Spec_gauss_pt_%d_Centrality_%d_EtaGap_%d_2nd",i,j,l);
                    FuncName_bw = Form("f_Spec_bw_pt_%d_Centrality_%d_EtaGap_%d_2nd",i,j,l);
                    FuncName_sig = Form("f_Spec_sig_pt_%d_Centrality_%d_EtaGap_%d_2nd",i,j,l);
                    FuncName_bg = Form("f_Spec_bg_pt_%d_Centrality_%d_EtaGap_%d_2nd",i,j,l);
                    FuncName_flow_sig = Form("f_flow_sig_pt_%d_Centrality_%d_EtaGap_%d_2nd",i,j,l);
                    FuncName_flow_bg = Form("f_flow_bg_pt_%d_Centrality_%d_EtaGap_%d_2nd",i,j,l);
                    FuncName_global = Form("f_flow_global_pt_%d_Centrality_%d_EtaGap_%d_2nd",i,j,l);
                    for(Int_t pt_bin = pt_new_bin_start[i]; pt_bin <= pt_new_bin_stop[i]; pt_bin++)
                    {
                        if(pt_bin == pt_new_bin_start[i])
                        {
                            h_mMass_Spec_SE_New_Bin[i][j][l] = (TH1F*)h_mMass_Spec_SE[pt_bin][j][l]->Clone();
                            h_mMass_Spec_ME_New_Bin[i][j][l] = (TH1F*)h_mMass_Spec_ME[pt_bin][j][l]->Clone();
                            p_mMass2_invMfit_New_Bin[i][j][l] = (TProfile*)p_mMass2_invMfit[pt_bin][j][l]->Clone();
                        }
                        else
                        {
                            h_mMass_Spec_SE_New_Bin[i][j][l] ->Add(h_mMass_Spec_SE[pt_bin][j][l],1.0);
                            h_mMass_Spec_ME_New_Bin[i][j][l] ->Add(h_mMass_Spec_ME[pt_bin][j][l],1.0);
                            p_mMass2_invMfit_New_Bin[i][j][l] ->Add(p_mMass2_invMfit[pt_bin][j][l],1.0);
                        }
                    }
                    h_mMass_Spec_SM_new_bin[i][j][l] = (TH1F*)h_mMass_Spec_SE_New_Bin[i][j][l]->Clone();
                    h_mMass_Spec_SM_New_Bin[i][j][l] = (TH1F*)h_mMass_Spec_SE_New_Bin[i][j][l]->Clone();
                    Int_t bin2_SE_start = h_mMass_Spec_SE_New_Bin[i][j][l]->FindBin(1.04);
                    Int_t bin2_SE_stop  = h_mMass_Spec_SE_New_Bin[i][j][l]->FindBin(1.05);
                    Float_t Inte2_SE = h_mMass_Spec_SE_New_Bin[i][j][l]->Integral(bin2_SE_start,bin2_SE_stop);
                    Int_t bin2_ME_start = h_mMass_Spec_ME_New_Bin[i][j][l]->FindBin(1.04);
                    Int_t bin2_ME_stop  = h_mMass_Spec_ME_New_Bin[i][j][l]->FindBin(1.05);
                    Float_t Inte2_ME = h_mMass_Spec_ME_New_Bin[i][j][l]->Integral(bin2_ME_start,bin2_ME_stop);
                    //h_mMass_Spec_Ana_New_Bin[i][j][1]->Scale(Inte2_SE/Inte2_ME);
                    h_mMass_Spec_ME_New_Bin[i][j][l]->Scale(Inte2_SE/Inte2_ME);
                    //h_mMass_Spec_SM_New_Bin[i][j][l]->Scale(Inte2_SE/Inte2_ME);
                    //h_mMass_Spec_SM_new_bin[i][j][l]->Scale(Inte2_SE/Inte2_ME);

                    h_mMass_Spec_SM_new_bin[i][j][l]->Add(h_mMass_Spec_SE_New_Bin[i][j][l],-1.0);
                    h_mMass_Spec_SM_New_Bin[i][j][l]->Add(h_mMass_Spec_ME_New_Bin[i][j][l],-1.0);

                    f_gauss_Spec_2[i][j][l] = new TF1(FuncName_g.Data(),Gaussion,BW_Start,BW_Stop,3); 
		    f_gauss_Spec_2[i][j][l]->SetParameter(0,1.019);
                    f_gauss_Spec_2[i][j][l]->SetParameter(1,0.07);
                    f_gauss_Spec_2[i][j][l]->SetParameter(2,1000);
                    f_gauss_Spec_2[i][j][l]->SetRange(1.019-gaussSig*0.07,1.019+gaussSig*0.07); // dchen
                    h_mMass_Spec_SM_New_Bin[i][j][l]->Fit(f_gauss_Spec_2[i][j][l],"MQNR","MQNR",1.,1.05);
                    //h_mMass_Spec_SM_New_Bin[i][j][l]->Fit(f_gauss_Spec_2[i][j][l],"MQR");
                    // seeds
                    Double_t d_seeds_mu     = f_gauss_Spec_2[i][j][l] -> GetParameter(0);
                    Double_t d_seeds_sigma  = f_gauss_Spec_2[i][j][l] -> GetParameter(1);
                    Double_t d_seeds_norm   = f_gauss_Spec_2[i][j][l] -> GetParameter(2);
                    //f_gauss_Spec_2[i][j][l]->SetRange(ParFit2_total[i][j][0]-gaussSig*ParFit2_total[i][j][1],ParFit2_total[i][j][0]+gaussSig*ParFit2_total[i][j][1]); // dchen
                    f_PolyBW_Spec_total[i][j][l] = new TF1(FuncName_bw.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
                    f_Background[i][j][l] = new TF1(FuncName_bg.Data(),Poly2,BW_Start,BW_Stop,5);
                    f_Background[i][j][l]->SetRange(d_seeds_mu-nSigmaPhi*d_seeds_sigma,d_seeds_mu+nSigmaPhi*d_seeds_sigma);
                    for(Int_t n_par = 0; n_par < 5; n_par++)
                    {
                        f_PolyBW_Spec_total[i][j][l]->ReleaseParameter(n_par);
                    }
                    f_PolyBW_Spec_total[i][j][l]->SetParameter(0,d_seeds_mu);
                    //f_PolyBW_Spec_total[i][j][l]->SetParameter(0,1.019);
                    f_PolyBW_Spec_total[i][j][l]->SetParLimits(0,1.014,1.024);
                    f_PolyBW_Spec_total[i][j][l]->SetParameter(1,d_seeds_sigma);
                    //f_PolyBW_Spec_total[i][j][l]->SetParameter(1,0.0055);
                    f_PolyBW_Spec_total[i][j][l]->SetParameter(2,d_seeds_norm);
                    //f_PolyBW_Spec_total[i][j][l]->SetParameter(2,10000);
                    f_PolyBW_Spec_total[i][j][l]->SetParameter(3,-6000);
                    f_PolyBW_Spec_total[i][j][l]->SetParameter(4,0.5);
                    f_PolyBW_Spec_total[i][j][l]->SetRange(d_seeds_mu-nSigmaPhi*d_seeds_sigma,d_seeds_mu+nSigmaPhi*d_seeds_sigma);
                    //h_mMass_Spec_SM_New_Bin[i][j][l]->Fit(f_gauss_Spec_2[i][j][l],"E+","QR",1.019-gaussSig*0.008,1.019+gaussSig*0.008);
                    h_mMass_Spec_ME_New_Bin[i][j][l]->Fit(f_Background[i][j][l],"E+","MQR",1.,1.05);
                    h_mMass_Spec_SM_New_Bin[i][j][l]->Fit(f_PolyBW_Spec_total[i][j][l],"E+","MQR",1.,1.05);
                    f_flow_2_Sig[i][j][l] = new TF1(FuncName_flow_sig.Data(),SignalFitting,BW_Start,BW_Stop,1);
                    f_flow_2_Bg[i][j][l] = new TF1(FuncName_flow_bg.Data(),BackgroundFitting,BW_Start,BW_Stop,5);
                    f_global[i][j][l] = new TF1(FuncName_global.Data(),TotalFitting,BW_Start,BW_Stop,4); 
		    dParSig[0] = f_PolyBW_Spec_total[i][j][l] -> GetParameter(0);
		    dParSig[1] = f_PolyBW_Spec_total[i][j][l] -> GetParameter(1);
		    dParSig[2] = f_PolyBW_Spec_total[i][j][l] -> GetParameter(2);
		    dParSig[3] = f_PolyBW_Spec_total[i][j][l] -> GetParameter(3);
		    dParSig[4] = f_PolyBW_Spec_total[i][j][l] -> GetParameter(4);
		    dParBg[0] = f_Background[i][j][l] -> GetParameter(0);
		    dParBg[1] = f_Background[i][j][l] -> GetParameter(1);
		    dParBg[2] = f_Background[i][j][l] -> GetParameter(2);

                    sprintf(name, "c1_invmass_Spec_Pt%d_Cent%d",i,j);
                    c1_invmass_Spec[i][j][l] = new TCanvas(name, name, 600, 500);
                    c1_invmass_Spec[i][j][l]->cd();
                    gPad->SetTopMargin(0.1);
                    gPad->SetRightMargin(0.1);
                    gPad->SetLeftMargin(0.2);
                    gPad->SetBottomMargin(0.2);
                    h_mMass_Spec_SE_New_Bin[i][j][l]->SetMarkerColor(1);
                    h_mMass_Spec_SE_New_Bin[i][j][l]->SetNdivisions(505,"X");
                    h_mMass_Spec_SE_New_Bin[i][j][l]->SetNdivisions(505,"Y");
                    h_mMass_Spec_SE_New_Bin[i][j][l]->SetYTitle("Counts");
                    h_mMass_Spec_SE_New_Bin[i][j][l]->SetXTitle("m_{inv} (GeV/c^{2})");
                    h_mMass_Spec_SE_New_Bin[i][j][l]->GetXaxis()->CenterTitle();
                    h_mMass_Spec_SE_New_Bin[i][j][l]->GetYaxis()->CenterTitle();
                    h_mMass_Spec_SE_New_Bin[i][j][l]->SetLineColor(1);//
                    //h_mMass_Spec_SE_New_Bin[i][j][l]->Draw("HIST");//
                    h_mMass_Spec_SE_New_Bin[i][j][l]->Draw();//
                    h_mMass_Spec_ME_New_Bin[i][j][l]->SetMarkerColor(2);
                    h_mMass_Spec_ME_New_Bin[i][j][l]->SetLineColor(2);
                    h_mMass_Spec_ME_New_Bin[i][j][l]->SetFillColor(3002);
                    //h_mMass_Spec_ME_New_Bin[i][j][l]->Draw("HIST SAME");
                    h_mMass_Spec_ME_New_Bin[i][j][l]->Draw("SAME");
                    h_mMass_Spec_SM_new_bin[i][j][l]->SetLineColor(4);
                    h_mMass_Spec_SM_New_Bin[i][j][l]->SetLineColor(4);
                    h_mMass_Spec_SM_New_Bin[i][j][l]->SetFillStyle(3003);
                    //h_mMass_Spec_SM_New_Bin[i][j][l]->Draw("HIST SAME");
                    //f_gauss_Spec_2[i][j][l]->Draw("SAME");
                    //h_mMass_Spec_SM_new_bin[i][j][l]->SetFillStyle(3003);
                    //h_mMass_Spec_SM_new_bin[i][j][l]->Draw("HIST SAME");
                    sprintf(name, "figures_phi/invmass_Pt%d_Cent%d.pdf",i,j);
                    c1_invmass_Spec[i][j][l]->SaveAs(name);
                    
                    sprintf(name, "c1_invmass_Spec_SM_Pt%d_Cent%d",i,j);
                    c1_invmass_Spec_SM[i][j][l] = new TCanvas(name, name, 600, 500);
                    c1_invmass_Spec_SM[i][j][l]->cd();
                    gPad->SetTopMargin(0.1);
                    gPad->SetRightMargin(0.1);
                    gPad->SetLeftMargin(0.2);
                    gPad->SetBottomMargin(0.2);
                    h_mMass_Spec_SM_New_Bin[i][j][l]->Draw();
                    sprintf(name, "figures_phi/invmass_SM_Pt%d_Cent%d.pdf",i,j);
                    c1_invmass_Spec_SM[i][j][l]->SaveAs(name);
  		    
		    sprintf(name, "c1_invmass_invMfit_Pt%d_Cent%d",i,j);
                    c1_invmass_invMfit[i][j][l] = new TCanvas(name, name, 600, 500);
                    c1_invmass_invMfit[i][j][l]->cd();
                    gPad->SetTopMargin(0.1);
                    gPad->SetRightMargin(0.1);
                    gPad->SetLeftMargin(0.2);
                    gPad->SetBottomMargin(0.2);
		    f_flow_2_Sig[i][j][l]->SetLineColor(kBlue);
		    f_flow_2_Bg[i][j][l]->SetLineColor(kRed);
                    p_mMass2_invMfit_New_Bin[i][j][l]->Fit(f_flow_2_Sig[i][j][l],"MNQE+","MQR",1.,1.05);
                    p_mMass2_invMfit_New_Bin[i][j][l]->Fit(f_flow_2_Bg[i][j][l],"MNQE+","MQR",1.,1.05);
                    //p_mMass2_invMfit_New_Bin[i][j][l]->Fit(f_flow_2_Sig[i][j][l],"QE+","MQR",1.,1.05);
		    double d_flow_seed = f_flow_2_Sig[i][j][l] -> GetParameter(0);
		    double d_flow_p0_seed = f_flow_2_Bg[i][j][l] -> GetParameter(0);
		    double d_flow_p1_seed = f_flow_2_Bg[i][j][l] -> GetParameter(1);
		    double d_flow_p2_seed = f_flow_2_Bg[i][j][l] -> GetParameter(2);
		    f_global[i][j][l] ->SetParameter(0, d_flow_p0_seed);
		    f_global[i][j][l] ->SetParameter(1, d_flow_p1_seed);
		    f_global[i][j][l] ->SetParameter(2, d_flow_p2_seed);
		    f_global[i][j][l] ->SetParameter(3, d_flow_seed);
                    p_mMass2_invMfit_New_Bin[i][j][l]->Fit(f_global[i][j][l],"QE+","MQR",1.,1.05);
    		    //txtFile << Form("%1.2f \t%f \t%f", (pt_low_phi[i]+pt_up_phi[i])/2.0, bin_content, bin_error) << endl;
                    double d_flow2_invm = f_global[i][j][l] -> GetParameter(3);
                    double d_err_2_invm = f_global[i][j][l] -> GetParError(3);
    		    txtFile_invm <<  Form("%1.2f \t%f \t%f", (pt_low_phi[i]+pt_up_phi[i])/2.0, d_flow2_invm, d_err_2_invm) << endl;
                    p_mMass2_invMfit_New_Bin[i][j][l]->Draw();//

                    sprintf(name, "figures_phi/flowinvmass_Pt%d_Cent%d.pdf",i,j);
                    c1_invmass_invMfit[i][j][l]->SaveAs(name);
                }
            //}
        }
    } // flowinvm

    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
            {
                for(Int_t i = 0; i < pt_total_New_phi; i++)
                {
                    for(Int_t pt_bin = pt_new_bin_start[i]; pt_bin <= pt_new_bin_stop[i]; pt_bin++)
                    {
                        if(pt_bin == pt_new_bin_start[i])
                        {
                            //h_mMass_Phi2_Ana_New_Bin[i][j][0][m] = (TH1F*)h_mMass_Phi2_Ana[pt_bin][j][0][m]->Clone();
                            //h_mMass_Phi2_Ana_New_Bin[i][j][1][m] = (TH1F*)h_mMass_Phi2_Ana[pt_bin][j][1][m]->Clone();
                            h_mMass_Phi2_SE_New_Bin[i][j][l][m] = (TH1F*)h_mMass_Phi2_SE[pt_bin][j][l][m]->Clone();
                            h_mMass_Phi2_ME_New_Bin[i][j][l][m] = (TH1F*)h_mMass_Phi2_ME[pt_bin][j][l][m]->Clone();
                        }
                        else
                        {
                            //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->Add(h_mMass_Phi2_Ana[pt_bin][j][0][m],1.0);
                            //h_mMass_Phi2_Ana_New_Bin[i][j][1][m]->Add(h_mMass_Phi2_Ana[pt_bin][j][1][m],1.0);
                            h_mMass_Phi2_SE_New_Bin[i][j][l][m]->Add(h_mMass_Phi2_SE[pt_bin][j][l][m],1.0);
                            h_mMass_Phi2_ME_New_Bin[i][j][l][m]->Add(h_mMass_Phi2_ME[pt_bin][j][l][m],1.0);
                        }
                    }
                    //h_mMass_Phi2_SM_new_bin[i][j][m] = (TH1F*)h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->Clone();
                    //h_mMass_Phi2_SM_new_bin[i][j][m] = (TH1F*)h_mMass_Phi2_SE_New_Bin[i][j][0][m]->Clone();
                    //Int_t bin2_SE_start = h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->FindBin(1.04);
                    //Int_t bin2_SE_stop  = h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->FindBin(1.05);
                    //Float_t Inte2_SE = h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->Integral(bin2_SE_start,bin2_SE_stop);
                    h_mMass_Phi2_SM_new_bin[i][j][l][m] = (TH1F*)h_mMass_Phi2_SE_New_Bin[i][j][l][m]->Clone();
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m] = (TH1F*)h_mMass_Phi2_SE_New_Bin[i][j][l][m]->Clone();
                    Int_t bin2_SE_start = h_mMass_Phi2_SE_New_Bin[i][j][l][m]->FindBin(1.04);
                    Int_t bin2_SE_stop  = h_mMass_Phi2_SE_New_Bin[i][j][l][m]->FindBin(1.05);
                    Float_t Inte2_SE = h_mMass_Phi2_SE_New_Bin[i][j][l][m]->Integral(bin2_SE_start,bin2_SE_stop);
                    //Int_t bin2_ME_start = h_mMass_Phi2_Ana_New_Bin[i][j][1][m]->FindBin(1.04);
                    //Int_t bin2_ME_stop  = h_mMass_Phi2_Ana_New_Bin[i][j][1][m]->FindBin(1.05);
                    //Float_t Inte2_ME = h_mMass_Phi2_Ana_New_Bin[i][j][1][m]->Integral(bin2_ME_start,bin2_ME_stop);
                    Int_t bin2_ME_start = h_mMass_Phi2_ME_New_Bin[i][j][l][m]->FindBin(1.04);
                    Int_t bin2_ME_stop  = h_mMass_Phi2_ME_New_Bin[i][j][l][m]->FindBin(1.05);
                    Float_t Inte2_ME = h_mMass_Phi2_ME_New_Bin[i][j][l][m]->Integral(bin2_ME_start,bin2_ME_stop);
                    //h_mMass_Phi2_Ana_New_Bin[i][j][1][m]->Scale(Inte2_SE/Inte2_ME);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Scale(Inte2_SE/Inte2_ME);
                    h_mMass_Phi2_SM_new_bin[i][j][l][m]->Scale(Inte2_SE/Inte2_ME);

                    h_mMass_Phi2_SM_new_bin[i][j][l][m]->Add(h_mMass_Phi2_SE_New_Bin[i][j][l][m],-1.0);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Add(h_mMass_Phi2_ME_New_Bin[i][j][l][m],-1.0);

                    sprintf(name, "c1_invmass_Pt%d_Cent%d_phi_Psi%d",i,j,m);
                    c1_invmass[i][j][l][m] = new TCanvas(name, name, 600, 500);
                    c1_invmass[i][j][l][m]->cd();
                    gPad->SetTopMargin(0.1);
                    gPad->SetRightMargin(0.1);
                    gPad->SetLeftMargin(0.2);
                    gPad->SetBottomMargin(0.2);
                    //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->SetMarkerColor(1);
                    //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->SetNdivisions(505,"X");
                    //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->SetNdivisions(505,"Y");
                    //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->SetYTitle("Counts");
                    //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->SetXTitle("m_{inv} (GeV/c^{2})");
                    //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->GetXaxis()->CenterTitle();
                    //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->GetYaxis()->CenterTitle();
                    //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->SetLineColor(1);//
                    //h_mMass_Phi2_Ana_New_Bin[i][j][0][m]->Draw("HIST");//
                    //h_mMass_Phi2_Ana_New_Bin[i][j][1][m]->SetMarkerColor(2);
                    //h_mMass_Phi2_Ana_New_Bin[i][j][1][m]->SetLineColor(2);
                    //h_mMass_Phi2_Ana_New_Bin[i][j][1][m]->SetFillColor(3002);
                    //h_mMass_Phi2_Ana_New_Bin[i][j][1][m]->Draw("HIST SAME");
                    //h_mMass_Phi2_SM_new_bin[i][j][m]->SetLineColor(4);
                    //h_mMass_Phi2_SM_new_bin[i][j][m]->SetFillStyle(3003);
                    //h_mMass_Phi2_SM_new_bin[i][j][m]->Draw("HIST SAME");
                    h_mMass_Phi2_SE_New_Bin[i][j][l][m]->SetMarkerColor(1);
                    h_mMass_Phi2_SE_New_Bin[i][j][l][m]->SetNdivisions(505,"X");
                    h_mMass_Phi2_SE_New_Bin[i][j][l][m]->SetNdivisions(505,"Y");
                    h_mMass_Phi2_SE_New_Bin[i][j][l][m]->SetYTitle("Counts");
                    h_mMass_Phi2_SE_New_Bin[i][j][l][m]->SetXTitle("m_{inv} (GeV/c^{2})");
                    h_mMass_Phi2_SE_New_Bin[i][j][l][m]->GetXaxis()->CenterTitle();
                    h_mMass_Phi2_SE_New_Bin[i][j][l][m]->GetYaxis()->CenterTitle();
                    h_mMass_Phi2_SE_New_Bin[i][j][l][m]->SetLineColor(1);//
                    h_mMass_Phi2_SE_New_Bin[i][j][l][m]->Draw("HIST");//
                    h_mMass_Phi2_ME_New_Bin[i][j][l][m]->SetMarkerColor(2);
                    h_mMass_Phi2_ME_New_Bin[i][j][l][m]->SetLineColor(2);
                    h_mMass_Phi2_ME_New_Bin[i][j][l][m]->SetFillColor(3002);
                    h_mMass_Phi2_ME_New_Bin[i][j][l][m]->Draw("HIST SAME");
                    h_mMass_Phi2_SM_new_bin[i][j][l][m]->SetLineColor(4);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetLineColor(4);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetFillStyle(3003);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Draw("HIST SAME");
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->SetFillStyle(3003);
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->Draw("HIST SAME");
                    //sprintf(name, "figures_phi/invmass_Pt%d_Cent%d_phi_Psi%d.pdf",i,j,m);
                    //c1_invmass[i][j][l][m]->SaveAs(name);
                }
            }
        }
    }
    cout << "test 2 ! " <<  endl;



    // merged original pT bins to get new pT bins
    //TH1F *h_mMass_Phi2_SM_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    //TH1F *h_mMass_Phi3_SM_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
            {
                for(Int_t i = 0; i < pt_total_New_phi; i++)
                {
                    for(Int_t pt_bin = pt_new_bin_start[i]; pt_bin <= pt_new_bin_stop[i]; pt_bin++)
                    {
                        if(pt_bin == pt_new_bin_start[i])
                        {
                            h_mMass_Phi2_SM_New_Bin[i][j][l][m] = (TH1F*)h_mMass_Phi2_SM[pt_bin][j][l][m]->Clone();
                            //h_mMass_Phi2_SM_new_bin[i][j][l][m] = (TH1F*)h_mMass_Phi2_SM[pt_bin][j][l][m]->Clone();
                            //h_mMass_Phi2_SM_new_bin[i][j][m] = (TH1F*)h_mMass_Phi2_SM[pt_bin][j][m]->Clone();
                            // h_mMass_Phi3_SM_New_Bin[i][j][m] = (TH1F*)h_mMass_Phi3_SM[pt_bin][j][m]->Clone();
                        }
                        else
                        {
                            h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Add(h_mMass_Phi2_SM[pt_bin][j][l][m],1.0);
                            //h_mMass_Phi2_SM_new_bin[i][j][l][m]->Add(h_mMass_Phi2_SM[pt_bin][j][l][m],1.0);
                            // h_mMass_Phi3_SM_New_Bin[i][j][m]->Add(h_mMass_Phi3_SM[pt_bin][j][m],1.0);
                        }
                    }
                    //if(i >= 6)
                    {
                        h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Rebin(4);
                        //h_mMass_Phi2_SM_new_bin[i][j][l][m]->Rebin(4);
                        // h_mMass_Phi3_SM_New_Bin[i][j][m]->Rebin(4);
                    }
                }
            }
        }
    }
    cout << "test 3 ! " <<  endl;

    // Poly+BW Fit for merged phi-psi bin
    TH1F *h_mMass_Phi2_SM_total[pt_total_New_phi][Centrality_total];
    // TH1F *h_mMass_Phi3_SM_total[pt_total_New_phi][Centrality_total];

    TF1  *f_PolyBW_Phi2_total[pt_total_New_phi][Centrality_total];
    // TF1  *f_PolyBW_Phi3_total[pt_total_New_phi][Centrality_total];

    Float_t ParFit2_total[pt_total_New_phi][Centrality_total][5];
    // Float_t ParFit3_total[pt_total_New_phi][Centrality_total][5];

    for(Int_t i = 0; i < pt_total_New_phi; i++) // pt bin
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    TString HistName;
                    if(m == 0)
                    {
                        HistName = Form("pt_%d_Centrality_%d_Total_2nd_Phi_SE",i,j);
                        h_mMass_Phi2_SM_total[i][j] = (TH1F*)h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Clone(HistName.Data());
                        //h_mMass_Phi2_SM_total[i][j] = (TH1F*)h_mMass_Phi2_SM_new_bin[i][j][l][m]->Clone(HistName.Data());
                        HistName = "f_"+HistName;
                        f_PolyBW_Phi2_total[i][j] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
                        for(Int_t n_par = 0; n_par < 5; n_par++)
                        {
                            f_PolyBW_Phi2_total[i][j]->ReleaseParameter(n_par);
                        }
                        f_PolyBW_Phi2_total[i][j]->SetParameter(0,1.019);
                        f_PolyBW_Phi2_total[i][j]->SetParLimits(0,1.014,1.024);
                        f_PolyBW_Phi2_total[i][j]->SetParameter(1,0.0055);
                        f_PolyBW_Phi2_total[i][j]->SetParameter(2,10000);
                        f_PolyBW_Phi2_total[i][j]->SetParameter(3,-6000);
                        f_PolyBW_Phi2_total[i][j]->SetParameter(4,0.5);
                        f_PolyBW_Phi2_total[i][j]->SetRange(BW_Start,BW_Stop);

                        /*
                        HistName = Form("pt_%d_Centrality_%d_Total_3rd_Phi_SE",i,j);
                        // h_mMass_Phi3_SM_total[i][j] = (TH1F*)h_mMass_Phi3_SM_New_Bin[i][j][m]->Clone(HistName.Data());
                        HistName = "f_"+HistName;
                        f_PolyBW_Phi3_total[i][j] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
                        for(Int_t n_par = 0; n_par < 4; n_par++)
                        {
                            f_PolyBW_Phi3_total[i][j]->ReleaseParameter(n_par);
                        }
                        f_PolyBW_Phi3_total[i][j]->SetParameter(0,1.019);
                        f_PolyBW_Phi3_total[i][j]->SetParLimits(0,1.014,1.024);
                        f_PolyBW_Phi3_total[i][j]->SetParameter(1,0.0055);
                        f_PolyBW_Phi3_total[i][j]->SetParameter(2,10000);
                        f_PolyBW_Phi3_total[i][j]->SetParameter(3,-6000);
                        f_PolyBW_Phi3_total[i][j]->SetParameter(4,0.5);
                        f_PolyBW_Phi3_total[i][j]->SetRange(BW_Start,BW_Stop);
                        */
                    }
                    else
                    {
                        h_mMass_Phi2_SM_total[i][j]->Add(h_mMass_Phi2_SM_New_Bin[i][j][l][m],1.0);
                        //h_mMass_Phi2_SM_total[i][j]->Add(h_mMass_Phi2_SM_new_bin[i][j][l][m],1.0);
                        // h_mMass_Phi3_SM_total[i][j]->Add(h_mMass_Phi3_SM_New_Bin[i][j][m],1.0);
                    }
                }
                h_mMass_Phi2_SM_total[i][j]->Fit(f_PolyBW_Phi2_total[i][j],"NQRN");
                // h_mMass_Phi3_SM_total[i][j]->Fit(f_PolyBW_Phi3_total[i][j],"NQRN");
                for(Int_t n_par = 0; n_par < 5; n_par++)
                {
                    ParFit2_total[i][j][n_par] = f_PolyBW_Phi2_total[i][j]->GetParameter(n_par);
                    // ParFit3_total[i][j][n_par] = f_PolyBW_Phi3_total[i][j]->GetParameter(n_par);
                }
            }
        }
    }
    cout << "test 4 ! " <<  endl;


    // Draw the merged histogram and fit line
    TCanvas *c_mMass_Phi2_SM_total[Centrality_total];
    //TCanvas *c_mMass_Phi3_SM_total[Centrality_total];
    TF1 *f_Poly_Phi2_total[pt_total_New_phi][Centrality_total];
    //TF1 *f_Poly_Phi3_total[pt_total_New_phi][Centrality_total];
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            TString CanName;

            CanName = Form("c_Phi2_SM_Total_Centrality_%d",j);
            c_mMass_Phi2_SM_total[j] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,600);
            c_mMass_Phi2_SM_total[j]->Divide(3,2);
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                c_mMass_Phi2_SM_total[j]->cd(i+1);
                c_mMass_Phi2_SM_total[j]->cd(i+1)->SetLeftMargin(0.2);
                c_mMass_Phi2_SM_total[j]->cd(i+1)->SetBottomMargin(0.2);
                c_mMass_Phi2_SM_total[j]->cd(i+1)->SetTicks(1,1);
                c_mMass_Phi2_SM_total[j]->cd(i+1)->SetGrid(0,0);
                h_mMass_Phi2_SM_total[i][j]->SetNdivisions(505,"X");
                h_mMass_Phi2_SM_total[i][j]->SetNdivisions(505,"Y");
                h_mMass_Phi2_SM_total[i][j]->SetYTitle("Counts");
                h_mMass_Phi2_SM_total[i][j]->SetXTitle("m_{inv} (GeV/c^{2})");
                h_mMass_Phi2_SM_total[i][j]->GetXaxis()->CenterTitle();
                h_mMass_Phi2_SM_total[i][j]->GetYaxis()->CenterTitle();
                h_mMass_Phi2_SM_total[i][j]->Draw("PE");
                f_PolyBW_Phi2_total[i][j]->SetLineColor(2);
                f_PolyBW_Phi2_total[i][j]->SetLineWidth(1);
                f_PolyBW_Phi2_total[i][j]->Draw("l same");
                CanName = Form("f_poly2_pt_%d_Centrality_%d",i,j);
                f_Poly_Phi2_total[i][j] = new TF1(CanName.Data(),Poly,BW_Start,BW_Stop,2);
                f_Poly_Phi2_total[i][j]->FixParameter(0,f_PolyBW_Phi2_total[i][j]->GetParameter(3));
                f_Poly_Phi2_total[i][j]->FixParameter(1,f_PolyBW_Phi2_total[i][j]->GetParameter(4));
                f_Poly_Phi2_total[i][j]->SetLineWidth(1);
                f_Poly_Phi2_total[i][j]->SetLineColor(4);
                f_Poly_Phi2_total[i][j]->Draw("l same");
                //h_mMass_Phi2_SM_total[i][j]->Add(f_Poly_Phi2_total[i][j],-1.0);
            }
            /*
            CanName = Form("c_Phi3_SM_Total_Centrality_%d",j);
            c_mMass_Phi3_SM_total[j] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,600);
            c_mMass_Phi3_SM_total[j]->Divide(3,2);
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                c_mMass_Phi3_SM_total[j]->cd(i+1);
                c_mMass_Phi3_SM_total[j]->cd(i+1)->SetLeftMargin(0.2);
                c_mMass_Phi3_SM_total[j]->cd(i+1)->SetBottomMargin(0.2);
                c_mMass_Phi3_SM_total[j]->cd(i+1)->SetTicks(1,1);
                c_mMass_Phi3_SM_total[j]->cd(i+1)->SetGrid(0,0);
                h_mMass_Phi3_SM_total[i][j]->SetNdivisions(505,"X");
                h_mMass_Phi3_SM_total[i][j]->SetNdivisions(505,"Y");
                h_mMass_Phi3_SM_total[i][j]->SetYTitle("Counts");
                h_mMass_Phi3_SM_total[i][j]->SetXTitle("m_{inv} (GeV/c^{2})");
                h_mMass_Phi3_SM_total[i][j]->GetXaxis()->CenterTitle();
                h_mMass_Phi3_SM_total[i][j]->GetYaxis()->CenterTitle();
                h_mMass_Phi3_SM_total[i][j]->Draw("PE");
                f_PolyBW_Phi3_total[i][j]->SetLineColor(2);
                f_PolyBW_Phi3_total[i][j]->SetLineWidth(1);
                f_PolyBW_Phi3_total[i][j]->Draw("l same");
                CanName = Form("f_poly2_pt_%d_Centrality_%d_3rd",i,j);
                f_Poly_Phi3_total[i][j] = new TF1(CanName.Data(),Poly,BW_Start,BW_Stop,2);
                f_Poly_Phi3_total[i][j]->FixParameter(0,f_PolyBW_Phi3_total[i][j]->GetParameter(3));
                f_Poly_Phi3_total[i][j]->FixParameter(1,f_PolyBW_Phi3_total[i][j]->GetParameter(4));
                f_Poly_Phi3_total[i][j]->SetLineWidth(1);
                f_Poly_Phi3_total[i][j]->SetLineColor(4);
                f_Poly_Phi3_total[i][j]->Draw("l same");
            }
            */
            CanName = Form("figures_phi/c_Phi2_SM_Total_Centrality_%d.pdf",j);
            //c_mMass_Phi2_SM_total[j]->SaveAs(CanName.Data());
            /*
            CanName = Form("figures_phi/c_Phi3_SM_Total_Centrality_%d.pdf",j);
            c_mMass_Phi3_SM_total[j]->SaveAs(CanName.Data());
            */
        }
    }
    cout << "test 5 ! " <<  endl;


    // Poly+BW Fit for phi-psi bin
    TF1 *f_PolyBW_Phi2[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    //TF1 *f_PolyBW_Phi3[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    for(Int_t i = 0; i < pt_total_New_phi; i++) // pt bin
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    TString HistName;

                    HistName = Form("f_pt_%d_Centrality_%d_phi_Psi_%d_2nd_Phi_SE",i,j,m);
                    f_PolyBW_Phi2[i][j][l][m] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
                    for(Int_t n_par = 0; n_par < 5; n_par++)
                    {
                        f_PolyBW_Phi2[i][j][l][m]->ReleaseParameter(n_par);
                    }
                    f_PolyBW_Phi2[i][j][l][m]->FixParameter(0,ParFit2_total[i][j][0]);
                    f_PolyBW_Phi2[i][j][l][m]->FixParameter(1,ParFit2_total[i][j][1]);
                    f_PolyBW_Phi2[i][j][l][m]->SetParameter(2,ParFit2_total[i][j][2]/7.0);
                    f_PolyBW_Phi2[i][j][l][m]->SetParameter(3,ParFit2_total[i][j][3]);
                    f_PolyBW_Phi2[i][j][l][m]->SetParameter(4,ParFit2_total[i][j][4]);
                    f_PolyBW_Phi2[i][j][l][m]->SetRange(BW_Start,BW_Stop);
                    /*
                    HistName = Form("f_pt_%d_Centrality_%d_phi_Psi_%d_3rd_Phi_SE",i,j,m);
                    f_PolyBW_Phi3[i][j][m] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
                    for(Int_t n_par = 0; n_par < 5; n_par++)
                    {
                        f_PolyBW_Phi3[i][j][m]->ReleaseParameter(n_par);
                    }
                    f_PolyBW_Phi3[i][j][m]->FixParameter(0,ParFit3_total[i][j][0]);
                    f_PolyBW_Phi3[i][j][m]->FixParameter(1,ParFit3_total[i][j][1]);
                    f_PolyBW_Phi3[i][j][m]->SetParameter(2,ParFit3_total[i][j][2]/7.0);
                    f_PolyBW_Phi3[i][j][m]->SetParameter(3,ParFit3_total[i][j][3]);
                    f_PolyBW_Phi3[i][j][m]->SetParameter(4,ParFit3_total[i][j][4]);
                    f_PolyBW_Phi3[i][j][m]->SetRange(BW_Start,BW_Stop);
                    */
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Fit(f_PolyBW_Phi2[i][j][l][m],"RMN");
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->Fit(f_PolyBW_Phi2[i][j][l][m],"RMN");
                    // h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Fit(f_PolyBW_Phi3[i][j][m],"RMN");
                }
            }
        }
    }
    cout << "test 6 ! " <<  endl;

    // subtract linear background from signal
    TF1 *f_Poly_Phi2[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    // TF1 *f_Poly_Phi3[pt_total_New_phi][Centrality_total][Phi_Psi_total];
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    TString FuncName;

                    FuncName = Form("f_poly2_pt_%d_Centrality_%d_phi_Psi_%d",i,j,m);
                    f_Poly_Phi2[i][j][l][m] = new TF1(FuncName.Data(),Poly,BW_Start,BW_Stop,2);
                    f_Poly_Phi2[i][j][l][m]->FixParameter(0,f_PolyBW_Phi2[i][j][l][m]->GetParameter(3));
                    f_Poly_Phi2[i][j][l][m]->FixParameter(1,f_PolyBW_Phi2[i][j][l][m]->GetParameter(4));
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Add(f_Poly_Phi2[i][j][l][m],-1.0);
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->Add(f_Poly_Phi2[i][j][l][m],-1.0);
                    /*
                    FuncName = Form("f_poly3_pt_%d_Centrality_%d_phi_Psi_%d",i,j,m);
                    f_Poly_Phi3[i][j][m] = new TF1(FuncName.Data(),Poly,BW_Start,BW_Stop,2);
                    f_Poly_Phi3[i][j][m]->FixParameter(0,f_PolyBW_Phi3[i][j][m]->GetParameter(3));
                    f_Poly_Phi3[i][j][m]->FixParameter(1,f_PolyBW_Phi3[i][j][m]->GetParameter(4));
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->Add(f_Poly_Phi3[i][j][m],-1.0);
                    */
                }
            }
        }
    }
    cout << "test 7 ! " <<  endl;

    // do gaussian fit and extract width of Phi from merged phi-Psi distribution
    // do breit wigner fit to get fit parameter for counting
    TH1F *h_mMass_Phi2_Sig_total[pt_total_New_phi][Centrality_total];
    // TH1F *h_mMass_Phi3_Sig_total[pt_total_New_phi][Centrality_total];

    TF1 *f_gauss_2[pt_total_New_phi][Centrality_total];
    TF1 *f_gauss_3[pt_total_New_phi][Centrality_total];

    //Float_t gaussSig = 0.60;
    //Float_t gaussSig = 1.5; // dchen
    //Float_t gaussSig = 3.; // dchen
    Float_t ParGaus_2[pt_total_New_phi][Centrality_total][2];
    Float_t ParGaus_3[pt_total_New_phi][Centrality_total][2];

    TF1 *f_bw_2[pt_total_New_phi][Centrality_total];
    TF1 *f_bw_3[pt_total_New_phi][Centrality_total];

    Float_t ParBW_2[pt_total_New_phi][Centrality_total][2];
    Float_t ParBW_3[pt_total_New_phi][Centrality_total][2];
    for(Int_t i = 0; i < pt_total_New_phi; i++) // pt bin
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                TString HistName;
                TString FuncName_g;
                TString FuncName_bw;

                HistName = Form("Sig_pt_%d_Centality_%d_2nd",i,j);
                FuncName_g = Form("f_gauss_pt_%d_Centrality_%d_2nd",i,j);
                FuncName_bw = Form("f_bw_pt_%d_Centrality_%d_2nd",i,j);
                for(Int_t m = 0; m < Phi_Psi_total; m++)
                {
                    if(m == 0)
                    {
                        //h_mMass_Phi2_Sig_total[i][j] = (TH1F*)h_mMass_Phi2_SM_new_bin[i][j][l][m]->Clone(HistName.Data());
                        h_mMass_Phi2_Sig_total[i][j] = (TH1F*)h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Clone(HistName.Data());
                        h_mMass_Phi2_Sig_total[i][j]->SetTitle("");
                        f_gauss_2[i][j] = new TF1(FuncName_g.Data(),Gaussion,BW_Start,BW_Stop,3); f_gauss_2[i][j]->SetParameter(0,ParFit2_total[i][j][0]);
                        //	    f_gauss_2[i][j]->SetParameter(1,ParFit2_total[i][j][1]);
                        f_gauss_2[i][j]->SetParameter(1,0.07);
                        f_gauss_2[i][j]->SetParameter(2,1000);
                        //f_gauss_2[i][j]->SetRange(ParFit2_total[i][j][0]-(gaussSig+2)*ParFit2_total[i][j][1],ParFit2_total[i][j][0]+(gaussSig+2)*ParFit2_total[i][j][1]);
                        f_gauss_2[i][j]->SetRange(ParFit2_total[i][j][0]-gaussSig*ParFit2_total[i][j][1],ParFit2_total[i][j][0]+gaussSig*ParFit2_total[i][j][1]); // dchen

                        f_bw_2[i][j] = new TF1(FuncName_bw.Data(),BreitWigner,BW_Start,BW_Stop,3);
                        f_bw_2[i][j]->SetParameter(0,ParFit2_total[i][j][0]);
                        f_bw_2[i][j]->SetParameter(1,ParFit2_total[i][j][1]);
                        f_bw_2[i][j]->SetParameter(2,1000);
                        f_bw_2[i][j]->SetRange(BW_Start,BW_Stop);
                    }
                    else
                    {
                        h_mMass_Phi2_Sig_total[i][j]->Add(h_mMass_Phi2_SM_New_Bin[i][j][l][m],1.0);
                        //h_mMass_Phi2_Sig_total[i][j]->Add(h_mMass_Phi2_SM_new_bin[i][j][l][m],1.0);
                    }
                }
                h_mMass_Phi2_Sig_total[i][j]->Fit(f_gauss_2[i][j],"MQR");
                if(i == 8 && j ==3) {h_mMass_Phi2_Sig_total[i][j]->Fit(f_gauss_2[i][j],"MQR");}
                ParGaus_2[i][j][0] = f_gauss_2[i][j]->GetParameter(0);
                ParGaus_2[i][j][1] = f_gauss_2[i][j]->GetParameter(1);
                if(i ==8 && j==3)
                {
                    cout << " ParGaus_2 1 = " << ParGaus_2[i][j][0] << endl;
                    cout << " ParGaus_2 2 = " << ParGaus_2[i][j][1] << endl;
                    //ParGaus_2[i][j][1] = -1*ParGaus_2[i][j][1];
                }

                h_mMass_Phi2_Sig_total[i][j]->Fit(f_bw_2[i][j],"MQNR");
                // f_bw_2[i][j]->SetLineColor()
                ParBW_2[i][j][0] = f_bw_2[i][j]->GetParameter(0);
                ParBW_2[i][j][1] = f_bw_2[i][j]->GetParameter(1);

                /*
                HistName = Form("Sig_pt_%d_Centality_%d_3rd",i,j);
                FuncName_g = Form("f_gauss_pt_%d_Centrality_%d_3rd",i,j);
                FuncName_bw = Form("f_bw_pt_%d_Centrality_%d_3rd",i,j);
                for(Int_t m = 0; m < Phi_Psi_total; m++)
                {
                    if(m == 0)
                    {
                        h_mMass_Phi3_Sig_total[i][j] = (TH1F*)h_mMass_Phi3_SM_New_Bin[i][j][m]->Clone(HistName.Data());
                        h_mMass_Phi3_Sig_total[i][j]->SetTitle("");
                        f_gauss_3[i][j] = new TF1(FuncName_g.Data(),Gaussion,BW_Start,BW_Stop,3);
                        f_gauss_3[i][j]->SetParameter(0,ParFit3_total[i][j][0]);
                        //	    f_gauss_3[i][j]->SetParameter(1,ParFit3_total[i][j][1]);
                        f_gauss_2[i][j]->SetParameter(1,0.07);
                        f_gauss_3[i][j]->SetParameter(2,1000);
                        f_gauss_3[i][j]->SetRange(ParFit3_total[i][j][0]-gaussSig*ParFit3_total[i][j][1],ParFit3_total[i][j][0]+gaussSig*ParFit3_total[i][j][1]);

                        f_bw_3[i][j] = new TF1(FuncName_bw.Data(),BreitWigner,BW_Start,BW_Stop,3);
                        f_bw_3[i][j]->SetParameter(0,ParFit3_total[i][j][0]);
                        f_bw_3[i][j]->SetParameter(1,ParFit3_total[i][j][1]);
                        f_bw_3[i][j]->SetParameter(2,1000);
                        f_bw_3[i][j]->SetRange(BW_Start,BW_Stop);
                    }
                    else
                    {
                        h_mMass_Phi3_Sig_total[i][j]->Add(h_mMass_Phi3_SM_New_Bin[i][j][m],1.0);
                    }
                }
                h_mMass_Phi3_Sig_total[i][j]->Fit(f_gauss_3[i][j],"MQNR");
                ParGaus_3[i][j][0] = f_gauss_3[i][j]->GetParameter(0);
                ParGaus_3[i][j][1] = f_gauss_3[i][j]->GetParameter(1);

                h_mMass_Phi3_Sig_total[i][j]->Fit(f_bw_3[i][j],"MQNR");
                ParBW_3[i][j][0] = f_bw_3[i][j]->GetParameter(0);
                ParBW_3[i][j][1] = f_bw_3[i][j]->GetParameter(1);
                */
            }
        }
    }
    cout << "test 8 ! " <<  endl;

    // Draw InvMass distribution for all pt bin with Gaussian fit
    TCanvas *c2_phi_Psi[pt_total_New_phi][Centrality_total];
    // TCanvas *c3_phi_Psi[pt_total_New_phi][Centrality_total];
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                TString CanName;

                CanName = Form("c2_phi_Psi_pt_%d_Centrality_%d",i,j);
                c2_phi_Psi[i][j] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
                c2_phi_Psi[i][j]->Divide(3,3);
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    c2_phi_Psi[i][j]->cd(m+1);
                    c2_phi_Psi[i][j]->cd(m+1)->SetLeftMargin(0.2);
                    c2_phi_Psi[i][j]->cd(m+1)->SetBottomMargin(0.2);
                    c2_phi_Psi[i][j]->cd(m+1)->SetTicks(1,1);
                    c2_phi_Psi[i][j]->cd(m+1)->SetGrid(0,0);
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->SetTitle("");
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->SetNdivisions(505,"X");
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->SetNdivisions(505,"Y");
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetXaxis()->SetLabelSize(0.05);
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetYaxis()->SetLabelSize(0.05);
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetYaxis()->SetTitle("Counts");
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetXaxis()->SetTitleSize(0.06);
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetYaxis()->SetTitleSize(0.06);
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetYaxis()->SetTitleOffset(1.7);
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetXaxis()->CenterTitle();
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetYaxis()->CenterTitle();
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->Draw("pE");
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetTitle("");
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetNdivisions(505,"X");
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetNdivisions(505,"Y");
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetXaxis()->SetLabelSize(0.05);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->SetLabelSize(0.05);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->SetTitle("Counts");
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetXaxis()->SetTitleSize(0.06);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->SetTitleSize(0.06);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->SetTitleOffset(1.7);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetXaxis()->CenterTitle();
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->CenterTitle();
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Draw("pE");
                    Float_t x1 = ParGaus_2[i][j][0] - nSigmaPhi*ParGaus_2[i][j][1];
                    Float_t x2 = ParGaus_2[i][j][0] + nSigmaPhi*ParGaus_2[i][j][1];


					cout << "shaowei  = " << ParGaus_2[i][j][1] << endl;

                    //Float_t y = h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetBinContent(h_mMass_Phi2_SM_new_bin[i][j][l][m]->FindBin(ParGaus_2[i][j][0]));
                    Float_t y = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinContent(h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_2[i][j][0]));
                    PlotLine(x1,x1,0,y,4,2,2);
                    PlotLine(x2,x2,0,y,4,2,2);
                    PlotLine(0.98,1.05,0,0,1,2,2);
                }
                c2_phi_Psi[i][j]->cd(8);
                c2_phi_Psi[i][j]->cd(8)->SetLeftMargin(0.2);
                c2_phi_Psi[i][j]->cd(8)->SetBottomMargin(0.2);
                c2_phi_Psi[i][j]->cd(8)->SetTicks(1,1);
                c2_phi_Psi[i][j]->cd(8)->SetGrid(0,0);
                h_mMass_Phi2_Sig_total[i][j]->SetTitle("Merged Inv. Mass distribution");
                h_mMass_Phi2_Sig_total[i][j]->SetNdivisions(505,"X");
                h_mMass_Phi2_Sig_total[i][j]->SetNdivisions(505,"Y");
                h_mMass_Phi2_Sig_total[i][j]->GetXaxis()->SetLabelSize(0.05);
                h_mMass_Phi2_Sig_total[i][j]->GetYaxis()->SetLabelSize(0.05);
                h_mMass_Phi2_Sig_total[i][j]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
                h_mMass_Phi2_Sig_total[i][j]->GetYaxis()->SetTitle("Counts");
                h_mMass_Phi2_Sig_total[i][j]->GetXaxis()->SetTitleSize(0.06);
                h_mMass_Phi2_Sig_total[i][j]->GetYaxis()->SetTitleSize(0.06);
                h_mMass_Phi2_Sig_total[i][j]->GetYaxis()->SetTitleOffset(1.7);
                h_mMass_Phi2_Sig_total[i][j]->GetXaxis()->CenterTitle();
                h_mMass_Phi2_Sig_total[i][j]->GetYaxis()->CenterTitle();
                h_mMass_Phi2_Sig_total[i][j]->Draw("pE");
                f_gauss_2[i][j]->SetLineColor(4);
                f_gauss_2[i][j]->Draw("l same");
                PlotLine(0.98,1.05,0,0,1,2,2);
                /*
                CanName = Form("c3_phi_Psi_pt_%d_Centrality_%d",i,j);
                c3_phi_Psi[i][j] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
                c3_phi_Psi[i][j]->Divide(3,3);
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    c3_phi_Psi[i][j]->cd(m+1);
                    c3_phi_Psi[i][j]->cd(m+1)->SetLeftMargin(0.2);
                    c3_phi_Psi[i][j]->cd(m+1)->SetBottomMargin(0.2);
                    c3_phi_Psi[i][j]->cd(m+1)->SetTicks(1,1);
                    c3_phi_Psi[i][j]->cd(m+1)->SetGrid(0,0);
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->SetTitle("");
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->SetNdivisions(505,"X");
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->SetNdivisions(505,"Y");
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->GetXaxis()->SetLabelSize(0.05);
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->GetYaxis()->SetLabelSize(0.05);
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->GetYaxis()->SetTitle("Counts");
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->GetXaxis()->SetTitleSize(0.06);
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->GetYaxis()->SetTitleSize(0.06);
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->GetYaxis()->SetTitleOffset(1.7);
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->GetXaxis()->CenterTitle();
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->GetYaxis()->CenterTitle();
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->Draw("pE");
                    Float_t x1 = ParGaus_3[i][j][0] - nSigmaPhi*ParGaus_3[i][j][1];
                    Float_t x2 = ParGaus_3[i][j][0] + nSigmaPhi*ParGaus_3[i][j][1];
                    Float_t y = h_mMass_Phi3_SM_New_Bin[i][j][m]->GetBinContent(h_mMass_Phi3_SM_New_Bin[i][j][m]->FindBin(ParGaus_3[i][j][0]));
                    PlotLine(x1,x1,0,y,4,2,2);
                    PlotLine(x2,x2,0,y,4,2,2);
                    PlotLine(0.98,1.05,0,0,1,2,2);
                }
                c3_phi_Psi[i][j]->cd(8);
                c3_phi_Psi[i][j]->cd(8)->SetLeftMargin(0.2);
                c3_phi_Psi[i][j]->cd(8)->SetBottomMargin(0.2);
                c3_phi_Psi[i][j]->cd(8)->SetTicks(1,1);
                c3_phi_Psi[i][j]->cd(8)->SetGrid(0,0);
                h_mMass_Phi3_Sig_total[i][j]->SetTitle("Merged Inv. Mass distribution");
                h_mMass_Phi3_Sig_total[i][j]->SetNdivisions(505,"X");
                h_mMass_Phi3_Sig_total[i][j]->SetNdivisions(505,"Y");
                h_mMass_Phi3_Sig_total[i][j]->GetXaxis()->SetLabelSize(0.05);
                h_mMass_Phi3_Sig_total[i][j]->GetYaxis()->SetLabelSize(0.05);
                h_mMass_Phi3_Sig_total[i][j]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
                h_mMass_Phi3_Sig_total[i][j]->GetYaxis()->SetTitle("Counts");
                h_mMass_Phi3_Sig_total[i][j]->GetXaxis()->SetTitleSize(0.06);
                h_mMass_Phi3_Sig_total[i][j]->GetYaxis()->SetTitleSize(0.06);
                h_mMass_Phi3_Sig_total[i][j]->GetYaxis()->SetTitleOffset(1.7);
                h_mMass_Phi3_Sig_total[i][j]->GetXaxis()->CenterTitle();
                h_mMass_Phi3_Sig_total[i][j]->GetYaxis()->CenterTitle();
                h_mMass_Phi3_Sig_total[i][j]->Draw("pE");
                f_gauss_3[i][j]->SetLineColor(2);
                //f_gauss_3[i][j]->Draw("l same");
                PlotLine(0.98,1.05,0,0,1,2,2);
                */
            }
        }
    }
    cout << "test 9 ! " <<  endl;

    // calculate total counts and errors for each phi-Psi bin by bin counting
    TH1F *h_Counts2[pt_total_New_phi][Centrality_total];
    // TH1F *h_Counts3[pt_total_New_phi][Centrality_total];
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                TString HistName;

                HistName = Form("Counts2_pt_%d_Centrality_%d",i,j);
                h_Counts2[i][j] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[0]);
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    Float_t counts = 0.0;
                    Float_t errors = 0.0;
                    Float_t bin_center = PI_max[0]/14.0+m*PI_max[0]/7.0;
                    //Int_t bin_start = h_mMass_Phi2_SM_new_bin[i][j][l][m]->FindBin(ParGaus_2[i][j][0]-nSigmaPhi*ParGaus_2[i][j][1]);
                    //Int_t bin_stop  = h_mMass_Phi2_SM_new_bin[i][j][l][m]->FindBin(ParGaus_2[i][j][0]+nSigmaPhi*ParGaus_2[i][j][1]);
                    Int_t bin_start = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_2[i][j][0]-nSigmaPhi*ParGaus_2[i][j][1]);
                    Int_t bin_stop  = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_2[i][j][0]+nSigmaPhi*ParGaus_2[i][j][1]);
                    if(i ==8 && j==3)
                    {
                        cout <<"binstart = " << bin_start << endl;
                        cout <<"binstop = " << bin_stop << endl;
                    }
                    for(Int_t bin = bin_start; bin <= bin_stop; bin++)
                    {
                        //counts += h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetBinContent(bin);
                        //errors += h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetBinError(bin)*h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetBinError(bin);
                        counts += h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinContent(bin);
                        errors += h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinError(bin)*h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinError(bin);
                    }
                        //cout << "cos = " << counts << endl;
                        //cout << "err = " << errors << endl;
                    h_Counts2[i][j]->SetBinContent(h_Counts2[i][j]->FindBin(bin_center),counts);
                    h_Counts2[i][j]->SetBinError(h_Counts2[i][j]->FindBin(bin_center),TMath::Sqrt(errors));
                }
                /*
                HistName = Form("Counts3_pt_%d_Centrality_%d",i,j);
                h_Counts3[i][j] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[1]);
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    Float_t counts = 0.0;
                    Float_t errors = 0.0;
                    Float_t bin_center = PI_max[1]/14.0+m*PI_max[1]/7.0;
                    Int_t bin_start = h_mMass_Phi3_SM_New_Bin[i][j][m]->FindBin(ParGaus_3[i][j][0]-nSigmaPhi*ParGaus_3[i][j][1]);
                    Int_t bin_stop  = h_mMass_Phi3_SM_New_Bin[i][j][m]->FindBin(ParGaus_3[i][j][0]+nSigmaPhi*ParGaus_3[i][j][1]);
                    for(Int_t bin = bin_start; bin <= bin_stop; bin++)
                    {
                        counts += h_mMass_Phi3_SM_New_Bin[i][j][m]->GetBinContent(bin);
                        errors += h_mMass_Phi3_SM_New_Bin[i][j][m]->GetBinError(bin)*h_mMass_Phi3_SM_New_Bin[i][j][m]->GetBinError(bin);
                    }
                    h_Counts3[i][j]->SetBinContent(h_Counts3[i][j]->FindBin(bin_center),counts);
                    h_Counts3[i][j]->SetBinError(h_Counts3[i][j]->FindBin(bin_center),TMath::Sqrt(errors));
                }
                */
           }
        }
    }
    cout << "test 10 ! " <<  endl;

    // Draw phi-Psi distribution with counting
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                c2_phi_Psi[i][j]->cd(9);
                c2_phi_Psi[i][j]->cd(9)->SetLeftMargin(0.2);
                c2_phi_Psi[i][j]->cd(9)->SetBottomMargin(0.2);
                c2_phi_Psi[i][j]->cd(9)->SetTicks(1,1);
                c2_phi_Psi[i][j]->cd(9)->SetGrid(0,0);
                TString Title = Form("%1.1f < p_{T} < %1.1f GeV/c",pt_low_phi[i],pt_up_phi[i]);
                h_Counts2[i][j]->SetTitle(Title.Data());
                h_Counts2[i][j]->SetNdivisions(505,"X");
                h_Counts2[i][j]->SetNdivisions(505,"Y");
                h_Counts2[i][j]->GetXaxis()->SetLabelSize(0.05);
                h_Counts2[i][j]->GetYaxis()->SetLabelSize(0.05);
                h_Counts2[i][j]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
                h_Counts2[i][j]->GetYaxis()->SetTitle("Counts");
                h_Counts2[i][j]->GetXaxis()->SetTitleSize(0.06);
                h_Counts2[i][j]->GetYaxis()->SetTitleSize(0.06);
                h_Counts2[i][j]->GetYaxis()->SetTitleOffset(1.7);
                h_Counts2[i][j]->GetXaxis()->CenterTitle();
                h_Counts2[i][j]->GetYaxis()->CenterTitle();
                Int_t bin_low = h_Counts2[i][j]->FindBin(PI_max[0]/14.0+6*PI_max[0]/7.0);
                Int_t bin_up  = h_Counts2[i][j]->FindBin(PI_max[0]/14.0+0*PI_max[0]/7.0);
                h_Counts2[i][j]->GetYaxis()->SetRangeUser(h_Counts2[i][j]->GetBinContent(bin_low)*0.5,h_Counts2[i][j]->GetBinContent(bin_up)*1.5);
                h_Counts2[i][j]->SetLineColor(2);
                h_Counts2[i][j]->SetMarkerColor(2);
                h_Counts2[i][j]->Draw("pE");
                TLegend *leg2 = new TLegend(0.5,0.7,0.88,0.8);
                leg2->SetFillColor(0);
                leg2->SetFillStyle(0);
                leg2->SetBorderSize(0);
                leg2->AddEntry(h_Counts2[i][j],"Counting","p");
                leg2->Draw("same");
                /*
                c3_phi_Psi[i][j]->cd(9);
                c3_phi_Psi[i][j]->cd(9)->SetLeftMargin(0.2);
                c3_phi_Psi[i][j]->cd(9)->SetBottomMargin(0.2);
                c3_phi_Psi[i][j]->cd(9)->SetTicks(1,1);
                c3_phi_Psi[i][j]->cd(9)->SetGrid(0,0);
                h_Counts3[i][j]->SetTitle(Title.Data());
                h_Counts3[i][j]->SetNdivisions(505,"X");
                h_Counts3[i][j]->SetNdivisions(505,"Y");
                h_Counts3[i][j]->GetXaxis()->SetLabelSize(0.05);
                h_Counts3[i][j]->GetYaxis()->SetLabelSize(0.05);
                h_Counts3[i][j]->GetXaxis()->SetTitle("#phi-#Psi_{3}");
                h_Counts3[i][j]->GetYaxis()->SetTitle("Counts");
                h_Counts3[i][j]->GetXaxis()->SetTitleSize(0.06);
                h_Counts3[i][j]->GetYaxis()->SetTitleSize(0.06);
                h_Counts3[i][j]->GetYaxis()->SetTitleOffset(1.7);
                h_Counts3[i][j]->GetXaxis()->CenterTitle();
                h_Counts3[i][j]->GetYaxis()->CenterTitle();
                bin_low = h_Counts3[i][j]->FindBin(PI_max[1]/14.0+6*PI_max[1]/7.0);
                bin_up  = h_Counts3[i][j]->FindBin(PI_max[1]/14.0+0*PI_max[1]/7.0);
                h_Counts3[i][j]->GetYaxis()->SetRangeUser(h_Counts3[i][j]->GetBinContent(bin_low)*0.8,h_Counts3[i][j]->GetBinContent(bin_up)*1.2);
                h_Counts3[i][j]->SetLineColor(2);
                h_Counts3[i][j]->SetMarkerColor(2);
                h_Counts3[i][j]->Draw("pE");
                TLegend *leg3 = new TLegend(0.5,0.7,0.88,0.8);
                leg3->SetFillColor(0);
                leg3->SetBorderSize(0);
                leg3->AddEntry(h_Counts3[i][j],"Counting","p");
                leg3->Draw("same");
                */
            }
        }
    }
    cout << "test 11 ! " <<  endl;

    // Breit Wigner fit for phi-Psi bin
    TF1 *f_bw_phi_2[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
    //TF1 *f_bw_phi_3[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

    // calculate total counts and errors for each phi-Psi bin by BW Integration
    TH1F *h_Counts2_bw[pt_total_New_phi][Centrality_total];
    // TH1F *h_Counts3_bw[pt_total_New_phi][Centrality_total];
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                TString HistName;

                HistName = Form("Counts2_pt_%d_Centrality_%d_bw",i,j);
                h_Counts2_bw[i][j] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[0]);
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    TString FuncName_bw = Form("f_bw_pt_%d_Centrality_%d_phi_Psi_%d_2nd",i,j,m);
                    f_bw_phi_2[i][j][l][m] = new TF1(FuncName_bw.Data(),BreitWigner,BW_Start,BW_Stop,3);;
                    f_bw_phi_2[i][j][l][m]->FixParameter(0,ParBW_2[i][j][0]);
                    f_bw_phi_2[i][j][l][m]->FixParameter(1,ParBW_2[i][j][1]);
                    f_bw_phi_2[i][j][l][m]->SetParameter(2,1000);
                    f_bw_phi_2[i][j][l][m]->SetRange(BW_Start,BW_Stop);
                    //h_mMass_Phi2_SM_new_bin[i][j][l][m]->Fit(f_bw_phi_2[i][j][l][m],"NMQR");
                    //Float_t bin_width = h_mMass_Phi2_SM_new_bin[i][j][l][m]->GetBinWidth(1);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Fit(f_bw_phi_2[i][j][l][m],"NMQR");
                    Float_t bin_width = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinWidth(1);
                    Float_t Inte_start = ParGaus_2[i][j][0]-nSigmaPhi*ParGaus_2[i][j][1]-0.5*bin_width;
                    Float_t Inte_stop  = ParGaus_2[i][j][0]+nSigmaPhi*ParGaus_2[i][j][1]+0.5*bin_width;
                    Float_t counts = f_bw_phi_2[i][j][l][m]->Integral(Inte_start,Inte_stop)/bin_width;
                    Float_t errors = f_bw_phi_2[i][j][l][m]->IntegralError(Inte_start,Inte_stop)/bin_width;
                    Float_t bin_center = PI_max[0]/14.0+m*PI_max[0]/7.0;
                    h_Counts2_bw[i][j]->SetBinContent(h_Counts2_bw[i][j]->FindBin(bin_center),counts);
                    h_Counts2_bw[i][j]->SetBinError(h_Counts2_bw[i][j]->FindBin(bin_center),errors);
                }
                /*
                HistName = Form("Counts3_pt_%d_Centrality_%d_bw",i,j);
                h_Counts3_bw[i][j] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[1]);
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    TString FuncName_bw = Form("f_bw_pt_%d_Centrality_%d_phi_Psi_%d_3rd",i,j,m);
                    f_bw_phi_3[i][j][m] = new TF1(FuncName_bw.Data(),BreitWigner,BW_Start,BW_Stop,3);;
                    f_bw_phi_3[i][j][m]->FixParameter(0,ParBW_3[i][j][0]);
                    f_bw_phi_3[i][j][m]->FixParameter(1,ParBW_3[i][j][1]);
                    f_bw_phi_3[i][j][m]->SetParameter(2,1000);
                    f_bw_phi_3[i][j][m]->SetRange(BW_Start,BW_Stop);
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->Fit(f_bw_phi_3[i][j][m],"NMQR");
                    Float_t bin_width = h_mMass_Phi3_SM_New_Bin[i][j][m]->GetBinWidth(1);
                    Float_t Inte_start = ParGaus_3[i][j][0]-nSigmaPhi*ParGaus_3[i][j][1]-0.5*bin_width;
                    Float_t Inte_stop  = ParGaus_3[i][j][0]+nSigmaPhi*ParGaus_3[i][j][1]+0.5*bin_width;
                    Float_t counts = f_bw_phi_3[i][j][m]->Integral(Inte_start,Inte_stop)/bin_width;
                    Float_t errors = f_bw_phi_3[i][j][m]->IntegralError(Inte_start,Inte_stop)/bin_width;
                    Float_t bin_center = PI_max[1]/14.0+m*PI_max[1]/7.0;
                    h_Counts3_bw[i][j]->SetBinContent(h_Counts3_bw[i][j]->FindBin(bin_center),counts);
                    h_Counts3_bw[i][j]->SetBinError(h_Counts3_bw[i][j]->FindBin(bin_center),errors);
                }
                */
            }
        }
    }
    cout << "test 12 ! " <<  endl;

    // Draw Breit Wigner fit for InvMass distribution of all pt bin
    TCanvas *c2_phi_Psi_bw[pt_total_New_phi][Centrality_total];
    // TCanvas *c3_phi_Psi_bw[pt_total_New_phi][Centrality_total];
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                TString CanName;

                CanName = Form("c2_phi_Psi_pt_%d_Centrality_%d_bw",i,j);
                c2_phi_Psi_bw[i][j] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
                c2_phi_Psi_bw[i][j]->Divide(3,3);
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    c2_phi_Psi_bw[i][j]->cd(m+1);
                    c2_phi_Psi_bw[i][j]->cd(m+1)->SetLeftMargin(0.2);
                    c2_phi_Psi_bw[i][j]->cd(m+1)->SetBottomMargin(0.2);
                    c2_phi_Psi_bw[i][j]->cd(m+1)->SetTicks(1,1);
                    c2_phi_Psi_bw[i][j]->cd(m+1)->SetGrid(0,0);
                    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Draw("pE");
                    f_bw_phi_2[i][j][l][m]->SetLineColor(2); // red fit on the invMass
                    f_bw_phi_2[i][j][l][m]->Draw("l Same");
                    Float_t x1 = ParGaus_2[i][j][0] - nSigmaPhi*ParGaus_2[i][j][1];
                    Float_t x2 = ParGaus_2[i][j][0] + nSigmaPhi*ParGaus_2[i][j][1];
                    Float_t y = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinContent(h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_2[i][j][0]));
                    PlotLine(x1,x1,0,y,4,2,2);
                    PlotLine(x2,x2,0,y,4,2,2);
                    PlotLine(0.98,1.05,0,0,1,2,2);
                }
                c2_phi_Psi_bw[i][j]->cd(8);
                c2_phi_Psi_bw[i][j]->cd(8)->SetLeftMargin(0.2);
                c2_phi_Psi_bw[i][j]->cd(8)->SetBottomMargin(0.2);
                c2_phi_Psi_bw[i][j]->cd(8)->SetTicks(1,1);
                c2_phi_Psi_bw[i][j]->cd(8)->SetGrid(0,0);
                h_mMass_Phi2_Sig_total[i][j]->Draw("pE");
                f_bw_2[i][j]->SetLineColor(2);
                f_bw_2[i][j]->Draw("l same");
                f_gauss_2[i][j]->SetLineColor(6);
                f_gauss_2[i][j]->Draw("l same");
                PlotLine(0.98,1.05,0,0,1,2,2);
                /*
                CanName = Form("c3_phi_Psi_pt_%d_Centrality_%d_bw",i,j);
                c3_phi_Psi_bw[i][j] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
                c3_phi_Psi_bw[i][j]->Divide(3,3);
                for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
                {
                    c3_phi_Psi_bw[i][j]->cd(m+1);
                    c3_phi_Psi_bw[i][j]->cd(m+1)->SetLeftMargin(0.2);
                    c3_phi_Psi_bw[i][j]->cd(m+1)->SetBottomMargin(0.2);
                    c3_phi_Psi_bw[i][j]->cd(m+1)->SetTicks(1,1);
                    c3_phi_Psi_bw[i][j]->cd(m+1)->SetGrid(0,0);
                    h_mMass_Phi3_SM_New_Bin[i][j][m]->Draw("pE");
                    f_bw_phi_3[i][j][m]->SetLineColor(2);
                    f_bw_phi_3[i][j][m]->Draw("l Same");
                    Float_t x1 = ParGaus_3[i][j][0] - nSigmaPhi*ParGaus_3[i][j][1];
                    Float_t x2 = ParGaus_3[i][j][0] + nSigmaPhi*ParGaus_3[i][j][1];

                    Float_t y = h_mMass_Phi3_SM_New_Bin[i][j][m]->GetBinContent(h_mMass_Phi3_SM_New_Bin[i][j][m]->FindBin(ParGaus_3[i][j][0]));
                    PlotLine(x1,x1,0,y,4,2,2);
                    PlotLine(x2,x2,0,y,4,2,2);
                    PlotLine(0.98,1.05,0,0,1,2,2);
                }
                c3_phi_Psi_bw[i][j]->cd(8);
                c3_phi_Psi_bw[i][j]->cd(8)->SetLeftMargin(0.2);
                c3_phi_Psi_bw[i][j]->cd(8)->SetBottomMargin(0.2);
                c3_phi_Psi_bw[i][j]->cd(8)->SetTicks(1,1);
                c3_phi_Psi_bw[i][j]->cd(8)->SetGrid(0,0);
                h_mMass_Phi3_Sig_total[i][j]->Draw("pE");
                f_bw_3[i][j]->SetLineColor(2);
                //f_bw_3[i][j]->Draw("l same");
                PlotLine(0.98,1.05,0,0,1,2,2);
                */
            }
        }
    }
    cout << "test 13 ! " <<  endl;

    // Draw phi-Psi distribution with Breit Wignar Fit
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                c2_phi_Psi_bw[i][j]->cd(9);
                c2_phi_Psi_bw[i][j]->cd(9)->SetLeftMargin(0.2);
                c2_phi_Psi_bw[i][j]->cd(9)->SetBottomMargin(0.2);
                c2_phi_Psi_bw[i][j]->cd(9)->SetTicks(1,1);
                c2_phi_Psi_bw[i][j]->cd(9)->SetGrid(0,0);
                TString Title = Form("%1.1f < p_{T} < %1.1f GeV/c",pt_low_phi[i],pt_up_phi[i]);
                h_Counts2_bw[i][j]->SetTitle(Title.Data());
                h_Counts2_bw[i][j]->SetNdivisions(505,"X");
                h_Counts2_bw[i][j]->SetNdivisions(505,"Y");
                h_Counts2_bw[i][j]->GetXaxis()->SetLabelSize(0.05);
                h_Counts2_bw[i][j]->GetYaxis()->SetLabelSize(0.05);
                h_Counts2_bw[i][j]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
                h_Counts2_bw[i][j]->GetYaxis()->SetTitle("Counts");
                h_Counts2_bw[i][j]->GetXaxis()->SetTitleSize(0.06);
                h_Counts2_bw[i][j]->GetYaxis()->SetTitleSize(0.06);
                h_Counts2_bw[i][j]->GetYaxis()->SetTitleOffset(1.7);
                h_Counts2_bw[i][j]->GetXaxis()->CenterTitle();
                h_Counts2_bw[i][j]->GetYaxis()->CenterTitle();
                Int_t bin_low = h_Counts2_bw[i][j]->FindBin(PI_max[0]/14.0+6*PI_max[0]/7.0);
                Int_t bin_up  = h_Counts2_bw[i][j]->FindBin(PI_max[0]/14.0+0*PI_max[0]/7.0);
                h_Counts2_bw[i][j]->GetYaxis()->SetRangeUser(h_Counts2_bw[i][j]->GetBinContent(bin_low)*0.7,h_Counts2_bw[i][j]->GetBinContent(bin_up)*1.3);
                h_Counts2_bw[i][j]->Draw("pE");
                TLegend *leg2 = new TLegend(0.5,0.7,0.88,0.8);
                leg2->SetFillColor(0);
                leg2->SetFillStyle(0);
                leg2->SetBorderSize(0);
                leg2->AddEntry(h_Counts2_bw[i][j],"Breit Wigner","p");
                leg2->Draw("same");
                /*
                c3_phi_Psi_bw[i][j]->cd(9);
                c3_phi_Psi_bw[i][j]->cd(9)->SetLeftMargin(0.2);
                c3_phi_Psi_bw[i][j]->cd(9)->SetBottomMargin(0.2);
                c3_phi_Psi_bw[i][j]->cd(9)->SetTicks(1,1);
                c3_phi_Psi_bw[i][j]->cd(9)->SetGrid(0,0);
                h_Counts3_bw[i][j]->SetTitle(Title.Data());
                h_Counts3_bw[i][j]->SetNdivisions(505,"X");
                h_Counts3_bw[i][j]->SetNdivisions(505,"Y");
                h_Counts3_bw[i][j]->GetXaxis()->SetLabelSize(0.05);
                h_Counts3_bw[i][j]->GetYaxis()->SetLabelSize(0.05);
                h_Counts3_bw[i][j]->GetXaxis()->SetTitle("#phi-#Psi_{3}");
                h_Counts3_bw[i][j]->GetYaxis()->SetTitle("Counts");
                h_Counts3_bw[i][j]->GetXaxis()->SetTitleSize(0.06);
                h_Counts3_bw[i][j]->GetYaxis()->SetTitleSize(0.06);
                h_Counts3_bw[i][j]->GetYaxis()->SetTitleOffset(1.7);
                h_Counts3_bw[i][j]->GetXaxis()->CenterTitle();
                h_Counts3_bw[i][j]->GetYaxis()->CenterTitle();
                bin_low = h_Counts3_bw[i][j]->FindBin(PI_max[1]/14.0+6*PI_max[1]/7.0);
                bin_up  = h_Counts3_bw[i][j]->FindBin(PI_max[1]/14.0+0*PI_max[1]/7.0);
                h_Counts3_bw[i][j]->GetYaxis()->SetRangeUser(h_Counts3_bw[i][j]->GetBinContent(bin_low)*0.8,h_Counts3_bw[i][j]->GetBinContent(bin_up)*1.2);
                h_Counts3_bw[i][j]->Draw("pE");
                TLegend *leg3 = new TLegend(0.5,0.7,0.88,0.8);
                leg3->SetFillColor(0);
                leg3->SetBorderSize(0);
                leg3->AddEntry(h_Counts3_bw[i][j],"Breit Wigner","p");
                leg3->Draw("same");
                */
            }
        }
    }
    cout << "test 14 ! " <<  endl;

    // do cos fit to extract raw v2 and v3
    TF1 *f_phi2[pt_total_New_phi][Centrality_total];
    TF1 *f_phi3[pt_total_New_phi][Centrality_total];

    TF1 *f_phi2_bw[pt_total_New_phi][Centrality_total];
    TF1 *f_phi3_bw[pt_total_New_phi][Centrality_total];
    /*
    */
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                TString Flow_phi;

                Flow_phi = Form("flow_pt_%d_Centrality_%d_2nd_phi",i,j);
                f_phi2[i][j] = new TF1(Flow_phi.Data(),flow_2,0.0,PI_max[0],2);
                f_phi2[i][j]->SetParameter(0,2.0);
                f_phi2[i][j]->SetParameter(1,1.0);
                h_Counts2[i][j]->Fit(f_phi2[i][j],"NQM");

                // Flow_phi = Form("flow_pt_%d_Centrality_%d_3rd_phi",i,j);
                // f_phi3[i][j] = new TF1(Flow_phi.Data(),flow_3,0.0,PI_max[1],2);
                // f_phi3[i][j]->SetParameter(0,2.0);
                // f_phi3[i][j]->SetParameter(1,1.0);
                // h_Counts3[i][j]->Fit(f_phi3[i][j],"NQM");

                Flow_phi = Form("flow_pt_%d_Centrality_%d_2nd_phi_bw",i,j);
                f_phi2_bw[i][j] = new TF1(Flow_phi.Data(),flow_2,0.0,PI_max[0],2);
                f_phi2_bw[i][j]->SetParameter(0,2.0);
                f_phi2_bw[i][j]->SetParameter(1,1.0);
                h_Counts2_bw[i][j]->Fit(f_phi2_bw[i][j],"NQM");

                // Flow_phi = Form("flow_pt_%d_Centrality_%d_3rd_phi_bw",i,j);
                // f_phi3_bw[i][j] = new TF1(Flow_phi.Data(),flow_3,0.0,PI_max[1],2);
                // f_phi3_bw[i][j]->SetParameter(0,2.0);
                // f_phi3_bw[i][j]->SetParameter(1,1.0);
                // h_Counts3_bw[i][j]->Fit(f_phi3_bw[i][j],"NQM");
            }
        }
    }
    cout << "test 15 ! " <<  endl;

    // Draw the fit line
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                c2_phi_Psi[i][j]->cd(9);
                f_phi2[i][j]->SetLineColor(2);
                f_phi2[i][j]->Draw("l Same");

                // c3_phi_Psi[i][j]->cd(9);
                // f_phi3[i][j]->SetLineColor(2);
                // f_phi3[i][j]->Draw("l Same");

                c2_phi_Psi_bw[i][j]->cd(9);
                f_phi2_bw[i][j]->SetLineColor(2);
                f_phi2_bw[i][j]->Draw("l Same");

                // c3_phi_Psi_bw[i][j]->cd(9);
                // f_phi3_bw[i][j]->SetLineColor(2);
                // f_phi3_bw[i][j]->Draw("l Same");
            }
        }
    }

    // Draw phi-Psi distribution from Counting and Breit Wignar Fit in the same plot with fit line
    TCanvas *c_phi_Psi2[Centrality_total];
    TCanvas *c_phi_Psi3[Centrality_total];
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            TString CanName;

            CanName = Form("c_phi_Psi2_Centrality_%d",j);
            c_phi_Psi2[j] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,600);
            c_phi_Psi2[j]->Divide(3,2);
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                c_phi_Psi2[j]->cd(i+1);
                c_phi_Psi2[j]->cd(i+1)->SetLeftMargin(0.2);
                c_phi_Psi2[j]->cd(i+1)->SetBottomMargin(0.2);
                c_phi_Psi2[j]->cd(i+1)->SetTicks(1,1);
                c_phi_Psi2[j]->cd(i+1)->SetGrid(0,0);

                h_Counts2[i][j]->SetMarkerStyle(20);
                h_Counts2[i][j]->SetMarkerSize(1.0);
                h_Counts2[i][j]->SetMarkerColor(1);
                h_Counts2[i][j]->Draw("pE");
                TLegend *leg2 = new TLegend(0.5,0.6,0.88,0.8);
                leg2->SetFillColor(0);
                leg2->SetFillStyle(0);
                leg2->SetBorderSize(0);
                leg2->AddEntry(h_Counts2[i][j],"Counting","p");
                f_phi2[i][j]->SetLineColor(1);
                f_phi2[i][j]->SetLineStyle(1);
                f_phi2[i][j]->Draw("l Same");

                h_Counts2_bw[i][j]->SetMarkerStyle(24);
                h_Counts2_bw[i][j]->SetMarkerSize(1.0);
                h_Counts2_bw[i][j]->SetMarkerColor(2);
                h_Counts2_bw[i][j]->Draw("pE same");
                leg2->AddEntry(h_Counts2_bw[i][j],"Breit Wigner","p");
                f_phi2_bw[i][j]->SetLineColor(2);
                f_phi2_bw[i][j]->SetLineStyle(2);
                //f_phi2_bw[i][j]->Draw("l Same");
                leg2->Draw("same");
            }

            /*
            CanName = Form("c_phi_Psi3_Centrality_%d",j);
            c_phi_Psi3[j] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,600);
            c_phi_Psi3[j]->Divide(3,2);
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                c_phi_Psi3[j]->cd(i+1);
                c_phi_Psi3[j]->cd(i+1)->SetLeftMargin(0.2);
                c_phi_Psi3[j]->cd(i+1)->SetBottomMargin(0.2);
                c_phi_Psi3[j]->cd(i+1)->SetTicks(1,1);
                c_phi_Psi3[j]->cd(i+1)->SetGrid(0,0);

                h_Counts3[i][j]->SetMarkerStyle(20);
                h_Counts3[i][j]->SetMarkerSize(1.0);
                h_Counts3[i][j]->SetMarkerColor(1);
                h_Counts3[i][j]->Draw("pE");
                TLegend *leg3 = new TLegend(0.5,0.6,0.88,0.8);
                leg3->SetFillColor(0);
                leg3->SetBorderSize(0);
                leg3->AddEntry(h_Counts3[i][j],"Counting","p");
                f_phi3[i][j]->SetLineColor(1);
                f_phi3[i][j]->SetLineStyle(1);
                f_phi3[i][j]->Draw("l Same");

                h_Counts3_bw[i][j]->SetMarkerStyle(24);
                h_Counts3_bw[i][j]->SetMarkerSize(1.0);
                h_Counts3_bw[i][j]->SetMarkerColor(2);
                h_Counts3_bw[i][j]->Draw("pE same");
                leg3->AddEntry(h_Counts3_bw[i][j],"Breit Wigner","p");
                f_phi3_bw[i][j]->SetLineColor(2);
                f_phi3_bw[i][j]->SetLineStyle(2);
                f_phi3_bw[i][j]->Draw("l Same");
                leg3->Draw("same");
            }
            */
        }
    }
    cout << "test 16 ! " <<  endl;

    // fill raw v2 and v3 into histogram for counting and Breit Wignar
    TH1F *h_flow_2[Centrality_total];
    // TH1F *h_flow_3[Centrality_total];

    TH1F *h_flow_2_bw[Centrality_total];
    // TH1F *h_flow_3_bw[Centrality_total];
    //txtFile << Form("##columns: pT ,  v2 , stat.error, syst.error") << endl;
    //txtFile << Form("#---------- phi-meson v2 ------------") << endl;
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
    cout << "test 16.1 ! " <<  endl;
    	//txtFile << Form("# phi-meson v2 (centrality bin %d) in 19.6 GeV ===> ", j) << endl;
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            TString HistName;

    cout << "test 16.2 ! " <<  endl;
            HistName = Form("flow_2_Centrality_%d_phi",j);
            h_flow_2[j] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.6);
    	    //txtFile << Form("# Counting method ===> ", j) << endl;
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                Float_t bin_center = h_flow_2[j]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
                Float_t bin_content = f_phi2[i][j]->GetParameter(1);
                Float_t bin_error = f_phi2[i][j]->GetParError(1);
                h_flow_2[j]->SetBinContent(bin_center,bin_content);
                h_flow_2[j]->SetBinError(bin_center,bin_error);
		cout << "Counting centrality: "<< j << "pt: "<<i<< " v2: "<< bin_content << " error: "<< bin_error << endl;
    		//txtFile << Form("%1.2f \t%f \t%f", (pt_low_phi[i]+pt_up_phi[i])/2.0, bin_content, bin_error) << endl;
    cout << "test 16.3 ! " <<  endl;
            }
            /*
            HistName = Form("flow_3_Centrality_%d_phi",j);
            h_flow_3[j] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.6);
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                Int_t bin_center = h_flow_3[j]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
                Float_t bin_content = f_phi3[i][j]->GetParameter(1);
                Float_t bin_error = f_phi3[i][j]->GetParError(1);
                h_flow_3[j]->SetBinContent(bin_center,bin_content);
                h_flow_3[j]->SetBinError(bin_center,bin_error);
            }
            */
            HistName = Form("flow_2_Centrality_%d_phi_bw",j);
            h_flow_2_bw[j] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.6);
    	    //txtFile << Form("# BW fitting method ===> ", j) << endl;
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                Float_t bin_center = h_flow_2_bw[j]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
                Float_t bin_content = f_phi2_bw[i][j]->GetParameter(1);
                Float_t bin_error = f_phi2_bw[i][j]->GetParError(1);
                h_flow_2_bw[j]->SetBinContent(bin_center,bin_content);
                h_flow_2_bw[j]->SetBinError(bin_center,bin_error);
		//cout << "BW centrality: "<< j << "pt: "<<i<< " v2: "<< bin_content << " error: "<< bin_error << endl;
    		//txtFile << Form("%1.2f \t%f \t%f", (pt_low_phi[i]+pt_up_phi[i])/2.0, bin_content, bin_error) << endl;
    cout << "test 16.4 ! " <<  endl;
            }
            /*
            HistName = Form("flow_3_Centrality_%d_phi_bw",j);
            h_flow_3_bw[j] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.6);
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                Int_t bin_center = h_flow_3_bw[j]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
                Float_t bin_content = f_phi3_bw[i][j]->GetParameter(1);
                Float_t bin_error = f_phi3_bw[i][j]->GetParError(1);
                h_flow_3_bw[j]->SetBinContent(bin_center,bin_content);
                h_flow_3_bw[j]->SetBinError(bin_center,bin_error);
            }
            */
        }
    }
    cout << "test 17 ! " <<  endl;

    // Draw raw v2 and v3 into phi-Psi distribution
    TCanvas *c_raw_v2[Centrality_total];
    TCanvas *c_raw_v3[Centrality_total];
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            TString CanName;
            CanName = Form("c_raw_v2_Centrality_%d",j);
            c_raw_v2[j] = new TCanvas(CanName.Data(),CanName.Data(),10,10,800,800);
            c_raw_v2[j]->cd();
            c_raw_v2[j]->cd()->SetLeftMargin(0.2);
            c_raw_v2[j]->cd()->SetBottomMargin(0.2);
            c_raw_v2[j]->cd()->SetTicks(1,1);
            c_raw_v2[j]->cd()->SetGrid(0,0);

            h_flow_2[j]->SetMarkerStyle(20);
            h_flow_2[j]->SetMarkerSize(1.0);
            h_flow_2[j]->SetMarkerColor(1);
            h_flow_2[j]->DrawCopy("pE");

            h_flow_2_bw[j]->SetMarkerStyle(24);
            h_flow_2_bw[j]->SetMarkerSize(1.0);
            h_flow_2_bw[j]->SetMarkerColor(2);
            h_flow_2_bw[j]->DrawCopy("pE same");
            /*
            CanName = Form("c_raw_v3_Centrality_%d",j);
            c_raw_v3[j] = new TCanvas(CanName.Data(),CanName.Data(),10,10,800,800);
            c_raw_v3[j]->cd();
            c_raw_v3[j]->cd()->SetLeftMargin(0.2);
            c_raw_v3[j]->cd()->SetBottomMargin(0.2);
            c_raw_v3[j]->cd()->SetTicks(1,1);
            c_raw_v3[j]->cd()->SetGrid(0,0);

            h_flow_3[j]->SetMarkerStyle(20);
            h_flow_3[j]->SetMarkerSize(1.0);
            h_flow_3[j]->SetMarkerColor(1);
            h_flow_3[j]->DrawCopy("pE");

            h_flow_3_bw[j]->SetMarkerStyle(24);
            h_flow_3_bw[j]->SetMarkerSize(1.0);
            h_flow_3_bw[j]->SetMarkerColor(2);
            h_flow_3_bw[j]->DrawCopy("pE same");
            */
        }
    }
    cout << "test 18 ! " <<  endl;

    //--------------------------------------------------------------------------------------------------
    // calculate the final resolution correction
    // read the file from Same Event and Mixed Event
    TH1F *h_Yields_Phi_SE[9];
    TH1F *h_Yields_Phi_ME[9];
    for(Int_t j = 0; j < 9; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            TString HistName;

            HistName = Form("Yields_Centrality_%d_EtaGap_0_Phi_SE",j);
            h_Yields_Phi_SE[j] = (TH1F*)File_Ana_SE->Get(HistName.Data())->Clone();

            HistName = Form("Yields_Centrality_%d_EtaGap_0_Phi_ME",j);
            h_Yields_Phi_ME[j] = (TH1F*)File_Ana_ME->Get(HistName.Data())->Clone();
        }
    }
    cout << "test 19 ! " <<  endl;
    
    TH1F *h_Yields_Phi_SM[9];
    for(Int_t j = 0; j < 9; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            Int_t bin_SE_start = h_Yields_Phi_SE[j]->FindBin(1.04);
            Int_t bin_SE_stop  = h_Yields_Phi_SE[j]->FindBin(1.05);
            Float_t Inte_SE = h_Yields_Phi_SE[j]->Integral(bin_SE_start,bin_SE_stop);
            h_Yields_Phi_SM[j] = (TH1F*)h_Yields_Phi_SE[j]->Clone();

            Int_t bin_ME_start = h_Yields_Phi_ME[j]->FindBin(1.04);
            Int_t bin_ME_stop  = h_Yields_Phi_ME[j]->FindBin(1.05);
            Float_t Inte_ME = h_Yields_Phi_ME[j]->Integral(bin_ME_start,bin_ME_stop);
            h_Yields_Phi_ME[j]->Scale(Inte_SE/Inte_ME);
            h_Yields_Phi_SM[j]->Add(h_Yields_Phi_ME[j],-1.0);
        }
    }
    cout << "test 20 ! " <<  endl;

    TF1  *f_PolyBW_Phi_Yields[9];

    Float_t ParFit_Yields[9][5];

    for(Int_t j = 0; j < 9; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            TString HistName;

            HistName = Form("f_Centrality_%d_Yields_Phi_SE",j);
            f_PolyBW_Phi_Yields[j] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
            for(Int_t n_par = 0; n_par < 5; n_par++)
            {
                f_PolyBW_Phi_Yields[j]->ReleaseParameter(n_par);
            }
            f_PolyBW_Phi_Yields[j]->SetParameter(0,1.019);
            f_PolyBW_Phi_Yields[j]->SetParLimits(0,1.014,1.024);
            f_PolyBW_Phi_Yields[j]->SetParameter(1,0.0055);
            f_PolyBW_Phi_Yields[j]->SetParameter(2,10000);
            f_PolyBW_Phi_Yields[j]->SetParameter(3,-6000);
            f_PolyBW_Phi_Yields[j]->SetParameter(4,0.5);
            f_PolyBW_Phi_Yields[j]->SetRange(BW_Start,BW_Stop);

            h_Yields_Phi_SM[j]->Fit(f_PolyBW_Phi_Yields[j],"RM");
            for(Int_t n_par = 0; n_par < 5; n_par++)
            {
                ParFit_Yields[j][n_par] = f_PolyBW_Phi_Yields[j]->GetParameter(n_par);
            }
        }
    }
    cout << "test 21 ! " <<  endl;

    TF1 *f_Poly_Phi_Yields[9];

    for(Int_t j = 0; j < 9; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            TString FuncName;

            FuncName = Form("f_poly_Centrality_%d_Yields",j);
            f_Poly_Phi_Yields[j] = new TF1(FuncName.Data(),Poly,BW_Start,BW_Stop,2);
            f_Poly_Phi_Yields[j]->FixParameter(0,ParFit_Yields[j][3]);
            f_Poly_Phi_Yields[j]->FixParameter(1,ParFit_Yields[j][4]);
            h_Yields_Phi_SM[j]->Add(f_Poly_Phi_Yields[j],-1.0);
        }
    }
    cout << "test 22 ! " <<  endl;

    // calculate total counts and for each centrality bin and eta bin
    Float_t Counts_Yields[9];
    for(Int_t j = 0; j < 9; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            Counts_Yields[j] = 0.0;
            Int_t bin_start = h_Yields_Phi_SM[j]->FindBin(ParFit_Yields[j][0]-3.0*ParFit_Yields[j][1]);
            Int_t bin_stop  = h_Yields_Phi_SM[j]->FindBin(ParFit_Yields[j][0]+3.0*ParFit_Yields[j][1]);
            for(Int_t bin = bin_start; bin <= bin_stop; bin++)
            {
                  
	        if( h_Yields_Phi_SM[j]->GetBinContent(bin) >=0) Counts_Yields[j] += h_Yields_Phi_SM[j]->GetBinContent(bin);
            }
        }
    }
    cout << "test 23 ! " <<  endl;

    // get resolution
    TString input_res = Form("data/19gev_shifted.root");
    TFile *input = TFile::Open(input_res.Data());
    TProfile *p_res2;
    TProfile *p_res3;
    Double_t mean_res_2_phi[Centrality_total];
    Double_t mean_res_3_phi[Centrality_total];
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++)
    {
        TString ProName;
        ProName = Form("p_r2_tpc_AB_sub0_1");
        p_res2 = (TProfile*)input->Get(ProName.Data());
        // ProName = Form("Res3_EP",l);
        // p_res3[l] = (TProfile*)input->Get(ProName.Data());
    }
    cout << "test 24 ! " <<  endl;

    for(Int_t j = Centrality_start; j < Centrality_stop; j++)
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            mean_res_2_phi[j] = 0.0;
            // mean_res_3_phi[j] = 0.0;

            //Double_t res_2[9]={0.135779, 0.219408, 0.248626, 0.208038, 0.132007,1,1,1,1};
            //Double_t res_2[9]={1,1,1,1, 0.132007, 0.208038, 0.248626, 0.219408, 0.135779};
            Double_t res_2[9]={0.0};// dchen
            Float_t phi2_total = 0.0;
            for(Int_t cent = cent_low[j]; cent <= cent_up[j]; cent++)
            {
                if( (p_res2->GetBinContent(cent+1)) > 0 )
		{
                 res_2[cent] = TMath::Sqrt(p_res2->GetBinContent(cent+1));
		} else
		{
                 res_2[cent] = -999.0;
		}
                //res_2[cent] = 0.3; // dchen test
                phi2_total += Counts_Yields[cent];
                cout <<  "centality_bin = " << cent << ", Counts_Yields = " << Counts_Yields[cent] << endl;
            }
	    cout << "phi2_total = " << phi2_total << endl;
            for(Int_t cent = cent_low[j]; cent <= cent_up[j]; cent++)
            {
                mean_res_2_phi[j] += Counts_Yields[cent]/(res_2[cent]*phi2_total);
		cout << "res2= " << res_2[cent] << endl;
            }
	    //mean_res_2_phi[j] = 5.0; // test dchen
            cout <<  "centality_bin = " << j << ", mean_res_2_phi = " << mean_res_2_phi[j] << endl;
            /*
            Float_t phi3_total = 0.0;
            Double_t res_3[9];
            for(Int_t i = 0; i < 9; i++)
            {
                res_3[i] = -999.9;
            }
            for(Int_t cent = cent_low[j]; cent <= cent_up[j]; cent++)
            {
                if(p_res3[l]->GetBinContent(cent+1) > 0)
                {
                    res_3[cent] = TMath::Sqrt(p_res3[l]->GetBinContent(cent+1));
                    phi3_total += Counts_Yields[cent][l];
                }
            }
            for(Int_t cent = cent_low[j]; cent <= cent_up[j]; cent++)
            {
                mean_res_3_phi[j] += Counts_Yields[cent][l]/(res_3[cent]*phi3_total);
            }
            cout <<  "centality_bin = " << j << ", eta_bin = " << l << ", mean_res_3_phi = " << mean_res_3_phi[j] << endl;
            */
        }
    }
    cout << "test 25 ! " <<  endl;

    // scale raw v2 and v3 with resoltuion correction
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            h_flow_2[j]->Scale(mean_res_2_phi[j]);
            // h_flow_3[j]->Scale(mean_res_3_phi[j]);
            h_flow_2_bw[j]->Scale(mean_res_2_phi[j]);
            // h_flow_3_bw[j]->Scale(mean_res_3_phi[j]);
        }
    }
    cout << "test 26 ! " <<  endl;

    // save v2 and v3 into TGraphAsymmErrors and into file
    TGraphAsymmErrors *g_flow_2[Centrality_total];
    // TGraphAsymmErrors *g_flow_3[Centrality_total];
    TGraphAsymmErrors *g_flow_2_bw[Centrality_total];
    // TGraphAsymmErrors *g_flow_3_bw[Centrality_total];
    txtFile << Form("##columns: pT ,  v2 , stat.error, syst.error") << endl;
    txtFile << Form("#---------- phi-meson v2 ------------") << endl;
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
    	txtFile << Form("# phi-meson v2 (centrality bin %d) in 19.6 GeV ===> ", j) << endl;
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
    	    txtFile << Form("# Counting method ===> ", j) << endl;
            TString g_Name_2 = Form("g_Centrality_%d_2nd",j);
            g_flow_2[j] = new TGraphAsymmErrors();
            g_flow_2[j]->SetName(g_Name_2.Data());
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                Float_t bin_center = h_flow_2[j]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
                Float_t bin_content = h_flow_2[j]->GetBinContent(bin_center);
                Float_t bin_error   = h_flow_2[j]->GetBinError(bin_center);
		cout << "Counting centrality: "<< j << ", pt: "<<i<< " v2: "<< bin_content << " error: "<< bin_error << endl;
                g_flow_2[j]->SetPoint(i,(pt_low_phi[i]+pt_up_phi[i])/2.0,bin_content);
                g_flow_2[j]->SetPointError(i,0.0,0.0,bin_error,bin_error);
    		txtFile << Form("%1.2f \t%f \t%f", (pt_low_phi[i]+pt_up_phi[i])/2.0, bin_content, bin_error) << endl;
            }
            /*
            TString g_Name_3 = Form("g_Centrality_%d_3rd",j);
            g_flow_3[j] = new TGraphAsymmErrors();
            g_flow_3[j]->SetName(g_Name_3.Data());
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                Int_t bin_center = h_flow_3[j]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
                Float_t bin_content = h_flow_3[j]->GetBinContent(bin_center);
                Float_t bin_error   = h_flow_3[j]->GetBinError(bin_center);
                g_flow_3[j]->SetPoint(i,(pt_low_phi[i]+pt_up_phi[i])/2.0,bin_content);
                g_flow_3[j]->SetPointError(i,0.0,0.0,bin_error,bin_error);
            }
            */
            TString g_Name_2_bw = Form("g_Centrality_%d_2nd_bw",j);
            g_flow_2_bw[j] = new TGraphAsymmErrors();
            g_flow_2_bw[j]->SetName(g_Name_2_bw.Data());
    	    txtFile << Form("# BW fitting method ===> ", j) << endl;
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                Float_t bin_center = h_flow_2_bw[j]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
                Float_t bin_content = h_flow_2_bw[j]->GetBinContent(bin_center);
                Float_t bin_error   = h_flow_2_bw[j]->GetBinError(bin_center);
		cout << "BW centrality: "<< j << ", pt: "<<i<< " v2: "<< bin_content << " error: "<< bin_error << endl;
                g_flow_2_bw[j]->SetPoint(i,(pt_low_phi[i]+pt_up_phi[i])/2.0+0.05,bin_content);
                g_flow_2_bw[j]->SetPointError(i,0.0,0.0,bin_error,bin_error);
    		txtFile << Form("%1.2f \t%f \t%f", (pt_low_phi[i]+pt_up_phi[i])/2.0, bin_content, bin_error) << endl;
            }
	    cout << endl;
            /*
            TString g_Name_3_bw = Form("g_Centrality_%d_3rd_bw",j);
            g_flow_3_bw[j] = new TGraphAsymmErrors();
            g_flow_3_bw[j]->SetName(g_Name_3_bw.Data());
            for(Int_t i = 0; i < pt_total_New_phi; i++)
            {
                Int_t bin_center = h_flow_3_bw[j]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
                Float_t bin_content = h_flow_3_bw[j]->GetBinContent(bin_center);
                Float_t bin_error   = h_flow_3_bw[j]->GetBinError(bin_center);
                g_flow_3_bw[j]->SetPoint(i,(pt_low_phi[i]+pt_up_phi[i])/2.0+0.05,bin_content);
                g_flow_3_bw[j]->SetPointError(i,0.0,0.0,bin_error,bin_error);
            }
            */
        }
    }
    cout << "test 27 ! " <<  endl;

    TFile *File_OutPut = new TFile("flow/Flow_Phi.root","RECREATE");
    File_OutPut->cd();
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            h_flow_2[j]->Write();
            g_flow_2[j]->Write();
            // g_flow_3[j]->Write();
            g_flow_2_bw[j]->Write();
            h_flow_2_bw[j]->Write();
            // g_flow_3_bw[j]->Write();
        }
    }
    cout << "test 28 ! " <<  endl;
    File_OutPut->Close();
    //  File_Ana->Close();
    //  File_Ana->Close();
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                TString CanName;

                CanName = Form("./figures_phi/c2_phi_Psi_pt_%d_Centrality_%d.pdf",i,j);
                //c2_phi_Psi[i][j]->SaveAs(CanName.Data());
                // CanName = Form("./figures_phi/c3_phi_Psi_pt_%d_Centrality_%d.pdf",i,j);
                // c3_phi_Psi[i][j]->SaveAs(CanName.Data());
                CanName = Form("./figures_phi/c2_phi_Psi_pt_%d_Centrality_%d_bw.pdf",i,j);
                //c2_phi_Psi_bw[i][j]->SaveAs(CanName.Data());
                // CanName = Form("./figures_phi/c3_phi_Psi_pt_%d_Centrality_%d_bw.pdf",i,j);
                // c3_phi_Psi_bw[i][j]->SaveAs(CanName.Data());
            }
        }
    }
    cout << "test 29 ! " <<  endl;
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
        for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
        {
            TString CanName;

            CanName = Form("./figures_phi/c_phi_Psi2_Centrality_%d.pdf",j);
            //c_phi_Psi2[j]->SaveAs(CanName.Data());
            // CanName = Form("./figures_phi/c_phi_Psi3_Centrality_%d.pdf",j);
            // c_phi_Psi3[j]->SaveAs(CanName.Data());
        }
    }
    cout << "test 30 ! " <<  endl;

    // draw poster fit figure
    TCanvas *c_play = new TCanvas("c_play","c_play",1400,10,900,900);
    c_play->SetLeftMargin(0.15);
    c_play->SetBottomMargin(0.15);
    c_play->SetTicks(1,1);
    c_play->SetGrid(0,0);
    c_play->SetFillColor(0);
    c_play->SetFillStyle(4000);
    c_play->SetFrameFillColor(18);
    /*
    //  h_Counts3[10][0][0]->SetTitle("");
    h_Counts3[2][0][0]->GetXaxis()->SetTitle("#phi-#Psi_{3}");
    h_Counts3[2][0][0]->GetYaxis()->SetTitle("dN/d(#phi-#Psi_{3})");
    h_Counts3[2][0][0]->GetXaxis()->SetTitleSize(0.06);
    h_Counts3[2][0][0]->GetYaxis()->SetTitleSize(0.06);
    h_Counts3[2][0][0]->GetXaxis()->SetTitleOffset(0.9);
    h_Counts3[2][0][0]->GetXaxis()->CenterTitle();
    h_Counts3[2][0][0]->GetYaxis()->CenterTitle();
    h_Counts3[2][0][0]->GetXaxis()->SetTitleColor(0);
    h_Counts3[2][0][0]->GetYaxis()->SetTitleColor(0);
    h_Counts3[2][0][0]->SetNdivisions(505,"X");
    h_Counts3[2][0][0]->SetNdivisions(505,"Y");
    h_Counts3[2][0][0]->GetXaxis()->SetLabelSize(0.04);
    h_Counts3[2][0][0]->GetYaxis()->SetLabelSize(0.04);
    h_Counts3[2][0][0]->GetXaxis()->SetLabelColor(0);
    h_Counts3[2][0][0]->GetYaxis()->SetLabelColor(0);
    h_Counts3[2][0][0]->Draw("pE");
    */
    TFile *file_counts = new TFile("flow/Counts.root","RECREATE");
    file_counts->cd();
    for(Int_t i = 0; i < pt_total_New_phi; i++)
    {
        for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
            {
                h_Counts2[i][j]->Write();
                // h_Counts3[i][j]->Write();
            }
        }
    }
    cout << "test end ! " <<  endl;
    file_counts->Close();
    txtFile.close();
    txtFile_invm.close();
}
/*
Double_t PolyBreitWigner(Double_t *x_val, Double_t *par)
{
    Double_t x = x_val[0];
    Double_t m0 = par[0];
    Double_t Gamma = par[1];
    Double_t Norm = par[2];

    Double_t denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
    Double_t BW = Norm*Gamma/denom;

    Double_t Poly = par[3] + par[4]*x;

    Double_t y = BW + Poly;

    return y;
}
*/
// Fitting functions for Flow VS Invariant Mass
Double_t BackgroundFitting(Double_t *x, Double_t *p)
{    
    Double_t m0 = dParSig[0];
    Double_t Gamma = dParSig[1];
    Double_t Norm = dParSig[2];

    Double_t denom = 2.0*TMath::Pi()*((x[0]-m0)*(x[0]-m0)+Gamma*Gamma/4.0);
    Double_t BW = Norm*Gamma/denom;

    Double_t Poly = dParSig[3] + dParSig[4]*x[0];

    Double_t Sig = BW + Poly;


  //if( x[0] >= (dParSig[1] - (nSigmaPhi*dParSig[2])) && x[0] <= (dParSig[1] + (nSigmaPhi*dParSig[2])))
   {
     return (((dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))
     /(Sig +(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))))
     * /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3)*/);
   }
   //else
   //{
    // return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3)*/;
   //}
}
Double_t SignalFitting(Double_t *x, Double_t *p)
{    
    Double_t m0 = dParSig[0];
    Double_t Gamma = dParSig[1];
    Double_t Norm = dParSig[2];

    Double_t denom = 2.0*TMath::Pi()*((x[0]-m0)*(x[0]-m0)+Gamma*Gamma/4.0);
    Double_t BW = Norm*Gamma/denom;

    Double_t Poly = dParSig[3] + dParSig[4]*x[0];

    Double_t Sig = BW + Poly;


  //if( x[0] >= (dParSig[1] - (nSigmaPhi*dParSig[2])) && x[0] <= (dParSig[1] + (nSigmaPhi*dParSig[2])))
   {
     return (( Sig
     /(Sig +(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))))
     * /*p[0]*//*(p[0] + p[1] * x[0])*/p[0] /*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3)*/);
   }
   //else
   //{
    // return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3)*/;
   //}
}
Double_t TotalFitting(Double_t *x, Double_t *p)
{   
    Double_t m0 = dParSig[0];
    Double_t Gamma = dParSig[1];
    Double_t Norm = dParSig[2];

    Double_t denom = 2.0*TMath::Pi()*((x[0]-m0)*(x[0]-m0)+Gamma*Gamma/4.0);
    Double_t BW = Norm*Gamma/denom;

    Double_t Poly = dParSig[3] + dParSig[4]*x[0];

    Double_t Sig = BW + Poly;


  //if( x[0] >= (dParSig[1] - (nSigmaPhi*dParSig[2])) && x[0] <= (dParSig[1] + (nSigmaPhi*dParSig[2])))
   //{
     return (
       (p[/*1*//*2*/3/*4*/]* Sig
     /( Sig +(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2)))))
     +
       (((dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))
     /( Sig +(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))))
     * /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3))*/);
   //}
   //else
   //{
     //return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3))*/;
   //}
}
