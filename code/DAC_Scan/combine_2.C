#include "sPhenixStyle.C"

// TGraphAsymmErrors * h1D_to_GRAE(TH1D * h1D_in) {
//     TGraphAsymmErrors * g = new TGraphAsymmErrors();

//     for (int i = 1; i <= h1D_in->GetNbinsX(); i++) {
//         double x = h1D_in->GetBinCenter(i);
//         double y = h1D_in->GetBinContent(i);
//         double y_err = h1D_in->GetBinError(i);
//         double x_err = 0.5 * h1D_in->GetBinWidth(i);

//         if (y <)

//         g->SetPoint(g->GetN(), x, y);
//         g->SetPointError(g->GetN() - 1, x_err, x_err, y_err, y_err);
//     }


//     return g;

// }

double Yaxis_max = 2200;
std::pair<double,double> fit_range_pair = {8, 140};
std::pair<double,double> signal_purity_range = {40, 210};
bool ShowReducedChi2 = true;

const std::vector<std::string> color_code = {
    "#9e0142",
    "#66c2a5",
    "#f46d43",
    "#3288bd",
    "#fee08b",
    "#5e4fa2",
    "#00A1FF",
    "#FF42A1",
    "#000000",
    
    
    "#abdda4",
    "#e6f598",
    "#fdae61",
    "#d53e4f"
};

const std::vector<std::string> syst_HighLow_color_code = {
    "#9e0142",
    "#FF42A1",

    "#005493",
    "#76D6FF",

    "#fee08b",
    "#f46d43",

    "#367C22",
    "#abdda4"
};

const std::vector<int> marker_code = {
    25,
    28,
    27,
    26
    // 21,
    // 23,
    // 27,
    // 28,
    // 29
};

struct RelativeVariationSet{
    TGraph * MaxVariationUp;
    TGraph * MaxVariationDown;

    double all_MaxVariationUp;
    double all_MaxVariationDown;
};

double Func_Langaus(double *x, double *par)
{
	// + Fit parameters:
	//------------------
	// * Signal part:
	// * - par[0]: Width (scale) parameter of Landau density
	// * - par[1]: Most Probable (MP, location) parameter of Landau density
	// * - par[2]: Total area (integral -inf to inf, normalization constant)
	// * - par[3]: Width (sigma) of convoluted Gaussian function

	// par 4 : size 
	// par 5 : width
	// par 6 : scale
	
	//Note1: In the Landau distribution (represented by the CERNLIB approximation),
	//      the maximum is located at x=-0.22278298 with the location parameter=0.
	//      This shift is corrected within this function, so that the actual
	//      maximum is identical to the MP parameter.
	
	//Note2: In a convolution, the variable in the integral run from -inf to +inf.
	//       We can replace the infinity by a number whose magnitude large enough
	//       that beyond its value, the contribution to the convolution is neglectable
	
	// + Numeric constants
	//--------------------
	double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
	double mpshift  = -0.22278298;       // Landau maximum location
	
	// + Control constants
	//--------------------
	double nStep  = 100.0;    // number of convolution steps
	double extSig = 5.0;      // convolution extends to +- [extSig * Gaussian sigmas]
	
	// + Variables
	//------------
	double xx;            // Variable of integration
	double meanLand;      // "Mean" value of Landau distribution used in ROOT
	double fland;         // Landau in the convolution integral
	double fgaus;         // Gaussian in the convolution integral
	double sum = 0.0;     // The sum that replace the integral
	double xlow,xupp;     // Lowest and highest boundary of the integration
	double step;          // step size
	double signal;        // Value of the convolution Land x Gaus at x[0]
	double background;    // Value of the background at x[0]
	
	// + MP shift correction
	//----------------------
	meanLand = par[1] - mpshift * par[0];
	
	// + Range of convolution integral
	//--------------------------------
	xlow = x[0] - extSig * par[3];
	xupp = x[0] + extSig * par[3];
	
	step = (xupp-xlow) / nStep;
	
	// + Convolution at x[0]: integral of Landau and Gaussian by sum
	//--------------------------------------------------------------
	for(double i=1.0; i<=nStep/2; i++)
	{
		xx = xlow + (i-0.5) * step;
		fland = TMath::Landau(xx,meanLand,par[0]) / par[0];
		fgaus = TMath::Gaus(x[0],xx,par[3]);
		sum += fland * fgaus;
		
		xx = xupp - (i-0.5) * step;
		fland = TMath::Landau(xx,meanLand,par[0]) / par[0];
		fgaus = TMath::Gaus(x[0],xx,par[3]);
		sum += fland * fgaus;
	}
	signal = par[2] * step * sum * invsq2pi / par[3];

	// double gaussian_eq  = par[4]*( 1/( par[5]* sqrt(2*TMath::Pi()) ) ) * TMath::Exp( -1*( pow(x[0]-par[6],2)/(2*pow(par[5],2)) ) );
	
	return signal;
}

double Func_Langaus_bkg(double *x, double *par)
{
	// + Fit parameters:
	//------------------
	// * Signal part:
	// * - par[0]: Width (scale) parameter of Landau density
	// * - par[1]: Most Probable (MP, location) parameter of Landau density
	// * - par[2]: Total area (integral -inf to inf, normalization constant)
	// * - par[3]: Width (sigma) of convoluted Gaussian function

	// par 4 : size 
	// par 5 : width
	// par 6 : scale
	
	//Note1: In the Landau distribution (represented by the CERNLIB approximation),
	//      the maximum is located at x=-0.22278298 with the location parameter=0.
	//      This shift is corrected within this function, so that the actual
	//      maximum is identical to the MP parameter.
	
	//Note2: In a convolution, the variable in the integral run from -inf to +inf.
	//       We can replace the infinity by a number whose magnitude large enough
	//       that beyond its value, the contribution to the convolution is neglectable
	
	// + Numeric constants
	//--------------------
	double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
	double mpshift  = -0.22278298;       // Landau maximum location
	
	// + Control constants
	//--------------------
	double nStep  = 100.0;    // number of convolution steps
	double extSig = 5.0;      // convolution extends to +- [extSig * Gaussian sigmas]
	
	// + Variables
	//------------
	double xx;            // Variable of integration
	double meanLand;      // "Mean" value of Landau distribution used in ROOT
	double fland;         // Landau in the convolution integral
	double fgaus;         // Gaussian in the convolution integral
	double sum = 0.0;     // The sum that replace the integral
	double xlow,xupp;     // Lowest and highest boundary of the integration
	double step;          // step size
	double signal;        // Value of the convolution Land x Gaus at x[0]
	double background;    // Value of the background at x[0]
	
	// + MP shift correction
	//----------------------
	meanLand = par[1] - mpshift * par[0];
	
	// + Range of convolution integral
	//--------------------------------
	xlow = x[0] - extSig * par[3];
	xupp = x[0] + extSig * par[3];
	
	step = (xupp-xlow) / nStep;
	
	// + Convolution at x[0]: integral of Landau and Gaussian by sum
	//--------------------------------------------------------------
	for(double i=1.0; i<=nStep/2; i++)
	{
		xx = xlow + (i-0.5) * step;
		fland = TMath::Landau(xx,meanLand,par[0]) / par[0];
		fgaus = TMath::Gaus(x[0],xx,par[3]);
		sum += fland * fgaus;
		
		xx = xupp - (i-0.5) * step;
		fland = TMath::Landau(xx,meanLand,par[0]) / par[0];
		fgaus = TMath::Gaus(x[0],xx,par[3]);
		sum += fland * fgaus;
	}
	signal = par[2] * step * sum * invsq2pi / par[3];

	// double gaussian_eq  = par[4]*( 1/( par[5]* sqrt(2*TMath::Pi()) ) ) * TMath::Exp( -1*( pow(x[0]-par[6],2)/(2*pow(par[5],2)) ) );
	
	double exp_bkg_eq = par[4] * TMath::Exp(par[5] * (x[0] + par[6]));

	double gaus_bkg_eq = par[7]*TMath::Gaus(x[0],0,par[8]);

	return signal + exp_bkg_eq + gaus_bkg_eq;
}

double SingleExp(double *x, double *par)
{
	double exp_bkg_eq = par[0] * TMath::Exp(par[1] * (x[0] + par[2]));	

	return exp_bkg_eq;
}

double SingleGaus(double *x, double *par)
{
	// double gaussian_eq  = par[4]*( 1/( par[5]* sqrt(2*TMath::Pi()) ) ) * TMath::Exp( -1*( pow(x[0]-par[6],2)/(2*pow(par[5],2)) ) );

	double gaus_bkg_eq = par[0]*TMath::Gaus(x[0],par[1],par[2]);

	return gaus_bkg_eq;
}

void GetVariation (
    std::vector<TH1D*> VariationSet_vec, // note : the first one has to be the baseline
    std::vector<std::string> file_title_vec,
    std::string leg_header_in,
    std::string final_output_directory,
    std::string plot_name_in
) {
    SetsPhenixStyle();

    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleAlign(22);
    gStyle->SetTitleX(0.5);
	gStyle->SetTitleY(0.95);
	gStyle->SetTitleFillColor(10);
	gStyle->SetPadRightMargin(0.09);
    gStyle->SetPadTopMargin(0.1);

    TCanvas * c1 = new TCanvas("", "", 950, 800);
    c1 -> cd();

    TPad * pad1 = new TPad("", "", 0, 0.31, 1.0, 0.85);
    pad1 -> SetRightMargin(0.05);
    pad1 -> Draw();
    pad1 -> cd();

    TLegend * leg_variation = new TLegend(0.21,0.83,0.41,0.98);
    leg_variation -> SetBorderSize(0);
    leg_variation -> SetTextSize(0.03);
    leg_variation -> SetHeader(leg_header_in.c_str());

    TLine * line = new TLine();
    line -> SetLineStyle(2);
    line -> SetLineWidth(2);
    line -> SetLineColor(28);

    TLatex * ltx_draw = new TLatex();
    ltx_draw->SetNDC();
    ltx_draw->SetTextSize(0.03);

    TLatex * ltx_draw_all = new TLatex();
    ltx_draw_all->SetNDC();
    ltx_draw_all->SetTextSize(0.13);

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    // ltx->SetTextAlign(31);

    VariationSet_vec[0] -> SetMarkerStyle(20);
    VariationSet_vec[0] -> SetMarkerSize(0.8);
    VariationSet_vec[0] -> SetMarkerColor(1);
    VariationSet_vec[0] -> SetLineColor(1);
    VariationSet_vec[0] -> SetLineWidth(1);
    VariationSet_vec[0] -> GetXaxis() -> SetTitle("DAC value");
    VariationSet_vec[0] -> GetYaxis() -> SetTitle("Entries (A.U.)");
    VariationSet_vec[0] -> GetXaxis() -> SetLabelOffset(999);
    // VariationSet_vec[0] -> GetXaxis() -> SetLabelSize(0);
    VariationSet_vec[0] -> GetXaxis() -> SetTitleOffset(999);

    // VariationSet_vec[0] -> GetYaxis() -> SetRangeUser(
    //     VariationSet_vec[0] -> GetYaxis() -> GetXmin() * 0.8,
    //     VariationSet_vec[0] -> GetYaxis() -> GetXmax() * 1.2
    // );

    VariationSet_vec[0] -> SetMaximum(Yaxis_max);
    VariationSet_vec[0] -> SetMinimum(0);

    // VariationSet_vec[0] -> GetXaxis() -> SetLimits(left_pos_edge, right_pos_edge);

    VariationSet_vec[0] -> Draw("ep");
    leg_variation -> AddEntry(VariationSet_vec[0], file_title_vec[0].c_str(), "ep");

    // note : i starts at 1, because the first one is the baseline
    for (int i = 1; i < VariationSet_vec.size(); i++) {
        VariationSet_vec[i] -> SetMarkerStyle(marker_code[i-1]);
        VariationSet_vec[i] -> SetMarkerSize(1);
        VariationSet_vec[i] -> SetMarkerColor(TColor::GetColor(color_code[i-1].c_str()));
        VariationSet_vec[i] -> SetLineColorAlpha(1,0);
        VariationSet_vec[i] -> Draw("p same");

        leg_variation -> AddEntry(VariationSet_vec[i], file_title_vec[i].c_str(), "p");
    }

    // line -> DrawLine(
    //     VariationSet_vec[0]->GetXaxis()->GetXmin(), 99., 
    //     VariationSet_vec[0]->GetXaxis()->GetXmax(), 99.
    // );

    ltx_draw -> DrawLatex(0.4, 0.8, "*Statiscial error bars of varied data points are removed for better visualization");

    // Division : -----------------------------------------------------------------------------------------------------------------------------------------------------

    c1 -> cd();
    TPad * pad2 = new TPad("", "", 0, 0., 1.0, 0.38);
    pad2 -> SetRightMargin(0.05);
    pad2 -> Draw();
    pad2 -> cd();

    for (int i = 1; i < VariationSet_vec.size(); i++) {
        std::string plot_draw_text = (i == 1) ? "p" : "p same";
        
        // TGraph * temp_ratio_h1D = GetRatioGr(
        //     VariationSet_vec[i],
        //     VariationSet_vec[0] // note : baseline 
        // );

        TH1D * temp_ratio_h1D = (TH1D*) VariationSet_vec[i] -> Clone(Form("temp_ratio_h1D_%d", i));
        temp_ratio_h1D -> Divide(VariationSet_vec[0]);

        temp_ratio_h1D -> SetMaximum(2);
        temp_ratio_h1D -> SetMinimum(0);
        temp_ratio_h1D -> GetXaxis() -> CenterTitle(true);
        temp_ratio_h1D -> GetXaxis() -> SetTitle("DAC value");
        temp_ratio_h1D -> GetYaxis() -> SetTitle("Ratio to baseline");
        // temp_ratio_h1D -> GetXaxis() -> SetLimits(left_pos_edge, right_pos_edge);
        temp_ratio_h1D -> Draw(plot_draw_text.c_str());
    }

    line -> DrawLine(VariationSet_vec[0]->GetXaxis()->GetXmin(), 1, VariationSet_vec[0]->GetXaxis()->GetXmax(), 1);

    c1 -> cd();
    leg_variation -> Draw("same");

   
    c1 -> Print(Form("%s/%s.pdf", final_output_directory.c_str(), plot_name_in.c_str()));
    // c1 -> Clear();

    // Division : -----------------------------------------------------------------------------------------------------------------------------------------------------

    TGraph * MaxVariationUp = new TGraph();
    TGraph * MaxVariationDown = new TGraph();

    TGraphErrors * MaxVariationUp_Ratio = new TGraphErrors();
    TGraphErrors * MaxVariationDown_Ratio = new TGraphErrors();

    MaxVariationUp_Ratio -> SetMarkerStyle(20);
    MaxVariationUp_Ratio -> SetMarkerSize(0.8);
    MaxVariationUp_Ratio -> SetMarkerColor(kBlack);
    MaxVariationUp_Ratio -> SetLineWidth(1);
    MaxVariationUp_Ratio -> SetLineColor(kBlack);
    
    MaxVariationUp_Ratio -> GetXaxis() -> SetTitle("DAC value");
    MaxVariationUp_Ratio -> GetYaxis() -> SetTitle("Maximal Variation / baseline");
    
    MaxVariationDown_Ratio -> SetMarkerStyle(26);
    MaxVariationDown_Ratio -> SetMarkerSize(0.8);
    MaxVariationDown_Ratio -> SetMarkerColor(kRed);
    MaxVariationDown_Ratio -> SetLineWidth(1);
    MaxVariationDown_Ratio -> SetLineColor(kRed);


    for (int i = 0; i < VariationSet_vec.size(); i++) {
        if (i == 0) {
            for (int j = 1; j <= VariationSet_vec[i] -> GetNbinsX(); j++) {
                double x = VariationSet_vec[i] -> GetBinCenter(j);
                double y = VariationSet_vec[i] -> GetBinContent(j);
                MaxVariationUp->SetPoint(MaxVariationUp->GetN(), x, y);
                MaxVariationDown->SetPoint(MaxVariationDown->GetN(), x, y);
                // std::cout<<"Index: "<<i<<", BinCenter: "<<x<<", BinContent: "<<y<<std::endl;
            }
        }

        else {
            for (int j = 1; j <= VariationSet_vec[i] -> GetNbinsX(); j++) {
                double x = VariationSet_vec[i] -> GetBinCenter(j);
                double y = VariationSet_vec[i] -> GetBinContent(j);

                if ( fabs(MaxVariationUp->GetPointX(j-1) - x) > 0.0001 ) {
                    std::cout<<"Error : the x value of the points are not the same"<<std::endl;
                    std::cout<<"x1 : "<<MaxVariationUp->GetPointX(j-1)<<", x2 : "<<x<<std::endl;
                    exit(1);
                }

                if (y > MaxVariationUp->GetPointY(j-1)) {
                    MaxVariationUp->SetPoint(j-1, x, y);
                }

                if (y < MaxVariationDown->GetPointY(j-1)) {
                    MaxVariationDown->SetPoint(j-1, x, y);
                }
            }
        }

    }

    for (int i = 0; i < MaxVariationUp->GetN(); i++) {
        
        double ratio_up =  MaxVariationUp->GetPointY(i) / VariationSet_vec[0]->GetBinContent(VariationSet_vec[0]->FindBin(MaxVariationUp->GetPointX(i)));
        double ratio_down =  MaxVariationDown->GetPointY(i) / VariationSet_vec[0]->GetBinContent(VariationSet_vec[0]->FindBin(MaxVariationDown->GetPointX(i)));

        std::cout<<"test, i: "<<i<<", ratio_up: "<<ratio_up<<", ratio_down: "<<ratio_down<<std::endl;

        if (!std::isnan(ratio_up)) {
            MaxVariationUp_Ratio -> SetPoint(
                MaxVariationUp_Ratio -> GetN(),
                MaxVariationUp->GetPointX(i),
                ratio_up
            );
            MaxVariationUp_Ratio -> SetPointError(
                MaxVariationUp_Ratio -> GetN() - 1,
                VariationSet_vec[0] -> GetBinWidth(1) / 2.,
                0
            );
        }

        if (!std::isnan(ratio_down)) {
            MaxVariationDown_Ratio -> SetPoint(
                MaxVariationDown_Ratio -> GetN(),
                MaxVariationDown->GetPointX(i),
                ratio_down
            );
            MaxVariationDown_Ratio -> SetPointError(
                MaxVariationDown_Ratio -> GetN() - 1,
                VariationSet_vec[0] -> GetBinWidth(1) / 2.,
                0
            );
        }
        

        MaxVariationUp->SetPoint(
            i, 
            MaxVariationUp->GetPointX(i), 
            fabs(
                MaxVariationUp->GetPointY(i) - 
                VariationSet_vec[0]->GetBinContent(VariationSet_vec[0]->FindBin(MaxVariationUp->GetPointX(i)))
            )
        ); // note : for example, 0.1

        MaxVariationDown->SetPoint(
            i, 
            MaxVariationDown->GetPointX(i), 
            fabs(
                MaxVariationDown->GetPointY(i) - 
                VariationSet_vec[0]->GetBinContent(VariationSet_vec[0]->FindBin(MaxVariationDown->GetPointX(i)))
            )
        ); // note : for example, 0.010

        std::cout<<"Index: "<<i<<", BinCenter: "<<MaxVariationUp->GetPointX(i)<<", MaxVariationUp: "<<MaxVariationUp->GetPointY(i)<<", MaxVariationDown: "<<MaxVariationDown->GetPointY(i)<<std::endl;

    }    

    TLegend * leg_maximal_variation = new TLegend(0.19,0.76,0.41,0.86);
    leg_maximal_variation -> SetBorderSize(0);
    leg_maximal_variation -> SetTextSize(0.03);
    leg_maximal_variation -> AddEntry(MaxVariationUp_Ratio, "Maximal Variation Up", "ep");
    leg_maximal_variation -> AddEntry(MaxVariationDown_Ratio, "Maximal Variation Down", "ep");

    c1 -> cd();
    MaxVariationUp_Ratio -> GetXaxis() -> SetLimits(VariationSet_vec[0]->GetXaxis()->GetXmin(), VariationSet_vec[0]->GetXaxis()->GetXmax());
    MaxVariationUp_Ratio -> GetYaxis() -> SetRangeUser(0,2.);
    MaxVariationUp_Ratio -> Draw("ap");
    MaxVariationDown_Ratio -> Draw("p same");
    line -> DrawLine(
        MaxVariationUp_Ratio->GetXaxis()->GetXmin(), 1, 
        MaxVariationUp_Ratio->GetXaxis()->GetXmax(), 1
    );
    leg_maximal_variation -> Draw("same");
    c1 -> Print(Form("%s/%s_MaximalRatio.pdf", final_output_directory.c_str(), plot_name_in.c_str()));
    // c1 -> Clear();

    // Division : -----------------------------------------------------------------------------------------------------------------------------------------------------

    TGraphAsymmErrors * FinalSyst_Asymmetry_gr = new TGraphAsymmErrors();
    for (int i = 0; i < MaxVariationUp->GetN(); i++){

        double x = MaxVariationUp->GetPointX(i);

        double y           = VariationSet_vec[0] -> GetBinContent(VariationSet_vec[0]->FindBin(x));
        double x_err = 0.5 * VariationSet_vec[0] -> GetBinWidth(1);

        double y_err_up = MaxVariationUp->GetPointY(i);
        double y_err_down = MaxVariationDown->GetPointY(i);


        FinalSyst_Asymmetry_gr -> SetPoint(FinalSyst_Asymmetry_gr->GetN(), x, y);
        FinalSyst_Asymmetry_gr -> SetPointError(FinalSyst_Asymmetry_gr->GetN() - 1, x_err, x_err, y_err_down, y_err_up);
    }

    FinalSyst_Asymmetry_gr -> SetMarkerStyle(20);
    FinalSyst_Asymmetry_gr -> SetMarkerSize(0);
    FinalSyst_Asymmetry_gr -> SetMarkerColorAlpha(1,0);
    FinalSyst_Asymmetry_gr -> SetFillColorAlpha(1,0.5);
    // FinalSyst_Asymmetry_gr -> SetLineColorAlpha(1,0);
    // FinalSyst_Asymmetry_gr -> SetLineWidth(1);
    FinalSyst_Asymmetry_gr->SetFillColorAlpha(92, 0.85);

    TLegend * leg_syst = new TLegend(0.4,0.7,0.8,0.78);
	leg_syst -> SetBorderSize(0);
	leg_syst -> SetMargin(0.15);
	leg_syst -> SetTextSize(0.03);

    TGraphErrors * gr_syst_style = new TGraphErrors();
    gr_syst_style->SetPoint(0,0,0);
    gr_syst_style->SetMarkerStyle(VariationSet_vec[0]->GetMarkerStyle());
    gr_syst_style->SetMarkerSize(VariationSet_vec[0]->GetMarkerSize());
    gr_syst_style->SetMarkerColor(VariationSet_vec[0]->GetMarkerColor());
    gr_syst_style->SetLineWidth(0);
    gr_syst_style->SetLineColorAlpha(1,0);
    gr_syst_style->SetFillColorAlpha(FinalSyst_Asymmetry_gr->GetFillColor(), 0.92);


    leg_syst->AddEntry(gr_syst_style, "Data", "epf");

    c1 -> cd();
    
    VariationSet_vec[0] -> GetXaxis() -> SetLabelOffset(0);
    VariationSet_vec[0] -> GetXaxis() -> SetTitleOffset(0);

    VariationSet_vec[0]->Draw("ep");
    FinalSyst_Asymmetry_gr->Draw("P2 same");
    VariationSet_vec[0]->Draw("ep same");

    leg_syst->Draw("same");
    ltx->DrawLatex(0.38, 0.83, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021")); 
    c1 ->Print(Form("%s/%s_StatSystUnc.pdf", final_output_directory.c_str(), plot_name_in.c_str() )); 
    c1 -> Clear();

    // Division : -----------------------------------------------------------------------------------------------------------------------------------------------------

    TGraphAsymmErrors * Final_Asymmetry_gr = new TGraphAsymmErrors();

    for (int i = 0; i < MaxVariationUp->GetN(); i++){

        double x = MaxVariationUp->GetPointX(i);

        double y           = VariationSet_vec[0] -> GetBinContent(VariationSet_vec[0]->FindBin(x));
        double x_err = 0.5 * VariationSet_vec[0] -> GetBinWidth(1);
        
        double y_err_stat = VariationSet_vec[0] -> GetBinError(VariationSet_vec[0]->FindBin(x));

        double y_err_syst_up = MaxVariationUp->GetPointY(i);
        double y_err_syst_down = MaxVariationDown->GetPointY(i);

        double y_err_final_up = sqrt(y_err_stat * y_err_stat + y_err_syst_up * y_err_syst_up);
        double y_err_final_down = sqrt(y_err_stat * y_err_stat + y_err_syst_down * y_err_syst_down);

        if (x < 8 || x > 176) {continue;}

        Final_Asymmetry_gr -> SetPoint(Final_Asymmetry_gr->GetN(), x, y);
        Final_Asymmetry_gr -> SetPointError(Final_Asymmetry_gr->GetN() - 1, x_err, x_err, y_err_final_down, y_err_final_up);

        std::cout<<"("<<x<<", "<<y<<"), ("<<x<<", "<<y - y_err_final_down<<"), ("<<x<<", "<<y + y_err_final_up<<")"<<std::endl;
    }

    
    Final_Asymmetry_gr->SetMarkerStyle(VariationSet_vec[0]->GetMarkerStyle());
    Final_Asymmetry_gr->SetMarkerSize(VariationSet_vec[0]->GetMarkerSize());
    Final_Asymmetry_gr->SetMarkerColor(VariationSet_vec[0]->GetMarkerColor());

    Final_Asymmetry_gr->GetXaxis()->SetTitle("DAC value");
    Final_Asymmetry_gr->GetYaxis()->SetTitle("Entries (A.U.)");
    Final_Asymmetry_gr->GetXaxis()->SetLimits(VariationSet_vec[0]->GetXaxis()->GetXmin(), VariationSet_vec[0]->GetXaxis()->GetXmax());
    Final_Asymmetry_gr->GetYaxis()->SetRangeUser(0, Yaxis_max);



    TF1 *LandGausComp = new TF1("LandGausComp", Func_Langaus, 0, 200,4);
	LandGausComp -> SetLineColor(65);
	LandGausComp -> SetLineStyle(2);
	LandGausComp -> SetLineWidth(2);

	TF1 * BkgComp = new TF1("BkgComp", SingleExp, 0, 200,3);
	BkgComp -> SetLineColor(93);
	BkgComp -> SetLineStyle(2);
	BkgComp -> SetLineWidth(2);

	TF1 * BkgComp2 = new TF1("BkgComp2", SingleGaus, 0, 200,3);
	BkgComp2 -> SetLineColor(30);
	BkgComp2 -> SetLineStyle(2);
	BkgComp2 -> SetLineWidth(2);

    double exp_par4 = (plot_name_in.find("L2") != plot_name_in.npos) ? 0.1 : 1.02518e+04;
    double exp_par5 = (plot_name_in.find("L2") != plot_name_in.npos) ? -0.64 : -8.05366e-01;
    double exp_par6 = (plot_name_in.find("L2") != plot_name_in.npos) ? -24.7 : -7.79397e+00;

    double gaus_par7 = (plot_name_in.find("L2") != plot_name_in.npos) ? 1.41845e+02 : 1.35020e+02;
    double gaus_par8 = (plot_name_in.find("L2") != plot_name_in.npos) ? 1.81768e+01 : 1.90138e+01;

	TF1 * land_gaus_bkg_fit = new TF1("land_gaus_bkg_fit", Func_Langaus_bkg, 0, 200, 9);
	land_gaus_bkg_fit -> SetLineColor(TColor::GetColor("#F5321D"));
	land_gaus_bkg_fit -> SetLineWidth(2);
	land_gaus_bkg_fit -> SetParameters(
		// 3.52106e+00, 7.06707e+01 ,  4.36460e+04, 1.28136e+01,
        3.25141e+00, 7.30421e+01, 4.39288e+04, 1.38576e+01,
		exp_par4, exp_par5, exp_par6,
		gaus_par7, gaus_par8
	);

	land_gaus_bkg_fit -> SetParLimits(7,0,100000);
	// land_gaus_bkg_fit -> SetParLimits(8,0,20);

    TLegend * leg_final = new TLegend(0.55,0.68,0.8,0.8);
	leg_final -> SetBorderSize(0);
	leg_final -> SetMargin(0.2);
	leg_final -> SetTextSize(0.025);

	leg_final->AddEntry(land_gaus_bkg_fit, "Total Fit", "l");
	leg_final->AddEntry(LandGausComp, "Signal (Landau #otimes Gaussian)", "l");
	leg_final->AddEntry(BkgComp, "Background (Exponential)", "l");
	leg_final->AddEntry(BkgComp2, "Background (Gaussian)", "l");

    TLatex * ltx_fit = new TLatex();
    ltx_fit->SetNDC();
    ltx_fit->SetTextSize(0.03);






    Final_Asymmetry_gr -> Fit(land_gaus_bkg_fit, "", "", fit_range_pair.first, fit_range_pair.second);

    LandGausComp->SetParameters(
        land_gaus_bkg_fit->GetParameter(0),
        land_gaus_bkg_fit->GetParameter(1),
        land_gaus_bkg_fit->GetParameter(2),
        land_gaus_bkg_fit->GetParameter(3)
    );

    BkgComp->SetParameters(
        land_gaus_bkg_fit->GetParameter(4),
        land_gaus_bkg_fit->GetParameter(5),
        land_gaus_bkg_fit->GetParameter(6)
    );

    BkgComp2->SetParameters(
        land_gaus_bkg_fit->GetParameter(7),
        0,
        land_gaus_bkg_fit->GetParameter(8)
        // land_gaus_bkg_fit->GetParameter(9)
    );

    c1 -> cd();
    Final_Asymmetry_gr -> Draw("ap");

    land_gaus_bkg_fit->Draw("lsame");
    LandGausComp->Draw("lsame");
    BkgComp->Draw("lsame");
    BkgComp2->Draw("lsame");
    


    double signal_purity = LandGausComp->Integral(signal_purity_range.first, signal_purity_range.second)/land_gaus_bkg_fit->Integral(signal_purity_range.first, signal_purity_range.second);
    cout<<"ChiSquare: "<<land_gaus_bkg_fit->GetChisquare()<<", NDF: "<<land_gaus_bkg_fit->GetNDF()<<", ReducedChiSquare: "<<land_gaus_bkg_fit->GetChisquare()/land_gaus_bkg_fit->GetNDF()<<endl;
    
    cout<<Form("Full Integral {%.0f,%.0f}: ",signal_purity_range.first, signal_purity_range.second)<<land_gaus_bkg_fit->Integral(signal_purity_range.first, signal_purity_range.second)
    <<Form(", LandGausComp Integral {%.0f, %.0f}: ",signal_purity_range.first,signal_purity_range.second)<<LandGausComp->Integral(signal_purity_range.first, signal_purity_range.second)
    <<", Signal purity: "<<signal_purity<<endl;

    cout<<Form("From 15 to 30, the integral: %.4f, From 15 to 210: %.4f", LandGausComp->Integral(15, 30), LandGausComp->Integral(15, 210))<<", ratio: "<<LandGausComp->Integral(15, 30)/LandGausComp->Integral(15, 210)<<endl;
    // cout<<Form("From 15 to 30, the integralError: %.4f, From 15 to 210: %.4f", LandGausComp->IntegralError(15, 30), LandGausComp->IntegralError(15, 210))<<endl;

    cout<<Form("From 0 to 15, the integral: %.4f, From 0 to 210: %.4f", LandGausComp->Integral(0, 15), LandGausComp->Integral(0, 210))<<", ratio: "<<LandGausComp->Integral(0, 15)/LandGausComp->Integral(0, 210)<<endl;
    // cout<<Form("From 0 to 15, the integralError: %.4f, From 0 to 210: %.4f", LandGausComp->IntegralError(0, 15), LandGausComp->IntegralError(0, 210))<<endl;
    
    cout<<endl;

    double ltx_fit_Y_offset = (ShowReducedChi2) ? 0 : 0.04;
    if (ShowReducedChi2){ltx_fit->DrawLatex(0.5, 0.62, Form("#chi^{2} / NDF : %.2f / %d = %.2f",land_gaus_bkg_fit->GetChisquare(),land_gaus_bkg_fit->GetNDF(), land_gaus_bkg_fit->GetChisquare()/land_gaus_bkg_fit->GetNDF()));}
    ltx_fit->DrawLatex(0.5, 0.58 + ltx_fit_Y_offset, Form("Landau MPV: %.2f #pm %.2f",land_gaus_bkg_fit->GetParameter(1), land_gaus_bkg_fit->GetParError(1)));
    ltx_fit->DrawLatex(0.5, 0.54 + ltx_fit_Y_offset, Form("Signal Purity (DAC #geq %.0f): %.2f%",signal_purity_range.first, signal_purity * 100.));

    leg_final->Draw("same");

    ltx->DrawLatex(0.38, 0.83, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021")); 

    c1 ->Print(Form("%s/%s_Final.pdf", final_output_directory.c_str(), plot_name_in.c_str() )); 
    // c1 -> Clear();











    // // note : for the "all" case
    // double all_MaxVariationUp   = VariationSet_vec[0].all_Final_Effi;
    // double all_MaxVariationDown = VariationSet_vec[0].all_Final_Effi;

    // for (int i = 1; i < VariationSet_vec.size(); i++) {
    //     if (VariationSet_vec[i].all_Final_Effi > all_MaxVariationUp) {
    //         all_MaxVariationUp = VariationSet_vec[i].all_Final_Effi;
    //     }
    //     if (VariationSet_vec[i].all_Final_Effi < all_MaxVariationDown) {
    //         all_MaxVariationDown = VariationSet_vec[i].all_Final_Effi;
    //     }
    // }

    // all_MaxVariationUp = fabs(all_MaxVariationUp / VariationSet_vec[0].all_Final_Effi - 1.0);
    // all_MaxVariationDown = fabs(all_MaxVariationDown / VariationSet_vec[0].all_Final_Effi - 1.0);

    // std::cout<<"leg_header_in : "<<leg_header_in<<std::endl;
    // std::cout<<"all_MaxVariationUp: "<<all_MaxVariationUp<<std::endl;
    // std::cout<<"all_MaxVariationDown: "<<all_MaxVariationDown<<std::endl;
    // std::cout<<"MaxVariationUp in Y axis: ";
    // for (int i = 0; i < MaxVariationUp->GetN(); i++) {
    //     std::cout<<Form("%.5f", MaxVariationUp->GetPointY(i))<<", ";
    // }
    // std::cout<<std::endl;

    // std::cout<<"MaxVariationDown in Y axis: ";
    // for (int i = 0; i < MaxVariationDown->GetN(); i++) {
    //     std::cout<<Form("%.5f", MaxVariationDown->GetPointY(i))<<", ";
    // }
    // std::cout<<std::endl;
    // std::cout<<std::endl;

    // return {
    //     MaxVariationUp,
    //     MaxVariationDown,

    //     all_MaxVariationUp,
    //     all_MaxVariationDown
    // };

}


int combine_2(){

    std::string input_directory_1 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out";
    std::string input_filename_1 = "DAC_scan_out.root";

    std::string input_directory_2 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out_Scan2Scan3Match1Bin";
    std::string input_filename_2 = "DAC_scan_out.root";

    std::string input_directory_3 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out_UsedOverflow";
    std::string input_filename_3 = "DAC_scan_out.root";

    std::string input_directory_4 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out_AllOnly6thBin";
    std::string input_filename_4 = "DAC_scan_out.root";

    // SetsPhenixStyle();

    // gStyle->SetTitleBorderSize(0);
    // gStyle->SetTitleAlign(22);
    // gStyle->SetTitleX(0.5);
	// gStyle->SetTitleY(0.95);
	// gStyle->SetTitleFillColor(10);
	// gStyle->SetPadRightMargin(0.09);
    // gStyle->SetPadTopMargin(0.1);

    TCanvas * c1 = new TCanvas("c1", "c1", 950, 800);

    TFile * file_in_1 = TFile::Open(Form("%s/%s", input_directory_1.c_str(), input_filename_1.c_str()));
    TFile * file_in_2 = TFile::Open(Form("%s/%s", input_directory_2.c_str(), input_filename_2.c_str()));
    TFile * file_in_3 = TFile::Open(Form("%s/%s", input_directory_3.c_str(), input_filename_3.c_str()));
    TFile * file_in_4 = TFile::Open(Form("%s/%s", input_directory_4.c_str(), input_filename_4.c_str()));

    TH1D * h1D_DAC_hist_all_1[3];
    TH1D * h1D_DAC_hist_all_2[3];
    TH1D * h1D_DAC_hist_all_3[3];
    TH1D * h1D_DAC_hist_all_4[3];

    h1D_DAC_hist_all_1[0] = (TH1D*) file_in_1 -> Get("DAC_hist_all_0"); // h1D_DAC_hist_all_1[0] -> SetTitle("DAC_hist_all_0");
    h1D_DAC_hist_all_1[1] = (TH1D*) file_in_1 -> Get("DAC_hist_all_1"); // h1D_DAC_hist_all_1[1] -> SetTitle("DAC_hist_all_1");
    h1D_DAC_hist_all_1[2] = (TH1D*) file_in_1 -> Get("DAC_hist_all_2"); // h1D_DAC_hist_all_1[2] -> SetTitle("DAC_hist_all_2");

    h1D_DAC_hist_all_2[0] = (TH1D*) file_in_2 -> Get("DAC_hist_all_0"); // h1D_DAC_hist_all_2[0] -> SetTitle("DAC_hist_all_0");
    h1D_DAC_hist_all_2[1] = (TH1D*) file_in_2 -> Get("DAC_hist_all_1"); // h1D_DAC_hist_all_2[1] -> SetTitle("DAC_hist_all_1");
    h1D_DAC_hist_all_2[2] = (TH1D*) file_in_2 -> Get("DAC_hist_all_2"); // h1D_DAC_hist_all_2[2] -> SetTitle("DAC_hist_all_2");

    h1D_DAC_hist_all_3[0] = (TH1D*) file_in_3 -> Get("DAC_hist_all_0"); // h1D_DAC_hist_all_3[0] -> SetTitle("DAC_hist_all_0");
    h1D_DAC_hist_all_3[1] = (TH1D*) file_in_3 -> Get("DAC_hist_all_1"); // h1D_DAC_hist_all_3[1] -> SetTitle("DAC_hist_all_1");
    h1D_DAC_hist_all_3[2] = (TH1D*) file_in_3 -> Get("DAC_hist_all_2"); // h1D_DAC_hist_all_3[2] -> SetTitle("DAC_hist_all_2");

    h1D_DAC_hist_all_4[0] = (TH1D*) file_in_4 -> Get("DAC_hist_all_0"); // h1D_DAC_hist_all_4[0] -> SetTitle("DAC_hist_all_0");
    h1D_DAC_hist_all_4[1] = (TH1D*) file_in_4 -> Get("DAC_hist_all_1"); // h1D_DAC_hist_all_4[1] -> SetTitle("DAC_hist_all_1");
    h1D_DAC_hist_all_4[2] = (TH1D*) file_in_4 -> Get("DAC_hist_all_2"); // h1D_DAC_hist_all_4[2] -> SetTitle("DAC_hist_all_2");


    GetVariation(
        {
            h1D_DAC_hist_all_1[0],
            h1D_DAC_hist_all_2[0],
            h1D_DAC_hist_all_3[0],
            h1D_DAC_hist_all_4[0]
        },
        {
            "Baseline (two bins for dist. matching)",
            "Scan2 matched by only the 6th bin",
            "Overflow bins included for matching",
            "Only first overlapped bin used"
        },
        "Different matching methods",
        input_directory_4,
        "Variation_match_method_L0"

        // std::vector<TH1D*> VariationSet_vec, // note : the first one has to be the baseline
        // std::vector<std::string> file_title_vec,
        // std::string leg_header_in,
        // std::string final_output_directory,
        // std::string plot_name_in
    );


    GetVariation(
        {
            h1D_DAC_hist_all_1[1],
            h1D_DAC_hist_all_2[1],
            h1D_DAC_hist_all_3[1],
            h1D_DAC_hist_all_4[1]
        },
        {
            "Baseline (two bins for dist. matching)",
            "Scan2 matched by only the 6th bin",
            "Overflow bins included for matching",
            "Only first overlapped bin used"
        },
        "Different matching methods",
        input_directory_4,
        "Variation_match_method_L1"

        // std::vector<TH1D*> VariationSet_vec, // note : the first one has to be the baseline
        // std::vector<std::string> file_title_vec,
        // std::string leg_header_in,
        // std::string final_output_directory,
        // std::string plot_name_in
    );

    GetVariation(
        {
            h1D_DAC_hist_all_1[2],
            h1D_DAC_hist_all_2[2],
            h1D_DAC_hist_all_3[2],
            h1D_DAC_hist_all_4[2]
        },
        {
            "Baseline (two bins for dist. matching)",
            "Scan2 matched by only the 6th bin",
            "Overflow bins included for matching",
            "Only first overlapped bin used"
        },
        "Different matching methods",
        input_directory_4,
        "Variation_match_method_L2"

        // std::vector<TH1D*> VariationSet_vec, // note : the first one has to be the baseline
        // std::vector<std::string> file_title_vec,
        // std::string leg_header_in,
        // std::string final_output_directory,
        // std::string plot_name_in
    );

    return 888;
}