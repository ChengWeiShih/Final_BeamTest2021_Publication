#include "sPhenixStyle.C"


double Gaus_function (double *x, double*par)
{
	double gaussian_eq  = TMath::Exp( -1*( pow(x[0],2)/(2*pow(par[0],2)) ) );
	return gaussian_eq;
}

double Exp_bkg_func (double *x, double *par)
{
	return par[0] * TMath::Exp(par[1] * (x[0] + par[2])) + par[3];  
	// TMath::Exp(par[0] * TMath(x[0]+par[1])) + par[2];
}

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
	
	// double exp_bkg_eq = par[4] * TMath::Exp(par[5] * (x[0] + par[6]));
	// double gaus_bkg_eq = par[7]*TMath::Gaus(x[0],0,par[8]);

	double gaus_bkg_eq1 = par[4]*TMath::Gaus(x[0],0,par[5]);
	double gaus_bkg_eq2 = par[6]*TMath::Gaus(x[0],0,par[7]);

	// return signal + exp_bkg_eq + gaus_bkg_eq;
	return signal + gaus_bkg_eq1 + gaus_bkg_eq2;
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

void SoNRatio_Final(
    double MPV,
    double MPV_error,
    TF1 * BkgComp1_in, // note : skinny
    TF1 * BkgComp2_in,  // note : the fat one 

    std::string output_directory_in
){
    TF1 * temp_BKgComp1 = new TF1("",SingleGaus,-200,200, 3);
    TF1 * temp_BKgComp2 = new TF1("",SingleGaus,-200,200, 3);

    TRandom3 * rand = new TRandom3(0);

    std::cout<<"MPV: "<<MPV<<", MPV_error: "<<MPV_error<<std::endl;
    
    std::cout<<"BkgComp1_in, par0: "<<BkgComp1_in->GetParameter(0)<<" #pm "<<BkgComp1_in->GetParError(0)<<std::endl;
    std::cout<<"BkgComp1_in, par1: "<<BkgComp1_in->GetParameter(1)<<" #pm "<<BkgComp1_in->GetParError(1)<<std::endl;
    std::cout<<"BkgComp1_in, par2: "<<BkgComp1_in->GetParameter(2)<<" #pm "<<BkgComp1_in->GetParError(2)<<std::endl;

    std::cout<<"BkgComp2_in, par0: "<<BkgComp2_in->GetParameter(0)<<" #pm "<<BkgComp2_in->GetParError(0)<<std::endl;
    std::cout<<"BkgComp2_in, par1: "<<BkgComp2_in->GetParameter(1)<<" #pm "<<BkgComp2_in->GetParError(1)<<std::endl;
    std::cout<<"BkgComp2_in, par2: "<<BkgComp2_in->GetParameter(2)<<" #pm "<<BkgComp2_in->GetParError(2)<<std::endl;

    TH1D * h1D_SoN_hist = new TH1D("h1D_SoN_hist",";Signal-to-Noise ratio (S/N);Entries",100, 0, 30);
    TH1D * h1D_BkgComp_width = new TH1D("h1D_BkgComp_width",";Weighted average noise width;Entries",100, 0, 10);


    for (int i = 0; i < 10000; i++) {

        double this_MPV = rand->Gaus(MPV, MPV_error);

        double this_BkgComp1_Height = rand -> Gaus(BkgComp1_in -> GetParameter(0), BkgComp1_in -> GetParError(0) );
        double this_BkgComp1_Mean = 0;
        double this_BkgCOmp1_Width = rand -> Gaus(BkgComp1_in -> GetParameter(2), BkgComp1_in -> GetParError(2) );

        double this_BkgComp2_Height = rand -> Gaus(BkgComp2_in -> GetParameter(0), BkgComp2_in -> GetParError(0) );
        double this_BkgComp2_Mean = 0;
        double this_BkgCOmp2_Width = BkgComp2_in -> GetParameter(2); // rand -> Gaus(BkgComp2_in -> GetParameter(2), BkgComp2_in -> GetParError(2) );

        temp_BKgComp1 -> SetParameters(
            this_BkgComp1_Height,
            this_BkgComp1_Mean,
            this_BkgCOmp1_Width
        );

        temp_BKgComp2 -> SetParameters(
            this_BkgComp2_Height,
            this_BkgComp2_Mean,
            this_BkgCOmp2_Width
        );

        
        double BkgComp1Width    = temp_BKgComp1 -> GetParameter(2);
        double BkgComp1Integral = temp_BKgComp1 -> Integral(-60,60);

        double BkgComp2Width    = temp_BKgComp2 -> GetParameter(2);
        double BkgComp2Integral = temp_BKgComp2 -> Integral(-60,60);

        double weight_average_bkg_GausWidth = (BkgComp1Width * BkgComp1Integral + BkgComp2Width * BkgComp2Integral) / (BkgComp1Integral + BkgComp2Integral);

        double SoN_ratio = this_MPV / weight_average_bkg_GausWidth;

        h1D_BkgComp_width -> Fill(weight_average_bkg_GausWidth);

        h1D_SoN_hist -> Fill(SoN_ratio);
    }

    TCanvas * c2 = new TCanvas("c2", "c2", 950, 800);
    c2 -> cd();

    TLatex * ltx_fit = new TLatex();
    ltx_fit->SetNDC();
    ltx_fit->SetTextSize(0.03);

    h1D_BkgComp_width -> SetMaximum( h1D_BkgComp_width -> GetBinContent(h1D_BkgComp_width -> GetMaximumBin()) * 1.5 );
    h1D_BkgComp_width -> Draw("hist");
    h1D_BkgComp_width -> Fit("gaus");
    h1D_BkgComp_width -> GetFunction("gaus") -> SetLineColor(kRed);
    h1D_BkgComp_width -> GetFunction("gaus") -> SetLineWidth(1);
    h1D_BkgComp_width -> GetFunction("gaus") -> SetNpx(1000);
    h1D_BkgComp_width -> GetFunction("gaus") -> Draw("lsame");
    ltx_fit -> DrawLatex(0.3, 0.8, Form("Gaus Mean: %.2f", h1D_BkgComp_width -> GetFunction("gaus") -> GetParameter(1) ));
    ltx_fit -> DrawLatex(0.3, 0.75, Form("Gaus Sigma: %.2f", h1D_BkgComp_width -> GetFunction("gaus") -> GetParameter(2) ));
    c2 -> Print(Form("%s/h1D_BkgComp_width.pdf", output_directory_in.c_str() ));
    c2 -> Clear();

    c2 -> cd();
    h1D_SoN_hist -> SetMaximum( h1D_SoN_hist -> GetBinContent(h1D_SoN_hist -> GetMaximumBin()) * 1.5 );
    h1D_SoN_hist -> Draw("hist");
    h1D_SoN_hist -> Fit("gaus");
    h1D_SoN_hist -> GetFunction("gaus") -> SetLineColor(kRed);
    h1D_SoN_hist -> GetFunction("gaus") -> SetLineWidth(1);
    h1D_SoN_hist -> GetFunction("gaus") -> SetNpx(1000);
    h1D_SoN_hist -> GetFunction("gaus") -> Draw("lsame");
    ltx_fit -> DrawLatex(0.3, 0.8, Form("Gaus Mean: %.2f", h1D_SoN_hist -> GetFunction("gaus") -> GetParameter(1) ));
    ltx_fit -> DrawLatex(0.3, 0.75, Form("Gaus Sigma: %.2f", h1D_SoN_hist -> GetFunction("gaus") -> GetParameter(2) ));
    c2 -> Print(Form("%s/h1D_SoN_hist.pdf", output_directory_in.c_str() ));
    c2 -> Clear();
}



int LbLVariation(){

	std::pair<double,double> fit_range_pair = {8, 160};
    double Yaxis_max = 2500;
	bool ShowReducedChi2 = false;


    std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out_AllOnly6thBin";
    std::string input_filename = "DAC_scan_out.root";

    SetsPhenixStyle();

    gStyle -> SetPadRightMargin(0.09);

    TLegend * legend = new TLegend(0.6, 0.74, 0.9, 0.84);
    legend -> SetNColumns(2);

    TLegend * leg_final = new TLegend(0.55,0.63,0.8,0.82);
	leg_final -> SetBorderSize(0);
	leg_final -> SetMargin(0.2);
	leg_final -> SetTextSize(0.025);

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    // ltx->SetTextAlign(31);

    TFile * file_in = TFile::Open(Form("%s/%s", input_directory.c_str(), input_filename.c_str()));
    
    TH1D * l0_hist = (TH1D*) file_in -> Get("DAC_hist_all_0");
    TH1D * l1_hist = (TH1D*) file_in -> Get("DAC_hist_all_1");
    TH1D * l2_hist = (TH1D*) file_in -> Get("DAC_hist_all_2");

    TF1 * l0_fit = (TF1*) file_in -> Get("land_gaus_bkg_fit_0");
    TF1 * l1_fit = (TF1*) file_in -> Get("land_gaus_bkg_fit_1");
    TF1 * l2_fit = (TF1*) file_in -> Get("land_gaus_bkg_fit_2");

    l0_hist -> SetLineColor(kRed);
    l0_hist -> SetMarkerColor(kRed);

    l2_hist -> SetLineColor(kBlue);
    l2_hist -> SetMarkerColor(kBlue);

    l0_fit -> SetLineColor(kRed-9);
    l1_fit -> SetLineColor(kGray);
    l2_fit -> SetLineColor(kBlue-9);

    legend -> AddEntry(l0_hist, "L0", "lep");
    legend -> AddEntry(l0_fit, "L0-fit", "l");

    // legend -> AddEntry(l1_hist, "L1", "lep");
    // legend -> AddEntry(l1_fit, "L1-fit", "l");

    legend -> AddEntry(l2_hist, "L2", "lep");
    legend -> AddEntry(l2_fit, "L2-fit", "l");

    TCanvas * c1 = new TCanvas("c1", "c1", 950, 800);
    c1 -> cd();

    l0_hist -> SetMaximum(Yaxis_max);

    l0_hist -> Draw("ep");
    // l1_hist -> Draw("epsame");
    l2_hist -> Draw("epsame");

    l0_fit -> Draw("same");
    // l1_fit -> Draw("same");
    l2_fit -> Draw("same");

    legend -> Draw();

    ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));

    c1 -> Print(Form("%s/LbLVariation_%.0f.pdf", input_directory.c_str(), Yaxis_max));
    c1 -> Clear();

    c1 -> cd();

    l0_hist -> SetMaximum(8000);

    l0_hist -> Draw("ep");
    // l1_hist -> Draw("epsame");
    l2_hist -> Draw("epsame");

    l0_fit -> Draw("same");
    // l1_fit -> Draw("same");
    l2_fit -> Draw("same");

    legend -> Draw();

    ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));

    c1 -> Print(Form("%s/%s", input_directory.c_str(), "LbLVariation_8000.pdf"));
    c1 -> Clear();


    //Division: --------------------------------------------------------------------------------------------------------------------------------------------
    TH1D * Final_hist = new TH1D(
        "Final_hist", 
        "Final_hist", 
        l0_hist -> GetNbinsX(),
        l0_hist -> GetXaxis() -> GetXmin(),
        l0_hist -> GetXaxis() -> GetXmax()
    );

    Final_hist -> SetMarkerStyle(20);
    Final_hist -> SetMarkerSize(0.8);
    Final_hist -> SetMarkerColor(1);
    Final_hist -> SetLineColor(1);
    Final_hist -> SetLineWidth(1);
    Final_hist -> GetXaxis() -> SetTitle("DAC value");
    Final_hist -> GetYaxis() -> SetTitle("Entries (A.U.)");

    TGraphAsymmErrors * Final_syst = new TGraphAsymmErrors();
    Final_syst -> SetMarkerStyle(20);
    Final_syst -> SetMarkerSize(0);
    Final_syst -> SetMarkerColorAlpha(1,0);
    // Final_syst -> SetLineColorAlpha(1,0);
    // Final_syst -> SetLineWidth(1);
    Final_syst->SetFillColorAlpha(92, 0.85);
    Final_syst -> GetXaxis() -> SetTitle("DAC value");
    Final_syst -> GetYaxis() -> SetTitle("Entries (A.U.)");


    for (int i = 1; i <= l0_hist -> GetNbinsX(); i++){
        
        double l0_center = l0_hist -> GetBinCenter(i);
        double l1_width = l1_hist -> GetBinWidth(i) / 2.;
        double l0_content = l0_hist -> GetBinContent(i);
        double l0_error = l0_hist -> GetBinError(i);
        double l0_weight = 1./pow(l0_error,2);


        double l2_content = l2_hist -> GetBinContent(i);
        double l2_error = l2_hist -> GetBinError(i);
        double l2_weight = 1./pow(l2_error,2);

        if (l0_center <= 8) {continue;}
        if (l0_center > 176) {continue;}

        double Final_content = (l0_content * l0_weight + l2_content * l2_weight) / (l0_weight + l2_weight);
        double Final_error = 1./sqrt(l0_weight + l2_weight);

        double upper_one = std::max(l0_content, l2_content);
        double lower_one = std::min(l0_content, l2_content);

        Final_hist -> SetBinContent(i, Final_content);
        Final_hist -> SetBinError(i, Final_error);

        std::cout<<"bin: "<<i<<" l0_center: "<<l0_center
                 <<" l0_content: "<<l0_content
                 <<" l2_content: "<<l2_content
                 <<" Final_content: "<<Final_content
                 <<" Final_error: "<<Final_error
                 <<" upper_one: "<<upper_one
                 <<" lower_one: "<<lower_one<<std::endl;

        Final_syst -> SetPoint(
            Final_syst->GetN(),
            l0_center,
            Final_content
        );

        Final_syst -> SetPointError(
            Final_syst->GetN()-1,
            l1_width,
            l1_width,
            (Final_content - lower_one),
            (upper_one - Final_content)
        );

    }

    Final_syst -> GetXaxis() -> SetLimits(
        Final_hist -> GetXaxis() -> GetXmin(),
        Final_hist -> GetXaxis() -> GetXmax()
    );
    Final_syst -> GetYaxis() -> SetRangeUser(
        0,
        Yaxis_max
    );

    c1 -> cd();
    Final_syst -> Draw("AP2");
    Final_hist -> Draw("ep same");


    //Division: --------------------------------------------------------------------------------------------------------------------------------------------

    

    TF1 *LandGausComp = new TF1("LandGausComp", Func_Langaus, 0, 200,4);
	LandGausComp -> SetLineColor(61);
	LandGausComp -> SetLineStyle(2);
	LandGausComp -> SetLineWidth(2);

	// TF1 * BkgComp = new TF1("BkgComp", SingleExp, 0, 200,3);
	// BkgComp -> SetLineColor(93);
	// BkgComp -> SetLineStyle(2);
	// BkgComp -> SetLineWidth(2);

	TF1 * BkgComp = new TF1("BkgComp", SingleGaus, 0, 200,3);
	BkgComp -> SetLineColor(30);
	BkgComp -> SetLineStyle(2);
	BkgComp -> SetLineWidth(2);

	TF1 * BkgComp2 = new TF1("BkgComp2", SingleGaus, 0, 200,3);
	BkgComp2 -> SetLineColor(6);
	BkgComp2 -> SetLineStyle(2);
	BkgComp2 -> SetLineWidth(2);

	TF1 * land_gaus_bkg_fit = new TF1("land_gaus_bkg_fit", Func_Langaus_bkg, 0, 200, 8);
	land_gaus_bkg_fit -> SetLineColor(TColor::GetColor("#F5321D"));
	land_gaus_bkg_fit -> SetLineWidth(2);
	land_gaus_bkg_fit -> SetParameters(
		3.47183e+00, 7.19476e+01,  4.14353e+04, 1.32139e+01,
		3.52018e+04, 3.94777e+00,
		1.32515e+02, 2.00000e+01
	);

    land_gaus_bkg_fit -> SetParLimits(4,0,100000000);
	land_gaus_bkg_fit -> SetParLimits(5,0,20);

	land_gaus_bkg_fit -> SetParLimits(6,0,100000000);
	land_gaus_bkg_fit -> SetParLimits(7,0,20);

    TGraphErrors * gr_syst_style = new TGraphErrors();
    gr_syst_style->SetPoint(0,0,0);
    gr_syst_style->SetMarkerStyle(Final_hist->GetMarkerStyle());
    gr_syst_style->SetMarkerSize(Final_hist->GetMarkerSize());
    gr_syst_style->SetMarkerColor(Final_hist->GetMarkerColor());
    gr_syst_style->SetLineWidth(0);
    gr_syst_style->SetLineColorAlpha(1,0);
    gr_syst_style->SetFillColorAlpha(Final_syst->GetFillColor(), 0.92);

    leg_final->AddEntry(gr_syst_style, "Data","epf");
	leg_final->AddEntry(land_gaus_bkg_fit, "Total fit", "l");
	leg_final->AddEntry(LandGausComp, "Signal (Landau #otimes Gaussian)", "l");
	leg_final->AddEntry(BkgComp, "Background (Gaussian1)", "l");
	leg_final->AddEntry(BkgComp2, "Background (Gaussian2)", "l");

    Final_hist -> Fit(land_gaus_bkg_fit, "", "", fit_range_pair.first, fit_range_pair.second);


    LandGausComp->SetParameters(
        land_gaus_bkg_fit->GetParameter(0),
        land_gaus_bkg_fit->GetParameter(1),
        land_gaus_bkg_fit->GetParameter(2),
        land_gaus_bkg_fit->GetParameter(3)
    );

    BkgComp->SetParameters(
        land_gaus_bkg_fit->GetParameter(4),
        0.,
        land_gaus_bkg_fit->GetParameter(5)
        // land_gaus_bkg_fit->GetParameter(9)
    );

    BkgComp->SetParError(0, land_gaus_bkg_fit->GetParError(4) );
    BkgComp->SetParError(1, 0.);
    BkgComp->SetParError(2, land_gaus_bkg_fit->GetParError(5) );

    BkgComp2->SetParameters(
        land_gaus_bkg_fit->GetParameter(6),
        0.,
        land_gaus_bkg_fit->GetParameter(7)
        // land_gaus_bkg_fit->GetParameter(9)
    );

    BkgComp2->SetParError(0, land_gaus_bkg_fit->GetParError(6) );
    BkgComp2->SetParError(1, 0.);
    BkgComp2->SetParError(2, land_gaus_bkg_fit->GetParError(7) );

    // land_gaus_bkg_fit -> Draw("lsame");
    LandGausComp->Draw("lsame");
    BkgComp->Draw("lsame");
    BkgComp2->Draw("lsame");

    // DAC_hist_all[i1]->Draw("ep same");

    leg_final->Draw("same");

    ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));


    double Gaus1Width    = BkgComp->GetParameter(2);
    double Gaus1Integral = BkgComp->Integral(-60,60);

    double Gaus2Width    = BkgComp2->GetParameter(2);
    double Gaus2Integral = BkgComp2->Integral(-60,60);

    double weight_average_GausWidth = (Gaus1Width * Gaus1Integral + Gaus2Width * Gaus2Integral) / (Gaus1Integral + Gaus2Integral);
    
    std::cout<<"Bkg Gaus1 width: "<<Gaus1Width<<", Integral: "<<Gaus1Integral<<std::endl;
    std::cout<<"Bkg Gaus2 width: "<<Gaus2Width<<", Integral: "<<Gaus2Integral<<std::endl;
    std::cout<<"weight_average_Gaus_width: "<<weight_average_GausWidth<<std::endl;

    // std::cout<<"Land MPV: "<<land_gaus_bkg_fit->GetParameter(1)<<std::endl;
    std::cout<<"Signal component, MPV: "<<land_gaus_bkg_fit->GetParameter(1)<<", S/N : "<<land_gaus_bkg_fit->GetParameter(1) / weight_average_GausWidth<<std::endl;

    cout<<Form("From 0 to 15, the integral: %.4f, From 0 to 210: %.4f", LandGausComp->Integral(0, 15), LandGausComp->Integral(0, 210))<<", ratio: "<<LandGausComp->Integral(0, 15)/LandGausComp->Integral(0, 210)<<endl;

    SoNRatio_Final(
        land_gaus_bkg_fit->GetParameter(1),
        land_gaus_bkg_fit->GetParError(1),
        BkgComp, // note : skinny
        BkgComp2,  // note : the fat one 

        input_directory
    );

    c1 -> cd();

    // std::cout<<"Gaus1 Integral error: "<< BkgComp->IntegralError(-60,60) << std::endl;
    // std::cout<<"Gaus2 Integral error: "<< BkgComp2->IntegralError(-60,60) << std::endl;

    double l0_fit_MPV = l0_fit -> GetParameter(1);
    double l2_fit_MPV = l2_fit -> GetParameter(1);

    double MPV_syst = std::max(
        std::abs(l0_fit_MPV - land_gaus_bkg_fit->GetParameter(1)),
        std::abs(l2_fit_MPV - land_gaus_bkg_fit->GetParameter(1))
    );

    std::cout<<"l0_fit_MPV: "<<l0_fit_MPV<<std::endl;
    std::cout<<"l2_fit_MPV: "<<l2_fit_MPV<<std::endl;

    std::cout<<Form("Landau MPV: %.2f #pm%.2f (stat.) #pm%.2f (syst.)",land_gaus_bkg_fit->GetParameter(1), land_gaus_bkg_fit->GetParError(1), MPV_syst)<<std::endl;

    TLatex * ltx_fit2 = new TLatex();
    ltx_fit2->SetNDC();
    ltx_fit2->SetTextSize(0.03);

    double ltx_fit_Y_offset = (ShowReducedChi2) ? 0 : 0.04;
    if (ShowReducedChi2){ltx_fit2->DrawLatex(0.4, 0.55, Form("#chi^{2} / NDF : %.2f / %d = %.2f",land_gaus_bkg_fit->GetChisquare(),land_gaus_bkg_fit->GetNDF(), land_gaus_bkg_fit->GetChisquare()/land_gaus_bkg_fit->GetNDF()));}
    ltx_fit2->DrawLatex(0.4, 0.51 + ltx_fit_Y_offset, Form("Landau MPV: %.2f #pm%.2f (stat.) #pm%.2f (syst.)",land_gaus_bkg_fit->GetParameter(1), land_gaus_bkg_fit->GetParError(1), MPV_syst));

    c1 -> Print(Form("%s/LbLVariation_Final_%.0f_%.0f.pdf", input_directory.c_str(), fit_range_pair.first, fit_range_pair.second));

    return 888;

}
