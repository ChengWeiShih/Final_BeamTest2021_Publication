#include "sPhenixStyle.C"

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

TGraphAsymmErrors * TH1D_converted_to_TGraphAsymmErrors(TH1D * hist, double NoYError = true, bool skipYZero = true) {

    TGraphAsymmErrors * graph = new TGraphAsymmErrors();

    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        double x = hist->GetBinCenter(i);
        double y = hist->GetBinContent(i);
        double x_err = hist->GetBinWidth(i) / 2;
        double y_err = (NoYError) ? 0 : hist->GetBinError(i);

        if (skipYZero && fabs(y) < 0.00001) {
            continue;
        }

        graph->SetPoint(graph->GetN(), x, y);
        graph->SetPointError(graph->GetN() - 1, x_err, x_err, y_err, y_err);

    }

    return graph;
}

int combine() {

    string file_directory_1 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out";
    string file_name_1 = "DAC_scan_out.root";

    string file_directory_2 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out_Scan2Scan3Match1Bin";
    string file_name_2 = "DAC_scan_out.root";

    string output_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/combine";

    double Yaxis_max = 2200;
    std::pair<double,double> fit_range_pair = {8, 140};
	std::pair<double,double> signal_purity_range = {40, 210};
    bool ShowReducedChi2 = false;

    TFile * file_in_1 = TFile::Open(Form("%s/%s", file_directory_1.c_str(), file_name_1.c_str() ) );
    TFile * file_in_2 = TFile::Open(Form("%s/%s", file_directory_2.c_str(), file_name_2.c_str() ) );

    std::vector<TH1D*> hist_vec_1; hist_vec_1.clear();
    std::vector<TH1D*> hist_vec_2; hist_vec_2.clear();

    std::vector<TGraphAsymmErrors*> TGAE_vec_1_syst; TGAE_vec_1_syst.clear();
    TGraphErrors * gr_syst_style = new TGraphErrors();

    std::vector<TGraphAsymmErrors*> TGAE_vec_ErrorComb; TGAE_vec_ErrorComb.clear();

    SetsPhenixStyle();

    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleAlign(22);
    gStyle->SetTitleX(0.5);
	gStyle->SetTitleY(0.95);
	gStyle->SetTitleFillColor(10);
	gStyle->SetPadTopMargin(0.09);
	gStyle->SetPadRightMargin(0.09);
    gStyle->SetPadTopMargin(0.05);

    TLegend * leg = new TLegend(0.4,0.72,0.8,0.84);
	leg -> SetBorderSize(0);
	leg -> SetMargin(0.15);
	leg -> SetTextSize(0.025);

    TLegend * leg_syst = new TLegend(0.4,0.78,0.8,0.84);
	leg_syst -> SetBorderSize(0);
	leg_syst -> SetMargin(0.15);
	leg_syst -> SetTextSize(0.03);

	TLegend * leg_final = new TLegend(0.55,0.72,0.8,0.84);
	leg_final -> SetBorderSize(0);
	leg_final -> SetMargin(0.2);
	leg_final -> SetTextSize(0.025);

	TLatex * ltx_matrix = new TLatex();
    ltx_matrix->SetNDC();
    ltx_matrix->SetTextSize(0.07);
    // ltx_matrix->SetTextAlign(31);

	TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    // ltx->SetTextAlign(31);

	TLatex * ltx_fit = new TLatex();
    ltx_fit->SetNDC();
    ltx_fit->SetTextSize(0.03);

    TCanvas * c1 = new TCanvas("c1", "c1", 950, 800);
    // c1 -> Print(Form("%s/comp.pdf(", output_directory.c_str()));

    TCanvas * c2 = new TCanvas("c2", "c2", 950, 800);
    // c2 -> Print(Form("%s/comp_syst.pdf(", output_directory.c_str()));

    for (int i = 0; i < 3; i++)
    {
        hist_vec_1.push_back((TH1D*)file_in_1->Get(Form("DAC_hist_all_%d", i)));
        hist_vec_2.push_back((TH1D*)file_in_2->Get(Form("DAC_hist_all_%d", i)));

        TGAE_vec_1_syst.push_back(TH1D_converted_to_TGraphAsymmErrors(hist_vec_1.back(), true));
        for (int point_i = 0; point_i < TGAE_vec_1_syst.back()->GetN(); point_i++){
            int the_bin_id = hist_vec_2.back()->FindBin(TGAE_vec_1_syst.back()->GetPointX(point_i));
            double diff_in_y = TGAE_vec_1_syst.back()->GetPointY(point_i) - hist_vec_2.back()->GetBinContent(the_bin_id);

            // std::cout<<"point_i : "<<point_i<<", Xpos: "<<TGAE_vec_1_syst.back()->GetPointX(point_i)<<", Ypos: "<<TGAE_vec_1_syst.back()->GetPointY(point_i)
            // <<", the_bin_id : "<<the_bin_id<<", Xpos: "<<hist_vec_2.back()->GetBinCenter(the_bin_id)<<", Ypos: "<<hist_vec_2.back()->GetBinContent(the_bin_id)
            // <<", diff_in_y : "<<diff_in_y<<std::endl;

            if (diff_in_y > 0){
                TGAE_vec_1_syst.back()->SetPointEYlow(point_i, fabs(diff_in_y));
            }
            else if (diff_in_y < 0){
                TGAE_vec_1_syst.back()->SetPointEYhigh(point_i, fabs(diff_in_y));
            }
        }
        TGAE_vec_1_syst.back() -> SetMarkerStyle(20);
        TGAE_vec_1_syst.back() -> SetMarkerSize(0);
        TGAE_vec_1_syst.back() -> SetMarkerColorAlpha(1,0);
        TGAE_vec_1_syst.back() -> SetFillColorAlpha(1,0.5);
        // TGAE_vec_1_syst.back() -> SetLineColorAlpha(1,0);
        // TGAE_vec_1_syst.back() -> SetLineWidth(1);
        TGAE_vec_1_syst.back()->SetFillColorAlpha(92, 0.85);

        TGAE_vec_ErrorComb.push_back(
            (TGraphAsymmErrors*) TGAE_vec_1_syst.back()->Clone( Form("TGAE_vec_ErrorComb_l%d", i) )
        );

        hist_vec_1.back()->SetLineColor(kBlack);
        hist_vec_1.back()->SetMarkerSize(0.8);
        hist_vec_1.back()->SetMarkerColor(kBlack);
        hist_vec_1.back()->SetMaximum(Yaxis_max);

        hist_vec_2.back()->SetLineColor(kRed);
        hist_vec_2.back()->SetMarkerSize(0.8);
        hist_vec_2.back()->SetMarkerColor(kRed);
        hist_vec_2.back()->SetMaximum(Yaxis_max);

        if (i == 0)
        {
            leg->AddEntry(hist_vec_1.back(), "Baseline (two bins for dist. matching)", "ep");
            leg->AddEntry(hist_vec_2.back(), "Variation (Scan2 matched by only the 6th bin)", "ep");

            gr_syst_style->SetPoint(0,0,0);
            gr_syst_style->SetMarkerStyle(hist_vec_1.back()->GetMarkerStyle());
            gr_syst_style->SetMarkerSize(hist_vec_1.back()->GetMarkerSize());
            gr_syst_style->SetMarkerColor(hist_vec_1.back()->GetMarkerColor());
            gr_syst_style->SetLineWidth(0);
            gr_syst_style->SetLineColorAlpha(1,0);
            gr_syst_style->SetFillColorAlpha(TGAE_vec_1_syst.back()->GetFillColor(), 0.92);


            leg_syst->AddEntry(gr_syst_style, "Data", "epf");
        }

        c1 -> cd();

        hist_vec_2.back()->Draw("ep");
        hist_vec_1.back()->Draw("ep same");

        ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021")); 
        leg ->Draw("same");

        c1 -> Print(Form("%s/l%d_comp.pdf", output_directory.c_str(),i));

        // Division : ------------------------------------------------------------------------------

        c2 -> cd();
        hist_vec_1.back()->Draw("ep");
        TGAE_vec_1_syst.back()->Draw("P2 same");
        hist_vec_1.back()->Draw("ep same");

        leg_syst->Draw("same");
        ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021")); 
        c2 ->Print(Form("%s/l%d_comp_syst.pdf", output_directory.c_str(),i));
        
    } 

    // c1 -> Print(Form("%s/comp.pdf)", output_directory.c_str()));
    // c2 -> Print(Form("%s/comp_syst.pdf)", output_directory.c_str()));

    c1 -> Clear();



    // Division : ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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

	TF1 * land_gaus_bkg_fit = new TF1("land_gaus_bkg_fit", Func_Langaus_bkg, 0, 200, 9);
	land_gaus_bkg_fit -> SetLineColor(TColor::GetColor("#F5321D"));
	land_gaus_bkg_fit -> SetLineWidth(2);
	land_gaus_bkg_fit -> SetParameters(
		8, 73.8,  26407, 12,
		1000,-0.5, -5,
		10, 20
	);

	land_gaus_bkg_fit -> SetParLimits(7,0,100000);
	land_gaus_bkg_fit -> SetParLimits(8,0,20);

	leg_final->AddEntry(land_gaus_bkg_fit, "Total Fit", "l");
	leg_final->AddEntry(LandGausComp, "Signal (Landau #otimes Gaussian)", "l");
	leg_final->AddEntry(BkgComp, "Background (Exponential)", "l");
	leg_final->AddEntry(BkgComp2, "Background (Gaussian)", "l");



    // c1 -> Print(Form("%s/comp_final.pdf(", output_directory.c_str()));
    for (int i = 0; i < 3; i++)
    {

        for (int point_i = 0; point_i < TGAE_vec_ErrorComb[i]->GetN(); point_i++){
            int the_bin_id = hist_vec_1[i]->FindBin(TGAE_vec_ErrorComb[i]->GetPointX(point_i));

            double stat_error = hist_vec_1[i]->GetBinError(the_bin_id);       
            double syst_error_high = TGAE_vec_ErrorComb[i]->GetErrorYhigh(point_i);
            double syst_error_low = TGAE_vec_ErrorComb[i]->GetErrorYlow(point_i);

            std::cout<<"layer: "<<i<<", point_i : "<<point_i<<", Xpos: "<<TGAE_vec_ErrorComb[i]->GetPointX(point_i)<<", Ypos: "<<TGAE_vec_ErrorComb[i]->GetPointY(point_i)<<", syst_error_high : "<<syst_error_high<<", syst_error_low : "<<syst_error_low<<std::endl;
            std::cout<<"layer: "<<i<<", bin_id : "<<the_bin_id<<", Xpos: "<<hist_vec_1[i]->GetBinCenter(the_bin_id)<<", Ypos: "<<hist_vec_1[i]->GetBinContent(the_bin_id)<<", stat_error : "<<stat_error<<std::endl;
            std::cout<<std::endl;

            double comb_error_high = sqrt(stat_error*stat_error + syst_error_high*syst_error_high);
            double comb_error_low = sqrt(stat_error*stat_error + syst_error_low*syst_error_low);
            
            TGAE_vec_ErrorComb[i]->SetPointEYhigh(point_i, comb_error_high);
            TGAE_vec_ErrorComb[i]->SetPointEYlow(point_i,  comb_error_low);
        }

        TGAE_vec_ErrorComb[i] -> SetMarkerStyle(20);
        TGAE_vec_ErrorComb[i] -> SetMarkerSize(0.8);
        TGAE_vec_ErrorComb[i] -> SetMarkerColor(1);
        TGAE_vec_ErrorComb[i] -> SetLineColor(1);
        TGAE_vec_ErrorComb[i] -> SetLineWidth(2);
        TGAE_vec_ErrorComb[i] -> GetYaxis() -> SetRangeUser(0, Yaxis_max);
        TGAE_vec_ErrorComb[i] -> GetXaxis() -> SetLimits(
            hist_vec_1[0]->GetXaxis()->GetXmin(),
            hist_vec_1[0]->GetXaxis()->GetXmax()
        );
        TGAE_vec_ErrorComb[i] -> GetXaxis() -> SetTitle(
            hist_vec_1[0]->GetXaxis()->GetTitle()
        );
        TGAE_vec_ErrorComb[i] -> GetYaxis() -> SetTitle(
            hist_vec_1[0]->GetYaxis()->GetTitle()
        );

        TGAE_vec_ErrorComb[i] -> Fit(land_gaus_bkg_fit, "", "", fit_range_pair.first, fit_range_pair.second);

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

        
        TGAE_vec_ErrorComb[i] -> Draw("ap");

        LandGausComp->Draw("lsame");
		BkgComp->Draw("lsame");
		BkgComp2->Draw("lsame");
        // land_gaus_bkg_fit->Draw("lsame");


        double signal_purity = LandGausComp->Integral(signal_purity_range.first, signal_purity_range.second)/land_gaus_bkg_fit->Integral(signal_purity_range.first, signal_purity_range.second);
		cout<<"Layer "<<i<<endl;
		cout<<"ChiSquare: "<<land_gaus_bkg_fit->GetChisquare()<<", NDF: "<<land_gaus_bkg_fit->GetNDF()<<", ReducedChiSquare: "<<land_gaus_bkg_fit->GetChisquare()/land_gaus_bkg_fit->GetNDF()<<endl;
		
		cout<<Form("Full Integral {%.0f,%.0f}: ",signal_purity_range.first, signal_purity_range.second)<<land_gaus_bkg_fit->Integral(signal_purity_range.first, signal_purity_range.second)
		<<Form(", LandGausComp Integral {%.0f, %.0f}: ",signal_purity_range.first,signal_purity_range.second)<<LandGausComp->Integral(signal_purity_range.first, signal_purity_range.second)
		<<", Signal purity: "<<signal_purity<<endl;
		
		cout<<endl;

		double ltx_fit_Y_offset = (ShowReducedChi2) ? 0 : 0.04;
		if (ShowReducedChi2){ltx_fit->DrawLatex(0.45, 0.65, Form("#chi^{2} / NDF : %.2f / %d = %.2f",land_gaus_bkg_fit->GetChisquare(),land_gaus_bkg_fit->GetNDF(), land_gaus_bkg_fit->GetChisquare()/land_gaus_bkg_fit->GetNDF()));}
		ltx_fit->DrawLatex(0.45, 0.61 + ltx_fit_Y_offset, Form("Landau MPV: %.2f #pm %.2f",land_gaus_bkg_fit->GetParameter(1), land_gaus_bkg_fit->GetParError(1)));
		ltx_fit->DrawLatex(0.45, 0.57 + ltx_fit_Y_offset, Form("Signal Purity (DAC #geq %.0f): %.2f%",signal_purity_range.first, signal_purity * 100.));

        leg_final->Draw("same");

        ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021")); 

        c1 -> Print(Form("%s/l%d_comp_final.pdf", output_directory.c_str(),i));
    }

    // c1 -> Print(Form("%s/comp_final.pdf)", output_directory.c_str()));

    return 888;
}