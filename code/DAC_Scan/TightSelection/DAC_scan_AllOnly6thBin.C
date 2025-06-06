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

	double gaus_bkg_eq = par[4]*TMath::Gaus(x[0],0,par[5]);

	return signal + gaus_bkg_eq;
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

void DAC_scan_AllOnly6thBin()
{
	std::pair<int,int> column_range_pair = {9, 11}; // note : use 1 to 13 coordination, column: {9, 10, 11}
	std::pair<double,double> fit_range_pair = {8, 116};
	std::pair<double,double> signal_purity_range = {15, 210};
	int number_of_runs = 5;
	bool DSE_evt_remove = true;
	double Yaxis_max = 2200;
	double sig_peak_content = 1000;
	bool SetFirstBinToZero = false;
	bool rough_event_selection = true;
	bool ShowReducedChi2 = true;

	TString folder_direction = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/TightSelection";
	TString output_directory = folder_direction + "/DACScan_out_AllOnly6thBin";
	system(Form("mkdir -p %s",output_directory.Data()));

	SetsPhenixStyle();

	gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleAlign(22);
    gStyle->SetTitleX(0.5);
	gStyle->SetTitleY(0.95);
	gStyle->SetTitleFillColor(10);
	gStyle->SetPadTopMargin(0.09);
	gStyle->SetPadRightMargin(0.09);
	
	vector<int> Scan_order_vec = {
		75, 76, 77, 78, 
		71, 72, 73, 74
	};

	// vector<TString>file_order; file_order.clear();
	
	// file_order.push_back("run75_no_clone_filter_all_clusters");//Run 75
	// file_order.push_back("run76_no_clone_filter_all_clusters");//Run 76
	// file_order.push_back("run77_no_clone_filter_all_clusters");//Run 77
	// file_order.push_back("run78_no_clone_filter_all_clusters");//Run 78  
	// file_order.push_back("run71_no_clone_filter_all_clusters");//Run 71
	// file_order.push_back("run72_no_clone_filter_all_clusters");//Run 72
	// file_order.push_back("run73_no_clone_filter_all_clusters");//Run 73
	// file_order.push_back("run74_no_clone_filter_all_clusters");//Run 74
	
	// file_order.push_back("BeamData_20211210-1526_0_filter_all");//Run 75
	// file_order.push_back("BeamData_20211210-1534_0_filter_all");//Run 76
	// file_order.push_back("BeamData_20211210-1541_0_filter_all");//Run 77
	// file_order.push_back("BeamData_20211210-1549_0_filter_all");//Run 78  
	// file_order.push_back("BeamData_20211210-1452_0_filter_all");//Run 71
	// file_order.push_back("BeamData_20211210-1501_0_filter_all");//Run 72
	// file_order.push_back("BeamData_20211210-1510_0_filter_all");//Run 73
	// file_order.push_back("BeamData_20211210-1519_0_filter_all");//Run 74

	
	TString color_code[8]={"#343434","#1A3947","#566575","#797983","#EFBD9D","#FCA26E","#F5751D","#F5321D"};
	TString color_code_2[8]={"#CC768D","#19768D","#DDA573","#009193","#6E9193","#941100","#A08144","#517E66"};
	// TH1D * DAC_hist_adc[3][8];
	TH1D * DAC_hist_bin[3][8];
	TH1D * DAC_hist_combine[3][8];
	TH2D * DAC_hist_adc_pos[3][8];
	TH2D * DAC_hist_adc_pos_raw[3][8];
	TH1D * DAC_hist_all[3];
	for (int i1=0; i1<3; i1++)
	{
		for (int i=0; i<8; i++)
		{
			// DAC_hist_adc[i1][i] = new TH1D ("",Form("DAC Scan l%i, %i",i1,i),8,8+20*i,40+20*i);//is not used.
			DAC_hist_bin[i1][i] = new TH1D("",Form("Scan-%d, Layer %d;HitAdc;Entries",i,i1),8,0,8);//consider bin (adc) only
			DAC_hist_bin[i1][i]->SetLineColor( TColor::GetColor(Form("%s",color_code_2[7].Data())) );
			DAC_hist_bin[i1][i]->SetLineWidth(2);

			// DAC_hist_combine[i1][i] = new TH1D ("",Form("DAC Scan combine l%i %i",i1,i),50,0,200);//include the real adc setting
			DAC_hist_combine[i1][i] = new TH1D("",Form(""),50,0,200);//include the real adc setting 
			DAC_hist_combine[i1][i]->SetLineColor(TColor::GetColor(Form("%s",color_code_2[i].Data())));
			DAC_hist_combine[i1][i]->GetYaxis()->SetRangeUser(0,Yaxis_max);
			DAC_hist_combine[i1][i]->GetYaxis()->SetTitle("Entries (A.U.)");
			DAC_hist_combine[i1][i]->SetLineWidth(2);
			DAC_hist_combine[i1][i]->SetMarkerSize(0.3);

			DAC_hist_combine[i1][i]->GetXaxis()->SetTitle("DAC value");

			DAC_hist_adc_pos[i1][i] = new TH2D("",Form("Scan-%d, Layer %d",i,i1),8,0,8,128,-9.984,9.984);
			DAC_hist_adc_pos_raw[i1][i] = new TH2D("",Form("(Raw) Scan-%d, Layer %d",i,i1),8,0,8,128,-9.984,9.984);

		}	

		DAC_hist_all[i1] = new TH1D ("",Form(""),50,0,200);//include the real adc setting, combine all hist into one hist
		DAC_hist_all[i1] -> SetMarkerSize(0.8);
		DAC_hist_all[i1] -> SetLineColor( TColor::GetColor(Form("#1A3947")) );
		DAC_hist_all[i1] -> GetYaxis()->SetRangeUser(0,Yaxis_max);
		DAC_hist_all[i1] -> GetYaxis()->SetTitle("Entries (A.U.)");
		DAC_hist_all[i1] -> GetXaxis()->SetTitle("DAC value");
	}
	

	

	TLegend * leg = new TLegend(0.62,0.51,0.8,0.84);
	leg -> SetBorderSize(0);
	leg -> SetMargin(0.15);
	leg -> SetTextSize(0.02);

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

	TCanvas * c1 = new TCanvas ("c1","c1",2000,1800); 
	c1->Divide(3,3);


	int adc_setting_run[8][9]=
	{	
		{8   ,12 ,16 ,20 ,24 ,28 ,32 ,36 ,40},
		{28  ,32 ,36 ,40 ,44 ,48 ,52 ,56 ,60},
		{48  ,52 ,56 ,60 ,64 ,68 ,72 ,76 ,80},
		{68  ,72 ,76 ,80 ,84 ,88 ,92 ,96 ,100},
		{88  ,92 ,96 ,100,104,108,112,116,120},
		{108 ,112,116,120,124,128,132,136,140},
		{128 ,132,136,140,144,148,152,156,160},
		{148 ,152,156,160,164,168,172,176,180}
	};

	std::map<int,int> nominal_setting_map = {
		{15, 0},
		{30, 1},
		{60, 2},
		{90, 3},
		{120, 4},
		{150, 5},
		{180, 6},
		{210, 7}
	};

	TFile *f1;
	TFile *f2;
	TFile *f3;
	TTree *tree_both_in;


	for (int i=0;i<8;i++)
	{
		f1 = TFile::Open(Form("%s/Data_DACScanTight_Run%d_Column9_BaseLine.root", folder_direction.Data(), Scan_order_vec[i]));
		f2 = TFile::Open(Form("%s/Data_DACScanTight_Run%d_Column10_BaseLine.root", folder_direction.Data(), Scan_order_vec[i]));
		f3 = TFile::Open(Form("%s/Data_DACScanTight_Run%d_Column11_BaseLine.root", folder_direction.Data(), Scan_order_vec[i]));

		TH1D * temp_L0_dist_f1 = (TH1D*)f1->Get("h1D_GoodTrack_1HitClusAdc_L0");
		TH1D * temp_L1_dist_f1 = (TH1D*)f1->Get("h1D_GoodL1Hit_1HitClusAdc_L1");
		TH1D * temp_L2_dist_f1 = (TH1D*)f1->Get("h1D_GoodTrack_1HitClusAdc_L2");

		TH1D * temp_L0_dist_f2 = (TH1D*)f2->Get("h1D_GoodTrack_1HitClusAdc_L0");
		TH1D * temp_L1_dist_f2 = (TH1D*)f2->Get("h1D_GoodL1Hit_1HitClusAdc_L1");
		TH1D * temp_L2_dist_f2 = (TH1D*)f2->Get("h1D_GoodTrack_1HitClusAdc_L2");

		TH1D * temp_L0_dist_f3 = (TH1D*)f3->Get("h1D_GoodTrack_1HitClusAdc_L0");
		TH1D * temp_L1_dist_f3 = (TH1D*)f3->Get("h1D_GoodL1Hit_1HitClusAdc_L1");
		TH1D * temp_L2_dist_f3 = (TH1D*)f3->Get("h1D_GoodTrack_1HitClusAdc_L2");


		TH2D * temp_L0_dist_ADCPos_f1 = (TH2D*)f1->Get("h2D_GoodTrack_1HitClusAdc_Pos_L0");
		TH2D * temp_L1_dist_ADCPos_f1 = (TH2D*)f1->Get("h2D_GoodL1Hit_1HitClusAdc_Pos_L1");
		TH2D * temp_L2_dist_ADCPos_f1 = (TH2D*)f1->Get("h2D_GoodTrack_1HitClusAdc_Pos_L2");

		TH2D * temp_L0_dist_ADCPos_f2 = (TH2D*)f2->Get("h2D_GoodTrack_1HitClusAdc_Pos_L0");
		TH2D * temp_L1_dist_ADCPos_f2 = (TH2D*)f2->Get("h2D_GoodL1Hit_1HitClusAdc_Pos_L1");
		TH2D * temp_L2_dist_ADCPos_f2 = (TH2D*)f2->Get("h2D_GoodTrack_1HitClusAdc_Pos_L2");

		TH2D * temp_L0_dist_ADCPos_f3 = (TH2D*)f3->Get("h2D_GoodTrack_1HitClusAdc_Pos_L0");
		TH2D * temp_L1_dist_ADCPos_f3 = (TH2D*)f3->Get("h2D_GoodL1Hit_1HitClusAdc_Pos_L1");
		TH2D * temp_L2_dist_ADCPos_f3 = (TH2D*)f3->Get("h2D_GoodTrack_1HitClusAdc_Pos_L2");

		
		TH2D * temp_L0_dist_ADCPos_Raw_f1 = (TH2D*)f1->Get("h2D_raw_1HitClusAdc_Pos_L0");
		TH2D * temp_L1_dist_ADCPos_Raw_f1 = (TH2D*)f1->Get("h2D_raw_1HitClusAdc_Pos_L1");
		TH2D * temp_L2_dist_ADCPos_Raw_f1 = (TH2D*)f1->Get("h2D_raw_1HitClusAdc_Pos_L2");

		TH2D * temp_L0_dist_ADCPos_Raw_f2 = (TH2D*)f2->Get("h2D_raw_1HitClusAdc_Pos_L0");
		TH2D * temp_L1_dist_ADCPos_Raw_f2 = (TH2D*)f2->Get("h2D_raw_1HitClusAdc_Pos_L1");
		TH2D * temp_L2_dist_ADCPos_Raw_f2 = (TH2D*)f2->Get("h2D_raw_1HitClusAdc_Pos_L2");

		TH2D * temp_L0_dist_ADCPos_Raw_f3 = (TH2D*)f3->Get("h2D_raw_1HitClusAdc_Pos_L0");
		TH2D * temp_L1_dist_ADCPos_Raw_f3 = (TH2D*)f3->Get("h2D_raw_1HitClusAdc_Pos_L1");
		TH2D * temp_L2_dist_ADCPos_Raw_f3 = (TH2D*)f3->Get("h2D_raw_1HitClusAdc_Pos_L2");


		DAC_hist_bin[0][i] -> Add(temp_L0_dist_f1);
		DAC_hist_bin[0][i] -> Add(temp_L0_dist_f2);
		DAC_hist_bin[0][i] -> Add(temp_L0_dist_f3);

		// DAC_hist_bin[0][i] -> Add(temp_L1_dist_f1);
		// DAC_hist_bin[0][i] -> Add(temp_L1_dist_f2);
		// DAC_hist_bin[0][i] -> Add(temp_L1_dist_f3);

		// DAC_hist_bin[0][i] -> Add(temp_L2_dist_f1);
		// DAC_hist_bin[0][i] -> Add(temp_L2_dist_f2);
		// DAC_hist_bin[0][i] -> Add(temp_L2_dist_f3);


		DAC_hist_bin[1][i] -> Add(temp_L1_dist_f1);
		DAC_hist_bin[1][i] -> Add(temp_L1_dist_f2);
		DAC_hist_bin[1][i] -> Add(temp_L1_dist_f3);

		DAC_hist_bin[2][i] -> Add(temp_L2_dist_f1);
		DAC_hist_bin[2][i] -> Add(temp_L2_dist_f2);
		DAC_hist_bin[2][i] -> Add(temp_L2_dist_f3);



		DAC_hist_adc_pos[0][i] -> Add(temp_L0_dist_ADCPos_f1);
		DAC_hist_adc_pos[0][i] -> Add(temp_L0_dist_ADCPos_f2);
		DAC_hist_adc_pos[0][i] -> Add(temp_L0_dist_ADCPos_f3);

		DAC_hist_adc_pos[1][i] -> Add(temp_L1_dist_ADCPos_f1);
		DAC_hist_adc_pos[1][i] -> Add(temp_L1_dist_ADCPos_f2);
		DAC_hist_adc_pos[1][i] -> Add(temp_L1_dist_ADCPos_f3);

		DAC_hist_adc_pos[2][i] -> Add(temp_L2_dist_ADCPos_f1);
		DAC_hist_adc_pos[2][i] -> Add(temp_L2_dist_ADCPos_f2);
		DAC_hist_adc_pos[2][i] -> Add(temp_L2_dist_ADCPos_f3);

		DAC_hist_adc_pos[0][i] -> GetXaxis() -> SetTitle( temp_L0_dist_ADCPos_f1 -> GetXaxis() -> GetTitle() );
		DAC_hist_adc_pos[0][i] -> GetYaxis() -> SetTitle( temp_L0_dist_ADCPos_f1 -> GetYaxis() -> GetTitle() );

		DAC_hist_adc_pos[1][i] -> GetXaxis() -> SetTitle( temp_L1_dist_ADCPos_f1 -> GetXaxis() -> GetTitle() );
		DAC_hist_adc_pos[1][i] -> GetYaxis() -> SetTitle( temp_L1_dist_ADCPos_f1 -> GetYaxis() -> GetTitle() );

		DAC_hist_adc_pos[2][i] -> GetXaxis() -> SetTitle( temp_L2_dist_ADCPos_f1 -> GetXaxis() -> GetTitle() );
		DAC_hist_adc_pos[2][i] -> GetYaxis() -> SetTitle( temp_L2_dist_ADCPos_f1 -> GetYaxis() -> GetTitle() );



		DAC_hist_adc_pos_raw[0][i] -> Add(temp_L0_dist_ADCPos_Raw_f1);
		DAC_hist_adc_pos_raw[0][i] -> Add(temp_L0_dist_ADCPos_Raw_f2);
		DAC_hist_adc_pos_raw[0][i] -> Add(temp_L0_dist_ADCPos_Raw_f3);

		DAC_hist_adc_pos_raw[1][i] -> Add(temp_L1_dist_ADCPos_Raw_f1);
		DAC_hist_adc_pos_raw[1][i] -> Add(temp_L1_dist_ADCPos_Raw_f2);
		DAC_hist_adc_pos_raw[1][i] -> Add(temp_L1_dist_ADCPos_Raw_f3);

		DAC_hist_adc_pos_raw[2][i] -> Add(temp_L2_dist_ADCPos_Raw_f1);
		DAC_hist_adc_pos_raw[2][i] -> Add(temp_L2_dist_ADCPos_Raw_f2);
		DAC_hist_adc_pos_raw[2][i] -> Add(temp_L2_dist_ADCPos_Raw_f3);

		DAC_hist_adc_pos_raw[0][i] -> GetXaxis() -> SetTitle( temp_L0_dist_ADCPos_Raw_f1 -> GetXaxis() -> GetTitle() );
		DAC_hist_adc_pos_raw[0][i] -> GetYaxis() -> SetTitle( temp_L0_dist_ADCPos_Raw_f1 -> GetYaxis() -> GetTitle() );

		DAC_hist_adc_pos_raw[1][i] -> GetXaxis() -> SetTitle( temp_L1_dist_ADCPos_Raw_f1 -> GetXaxis() -> GetTitle() );
		DAC_hist_adc_pos_raw[1][i] -> GetYaxis() -> SetTitle( temp_L1_dist_ADCPos_Raw_f1 -> GetYaxis() -> GetTitle() );

		DAC_hist_adc_pos_raw[2][i] -> GetXaxis() -> SetTitle( temp_L2_dist_ADCPos_Raw_f1 -> GetXaxis() -> GetTitle() );
		DAC_hist_adc_pos_raw[2][i] -> GetYaxis() -> SetTitle( temp_L2_dist_ADCPos_Raw_f1 -> GetYaxis() -> GetTitle() );



		for (int i2 = 0; i2 < 3; i2++){
			DAC_hist_bin[i2][i] -> SetBinContent(8,0); // note : set overflow bin to zero
			if (i < number_of_runs - 1) {DAC_hist_bin[i2][i] -> SetBinContent(7,0);} // note : set second overlapped bin to zero

			for (int bin_i = 1; bin_i <= 7; bin_i++){ // note : remove overflow bin
				double bin_content = DAC_hist_bin[i2][i]->GetBinContent(bin_i);
				
				DAC_hist_combine[i2][i]	-> SetBinContent(
					DAC_hist_combine[i2][i] -> FindBin(adc_setting_run[i][bin_i-1]),
					bin_content
				);
			}
		}

	}// end of i, the run

	if (SetFirstBinToZero){
		DAC_hist_bin[0][0]->SetBinContent(1,0);
		DAC_hist_bin[0][0]->SetBinError(1,0);
		
		DAC_hist_bin[1][0]->SetBinContent(1,0);
		DAC_hist_bin[1][0]->SetBinError(1,0);

		DAC_hist_bin[2][0]->SetBinContent(1,0);
		DAC_hist_bin[2][0]->SetBinError(1,0);

		DAC_hist_combine[0][0]->SetBinContent(3,0);
		DAC_hist_combine[0][0]->SetBinError(3,0);

		DAC_hist_combine[1][0]->SetBinContent(3,0);
		DAC_hist_combine[1][0]->SetBinError(3,0);

		DAC_hist_combine[2][0]->SetBinContent(3,0);
		DAC_hist_combine[2][0]->SetBinError(3,0);

	}

	// c1->Print( Form("%s/DAC_scan_matrix.pdf(",output_directory.Data()) );
	for (int i=0; i<3; i++)
	{	
		// c1->Divide(3,3);
		for (int i1=0; i1<8; i1++)
		{
			c1->cd(i1+1);
			DAC_hist_bin[i][i1]->SetStats(0);
			DAC_hist_bin[i][i1]->Draw("ep");	

		}
		c1->cd(9);
		ltx_matrix->DrawLatex(0.0, 0.5, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));

		c1->Print( Form("%s/l%d_DAC_scan_matrix.pdf",output_directory.Data(),i) );
		// c1->Clear();
	}
	// c1->Print( Form("%s/DAC_scan_matrix.pdf)",output_directory.Data()) );
	c1->Clear();

	gStyle->SetPadRightMargin(0.16);
	c1 = new TCanvas ("c1","c1",2000,1800);
	c1 -> Divide(3,3);
	for (int i = 0; i < 3; i++) {
		for (int i1 = 0; i1 < 8; i1++){
			c1->cd(i1+1);
			DAC_hist_adc_pos[i][i1]->SetStats(0);
			DAC_hist_adc_pos[i][i1]->Draw("colz0");
		}
		c1->cd(9);
		ltx_matrix->DrawLatex(0.0, 0.5, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));

		c1->Print( Form("%s/l%d_DACScan_AdcPos_matrix.pdf",output_directory.Data(), i) );
	}
	c1 -> Clear();

	c1 -> Divide(3,3);
	for (int i = 0; i < 3; i++) {
		for (int i1 = 0; i1 < 8; i1++){
			c1->cd(i1+1);
			DAC_hist_adc_pos_raw[i][i1]->SetStats(0);
			DAC_hist_adc_pos_raw[i][i1]->Draw("colz0");
		}
		c1->cd(9);
		ltx_matrix->DrawLatex(0.0, 0.5, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));

		c1->Print( Form("%s/l%d_DACScan_AdcPos_Raw_matrix.pdf",output_directory.Data(), i) );
	}


	// for (int i1=0; i1<3; i1++)//normalization
	// {	
	// 	for (int i=0; i<8; i++)
	// 	{
	// 		DAC_hist_combine[i1][i]->Scale(1./DAC_hist_combine[i1][i]->Integral(-1,-1)); 
	// 	}	
	// }
	double previous_content;
	double next_content;
	double scale_weight;

	gStyle->SetPadTopMargin(0.05);
	// gStyle->SetPadRightMargin(0.16);
	// gStyle->SetPadLeftMargin(0.16);

	gStyle->SetPadRightMargin(0.09);

	TCanvas * c2 = new TCanvas ("c2","c2",950, 800); 
	c2->cd();

	// for (int i1=0; i1<3; i1++)//combine all hist into one hist
	// {	
	// 	for (int i=0; i<8; i++)
	// 	{
	// 		DAC_hist_adc_pos[i1][i]->SetStats(0);
	// 		DAC_hist_adc_pos[i1][i]->Draw("colz0");
	// 		c2 -> Print( Form("%s/l%d_DACScan_AdcPos_Scan%d.pdf",output_directory.Data(), i1, i) );
	// 	}
	// }

	// gStyle->SetPadTopMargin(0.05);
	// gStyle->SetPadRightMargin(0.09);
	// gStyle->SetPadLeftMargin(0.16);
	// c2 = new TCanvas ("c2","c2",950, 800); 
	// c2 -> cd();


	// TPad *pad = new TPad(Form("pad1"), "", 0.0, 0.0, 1.0, 1.0);
	// //pad->SetTopMargin(0.12);
	// //pad->SetBottomMargin(0.120);
	// pad->SetLeftMargin(0.15);
	// pad->SetRightMargin(0.05);
	// //pad->SetGrid(1, 1);
	// pad->Draw("same");
	// pad->cd();

	// |_|_|_|_|_|_|_|_| first run, the 8th bin is negelected because it is the overflow bin. use the 6th and 7th bin to do the hist size match
	//           |_|_|_|_|_|_|_|_| second run
	//                     |_|_|_|_|_|_|_|_| third run, vice versa

	// c2->Print( Form("%s/DAC_scan_overlap.pdf(",output_directory.Data()) );
	
	for (int i = 0; i < number_of_runs; i++)
	{
		leg->AddEntry(DAC_hist_combine[0][i],Form("Scan %d, range: {%d-%d} DAC",i, adc_setting_run[i][0], adc_setting_run[i][7]),"epl");
	}

	int selected_scan = 3; // todo: change the selected scan
	for (int i1 = 0; i1 < 3; i1++)
	{
		double scan3_max_bin_content = DAC_hist_combine[i1][selected_scan] -> GetBinContent(DAC_hist_combine[i1][selected_scan]->GetMaximumBin()); 
		DAC_hist_combine[i1][selected_scan]->Scale(sig_peak_content/scan3_max_bin_content);

		for (int scan_i = selected_scan; scan_i > 0; scan_i--)
		{
			previous_content = DAC_hist_combine[i1][scan_i-1]->GetBinContent(3+scan_i*5);// note : adc5 and adc6 bin of the previous hist
			next_content     = DAC_hist_combine[i1][scan_i]  ->GetBinContent(3+scan_i*5);// note : adc0 and adc1 bin of the next hist  --> this 

			std::cout<<"layer: "<<i1<<", scan_i: "<<scan_i<<", previous_content: "<<previous_content<<", next_content: "<<next_content<<std::endl;

			if (next_content == 0) { scale_weight = 1; }
			else
			{
				if (previous_content==0) { previous_content = 1.; }
				scale_weight = next_content/previous_content;
			}

			DAC_hist_combine[i1][scan_i - 1] -> Scale(scale_weight);

			previous_content = 0;
			next_content = 0;
		}

		for (int scan_i = selected_scan; scan_i < number_of_runs-1; scan_i++)
		{
			previous_content = DAC_hist_combine[i1][scan_i]->GetBinContent(8+scan_i*5); //note : adc5 and adc6 bin of the previous hist ---> this
			next_content     = DAC_hist_combine[i1][scan_i+1]->GetBinContent(8+scan_i*5); //note : adc0 and adc1 bin of the next hist  

			std::cout<<"layer: "<<i1<<", scan_i: "<<scan_i<<", previous_content: "<<previous_content<<", next_content: "<<next_content<<std::endl;

			if (next_content == 0) { scale_weight = 1; }
			else
			{
				if (previous_content==0) { previous_content = 1.; }
				scale_weight = previous_content/next_content;
			}

			DAC_hist_combine[i1][scan_i + 1] -> Scale(scale_weight);

			previous_content = 0;
			next_content = 0;
		}


		for (int scan_i = 0; scan_i < number_of_runs; scan_i++)
		{
			DAC_hist_combine[i1][scan_i]->GetYaxis()->SetRangeUser(0,Yaxis_max);
			DAC_hist_combine[i1][scan_i]->SetStats(0);

			if (scan_i == 0){
				DAC_hist_combine[i1][scan_i]->Draw("hist");
				DAC_hist_combine[i1][scan_i]->Draw("ep same");
			}
			else
			{
				DAC_hist_combine[i1][scan_i]->Draw("hist same");	
				DAC_hist_combine[i1][scan_i]->Draw("ep same");	
			}
		}

		leg->Draw("same");
		ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021")); 

		c2->Print( Form("%s/l%d_DAC_scan_overlap.pdf",output_directory.Data(), i1) );		
		c2 -> Clear();
	}

	




	TF1 *LandGausComp = new TF1("LandGausComp", Func_Langaus, 0, 200,4);
	LandGausComp -> SetLineColor(61);
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

	TF1 * land_gaus_bkg_fit = new TF1("land_gaus_bkg_fit", Func_Langaus_bkg, 0, 200, 6);
	land_gaus_bkg_fit -> SetLineColor(TColor::GetColor("#F5321D"));
	land_gaus_bkg_fit -> SetLineWidth(2);
	land_gaus_bkg_fit -> SetParameters(
		8, 73.8,  26407, 12,
		// 1000,-0.5, -5,
		10, 20
	);

	land_gaus_bkg_fit -> SetParLimits(7,0,100000);
	land_gaus_bkg_fit -> SetParLimits(8,0,20);

	leg_final->AddEntry(land_gaus_bkg_fit, "Total Fit", "l");
	leg_final->AddEntry(LandGausComp, "Signal (Landau #otimes Gaussian)", "l");
	// leg_final->AddEntry(BkgComp, "Background (Exponential)", "l");
	leg_final->AddEntry(BkgComp2, "Background (Gaussian)", "l");

	// double width_v, width_e, mip_v, mip_e, area1_v, area1_e, gsig1_v, gsig1_e;
	// float  chi1, con1;
	// int    ndf1;
	
	TFile * DAC_scan_out = new TFile (Form("%s/DAC_scan_out.root", output_directory.Data()), "RECREATE");
	

	TLatex *tex11 = new TLatex();
	tex11 -> SetNDC();
	tex11 -> SetTextFont(42);
	tex11 -> SetTextSize(0.028);
	tex11 -> SetTextAlign(13);
	
	// merge all 8 hist into 1 hist
	for (int i1=0; i1<3; i1++)
	{
		DAC_hist_all[i1]->SetStats(0);
		for (int i=0; i<number_of_runs; i++)
		{	
			for (int i2=1; i2<51; i2++)
			{
				if (DAC_hist_combine[i1][i]->GetBinContent(i2)!=0)
				{
					if (DAC_hist_all[i1]->GetBinContent(i2)==0)
					{
						DAC_hist_all[i1]->SetBinContent(i2,DAC_hist_combine[i1][i]->GetBinContent(i2));
						DAC_hist_all[i1]->SetBinError(i2,DAC_hist_combine[i1][i]->GetBinError(i2));
					}
					else 
					{	
						double temp_original_content = DAC_hist_all[i1]->GetBinContent(i2);
						double temp_original_error   = DAC_hist_all[i1]->GetBinError(i2);
						double temp_original_weight  = 1./pow(temp_original_error,2);

						double temp_coming_content = DAC_hist_combine[i1][i]->GetBinContent(i2);
						double temp_coming_error   = DAC_hist_combine[i1][i]->GetBinError(i2);
						double temp_coming_weight  = 1./pow(temp_coming_error,2);

						double weight_average = (temp_original_content*temp_original_weight + temp_coming_content*temp_coming_weight)/(temp_original_weight + temp_coming_weight);
						double weight_average_error = sqrt(1./(temp_original_weight + temp_coming_weight));

						double ChiSquare = pow(temp_original_content - weight_average,2)/pow(temp_original_error,2) + pow(temp_coming_content - weight_average,2)/pow(temp_coming_error,2);
						double ndf = 1.;
						double ReducedChiSquare = ChiSquare/ndf;

						// double final_weight_average_error = (ReducedChiSquare > 1) ? sqrt(ReducedChiSquare)*weight_average_error : weight_average_error;
						double final_weight_average_error = (temp_original_error + temp_coming_error)/2.;

						std::cout<<"layer: "<<i1<<", scan_i: "<<i<<", bin: "<<i2<<", temp_original_content: "<<temp_original_content<<", temp_coming_content: "<<temp_coming_content<<", weight_average: "<<weight_average<<std::endl;
						std::cout<<"layer: "<<i1<<", scan_i: "<<i<<", bin: "<<i2<<", temp_original_error: "<<temp_original_error<<", temp_coming_error: "<<temp_coming_error<<", weight_average_error: "<<weight_average_error<<std::endl;
						std::cout<<"layer: "<<i1<<", scan_i: "<<i<<", bin: "<<i2<<", temp_original_weight: "<<temp_original_weight<<", temp_coming_weight: "<<temp_coming_weight<<std::endl;
						std::cout<<"layer: "<<i1<<", scan_i: "<<i<<", bin: "<<i2<<", ChiSquare: "<<ChiSquare<<", ndf: "<<ndf<<", ReducedChiSquare: "<<ReducedChiSquare<<", final_weight_average_error: "<<final_weight_average_error<<std::endl;
						std::cout<<std::endl;



						DAC_hist_all[i1]->SetBinContent(i2, weight_average);
						DAC_hist_all[i1]->SetBinError(i2, final_weight_average_error);
					}
				}
				
			}
			
		}
		DAC_hist_all[i1]->SetMarkerStyle(20);
		DAC_hist_all[i1]->SetMarkerSize(0.8);
		DAC_hist_all[i1]->SetMarkerColor( TColor::GetColor(Form("#1A3947")) );
		DAC_hist_all[i1]->Draw("ep");

		ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));
		
		DAC_scan_out -> cd();
		DAC_hist_all[i1] -> Write(Form("DAC_hist_all_%d",i1));

		// for (int i4 = 0; i4 < 50; i4++)
		// {
		// 	cout<<" layer : "<<i1<<" Bin : "<<i4+1<<" Center : "<<DAC_hist_all[i1]->GetBinCenter(i4+1)<<" Width : "<<DAC_hist_all[i1]->GetBinWidth(i4+1)<<" Entry : "<<DAC_hist_all[i1]->GetBinContent(i4+1)<<endl;
		// }

		if (i1==0)
		{
			// DAC_hist_all[i1] -> Fit(func_fit, "L0", "", 60, 140);
			// DAC_hist_all[i1] -> Fit(Land_Gaus_fit, "L0", "", 60, 140);

			DAC_hist_all[i1] -> Fit(land_gaus_bkg_fit, "", "", fit_range_pair.first, fit_range_pair.second);
		}
		else 
		{
			// func_fit -> SetParameters(4.5, 73.8,  26407, 12);
			// DAC_hist_all[i1] -> Fit(func_fit, "L0", "", 40, 130);	
			// DAC_hist_all[i1] -> Fit(Land_Gaus_fit, "L0", "", 40, 130);

			DAC_hist_all[i1] -> Fit(land_gaus_bkg_fit, "", "", fit_range_pair.first, fit_range_pair.second);
		}

		LandGausComp->SetParameters(
			land_gaus_bkg_fit->GetParameter(0),
			land_gaus_bkg_fit->GetParameter(1),
			land_gaus_bkg_fit->GetParameter(2),
			land_gaus_bkg_fit->GetParameter(3)
		);

		// BkgComp->SetParameters(
		// 	land_gaus_bkg_fit->GetParameter(4),
		// 	land_gaus_bkg_fit->GetParameter(5),
		// 	land_gaus_bkg_fit->GetParameter(6)
		// );

		BkgComp2->SetParameters(
			land_gaus_bkg_fit->GetParameter(4),
			0,
			land_gaus_bkg_fit->GetParameter(5)
			// land_gaus_bkg_fit->GetParameter(9)
		);

		std::cout<<"Gaussian noise, height: "<<land_gaus_bkg_fit->GetParameter(4)<<", width: "<<land_gaus_bkg_fit->GetParameter(5)<<std::endl;
		std::cout<<"Signal component, MPV: "<<land_gaus_bkg_fit->GetParameter(1)<<", S/N : "<<land_gaus_bkg_fit->GetParameter(1) / land_gaus_bkg_fit->GetParameter(5)<<std::endl;

		LandGausComp->Draw("lsame");
		BkgComp->Draw("lsame");
		BkgComp2->Draw("lsame");

		DAC_hist_all[i1]->Draw("ep same");

		leg_final->Draw("same");


		double signal_purity = LandGausComp->Integral(signal_purity_range.first, signal_purity_range.second)/land_gaus_bkg_fit->Integral(signal_purity_range.first, signal_purity_range.second);
		cout<<"Layer "<<i1<<endl;
		cout<<"ChiSquare: "<<land_gaus_bkg_fit->GetChisquare()<<", NDF: "<<land_gaus_bkg_fit->GetNDF()<<", ReducedChiSquare: "<<land_gaus_bkg_fit->GetChisquare()/land_gaus_bkg_fit->GetNDF()<<endl;
		
		cout<<Form("Full Integral {%.0f,%.0f}: ",signal_purity_range.first, signal_purity_range.second)<<land_gaus_bkg_fit->Integral(signal_purity_range.first, signal_purity_range.second)
		<<Form(", LandGausComp Integral {%.0f, %.0f}: ",signal_purity_range.first,signal_purity_range.second)<<LandGausComp->Integral(signal_purity_range.first, signal_purity_range.second)
		<<", Signal purity: "<<signal_purity<<endl;
		
		cout<<endl;

		double ltx_fit_Y_offset = (ShowReducedChi2) ? 0 : 0.04;
		if (ShowReducedChi2){ltx_fit->DrawLatex(0.45, 0.65, Form("#chi^{2} / NDF : %.2f / %d = %.2f",land_gaus_bkg_fit->GetChisquare(),land_gaus_bkg_fit->GetNDF(), land_gaus_bkg_fit->GetChisquare()/land_gaus_bkg_fit->GetNDF()));}
		ltx_fit->DrawLatex(0.45, 0.61 + ltx_fit_Y_offset, Form("Landau MPV: %.2f #pm %.2f",land_gaus_bkg_fit->GetParameter(1), land_gaus_bkg_fit->GetParError(1)));
		ltx_fit->DrawLatex(0.45, 0.57 + ltx_fit_Y_offset, Form("Signal Purity (DAC #geq %.0f): %.2f%",signal_purity_range.first, signal_purity * 100.));
		

		

		// * Extract
		// * - Extract signal width
		// width_v = func_fit -> GetParameter(0);  
		// width_e = func_fit -> GetParError(0);   
		// // * - Extract signal MPV
		// mip_v = func_fit -> GetParameter(1);    
		// mip_e = func_fit -> GetParError(1);     
		// // * - Extract signal area
		// area1_v = func_fit -> GetParameter(2);   
		// area1_e = func_fit -> GetParError(2);    
		// // * - Extract signal Gaussian sigma
		// gsig1_v = func_fit -> GetParameter(3);   
		// gsig1_e = func_fit -> GetParError(3);    
		// // * - Extract chi2/ndf
		// chi1 = func_fit -> GetChisquare();
		// ndf1 = func_fit -> GetNDF();

		// Land_Gaus_fit -> Draw("lsame");
		// func_fit->Draw("lsame");
		// land_gaus_bkg_fit->Draw("lsame");

		// tex11 -> DrawLatex(0.2, 0.800+0.05, Form("#bf{Hoa method}"));
		// tex11 -> DrawLatex(0.2, 0.750+0.05, Form("#bullet Width = %.3f  #pm %.3f", width_v , width_e));
		// tex11 -> DrawLatex(0.2, 0.700+0.05, Form("#bullet MIP = %.3f  #pm %.3f", mip_v , mip_e));
		// tex11 -> DrawLatex(0.2, 0.650+0.05, Form("#bullet Area = %.2f  #pm %.2f", area1_v , area1_e));
		// tex11 -> DrawLatex(0.2, 0.600+0.05, Form("#bullet sigma = %.3f  #pm %.3f", gsig1_v, gsig1_e));
		// tex11 -> DrawLatex(0.2, 0.550+0.05, Form("#bullet #chi^{2}/ndf = %.3f/%d = %.3f", chi1, ndf1, chi1/double(ndf1)));


		// tex11 -> DrawLatex(0.58, 0.800+0.05, Form("#bf{root con method}"));
		// tex11 -> DrawLatex(0.58, 0.750+0.05, Form("#bullet Size = %.2f  #pm %.2f", Land_Gaus_fit->GetParameter(0), Land_Gaus_fit->GetParError(0)));
		// tex11 -> DrawLatex(0.58, 0.700+0.05, Form("#bullet MPV = %.2f  #pm %.2f", Land_Gaus_fit->GetParameter(1), Land_Gaus_fit->GetParError(1)));
		// tex11 -> DrawLatex(0.58, 0.650+0.05, Form("#bullet Land width = %.2f #pm %.2f", Land_Gaus_fit->GetParameter(2), Land_Gaus_fit->GetParError(2)));
		// tex11 -> DrawLatex(0.58, 0.600+0.05, Form("#bullet Gaus width = %.2f #pm %.2f", Land_Gaus_fit->GetParameter(3), Land_Gaus_fit->GetParError(3)));
		// tex11 -> DrawLatex(0.58, 0.550+0.05, Form("#bullet #chi^{2}/ndf = %.1f/%d = %.2f", Land_Gaus_fit->GetChisquare(), Land_Gaus_fit->GetNDF(), Land_Gaus_fit->GetChisquare()/Land_Gaus_fit->GetNDF()));

		// cout<<"Layer"<<i1<<" area ratio above threshold 15  : "<<func_fit->Integral(15,1000)/func_fit->Integral(0,1000)<<endl;

		c2->Print( Form("%s/l%d_DAC_scan_all.pdf",output_directory.Data(), i1) ); 
		c2 -> Clear();
	}
	// c2->Print( Form("%s/DAC_scan_all.pdf)",output_directory.Data()) );

	DAC_scan_out->Close();
}