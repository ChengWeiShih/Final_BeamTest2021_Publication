#include "sPhenixStyle.C"

int CheckColumn(){

    std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan";
    std::string input_filename = "run77_no_clone_filter_all_clusters.root";

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

    TFile * file_in = TFile::Open(Form("%s/%s", input_directory.c_str(), input_filename.c_str()));
    TTree * tree_in = (TTree*) file_in -> Get("cluster_info");

    int DSE;
	int eID;
	vector<int>* layer = 0;
	vector<int>* chip = 0;
	vector<int>* Nhit = 0;
	vector<double>* position = 0;
	vector<double>* raw_data = 0;
	vector<double>* cluster_adc = 0;

    tree_in -> SetBranchAddress("layer", &layer);
    tree_in -> SetBranchAddress("chip", &chip);
    tree_in -> SetBranchAddress("Nhit", &Nhit);
    tree_in -> SetBranchAddress("position", &position);
    tree_in -> SetBranchAddress("raw_data", &raw_data);
    tree_in -> SetBranchAddress("cluster_adc", &cluster_adc);
    tree_in -> SetBranchAddress("eID",&eID);
    tree_in -> SetBranchAddress("DSE",&DSE);

    TH1D * h1D_1HitClusAdc_U9 = new TH1D("h1D_1HitClusAdc_U9", "h1D_1HitClusAdc_U9;Hit adc;Entries (A. U.)", 8, -0.5, 7.5);
    TH1D * h1D_1HitClusAdc_U10 = new TH1D("h1D_1HitClusAdc_U10", "h1D_1HitClusAdc_U10;Hit adc;Entries (A. U.)", 8, -0.5, 7.5);
    TH1D * h1D_1HitClusAdc_U11 = new TH1D("h1D_1HitClusAdc_U11", "h1D_1HitClusAdc_U11;Hit adc;Entries (A. U.)", 8, -0.5, 7.5);

    h1D_1HitClusAdc_U9 -> SetLineColor(kRed);
    h1D_1HitClusAdc_U10 -> SetLineColor(kBlue);
    h1D_1HitClusAdc_U11 -> SetLineColor(kGreen);

    h1D_1HitClusAdc_U9 -> SetLineWidth(2);
    h1D_1HitClusAdc_U10 -> SetLineWidth(2);
    h1D_1HitClusAdc_U11 -> SetLineWidth(2);

    for (int i = 0; i < tree_in -> GetEntries(); i++) {
        tree_in -> GetEntry(i);
        
        if (DSE != 0){continue;}

        if (true)
        {
            int l0_Clus_count = 0;
            int l1_Clus_count = 0;
            int l2_Clus_count = 0;

            for (int i3=0; i3<layer->size(); i3++)
            {
                if (layer->at(i3) == 0){l0_Clus_count++;}
                else if (layer->at(i3) == 1){l1_Clus_count++;}
                else if (layer->at(i3) == 2){l2_Clus_count++;}
            }

            if (l0_Clus_count > 1 || l1_Clus_count > 1 || l2_Clus_count > 1)
            {
                continue;
            }
        }

        for (int clus_i = 0; clus_i < layer->size(); clus_i++) {
            if (Nhit->at(clus_i) != 1){continue;}

            if (chip->at(clus_i) == 9)
            {
                h1D_1HitClusAdc_U9 -> Fill(nominal_setting_map[cluster_adc->at(clus_i)]);
            }
            else if (chip->at(clus_i) == 10)
            {
                h1D_1HitClusAdc_U10 -> Fill(nominal_setting_map[cluster_adc->at(clus_i)]);
            }
            else if (chip->at(clus_i) == 11)
            {
                h1D_1HitClusAdc_U11 -> Fill(nominal_setting_map[cluster_adc->at(clus_i)]);
            }

        }
    }

    h1D_1HitClusAdc_U9 -> SetBinContent(8,0);
    h1D_1HitClusAdc_U10 -> SetBinContent(8,0);
    h1D_1HitClusAdc_U11 -> SetBinContent(8,0);

    SetsPhenixStyle();
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleAlign(22);
    gStyle->SetTitleX(0.5);
	gStyle->SetTitleY(0.95);
	gStyle->SetTitleFillColor(10);
	gStyle->SetPadTopMargin(0.09);

    TCanvas * c1 = new TCanvas("c1", "c1", 950, 800);
    h1D_1HitClusAdc_U9 -> Scale(1./h1D_1HitClusAdc_U9 -> Integral());
    h1D_1HitClusAdc_U10 -> Scale(1./h1D_1HitClusAdc_U10 -> Integral());
    h1D_1HitClusAdc_U11 -> Scale(1./h1D_1HitClusAdc_U11 -> Integral());

    TLegend * leg = new TLegend(0.6, 0.25, 0.75, 0.45);
    leg -> AddEntry(h1D_1HitClusAdc_U9, "U9", "l");
    leg -> AddEntry(h1D_1HitClusAdc_U10, "U10", "l");
    leg -> AddEntry(h1D_1HitClusAdc_U11, "U11", "l");

    h1D_1HitClusAdc_U9 -> SetTitle("run77_h1D_1HitClusAdc_U9_U10_U11");
    h1D_1HitClusAdc_U9 -> Draw("hist");
    h1D_1HitClusAdc_U10 -> Draw("hist same");
    h1D_1HitClusAdc_U11 -> Draw("hist same");
    leg -> Draw("same");

    c1 -> Print(Form("%s/run77_h1D_1HitClusAdc_U9_U10_U11.pdf", input_directory.c_str()));

    return 888;
}