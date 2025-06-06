#include "sPhenixStyle.C"

void Dist_comp() {
    std::string input_directory_1 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out";
    std::string input_filename_1 = "DAC_scan_out.root";

    std::string input_directory_2 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out_Scan2Scan3Match1Bin";
    std::string input_filename_2 = "DAC_scan_out.root";

    std::string input_directory_3 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out_UsedOverflow";
    std::string input_filename_3 = "DAC_scan_out.root";

    std::string input_directory_4 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/DACScan_out_AllOnly6thBin";
    std::string input_filename_4 = "DAC_scan_out.root";

    SetsPhenixStyle();

    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleAlign(22);
    gStyle->SetTitleX(0.5);
	gStyle->SetTitleY(0.95);
	gStyle->SetTitleFillColor(10);
	gStyle->SetPadRightMargin(0.09);
    gStyle->SetPadTopMargin(0.1);

    TCanvas * c1 = new TCanvas("c1", "c1", 950, 800);

    

    TLegend * leg = new TLegend(0.4,0.72,0.8,0.84);
	leg -> SetBorderSize(0);
	leg -> SetMargin(0.15);
	leg -> SetTextSize(0.025);

    TFile * file_in_1 = TFile::Open(Form("%s/%s", input_directory_1.c_str(), input_filename_1.c_str()));
    TFile * file_in_2 = TFile::Open(Form("%s/%s", input_directory_2.c_str(), input_filename_2.c_str()));
    TFile * file_in_3 = TFile::Open(Form("%s/%s", input_directory_3.c_str(), input_filename_3.c_str()));
    TFile * file_in_4 = TFile::Open(Form("%s/%s", input_directory_4.c_str(), input_filename_4.c_str()));

    TH1D * h1D_DAC_hist_all_1[3];
    TH1D * h1D_DAC_hist_all_2[3];
    TH1D * h1D_DAC_hist_all_3[3];
    TH1D * h1D_DAC_hist_all_4[3];

    h1D_DAC_hist_all_1[0] = (TH1D*) file_in_1 -> Get("DAC_hist_all_0");
    h1D_DAC_hist_all_1[1] = (TH1D*) file_in_1 -> Get("DAC_hist_all_1");
    h1D_DAC_hist_all_1[2] = (TH1D*) file_in_1 -> Get("DAC_hist_all_2"); 

    h1D_DAC_hist_all_2[0] = (TH1D*) file_in_2 -> Get("DAC_hist_all_0");
    h1D_DAC_hist_all_2[1] = (TH1D*) file_in_2 -> Get("DAC_hist_all_1");
    h1D_DAC_hist_all_2[2] = (TH1D*) file_in_2 -> Get("DAC_hist_all_2"); 

    h1D_DAC_hist_all_3[0] = (TH1D*) file_in_3 -> Get("DAC_hist_all_0");
    h1D_DAC_hist_all_3[1] = (TH1D*) file_in_3 -> Get("DAC_hist_all_1");
    h1D_DAC_hist_all_3[2] = (TH1D*) file_in_3 -> Get("DAC_hist_all_2");

    h1D_DAC_hist_all_4[0] = (TH1D*) file_in_4 -> Get("DAC_hist_all_0");
    h1D_DAC_hist_all_4[1] = (TH1D*) file_in_4 -> Get("DAC_hist_all_1");
    h1D_DAC_hist_all_4[2] = (TH1D*) file_in_4 -> Get("DAC_hist_all_2"); 

    leg->AddEntry(h1D_DAC_hist_all_1[0], "Baseline (two bins for dist. matching)", "ep");
    leg->AddEntry(h1D_DAC_hist_all_2[0], "Scan2 matched by only the 6th bin", "ep");
    leg->AddEntry(h1D_DAC_hist_all_3[0], "Overflow bins included for matching", "ep");
    leg->AddEntry(h1D_DAC_hist_all_4[0], "Only first overlapped bin used", "ep");

    for (int i = 0; i < 3; i++) {
        c1 -> cd();

        h1D_DAC_hist_all_1[i] -> SetMarkerStyle(20);
        h1D_DAC_hist_all_1[i] -> SetMarkerSize(0.8);
        h1D_DAC_hist_all_1[i] -> SetMarkerColor(1);
        h1D_DAC_hist_all_1[i] -> SetLineColor(1);
        h1D_DAC_hist_all_1[i] -> SetMaximum(2200);
        h1D_DAC_hist_all_1[i] -> SetTitle(Form("DAC_hist_all_%d", i));

        h1D_DAC_hist_all_2[i] -> SetMarkerStyle(22);
        h1D_DAC_hist_all_2[i] -> SetMarkerSize(0.8);
        h1D_DAC_hist_all_2[i] -> SetMarkerColor(2);
        h1D_DAC_hist_all_2[i] -> SetLineColor(2);

        h1D_DAC_hist_all_3[i] -> SetMarkerStyle(21);
        h1D_DAC_hist_all_3[i] -> SetMarkerSize(0.8);
        h1D_DAC_hist_all_3[i] -> SetMarkerColor(8);
        h1D_DAC_hist_all_3[i] -> SetLineColor(8);
        h1D_DAC_hist_all_3[i] -> SetFillColorAlpha(8,0.5);

        h1D_DAC_hist_all_4[i] -> SetMarkerStyle(28);
        h1D_DAC_hist_all_4[i] -> SetMarkerSize(0.8);
        h1D_DAC_hist_all_4[i] -> SetMarkerColor(kBlue);
        h1D_DAC_hist_all_4[i] -> SetLineColor(kBlue);
        // h1D_DAC_hist_all_3[i] -> SetFillColorAlpha(8,0.5);

        h1D_DAC_hist_all_1[i] -> Draw("ep");
        h1D_DAC_hist_all_2[i] -> Draw("ep same");
        h1D_DAC_hist_all_3[i] -> Draw("ep same");
        h1D_DAC_hist_all_4[i] -> Draw("ep same");
        leg -> Draw("same");

        c1 -> Print(Form("%s/Dist_comp_%d.pdf", input_directory_3.c_str(), i));
        c1 -> Clear();
    }

}