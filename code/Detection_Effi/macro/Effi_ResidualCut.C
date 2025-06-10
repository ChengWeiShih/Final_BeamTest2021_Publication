#include "sPhenixStyle.C"

int Effi_ResidualCut() {
    
    string data_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/DetectionEffi_test";
    string data_input_filename = "Data_DetectionEffi_Run52_Column8_BaseLine.root";

    string MC_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/DetectionEffi_test/MC/DetectionEffi_test";
    string MC_input_filename = "MC_DetectionEffi_Run52_Column10_BaseLine.root";

    SetsPhenixStyle();
    TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);

    TLegend *legend1 = new TLegend(0.7, 0.45, 0.85, 0.6);
    legend1 -> SetBorderSize(0);
	legend1 -> SetTextSize(0.050);
	// legend1 -> SetNColumns (4);

    std::pair<double, double> boundary_cut = {-6, 6}; // note : unit: mm, for run 52

    TFile * data_file_in = TFile::Open(Form("%s/%s", data_input_directory.c_str(), data_input_filename.c_str()));
    TTree * data_tree_in = (TTree *)data_file_in->Get("tree_effi");
    double data_L0L2Interpolation;
    double data_L1Residual;
    bool data_L1Good;

    data_tree_in->SetBranchAddress("L0L2Interpolation", &data_L0L2Interpolation);
    data_tree_in->SetBranchAddress("L1Residual", &data_L1Residual);
    data_tree_in->SetBranchAddress("L1Good", &data_L1Good);


    TFile * MC_file_in = TFile::Open(Form("%s/%s", MC_input_directory.c_str(), MC_input_filename.c_str()));
    TTree * MC_tree_in = (TTree *)MC_file_in->Get("tree_effi");
    double MC_L0L2Interpolation;
    double MC_L1Residual;
    bool MC_L1Good;

    MC_tree_in->SetBranchAddress("L0L2Interpolation", &MC_L0L2Interpolation);
    MC_tree_in->SetBranchAddress("L1Residual", &MC_L1Residual);
    MC_tree_in->SetBranchAddress("L1Good", &MC_L1Good);


    TGraphAsymmErrors * data_gr = new TGraphAsymmErrors();
    TGraphAsymmErrors * MC_gr = new TGraphAsymmErrors();

    TH1D * h1D_data_denominator = new TH1D("h1D_data_denominator", "h1D_data_denominator", 1,0,1);
    TH1D * h1D_data_numerator = new TH1D("h1D_data_numerator", "h1D_data_numerator", 1,0,1);
    
    TH1D * h1D_MC_denominator = new TH1D("h1D_MC_denominator", "h1D_MC_denominator", 1,0,1);
    TH1D * h1D_MC_numerator = new TH1D("h1D_MC_numerator", "h1D_MC_numerator", 1,0,1);

    double data_denominator = 0;
    double MC_denominator = 0;

    for (int i = 0; i < data_tree_in -> GetEntries(); i++) {
        data_tree_in -> GetEntry(i);
        if (data_L0L2Interpolation > boundary_cut.first && data_L0L2Interpolation < boundary_cut.second) {
            data_denominator++;
        }
    }

    for (int i = 0; i < MC_tree_in -> GetEntries(); i++) {
        MC_tree_in -> GetEntry(i);
        if (MC_L0L2Interpolation > boundary_cut.first && MC_L0L2Interpolation < boundary_cut.second) {
            MC_denominator++;
        }
    }

    h1D_data_denominator -> SetBinContent(1, data_denominator);
    h1D_MC_denominator -> SetBinContent(1, MC_denominator);    

    for (int i = 0; i < 34; i++) {
        double cut = 0.01 + i * 0.03;// note : mm

        double data_numerator = 0;
        double MC_numerator = 0;

        for (int j = 0; j < data_tree_in -> GetEntries(); j++) {
            data_tree_in -> GetEntry(j);
            if (fabs(data_L1Residual) < cut && data_L0L2Interpolation > boundary_cut.first && data_L0L2Interpolation < boundary_cut.second ) {
                data_numerator++;
            }
        }

        for (int j = 0; j < MC_tree_in -> GetEntries(); j++) {
            MC_tree_in -> GetEntry(j);
            if (fabs(MC_L1Residual) < cut && MC_L0L2Interpolation > boundary_cut.first && MC_L0L2Interpolation < boundary_cut.second) {
                MC_numerator++;
            }
        }

        h1D_data_numerator -> SetBinContent(1, data_numerator);
        h1D_MC_numerator -> SetBinContent(1, MC_numerator);

        TEfficiency * data_detection_effi = new TEfficiency (*h1D_data_numerator,*h1D_data_denominator);
        TEfficiency * MC_detection_effi = new TEfficiency (*h1D_MC_numerator,*h1D_MC_denominator);


        data_gr -> SetPoint(data_gr->GetN(), cut, data_detection_effi -> GetEfficiency(1) * 100.);
        data_gr -> SetPointError(
            data_gr->GetN()-1, 
            0, 0,
            data_detection_effi -> GetEfficiencyErrorLow(1) * 100.,
            data_detection_effi -> GetEfficiencyErrorUp(1) * 100.
        );

        MC_gr -> SetPoint(MC_gr->GetN(), cut, MC_detection_effi -> GetEfficiency(1) * 100.);
        MC_gr -> SetPointError(
            MC_gr->GetN()-1, 
            0, 0,
            MC_detection_effi -> GetEfficiencyErrorLow(1) * 100.,
            MC_detection_effi -> GetEfficiencyErrorUp(1) * 100.
        );
    }

    data_gr -> SetName("data_gr");
    data_gr -> SetMarkerStyle(20);
    data_gr -> SetMarkerSize(0.5);
    data_gr -> SetMarkerColor(kBlack);
    data_gr -> SetLineColor(kBlack);
    data_gr -> GetXaxis() -> SetTitle("L1Residual cut [mm]");
    data_gr -> GetYaxis() -> SetTitle("Detection Efficiency (%)");
    data_gr -> GetYaxis() -> SetRangeUser(75, 101);


    MC_gr -> SetName("MC_gr");
    MC_gr -> SetMarkerStyle(20);
    MC_gr -> SetMarkerSize(0.5);
    MC_gr -> SetMarkerColor(kBlue);
    MC_gr -> SetLineColor(kBlack);
    MC_gr -> GetXaxis() -> SetTitle("L1Residual cut [mm]");
    MC_gr -> GetYaxis() -> SetTitle("Detection Efficiency (%)");
    MC_gr -> GetYaxis() -> SetRangeUser(75, 101);


    for (int i = 0; i < data_gr -> GetN(); i++) {
        double x, y;
        double exl, exh, eyl, eyh;
        data_gr -> GetPoint(i, x, y);
        exl = data_gr -> GetErrorXlow(i);
        exh = data_gr -> GetErrorXhigh(i);
        eyl = data_gr -> GetErrorYlow(i);
        eyh = data_gr -> GetErrorYhigh(i);
        
        std::cout<<"data_gr, residual cut: "<<x<<", effi: "<<y<<", error: "<<eyl<<", "<<eyh<<std::endl;
    }

    for (int i = 0; i < MC_gr -> GetN(); i++) {
        double x, y;
        double exl, exh, eyl, eyh;
        MC_gr -> GetPoint(i, x, y);
        exl = MC_gr -> GetErrorXlow(i);
        exh = MC_gr -> GetErrorXhigh(i);
        eyl = MC_gr -> GetErrorYlow(i);
        eyh = MC_gr -> GetErrorYhigh(i);
        
        std::cout<<"MC_gr, residual cut: "<<x<<", effi: "<<y<<", error: "<<eyl<<", "<<eyh<<std::endl;
    }

    TGraphAsymmErrors * data_MC_correlation_gr = new TGraphAsymmErrors();
    data_MC_correlation_gr -> SetName("data_MC_correlation_gr");
    if (data_gr -> GetN() != MC_gr -> GetN()) {
        std::cout<<"data_gr and MC_gr have different number of points!"<<std::endl;
        return 0;
    }
    for (int i = 0; i < data_gr -> GetN(); i++) {
        double x, y;
        double exl, exh, eyl, eyh;
        data_gr -> GetPoint(i, x, y);
        exl = data_gr -> GetErrorXlow(i);
        exh = data_gr -> GetErrorXhigh(i);
        eyl = data_gr -> GetErrorYlow(i);
        eyh = data_gr -> GetErrorYhigh(i);

        double MC_x, MC_y;
        double MC_exl, MC_exh, MC_eyl, MC_eyh;
        MC_gr -> GetPoint(i, MC_x, MC_y);
        MC_exl = MC_gr -> GetErrorXlow(i);
        MC_exh = MC_gr -> GetErrorXhigh(i);
        MC_eyl = MC_gr -> GetErrorYlow(i);
        MC_eyh = MC_gr -> GetErrorYhigh(i);

        // std::cout<<y/MC_y<<std::endl;

        data_MC_correlation_gr -> SetPoint(data_MC_correlation_gr->GetN(), y, MC_y);
        data_MC_correlation_gr -> SetPointError(
            data_MC_correlation_gr->GetN()-1,
            eyl, eyh,
            MC_eyl, MC_eyh 
        );

        std::cout<<"data_MC_correlation_gr, data effi: "<<y<<", MC effi: "<<MC_y<<", error: "<<eyl<<", "<<eyh<<", MC error: "<<MC_eyl<<", "<<MC_eyh<<std::endl;
        // std::cout<<"data_MC_correlation_gr, errXL: "
        // <<data_MC_correlation_gr->GetErrorXlow(data_MC_correlation_gr->GetN()-1)
        // <<", errXH: "<<data_MC_correlation_gr->GetErrorXhigh(data_MC_correlation_gr->GetN()-1)
        // <<", errYL: "<<data_MC_correlation_gr->GetErrorYlow(data_MC_correlation_gr->GetN()-1)
        // <<", errYH: "<<data_MC_correlation_gr->GetErrorYhigh(data_MC_correlation_gr->GetN()-1)
        // <<std::endl;
    }

    
    data_MC_correlation_gr -> SetMarkerStyle(20);
    data_MC_correlation_gr -> SetMarkerSize(0.8);
    data_MC_correlation_gr -> SetMarkerColor(kBlack);
    data_MC_correlation_gr -> SetLineColor(kBlack);
    data_MC_correlation_gr -> GetXaxis() -> SetTitle("Data Detection Efficiency (%)");
    data_MC_correlation_gr -> GetYaxis() -> SetTitle("MC Detection Efficiency (%)");

    legend1->AddEntry(data_gr, "Data", "p");
    legend1->AddEntry(MC_gr, "MC", "p");

    TFile * file_out = new TFile(Form("%s/Effi_ResidualCut.root", MC_input_directory.c_str()), "RECREATE");
    data_gr -> Write();
    MC_gr -> Write();
    data_MC_correlation_gr -> Write();

    c1 -> cd();
    data_gr -> Draw("AP");
    MC_gr -> Draw("P same");
    c1 -> Write("c1_Effi_ResidualCut");
    
    legend1 -> Draw("same");
    c1 -> Print(Form("%s/c1_Effi_ResidualCut.pdf", MC_input_directory.c_str()));
    
    data_gr -> GetYaxis() -> SetRangeUser(99.6, 100.1);
    data_gr -> GetXaxis() -> SetRangeUser(0.7, 1.1);
    c1 -> Print(Form("%s/c1_Effi_ResidualCut_focus.pdf", MC_input_directory.c_str()));

    c1 -> Clear();

    c1 -> cd();
    data_MC_correlation_gr -> Draw("AP");
    c1 -> Write("c1_data_MC_correlation_gr");


    file_out -> Close();

    return 888;
}