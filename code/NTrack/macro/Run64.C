#include "../NTrackAna.h"
#include "draw_style.h"


R__LOAD_LIBRARY(../libNTrackAna.so)

int Run64_temp(
    double noise_hit_distance,
    bool isColumnRejected
) {

    // double noise_hit_distance = 0.75; // note : unit : mm
    // bool isColumnRejected = false;
    std::vector<int> skip_column_vec = {4, 6};

    std::string output_plot_name_suffix = "";
    output_plot_name_suffix += (isColumnRejected) ? "_ColumnRejected" : "_AllColumn";
    output_plot_name_suffix += Form("_noiseHitDistance_%.0fum", noise_hit_distance * 1000);

    std::string data_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run64";
    std::string data_input_filename = "run64_no_clone_filter_all_clusters.root";
    std::string data_output_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run64/NTrack";
    std::string data_output_file_name_suffix = "test1";

    NTrackAna * data_run64 = new NTrackAna(
        true, 
        64,
        data_input_directory,
        data_input_filename,
        data_output_directory,
        data_output_file_name_suffix
    );
    if (isColumnRejected) {
        data_run64->SetSkipColumn(
            skip_column_vec
        );
    }

    data_run64->Get_l1_alignment();
    data_run64->GetNTrack(noise_hit_distance);

    TH1D * data_h1D_Column_NTrack = data_run64 -> Get_h1D_Column_NTrack();
    TH1D * data_h1D_NTrack_NoZero = data_run64 -> Get_h1D_NTrack_NoZero();
    TH1D * data_h1D_NTrackWithZero = data_run64 -> Get_h1D_NTrackWithZero();
    TH1D * data_h1D_l0_ClusPosFitDiff = data_run64 -> Get_h1D_l0_ClusPosFitDiff();
    TH1D * data_h1D_l1_ClusPosFitDiff = data_run64 -> Get_h1D_l1_ClusPosFitDiff();
    TH1D * data_h1D_l2_ClusPosFitDiff = data_run64 -> Get_h1D_l2_ClusPosFitDiff();
    TH1D * data_h1D_l1_residual = data_run64 -> Get_h1D_l1_residual();
    



    std::string MC_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run64/MC";
    std::string MC_input_filename = "cluster_information_offset-0.0000_adcinfo_SingleTrigger.root";
    std::string MC_output_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run64/MC/NTrack";
    std::string MC_output_file_name_suffix = "test1";

    NTrackAna * MC_run64 = new NTrackAna(
        false, 
        64,
        MC_input_directory,
        MC_input_filename,
        MC_output_directory,
        MC_output_file_name_suffix
    );
    MC_run64 -> SetBeamSpotColumn(10);
    if (isColumnRejected) {
        MC_run64 -> SetSkipColumn(
            skip_column_vec
        );
    }

    MC_run64->Get_l1_alignment();
    MC_run64->GetNTrack(noise_hit_distance);

    TH1D * MC_h1D_Column_NTrack = MC_run64 -> Get_h1D_Column_NTrack();
    TH1D * MC_h1D_NTrack_NoZero = MC_run64 -> Get_h1D_NTrack_NoZero();
    TH1D * MC_h1D_NTrackWithZero = MC_run64 -> Get_h1D_NTrackWithZero();
    TH1D * MC_h1D_l0_ClusPosFitDiff = MC_run64 -> Get_h1D_l0_ClusPosFitDiff();
    TH1D * MC_h1D_l1_ClusPosFitDiff = MC_run64 -> Get_h1D_l1_ClusPosFitDiff();
    TH1D * MC_h1D_l2_ClusPosFitDiff = MC_run64 -> Get_h1D_l2_ClusPosFitDiff();
    TH1D * MC_h1D_l1_residual = MC_run64 -> Get_h1D_l1_residual();






    TFile * file_out = new TFile(Form("%s/Run64_dataMC_comp.root", MC_output_directory.c_str()), "RECREATE");
    TCanvas * c2 = new TCanvas("c2", "c2", 800, 600);
    c2->cd();

    c2 -> cd();
    data_h1D_Column_NTrack -> SetMarkerStyle(20);
    data_h1D_Column_NTrack -> SetMarkerSize(0.8);
    data_h1D_Column_NTrack -> SetMarkerColor(kBlack);
    data_h1D_Column_NTrack -> Scale(1. / data_h1D_Column_NTrack -> Integral());
    MC_h1D_Column_NTrack -> Scale(1. / MC_h1D_Column_NTrack -> Integral());
    data_h1D_Column_NTrack -> Draw("ep");
    MC_h1D_Column_NTrack -> Draw("hist same");
    c2 -> Write("c1_h1D_Column_NTrack");
    c2 -> Clear();

    c2 -> cd();
    c2 -> SetLogy();
    data_h1D_NTrack_NoZero -> SetMarkerStyle(20);
    data_h1D_NTrack_NoZero -> SetMarkerSize(0.8);
    data_h1D_NTrack_NoZero -> SetMarkerColor(kBlack);
    data_h1D_NTrack_NoZero -> Scale(1. / data_h1D_NTrack_NoZero -> Integral());
    MC_h1D_NTrack_NoZero -> Scale(1. / MC_h1D_NTrack_NoZero -> Integral());
    data_h1D_NTrack_NoZero -> Draw("ep");
    MC_h1D_NTrack_NoZero -> Draw("hist same");
    c2 -> Write("c1_h1D_NTrack_NoZero");
    c2 -> Clear();

    c2 -> cd();
    c2 -> SetLogy();
    data_h1D_NTrackWithZero -> SetMarkerStyle(20);
    data_h1D_NTrackWithZero -> SetMarkerSize(0.8);
    data_h1D_NTrackWithZero -> SetMarkerColor(kBlack);
    data_h1D_NTrackWithZero -> Scale(1. / data_h1D_NTrackWithZero -> Integral());
    MC_h1D_NTrackWithZero -> Scale(1. / MC_h1D_NTrackWithZero -> Integral());
    data_h1D_NTrackWithZero -> Draw("ep");
    MC_h1D_NTrackWithZero -> Draw("hist same");
    c2 -> Write("c1_h1D_NTrackWithZero");
    c2 -> Clear();

    c2 -> cd();
    c2 -> SetLogy();
    data_h1D_l0_ClusPosFitDiff -> SetMarkerStyle(20);
    data_h1D_l0_ClusPosFitDiff -> SetMarkerSize(0.8);
    data_h1D_l0_ClusPosFitDiff -> SetMarkerColor(kBlack);
    data_h1D_l0_ClusPosFitDiff -> Scale(1. / data_h1D_l0_ClusPosFitDiff -> Integral());
    MC_h1D_l0_ClusPosFitDiff -> Scale(1. / MC_h1D_l0_ClusPosFitDiff -> Integral());
    data_h1D_l0_ClusPosFitDiff -> Draw("ep");
    MC_h1D_l0_ClusPosFitDiff -> Draw("hist same");
    c2 -> Write("c1_h1D_l0_ClusPosFitDiff");
    c2 -> Clear();

    c2 -> cd();
    c2 -> SetLogy();
    data_h1D_l1_ClusPosFitDiff -> SetMarkerStyle(20);
    data_h1D_l1_ClusPosFitDiff -> SetMarkerSize(0.8);
    data_h1D_l1_ClusPosFitDiff -> SetMarkerColor(kBlack);
    data_h1D_l1_ClusPosFitDiff -> Scale(1. / data_h1D_l1_ClusPosFitDiff -> Integral());
    MC_h1D_l1_ClusPosFitDiff -> Scale(1. / MC_h1D_l1_ClusPosFitDiff -> Integral());
    data_h1D_l1_ClusPosFitDiff -> Draw("ep");
    MC_h1D_l1_ClusPosFitDiff -> Draw("hist same");
    c2 -> Write("c1_h1D_l1_ClusPosFitDiff");
    c2 -> Clear();

    c2 -> cd();
    c2 -> SetLogy();
    data_h1D_l2_ClusPosFitDiff -> SetMarkerStyle(20);
    data_h1D_l2_ClusPosFitDiff -> SetMarkerSize(0.8);
    data_h1D_l2_ClusPosFitDiff -> SetMarkerColor(kBlack);
    data_h1D_l2_ClusPosFitDiff -> Scale(1. / data_h1D_l2_ClusPosFitDiff -> Integral());
    MC_h1D_l2_ClusPosFitDiff -> Scale(1. / MC_h1D_l2_ClusPosFitDiff -> Integral());
    data_h1D_l2_ClusPosFitDiff -> Draw("ep");
    MC_h1D_l2_ClusPosFitDiff -> Draw("hist same");
    c2 -> Write("c1_h1D_l2_ClusPosFitDiff");
    c2 -> Clear();
    
    file_out->Close();


    dataMC_comp(data_h1D_Column_NTrack, MC_h1D_Column_NTrack, MC_output_directory, Form("c1_h1D_Column_NTrack_%s",output_plot_name_suffix.c_str()), {"Column ID", "Number of reco. tracks (A.U.)"}, false, false);
    dataMC_comp(data_h1D_NTrack_NoZero, MC_h1D_NTrack_NoZero, MC_output_directory, Form("c1_h1D_NTrack_NoZero_%s",output_plot_name_suffix.c_str()), {"Number of reco. tracks","Entries (A.U.)"}, true, false);
    // dataMC_comp(data_h1D_NTrackWithZero, MC_h1D_NTrackWithZero, MC_output_directory, Form("c1_h1D_NTrackWithZero_%s",output_plot_name_suffix.c_str()), {"Number of reco. tracks","Entries (A.U.)"}, true, false);
    // dataMC_comp(data_h1D_l0_ClusPosFitDiff, MC_h1D_l0_ClusPosFitDiff, MC_output_directory, Form("c1_h1D_l0_ClusPosFitDiff_%s",output_plot_name_suffix.c_str()), {"ClusPos - Fit [mm] (L0)","Entries (A.U.)"}, true, false);
    // dataMC_comp(data_h1D_l1_ClusPosFitDiff, MC_h1D_l1_ClusPosFitDiff, MC_output_directory, Form("c1_h1D_l1_ClusPosFitDiff_%s",output_plot_name_suffix.c_str()), {"ClusPos - Fit [mm] (L1)","Entries (A.U.)"}, true, false);
    // dataMC_comp(data_h1D_l2_ClusPosFitDiff, MC_h1D_l2_ClusPosFitDiff, MC_output_directory, Form("c1_h1D_l2_ClusPosFitDiff_%s",output_plot_name_suffix.c_str()), {"ClusPos - Fit [mm] (L2)","Entries (A.U.)"}, true, false);
    dataMC_comp(data_h1D_l1_residual, MC_h1D_l1_residual, MC_output_directory, Form("c1_h1D_l1_residual_%s",output_plot_name_suffix.c_str()), {"L1 - (L2L0 interpolation) [mm]","Entries (A.U.)"}, true, false);

    data_run64->EndRun();
    MC_run64->EndRun();

    return 111;
}

void Run64(){
    // Run64_temp(0.75, false);
    // Run64_temp(0.75, true);

    // Run64_temp(0.234, false);
    // sleep(5);

    // Run64_temp(0.234, true);
    // sleep(5);

    // Run64_temp(0.5, false);
    // sleep(5);

    // Run64_temp(0.5, true);
    // sleep(5);

    // Run64_temp(0.75, false);
    // sleep(5);

    // Run64_temp(0.75, true);
    // sleep(5);

    // Run64_temp(0.65, false);
    // sleep(5);

    // Run64_temp(0.65, true);
    // sleep(5);

    // Run64_temp(0.85, false);
    // sleep(5);

    // Run64_temp(0.85, true);
    // sleep(5);

    Run64_temp(0.70, true);
    sleep(5);

    Run64_temp(0.70, false);
    // sleep(5);
}