#include "draw_style.h"

#include "../ResidualCompAna.h"
R__LOAD_LIBRARY(../libResidualCompAna.so)

int Run52_DataMC_comp()
{   
    // Division : ---------------------------------------------------------------------------------------------------------------------------

    std::string data_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52";
    std::string data_input_filename = "run52_no_clone_filter_all_clusters.root";
    std::string data_output_directory = data_input_directory + "/ResidualComp//Output";
    std::string data_output_file_name_suffix = "test1";

    ResidualCompAna * Run52Data = new ResidualCompAna(
        true,
        52,
        data_input_directory,
        data_input_filename,
        data_output_directory,
        data_output_file_name_suffix,

        8 // Selected_Column_in, // note : 1 to 13;
         // l1_alignment_correction_in,
         // l0l1_slope_correction_in,
         // ClusPhiSize_cut_in,
         // ClusAdc_cut_in
    );

    Run52Data -> Get_l1_alignment();
    Run52Data -> Get_l0l1_slope();
    Run52Data -> GetHistsForComp(
        0.01, // note : slope cut
        5.0   // note : pos cut
    );

    TH1D * data_h1D_l1_residual = (TH1D*) (Run52Data -> Get_h1D_l1_residual()) -> Clone("data_h1D_l1_residual");
    TH1D * data_h1D_scattering = (TH1D*) (Run52Data -> Get_h1D_scattering()) -> Clone("data_h1D_scattering");

    


    // Division : ---------------------------------------------------------------------------------------------------------------------------

    std::string MC_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/ResidualComp/MC";
    std::string MC_input_filename = "cluster_information_offset-0.0000_adcinfo_SingleTrigger.root";
    std::string MC_output_directory = MC_input_directory ;
    std::string MC_output_file_name_suffix = "test1";

    ResidualCompAna * Run52MC = new ResidualCompAna(
        false,
        52,
        MC_input_directory,
        MC_input_filename,
        MC_output_directory,
        MC_output_file_name_suffix,

        11 // Selected_Column_in, // note : 1 to 13;
         // l1_alignment_correction_in,
         // l0l1_slope_correction_in,
         // ClusPhiSize_cut_in,
         // ClusAdc_cut_in
    );

    Run52MC -> Get_l1_alignment();
    Run52MC -> Get_l0l1_slope();
    Run52MC -> GetHistsForComp(
        0.01, // note : slope cut
        5.0   // note : pos cut
    );

    TH1D * MC_h1D_l1_residual = (TH1D*) (Run52MC -> Get_h1D_l1_residual()) -> Clone("MC_h1D_l1_residual");
    TH1D * MC_h1D_scattering = (TH1D*) (Run52MC -> Get_h1D_scattering()) -> Clone("MC_h1D_scattering");

    


    // Division : ---------------------------------------------------------------------------------------------------------------------------    


    pair<double,double> ratio_Y_range_pair = {-3. + 1, 3 + 1};
    TString MC_selected_Column = "U11";

    // vector<TString> titles_scattering_vec = {"Scattering [slope_{l2l1} - slope_{l1l0}]","A.U."};
    // TString plot_scattering_name = Form("DataRun52U8_MC%s_scattering_",MC_selected_Column.Data());
    // dataMC_comp (data_scattering_hist, MC_scattering_hist, MC_folder_direction , plot_scattering_name, titles_scattering_vec, ratio_Y_range_pair,true, false);
    
    vector<TString> titles_residual_vec = {"L1 - (L2L0 interpolation) [mm]","Entries (A.U.)"};
    TString plot_residual_name = Form("DataRun52U8_MC%s_residual_",MC_selected_Column.Data());
    dataMC_comp (data_h1D_l1_residual, MC_h1D_l1_residual, MC_output_directory , plot_residual_name, titles_residual_vec, ratio_Y_range_pair,true, false);

    Run52Data -> EndRun();
    Run52MC -> EndRun();

    return 888;
}