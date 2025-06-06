#include "../DACScanTight.h"
R__LOAD_LIBRARY(../libDACScanTight.so)

// std::string data_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run72";
// std::string data_input_filename = "run72_no_clone_filter_all_clusters.root";
// std::string data_output_directory = data_input_directory + "/DetectionEffi_baseline700um";

std::tuple<std::string, double, double, double> RunDACScan_mother(
    int runnumber,
    std::string data_input_directory,
    std::string data_input_filename,
    std::string data_output_directory,

    double L1Residual_Alignment_corr = 0.,
    double L0L2Slope_correction = 0.,

    std::string data_output_file_name_suffix = "BaseLine",    
    int Selected_Column = 8, // note : 1 to 13,
    double Effi_slope_cut = 0.01, // note : slope cut
    double Effi_pos_cut = 0.4, // note : noise hit distance
    int    Effi_boundary_cut = 5, // note : boundary cut
    double Effi_cluster_adc_value_requirement = 15, // note : cluster adc value requirement

    double ClusDist_slope_cut = 0.01, // note : slope cut
    double ClusDist_pos_cut = 0.234, // note : noise hit distance

    bool prepare_ClusDist = false
){

    DACScanTight * RunData = new DACScanTight(
        true,
        runnumber,
        data_input_directory,
        data_input_filename,
        data_output_directory,
        data_output_file_name_suffix,

        Selected_Column, // Selected_Column_in, // note : 1 to 13;
        L1Residual_Alignment_corr, // l1_alignment_correction_in,
        L0L2Slope_correction // l0l2_slope_correction_in
    );

    RunData -> Get_l1_alignment();
    RunData -> Get_l0l2_slope();
    RunData -> GetDetectionEffi(
        Effi_slope_cut, // note : slope cut
        Effi_pos_cut, // note : noise hit distance
        Effi_boundary_cut, // note : boundary cut
        Effi_cluster_adc_value_requirement // note : cluster adc value requirement
    );

    if (prepare_ClusDist) {
        RunData -> GetControlClusDist(
            ClusDist_slope_cut, // note : slope cut
            ClusDist_pos_cut // note : noise hit distance
        );
    }

    RunData -> EndRun();

    return {
        RunData -> GetOutputFileName(),
        RunData -> Get_Final_Effi()[0],
        RunData -> Get_Final_Effi()[1],
        RunData -> Get_Final_Effi()[2]
    };
}



int RunDACScan() {

    string BaseLine_data_output_file_name_suffix = "BaseLine";
    
    int    BaseLine_Selected_Column1 = 9; // note : 1 to 13,
    int    BaseLine_Selected_Column2 = 10; // note : 1 to 13,
    int    BaseLine_Selected_Column3 = 11; // note : 1 to 13,

    double BaseLine_Effi_slope_cut = 0.01; // note : slope cut
    double BaseLine_Effi_pos_cut = 0.7; // note : noise hit distance
    int    BaseLine_Effi_boundary_cut = 8;// note : boundary cut
    double BaseLine_Effi_cluster_adc_value_requirement = 0; // note : cluster adc value requirement

    double BaseLine_ClusDist_slope_cut = 0.01; // note : slope cut
    double BaseLine_ClusDist_pos_cut = 0.7; // note : noise hit distance

    bool   BaseLine_prepare_ClusDist = true;

    // todo: numbers are from Run76
    // note : Column9
    double Column1_L1Residual_Alignment_corr = -0.31483618;
    double Column1_L0L2Slope_correction = 0.0048631790;

    // note : Column10
    double Column2_L1Residual_Alignment_corr = -0.29968099;
    double Column2_L0L2Slope_correction = 0.0046633715;

    // note : Column11
    double Column3_L1Residual_Alignment_corr = -0.26356519;
    double Column3_L0L2Slope_correction = 0.0045548639;

    for (int i = 71; i < 79; i++) {
        int runnumber = i;
        std::string data_input_directory = Form("/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan");
        std::string data_input_filename = Form("run%d_no_clone_filter_all_clusters.root",i);
        std::string data_output_directory = data_input_directory + "/TightSelection";


        // Division: BaseLine ---------------------------------------------------------------------------------
        RunDACScan_mother(
            runnumber,
            data_input_directory,
            data_input_filename,
            data_output_directory,
            Column1_L1Residual_Alignment_corr,
            Column1_L0L2Slope_correction,

            BaseLine_data_output_file_name_suffix,
            BaseLine_Selected_Column1,
            BaseLine_Effi_slope_cut,
            BaseLine_Effi_pos_cut,
            BaseLine_Effi_boundary_cut,
            BaseLine_Effi_cluster_adc_value_requirement,

            BaseLine_ClusDist_slope_cut,
            BaseLine_ClusDist_pos_cut,

            BaseLine_prepare_ClusDist
        );

        RunDACScan_mother(
            runnumber,
            data_input_directory,
            data_input_filename,
            data_output_directory,
            Column2_L1Residual_Alignment_corr,
            Column2_L0L2Slope_correction,

            BaseLine_data_output_file_name_suffix,
            BaseLine_Selected_Column2,
            BaseLine_Effi_slope_cut,
            BaseLine_Effi_pos_cut,
            BaseLine_Effi_boundary_cut,
            BaseLine_Effi_cluster_adc_value_requirement,

            BaseLine_ClusDist_slope_cut,
            BaseLine_ClusDist_pos_cut,

            BaseLine_prepare_ClusDist
        );

        RunDACScan_mother(
            runnumber,
            data_input_directory,
            data_input_filename,
            data_output_directory,
            Column3_L1Residual_Alignment_corr,
            Column3_L0L2Slope_correction,

            BaseLine_data_output_file_name_suffix,
            BaseLine_Selected_Column3,
            BaseLine_Effi_slope_cut,
            BaseLine_Effi_pos_cut,
            BaseLine_Effi_boundary_cut,
            BaseLine_Effi_cluster_adc_value_requirement,

            BaseLine_ClusDist_slope_cut,
            BaseLine_ClusDist_pos_cut,

            BaseLine_prepare_ClusDist
        );

    }

    return 888;
}