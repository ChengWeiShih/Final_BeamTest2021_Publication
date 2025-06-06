#include "../DetectionEffiAna.h"
R__LOAD_LIBRARY(../libDetectionEffiAna.so)

std::string data_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run72";
std::string data_input_filename = "run72_no_clone_filter_all_clusters.root";
std::string data_output_directory = data_input_directory + "/DetectionEffi_baseline700um";

std::tuple<std::string, double, double, double> Run72_mother(
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

    DetectionEffiAna * Run72Data = new DetectionEffiAna(
        true,
        72,
        data_input_directory,
        data_input_filename,
        data_output_directory,
        data_output_file_name_suffix,

        Selected_Column // Selected_Column_in, // note : 1 to 13;
         // l1_alignment_correction_in,
         // l0l2_slope_correction_in
    );

    Run72Data -> Get_l1_alignment();
    Run72Data -> Get_l0l2_slope();
    Run72Data -> GetDetectionEffi(
        Effi_slope_cut, // note : slope cut
        Effi_pos_cut, // note : noise hit distance
        Effi_boundary_cut, // note : boundary cut
        Effi_cluster_adc_value_requirement // note : cluster adc value requirement
    );

    if (prepare_ClusDist) {
        Run72Data -> GetControlClusDist(
            ClusDist_slope_cut, // note : slope cut
            ClusDist_pos_cut // note : noise hit distance
        );
    }

    Run72Data -> EndRun();

    return {
        Run72Data -> GetOutputFileName(),
        Run72Data -> Get_Final_Effi()[0],
        Run72Data -> Get_Final_Effi()[1],
        Run72Data -> Get_Final_Effi()[2]
    };
}



int Run72() {

    string BaseLine_data_output_file_name_suffix = "BaseLine";
    int    BaseLine_Selected_Column = 10; // note : 1 to 13,
    double BaseLine_Effi_slope_cut = 0.01; // note : slope cut
    double BaseLine_Effi_pos_cut = 0.7; // note : noise hit distance
    int    BaseLine_Effi_boundary_cut = 8;// note : boundary cut
    double BaseLine_Effi_cluster_adc_value_requirement = 0; // note : cluster adc value requirement

    double BaseLine_ClusDist_slope_cut = 0.01; // note : slope cut
    double BaseLine_ClusDist_pos_cut = 0.7; // note : noise hit distance

    bool   BaseLine_prepare_ClusDist = true;
    
    // Division: BaseLine ---------------------------------------------------------------------------------
    std::tuple<std::string, double, double, double> tuple_BaseLine = Run72_mother(
        BaseLine_data_output_file_name_suffix,
        BaseLine_Selected_Column,
        BaseLine_Effi_slope_cut,
        BaseLine_Effi_pos_cut,
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        BaseLine_prepare_ClusDist
    );

    return 888;
}