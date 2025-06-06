#include "../DetectionEffiAna.h"
R__LOAD_LIBRARY(../libDetectionEffiAna.so)

std::string data_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run77";
std::string data_input_filename = "run77_no_clone_filter_all_clusters.root";
std::string data_output_directory = data_input_directory + "/DetectionEffi_baseline700um";

std::tuple<std::string, double, double, double> Run77_mother(
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

    DetectionEffiAna * Run77Data = new DetectionEffiAna(
        true,
        77,
        data_input_directory,
        data_input_filename,
        data_output_directory,
        data_output_file_name_suffix,

        Selected_Column // Selected_Column_in, // note : 1 to 13;
         // l1_alignment_correction_in,
         // l0l2_slope_correction_in
    );

    Run77Data -> Get_l1_alignment();
    Run77Data -> Get_l0l2_slope();
    Run77Data -> GetDetectionEffi(
        Effi_slope_cut, // note : slope cut
        Effi_pos_cut, // note : noise hit distance
        Effi_boundary_cut, // note : boundary cut
        Effi_cluster_adc_value_requirement // note : cluster adc value requirement
    );

    if (prepare_ClusDist) {
        Run77Data -> GetControlClusDist(
            ClusDist_slope_cut, // note : slope cut
            ClusDist_pos_cut // note : noise hit distance
        );
    }

    Run77Data -> EndRun();

    return {
        Run77Data -> GetOutputFileName(),
        Run77Data -> Get_Final_Effi()[0],
        Run77Data -> Get_Final_Effi()[1],
        Run77Data -> Get_Final_Effi()[2]
    };
}



int Run77() {

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
    std::tuple<std::string, double, double, double> tuple_BaseLine = Run77_mother(
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