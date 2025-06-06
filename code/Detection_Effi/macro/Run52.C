#include "../DetectionEffiAna.h"
R__LOAD_LIBRARY(../libDetectionEffiAna.so)

// int Run52_Original(){
//     std::string data_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52";
//     std::string data_input_filename = "run52_no_clone_filter_all_clusters.root";
//     std::string data_output_directory = data_input_directory + "/DetectionEffi";

//     std::string data_output_file_name_suffix = "test1";    
//     int Selected_Column = 8; // note : 1 to 13;
//     double Effi_slope_cut = 0.01; // note : slope cut
//     double Effi_pos_cut = 0.4; // note : noise hit distance
//     int    Effi_boundary_cut = 5; // note : boundary cut
//     double Effi_cluster_adc_value_requirement = 15; // note : cluster adc value requirement

//     double ClusDist_slope_cut = 0.01; // note : slope cut
//     double ClusDist_pos_cut = 0.234; // note : noise hit distance


//     DetectionEffiAna * Run52Data = new DetectionEffiAna(
//         true,
//         52,
//         data_input_directory,
//         data_input_filename,
//         data_output_directory,
//         data_output_file_name_suffix,

//         Selected_Column // Selected_Column_in, // note : 1 to 13;
//          // l1_alignment_correction_in,
//          // l0l2_slope_correction_in
//     );

//     Run52Data -> Get_l1_alignment();
//     Run52Data -> Get_l0l2_slope();
//     Run52Data -> GetDetectionEffi(
//         Effi_slope_cut, // note : slope cut
//         Effi_pos_cut, // note : noise hit distance
//         Effi_boundary_cut, // note : boundary cut
//         Effi_cluster_adc_value_requirement // note : cluster adc value requirement
//     );
//     Run52Data -> GetControlClusDist(
//         ClusDist_slope_cut, // note : slope cut
//         ClusDist_pos_cut // note : noise hit distance
//     );

//     Run52Data -> EndRun();

//     return 888;
// }

std::string data_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52";
std::string data_input_filename = "run52_no_clone_filter_all_clusters.root";
std::string data_output_directory = data_input_directory + "/DetectionEffi_baseline700um";

std::tuple<std::string, double, double, double> Run52_mother(
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

    DetectionEffiAna * Run52Data = new DetectionEffiAna(
        true,
        52,
        data_input_directory,
        data_input_filename,
        data_output_directory,
        data_output_file_name_suffix,

        Selected_Column // Selected_Column_in, // note : 1 to 13;
         // l1_alignment_correction_in,
         // l0l2_slope_correction_in
    );

    Run52Data -> Get_l1_alignment();
    Run52Data -> Get_l0l2_slope();
    Run52Data -> GetDetectionEffi(
        Effi_slope_cut, // note : slope cut
        Effi_pos_cut, // note : noise hit distance
        Effi_boundary_cut, // note : boundary cut
        Effi_cluster_adc_value_requirement // note : cluster adc value requirement
    );

    if (prepare_ClusDist) {
        Run52Data -> GetControlClusDist(
            ClusDist_slope_cut, // note : slope cut
            ClusDist_pos_cut // note : noise hit distance
        );
    }

    Run52Data -> EndRun();

    return {
        Run52Data -> GetOutputFileName(),
        Run52Data -> Get_Final_Effi()[0],
        Run52Data -> Get_Final_Effi()[1],
        Run52Data -> Get_Final_Effi()[2]
    };
}



int Run52() {

    string BaseLine_data_output_file_name_suffix = "BaseLine";    
    int    BaseLine_Selected_Column = 8; // note : 1 to 13,
    double BaseLine_Effi_slope_cut = 0.01; // note : slope cut
    double BaseLine_Effi_pos_cut = 0.7; // note : noise hit distance
    int    BaseLine_Effi_boundary_cut = 8;// note : boundary cut
    double BaseLine_Effi_cluster_adc_value_requirement = 15; // note : cluster adc value requirement

    double BaseLine_ClusDist_slope_cut = 0.01; // note : slope cut
    double BaseLine_ClusDist_pos_cut = 0.7; // note : noise hit distance

    bool   BaseLine_prepare_ClusDist = true;
    
    // Division: BaseLine ---------------------------------------------------------------------------------
    std::tuple<std::string, double, double, double> tuple_BaseLine = Run52_mother(
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

    // Division: Vary Column ---------------------------------------------------------------------------------
    string VaryColumn_data_output_file_name_suffix = "VaryColumn";
    int    VaryColumn_Selected_Column = 9; // note : 1 to 13,

    std::tuple<std::string, double, double, double> tuple_VaryColumn = Run52_mother(
        VaryColumn_data_output_file_name_suffix, // note : changed
        VaryColumn_Selected_Column,              // note : changed
        BaseLine_Effi_slope_cut,
        BaseLine_Effi_pos_cut,
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );

    // Division: Vary Pos Cut Large ---------------------------------------------------------------------------------
    string VaryPosCutLarge_data_output_file_name_suffix = "VaryPosCutLarge";
    double VaryPosCutLarge_Effi_pos_cut = 1.0;

    std::tuple<std::string, double, double, double> tuple_VaryPosCutLarge = Run52_mother(
        VaryPosCutLarge_data_output_file_name_suffix, // note : changed
        BaseLine_Selected_Column,
        BaseLine_Effi_slope_cut,
        VaryPosCutLarge_Effi_pos_cut, // note : changed
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );

    // Division: Vary Pos Cut Small ---------------------------------------------------------------------------------
    string VaryPosCutSmall_data_output_file_name_suffix = "VaryPosCutSmall";
    double VaryPosCutSmall_Effi_pos_cut = 0.85;

    std::tuple<std::string, double, double, double> tuple_VaryPosCutSmall = Run52_mother(
        VaryPosCutSmall_data_output_file_name_suffix, // note : changed
        BaseLine_Selected_Column,
        BaseLine_Effi_slope_cut,
        VaryPosCutSmall_Effi_pos_cut, // note : changed
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );


    // Division: Vary SlopeCut Large ---------------------------------------------------------------------------------
    string VarySlopeCutLarge_data_output_file_name_suffix = "VarySlopeCutLarge";
    double VarySlopeCutLarge_Effi_slope_cut = 0.013;

    std::tuple<std::string, double, double, double> tuple_VarySlopeCutLarge = Run52_mother(
        VarySlopeCutLarge_data_output_file_name_suffix, // note : changed
        BaseLine_Selected_Column,
        VarySlopeCutLarge_Effi_slope_cut, // note : changed
        BaseLine_Effi_pos_cut,
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );

    // Division: Vary SlopeCut Small ---------------------------------------------------------------------------------
    string VarySlopeCutSmall_data_output_file_name_suffix = "VarySlopeCutSmall";
    double VarySlopeCutSmall_Effi_slope_cut = 0.007;

    std::tuple<std::string, double, double, double> tuple_VarySlopeCutSmall = Run52_mother(
        VarySlopeCutSmall_data_output_file_name_suffix, // note : changed
        BaseLine_Selected_Column,
        VarySlopeCutSmall_Effi_slope_cut, // note : changed
        BaseLine_Effi_pos_cut,
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );

    // Division: Vary Boundary Cut Large ---------------------------------------------------------------------------------
    string VaryBoundaryCutLarge_data_output_file_name_suffix = "VaryBoundaryCutLarge";
    int VaryBoundaryCutLarge_Effi_boundary_cut = 11;

    std::tuple<std::string, double, double, double> tuple_VaryBoundaryCutLarge = Run52_mother(
        VaryBoundaryCutLarge_data_output_file_name_suffix, // note : changed
        BaseLine_Selected_Column,
        BaseLine_Effi_slope_cut,
        BaseLine_Effi_pos_cut,
        VaryBoundaryCutLarge_Effi_boundary_cut, // note : changed
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );

    // Division: Vary Boundary Cut Small ---------------------------------------------------------------------------------
    string VaryBoundaryCutSmall_data_output_file_name_suffix = "VaryBoundaryCutSmall";
    int VaryBoundaryCutSmall_Effi_boundary_cut = 5;

    std::tuple<std::string, double, double, double> tuple_VaryBoundaryCutSmall = Run52_mother(
        VaryBoundaryCutSmall_data_output_file_name_suffix, // note : changed
        BaseLine_Selected_Column,
        BaseLine_Effi_slope_cut,
        BaseLine_Effi_pos_cut,
        VaryBoundaryCutSmall_Effi_boundary_cut, // note : changed
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );


    // Division: Vary SlopeCut Small ---------------------------------------------------------------------------------
    string VarySlopeCutTight_data_output_file_name_suffix = "VarySlopeCutTight";
    double VarySlopeCutTight_Effi_slope_cut = 0.004;

    std::tuple<std::string, double, double, double> tuple_VarySlopeCutTight = Run52_mother(
        VarySlopeCutTight_data_output_file_name_suffix, // note : changed
        BaseLine_Selected_Column,
        VarySlopeCutTight_Effi_slope_cut, // note : changed
        BaseLine_Effi_pos_cut,
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );



    std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<std::get<0>(tuple_BaseLine)<<", "<<std::get<1>(tuple_BaseLine)<<", "<<std::get<2>(tuple_BaseLine)<<", "<<std::get<3>(tuple_BaseLine)<<std::endl;
    std::cout<<std::endl;

    std::cout<<std::get<0>(tuple_VaryColumn)<<", "<<std::get<1>(tuple_VaryColumn)<<", "<<std::get<2>(tuple_VaryColumn)<<", "<<std::get<3>(tuple_VaryColumn)<<std::endl;
    std::cout<<std::endl;
    
    std::cout<<std::get<0>(tuple_VaryPosCutLarge)<<", "<<std::get<1>(tuple_VaryPosCutLarge)<<", "<<std::get<2>(tuple_VaryPosCutLarge)<<", "<<std::get<3>(tuple_VaryPosCutLarge)<<std::endl;
    std::cout<<std::get<0>(tuple_VaryPosCutSmall)<<", "<<std::get<1>(tuple_VaryPosCutSmall)<<", "<<std::get<2>(tuple_VaryPosCutSmall)<<", "<<std::get<3>(tuple_VaryPosCutSmall)<<std::endl;
    std::cout<<std::endl;
    
    std::cout<<std::get<0>(tuple_VarySlopeCutLarge)<<", "<<std::get<1>(tuple_VarySlopeCutLarge)<<", "<<std::get<2>(tuple_VarySlopeCutLarge)<<", "<<std::get<3>(tuple_VarySlopeCutLarge)<<std::endl;
    std::cout<<std::get<0>(tuple_VarySlopeCutSmall)<<", "<<std::get<1>(tuple_VarySlopeCutSmall)<<", "<<std::get<2>(tuple_VarySlopeCutSmall)<<", "<<std::get<3>(tuple_VarySlopeCutSmall)<<std::endl;
    std::cout<<std::endl;
    
    std::cout<<std::get<0>(tuple_VaryBoundaryCutLarge)<<", "<<std::get<1>(tuple_VaryBoundaryCutLarge)<<", "<<std::get<2>(tuple_VaryBoundaryCutLarge)<<", "<<std::get<3>(tuple_VaryBoundaryCutLarge)<<std::endl;
    std::cout<<std::get<0>(tuple_VaryBoundaryCutSmall)<<", "<<std::get<1>(tuple_VaryBoundaryCutSmall)<<", "<<std::get<2>(tuple_VaryBoundaryCutSmall)<<", "<<std::get<3>(tuple_VaryBoundaryCutSmall)<<std::endl;
    std::cout<<std::endl;

    std::cout<<std::get<0>(tuple_VarySlopeCutTight)<<", "<<std::get<1>(tuple_VarySlopeCutTight)<<", "<<std::get<2>(tuple_VarySlopeCutTight)<<", "<<std::get<3>(tuple_VarySlopeCutTight)<<std::endl;

    return 888;
}