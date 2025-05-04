#include "DetectionEffiAna.h"

DetectionEffiAna::DetectionEffiAna(
    bool isData_in,
    int runnumber_in,
    std::string input_directory_in,
    std::string input_filename_in,
    std::string output_directory_in,
    std::string output_file_name_suffix_in,

    int Selected_Column_in, // note : 1 to 13;
    double l1_alignment_correction_in,
    double l0l2_slope_correction_in
) : 
    isData(isData_in),
    runnumber(runnumber_in),
    input_directory(input_directory_in),
    input_filename(input_filename_in),
    output_directory(output_directory_in),
    output_file_name_suffix(output_file_name_suffix_in),

    Selected_Column(Selected_Column_in),
    l1_alignment_correction(l1_alignment_correction_in),
    l0l2_slope_correction(l0l2_slope_correction_in)
{
    if (Selected_Column < 1 || Selected_Column > 13) {
        std::cerr << "Error: Selected_Column must be between 1 and 13." << std::endl;
        exit(1);
    }

    system(Form("mkdir -p %s", output_directory.c_str()) );

    cluster_info_NoL1Aligned_vec.clear();

    PrepareInputData();
    PrepareOutputFileName();
    PrepareOutputFile();
    PrepareHistFit();
}

void DetectionEffiAna::PrepareInputData() {

    file_in = TFile::Open(Form("%s/%s", input_directory.c_str(), input_filename.c_str() ));
    if (!file_in || file_in->IsZombie()) {
        std::cerr << "Error: Unable to open input file: " << input_filename << std::endl;
        exit(1);
    }

    tree_in = (TTree*)file_in->Get("cluster_info");
    if (!tree_in) {
        std::cerr << "Error: Unable to find tree 'cluster_info' in file: " << input_filename << std::endl;
        exit(1);
    }

    std::cout<<"Number of events: "<<tree_in->GetEntries()<<std::endl;

    layer = 0;
    chip = 0;
    Nhit = 0;
    position = 0;
    cluster_adc = 0;

    tree_in -> SetBranchAddress("eID", &eID);
    if (isData) {tree_in -> SetBranchAddress("DSE", &DSE);}
    tree_in -> SetBranchAddress("layer", &layer);
    tree_in -> SetBranchAddress("chip", &chip);
    tree_in -> SetBranchAddress("Nhit", &Nhit);
    tree_in -> SetBranchAddress("position", &position);
    tree_in -> SetBranchAddress("cluster_adc", &cluster_adc);


    // std::vector<std::vector<std::vector<std::vector<ClusInfo>>>> cluster_info_NoL1Aligned_vec; // note : [event][chip][layer][clusters]

    // note : [layer][clusters]
    std::vector<std::vector<ClusInfo>> empty_SingleEvt_ClusInfo_vec = std::vector<std::vector<ClusInfo>>(3); 
    empty_SingleEvt_ClusInfo_vec[0].clear();
    empty_SingleEvt_ClusInfo_vec[1].clear();
    empty_SingleEvt_ClusInfo_vec[2].clear();
    std::cout<<"empty_SingleEvt_ClusInfo_vec.size() : "<<empty_SingleEvt_ClusInfo_vec.size()<<std::endl;
    std::cout<<"empty_SingleEvt_ClusInfo_vec[0].size() : "<<empty_SingleEvt_ClusInfo_vec[0].size()<<std::endl;
    std::cout<<"empty_SingleEvt_ClusInfo_vec[1].size() : "<<empty_SingleEvt_ClusInfo_vec[1].size()<<std::endl;
    std::cout<<"empty_SingleEvt_ClusInfo_vec[2].size() : "<<empty_SingleEvt_ClusInfo_vec[2].size()<<std::endl;

    // note : [chip][layer][clusters]
    std::vector<std::vector<std::vector<ClusInfo>>> empty_SingleEvt_ClusInfo_chip_vec = std::vector<std::vector<std::vector<ClusInfo>>>(13,empty_SingleEvt_ClusInfo_vec); 
    std::cout<<"empty_SingleEvt_ClusInfo_chip_vec.size() : "<<empty_SingleEvt_ClusInfo_chip_vec.size()<<std::endl;
    std::cout<<"empty_SingleEvt_ClusInfo_chip_vec[0].size() : "<<empty_SingleEvt_ClusInfo_chip_vec[0].size()<<std::endl;

    for (int i = 0; i < tree_in->GetEntries(); i++) {
        tree_in->GetEntry(i);

        if (i%5000 == 0){
            std::cout<<"Rading data, Event : "<<i<<std::endl;
        }

        if (layer->size() != chip->size() || layer->size() != Nhit->size() || layer->size() != position->size() || layer->size() != cluster_adc->size()) {
            std::cerr << "Error: Mismatched vector sizes in event " << i <<", eID: "<<eID<< std::endl;
            std::cerr << "layer size: " << layer->size() << ", chip size: " << chip->size() << ", Nhit size: " << Nhit->size() << ", position size: " << position->size() << ", cluster_adc size: " << cluster_adc->size() << std::endl;
            exit(1);
        }
        
        if (isData && DSE != 0) {continue;} // todo: DSEs are removed from the analysis
        if (layer->size() == 0) {continue;}

        cluster_info_NoL1Aligned_vec.push_back(empty_SingleEvt_ClusInfo_chip_vec);

        for (int clu_i = 0; clu_i < layer->size(); clu_i++) {

            ClusInfo clus_info;
            clus_info.adc = cluster_adc->at(clu_i);
            clus_info.size = Nhit->at(clu_i);
            clus_info.pos = position->at(clu_i);

            cluster_info_NoL1Aligned_vec.back()[ chip->at(clu_i) - 1 ][ layer->at(clu_i) ].push_back(clus_info);
        }
    }
    
}

void DetectionEffiAna::PrepareOutputFileName() {
    
    if (output_file_name_suffix.size() > 0 && output_file_name_suffix[0] != '_') {
        output_file_name_suffix = "_" + output_file_name_suffix;
    }

    output_filename = (isData) ? "Data" : "MC";
    output_filename += "_DetectionEffi";
    output_filename += Form("_Run%d", runnumber);
    output_filename += Form("_Column%d", Selected_Column);

    output_filename += output_file_name_suffix;

    output_filename += ".root";
}

void DetectionEffiAna::PrepareOutputFile() {

    file_out = new TFile(Form("%s/%s", output_directory.c_str(), output_filename.c_str()), "RECREATE");

    tree_out = new TTree("tree","tree");

    tree_out -> Branch("l1_alignment_correction", &l1_alignment_correction);
    tree_out -> Branch("l0l2_slope_correction", &l0l2_slope_correction);

    tree_out -> Branch("slope_cut", &slope_cut);
    tree_out -> Branch("noise_hit_distance", &noise_hit_distance);
    tree_out -> Branch("boundary_cut", &boundary_cut);
    tree_out -> Branch("cluster_adc_value_requirement", &cluster_adc_value_requirement);


    tree_out_effi = new TTree("tree_effi","tree_effi");
    tree_out_effi -> Branch("L0L2Interpolation", &L0L2Interpolation);
    tree_out_effi -> Branch("L1Good", &L1Good);
}

void DetectionEffiAna::PrepareHistFit() {

    // h1D_l1_alignment_before = new TH1D("h1D_l1_alignment_before","h1D_l1_alignment_before;L1 - (L2L0 interpolation) [mm];Entries",50,-0.4,1);
    // h1D_l1_alignment_after = new TH1D("h1D_l1_alignment_after","h1D_l1_alignment_after;L1 - (L2L0 interpolation) [mm];Entries",50,-0.4,1);
    
    // h1D_l0l1_slope_before = new TH1D("h1D_l0l1_slope_before","h1D_l0l1_slope_before;Slope (L1 - L0);Entries",50,-0.05,0.05);
    // h1D_l0l1_slope_after = new TH1D("h1D_l0l1_slope_after","h1D_l0l1_slope_after;Slope (L1 - L0);Entries",50,-0.05,0.05);

    fit_l1_alignment_before = new TF1("fit_l1_alignment_before","gaus",-5,5);
    fit_l1_alignment_before -> SetNpx(1000);

    fit_l1_alignment_after = new TF1("fit_l1_alignment_after","gaus",-5,5);
    fit_l1_alignment_after -> SetNpx(1000);

    fit_l0l2_slope_before = new TF1("fit_l0l2_slope_before","gaus",-5,5);
    fit_l0l2_slope_before -> SetNpx(1000);

    fit_l0l2_slope_after = new TF1("fit_l0l2_slope_after","gaus",-5,5);
    fit_l0l2_slope_after -> SetNpx(1000);


    h1D_effi_l1_residual = new TH1D("h1D_effi_l1_residual","h1D_effi_l1_residual;L1 - (L2L0 interpolation) [mm];Entries",50,-1,1);
    h1D_effi_l0l2_slope = new TH1D("h1D_effi_l0l2_slope","h1D_effi_l0l2_slope;Slope (L2 - L0);Entries",50,-0.05,0.05);


    // Division : --------------------------------------------------------------------------------
    // note : for other study

    h1D_ClusPos_Raw[0] = new TH1D("h1D_l0_ClusPos_Raw", "h1D_l0_ClusPos_Raw;ClusPos (L0) [mm];Entries",128,-9.984,9.984);
    h1D_ClusPos_Raw[1] = new TH1D("h1D_l1_ClusPos_Raw", "h1D_l1_ClusPos_Raw;ClusPos (L1) [mm];Entries",128,-9.984,9.984);
    h1D_ClusPos_Raw[2] = new TH1D("h1D_l2_ClusPos_Raw", "h1D_l2_ClusPos_Raw;ClusPos (L2) [mm];Entries",128,-9.984,9.984);

    h1D_ClusPos_GoodEvt[0] = new TH1D("h1D_l0_ClusPos_GoodEvt", "h1D_l0_ClusPos_GoodEvt;ClusPos (L0) [mm];Entries",128,-9.984,9.984);
    h1D_ClusPos_GoodEvt[1] = new TH1D("h1D_l1_ClusPos_GoodEvt", "h1D_l1_ClusPos_GoodEvt;ClusPos (L1) [mm];Entries",128,-9.984,9.984);
    h1D_ClusPos_GoodEvt[2] = new TH1D("h1D_l2_ClusPos_GoodEvt", "h1D_l2_ClusPos_GoodEvt;ClusPos (L2) [mm];Entries",128,-9.984,9.984);
    
    h1D_ClusPos_3Hit[0] = new TH1D("h1D_l0_ClusPos_3Hit", "h1D_l0_ClusPos_3Hit;ClusPos (L0) [mm];Entries",128,-9.984,9.984);
    h1D_ClusPos_3Hit[1] = new TH1D("h1D_l1_ClusPos_3Hit", "h1D_l1_ClusPos_3Hit;ClusPos (L1) [mm];Entries",128,-9.984,9.984);
    h1D_ClusPos_3Hit[2] = new TH1D("h1D_l2_ClusPos_3Hit", "h1D_l2_ClusPos_3Hit;ClusPos (L2) [mm];Entries",128,-9.984,9.984);
    
    h1D_ClusAdc[0] = new TH1D("h1D_l0_ClusAdc", "h1D_l0_ClusAdc;ClusAdc (L0);Entries",100,0,2000);
    h1D_ClusAdc[1] = new TH1D("h1D_l1_ClusAdc", "h1D_l1_ClusAdc;ClusAdc (L1);Entries",100,0,2000);
    h1D_ClusAdc[2] = new TH1D("h1D_l2_ClusAdc", "h1D_l2_ClusAdc;ClusAdc (L2);Entries",100,0,2000);
    
    h1D_ClusSize[0] = new TH1D("h1D_l0_ClusSize", "h1D_l0_ClusSize;ClusSize (L0);Entries",40,0,40);
    h1D_ClusSize[1] = new TH1D("h1D_l1_ClusSize", "h1D_l1_ClusSize;ClusSize (L1);Entries",40,0,40);
    h1D_ClusSize[2] = new TH1D("h1D_l2_ClusSize", "h1D_l2_ClusSize;ClusSize (L2);Entries",40,0,40);

}

double DetectionEffiAna::Get_l1_alignment() {
    if(l1_alignment_correction == 0.) {
        std::cout<<"l1_alignment_correction is 0, then measure it!"<<std::endl;

        h1D_l1_alignment_before = (TH1D *) (Sub_Get_l1_alignment_hist()) -> Clone("h1D_l1_alignment_before");
        h1D_l1_alignment_before -> SetTitle("h1D_l1_alignment_before");
        h1D_l1_alignment_before -> Fit(fit_l1_alignment_before,"N");

        l1_alignment_correction = -1 * fit_l1_alignment_before -> GetParameter(1);
        std::cout<<"l1_alignment_correction : "<<l1_alignment_correction<<std::endl;
    }
    else{
        std::cout<<"l1_alignment_correction is already set to "<<l1_alignment_correction<<std::endl;
    }

    h1D_l1_alignment_after = (TH1D *) (Sub_Get_l1_alignment_hist()) -> Clone("h1D_l1_alignment_after");
    h1D_l1_alignment_after -> SetTitle("h1D_l1_alignment_after");
    h1D_l1_alignment_after -> Fit(fit_l1_alignment_after,"N");
    std::cout<<"post correction, l1_alignment, new peak: "<<fit_l1_alignment_after -> GetParameter(1)<<std::endl;

    return l1_alignment_correction;
}

TH1D * DetectionEffiAna::Sub_Get_l1_alignment_hist(){
    TH1D * h1D_l1_alignment_temp = new TH1D("",";L1 - (L2L0 interpolation) [mm];Entries",50,-0.4,1);
    h1D_l1_alignment_temp -> Reset("ICESM");

    // double edge_exclusion_bottom = (lower_section_initial - INTT_strip_width / 2.) + INTT_strip_width * double(boundary_cut);
	// double edge_exclusion_upper = ( INTT_strip_width * 128. ) - INTT_strip_width * double(boundary_cut);

    // note : [event][chip][layer][clusters]
    for (int i = 0; i < cluster_info_NoL1Aligned_vec.size(); i++) {
        if (i % 5000 == 0){
            std::cout<<"Doing l1_alignment, Event : "<<i<<std::endl;
        }

        std::vector<ClusInfo> empty_ClusInfo_vec; empty_ClusInfo_vec.clear();

        std::vector<ClusInfo> CheckedCol_L0_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][0];
        std::vector<ClusInfo> CheckedCol_L1_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][1];
        std::vector<ClusInfo> CheckedCol_L2_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][2];

        std::vector<ClusInfo> CheckedColM1_L0_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L1_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L2_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][2] : empty_ClusInfo_vec;

        std::vector<ClusInfo> CheckedColP1_L0_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L1_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L2_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][2] : empty_ClusInfo_vec;

        // note : zero cluster in adjacent chips
        if ( 
            (
                CheckedColM1_L0_vec.size() + 
                CheckedColM1_L1_vec.size() + 
                CheckedColM1_L2_vec.size() + 
                
                CheckedColP1_L0_vec.size() + 
                CheckedColP1_L1_vec.size() + 
                CheckedColP1_L2_vec.size()
            ) != 0 
        ) {continue;}
        

        // note : single cluster in l0 and l2
        if (CheckedCol_L0_vec.size() != 1 || CheckedCol_L2_vec.size() != 1) {continue;}

        // note : single cluster in l1
        if (CheckedCol_L1_vec.size() != 1) {continue;}

        double expected_middle_pos = (CheckedCol_L0_vec[0].pos + CheckedCol_L2_vec[0].pos) / 2.;
        double diff = CheckedCol_L1_vec[0].pos - expected_middle_pos + l1_alignment_correction;

        h1D_l1_alignment_temp -> Fill(diff);

        // note : make sure middle layer has hits
        // todo : new method to find the closest cluster in the middle layer
        // todo : if I do the loop, and find the one closest to the expected position, I can be biased. Since the misalignment is still there
        // if (CheckedCol_L1_vec.size() != 0) {
        //     double expected_middle_pos = (CheckedCol_L0_vec[0] + CheckedCol_L2_vec[0]) / 2.;

        //     for (int clu_i = 0; clu_i < CheckedCol_L1_vec.size(); clu_i++){
        //         if (clu_i == 0) {
        //             least_diff = CheckedCol_L1_vec[clu_i].pos - expected_middle_pos;
        //         }
        //         else if (std::abs(least_diff) < std::abs(CheckedCol_L1_vec[clu_i].pos - expected_middle_pos)) {
        //             least_diff = CheckedCol_L1_vec[clu_i].pos - expected_middle_pos;
        //         }
        //     }
        // }
    
    } // note : loop of cluster_info_NoL1Aligned_vec

    return h1D_l1_alignment_temp;
}

double DetectionEffiAna::Get_l0l2_slope(){
    if(l0l2_slope_correction == 0.) {
        std::cout<<"l0l2_slope_correction is 0, then measure it!"<<std::endl;

        h1D_l0l2_slope_before = (TH1D *) (Sub_Get_l0l2_slope_hist()) -> Clone("h1D_l0l2_slope_before");
        h1D_l0l2_slope_before -> SetTitle("h1D_l0l2_slope_before");
        h1D_l0l2_slope_before -> Fit(fit_l0l2_slope_before,"N");

        l0l2_slope_correction = -1 * fit_l0l2_slope_before -> GetParameter(1);
        std::cout<<"l0l2_slope_correction : "<<l0l2_slope_correction<<std::endl;
    }
    else{
        std::cout<<"l0l2_slope_correction is already set to "<<l0l2_slope_correction<<std::endl;
    }

    h1D_l0l2_slope_after = (TH1D *) (Sub_Get_l0l2_slope_hist()) -> Clone("h1D_l0l2_slope_after");
    h1D_l0l2_slope_after -> SetTitle("h1D_l0l2_slope_after");
    h1D_l0l2_slope_after -> Fit(fit_l0l2_slope_after,"N");
    std::cout<<"post correction, l0l2 slope, new peak: "<<fit_l0l2_slope_after -> GetParameter(1)<<std::endl;

    return l0l2_slope_correction;
}

TH1D* DetectionEffiAna::Sub_Get_l0l2_slope_hist() {
    TH1D * h1D_l0l2_slope_temp = new TH1D("",";L2L0 slope;Entries",50,-0.05,0.05);
    h1D_l0l2_slope_temp -> Reset("ICESM");

    // double edge_exclusion_bottom = (lower_section_initial - INTT_strip_width / 2.) + INTT_strip_width * double(boundary_cut);
	// double edge_exclusion_upper = ( INTT_strip_width * 128. ) - INTT_strip_width * double(boundary_cut);

    // note : [event][chip][layer][clusters]
    for (int i = 0; i < cluster_info_NoL1Aligned_vec.size(); i++) {
        if (i % 5000 == 0){
            std::cout<<"Doing l0l2_slope, Event : "<<i<<std::endl;
        }

        std::vector<ClusInfo> empty_ClusInfo_vec; empty_ClusInfo_vec.clear();

        std::vector<ClusInfo> CheckedCol_L0_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][0];
        std::vector<ClusInfo> CheckedCol_L1_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][1];
        std::vector<ClusInfo> CheckedCol_L2_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][2];

        std::vector<ClusInfo> CheckedColM1_L0_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L1_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L2_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][2] : empty_ClusInfo_vec;

        std::vector<ClusInfo> CheckedColP1_L0_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L1_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L2_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][2] : empty_ClusInfo_vec;

        // note : zero cluster in adjacent chips
        if ( 
            (
                CheckedColM1_L0_vec.size() + 
                CheckedColM1_L1_vec.size() + 
                CheckedColM1_L2_vec.size() + 
                
                CheckedColP1_L0_vec.size() + 
                CheckedColP1_L1_vec.size() + 
                CheckedColP1_L2_vec.size()
            ) != 0 
        ) {continue;}
        

        // note : single cluster in l0 and l2
        if (CheckedCol_L0_vec.size() != 1 || CheckedCol_L2_vec.size() != 1) {continue;}

        // todo: this is masked, to be in line with what we did before
        // note : single cluster in l1
        // if (CheckedCol_L1_vec.size() != 1) {continue;}

        double l0l2_slope = (CheckedCol_L2_vec[0].pos - CheckedCol_L0_vec[0].pos) / actual_xpos[2] + l0l2_slope_correction;

        h1D_l0l2_slope_temp -> Fill(l0l2_slope);
    
    } // note : loop of cluster_info_NoL1Aligned_vec

    return h1D_l0l2_slope_temp;
}

void DetectionEffiAna::GetDetectionEffi(
    double slope_cut_in, 
    double noise_hit_distance_in, // note : unit : mm
    int boundary_cut_in, 
    double cluster_adc_value_requirement_in
){
    slope_cut = slope_cut_in;
    noise_hit_distance = noise_hit_distance_in;
    boundary_cut = boundary_cut_in;
    cluster_adc_value_requirement = cluster_adc_value_requirement_in;

    double edge_exclusion_bottom = (lower_section_initial - INTT_strip_width / 2.) + INTT_strip_width * double(boundary_cut);
	double edge_exclusion_upper = ( INTT_strip_width * 128. ) - INTT_strip_width * double(boundary_cut);

    int count_1 = 0;
    int count_2 = 0;
    int count_3 = 0;
    int count_4 = 0;
    int count_5 = 0;
    int count_6 = 0;
    int count_7 = 0;
    int count_8 = 0;

    for (int i = 0; i < cluster_info_NoL1Aligned_vec.size(); i++) {
        if (i % 5000 == 0){
            std::cout<<"Doing Detection Efficiency, Event : "<<i<<std::endl;
        }

        L1Good = 0;

        std::vector<ClusInfo> empty_ClusInfo_vec; empty_ClusInfo_vec.clear();

        std::vector<ClusInfo> CheckedCol_L0_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][0];
        std::vector<ClusInfo> CheckedCol_L1_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][1];
        std::vector<ClusInfo> CheckedCol_L2_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][2];

        std::vector<ClusInfo> CheckedColM1_L0_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L1_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L2_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][2] : empty_ClusInfo_vec;

        std::vector<ClusInfo> CheckedColP1_L0_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L1_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L2_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][2] : empty_ClusInfo_vec;

        count_1++;

        // note : zero cluster in adjacent chips
        if ( 
            (
                CheckedColM1_L0_vec.size() + 
                CheckedColM1_L1_vec.size() + 
                CheckedColM1_L2_vec.size() + 
                
                CheckedColP1_L0_vec.size() + 
                CheckedColP1_L1_vec.size() + 
                CheckedColP1_L2_vec.size()
            ) != 0 
        ) {continue;}
        
        count_2++;

        // note : single cluster in l0 and l2
        if (CheckedCol_L0_vec.size() != 1 || CheckedCol_L2_vec.size() != 1) {continue;}

        count_3++;

        // note : standalone cluster adc cut
        if ( CheckedCol_L0_vec[0].adc <= cluster_adc_value_requirement || CheckedCol_L2_vec[0].adc <= cluster_adc_value_requirement ) {continue;}

        count_4++;


        // note : edge exclusion cut of the l0
        if ( CheckedCol_L0_vec[0].pos <= edge_exclusion_bottom || CheckedCol_L0_vec[0].pos >= edge_exclusion_upper ) {continue;}

        count_5++;

        // note : edge exclusion cut of the l2
        if ( CheckedCol_L2_vec[0].pos <= edge_exclusion_bottom || CheckedCol_L2_vec[0].pos >= edge_exclusion_upper ) {continue;}

        count_6++;

        double l0l2_slope = (CheckedCol_L2_vec[0].pos - CheckedCol_L0_vec[0].pos) / actual_xpos[2] + l0l2_slope_correction;
        h1D_effi_l0l2_slope -> Fill(l0l2_slope);

        // note : the slope cut
        if ( fabs(l0l2_slope) >= slope_cut ) {continue;}

        count_7++; // note : denominator

        L0L2Interpolation = ( CheckedCol_L0_vec[0].pos + CheckedCol_L2_vec[0].pos ) / 2.;

        if (CheckedCol_L1_vec.size() != 0) {
            
            double diff;

            for (int clu_i = 0; clu_i < CheckedCol_L1_vec.size(); clu_i++){
                if (clu_i == 0) {
                    diff = (CheckedCol_L1_vec[clu_i].pos + l1_alignment_correction) - L0L2Interpolation;
                }
                else if (std::abs(diff) > std::abs((CheckedCol_L1_vec[clu_i].pos + l1_alignment_correction) - L0L2Interpolation)) {
                    diff = (CheckedCol_L1_vec[clu_i].pos + l1_alignment_correction) - L0L2Interpolation;
                }
            }

            h1D_effi_l1_residual -> Fill(diff);

            // note : the distance cut
            if ( fabs(diff) >= noise_hit_distance ) {
                L1Good = 0;
            }
            else {
                count_8++;
                L1Good = 1;
            }

        }
        else {
            L1Good = 0;
        }

        tree_out_effi -> Fill();
    }    

    std::cout<<"cluster_info_NoL1Aligned_vec.size() : "<<count_1<<std::endl;
    std::cout<<"zero cluster in adjacent chips      : "<<count_2<<std::endl;
    std::cout<<"single cluster in l0 and l2         : "<<count_3<<std::endl;
    std::cout<<"standalone cluster adc cut          : "<<count_4<<std::endl;
    std::cout<<"edge exclusion cut of the l0        : "<<count_5<<std::endl;
    std::cout<<"edge exclusion cut of the l2        : "<<count_6<<std::endl;
    std::cout<<"the slope cut                       : "<<count_7<<", the denominator"<<std::endl;
    std::cout<<"L1Good                              : "<<count_8<<", the numerator"<<std::endl;
    std::cout<<"Efficiency                          : "<<double(count_8) / double(count_7)<<std::endl;


    TH1D * total_hist = new TH1D ("","",1,0,1);
    TH1D * pass_hist = new TH1D ("","",1,0,1);
    
    // note : layer - 1 efficiency
    total_hist->SetBinContent(1,count_7); 
    pass_hist ->SetBinContent(1,count_8);

    TEfficiency * detection_effi = new TEfficiency (*pass_hist,*total_hist);
    printf("Efficiency by TEfficiency \n");
    printf("%.5f\t+%.5f\t-%.5f \n", detection_effi->GetEfficiency(1)*100.,detection_effi->GetEfficiencyErrorUp(1)*100.,detection_effi->GetEfficiencyErrorLow(1)*100.);

    Final_Effi = detection_effi->GetEfficiency(1)*100.;
    Final_Effi_StatErrorUp = detection_effi->GetEfficiencyErrorUp(1)*100.;
    Final_Effi_StatErrorDown = detection_effi->GetEfficiencyErrorLow(1)*100.;
}

void DetectionEffiAna:: GetControlClusDist(
    double slope_cut_in, 
    double noise_hit_distance_in // note : unit : mm
){

    for (int i = 0; i < cluster_info_NoL1Aligned_vec.size(); i++) {
        if (i % 50000 == 0){
            std::cout<<"Doing Detection Efficiency, Event : "<<i<<std::endl;
        }

        std::vector<ClusInfo> empty_ClusInfo_vec; empty_ClusInfo_vec.clear();

        std::vector<ClusInfo> CheckedCol_L0_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][0];
        std::vector<ClusInfo> CheckedCol_L1_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][1];
        std::vector<ClusInfo> CheckedCol_L2_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][2];

        std::vector<ClusInfo> CheckedColM1_L0_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L1_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L2_vec = (Selected_Column-1-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1-1][2] : empty_ClusInfo_vec;

        std::vector<ClusInfo> CheckedColP1_L0_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L1_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L2_vec = (Selected_Column-1+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1+1][2] : empty_ClusInfo_vec;

        for (int li = 0; li < 3; li++){
            for (int clu_i = 0; clu_i < cluster_info_NoL1Aligned_vec[i][Selected_Column-1][li].size(); clu_i++){
                double this_clus_pos = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][li][clu_i].pos;
                this_clus_pos = (li == 1) ? (this_clus_pos + l1_alignment_correction) : this_clus_pos;
                h1D_ClusPos_Raw[li] -> Fill(this_clus_pos);
            }
        }

        // note : zero cluster in adjacent chips
        if ( 
            (
                CheckedColM1_L0_vec.size() + 
                CheckedColM1_L1_vec.size() + 
                CheckedColM1_L2_vec.size() + 
                
                CheckedColP1_L0_vec.size() + 
                CheckedColP1_L1_vec.size() + 
                CheckedColP1_L2_vec.size()
            ) != 0 
        ) {continue;}
        

        // note : single cluster in l0 and l2
        if (CheckedCol_L0_vec.size() != 1 || CheckedCol_L2_vec.size() != 1) {continue;}


        double this_l0l2_slope = (CheckedCol_L2_vec[0].pos - CheckedCol_L0_vec[0].pos) / actual_xpos[2] + l0l2_slope_correction;

        // note : the slope cut
        if ( fabs(this_l0l2_slope) >= slope_cut_in ) {continue;}


        for (int li = 0; li < 3; li++){
            for (int clu_i = 0; clu_i < cluster_info_NoL1Aligned_vec[i][Selected_Column-1][li].size(); clu_i++){
                double this_clus_pos = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][li][clu_i].pos;
                this_clus_pos = (li == 1) ? (this_clus_pos + l1_alignment_correction) : this_clus_pos;
                h1D_ClusPos_GoodEvt[li] -> Fill(this_clus_pos);
            }
        }

        if (CheckedCol_L1_vec.size() != 1) {continue;}

        double diff = (CheckedCol_L1_vec[0].pos + l1_alignment_correction) - (CheckedCol_L0_vec[0].pos + CheckedCol_L2_vec[0].pos) / 2.;

        if (std::abs(diff) < noise_hit_distance_in) {
            for (int li = 0; li < 3; li++){
                for (int clu_i = 0; clu_i < cluster_info_NoL1Aligned_vec[i][Selected_Column-1][li].size(); clu_i++){
                    
                    double this_clus_pos = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][li][clu_i].pos;
                    this_clus_pos = (li == 1) ? (this_clus_pos + l1_alignment_correction) : this_clus_pos;
                    h1D_ClusPos_3Hit[li] -> Fill(this_clus_pos);
                
                    h1D_ClusAdc[li] -> Fill(cluster_info_NoL1Aligned_vec[i][Selected_Column-1][li][clu_i].adc);
                    h1D_ClusSize[li] -> Fill(cluster_info_NoL1Aligned_vec[i][Selected_Column-1][li][clu_i].size);
                }
            }
        }

    }
}


void DetectionEffiAna::EndRun(){
    file_out -> cd();
    if (h1D_l1_alignment_before != nullptr) {h1D_l1_alignment_before -> Write();}
    if (h1D_l1_alignment_after != nullptr) {h1D_l1_alignment_after -> Write();}
    if (h1D_l0l2_slope_before != nullptr) {h1D_l0l2_slope_before -> Write();}
    if (h1D_l0l2_slope_after != nullptr) {h1D_l0l2_slope_after -> Write();}

    if (h1D_effi_l1_residual != nullptr) {h1D_effi_l1_residual -> Write();}
    if (h1D_effi_l0l2_slope != nullptr) {h1D_effi_l0l2_slope -> Write();}

    tree_out -> Fill();
    tree_out -> Write();

    tree_out_effi -> Write();

    TCanvas * c1 = new TCanvas("c1","c1",800,600);
    c1 -> cd();

    if (h1D_l1_alignment_before != nullptr) {
        c1 -> cd();
        h1D_l1_alignment_before -> Draw();
        fit_l1_alignment_before -> Draw("lsame");
        c1 -> Write("c1_h1D_l1_alignment_before");
    }

    if (h1D_l1_alignment_after != nullptr) {
        c1 -> cd();
        h1D_l1_alignment_after -> Draw();
        fit_l1_alignment_after -> Draw("lsame");
        c1 -> Write("c1_h1D_l1_alignment_after");
    }

    if (h1D_l0l2_slope_before != nullptr) {
        c1 -> cd();
        h1D_l0l2_slope_before -> Draw();
        fit_l0l2_slope_before -> Draw("lsame");
        c1 -> Write("c1_h1D_l0l2_slope_before");
    }

    if (h1D_l0l2_slope_after != nullptr) {
        c1 -> cd();
        h1D_l0l2_slope_after -> Draw();
        fit_l0l2_slope_after -> Draw("lsame");
        c1 -> Write("c1_h1D_l0l2_slope_after");
    }

    for (int i = 0; i < 3; i++){
        h1D_ClusPos_Raw[i] -> Write();
    }

    for (int i = 0; i < 3; i++){
        h1D_ClusPos_GoodEvt[i] -> Write();
    }

    for (int i = 0; i < 3; i++){
        h1D_ClusPos_3Hit[i] -> Write();
    }

    for (int i = 0; i < 3; i++){
        h1D_ClusAdc[i] -> Write();
    }

    for (int i = 0; i < 3; i++){
        h1D_ClusSize[i] -> Write();
    }

    file_out -> Close();
    file_in -> Close();

    std::cout<<"EndRun, "<<output_filename<<std::endl;
}