#include "ResidualCompAna.h"

ResidualCompAna::ResidualCompAna(
    bool isData_in,
    int runnumber_in,
    std::string input_directory_in,
    std::string input_filename_in,
    std::string output_directory_in,
    std::string output_file_name_suffix_in,

    int Selected_Column_in, // note : 1 to 13;
    double l1_alignment_correction_in,
    double l0l1_slope_correction_in,
    int ClusPhiSize_cut_in,
    int ClusAdc_cut_in
 ) : 
    isData(isData_in),
    runnumber(runnumber_in),
    input_directory(input_directory_in),
    input_filename(input_filename_in),
    output_directory(output_directory_in),
    output_file_name_suffix(output_file_name_suffix_in),
    Selected_Column(Selected_Column_in),
    l1_alignment_correction(l1_alignment_correction_in),
    l0l1_slope_correction(l0l1_slope_correction_in),
    ClusPhiSize_cut(ClusPhiSize_cut_in),
    ClusAdc_cut(ClusAdc_cut_in),

    slope_cut(std::nan("")),
    pos_cut(std::nan(""))
{
    if (Selected_Column < 1 || Selected_Column > 13) {
        std::cerr << "Error: Selected_Column must be between 1 and 13." << std::endl;
        exit(1);
    }

    system(Form("mkdir -p %s", output_directory.c_str()) );

    cluster_info_NoL1Aligned_vec.clear();
    // cluster_info_PostL1Aligned_vec.clear();   

    PrepareInputData();
    PrepareOutputFileName();
    PrepareOutputFile();
    PrepareHistFit();

    ThreePointGrr = new TGraph();
    fit_linear_fit = new TF1("fit_linear_fit","pol1",-1,actual_xpos[2]+2);
}

void ResidualCompAna::PrepareInputData() {

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

void ResidualCompAna::PrepareOutputFileName() {
    
    if (output_file_name_suffix.size() > 0 && output_file_name_suffix[0] != '_') {
        output_file_name_suffix = "_" + output_file_name_suffix;
    }

    output_filename = (isData) ? "Data" : "MC";
    output_filename += "_Residual";
    output_filename += Form("_Run%d", runnumber);
    output_filename += Form("_Column%d", Selected_Column);

    output_filename += output_file_name_suffix;

    output_filename += ".root";
} 

void ResidualCompAna::PrepareOutputFile() {

    file_out = new TFile(Form("%s/%s", output_directory.c_str(), output_filename.c_str()), "RECREATE");

    tree_out = new TTree("tree","tree");

    tree_out->Branch("l1_alignment_correction", &l1_alignment_correction);
    tree_out->Branch("l0l1_slope_correction", &l0l1_slope_correction);
    tree_out->Branch("ClusPhiSize_cut", &ClusPhiSize_cut);
    tree_out->Branch("ClusAdc_cut", &ClusAdc_cut);
    tree_out->Branch("slope_cut", &slope_cut);
    tree_out->Branch("pos_cut", &pos_cut);
}

void ResidualCompAna::PrepareHistFit() {

    // h1D_l1_alignment_before = new TH1D("h1D_l1_alignment_before","h1D_l1_alignment_before;L1 - (L2L0 interpolation) [mm];Entries",50,-0.4,1);
    // h1D_l1_alignment_after = new TH1D("h1D_l1_alignment_after","h1D_l1_alignment_after;L1 - (L2L0 interpolation) [mm];Entries",50,-0.4,1);
    
    // h1D_l0l1_slope_before = new TH1D("h1D_l0l1_slope_before","h1D_l0l1_slope_before;Slope (L1 - L0);Entries",50,-0.05,0.05);
    // h1D_l0l1_slope_after = new TH1D("h1D_l0l1_slope_after","h1D_l0l1_slope_after;Slope (L1 - L0);Entries",50,-0.05,0.05);

    fit_l1_alignment_before = new TF1("fit_l1_alignment_before","gaus",-5,5);
    fit_l1_alignment_before -> SetNpx(1000);

    fit_l1_alignment_after = new TF1("fit_l1_alignment_after","gaus",-5,5);
    fit_l1_alignment_after -> SetNpx(1000);

    fit_l0l1_slope_before = new TF1("fit_l0l1_slope_before","gaus",-5,5);
    fit_l0l1_slope_before -> SetNpx(1000);

    fit_l0l1_slope_after = new TF1("fit_l0l1_slope_after","gaus",-5,5);
    fit_l0l1_slope_after -> SetNpx(1000);


    h1D_l1_residual = new TH1D("h1D_l1_residual","h1D_l1_residual;L1 - (L2L0 interpolation) [mm];Entries",50,-1,1);
    h1D_l2_residual = new TH1D("h1D_l2_residual","h1D_l2_residual;L2 - (L1L0 extrapolation) [mm];Entries",50,-1,1);
    h1D_scattering = new TH1D("h1D_scattering","h1D_scattering; (L2L1 slope) - (L1L0 slope)",50,-0.05,0.05);

}

double ResidualCompAna::Get_l1_alignment() {
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
    std::cout<<"post correction, l1_alignment, new peak location: "<<fit_l1_alignment_after -> GetParameter(1)<<std::endl;

    return l1_alignment_correction;
}

TH1D * ResidualCompAna::Sub_Get_l1_alignment_hist(){
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

double ResidualCompAna::Get_l0l1_slope(){
    if(l0l1_slope_correction == 0.) {
        std::cout<<"l0l1_slope_correction is 0, then measure it!"<<std::endl;

        h1D_l0l1_slope_before = (TH1D *) (Sub_Get_l0l1_slope_hist()) -> Clone("h1D_l0l1_slope_before");
        h1D_l0l1_slope_before -> SetTitle("h1D_l0l1_slope_before");
        h1D_l0l1_slope_before -> Fit(fit_l0l1_slope_before,"N");

        l0l1_slope_correction = -1 * fit_l0l1_slope_before -> GetParameter(1);
        std::cout<<"l0l1_slope_correction : "<<l0l1_slope_correction<<std::endl;
    }
    else {
        std::cout<<"l0l1_slope_correction is already set to "<<l0l1_slope_correction<<std::endl;
    }

    h1D_l0l1_slope_after = (TH1D *) (Sub_Get_l0l1_slope_hist()) -> Clone("h1D_l0l1_slope_after");
    h1D_l0l1_slope_after -> SetTitle("h1D_l0l1_slope_after");
    h1D_l0l1_slope_after -> Fit(fit_l0l1_slope_after,"N");
    std::cout<<"post correction, l0l1 slope, new peak: "<<fit_l0l1_slope_after -> GetParameter(1)<<std::endl;

    return l0l1_slope_correction;
}

TH1D * ResidualCompAna::Sub_Get_l0l1_slope_hist(){
    TH1D * h1D_l0l1_slope_temp = new TH1D("",";Slope (L1 - L0);Entries",50,-0.05,0.05);
    h1D_l0l1_slope_temp -> Reset("ICESM");

    int count_1 = 0;
    int count_2 = 0;
    int count_3 = 0;
    int count_4 = 0;
    int count_5 = 0;

    // note : [event][chip][layer][clusters]
    for (int i = 0; i < cluster_info_NoL1Aligned_vec.size(); i++) {
        if (i % 5000 == 0){
            std::cout<<"Doing l0l1_slope, Event : "<<i<<std::endl;
        }

        count_1++;

        std::vector<ClusInfo> CheckedCol_L0_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][0];
        std::vector<ClusInfo> CheckedCol_L1_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][1];
        std::vector<ClusInfo> CheckedCol_L2_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][2];

        int N_cluster_l0 = 0;
        int N_cluster_l1 = 0;
        int N_cluster_l2 = 0;

        // note : calculate the total_N_cluster in each layer
        for (int col_i = 0; col_i < 13; col_i++)
        {
            N_cluster_l0 += cluster_info_NoL1Aligned_vec[i][col_i][0].size();
            N_cluster_l1 += cluster_info_NoL1Aligned_vec[i][col_i][1].size();
            N_cluster_l2 += cluster_info_NoL1Aligned_vec[i][col_i][2].size();
        } 

        // note : only one cluster in each layer
        if ( N_cluster_l0 != 1 || N_cluster_l1 != 1 || N_cluster_l2 != 1 ) {continue;}

        count_2++;
        

        // note : single cluster in l0 and l2
        if (CheckedCol_L0_vec.size() != 1 || CheckedCol_L1_vec.size() != 1 || CheckedCol_L2_vec.size() != 1) {continue;}

        count_3++;

        // todo : the cluster adc cut
        // if (CheckedCol_L0_vec[0].adc <= ClusAdc_cut || CheckedCol_L1_vec[0].adc <= ClusAdc_cut || CheckedCol_L2_vec[0].adc <= ClusAdc_cut) {continue;}    
    
        double l0l1_slope = ((CheckedCol_L1_vec[0].pos + l1_alignment_correction) - CheckedCol_L0_vec[0].pos) / actual_xpos[1] + l0l1_slope_correction;

        h1D_l0l1_slope_temp -> Fill(l0l1_slope);

    } // note : loop of cluster_info_NoL1Aligned_vec

    std::cout<<"count_1 : "<<count_1<<std::endl;
    std::cout<<"count_2 : "<<count_2<<std::endl;
    std::cout<<"count_3 : "<<count_3<<std::endl;

    return h1D_l0l1_slope_temp;
}

void ResidualCompAna::GetHistsForComp(double slope_cut_in, double pos_cut_in){
    slope_cut = slope_cut_in;
    pos_cut = pos_cut_in;

    // note : [event][chip][layer][clusters]
    for (int i = 0; i < cluster_info_NoL1Aligned_vec.size(); i++) {
        if (i % 5000 == 0){
            std::cout<<"Doing l1_alignment, Event : "<<i<<std::endl;
        }

        std::vector<ClusInfo> CheckedCol_L0_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][0];
        std::vector<ClusInfo> CheckedCol_L1_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][1];
        std::vector<ClusInfo> CheckedCol_L2_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column-1][2];

        int N_cluster_l0 = 0;
        int N_cluster_l1 = 0;
        int N_cluster_l2 = 0;

        // note : calculate the total_N_cluster in each layer
        for (int col_i = 0; col_i < 13; col_i++)
        {
            N_cluster_l0 += cluster_info_NoL1Aligned_vec[i][col_i][0].size();
            N_cluster_l1 += cluster_info_NoL1Aligned_vec[i][col_i][1].size();
            N_cluster_l2 += cluster_info_NoL1Aligned_vec[i][col_i][2].size();
        } 

        // note : only one cluster in each layer
        if ( N_cluster_l0 != 1 || N_cluster_l1 != 1 || N_cluster_l2 != 1 ) {continue;}
        

        // note : single cluster in l0 and l2
        if (CheckedCol_L0_vec.size() != 1 || CheckedCol_L1_vec.size() != 1 || CheckedCol_L2_vec.size() != 1) {continue;}

        // todo : the cluster adc cut
        // if (CheckedCol_L0_vec[0].adc <= ClusAdc_cut || CheckedCol_L1_vec[0].adc <= ClusAdc_cut || CheckedCol_L2_vec[0].adc <= ClusAdc_cut) {continue;}    
    
        double l0l1_slope = ((CheckedCol_L1_vec[0].pos + l1_alignment_correction) - CheckedCol_L0_vec[0].pos) / actual_xpos[1] + l0l1_slope_correction;
        double l0l1_pos = ((CheckedCol_L1_vec[0].pos + l1_alignment_correction) + CheckedCol_L0_vec[0].pos) / 2.;

        // note : the slope cut
        if ( fabs(l0l1_slope) >= slope_cut ) {continue;}

        // note : the position cut
        if ( fabs(l0l1_pos) >= pos_cut ) {continue;}

        h1D_l1_residual -> Fill(
            (CheckedCol_L1_vec[0].pos + l1_alignment_correction) - ( (CheckedCol_L2_vec[0].pos + CheckedCol_L0_vec[0].pos) ) / 2.
        ); 
        // h1D_l2_residual
        h1D_scattering -> Fill(
            (( CheckedCol_L2_vec[0].pos - (CheckedCol_L1_vec[0].pos + l1_alignment_correction) )/actual_xpos[1]) - (( (CheckedCol_L1_vec[0].pos + l1_alignment_correction) - CheckedCol_L0_vec[0].pos )/actual_xpos[1])
        );

    } // note : loop of cluster_info_NoL1Aligned_vec
}

void ResidualCompAna::EndRun(){
    file_out->cd();
    if (h1D_l1_alignment_before != nullptr) {h1D_l1_alignment_before->Write();}
    if (h1D_l1_alignment_after != nullptr) {h1D_l1_alignment_after->Write();}
    if (h1D_l0l1_slope_before != nullptr) {h1D_l0l1_slope_before->Write();}
    if (h1D_l0l1_slope_after != nullptr) {h1D_l0l1_slope_after->Write();}

    h1D_l1_residual->Write();
    h1D_l2_residual->Write();
    h1D_scattering->Write();

    TCanvas * c1 = new TCanvas("c1","c1",800,600);
    c1 -> cd();

    if (h1D_l1_alignment_before != nullptr) {
        h1D_l1_alignment_before->Draw("hist");
        fit_l1_alignment_before -> Draw("lsame");
        c1 -> Write("c1_h1D_l1_alignment_before");
        c1 -> Clear();
    }

    if (h1D_l1_alignment_after != nullptr) {
        h1D_l1_alignment_after->Draw("hist");
        fit_l1_alignment_after -> Draw("lsame");
        c1 -> Write("c1_h1D_l1_alignment_after");
        c1 -> Clear();
    }

    if (h1D_l0l1_slope_before != nullptr) {
        h1D_l0l1_slope_before->Draw("hist");
        fit_l0l1_slope_before -> SetNpx(1000);
        fit_l0l1_slope_before -> Draw("lsame");
        c1 -> Write("c1_h1D_l0l1_slope_before");
        c1 -> Clear();
    }

    if (h1D_l0l1_slope_after != nullptr) {
        h1D_l0l1_slope_after->Draw("hist");
        fit_l0l1_slope_after -> SetNpx(1000);
        fit_l0l1_slope_after -> Draw("lsame");
        c1 -> Write("c1_h1D_l0l1_slope_after");
        c1 -> Clear();
    }


    tree_out->Fill();
    tree_out->Write();

    file_out->Close();
}