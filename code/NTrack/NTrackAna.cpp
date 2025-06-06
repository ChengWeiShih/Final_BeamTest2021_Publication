#include "NTrackAna.h"

NTrackAna::NTrackAna(
    bool isData_in,
    int runnumber_in,
    std::string input_directory_in,
    std::string input_filename_in,
    std::string output_directory_in,
    std::string output_file_name_suffix_in
) : isData(isData_in),
    runnumber(runnumber_in),
    input_directory(input_directory_in),
    input_filename(input_filename_in),
    output_directory(output_directory_in),
    output_file_name_suffix(output_file_name_suffix_in)
{
    h1D_l1_alignment_before_vec.clear();
    h1D_l1_alignment_after_vec.clear();

    fit_l1_alignment_before_vec.clear();
    fit_l1_alignment_after_vec.clear();

    system(Form("mkdir -p %s", output_directory.c_str()) );

    cluster_info_NoL1Aligned_vec.clear();

    h2D_TrackLines_vec.clear();

    PrepareInputData();
    PrepareOutputFileName();
    PrepareOutputFile();
    PrepareHistFit();

    l1_alignment_correction_vec = std::vector<double>(13,0);
    l1_alignment_correction_error_vec = std::vector<double>(13,0);
    l1_alignment_correction_final_vec = std::vector<double>(13,0);

    skip_column.clear();
}

void NTrackAna::PrepareInputData() {

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

void NTrackAna::PrepareOutputFileName() {
    
    if (output_file_name_suffix.size() > 0 && output_file_name_suffix[0] != '_') {
        output_file_name_suffix = "_" + output_file_name_suffix;
    }

    output_filename = (isData) ? "Data" : "MC";
    output_filename += "_NTrack";
    output_filename += Form("_Run%d", runnumber);

    output_filename += output_file_name_suffix;

    output_filename += ".root";
}

void NTrackAna::PrepareOutputFile() {

    file_out = new TFile(Form("%s/%s", output_directory.c_str(), output_filename.c_str()), "RECREATE");

    tree_out = new TTree("tree","tree");

    tree_out -> Branch("Column_ID", &Column_ID_out);
    tree_out -> Branch("l1_alignment_correction_column", &l1_alignment_correction_column_out);
    tree_out -> Branch("l1_alignment_correction_column_final", &l1_alignment_correction_column_final_out);
}

void NTrackAna::PrepareHistFit() {
    
    for (int i = 0; i < 13; i++){
        fit_l1_alignment_before_vec.push_back(new TF1(Form("fit_l1_alignment_Column%d_before", i+1), "gaus", -5, 5));
        fit_l1_alignment_after_vec.push_back(new TF1(Form("fit_l1_alignment_Column%d_after", i+1), "gaus", -5, 5));

        fit_l1_alignment_before_vec.back() -> SetNpx(1000);
        fit_l1_alignment_after_vec.back() -> SetNpx(1000);
    }

    gr_typeA_l1_alignment = new TGraphErrors();
    gr_typeB_l1_alignment = new TGraphErrors();

    fit_typeA_l1_alignment = new TF1("fit_typeA_l1_alignment", "pol1", 0, 14);
    fit_typeB_l1_alignment = new TF1("fit_typeB_l1_alignment", "pol1", 0, 14);

    f1_pol1_unity = new TF1("f1_pol1_unity", "pol1", -10, 10);
    f1_pol1_unity -> SetParameters(0, 1);
    f1_pol1_unity -> SetLineColor(kRed);


    h1D_Column_NTrack = new TH1D("h1D_Column_NTrack","h1D_Column_NTrack;Column ID;Number of tracks",13,0.5,13.5);

    h1D_NTrack_NoZero = new TH1D("h1D_NTrack_NoZero","h1D_NTrack_NoZero;Number of reco. tracks;Entries",13, -3.5, 9.5);

    h1D_NTrackWithZero = new TH1D("h1D_NTrackWithZero","h1D_NTrackWithZero;Number of reco. tracks;Entries",13, -3.5, 9.5);

    h2D_TrackPos = new TH2D("h2D_TrackPos","h2D_TrackPos;Column ID; Y axis [mm]",13, 0.5, 13.5, 128, -9.984, 9.984);

    h1D_l0_ClusPosFitDiff = new TH1D("h1D_l0_ClusPosFitDiff", "h1D_l0_ClusPosFitDiff;ClusPos - Fit [mm] (L0);Entries", 25., -1, 1);
    h1D_l1_ClusPosFitDiff = new TH1D("h1D_l1_ClusPosFitDiff", "h1D_l1_ClusPosFitDiff;ClusPos - Fit [mm] (L1);Entries", 25., -1, 1);
    h1D_l2_ClusPosFitDiff = new TH1D("h1D_l2_ClusPosFitDiff", "h1D_l2_ClusPosFitDiff;ClusPos - Fit [mm] (L2);Entries", 25., -1, 1);

    h1D_l1_residual = new TH1D("h1D_l1_residual","h1D_l1_residual;L1 - (L2L0 interpolation) [mm];Entries", 50., -1, 1);

    // h2D_TrackLines = new TH2D("h2D_TrackLines", "h2D_TrackLines;Z axis (beam direction) [mm];Y axis [mm]", 400, -200, 65, 300, -15, 15);

}

void NTrackAna::Get_l1_alignment() {
    
    // note : N columns
    for (int i = 0; i < 13; i++) {

        if(l1_alignment_correction_vec[i] == 0.) {
            std::cout<<Form("l1_alignment_correction_vec[%d] is 0, then measure it!", i)<<std::endl;

            h1D_l1_alignment_before_vec.push_back( (TH1D *) (Sub_Get_l1_alignment_hist(i)) -> Clone(Form("h1D_l1_alignment_Column%d_before",i+1)) );
            h1D_l1_alignment_before_vec.back() -> SetTitle( Form("h1D_l1_alignment_Column%d_before",i+1) );
            h1D_l1_alignment_before_vec.back() -> Fit(fit_l1_alignment_before_vec[i],"N");

            l1_alignment_correction_vec[i] = -1 * fit_l1_alignment_before_vec[i] -> GetParameter(1);
            l1_alignment_correction_error_vec[i] = fit_l1_alignment_before_vec[i] -> GetParError(1);
            std::cout<<Form("l1_alignment_correction_vec[%d] : ", i)<<l1_alignment_correction_vec[i]<<std::endl;
        }
        else{
            std::cout<<Form("l1_alignment_correction_vec[%d] is already set to ", i)<<l1_alignment_correction_vec[i]<<std::endl;
        }

        h1D_l1_alignment_after_vec.push_back( (TH1D *) (Sub_Get_l1_alignment_hist(i)) -> Clone( Form("h1D_l1_alignment_Column%d_after",i+1) ) );
        h1D_l1_alignment_after_vec.back() -> SetTitle( Form("h1D_l1_alignment_Column%d_after",i+1) );
        h1D_l1_alignment_after_vec.back() -> Fit(fit_l1_alignment_after_vec[i],"N");
        std::cout<<"post correction, l1_alignment, new peak: "<<fit_l1_alignment_after_vec[i] -> GetParameter(1)<<std::endl;
    }

    // todo : columns used to fit the l1 alignment correction are here
    // note : type B: 0 to 4
    for (int i = 1; i < 5; i++) { // note : i = 1 to 4 (Column: 2 to 5)
        gr_typeB_l1_alignment -> SetPoint(gr_typeB_l1_alignment -> GetN(), i+1, l1_alignment_correction_vec[i]);
        gr_typeB_l1_alignment -> SetPointError(gr_typeB_l1_alignment -> GetN()-1, 0, l1_alignment_correction_error_vec[i]);
    }

    // note : type A: 5 to 12
    for (int i = 8; i < 12; i++) { // note : i = 8 to 11 (Column: 9 to 12)
        gr_typeA_l1_alignment -> SetPoint(gr_typeA_l1_alignment -> GetN(), i+1, l1_alignment_correction_vec[i]);
        gr_typeA_l1_alignment -> SetPointError(gr_typeA_l1_alignment -> GetN()-1, 0, l1_alignment_correction_error_vec[i]);
    }

    gr_typeB_l1_alignment -> SetMarkerStyle(20);
    gr_typeB_l1_alignment -> SetMarkerSize(0.8);
    gr_typeB_l1_alignment -> SetMarkerColor(kBlack);
    gr_typeB_l1_alignment -> GetXaxis() -> SetTitle("Column ID");
    gr_typeB_l1_alignment -> GetXaxis() -> SetLimits(0, 6);
    gr_typeB_l1_alignment -> GetYaxis() -> SetTitle("l1 alignment correction (TypeB) [mm]");

    gr_typeA_l1_alignment -> SetMarkerStyle(20);
    gr_typeA_l1_alignment -> SetMarkerSize(0.8);
    gr_typeA_l1_alignment -> SetMarkerColor(kBlack);
    gr_typeA_l1_alignment -> GetXaxis() -> SetTitle("Column ID");
    gr_typeA_l1_alignment -> GetXaxis() -> SetLimits(5, 14);
    gr_typeA_l1_alignment -> GetYaxis() -> SetTitle("l1 alignment correction (TypeA) [mm]");


    gr_typeB_l1_alignment -> Fit(fit_typeB_l1_alignment, "N");
    gr_typeA_l1_alignment -> Fit(fit_typeA_l1_alignment, "N");
    std::cout<<"Type B, Column: 2 to 5, l1_alignment_correction function: y = "<<fit_typeB_l1_alignment->GetParameter(0)<<" + "<<fit_typeB_l1_alignment->GetParameter(1)<<" * x"<<std::endl;
    std::cout<<"Type A, Column: 9 to 12, l1_alignment_correction function: y = "<<fit_typeA_l1_alignment->GetParameter(0)<<" + "<<fit_typeA_l1_alignment->GetParameter(1)<<" * x"<<std::endl;

    if (isData) {
        for (int i = 0; i < 5; i++) {// note : i = 0 to 4 (Column: 1 to 5)
            l1_alignment_correction_final_vec[i] = fit_typeB_l1_alignment->GetParameter(0) + fit_typeB_l1_alignment->GetParameter(1) * (i+1);
        }

        for (int i = 5; i < 13; i++) {// note : i = 5 to 12 (Column: 6 to 13)
            l1_alignment_correction_final_vec[i] = fit_typeA_l1_alignment->GetParameter(0) + fit_typeA_l1_alignment->GetParameter(1) * (i+1);
        }
    }
    else {
        for (int i = 0; i < 13; i++) {
            l1_alignment_correction_final_vec[i] = (
                l1_alignment_correction_vec[BeamSpotColumn - 1 - 1] + 
                l1_alignment_correction_vec[BeamSpotColumn - 1] +
                l1_alignment_correction_vec[BeamSpotColumn - 1 + 1]
            ) / 3.;
        }
    }
}

TH1D * NTrackAna::Sub_Get_l1_alignment_hist(int Selected_Column){
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

        std::vector<ClusInfo> CheckedCol_L0_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column][0];
        std::vector<ClusInfo> CheckedCol_L1_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column][1];
        std::vector<ClusInfo> CheckedCol_L2_vec = cluster_info_NoL1Aligned_vec[i][Selected_Column][2];

        std::vector<ClusInfo> CheckedColM1_L0_vec = (Selected_Column-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L1_vec = (Selected_Column-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColM1_L2_vec = (Selected_Column-1 >= 0) ? cluster_info_NoL1Aligned_vec[i][Selected_Column-1][2] : empty_ClusInfo_vec;

        std::vector<ClusInfo> CheckedColP1_L0_vec = (Selected_Column+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column+1][0] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L1_vec = (Selected_Column+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column+1][1] : empty_ClusInfo_vec;
        std::vector<ClusInfo> CheckedColP1_L2_vec = (Selected_Column+1 <= 12) ? cluster_info_NoL1Aligned_vec[i][Selected_Column+1][2] : empty_ClusInfo_vec;

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
        double diff = CheckedCol_L1_vec[0].pos - expected_middle_pos + l1_alignment_correction_vec[Selected_Column];

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

void NTrackAna::GetNTrack(double  noise_hit_distance){

    TGraph * gr_temp = new TGraph();
    gr_temp -> Set(0);

    TF1 * linear_fit = new TF1("linear_fit","pol1",-1,actual_xpos[2]+2);

    std::vector<std::vector<double>> evt_trackline_pos_vec; evt_trackline_pos_vec.clear();

    double chi2_register;
    int cluster_register_l0;
    int cluster_register_l1;
    int cluster_register_l2;

    double hit3_best_fit_picker_info[10];

    int NTrack_Column = 0; // note : number of reco. tracks per event 

    for (int i = 0; i < cluster_info_NoL1Aligned_vec.size(); i++) {

        NTrack_Column = 0;
        evt_trackline_pos_vec.clear();

        // note : [event][chip][layer][clusters]
        if (i % 5000 == 0){
            std::cout<<"Doing NTrack, Event : "<<i<<std::endl;
        }

        for (int col_i = 0; col_i < 13; col_i++) {
            bool find_track = true;

            std::vector<std::vector<ClusInfo>> One_Column_Clus_vec = cluster_info_NoL1Aligned_vec[i][col_i];

            if ( std::find(skip_column.begin(), skip_column.end(), col_i+1) != skip_column.end()) {
                continue;
            }

            while ( One_Column_Clus_vec[0].size() > 0 && One_Column_Clus_vec[1].size() > 0 && One_Column_Clus_vec[2].size() > 0 && find_track == true ){

                // note : try to find the good track in the loops
                for (int l0 = 0; l0 < One_Column_Clus_vec[0].size(); l0 ++ ) { // note : # of clusters of the column of the layer 0

                    for (int l1 = 0; l1 < One_Column_Clus_vec[1].size(); l1 ++ ) { // note : # of clusters of the column of the layer 1

                        for (int l2 = 0; l2 < One_Column_Clus_vec[2].size(); l2 ++ ) { // note : # of clusters of the column of the layer 2

                            double hit3_Y_data[3] = {
                                One_Column_Clus_vec[0][l0].pos,
                                One_Column_Clus_vec[1][l1].pos + l1_alignment_correction_final_vec[col_i],
                                One_Column_Clus_vec[2][l2].pos
                            };

                            gr_temp -> Set(0);
                            gr_temp -> SetPoint(0, actual_xpos[0], hit3_Y_data[0]);
                            gr_temp -> SetPoint(1, actual_xpos[1], hit3_Y_data[1]);
                            gr_temp -> SetPoint(2, actual_xpos[2], hit3_Y_data[2]);
                            gr_temp -> Fit(linear_fit, "QN");


                            if (l0+l1+l2 == 0) {
                                chi2_register = ( linear_fit->GetChisquare()/double (linear_fit->GetNDF()) );
                                cluster_register_l0 = l0;
                                cluster_register_l1 = l1;
                                cluster_register_l2 = l2;


                                hit3_best_fit_picker_info[0] = ( hit3_Y_data[0] - ( linear_fit -> GetParameter(1) * actual_xpos[0] + linear_fit -> GetParameter(0) ) );
                                hit3_best_fit_picker_info[1] = ( hit3_Y_data[1] - ( linear_fit -> GetParameter(1) * actual_xpos[1] + linear_fit -> GetParameter(0) ) );
                                hit3_best_fit_picker_info[2] = ( hit3_Y_data[2] - ( linear_fit -> GetParameter(1) * actual_xpos[2] + linear_fit -> GetParameter(0) ) );
                                
                                hit3_best_fit_picker_info[3] = hit3_Y_data[0];
                                hit3_best_fit_picker_info[4] = hit3_Y_data[1];
                                hit3_best_fit_picker_info[5] = hit3_Y_data[2];

                                hit3_best_fit_picker_info[6] = linear_fit -> GetParameter(0); // note : offset
                                hit3_best_fit_picker_info[7] = linear_fit -> GetParameter(1); // note : slope

                                hit3_best_fit_picker_info[8] = linear_fit -> Eval(actual_xpos[0]);
                                hit3_best_fit_picker_info[9] = linear_fit -> Eval(actual_xpos[2]);

                            }
                            else {
                                if ( linear_fit->GetChisquare()/double (linear_fit->GetNDF()) < chi2_register ) {

                                    chi2_register = ( linear_fit->GetChisquare()/double (linear_fit->GetNDF()) );
                                    cluster_register_l0 = l0;
                                    cluster_register_l1 = l1;
                                    cluster_register_l2 = l2;

                                    hit3_best_fit_picker_info[0] = ( hit3_Y_data[0] - ( linear_fit -> GetParameter(1) * actual_xpos[0] + linear_fit -> GetParameter(0) ) );
                                    hit3_best_fit_picker_info[1] = ( hit3_Y_data[1] - ( linear_fit -> GetParameter(1) * actual_xpos[1] + linear_fit -> GetParameter(0) ) );
                                    hit3_best_fit_picker_info[2] = ( hit3_Y_data[2] - ( linear_fit -> GetParameter(1) * actual_xpos[2] + linear_fit -> GetParameter(0) ) );
                                    
                                    hit3_best_fit_picker_info[3] = hit3_Y_data[0];
                                    hit3_best_fit_picker_info[4] = hit3_Y_data[1];
                                    hit3_best_fit_picker_info[5] = hit3_Y_data[2];

                                    hit3_best_fit_picker_info[6] = linear_fit -> GetParameter(0); // note : offset
                                    hit3_best_fit_picker_info[7] = linear_fit -> GetParameter(1); // note : slope

                                    hit3_best_fit_picker_info[8] = linear_fit -> Eval(actual_xpos[0]);
                                    hit3_best_fit_picker_info[9] = linear_fit -> Eval(actual_xpos[2]);

                                }
                            }

                        } // note : loop of l2
                    } // note : loop of l1
                } // note : loop of l0


                h1D_l0_ClusPosFitDiff -> Fill(hit3_best_fit_picker_info[0]);
                h1D_l1_ClusPosFitDiff -> Fill(hit3_best_fit_picker_info[1]);
                h1D_l2_ClusPosFitDiff -> Fill(hit3_best_fit_picker_info[2]);
                h1D_l1_residual -> Fill(
                    hit3_best_fit_picker_info[4] - (hit3_best_fit_picker_info[3]+hit3_best_fit_picker_info[5])/2.
                );

                if ( 
                        // fabs(hit3_best_fit_picker_info[0]) > noise_hit_distance || 
                        // fabs(hit3_best_fit_picker_info[1]) > noise_hit_distance || 
                        // fabs(hit3_best_fit_picker_info[2]) > noise_hit_distance 

                        fabs(hit3_best_fit_picker_info[4] - (hit3_best_fit_picker_info[3]+hit3_best_fit_picker_info[5])/2.) > noise_hit_distance

                    ) {
                    find_track = false; // note : which means, even if we have at least one cluster in each layer in this event. they can not form a good track.
                }
                else {
                    NTrack_Column += 1;

                    h1D_Column_NTrack -> Fill(col_i + 1);

                    h2D_TrackPos -> Fill(
                        col_i + 1, 
                        (hit3_best_fit_picker_info[3] + hit3_best_fit_picker_info[4] + hit3_best_fit_picker_info[5]) / 3.
                    );

                    evt_trackline_pos_vec.push_back(
                        {
                            actual_xpos[0], hit3_best_fit_picker_info[8],
                            actual_xpos[2], hit3_best_fit_picker_info[9]
                        }
                    );

                    // TH2DSampleLineFill(
                    //     h2D_TrackLines, 
                    //     0.1,
                    //     {actual_xpos[0], hit3_best_fit_picker_info[8]},
                    //     {actual_xpos[2], hit3_best_fit_picker_info[9]}
                    // );


                    One_Column_Clus_vec[0].erase( One_Column_Clus_vec[0].begin() + cluster_register_l0 );
                    One_Column_Clus_vec[1].erase( One_Column_Clus_vec[1].begin() + cluster_register_l1 );
                    One_Column_Clus_vec[2].erase( One_Column_Clus_vec[2].begin() + cluster_register_l2 );

                }

            } // note : while loop

        } // note : loop of column

        if (NTrack_Column != 0) {h1D_NTrack_NoZero -> Fill(NTrack_Column);}

        if (evt_trackline_pos_vec.size() >= 4) {
            h2D_TrackLines_vec.push_back(
                new TH2D(
                    Form("h2D_TrackLines_evt%d", i),
                    Form("h2D_TrackLines_evt%d;Z axis (beam direction) [mm];Y axis [mm]", i),
                    400, -300, 65, 
                    300, -15, 15
                )
            );

            for (int track_i = 0; track_i < evt_trackline_pos_vec.size(); track_i++) {
                TH2DSampleLineFill(
                    h2D_TrackLines_vec.back(), 
                    0.1,
                    {evt_trackline_pos_vec[track_i][0], evt_trackline_pos_vec[track_i][1]},
                    {evt_trackline_pos_vec[track_i][2], evt_trackline_pos_vec[track_i][3]}
                );
            }
        }
        
        h1D_NTrackWithZero -> Fill(NTrack_Column);

    } // note : loop of event


}

void NTrackAna::TH2DSampleLineFill(TH2D * hist_in, double segmentation, std::pair<double,double> inner_clu, std::pair<double,double> outer_clu) {
    double x_min = hist_in -> GetXaxis() -> GetXmin();
    double x_max = hist_in -> GetXaxis() -> GetXmax();
    double y_min = hist_in -> GetYaxis() -> GetXmin();
    double y_max = hist_in -> GetYaxis() -> GetXmax();

    double seg_x, seg_y;
    double angle;
    int n_seg = 0;

    while (true) {
        angle = atan2(inner_clu.second-outer_clu.second, inner_clu.first-outer_clu.first);
        seg_x = (n_seg * segmentation) * cos(angle) + outer_clu.first; // note : atan2(y,x), point.first is the radius
        seg_y = (n_seg * segmentation) * sin(angle) + outer_clu.second;
        
        if ( (seg_x > x_min && seg_x < x_max && seg_y > y_min && seg_y < y_max) != true ) {break;}
        hist_in -> Fill(seg_x, seg_y);
        n_seg += 1;
    }

    n_seg = 1;
    while (true) {
        angle = atan2(inner_clu.second-outer_clu.second, inner_clu.first-outer_clu.first);
        seg_x = (-1 * n_seg * segmentation) * cos(angle) + outer_clu.first; // note : atan2(y,x), point.first is the radius
        seg_y = (-1 * n_seg * segmentation) * sin(angle) + outer_clu.second;
        
        if ( (seg_x > x_min && seg_x < x_max && seg_y > y_min && seg_y < y_max) != true ) {break;}
        hist_in -> Fill(seg_x, seg_y);
        n_seg += 1;
    }
}



void NTrackAna::EndRun() {
    std::cout<<"EndRun()"<<std::endl;

    file_out -> cd();

    for (int i = 0; i < 13; i++) {
        h1D_l1_alignment_before_vec[i] -> Write();
        h1D_l1_alignment_after_vec[i] -> Write();
    }

    TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);
    c1 -> cd();
    gr_typeB_l1_alignment -> Draw("ap");
    fit_typeB_l1_alignment -> Draw("same");
    c1 -> Write("c1_typeB_l1_alignment");

    c1 -> Clear();
    gr_typeA_l1_alignment -> Draw("ap");
    fit_typeA_l1_alignment -> Draw("same");
    c1 -> Write("c1_typeA_l1_alignment");

    // Division : --------------------------------------------------------------------------------
    c1 -> Clear();
    TGraphErrors * gr_l1_alignment_correlation = new TGraphErrors();
    for (int i = 0; i < 13; i++) {
        gr_l1_alignment_correlation -> SetPoint(gr_l1_alignment_correlation->GetN(), l1_alignment_correction_vec[i], l1_alignment_correction_final_vec[i]);
        gr_l1_alignment_correlation -> SetPointError(gr_l1_alignment_correlation->GetN()-1, l1_alignment_correction_error_vec[i], 0);

        Column_ID_out = i + 1;
        l1_alignment_correction_column_out = l1_alignment_correction_vec[i];
        l1_alignment_correction_column_final_out = l1_alignment_correction_final_vec[i];
        tree_out -> Fill();
    }

    gr_l1_alignment_correlation -> SetMarkerStyle(20);
    gr_l1_alignment_correlation -> SetMarkerSize(0.8);
    gr_l1_alignment_correlation -> SetMarkerColor(kBlack);
    gr_l1_alignment_correlation -> GetXaxis() -> SetTitle("l1 alignment correction (Each) [mm]");
    gr_l1_alignment_correlation -> GetYaxis() -> SetTitle("l1 alignment correction (Fit) [mm]");

    gr_l1_alignment_correlation -> Draw("ap");
    f1_pol1_unity -> Draw("lsame");

    c1 -> Write("c1_l1_alignment_correlation");

    //Division : --------------------------------------------------------------------------------

    // h2D_TrackLines -> Write();


    h1D_Column_NTrack -> Write();
    h1D_NTrack_NoZero -> Write();
    h1D_NTrackWithZero -> Write();
    h2D_TrackPos -> Write();

    h1D_l0_ClusPosFitDiff -> Write();
    h1D_l1_ClusPosFitDiff -> Write();
    h1D_l2_ClusPosFitDiff -> Write();

    h1D_l1_residual -> Write();

    for (int i = 0; i < h2D_TrackLines_vec.size(); i++) {
        h2D_TrackLines_vec[i] -> Write();
    }

    tree_out -> Write();
    file_out -> Close();

    file_in -> Close();
}