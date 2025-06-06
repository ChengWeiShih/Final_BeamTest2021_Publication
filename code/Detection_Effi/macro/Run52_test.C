#include "../DetectionEffiAna.h"
R__LOAD_LIBRARY(../libDetectionEffiAna.so)

void Characterize_Rate1D (TH1D *hist,  int hcolour)
{
	float ratio = 7.0/3.0;
	hist -> SetTitle("");
	hist -> SetLineColor(hcolour);
	hist -> SetMarkerColor(hcolour);
	hist -> SetMarkerStyle(20);
	hist -> SetMarkerSize(0.8);
	hist -> SetFillColor(hcolour);

	hist -> GetXaxis() -> SetTitleSize   (0.052*ratio);
	hist -> GetXaxis() -> SetTitleOffset (0.940);
	hist -> GetXaxis() -> SetLabelSize   (0.042*ratio);
	hist -> GetXaxis() -> SetLabelOffset (0.004*ratio);

	hist -> GetYaxis() -> SetTitle("Data/MC");
	hist -> GetYaxis() -> SetTitleSize(0.052*ratio);
	hist -> GetYaxis() -> SetTitleOffset(0.360);
	hist -> GetYaxis() -> SetLabelSize(0.042*ratio);
	hist -> GetYaxis() -> SetLabelOffset(0.006);
	hist -> GetYaxis() -> SetRangeUser(0.4, 3.2);
	hist -> GetYaxis() -> SetNdivisions(505);

    hist -> SetStats(0);
}

void Characterize_Hist1D (TH1D *hist,  int statsbox = 0, TString color_ID = "#1A3947", bool set_logY = false, bool do_comp = false, TString title = "", TString X_title = "", TString Y_title = "", bool X_label_less = false, bool Y_label_far = false )
{
	float ratio = 1.;
	hist -> SetTitle       ("");

	hist -> SetLineColor   (TColor::GetColor(color_ID));
    hist -> SetLineWidth   (3);

	hist -> SetMarkerColor (TColor::GetColor(color_ID));
	hist -> SetMarkerStyle (20);
	hist -> SetMarkerSize  (0.8);
	// hist -> SetFillColor   (TColor::GetColor(color_ID));

	hist -> GetXaxis() -> SetTitleSize   (0.05);
	hist -> GetXaxis() -> SetTitleOffset (0.8);

	hist -> GetXaxis() -> SetLabelSize   (0.042*ratio);

    double X_label_offset = (do_comp) ? 0.01 : 0.004;
    hist -> GetXaxis() -> SetLabelOffset (X_label_offset*ratio);
    

	// hist -> GetYaxis() -> SetTitle       ("Data/MC");
	hist -> GetYaxis() -> SetTitleSize   (0.05);
    
    // double amount_of_y_offset = (set_logY) ? 1.5 : 0.94;
    double amount_of_y_offset = (set_logY) ? 0.94 : 0.94;
	hist -> GetYaxis() -> SetTitleOffset (amount_of_y_offset);

    if (Y_label_far == true)
    {
        hist -> GetYaxis() -> SetTitleOffset (1.5);
    }

	hist -> GetYaxis() -> SetLabelSize   (0.042*ratio);
	hist -> GetYaxis() -> SetLabelOffset (0.006);
	// hist -> GetYaxis() -> SetRangeUser   (0.4, 3);

    hist -> SetTitle(title);
    hist -> GetXaxis() -> SetTitle(X_title);
    hist -> GetYaxis() -> SetTitle(Y_title);

	hist -> GetYaxis() -> SetNdivisions  (505);

    if (X_label_less == true)
    {
        hist -> GetXaxis() -> SetNdivisions  (505);
    }


    if (statsbox == 0) {hist -> SetStats(0);}
}

void Characterize_Pad (TPad *pad, float left = 0.15, float right = 0.1, float top = 0.1, float bottom = 0.12, bool set_logY = false, int setgrid_bool = 0)
{
	if (setgrid_bool == true) {pad -> SetGrid (1, 1);}
	pad -> SetLeftMargin   (left);
	pad -> SetRightMargin  (right);
	pad -> SetTopMargin    (top);
	pad -> SetBottomMargin (bottom);
    pad -> SetTicks(1,1);
    if (set_logY == true)
    {
        pad -> SetLogy (1);
    }
	
}

// note : for the comparison of data and MC
// note : titles_vec[0] : x
// note : titles_vec[1] : y
void dataMC_comp (TH1D* hist_data, TH1D* hist_MC, TString folder_direction, TString plot_name, vector<TString> titles_vec, vector<TString> legend_vec, bool linear_or_log = false, bool statsbox_bool = false)
{
    TCanvas * c3 = new TCanvas("c3","c3",850,800);
    c3 -> cd();

    // note : 0) is for SetGrid
    // note : 0,0) is for logY
    TPad *pad_obj = new TPad(Form("pad_obj"), "", 0.0, 0.30, 1.0, 1.0);
    Characterize_Pad(pad_obj, 0.15, 0.1, 0.1, 0.016 , 0, 0);
    pad_obj -> Draw();

    TPad *pad_ratio = new TPad(Form("pad_ratio"), "", 0.0, 0.0, 1.0, 0.30);
    Characterize_Pad(pad_ratio, 0.15, 0.1, 0.04, 0.350, 0, 1);
    pad_ratio -> Draw();

    TLegend *legend1 = new TLegend (0.7, 0.55, 0.85, 0.7);
    legend1 -> SetBorderSize(0);
	legend1 -> SetTextSize (0.050);
	// legend1 -> SetNColumns (4);

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    // ltx->SetTextAlign(31);

    if (linear_or_log == false) // note : linear
    {
        plot_name += "_Linear";
    }
    else if (linear_or_log == true) // note : log
    {   
        pad_obj -> SetLogy();
        plot_name += "_Log";
    }

    Characterize_Hist1D(hist_data,0,"#1A3947",linear_or_log,true ,"",titles_vec[0],titles_vec[1]);
    // hist_data -> GetXaxis() -> SetTitle(titles_vec[0]);
    // hist_data -> GetYaxis() -> SetTitle(titles_vec[1]);
    // hist_data -> GetYaxis() -> SetRangeUser(0,1);

    Characterize_Hist1D(hist_MC,0,"#A08144",linear_or_log,true, "",titles_vec[0],titles_vec[1]);
    // hist_MC -> GetXaxis() -> SetTitle(titles_vec[0]);
    // hist_MC -> GetYaxis() -> SetTitle(titles_vec[1]);
    // hist_MC -> SetTitle(titles_vec[2]);
    // hist_MC -> GetYaxis() -> SetRangeUser(0,1);

    if (statsbox_bool == false)
    {
        hist_MC -> SetStats(0); // note : remove the box
        hist_data -> SetStats(0);
    } 

    // note : normalize
    hist_data->Scale(1./hist_data->Integral(-1,-1));
    hist_MC->Scale(1./hist_MC->Integral(-1,-1));

    legend1 -> AddEntry(hist_data, Form("%s",legend_vec[0].Data()),  "pl");
    legend1 -> AddEntry(hist_MC, Form("%s",legend_vec[1].Data()),  "f");
    
     
    double Y_axis_max = (linear_or_log) ? 10 : 1;
    hist_MC -> SetMaximum(Y_axis_max);


    TH1D * hist_ratio = (TH1D*)hist_data -> Clone(); 
    hist_ratio -> Divide (hist_MC);

    Characterize_Rate1D(hist_ratio,1);
    hist_ratio -> SetMaximum(4.0);
    hist_ratio -> SetMinimum(0.0);
    
    pad_obj -> cd();
    hist_MC -> Draw("hist");
    hist_data -> Draw("ep same");
    legend1 -> Draw("same");

    ltx->DrawLatex(0.45, 0.81, Form("#it{sPHENIX INTT} #bf{Beam Test 2021}"));

    pad_ratio -> cd();
    hist_ratio -> Draw("ep");

    TLine * ratio_1_line = new TLine(hist_ratio -> GetXaxis() -> GetXmin(), 1., hist_ratio->GetXaxis()->GetXmax(), 1.);
    ratio_1_line -> SetLineStyle(9);
    ratio_1_line -> SetLineWidth(2);
    ratio_1_line -> SetLineColor(2);

    ratio_1_line -> Draw("lsame");
    

    // TString output_plot_name = plot_name.ReplaceAll(" ","_");
    // output_plot_name = output_plot_name.ReplaceAll(",","");

    c3 -> Print( Form("%s/%s.pdf",folder_direction.Data(),plot_name.Data()) );
    c3 -> Clear();
    
}


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
std::string data_output_directory = data_input_directory + "/DetectionEffi_test";

std::string MC_input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/DetectionEffi_test/MC";
std::string MC_input_filename = "cluster_information_offset-0.0000_adcinfo_SingleTrigger.root";
std::string MC_output_directory = MC_input_directory + "/DetectionEffi_test";

std::tuple<std::string, double, double, double, DetectionEffiAna*> Run52_mother(
    bool isData = true,
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
 
    std::string final_input_directory = (isData) ? data_input_directory : MC_input_directory;
    std::string final_input_filename = (isData) ? data_input_filename : MC_input_filename;
    std::string final_output_directory = (isData) ? data_output_directory : MC_output_directory;

    DetectionEffiAna * Run52Data = new DetectionEffiAna(
        isData,
        52,
        final_input_directory,
        final_input_filename,
        final_output_directory,
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

    // Run52Data -> EndRun();

    return {
        Run52Data -> GetOutputFileName(),
        Run52Data -> Get_Final_Effi()[0],
        Run52Data -> Get_Final_Effi()[1],
        Run52Data -> Get_Final_Effi()[2],
        Run52Data
    };
}



int Run52_test() {

    bool data_isData = true;
    bool MC_isData = false;

    string BaseLine_data_output_file_name_suffix = "BaseLine";
    int    BaseLine_Selected_Column = 8; // note : 1 to 13,
    double BaseLine_Effi_slope_cut = 0.01; // note : slope cut
    double BaseLine_Effi_pos_cut = 1.0; // note : noise hit distance
    int    BaseLine_Effi_boundary_cut = 8;// note : boundary cut
    double BaseLine_Effi_cluster_adc_value_requirement = 15; // note : cluster adc value requirement

    double BaseLine_ClusDist_slope_cut = 0.01; // note : slope cut
    double BaseLine_ClusDist_pos_cut = 0.234; // note : noise hit distance

    bool   BaseLine_prepare_ClusDist = true;
    
    // Division: BaseLine ---------------------------------------------------------------------------------
    std::tuple<std::string, double, double, double, DetectionEffiAna*> tuple_BaseLine = Run52_mother(
        data_isData,
        BaseLine_data_output_file_name_suffix,
        BaseLine_Selected_Column,
        BaseLine_Effi_slope_cut,
        BaseLine_Effi_pos_cut,
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );

    // Division: Vary Column ---------------------------------------------------------------------------------
    string VaryColumn_data_output_file_name_suffix = "VaryColumn";
    int    VaryColumn_Selected_Column = 9; // note : 1 to 13,

    std::tuple<std::string, double, double, double, DetectionEffiAna*> tuple_VaryColumn = Run52_mother(
        data_isData,
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

    std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<std::get<0>(tuple_BaseLine)<<", "<<std::get<1>(tuple_BaseLine)<<", "<<std::get<2>(tuple_BaseLine)<<", "<<std::get<3>(tuple_BaseLine)<<std::endl;
    std::cout<<std::endl;

    std::cout<<std::get<0>(tuple_VaryColumn)<<", "<<std::get<1>(tuple_VaryColumn)<<", "<<std::get<2>(tuple_VaryColumn)<<", "<<std::get<3>(tuple_VaryColumn)<<std::endl;
    std::cout<<std::endl;



    // Division: BaseLine ---------------------------------------------------------------------------------
    std::tuple<std::string, double, double, double, DetectionEffiAna*> MC_tuple_BaseLine = Run52_mother(
        MC_isData,
        BaseLine_data_output_file_name_suffix,
        10,
        BaseLine_Effi_slope_cut,
        BaseLine_Effi_pos_cut,
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );

    // Division: Vary Column ---------------------------------------------------------------------------------

    std::tuple<std::string, double, double, double, DetectionEffiAna*> MC_tuple_VaryColumn = Run52_mother(
        MC_isData,
        VaryColumn_data_output_file_name_suffix, // note : changed
        11,              // note : changed
        BaseLine_Effi_slope_cut,
        BaseLine_Effi_pos_cut,
        BaseLine_Effi_boundary_cut,
        BaseLine_Effi_cluster_adc_value_requirement,

        BaseLine_ClusDist_slope_cut,
        BaseLine_ClusDist_pos_cut,

        false
    );

    std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<std::get<0>(MC_tuple_BaseLine)<<", "<<std::get<1>(MC_tuple_BaseLine)<<", "<<std::get<2>(MC_tuple_BaseLine)<<", "<<std::get<3>(MC_tuple_BaseLine)<<std::endl;
    std::cout<<std::endl;

    std::cout<<std::get<0>(MC_tuple_VaryColumn)<<", "<<std::get<1>(MC_tuple_VaryColumn)<<", "<<std::get<2>(MC_tuple_VaryColumn)<<", "<<std::get<3>(MC_tuple_VaryColumn)<<std::endl;
    std::cout<<std::endl;

    
    TH1D * data_baseline_effi_l1_residual = (std::get<4>(tuple_BaseLine)) -> Get_h1D_effi_l1_residual();
    TH1D * data_baseline_effi_l1_residual_wide = (std::get<4>(tuple_BaseLine)) -> Get_h1D_effi_l1_residual_wide();

    TH1D * data_varycolumn_effi_l1_residual = (std::get<4>(tuple_VaryColumn)) -> Get_h1D_effi_l1_residual();
    TH1D * data_varycolumn_effi_l1_residual_wide = (std::get<4>(tuple_VaryColumn)) -> Get_h1D_effi_l1_residual_wide();

    TH1D * MC_baseline_effi_l1_residual = (std::get<4>(MC_tuple_BaseLine)) -> Get_h1D_effi_l1_residual();
    TH1D * MC_baseline_effi_l1_residual_wide = (std::get<4>(MC_tuple_BaseLine)) -> Get_h1D_effi_l1_residual_wide();

    TH1D * MC_varycolumn_effi_l1_residual = (std::get<4>(MC_tuple_VaryColumn)) -> Get_h1D_effi_l1_residual();
    TH1D * MC_varycolumn_effi_l1_residual_wide = (std::get<4>(MC_tuple_VaryColumn)) -> Get_h1D_effi_l1_residual_wide();


    dataMC_comp(data_baseline_effi_l1_residual, data_varycolumn_effi_l1_residual, MC_output_directory, Form("c1_data_U8U9_l1_residual"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data U8", "Data U9"}, true, false);
    dataMC_comp(data_baseline_effi_l1_residual_wide, data_varycolumn_effi_l1_residual_wide, MC_output_directory, Form("c1_data_U8U9_l1_residual_wide"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data U8", "Data U9"}, true, false);

    dataMC_comp(MC_baseline_effi_l1_residual, MC_varycolumn_effi_l1_residual, MC_output_directory, Form("c1_MC_U10U11_l1_residual"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"MC U10", "MC U11"}, true, false);
    dataMC_comp(MC_baseline_effi_l1_residual_wide, MC_varycolumn_effi_l1_residual_wide, MC_output_directory, Form("c1_MC_U10U11_l1_residual_wide"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"MC U10", "MC U11"}, true, false);

    dataMC_comp(data_baseline_effi_l1_residual, MC_baseline_effi_l1_residual, MC_output_directory, Form("c1_dataU8_MCU10_l1_residual"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data U8", "MC U10"}, true, false);
    dataMC_comp(data_baseline_effi_l1_residual_wide, MC_baseline_effi_l1_residual_wide, MC_output_directory, Form("c1_dataU8_MCU10_l1_residual_wide"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data U8", "MC U10"}, true, false);

    dataMC_comp(data_varycolumn_effi_l1_residual, MC_varycolumn_effi_l1_residual, MC_output_directory, Form("c1_dataU9_MCU11_l1_residual"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data U9", "MC U11"}, true, false);
    dataMC_comp(data_varycolumn_effi_l1_residual_wide, MC_varycolumn_effi_l1_residual_wide, MC_output_directory, Form("c1_dataU9_MCU11_l1_residual_wide"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data U9", "MC U11"}, true, false);

    // dataMC_comp(data_baseline_effi_l1_residual, MC_baseline_effi_l1_residual, MC_output_directory, Form("c1_dataU8_MCU10_l1_residual"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data", "MC"}, true, false);
    // dataMC_comp(data_baseline_effi_l1_residual_wide, MC_baseline_effi_l1_residual_wide, MC_output_directory, Form("c1_dataU8_MCU10_l1_residual_wide"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data", "MC"}, true, false);

    // dataMC_comp(data_varycolumn_effi_l1_residual, MC_varycolumn_effi_l1_residual, MC_output_directory, Form("c1_dataU9_MCU11_l1_residual"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data", "MC"}, true, false);
    // dataMC_comp(data_varycolumn_effi_l1_residual_wide, MC_varycolumn_effi_l1_residual_wide, MC_output_directory, Form("c1_dataU9_MCU11_l1_residual_wide"), {"L1 - (L2L0 interpolation)", "Entries (A.U.)"}, {"Data", "MC"}, true, false);

    (std::get<4>(tuple_BaseLine)) -> EndRun();
    (std::get<4>(tuple_VaryColumn)) -> EndRun();
    (std::get<4>(MC_tuple_BaseLine)) -> EndRun();
    (std::get<4>(MC_tuple_VaryColumn)) -> EndRun();


    return 888;
}