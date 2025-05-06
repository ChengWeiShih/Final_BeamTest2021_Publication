#include "sPhenixStyle.C"

struct VariationSet{
    TGraphAsymmErrors * effi_pos_gr;
    double all_Final_Effi;
    double all_Final_Effi_StatErrorUp;
    double all_Final_Effi_StatErrorDown;
};

struct RelativeVariationSet{
    TGraph * MaxVariationUp;
    TGraph * MaxVariationDown;

    double all_MaxVariationUp;
    double all_MaxVariationDown;
};

const std::vector<std::string> color_code = {
    "#9e0142",
    "#66c2a5",
    "#f46d43",
    "#3288bd",
    "#fee08b",
    "#5e4fa2",
    "#00A1FF",
    "#FF42A1",
    "#000000",
    
    
    "#abdda4",
    "#e6f598",
    "#fdae61",
    "#d53e4f"
};

const std::vector<std::string> syst_HighLow_color_code = {
    "#9e0142",
    "#FF42A1",

    "#005493",
    "#76D6FF",

    "#fee08b",
    "#f46d43",

    "#367C22",
    "#abdda4"
};

const std::vector<int> marker_code = {
    25,
    28,
    27,
    26
    // 21,
    // 23,
    // 27,
    // 28,
    // 29
};

// note : Run52 single
const double Run52_BeamSpot_region_High = 6.0;
const double Run52_BeamSpot_region_Low = -6.0;

// note : Run89 single 
const double Run89_BeamSpot_region_High = 2.0;
const double Run89_BeamSpot_region_Low = -10.0;


// note : plot config.
double left_pos_edge = -10.;
double right_pos_edge = 10;
int N_bins_pos_hist = 10; // todo : -----> A

// todo : change here, for different runs // todo : -----> B
// note : for the data points 
// note : Run52 single
// double Pos_BeamSpot_region_High = 6.0;
// double Pos_BeamSpot_region_Low = -6.0;

// note : Run89 single 
// double Pos_BeamSpot_region_High = 2.0;
// double Pos_BeamSpot_region_Low = -10.0;

// note : Run52 + Run 89 combined 
double Pos_BeamSpot_region_High = 6.0;
double Pos_BeamSpot_region_Low = -10.0;

TGraph * GetRatioGr(TGraphAsymmErrors * gr_numerator, TGraphAsymmErrors * gr_denominator)
{
    TGraph * gr_out = new TGraph();
    
    if (gr_numerator -> GetN() != gr_denominator -> GetN()){
        std::cout<<"Error : the number of points in the numerator and denominator are not the same"<<std::endl;
        exit(1);
    }

    for (int i = 0; i < gr_numerator -> GetN(); i++){
        double x = 0;
        double y = 0;
        gr_numerator -> GetPoint(i, x, y);

        double x2 = 0;
        double y2 = 0;
        gr_denominator -> GetPoint(i, x2, y2);

        if (fabs(x - x2) > 0.0001) {
            std::cout<<"Error : the x value of the points are not the same"<<std::endl;
            std::cout<<"x1 : "<<x<<", x2 : "<<x2<<std::endl;
            exit(1);
        }

        if (y2 == 0) {continue;}

        gr_out -> SetPoint(gr_out -> GetN(), x, y / y2);
    }

    gr_out -> SetMarkerStyle(gr_numerator -> GetMarkerStyle());
    gr_out -> SetMarkerSize(gr_numerator -> GetMarkerSize());
    gr_out -> SetMarkerColor(gr_numerator -> GetMarkerColor());

    return gr_out;
}

VariationSet MakeEffiPos(std::vector<std::string> input_file_vec) { // note : for combining runs
    std::vector<TFile *> file_vec; file_vec.clear();
    std::vector<TTree *> tree_vec; tree_vec.clear();

    double L0L2Interpolation;
    bool L1Good;

    TH1D * good_event_hist = new TH1D ("","",N_bins_pos_hist,left_pos_edge,right_pos_edge);
    TH1D * total_event_hist = new TH1D ("","",N_bins_pos_hist,left_pos_edge,right_pos_edge);

    TH1D * all_good_event_hist = new TH1D ("","",1,0,1);
    TH1D * all_total_event_hist = new TH1D ("","",1,0,1);

    int all_total_event_count = 0;
    int all_good_event_count = 0;

    good_event_hist -> Reset("ICESM");
    total_event_hist -> Reset("ICESM");

    for (int i = 0; i < input_file_vec.size(); i++) {
        file_vec.push_back(new TFile(input_file_vec[i].c_str()));
        tree_vec.push_back((TTree *)file_vec[i]->Get("tree_effi"));

        tree_vec.back() -> SetBranchAddress("L0L2Interpolation", &L0L2Interpolation);
        tree_vec.back() -> SetBranchAddress("L1Good", &L1Good);

        for ( int evt_i = 0; evt_i < tree_vec.back()->GetEntries(); evt_i++ ) {
            tree_vec.back() -> GetEntry(evt_i);

            // todo : the evnet rejection if it outside the defined beam spot region
            if (input_file_vec[i].find("Run52") != std::string::npos) {
                if (L0L2Interpolation > Run52_BeamSpot_region_High || L0L2Interpolation < Run52_BeamSpot_region_Low) {
                    continue;
                }
            }
            else if (input_file_vec[i].find("Run89") != std::string::npos) {
                if (L0L2Interpolation > Run89_BeamSpot_region_High || L0L2Interpolation < Run89_BeamSpot_region_Low) {
                    continue;
                }
            }
            else {
                std::cout<<"!!!!!!!!! It's nither Run52 nor Run89 !!!!!!!!!"<<std::endl;
                std::cout<<"input_file_vec[i] : "<<input_file_vec[i].c_str()<<std::endl;
                std::cout<<"All region is used"<<std::endl;
            }

            all_total_event_count++;

            total_event_hist -> Fill(L0L2Interpolation);
            if (L1Good == 1) {
                good_event_hist -> Fill(L0L2Interpolation);
                all_good_event_count++;
            }
        }
    }

    TEfficiency * detection_effi_pos = new TEfficiency(*good_event_hist,*total_event_hist);

    all_good_event_hist ->SetBinContent(1, all_good_event_count);
    all_total_event_hist->SetBinContent(1, all_total_event_count);
    TEfficiency * all_detection_effi = new TEfficiency(*all_good_event_hist,*all_total_event_hist);

    double all_Final_Effi = all_detection_effi->GetEfficiency(1)*100.;
    double all_Final_Effi_StatErrorUp = all_detection_effi->GetEfficiencyErrorUp(1)*100.;
    double all_Final_Effi_StatErrorDown = all_detection_effi->GetEfficiencyErrorLow(1)*100.;

    TGraphAsymmErrors * effi_pos_gr = new  TGraphAsymmErrors();
    for (int i = 1; i <= N_bins_pos_hist; i++) {
        double bin_center = good_event_hist -> GetBinCenter(i);
        double bin_width  = good_event_hist -> GetBinWidth(i)/2.;

        double bin_eff            = detection_effi_pos -> GetEfficiency(i) * 100;
        double bin_eff_error_up   = detection_effi_pos -> GetEfficiencyErrorUp(i) * 100;
        double bin_eff_error_down = detection_effi_pos -> GetEfficiencyErrorLow(i) * 100;


        if (bin_center > Pos_BeamSpot_region_Low && bin_center < Pos_BeamSpot_region_High) {
            effi_pos_gr -> SetPoint(effi_pos_gr -> GetN(), bin_center, bin_eff);
            effi_pos_gr -> SetPointError(effi_pos_gr -> GetN() - 1, bin_width, bin_width, bin_eff_error_down, bin_eff_error_up);
        }
    }

    std::cout<<"input_file_vec.size() : "<<input_file_vec.size()<<std::endl;
    for (int i = 0; i < input_file_vec.size(); i++) {
        std::cout<<input_file_vec[i].c_str()<<std::endl;
    }
    std::cout<<"all_good_event_count : "<<all_good_event_count<<std::endl;
    std::cout<<"all_total_event_count : "<<all_total_event_count<<std::endl;
    std::cout<<"all_Final_Effi : "<<all_Final_Effi<<", all_Final_Effi_StatErrorUp : "<<all_Final_Effi_StatErrorUp<<", all_Final_Effi_StatErrorDown : "<<all_Final_Effi_StatErrorDown<<std::endl;
    std::cout<<std::endl;

    return {
        effi_pos_gr,
        all_Final_Effi,
        all_Final_Effi_StatErrorUp,
        all_Final_Effi_StatErrorDown
    };

}

RelativeVariationSet GetVariation (
    std::vector<VariationSet> VariationSet_vec, // note : the first one has to be the baseline
    std::vector<std::string> file_title_vec,
    std::string leg_header_in,
    std::string final_output_directory,
    std::string plot_name_in
) {
    SetsPhenixStyle();
    TCanvas * c1 = new TCanvas("", "", 950, 800);
    c1 -> cd();

    TPad * pad1 = new TPad("", "", 0, 0.31, 0.85, 0.85);
    pad1 -> SetRightMargin(0.01);
    pad1 -> Draw();
    pad1 -> cd();

    TLegend * leg_variation = new TLegend(0.21,0.85,0.41,1);
    leg_variation -> SetBorderSize(0);
    leg_variation -> SetTextSize(0.03);
    leg_variation -> SetHeader(leg_header_in.c_str());

    TLine * line = new TLine();
    line -> SetLineStyle(2);
    line -> SetLineWidth(2);
    line -> SetLineColor(28);

    TLatex * ltx_draw = new TLatex();
    ltx_draw->SetNDC();
    ltx_draw->SetTextSize(0.03);

    TLatex * ltx_draw_all = new TLatex();
    ltx_draw_all->SetNDC();
    ltx_draw_all->SetTextSize(0.13);

    VariationSet_vec[0].effi_pos_gr -> SetMarkerStyle(20);
    VariationSet_vec[0].effi_pos_gr -> SetMarkerSize(1);
    VariationSet_vec[0].effi_pos_gr -> SetMarkerColor(1);
    VariationSet_vec[0].effi_pos_gr -> GetXaxis() -> SetTitle("Y axis [mm]");
    VariationSet_vec[0].effi_pos_gr -> GetYaxis() -> SetTitle("Hit Detection Efficiency (%)");
    VariationSet_vec[0].effi_pos_gr -> GetXaxis() -> SetLabelOffset(999);
    VariationSet_vec[0].effi_pos_gr -> GetXaxis() -> SetLabelSize(0);
    VariationSet_vec[0].effi_pos_gr -> GetXaxis() -> SetTitleOffset(999);

    // VariationSet_vec[0].effi_pos_gr -> GetYaxis() -> SetRangeUser(
    //     VariationSet_vec[0].effi_pos_gr -> GetYaxis() -> GetXmin() * 0.8,
    //     VariationSet_vec[0].effi_pos_gr -> GetYaxis() -> GetXmax() * 1.2
    // );

    VariationSet_vec[0].effi_pos_gr -> GetYaxis() -> SetRangeUser(98, 101);

    VariationSet_vec[0].effi_pos_gr -> GetXaxis() -> SetLimits(left_pos_edge, right_pos_edge);

    VariationSet_vec[0].effi_pos_gr -> Draw("AP");
    leg_variation -> AddEntry(VariationSet_vec[0].effi_pos_gr, file_title_vec[0].c_str(), "p");

    // note : i starts at 1, because the first one is the baseline
    for (int i = 1; i < VariationSet_vec.size(); i++) {
        VariationSet_vec[i].effi_pos_gr -> SetMarkerStyle(marker_code[i-1]);
        VariationSet_vec[i].effi_pos_gr -> SetMarkerSize(1);
        VariationSet_vec[i].effi_pos_gr -> SetMarkerColor(TColor::GetColor(color_code[i-1].c_str()));
        VariationSet_vec[i].effi_pos_gr -> SetLineColorAlpha(1,0);
        VariationSet_vec[i].effi_pos_gr -> Draw("p same");

        leg_variation -> AddEntry(VariationSet_vec[i].effi_pos_gr, file_title_vec[i].c_str(), "p");
    }

    line -> DrawLine(
        VariationSet_vec[0].effi_pos_gr->GetXaxis()->GetXmin(), 99., 
        VariationSet_vec[0].effi_pos_gr->GetXaxis()->GetXmax(), 99.
    );

    ltx_draw -> DrawLatex(0.21, 0.25, "*Statiscial error bars of varied data points are removed for better visualization");

    // Division : -----------------------------------------------------------------------------------------------------------------------------------------------------

    c1 -> cd();
    TPad * pad2 = new TPad("", "", 0, 0., 0.85, 0.38);
    pad2 -> SetRightMargin(0.01);
    pad2 -> Draw();
    pad2 -> cd();

    for (int i = 1; i < VariationSet_vec.size(); i++) {
        std::string plot_draw_text = (i == 1) ? "ap" : "p same";
        
        TGraph * temp_gr = GetRatioGr(
            VariationSet_vec[i].effi_pos_gr,
            VariationSet_vec[0].effi_pos_gr // note : baseline 
        );
        temp_gr -> GetYaxis() -> SetRangeUser(0.99,1.01);
        temp_gr -> GetXaxis() -> CenterTitle(true);
        temp_gr -> GetXaxis() -> SetTitle("Y axis [mm]");
        temp_gr -> GetYaxis() -> SetTitle("Ratio to baseline");
        temp_gr -> GetXaxis() -> SetLimits(left_pos_edge, right_pos_edge);
        temp_gr -> Draw(plot_draw_text.c_str());
    }

    line -> DrawLine(VariationSet_vec[0].effi_pos_gr->GetXaxis()->GetXmin(), 1, VariationSet_vec[0].effi_pos_gr->GetXaxis()->GetXmax(), 1);

    // Division : -----------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> cd();
    TPad * pad3 = new TPad("", "", 0.85, 0.31, 1., 0.85);
    pad3 -> SetLeftMargin(0.08);
    pad3 -> Draw();
    pad3 -> cd();

    TGraphAsymmErrors * all_baseline_gr = new TGraphAsymmErrors();
    all_baseline_gr -> SetPoint(0, 0, VariationSet_vec[0].all_Final_Effi);
    all_baseline_gr -> SetPointError(
        0, 
        0, 0, 
        VariationSet_vec[0].all_Final_Effi_StatErrorDown, VariationSet_vec[0].all_Final_Effi_StatErrorUp
    );
    all_baseline_gr -> SetMarkerStyle(20);
    all_baseline_gr -> SetMarkerSize(1);
    all_baseline_gr -> SetMarkerColor(1);
    all_baseline_gr -> GetXaxis() -> SetLabelOffset(999);
    all_baseline_gr -> GetXaxis() -> SetNdivisions(2);
    all_baseline_gr -> GetXaxis() -> SetTickLength(0.06);
    all_baseline_gr -> GetXaxis() -> SetLimits(left_pos_edge, right_pos_edge);

    all_baseline_gr -> GetYaxis() -> SetRangeUser(98, 101);
    all_baseline_gr -> GetYaxis() -> SetLabelSize(0);
    all_baseline_gr -> GetYaxis() -> SetTickLength(0.13);

    all_baseline_gr -> Draw("ap");

    // ltx_draw_all -> DrawLatex(0.2, 0.7, "*Combined");
    line -> DrawLine(
        all_baseline_gr->GetXaxis()->GetXmin(), 99., 
        all_baseline_gr->GetXaxis()->GetXmax(), 99.
    );

    for (int i = 1; i < VariationSet_vec.size(); i++) {
        TGraph * all_variation_gr = new TGraph();
        all_variation_gr -> SetPoint(0, 0, VariationSet_vec[i].all_Final_Effi);

        all_variation_gr -> SetMarkerStyle(marker_code[i-1]);
        all_variation_gr -> SetMarkerSize(1);
        all_variation_gr -> SetMarkerColor(TColor::GetColor(color_code[i-1].c_str()));
        all_variation_gr -> Draw("p same");
    }

    // Division : -----------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> cd();
    TPad * pad4 = new TPad("", "", 0.85, 0., 1., 0.38);
    pad4 -> SetLeftMargin(0.08);
    pad4 -> Draw();
    pad4 -> cd();

    for (int i = 1; i < VariationSet_vec.size(); i++) {
        std::string plot_draw_text = (i == 1) ? "ap" : "p same";

        TGraph * all_variation_ratio_gr = new TGraph();
        all_variation_ratio_gr -> SetPoint(0, 0, VariationSet_vec[i].all_Final_Effi / VariationSet_vec[0].all_Final_Effi);

        all_variation_ratio_gr -> SetMarkerStyle(marker_code[i-1]);
        all_variation_ratio_gr -> SetMarkerSize(1);
        all_variation_ratio_gr -> SetMarkerColor(TColor::GetColor(color_code[i-1].c_str()));

        all_variation_ratio_gr -> GetXaxis() -> SetLimits(left_pos_edge, right_pos_edge);
        all_variation_ratio_gr -> GetYaxis() -> SetRangeUser(0.99, 1.01);

        all_variation_ratio_gr -> GetXaxis() -> SetTitle("(Combined)");
        all_variation_ratio_gr -> GetXaxis() -> CenterTitle(true);
        all_variation_ratio_gr -> GetXaxis() -> SetTitleSize(0.11);
        all_variation_ratio_gr -> GetXaxis() -> SetTitleOffset(0.2);
        // all_variation_ratio_gr -> GetXaxis() -> SetLabelOffset(-0.04);
        all_variation_ratio_gr -> GetXaxis() -> SetLabelOffset(999);
        all_variation_ratio_gr -> GetXaxis() -> SetLabelSize(0.11);
        all_variation_ratio_gr -> GetXaxis() -> SetNdivisions(2);
        all_variation_ratio_gr -> GetXaxis() -> SetTickLength(0.04);

        all_variation_ratio_gr -> GetYaxis() -> SetLabelSize(0);
        all_variation_ratio_gr -> GetYaxis() -> SetTickLength(0.13);

        all_variation_ratio_gr -> Draw(plot_draw_text.c_str());
    }

    line -> DrawLine(
        all_baseline_gr->GetXaxis()->GetXmin(), 1, 
        all_baseline_gr->GetXaxis()->GetXmax(), 1
    );

    c1 -> cd();
    leg_variation -> Draw("same");
    c1 -> Print(Form("%s/%s.pdf", final_output_directory.c_str(), plot_name_in.c_str()));
    c1 -> Clear();

    // Division : -----------------------------------------------------------------------------------------------------------------------------------------------------

    TGraph * MaxVariationUp = new TGraph();
    TGraph * MaxVariationDown = new TGraph();

    for (int i = 0; i < VariationSet_vec.size(); i++) {
        if (i == 0) {
            for (int j = 0; j < VariationSet_vec[i].effi_pos_gr -> GetN(); j++) {
                double x = VariationSet_vec[i].effi_pos_gr -> GetPointX(j);
                double y = VariationSet_vec[i].effi_pos_gr -> GetPointY(j);
                MaxVariationUp->SetPoint(MaxVariationUp->GetN(), x, y);
                MaxVariationDown->SetPoint(MaxVariationDown->GetN(), x, y);
            }
        }

        else {
            for (int j = 0; j < VariationSet_vec[i].effi_pos_gr -> GetN(); j++) {
                double x = VariationSet_vec[i].effi_pos_gr -> GetPointX(j);
                double y = VariationSet_vec[i].effi_pos_gr -> GetPointY(j);

                if ( fabs(MaxVariationUp->GetPointX(j) - x) > 0.0001 ) {
                    std::cout<<"Error : the x value of the points are not the same"<<std::endl;
                    std::cout<<"x1 : "<<MaxVariationUp->GetPointX(j)<<", x2 : "<<x<<std::endl;
                    exit(1);
                }

                if (y > MaxVariationUp->GetPointY(j)) {
                    MaxVariationUp->SetPoint(j, x, y);
                }

                if (y < MaxVariationDown->GetPointY(j)) {
                    MaxVariationDown->SetPoint(j, x, y);
                }
            }
        }

    }

    for (int i = 0; i < MaxVariationUp->GetN(); i++) {
        MaxVariationUp->SetPoint(i, MaxVariationUp->GetPointX(i), fabs(MaxVariationUp->GetPointY(i) / VariationSet_vec[0].effi_pos_gr->GetPointY(i) - 1.0)); // note : for example, 0.1
        MaxVariationDown->SetPoint(i, MaxVariationDown->GetPointX(i), fabs(MaxVariationDown->GetPointY(i) / VariationSet_vec[0].effi_pos_gr->GetPointY(i) - 1.0)); // note : for example, 0.010
    }    

    // note : for the "all" case
    double all_MaxVariationUp   = VariationSet_vec[0].all_Final_Effi;
    double all_MaxVariationDown = VariationSet_vec[0].all_Final_Effi;

    for (int i = 1; i < VariationSet_vec.size(); i++) {
        if (VariationSet_vec[i].all_Final_Effi > all_MaxVariationUp) {
            all_MaxVariationUp = VariationSet_vec[i].all_Final_Effi;
        }
        if (VariationSet_vec[i].all_Final_Effi < all_MaxVariationDown) {
            all_MaxVariationDown = VariationSet_vec[i].all_Final_Effi;
        }
    }

    all_MaxVariationUp = fabs(all_MaxVariationUp / VariationSet_vec[0].all_Final_Effi - 1.0);
    all_MaxVariationDown = fabs(all_MaxVariationDown / VariationSet_vec[0].all_Final_Effi - 1.0);

    std::cout<<"leg_header_in : "<<leg_header_in<<std::endl;
    std::cout<<"all_MaxVariationUp: "<<all_MaxVariationUp<<std::endl;
    std::cout<<"all_MaxVariationDown: "<<all_MaxVariationDown<<std::endl;
    std::cout<<"MaxVariationUp in Y axis: ";
    for (int i = 0; i < MaxVariationUp->GetN(); i++) {
        std::cout<<Form("%.5f", MaxVariationUp->GetPointY(i))<<", ";
    }
    std::cout<<std::endl;

    std::cout<<"MaxVariationDown in Y axis: ";
    for (int i = 0; i < MaxVariationDown->GetN(); i++) {
        std::cout<<Form("%.5f", MaxVariationDown->GetPointY(i))<<", ";
    }
    std::cout<<std::endl;
    std::cout<<std::endl;

    return {
        MaxVariationUp,
        MaxVariationDown,

        all_MaxVariationUp,
        all_MaxVariationDown
    };

}

std::pair<double,double> CombineSyst(
    std::vector<RelativeVariationSet> RelativeVariationSet_vec,
    std::vector<std::string> Variation_title_vec,

    TGraphAsymmErrors * effi_pos_gr_baseline,
    std::vector<double> all_baseline_vec, // note : {efficiency, up_stat_error, down_stat_error}

    std::string output_directory
) {
    SetsPhenixStyle();
    TCanvas * c1 = new TCanvas("", "", 950, 800);
    c1 -> cd();

    TLegend * leg_errors = new TLegend(0.21,0.73,0.41,0.87);
    leg_errors -> SetBorderSize(0);
    leg_errors -> SetTextSize(0.03);

    TLegend * leg_errors_Low = new TLegend(0.21,0.19,0.41,0.33);
    leg_errors_Low -> SetBorderSize(0);
    leg_errors_Low -> SetTextSize(0.03);

    TLegend * leg_errors_Final = new TLegend(0.21,0.39,0.41,0.53);
    leg_errors_Final -> SetBorderSize(0);
    leg_errors_Final -> SetTextSize(0.03);

    TH1D * h1D_ErrorSum_Up   = new TH1D("h1D_ErrorSum_Up","h1D_ErrorSum_Up;Y axis [mm];Relative Variation",    N_bins_pos_hist,left_pos_edge,right_pos_edge);
    TH1D * h1D_ErrorSum_Down = new TH1D("h1D_ErrorSum_Down","h1D_ErrorSum_Down;Y axis [mm];Relative Variation",N_bins_pos_hist,left_pos_edge,right_pos_edge);
    h1D_ErrorSum_Up -> Reset("ICESM");
    h1D_ErrorSum_Down -> Reset("ICESM");

    std::vector<TGraph*> Variation_Ratio_up_vec; Variation_Ratio_up_vec.clear();
    std::vector<TGraph*> Variation_Ratio_down_vec; Variation_Ratio_down_vec.clear();

    TLine * line = new TLine();
    line -> SetLineStyle(2);
    line -> SetLineWidth(2);
    line -> SetLineColor(28);

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    // ltx->SetTextAlign(31);

    TLatex * ltx_draw = new TLatex();
    ltx_draw->SetNDC();
    ltx_draw->SetTextSize(0.04);

    double all_sum_MaxVariationUp = 0;
    double all_sum_MaxVariationDown = 0;

    leg_errors -> AddEntry(h1D_ErrorSum_Up,"Total Syst. High","l");
    leg_errors_Low -> AddEntry(h1D_ErrorSum_Down, "Total Syst. Low","l");    

    for (int i = 0; i < RelativeVariationSet_vec.size(); i++) {

        for (int j = 0; j < RelativeVariationSet_vec[i].MaxVariationUp -> GetN(); j++) {
            double x = RelativeVariationSet_vec[i].MaxVariationUp -> GetPointX(j);
            double y = RelativeVariationSet_vec[i].MaxVariationUp -> GetPointY(j);

            int X_bin_index = h1D_ErrorSum_Up -> GetXaxis() -> FindBin(x);

            h1D_ErrorSum_Up -> SetBinContent(
                X_bin_index, 
                h1D_ErrorSum_Up -> GetBinContent(X_bin_index) + pow(y, 2)
            );
        }

        for (int j = 0; j < RelativeVariationSet_vec[i].MaxVariationDown -> GetN(); j++) {
            double x = RelativeVariationSet_vec[i].MaxVariationDown -> GetPointX(j);
            double y = RelativeVariationSet_vec[i].MaxVariationDown -> GetPointY(j);

            int X_bin_index = h1D_ErrorSum_Down -> GetXaxis() -> FindBin(x);

            h1D_ErrorSum_Down -> SetBinContent(
                X_bin_index, 
                h1D_ErrorSum_Down -> GetBinContent(X_bin_index) + pow(y, 2)
            );
        }

        // Division : --------------------------------------------------------------------------------
        Variation_Ratio_up_vec.push_back(
            (TGraph*) RelativeVariationSet_vec[i].MaxVariationUp -> Clone(Form("Variation_Ratio_up_%d", i))
        );
        Variation_Ratio_down_vec.push_back(
            (TGraph*) RelativeVariationSet_vec[i].MaxVariationDown -> Clone(Form("Variation_Ratio_down_%d", i))
        );

        Variation_Ratio_up_vec.back() -> SetMarkerStyle(marker_code[i]);
        Variation_Ratio_up_vec.back() -> SetMarkerSize(1);
        Variation_Ratio_up_vec.back() -> SetMarkerColor(TColor::GetColor(syst_HighLow_color_code[i*2].c_str()));

        Variation_Ratio_down_vec.back() -> SetMarkerStyle(marker_code[i]);
        Variation_Ratio_down_vec.back() -> SetMarkerSize(1);
        Variation_Ratio_down_vec.back() -> SetMarkerColor(TColor::GetColor(syst_HighLow_color_code[i*2+1].c_str()));

        leg_errors -> AddEntry(
            Variation_Ratio_up_vec.back(), 
            Form("Syst. High, %s", Variation_title_vec[i].c_str()), 
            "p"
        );

        leg_errors_Low -> AddEntry(
            Variation_Ratio_down_vec.back(), 
            Form("Syst. Low, %s", Variation_title_vec[i].c_str()), 
            "p"
        );

        for (int j = 0; j < Variation_Ratio_up_vec.back() -> GetN(); j++) {
            Variation_Ratio_up_vec.back() -> SetPoint(
                j, 
                Variation_Ratio_up_vec.back() -> GetPointX(j), 
                Variation_Ratio_up_vec.back() -> GetPointY(j) + 1.0
            );   

            Variation_Ratio_down_vec.back() -> SetPoint(
                j, 
                Variation_Ratio_down_vec.back() -> GetPointX(j), 
                1.0 - Variation_Ratio_down_vec.back() -> GetPointY(j)
            );
        }

        // Division : --------------------------------------------------------------------------------

        all_sum_MaxVariationUp += pow(RelativeVariationSet_vec[i].all_MaxVariationUp, 2.);
        all_sum_MaxVariationDown += pow(RelativeVariationSet_vec[i].all_MaxVariationDown, 2.);
    }

    all_sum_MaxVariationUp = sqrt(all_sum_MaxVariationUp);
    all_sum_MaxVariationDown = sqrt(all_sum_MaxVariationDown);

    for (int i = 1; i <= h1D_ErrorSum_Up -> GetNbinsX(); i++) {
        double bin_content_up   = h1D_ErrorSum_Up   -> GetBinContent(i);
        double bin_content_down = h1D_ErrorSum_Down -> GetBinContent(i);

        if (bin_content_up > 0) {
            bin_content_up   = sqrt(bin_content_up);
        }
        else {
            bin_content_up   = 0;
        }

        if (bin_content_down > 0) {
            bin_content_down = sqrt(bin_content_down);
        }
        else {
            bin_content_down = 0;
        }

        h1D_ErrorSum_Up   -> SetBinContent(i, bin_content_up);
        h1D_ErrorSum_Down -> SetBinContent(i, bin_content_down);
    }

    
    TH1D * h1D_ErrorSum_Ratio_Up = (TH1D*) h1D_ErrorSum_Up -> Clone("h1D_ErrorSum_Ratio_Up");
    TH1D * h1D_ErrorSum_Ratio_Down = (TH1D*) h1D_ErrorSum_Down -> Clone("h1D_ErrorSum_Ratio_Down");

    h1D_ErrorSum_Ratio_Up -> SetLineColor(1);
    h1D_ErrorSum_Ratio_Down -> SetLineColor(1);

    h1D_ErrorSum_Ratio_Up -> GetXaxis() -> CenterTitle(true);
    h1D_ErrorSum_Ratio_Down -> GetXaxis() -> CenterTitle(true);

    h1D_ErrorSum_Ratio_Up -> GetYaxis() -> SetTitle("Variation / Baseline");
    h1D_ErrorSum_Ratio_Down -> GetYaxis() -> SetTitle("Variation / Baseline");

    h1D_ErrorSum_Ratio_Up -> GetYaxis() -> SetRangeUser(0.99, 1.01);
    h1D_ErrorSum_Ratio_Up -> GetYaxis() -> SetTitleOffset(2.0);

    for (int i = 1; i <= h1D_ErrorSum_Ratio_Up -> GetNbinsX(); i++) {
        double y_up = h1D_ErrorSum_Ratio_Up -> GetBinContent(i);
        double y_down = h1D_ErrorSum_Ratio_Down -> GetBinContent(i);

        h1D_ErrorSum_Ratio_Up -> SetBinContent(i, y_up + 1);
        h1D_ErrorSum_Ratio_Down -> SetBinContent(i, 1 - y_down);
    }

    c1 -> cd();
    TPad * pad1 = new TPad("", "", 0, 0.0, 0.85, 1);
    pad1 -> SetRightMargin(0.025);
    pad1 -> SetLeftMargin(0.2);
    pad1 -> Draw();
    pad1 -> cd();

    h1D_ErrorSum_Ratio_Up -> Draw("hist");
    h1D_ErrorSum_Ratio_Down -> Draw("hist same");

    for (int i = 0; i < Variation_Ratio_up_vec.size(); i++) {
        Variation_Ratio_up_vec[i] -> Draw("p same");
        Variation_Ratio_down_vec[i] -> Draw("p same");
    }
    leg_errors -> Draw("same");
    leg_errors_Low -> Draw("same");
    line -> DrawLine(
        h1D_ErrorSum_Ratio_Up -> GetXaxis() -> GetXmin(), 1, 
        h1D_ErrorSum_Ratio_Up -> GetXaxis() -> GetXmax(), 1
    );
    ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));


    c1 -> cd();
    TPad * pad3 = new TPad("", "", 0.85, 0.0, 1., 1.0);
    pad3 -> SetLeftMargin(0.00);
    // pad3 -> SetRightMargin(0.00);
    pad3 -> Draw();
    pad3 -> cd();

    TH1D * all_h1D_ErrorSum_Ratio_Up = new TH1D("all_h1D_ErrorSum_Ratio_Up","all_h1D_ErrorSum_Ratio_Up",1,left_pos_edge,right_pos_edge);
    TH1D * all_h1D_ErrorSum_Ratio_Down = new TH1D("all_h1D_ErrorSum_Ratio_Down","all_h1D_ErrorSum_Ratio_Down",1,left_pos_edge,right_pos_edge);

    all_h1D_ErrorSum_Ratio_Up -> SetBinContent(1, all_sum_MaxVariationUp + 1);
    all_h1D_ErrorSum_Ratio_Up -> SetLineColor(1);
    all_h1D_ErrorSum_Ratio_Down -> SetBinContent(1, 1 - all_sum_MaxVariationDown);
    all_h1D_ErrorSum_Ratio_Down -> SetLineColor(1);

    all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> SetNdivisions(2);
    all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> SetLabelOffset(999);
    all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> SetLabelSize(0.2);
    all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> SetTickLength(0.03);
    all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> SetTitle("(Combined)");
    all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> SetTitleSize(0.17);
    all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> SetTitleOffset(0.1);
    all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> CenterTitle(true);



    all_h1D_ErrorSum_Ratio_Up -> GetYaxis() -> SetLabelSize(0);
    all_h1D_ErrorSum_Ratio_Up -> GetYaxis() -> SetTickLength(0.13);
    all_h1D_ErrorSum_Ratio_Up -> GetYaxis() -> SetRangeUser(0.99, 1.01);

    all_h1D_ErrorSum_Ratio_Up -> Draw("hist");
    all_h1D_ErrorSum_Ratio_Down -> Draw("hist same");

    for (int i = 0; i < RelativeVariationSet_vec.size(); i++){
        TGraph * all_variation_Ratio_up_gr = new TGraph();
        TGraph * all_variation_Ratio_down_gr = new TGraph();

        all_variation_Ratio_up_gr -> SetPoint(0, 0, RelativeVariationSet_vec[i].all_MaxVariationUp + 1);
        all_variation_Ratio_down_gr -> SetPoint(0, 0, 1 - RelativeVariationSet_vec[i].all_MaxVariationDown);

        all_variation_Ratio_up_gr -> SetMarkerStyle(marker_code[i]);
        all_variation_Ratio_up_gr -> SetMarkerSize(1);
        all_variation_Ratio_up_gr -> SetMarkerColor(TColor::GetColor(syst_HighLow_color_code[i*2].c_str()));

        all_variation_Ratio_down_gr -> SetMarkerStyle(marker_code[i]);
        all_variation_Ratio_down_gr -> SetMarkerSize(1);
        all_variation_Ratio_down_gr -> SetMarkerColor(TColor::GetColor(syst_HighLow_color_code[i*2+1].c_str()));

        all_variation_Ratio_up_gr -> Draw("p same");
        all_variation_Ratio_down_gr -> Draw("p same");
    }

    line -> DrawLine(
        all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> GetXmin(), 1, 
        all_h1D_ErrorSum_Ratio_Up -> GetXaxis() -> GetXmax(), 1
    );


    c1 -> Print(Form("%s/SummaryRatio_VariationToBaseline.pdf", output_directory.c_str()));
    c1 -> Clear();

    // Division : -----------------------------------------------------------------------------------------------------------------------------------------------------

    TLine * line_99 = new TLine(left_pos_edge,99,right_pos_edge,99);
    line_99 -> SetLineColor(TColor::GetColor("#941100"));
    line_99 -> SetLineWidth(5);
    line_99 -> SetLineStyle(7);

    effi_pos_gr_baseline -> SetMarkerStyle(20);
    effi_pos_gr_baseline -> SetMarkerSize(1);
    effi_pos_gr_baseline -> SetMarkerColor(1);
    effi_pos_gr_baseline -> GetXaxis() -> CenterTitle(true);
    effi_pos_gr_baseline -> GetXaxis() -> SetTitle("Y axis [mm]");
    effi_pos_gr_baseline -> GetYaxis() -> SetTitle("Hit Detection Efficiency (%)");
    effi_pos_gr_baseline -> GetYaxis() -> SetRangeUser(95, 102);
    effi_pos_gr_baseline -> GetXaxis() -> SetLimits(-10, 10);

    TGraphAsymmErrors * effi_pos_gr_syst = (TGraphAsymmErrors*) effi_pos_gr_baseline -> Clone("effi_pos_gr_syst");
    effi_pos_gr_syst -> SetMarkerStyle(20);
    effi_pos_gr_syst -> SetMarkerSize(0);
    effi_pos_gr_syst -> SetMarkerColorAlpha(1,0);
    effi_pos_gr_syst -> SetFillColorAlpha(1,0.5);
    effi_pos_gr_syst -> SetFillColorAlpha(92, 0.85);

    for (int i = 0; i < effi_pos_gr_syst -> GetN(); i++) {
        double x = effi_pos_gr_syst -> GetPointX(i);
        int x_bin_index = h1D_ErrorSum_Up -> GetXaxis() -> FindBin(x);

        double y = effi_pos_gr_syst -> GetPointY(i);
        double y_error_up = h1D_ErrorSum_Up -> GetBinContent(x_bin_index) * y;
        double y_error_down = h1D_ErrorSum_Down -> GetBinContent(x_bin_index) * y;

        // todo : no check for y_error_down 
        y_error_up = (y + y_error_up > 100) ? 100 - y : y_error_up;
        // y_error_down = (y - y_error_down > 100) ? 100 : y_error_down;

        effi_pos_gr_syst -> SetPointError(
            i, 
            h1D_ErrorSum_Down->GetBinWidth(x_bin_index)/2., h1D_ErrorSum_Down->GetBinWidth(x_bin_index)/2., 
            y_error_down, y_error_up
        );
    }
    
    TGraphErrors * gr_syst_template = new TGraphErrors();
    gr_syst_template->SetPoint(0,0,0);
    gr_syst_template->SetMarkerStyle(20);
    gr_syst_template->SetMarkerSize(0);
    gr_syst_template->SetMarkerColorAlpha(1,0);
    gr_syst_template->GetXaxis()->SetTitle("Y axis [mm]");
    gr_syst_template->GetYaxis()->SetTitle("Hit Detection Efficiency (%)");
    gr_syst_template->GetXaxis()->CenterTitle(true);
    gr_syst_template->GetYaxis()->CenterTitle(true);
    gr_syst_template->GetYaxis()->SetRangeUser(95, 102);
    gr_syst_template->GetXaxis()->SetLimits(left_pos_edge, right_pos_edge);

    TGraphErrors * gr_syst_style = new TGraphErrors();
    gr_syst_style->SetMarkerStyle(effi_pos_gr_baseline->GetMarkerStyle());
    gr_syst_style->SetMarkerSize(effi_pos_gr_baseline->GetMarkerSize());
    gr_syst_style->SetMarkerColor(effi_pos_gr_baseline->GetMarkerColor());
    gr_syst_style->SetLineWidth(0);
    gr_syst_style->SetLineColorAlpha(1,0);
    gr_syst_style->SetFillColorAlpha(effi_pos_gr_syst->GetFillColor(), 0.85);


    leg_errors_Final -> AddEntry(
        gr_syst_style, 
        "Data", 
        "epf"
    );
    leg_errors_Final -> AddEntry(
        line_99,
        "99% efficiency line",
        "l"
    );

    c1 -> cd();
    pad1 = new TPad("", "", 0, 0.0, 0.85, 1);
    pad1 -> SetRightMargin(0.026);
    pad1 -> Draw();
    pad1 -> cd();

    gr_syst_template -> Draw("AP");
    effi_pos_gr_syst -> Draw("P2 same");
    // effi_pos_gr_syst -> Draw("P2");
    effi_pos_gr_baseline -> Draw("ep same");

    line_99 -> Draw("same");
    ltx->DrawLatex(0.38, 0.88, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));
    ltx_draw->DrawLatex(0.44, 0.21, "*Beam spot region shown only");

    leg_errors_Final -> Draw("same");

    ltx_draw->SetTextSize(0.025);
    ltx_draw->DrawLatex(0.2, 0.3, Form("INTT Ladder Hit Detection Effi. (Comb.): %.2f ^{+%.2f}_{-%.2f} (stat.) ^{+%.2f}_{-%.2f} (syst.) %%",
            all_baseline_vec[0], 
            all_baseline_vec[1], all_baseline_vec[2], 
            all_sum_MaxVariationUp * all_baseline_vec[0], all_sum_MaxVariationDown * all_baseline_vec[0]
        )
    );

    c1 -> cd();
    pad3 = new TPad("", "", 0.85, 0.0, 1., 1.0);
    pad3 -> SetLeftMargin(0.00);
    // pad3 -> SetRightMargin(0.00);
    pad3 -> Draw();
    pad3 -> cd();

    TGraphAsymmErrors * all_effi_baseline_gr = new TGraphAsymmErrors();
    all_effi_baseline_gr -> SetPoint(0, 0, all_baseline_vec[0]);
    all_effi_baseline_gr -> SetPointError(
        0, 
        5, 5, 
        all_baseline_vec[2], all_baseline_vec[1]
    );

    all_effi_baseline_gr -> SetMarkerStyle(20);
    all_effi_baseline_gr -> SetMarkerSize(1);
    all_effi_baseline_gr -> SetMarkerColor(1);
    
    all_effi_baseline_gr -> GetXaxis() -> CenterTitle(true);
    all_effi_baseline_gr -> GetXaxis() -> SetTitle("(Combined)");
    all_effi_baseline_gr -> GetXaxis() -> SetTitleSize(0.17);
    all_effi_baseline_gr -> GetXaxis() -> SetTitleOffset(0.1);
    all_effi_baseline_gr -> GetXaxis() -> SetNdivisions(2);
    all_effi_baseline_gr -> GetXaxis() -> SetLabelOffset(999);
    all_effi_baseline_gr -> GetXaxis() -> SetTickLength(0.03);
    all_effi_baseline_gr -> GetXaxis() -> SetLimits(-10, 10);

    all_effi_baseline_gr -> GetYaxis() -> SetTickLength(0.13);
    all_effi_baseline_gr -> GetYaxis() -> SetLabelSize(0);
    all_effi_baseline_gr -> GetYaxis() -> SetTitleSize(0);
    all_effi_baseline_gr -> GetYaxis() -> SetRangeUser(95, 102);
    

    TGraphAsymmErrors * all_effi_syst_gr = new TGraphAsymmErrors();
    all_effi_syst_gr -> SetPoint(0, 0, all_baseline_vec[0]);
    all_effi_syst_gr -> SetPointError(
        0, 
        5, 5, 
        all_sum_MaxVariationDown * all_baseline_vec[0], all_sum_MaxVariationUp * all_baseline_vec[0]
    );
    all_effi_syst_gr -> SetMarkerStyle(20);
    all_effi_syst_gr -> SetMarkerSize(0);
    all_effi_syst_gr -> SetMarkerColorAlpha(1,0);
    all_effi_syst_gr -> SetFillColorAlpha(1,0.5);
    all_effi_syst_gr -> SetFillColorAlpha(92, 0.85);

    all_effi_baseline_gr -> Draw("AP");
    all_effi_syst_gr -> Draw("P2 same");
    all_effi_baseline_gr -> Draw("ep same");

    line_99 -> DrawLine(
        all_effi_baseline_gr->GetXaxis()->GetXmin(), 99, 
        all_effi_baseline_gr->GetXaxis()->GetXmax(), 99
    );

    c1 -> Print(Form("%s/Final_Detection_Efficiency.pdf", output_directory.c_str()));

    return std::make_pair(
        all_sum_MaxVariationUp,
        all_sum_MaxVariationDown
    );
}

// note : Run52 single 
int FinalResult_Run52() {
    std::string mother_directory_Run52 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/DetectionEffi"; 
    std::string final_output_directory = mother_directory_Run52 + "/FinalResult";
    system(Form("mkdir -p %s", final_output_directory.c_str()));

    VariationSet BaseLine_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_BaseLine.root", mother_directory_Run52.c_str()),
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VaryColumn_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column9_VaryColumn.root", mother_directory_Run52.c_str()),
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VaryPosCutLarge_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VaryPosCutLarge.root", mother_directory_Run52.c_str()),
        }
    );

    VariationSet VaryPosCutSmall_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VaryPosCutSmall.root", mother_directory_Run52.c_str()),
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VarySlopeCutLarge_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VarySlopeCutLarge.root", mother_directory_Run52.c_str()),
        }
    );

    VariationSet VarySlopeCutSmall_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VarySlopeCutSmall.root", mother_directory_Run52.c_str()),
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VaryBoundaryCutLarge_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VaryBoundaryCutLarge.root", mother_directory_Run52.c_str()),
        }
    );

    VariationSet VaryBoundaryCutSmall_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VaryBoundaryCutSmall.root", mother_directory_Run52.c_str()),
        }
    );    

    // Division : ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    RelativeVariationSet VaryColumn_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VaryColumn_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (Column 8)",
            "Column 9"
        }, // note : file_title_vec,
        "Column variation, Run52", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VaryColumn"    // note : plot_name_in
    );

    RelativeVariationSet VaryPosCut_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VaryPosCutLarge_EffiSet,
            VaryPosCutSmall_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (|Residual| < 0.5 mm)",
            "|Residual| < 0.7 mm",
            "|Residual| < 0.3 mm"
        }, // note : file_title_vec,
        "Residual variation, Run52", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VaryPosCut"    // note : plot_name_in
    );

    RelativeVariationSet VarySlopeCut_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VarySlopeCutLarge_EffiSet,
            VarySlopeCutSmall_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (|slope| < 0.01)",
            "|slope| < 0.013",
            "|slope| < 0.007"
        }, // note : file_title_vec,
        "Slope cut variation, Run52", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VarySlopeCut"    // note : plot_name_in
    );

    RelativeVariationSet VaryBoundaryCut_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VaryBoundaryCutLarge_EffiSet,
            VaryBoundaryCutSmall_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (Hits of L0 & L2 to edge < 8 channels)",
            "Hits of L0 & L2 to edge < 11 channels",
            "Hits of L0 & L2 to edge < 5 channels"
        }, // note : file_title_vec,
        "Boundary cut variation, Run52", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VaryBoundaryCut"    // note : plot_name_in
    );

    // Division : ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    std::pair<double,double> all_syst_unc_pair = CombineSyst(
        {
            VaryColumn_variation,
            VaryPosCut_variation,
            VarySlopeCut_variation,
            VaryBoundaryCut_variation
        }, // note : std::vector<RelativeVariationSet> RelativeVariationSet_vec
        
        {
            "Column variation",
            "Residual cut variation",
            "Slope cut variation",
            "Boundary cut variation"
        }, // note : std::vector<std::string> Variation_title_vec
         
        BaseLine_EffiSet.effi_pos_gr, // note : TGraphAsymmErrors * effi_pos_gr_baseline

        {
            BaseLine_EffiSet.all_Final_Effi,
            BaseLine_EffiSet.all_Final_Effi_StatErrorUp,
            BaseLine_EffiSet.all_Final_Effi_StatErrorDown
        }, // note : {efficiency, up_stat_error, down_stat_error}

        final_output_directory // note : std::string output_directory
    );

    std::cout<<"Final Result : "<<std::endl;
    std::cout<<"Detection Efficiency : "<<BaseLine_EffiSet.all_Final_Effi<<std::endl;
    std::cout<<"Statistical Uncertainty (Up) : "<<BaseLine_EffiSet.all_Final_Effi_StatErrorUp<<" - (Down): "<<BaseLine_EffiSet.all_Final_Effi_StatErrorDown<<std::endl;
    std::cout<<"Systematic Uncertainty  (Up): "<<all_syst_unc_pair.first * BaseLine_EffiSet.all_Final_Effi<<" - (Down): "<<all_syst_unc_pair.second * BaseLine_EffiSet.all_Final_Effi<<std::endl;

    return 888;
}

// note : Run89 single 
int FinalResult_Run89() {
    std::string mother_directory_Run89 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run89/DetectionEffi"; 
    std::string final_output_directory = mother_directory_Run89 + "/FinalResult";
    system(Form("mkdir -p %s", final_output_directory.c_str()));

    VariationSet BaseLine_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run89_Column10_BaseLine.root", mother_directory_Run89.c_str()),
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VaryColumn_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run89_Column11_VaryColumn.root", mother_directory_Run89.c_str()),
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VaryPosCutLarge_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run89_Column10_VaryPosCutLarge.root", mother_directory_Run89.c_str()),
        }
    );

    VariationSet VaryPosCutSmall_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run89_Column10_VaryPosCutSmall.root", mother_directory_Run89.c_str()),
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VarySlopeCutLarge_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run89_Column10_VarySlopeCutLarge.root", mother_directory_Run89.c_str()),
        }
    );

    VariationSet VarySlopeCutSmall_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run89_Column10_VarySlopeCutSmall.root", mother_directory_Run89.c_str()),
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VaryBoundaryCutLarge_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run89_Column10_VaryBoundaryCutLarge.root", mother_directory_Run89.c_str()),
        }
    );

    VariationSet VaryBoundaryCutSmall_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run89_Column10_VaryBoundaryCutSmall.root", mother_directory_Run89.c_str()),
        }
    );    

    // Division : ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    RelativeVariationSet VaryColumn_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VaryColumn_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (Column 10)",
            "Column 11"
        }, // note : file_title_vec,
        "Column variation, Run89", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VaryColumn"    // note : plot_name_in
    );

    RelativeVariationSet VaryPosCut_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VaryPosCutLarge_EffiSet,
            VaryPosCutSmall_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (|Residual| < 0.5 mm)",
            "|Residual| < 0.7 mm",
            "|Residual| < 0.3 mm"
        }, // note : file_title_vec,
        "Residual variation, Run89", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VaryPosCut"    // note : plot_name_in
    );

    RelativeVariationSet VarySlopeCut_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VarySlopeCutLarge_EffiSet,
            VarySlopeCutSmall_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (|slope| < 0.01)",
            "|slope| < 0.013",
            "|slope| < 0.007"
        }, // note : file_title_vec,
        "Slope cut variation, Run89", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VarySlopeCut"    // note : plot_name_in
    );

    RelativeVariationSet VaryBoundaryCut_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VaryBoundaryCutLarge_EffiSet,
            VaryBoundaryCutSmall_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (Hits of L0 & L2 to edge < 8 channels)",
            "Hits of L0 & L2 to edge < 11 channels",
            "Hits of L0 & L2 to edge < 5 channels"
        }, // note : file_title_vec,
        "Boundary cut variation, Run89", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VaryBoundaryCut"    // note : plot_name_in
    );

    // Division : ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    std::pair<double,double> all_syst_unc_pair = CombineSyst(
        {
            VaryColumn_variation,
            VaryPosCut_variation,
            VarySlopeCut_variation,
            VaryBoundaryCut_variation
        }, // note : std::vector<RelativeVariationSet> RelativeVariationSet_vec
        
        {
            "Column variation",
            "Residual cut variation",
            "Slope cut variation",
            "Boundary cut variation"
        }, // note : std::vector<std::string> Variation_title_vec
         
        BaseLine_EffiSet.effi_pos_gr, // note : TGraphAsymmErrors * effi_pos_gr_baseline

        {
            BaseLine_EffiSet.all_Final_Effi,
            BaseLine_EffiSet.all_Final_Effi_StatErrorUp,
            BaseLine_EffiSet.all_Final_Effi_StatErrorDown
        }, // note : {efficiency, up_stat_error, down_stat_error}

        final_output_directory // note : std::string output_directory
    );

    std::cout<<"Final Result : "<<std::endl;
    std::cout<<"Detection Efficiency : "<<BaseLine_EffiSet.all_Final_Effi<<std::endl;
    std::cout<<"Statistical Uncertainty (Up) : "<<BaseLine_EffiSet.all_Final_Effi_StatErrorUp<<" - (Down): "<<BaseLine_EffiSet.all_Final_Effi_StatErrorDown<<std::endl;
    std::cout<<"Systematic Uncertainty (Up)  : "<<all_syst_unc_pair.first * BaseLine_EffiSet.all_Final_Effi<<" - (Down): "<<all_syst_unc_pair.second * BaseLine_EffiSet.all_Final_Effi<<std::endl;

    return 888;
}


// note : Run52 + Run89 combined 
int FinalResult_Run52Run89Comb() {
    std::string mother_directory_Run52 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/DetectionEffi";
    std::string mother_directory_Run89 = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run89/DetectionEffi";
    std::string final_output_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52Run89Comb/DetectionEffi";
    system(Form("mkdir -p %s", final_output_directory.c_str()));

    VariationSet BaseLine_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_BaseLine.root", mother_directory_Run52.c_str()),
            Form("%s/Data_DetectionEffi_Run89_Column10_BaseLine.root", mother_directory_Run89.c_str())
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VaryColumn_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column9_VaryColumn.root", mother_directory_Run52.c_str()),
            Form("%s/Data_DetectionEffi_Run89_Column11_VaryColumn.root", mother_directory_Run89.c_str())
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VaryPosCutLarge_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VaryPosCutLarge.root", mother_directory_Run52.c_str()),
            Form("%s/Data_DetectionEffi_Run89_Column10_VaryPosCutLarge.root", mother_directory_Run89.c_str())
        }
    );

    VariationSet VaryPosCutSmall_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VaryPosCutSmall.root", mother_directory_Run52.c_str()),
            Form("%s/Data_DetectionEffi_Run89_Column10_VaryPosCutSmall.root", mother_directory_Run89.c_str())
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VarySlopeCutLarge_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VarySlopeCutLarge.root", mother_directory_Run52.c_str()),
            Form("%s/Data_DetectionEffi_Run89_Column10_VarySlopeCutLarge.root", mother_directory_Run89.c_str())
        }
    );

    VariationSet VarySlopeCutSmall_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VarySlopeCutSmall.root", mother_directory_Run52.c_str()),
            Form("%s/Data_DetectionEffi_Run89_Column10_VarySlopeCutSmall.root", mother_directory_Run89.c_str())
        }
    );

    // Division : -----------------------------------------------------------------------------------------------------
    VariationSet VaryBoundaryCutLarge_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VaryBoundaryCutLarge.root", mother_directory_Run52.c_str()),
            Form("%s/Data_DetectionEffi_Run89_Column10_VaryBoundaryCutLarge.root", mother_directory_Run89.c_str())
        }
    );

    VariationSet VaryBoundaryCutSmall_EffiSet = MakeEffiPos(
        {
            Form("%s/Data_DetectionEffi_Run52_Column8_VaryBoundaryCutSmall.root", mother_directory_Run52.c_str()),
            Form("%s/Data_DetectionEffi_Run89_Column10_VaryBoundaryCutSmall.root", mother_directory_Run89.c_str())
        }
    );    

    // Division : ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    RelativeVariationSet VaryColumn_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VaryColumn_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (Run52: Column 8, Run89: Column 10)",
            "Run52: Column 9, Run89: Column 11"
        }, // note : file_title_vec,
        "Column variation, Run52+Run89", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VaryColumn"    // note : plot_name_in
    );

    RelativeVariationSet VaryPosCut_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VaryPosCutLarge_EffiSet,
            VaryPosCutSmall_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (|Residual| < 0.5 mm)",
            "|Residual| < 0.7 mm",
            "|Residual| < 0.3 mm"
        }, // note : file_title_vec,
        "Residual variation, Run52+Run89", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VaryPosCut"    // note : plot_name_in
    );

    RelativeVariationSet VarySlopeCut_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VarySlopeCutLarge_EffiSet,
            VarySlopeCutSmall_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (|slope| < 0.01)",
            "|slope| < 0.013",
            "|slope| < 0.007"
        }, // note : file_title_vec,
        "Slope cut variation, Run52+Run89", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VarySlopeCut"    // note : plot_name_in
    );

    RelativeVariationSet VaryBoundaryCut_variation = GetVariation(
        {
            BaseLine_EffiSet,
            VaryBoundaryCutLarge_EffiSet,
            VaryBoundaryCutSmall_EffiSet
        }, // note : the first one has to be the baseline
        {
            "Baseline (Hits of L0 & L2 to edge < 8 channels)",
            "Hits of L0 & L2 to edge < 11 channels",
            "Hits of L0 & L2 to edge < 5 channels"
        }, // note : file_title_vec,
        "Boundary cut variation, Run52+Run89", // note : leg_header_in,
        final_output_directory,   // note : final_output_directory,
        "Variation_VaryBoundaryCut"    // note : plot_name_in
    );

    // Division : ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    std::pair<double,double> all_syst_unc_pair = CombineSyst(
        {
            VaryColumn_variation,
            VaryPosCut_variation,
            VarySlopeCut_variation,
            VaryBoundaryCut_variation
        }, // note : std::vector<RelativeVariationSet> RelativeVariationSet_vec
        
        {
            "Column variation",
            "Residual cut variation",
            "Slope cut variation",
            "Boundary cut variation"
        }, // note : std::vector<std::string> Variation_title_vec
         
        BaseLine_EffiSet.effi_pos_gr, // note : TGraphAsymmErrors * effi_pos_gr_baseline

        {
            BaseLine_EffiSet.all_Final_Effi,
            BaseLine_EffiSet.all_Final_Effi_StatErrorUp,
            BaseLine_EffiSet.all_Final_Effi_StatErrorDown
        }, // note : {efficiency, up_stat_error, down_stat_error}

        final_output_directory // note : std::string output_directory
    );

    std::cout<<"Final Result : "<<std::endl;
    std::cout<<"Detection Efficiency : "<<BaseLine_EffiSet.all_Final_Effi<<std::endl;
    std::cout<<"Statistical Uncertainty (Up) : "<<BaseLine_EffiSet.all_Final_Effi_StatErrorUp<<" -  (Down): "<<BaseLine_EffiSet.all_Final_Effi_StatErrorDown<<std::endl;
    std::cout<<"Systematic Uncertainty (Up)  : "<<all_syst_unc_pair.first * BaseLine_EffiSet.all_Final_Effi<<" -  (Down): "<<all_syst_unc_pair.second * BaseLine_EffiSet.all_Final_Effi<<std::endl;

    return 888;
}