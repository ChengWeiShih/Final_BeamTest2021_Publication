#ifndef NTRACKANA_H
#define NTRACKANA_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <TGraph.h>

// #include "sPhenixStyle.h"

class NTrackAna{
    public:
        NTrackAna(
            bool isData_in,
            int runnumber_in,
            std::string input_directory_in,
            std::string input_filename_in,
            std::string output_directory_in,
            std::string output_file_name_suffix_in
        );

        void SetBeamSpotColumn(int column_in) {
            BeamSpotColumn = column_in;
        };

        void SetSkipColumn(std::vector<int> skip_column_in) {
            skip_column = skip_column_in;
        };

        std::string GetOutputFileName(){
            return output_filename;
        };

        void Get_l1_alignment(); // note : you may just do it twice in one run, before and after, if the input correction is nan

        void GetNTrack(double  noise_hit_distance);

        TH1D * Get_h1D_Column_NTrack() {
            return h1D_Column_NTrack;
        };

        TH1D * Get_h1D_NTrack_NoZero() {
            return h1D_NTrack_NoZero;
        };

        TH1D * Get_h1D_NTrackWithZero() {
            return h1D_NTrackWithZero;
        };

        TH1D * Get_h1D_l0_ClusPosFitDiff() {
            return h1D_l0_ClusPosFitDiff;
        };

        TH1D * Get_h1D_l1_ClusPosFitDiff() {
            return h1D_l1_ClusPosFitDiff;
        };

        TH1D * Get_h1D_l2_ClusPosFitDiff() {
            return h1D_l2_ClusPosFitDiff;
        };

        TH1D * Get_h1D_l1_residual() {
            return h1D_l1_residual;
        };


        void EndRun();

    protected:
        // Division : --------------------------------------------------------------------------------
        struct ClusInfo{
            int adc;
            int size;
            double pos;
        };

        // Division : --------------------------------------------------------------------------------
        bool isData;
        int runnumber;
        std::string input_directory;
        std::string input_filename;
        std::string output_directory;
        std::string output_file_name_suffix;

        std::string output_filename;
        std::vector<double> l1_alignment_correction_vec;
        std::vector<double> l1_alignment_correction_error_vec;
        std::vector<double> l1_alignment_correction_final_vec;

        int BeamSpotColumn; // note : 1 to 13;
        std::vector<int> skip_column; // note : 1 to 13;

        TFile * file_in;
        TTree * tree_in;
        void PrepareInputData();

        int eID;
        int DSE;
        std::vector<int>* layer;
        std::vector<int>* chip;
        std::vector<int>* Nhit;
        std::vector<double>* position;
        std::vector<double>* cluster_adc;

        std::vector<std::vector<std::vector<std::vector<ClusInfo>>>> cluster_info_NoL1Aligned_vec; // note : [event][chip][layer][clusters]
        TH1D * Sub_Get_l1_alignment_hist(int Selected_Column); // note : Selected_Column: 0 to 12;

        // Division : --------------------------------------------------------------------------------
        void TH2DSampleLineFill(TH2D * hist_in, double segmentation, std::pair<double,double> inner_clu, std::pair<double,double> outer_clu);
        std::vector<TH2D *> h2D_TrackLines_vec;

        // Division : --------------------------------------------------------------------------------
        void PrepareOutputFileName();
        void PrepareOutputFile();
        void PrepareHistFit();
        TFile * file_out;
        TTree * tree_out;

        int Column_ID_out; // note : 1 to 13;
        double l1_alignment_correction_column_out;
        double l1_alignment_correction_column_final_out;

        std::vector<TF1*> fit_l1_alignment_before_vec;
        std::vector<TF1*> fit_l1_alignment_after_vec;
        std::vector<TH1D*> h1D_l1_alignment_before_vec;
        std::vector<TH1D*> h1D_l1_alignment_after_vec;

        TH1D * h1D_Column_NTrack;
        TH1D * h1D_NTrack_NoZero;
        TH1D * h1D_NTrackWithZero;
        TH2D * h2D_TrackPos;

        TH1D * h1D_l0_ClusPosFitDiff;
        TH1D * h1D_l1_ClusPosFitDiff;
        TH1D * h1D_l2_ClusPosFitDiff;

        TH1D * h1D_l1_residual;

        // Division : --------------------------------------------------------------------------------
        // note : for preparing the final 
        TGraphErrors * gr_typeA_l1_alignment;
        TGraphErrors * gr_typeB_l1_alignment;
        TF1 * fit_typeA_l1_alignment;
        TF1 * fit_typeB_l1_alignment;

        TF1 * f1_pol1_unity;

        // Division : --------------------------------------------------------------------------------
        std::string color_code[8]   = {"#343434","#1A3947","#566575","#797983","#EFBD9D","#FCA26E","#F5751D","#F5321D"};
        std::string color_code_2[8] = {"#CC768D","#19768D","#DDA573","#009193","#6E9193","#941100","#A08144","#517E66"};

        double INTT_strip_width = 0.078;
        double lower_section_initial = -9.945;
        double upper_section_initial = 0.039;

        // note : the actual ladder position (unit : mm), the value is no longer to be 26.1, 
        // note : since 26.1 is the gap between ladders, without the consideration of ladder thickness
        double actual_xpos[3] = {0,29.552,59.104};

};

#endif