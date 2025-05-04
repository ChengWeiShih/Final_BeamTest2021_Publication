#ifndef RESIDUALCOMPANA_H
#define RESIDUALCOMPANA_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TGraph.h>

#include "sPhenixStyle.h"

class ResidualCompAna{
    public:
        ResidualCompAna(
            bool isData_in,
            int runnumber_in,
            std::string input_directory_in,
            std::string input_filename_in,
            std::string output_directory_in,
            std::string output_file_name_suffix_in,

            int Selected_Column_in, // note : 1 to 13;
            double l1_alignment_correction_in = 0.,
            double l0l1_slope_correction_in = 0.,
            int ClusPhiSize_cut_in = -999,
            int ClusAdc_cut_in = -999
        );

        std::string GetOutputFileName(){
            return output_filename;
        };
        double Get_l1_alignment(); // note : you may just do it twice in one run, before and after, if the input correction is nan
        double Get_l0l1_slope(); // note : you may just do it twice in one run, before and after, if the input correction is nan
        void GetHistsForComp(
            double slope_cut_in,
            double pos_cut_in
        );

        TH1D * Get_h1D_l1_residual() {
            return h1D_l1_residual;
        }
        TH1D * Get_h1D_l2_residual() {
            return h1D_l2_residual;
        }
        TH1D * Get_h1D_scattering() {
            return h1D_scattering;
        }

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
        int Selected_Column; // note : 1 to 13;
        double l1_alignment_correction;
        double l0l1_slope_correction;
        int ClusPhiSize_cut;
        int ClusAdc_cut;

        std::string output_filename;
        double slope_cut;
        double pos_cut;


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

        // Division : --------------------------------------------------------------------------------
        std::vector<std::vector<std::vector<std::vector<ClusInfo>>>> cluster_info_NoL1Aligned_vec; // note : [event][chip][layer][clusters]
        // std::vector<std::vector<std::vector<ClusInfo>>> cluster_info_PostL1Aligned_vec; // note : [event][layer][clusters]

        TF1 * fit_l1_alignment_before;
        TF1 * fit_l1_alignment_after;
        TF1 * fit_l0l1_slope_before;
        TF1 * fit_l0l1_slope_after;

        TH1D * Sub_Get_l1_alignment_hist();
        TH1D * Sub_Get_l0l1_slope_hist();

        TGraph * ThreePointGrr;
        TF1 * fit_linear_fit;

        // Division : --------------------------------------------------------------------------------
        void PrepareOutputFileName();
        void PrepareOutputFile();
        void PrepareHistFit();
        TFile * file_out;
        TTree * tree_out;
        
        TH1D * h1D_l1_alignment_before = nullptr;
        TH1D * h1D_l1_alignment_after = nullptr; 
        TH1D * h1D_l0l1_slope_before = nullptr;
        TH1D * h1D_l0l1_slope_after = nullptr;

        TH1D * h1D_l1_residual = nullptr; // note : original, use l2 and l0 to check l1
        TH1D * h1D_l2_residual = nullptr; // note : new one, used l1 and l0 to check l2
        TH1D * h1D_scattering = nullptr;  // note : angle between (l2-l1) and (l1-l0)

        

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