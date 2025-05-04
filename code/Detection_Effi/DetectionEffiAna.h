#ifndef DETECTIONEFFIANA_H
#define DETECTIONEFFIANA_H

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
#include <TEfficiency.h>

#include "sPhenixStyle.h"

class DetectionEffiAna{
    public:
        DetectionEffiAna(
            bool isData_in,
            int runnumber_in,
            std::string input_directory_in,
            std::string input_filename_in,
            std::string output_directory_in,
            std::string output_file_name_suffix_in,

            int Selected_Column_in, // note : 1 to 13;
            double l1_alignment_correction_in = 0.,
            double l0l2_slope_correction_in = 0.
        );

        std::string GetOutputFileName(){
            return output_filename;
        };


        double Get_l1_alignment(); // note : you may just do it twice in one run, before and after, if the input correction is nan
        double Get_l0l2_slope(); 

        void GetDetectionEffi(
            double slope_cut_in = 0.01, 
            double noise_hit_distance_in = 0.234, // note : unit : mm
            int boundary_cut_in = 5, 
            double cluster_adc_value_requirement_in = 22.5
        );

        void GetControlClusDist(
            double slope_cut_in = 0.01, 
            double noise_hit_distance_in = 0.234 // note : unit : mm
        );

        std::vector<double> Get_Final_Effi() {
            return {Final_Effi, Final_Effi_StatErrorUp, Final_Effi_StatErrorDown};
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
        double l0l2_slope_correction;

        std::string output_filename;

        double slope_cut;
        double noise_hit_distance;
        int boundary_cut;
        double cluster_adc_value_requirement;

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
        TH1D * h1D_ClusPos_Raw[3];

        TH1D * h1D_ClusPos_GoodEvt[3];
        
        TH1D * h1D_ClusPos_3Hit[3];
        
        TH1D * h1D_ClusAdc[3];
    
        TH1D * h1D_ClusSize[3];
        

        // Division : --------------------------------------------------------------------------------
        std::vector<std::vector<std::vector<std::vector<ClusInfo>>>> cluster_info_NoL1Aligned_vec; // note : [event][chip][layer][clusters]

        TF1 * fit_l1_alignment_before;
        TF1 * fit_l1_alignment_after;
        TF1 * fit_l0l2_slope_before;
        TF1 * fit_l0l2_slope_after;

        TH1D * Sub_Get_l1_alignment_hist();
        TH1D * Sub_Get_l0l2_slope_hist();

        // Division : --------------------------------------------------------------------------------
        void PrepareOutputFileName();
        void PrepareOutputFile();
        void PrepareHistFit();
        TFile * file_out;
        TTree * tree_out;
        TTree * tree_out_effi;
        
        TH1D * h1D_l1_alignment_before = nullptr;
        TH1D * h1D_l1_alignment_after = nullptr; 
        TH1D * h1D_l0l2_slope_before = nullptr;
        TH1D * h1D_l0l2_slope_after = nullptr;


        TH1D * h1D_effi_l1_residual = nullptr; // note : original, use l2 and l0 to check l1
        TH1D * h1D_effi_l0l2_slope = nullptr; // note : original, use l2 and l0 to check l1

        double L0L2Interpolation;
        bool L1Good;

        double Final_Effi;
        double Final_Effi_StatErrorUp;
        double Final_Effi_StatErrorDown;

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