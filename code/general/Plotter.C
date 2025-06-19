#include "sPhenixStyle.C"

void h1D_Plotter(bool isLogY = false){

    // std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/ResidualComp/Output";
    // std::string input_filename = "Data_Residual_Run52_Column8_test1.root";

    // std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/ResidualComp/MC";
    // std::string input_filename = "MC_Residual_Run52_Column11_test1.root";

    // std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run89/DetectionEffi_baseline700um";
    // std::string input_filename = "Data_DetectionEffi_Run89_Column10_BaseLine.root";

    std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/DetectionEffi_baseline700um";
    std::string input_filename = "Data_DetectionEffi_Run52_Column8_BaseLine.root";

    std::string output_name_flag = (isLogY) ? "LogY_" : "";

    // std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run64/NTrack";
    // std::string input_filename = "Data_NTrack_Run64_test1.root";


    SetsPhenixStyle();

    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleAlign(22);
    gStyle->SetTitleX(0.5);
	gStyle->SetTitleY(0.95);
	gStyle->SetTitleFillColor(10);
	gStyle->SetPadTopMargin(0.09);


    TCanvas * c1 = new TCanvas("c1", "c1", 950, 800);
    c1 -> cd();
    c1 -> SetLogy(isLogY);

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    // ltx->SetTextAlign(31);

    TFile * file_in = TFile::Open(Form("%s/%s", input_directory.c_str(), input_filename.c_str()));

    for (TObject* keyAsObj : *file_in->GetListOfKeys())
    {
        auto key = dynamic_cast<TKey*>(keyAsObj);
        std::string hist_name  = key->GetName();
        std::string class_name = key->GetClassName();

        if (class_name == "TH1D")
        {
            c1 -> cd();
            
            TH1D * hist = (TH1D*) file_in -> Get( hist_name.c_str() );

            double Y_axis_max_factor = (isLogY) ? 1000. : 1.5;

            hist -> SetMaximum( hist -> GetBinContent(hist -> GetMaximumBin()) * Y_axis_max_factor );
            hist -> GetXaxis() -> SetNdivisions(505);
            if (TString(hist->GetYaxis()->GetTitle()).IsNull())
            {
                hist -> GetYaxis() -> SetTitle("Entries");
            }

            hist -> Draw("hist");

            ltx->DrawLatex(0.38, 0.8, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));

            c1 -> Print(Form("%s/%s%s.pdf", input_directory.c_str(), output_name_flag.c_str(), hist_name.c_str()));
            c1 -> Clear();
        }
    }

}


void c1_Plotter(){

    // std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/ResidualComp/Output";
    // std::string input_filename = "Data_Residual_Run52_Column8_test1.root";

    // std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/ResidualComp/MC";
    // std::string input_filename = "MC_Residual_Run52_Column11_test1.root";

    // std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run89/DetectionEffi";
    // std::string input_filename = "Data_DetectionEffi_Run89_Column10_BaseLine.root";

    std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run64/NTrack_HotChannelMask";
    std::string input_filename = "Data_NTrack_Run64_test1.root";


    SetsPhenixStyle();

    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleAlign(22);
    gStyle->SetTitleX(0.5);
	gStyle->SetTitleY(1.1);
    gStyle->SetTitleYOffset(1.55);
	gStyle->SetTitleFillColor(10);
	gStyle->SetPadTopMargin(0.09);
    gStyle->SetPadLeftMargin(0.18);


    // TCanvas * c1 = new TCanvas("c1", "c1", 950, 800);
    // c1 -> cd();

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    // ltx->SetTextAlign(31);

    TFile * file_in = TFile::Open(Form("%s/%s", input_directory.c_str(), input_filename.c_str()));

    for (TObject* keyAsObj : *file_in->GetListOfKeys())
    {
        auto key = dynamic_cast<TKey*>(keyAsObj);
        std::string hist_name  = key->GetName();
        std::string class_name = key->GetClassName();

        if (class_name == "TCanvas")
        {
            // c1 -> cd();
            
            TCanvas * c1 = (TCanvas*) file_in -> Get( hist_name.c_str() );

            // hist -> SetMaximum( hist -> GetBinContent(hist -> GetMaximumBin()) * 1.5 );
            // hist -> GetXaxis() -> SetNdivisions(505);
            // if (TString(hist->GetYaxis()->GetTitle()).IsNull())
            // {
            //     hist -> GetYaxis() -> SetTitle("Entries");
            // }

            // hist -> Draw("hist");

            c1 -> Draw();

            ltx->DrawLatex(0.42, 0.92, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));

            c1 -> Print(Form("%s/%s.pdf", input_directory.c_str(), hist_name.c_str()));
            c1 -> Clear();
        }
    }

}

void Plot_L1AlignmentBefore(){

    std::string input_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/ResidualComp/Output";
    std::string input_filename = "Data_Residual_Run52_Column8_test1.root";
    std::string input_object = "c1_h1D_l1_alignment_before";
    
    TFile * file_in = TFile::Open(Form("%s/%s", input_directory.c_str(), input_filename.c_str()));

    SetsPhenixStyle();

    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleAlign(22);
    gStyle->SetTitleX(0.5);
	gStyle->SetTitleY(0.95);
	gStyle->SetTitleFillColor(10);
	gStyle->SetPadTopMargin(0.09);


    TCanvas * c2 = new TCanvas("c2", "c2", 950, 800);
    TCanvas * c1 = (TCanvas*) file_in -> Get(input_object.c_str());

    TH1D * hist = (TH1D*) c1 -> GetListOfPrimitives()->FindObject("h1D_l1_alignment_before");
    TF1 * fit   = (TF1*) c1 -> GetListOfPrimitives()->FindObject("fit_l1_alignment_before");

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    // ltx->SetTextAlign(31);

    c2 -> cd();
    hist -> SetMaximum( hist -> GetBinContent(hist -> GetMaximumBin()) * 1.5 );
    hist -> SetTitle("");
    hist -> GetXaxis() -> SetTitle("L1 misalignment [mm]");
    hist -> GetYaxis() -> SetMaxDigits(3);
    hist -> Draw("hist");
    fit -> SetLineColor(kRed);
    fit -> SetLineWidth(2);
    fit -> Draw("same");

    ltx->DrawLatex(0.38, 0.85, Form("#it{#bf{sPHENIX INTT}} Beam Test 2021"));
    ltx->DrawLatex(0.32, 0.77, Form("Run52, Column 8") );
    ltx->DrawLatex(0.32, 0.72, Form("Fit Gaus. mean: %.3f mm", fit->GetParameter(1)) );

    string input_filename_noSuffix = input_filename.substr(0, input_filename.find_last_of('.'));

    c2 -> Print(Form("%s/%s_%s.pdf", input_directory.c_str(), input_filename_noSuffix.c_str(), input_object.c_str()));
}