#include "sPhenixStyle.C"

int test_1() {

    SetsPhenixStyle();
    TCanvas * c1 = new TCanvas("", "", 950, 800);
    c1 -> cd();

    // Division: ------------------------------------------------------------------------------------------------------
    TPad * pad1 = new TPad("", "", 0, 0.31, 0.85, 0.85);
    pad1 -> SetRightMargin(0.01);
    pad1 -> Draw();
    pad1 -> cd();

    TH1D * aaa = new TH1D("aaa", "aaa", 100, -5, 5);
    aaa -> GetXaxis() -> SetTitle("Y axis [mm]");
    aaa -> GetYaxis() -> SetTitle("Hit Detection Efficiency (%)");
    aaa -> GetXaxis() -> SetLabelOffset(999);
    aaa -> GetXaxis() -> SetLabelSize(0);
    aaa -> GetXaxis() -> SetTitleOffset(999);

    aaa -> Draw("hist");

    // Division: ------------------------------------------------------------------------------------------------------
    c1 -> cd();
    TPad * pad2 = new TPad("", "", 0, 0., 0.85, 0.38);
    pad2 -> SetRightMargin(0.01);
    pad2 -> Draw();
    pad2 -> cd();

    TH1D * ddd = new TH1D("ddd", "ddd", 100, -5, 5);
    ddd -> GetXaxis() -> SetTitle("Y axis [mm]");
    ddd -> GetYaxis() -> SetTitle("Hit Detection Efficiency (%)");
    ddd -> Draw("hist");

    // Division: ------------------------------------------------------------------------------------------------------
    c1 -> cd();
    TPad * pad3 = new TPad("", "", 0.85, 0.31, 1., 0.85);
    pad3 -> SetLeftMargin(0.05);
    pad3 -> Draw();
    pad3 -> cd();

    TH1D * his_p3 = new TH1D("his_p3", "his_p3", 1, -5, 5);   
    his_p3 -> GetXaxis() -> SetLabelOffset(999);
    his_p3 -> GetXaxis() -> SetNdivisions(200);
    his_p3 -> GetXaxis() -> SetTickLength(0.06);

    his_p3 -> GetYaxis() -> SetLabelSize(0);
    his_p3 -> GetYaxis() -> SetTickLength(0.13);

    his_p3 -> Draw("hist");

    // Division: ------------------------------------------------------------------------------------------------------
    c1 -> cd();
    TPad * pad4 = new TPad("", "", 0.85, 0., 1., 0.38);
    pad4 -> SetLeftMargin(0.05);
    pad4 -> Draw();
    pad4 -> cd();

    TH1D * ddd2 = new TH1D("ddd2", "ddd2", 1, -5, 5);
    ddd2 -> GetXaxis() -> SetTitle("Inclusive");
    ddd2 -> GetXaxis() -> SetLabelOffset(-0.04);
    ddd2 -> GetXaxis() -> SetLabelSize(0.11);
    ddd2 -> GetXaxis() -> SetNdivisions(2);
    ddd2 -> GetXaxis() -> SetTickLength(0.04);

    ddd2 -> GetYaxis() -> SetLabelSize(0);
    ddd2 -> GetYaxis() -> SetTickLength(0.13);

    ddd2 -> Draw("hist");


    // Division: ------------------------------------------------------------------------------------------------------
    c1 -> Print("test.pdf");

    return 111;
}