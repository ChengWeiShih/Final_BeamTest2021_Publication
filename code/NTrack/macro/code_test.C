#include "draw_style.h"

int code_test() {

    TH1D * aaa = new TH1D("aaa", "Test Histogram;Xaxis;Yaxis", 100, -10,10);
    aaa -> FillRandom("gaus",1000);

    TH1D * bbb = new TH1D("bbb", "Test Histogram 2;Xaxis;Yaxis", 100, -10,10);
    bbb -> FillRandom("gaus",5000);

    std::string output_directory = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/code/NTrack/macro";

    dataMC_comp(aaa, bbb, output_directory, Form("c1_h1D_l1_residual_test"), {"L1 - (L2L0 interpolation) [mm]","Entries (A.U.)"}, false, false);

    return 111;
}