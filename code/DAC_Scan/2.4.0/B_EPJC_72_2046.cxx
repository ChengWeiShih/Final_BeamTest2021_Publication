//------------------------------------------------------------------------------
//
// BLUE: A ROOT class implementing the Best Linear Unbiased Estimate method.
// 
// Copyright (C) 2012-2019, Richard.Nisius@mpp.mpg.de
// All rights reserved
//
// This file is part of BLUE - Version 2.2.0.
//
// BLUE is free software: you can redistribute it and/or modify it under the 
// terms of the GNU Lesser General Public License as published by the Free 
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// For the licensing terms see the file COPYING or http://www.gnu.org/licenses.
//
//------------------------------------------------------------------------------
#include "Blue.h"

void B_EPJC_72_2046(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0:  mtop (1d + 2d full combination)
  //  1:  mtop (1d independent combination)
  //     +mtop (2d independent combination)
  //  2:  mtop(1d), mtop(2d) (treat 1d + 2d as two correlated observables)
  //  3:  mtop(el), mtop(mu) (treat el + mu as two correlated observables)
  //      [Shows how to use InspectPair()]
  //  4:   same as 3, but mtop(el) and mtop(mu) are independent
  //      [Shows how to de-activate uncertainties with SetInActiveUnc()]
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  4;
  static const Int_t NumUnc = 20;
  static const Int_t MaxObs = 2;
  Int_t NumObs =  2;

  TString NamEst[NumEst] = {"  1d-el", "  1d-mu", "  2d-el", "  2d-mu"};

  //                               0          1          2          3
  TString NamUnc[NumUnc] = {"   Stat", "   Meth", "  MCgen", "    Had",
  //                               4          5          6          7
			    "  Pile", "     UE",  "     CR", "   IFSR",
  //                               8          9         10         11
			    "    PDF", "  Wnorm", "  Wsha", "  Qnorm", 
  //                              12         13         14         15
			    "   Qsha", "   iJES", "    JES", "   bJES",
  //                              16         17         18         19
			    "  btag", "    JER",  "   JEFF", "    MET"};

  TString NamObs[MaxObs] = {"   Mtop", "   Mtop"};

  // Index for which estimates determines which observable
  Int_t IWhichObs[] = {0, 0, 0, 0};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_EPJC_72_2046: I use Flag = %2i \n",Flag);
    NumObs = 1;
  }else if(Flag == 1){
    printf("... B_EPJC_72_2046: I use Flag = %2i \n",Flag);
    NumObs = 1;
    NamObs[0] = "Mtop-xd";
  }else if(Flag == 2){
    IWhichObs[2] = 1;
    IWhichObs[3] = 1;
    NamObs[0] = "Mtop-1d";
    NamObs[1] = "Mtop-2d";
  }else if(Flag == 3 || Flag == 4){
    IWhichObs[1] = 1;
    IWhichObs[3] = 1;
    if(Flag == 3){
      NamObs[0] = "Mtop-el";
      NamObs[1] = "Mtop-mu";
    }else{
      NamObs[0] = "Mtop-xx";
    }
  }else{
    printf("... B_EPJC72_2046: Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  // Estimates 0-3 == 1d-el, 1d-mu, 2d-el, 2d-mu
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //         0     1     2     3     4     5     6     7     8     9    10
    172.93, 1.46, 0.00, 0.07, 0.81, 0.33, 0.05,  0.06, 0.47, 1.45, 0.22, 0.16,
    //11      12    13    14    15    16    17     18    19
    0.11,   0.07, 0.14, 1.21, 1.09, 0.21, 0.34,  0.08, 0.05,
    //         0     1     2     3     4     5      6     7     8     9    10
    175.54, 1.13, 0.00, 0.04, 0.69, 0.52, 0.05,  0.10, 0.74, 1.40, 0.09, 0.19,
    //11      12    13    14    15    16    17     18    19
    0.18,   0.05, 0.12, 1.25, 1.21, 0.13, 0.38,  0.11, 0.05,
    //         0     1     2     3     4     5      6     7     8     9    10
    174.30, 0.83, 0.59, 0.10, 0.39, 0.20, 0.01,  0.42, 0.32, 1.04, 0.10, 0.34,
    //11      12    13    14    15    16    17     18    19
    0.07,   0.25, 0.38, 0.63, 1.61, 0.31, 0.07, 0.005, 0.12,
    //         0     1     2     3     4     5      6     7     8     9    10
    175.01, 0.74, 0.51, 0.03, 0.22, 0.06, 0.01,  0.96, 1.04, 0.95, 0.10, 0.44,
    //11      12    13    14    15    16    17     18    19
    0.22,   0.33, 0.30, 0.71, 1.53, 0.26, 0.07,  0.01, 0.16
  };

  static const Int_t LenCor = NumEst * NumEst;
  // Statistical Correlation = Uncertainty 0
  Double_t CorStat[LenCor] = {
    1.00, 0.00, 0.15, 0.00,
    0.00, 1.00, 0.00, 0.16,
    0.15, 0.00, 1.00, 0.00,
    0.00, 0.16, 0.00, 1.00
  };
  
  // Qnorm, Qshape = Uncertainty 12, 13
  // Dense format, use only the off diagonal elments of one half of the matrix
  static const Int_t LenCorSig = (NumEst * NumEst - NumEst)/2;
  Double_t CorQCD[LenCorSig] = {
    0.00, 1.00, 0.00, 0.00, 1.00, 0.00,
  };

  // Uncorrelated: Meth/iJES = Uncertainty 1/2
  // The rest uses full correlation: 
  // MCGen/Had/Pile/UE/CR/IFSR/PDF/Wnorm/Wsha = Unc 3-11
  // JES/bJES/btag/JER/JEFF/MET = Unc 14-19

  // Define formats for figures and latex file
  const TString FilBas = "B_EPJC_72_2046";
  TString FilNam;
  char Buffer[100];
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%4.2f";
  const TString ForRho = ForUnc;
  const TString ForPul = ForRho;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";

  // Construct Object 
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);

  // Fill estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }

  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(k ==0){
      myBlue->FillCor(k,&CorStat[0]);
    }else if(k == 1 || k == 2){
      myBlue->FillCor(k,0.0);
    }else if(k == 12 || k == 13){
      myBlue->FillCor(-k,&CorQCD[0]);
    }else{
      myBlue->FillCor(k,1.0);
    }
  }

  // Fix input, solve according to Flag and finally delete
  if(Flag == 0){
    // Fix and solve
    myBlue->FixInp(); 
    myBlue->PrintEst();
    myBlue->Solve();
    printf("... B_EPJC_72_2046: mtop(1d + 2d full combination) \n");
    myBlue->PrintResult();
    myBlue->PrintCompatEst("B_EPJC_72_2046");
    myBlue->PrintPull();
    myBlue->LatexResult("B_EPJC_72_2046");
    myBlue->DisplayResult(0,"B_EPJC_72_2046");
  }else if(Flag == 1){
    // Set more detailed print level 
    myBlue->SetPrintLevel(1);

    // Change, fix and solve
    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(3);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_EPJC_72_2046: mtop(1d independent combination) \n");
    myBlue->PrintResult();
    myBlue->PrintPull();

    // Reset, change, fix and solve
    myBlue->ResetInp();
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(1);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_EPJC_72_2046: mtop(2d independent combination) \n");
    myBlue->PrintResult();
    myBlue->PrintPull();
  }else if(Flag == 2){
    myBlue->FixInp();
    myBlue->PrintCompatEst();
    myBlue->Solve();
    myBlue->PrintWeight();
    myBlue->PrintRhoRes();
    printf("... B_EPJC_72_2046:  mtop(1d), mtop(2d) (treat 1d + 2d");
    printf(" as two correlated Observables) \n");
    myBlue->PrintResult();
    myBlue->PrintPull();
    myBlue->PrintCompatObs();
  }else if(Flag == 3){
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_EPJC_72_2046: mtop(el), mtop(mu) (treat xd-el + xd-mu");
    printf(" as two correlated Observables) \n");
    myBlue->PrintResult();
    myBlue->PrintPull();

    // Inspect a pair
    myBlue->PrintEst(1);
    myBlue->PrintEst(3);   
    myBlue->InspectPair(1,3,"B_EPJC_72_2046",1);
  }else if(Flag == 4){
    myBlue->SetInActiveEst(1);
    myBlue->SetInActiveEst(3);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_EPJC_72_2046: mtop(el) from 1d(el) and 2d(el) analyses \n");
    myBlue->PrintResult();
    printf("... B_EPJC_72_2046: Pull(%2i) = %5.3f \n", 0, myBlue->GetPull(0));
    printf("... B_EPJC_72_2046: Pull(%2i) = %5.3f \n", 2, myBlue->GetPull(2));
    myBlue->PrintPull();

    myBlue->ResetInp();
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(2);
    myBlue->FixInp();
    myBlue->PrintListEst();
    myBlue->PrintListUnc();
    myBlue->PrintListObs();
    myBlue->Solve();
    printf("... B_EPJC_72_2046: mtop(mu) from 1d(mu) and 2d(mu) analyses \n");
    myBlue->PrintResult();
    myBlue->PrintPull();

    // De-activate some small uncertainties
    myBlue->ReleaseInp();
    myBlue->SetInActiveUnc(2);
    myBlue->SetInActiveUnc(5);
    myBlue->SetInActiveUnc(9);
    myBlue->SetInActiveUnc(18);
    myBlue->SetInActiveUnc(19);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_EPJC_72_2046: mtop(mu) from 1d(mu) and 2d(mu) analyses \n");
    printf("... B_EPJC_72_2046: after removing some unimportant");
    printf(" sources of uncertainty. The rest of the uncertainties \n");
    printf("... B_EPJC_72_2046: stays almost unchanged \n");
    myBlue->PrintResult();

    // De-activate a big one
    myBlue->ReleaseInp();
    myBlue->SetInActiveUnc(8);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_EPJC_72_2046: mtop(mu) from 1d(mu) and 2d(mu) analyses \n");
    printf("... B_EPJC_72_2046: after in addition removing");
    printf(" ISR/FSR uncertainty. \n");
    printf("... B_EPJC_72_2046: Now the 1d(mu) gets");
    printf(" more precise than 2d(mu) \n");
    myBlue->PrintResult();

    // Artificially rescale the correlation for the JES(14)
    myBlue->ReleaseInp();
    myBlue->SetRhoFacUnc(14, 0.5);
    myBlue->FixInp();
    myBlue->PrintCor(14);
    myBlue->Solve();
    printf("... B_EPJC_72_2046: mtop(mu) from 1d(mu) and 2d(mu) analyses \n");
    printf("... B_EPJC_72_2046: after also artificially setting");
    printf(" rho(JES)=0.5.\n");
    myBlue->PrintResult();
  }
  // Delete Object
  delete myBlue;
  myBlue = NULL;
  return;
}
