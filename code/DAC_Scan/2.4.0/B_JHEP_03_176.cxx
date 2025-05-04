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

void B_JHEP_03_176(Int_t Flag = 0){
  //----------------------------------------------------------------------------
  //  Flag == 0: LHCb combination
  //----------------------------------------------------------------------------
  
  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  8;
  static const Int_t NumUnc =  9;
  static const Int_t NumObs =  4;
  
  // Array of names of estimates
  const TString NamEst[NumEst] = {
    "Rat0_{do}", "Rat0_{up}", "Rat1_{do}", "Rat1_{up}",
    "Rat2_{do}", "Rat2_{up}", "Rat3_{do}", "Rat3_{up}",
  };
  
  // Array of names of uncertainties
  TString NamUnc[NumUnc] = {
    //    0         1         2         3         4         5
    "  Stat", "SizSim", "   PID", "Tracks", " Trigg", "MatDes", 
    //    6         7         8         9        10
    "FitMod", "SecDec", " DPBin"
  };
  
  // Array of names of observables
  const TString NamObs[NumObs] = {
    "Ratio0", "Ratio1", "Ratio2", "Ratio3"
  };
  
  // Index for which estimates determines which observable
  const Int_t IWhichObs[NumEst] = {0, 0, 1, 1, 2, 2, 3, 3};
  
  // The steering flag
  if(Flag == 0){
    printf("... B_JHEP_03_176: ---------------------------- \n");
    printf("... B_JHEP_03_176: LHCb combination, Flag = %2i \n",
	   Flag);
  }else{
    printf("... B_JHEP_03_176: ---------------------------- \n");
    printf("... B_JHEP_03_176:   Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  // Fill estimates and percentage uncertainties
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    // Ratio 1
    0.653, 0.54, 0.34, 0.022, 0.22, 0.011, 0.53, 0.14, 0.18, 0.09,
    0.655, 0.54, 0.32, 0.030, 0.22, 0.011, 0.53, 0.13, 0.18, 0.11,
    // Ratio 2
    5.220, 0.25, 0.47, 0.019, 0.069, 0.0025, 0.0, 0.03, 0.25, 0.05,
    5.244, 0.25, 0.52, 0.020, 0.070, 0.0024, 0.0, 0.07, 0.24, 0.03,
    //Ratio 3
    2.333, 1.4, 1.0, 0.022, 0.079, 0.0050, 0.0, 0.64, 0.31, 0.30,
    2.419, 1.4, 1.2, 0.023, 0.080, 0.0057, 0.0, 0.54, 0.38, 0.07,
    //Ratio 4
    103.00, 0.03, 0.75, 0.013, 0.11, 0.0057, 0.27, 0.06, 0.11, 0.13,
    102.59, 0.03, 0.81, 0.021, 0.10, 0.0060, 0.27, 0.06, 0.09, 0.28
  };
  
  // Multiply percentage uncertainties with measured value
  for(Int_t i = 0; i<NumEst; i++){
    for(Int_t k = 0; k<NumUnc; k++){
      XEst[i*(NumUnc+1)+k+1] = XEst[i*(NumUnc+1)] * XEst[i*(NumUnc+1)+k+1]/100.;
    }
  }
  
  // Define estimator correlations
  Double_t RhoVal[NumUnc] = {
    //0    1    2    3    4    5    6    7    8
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0
  };
  
  //-- Local Structures for Blue output
  // TMatrices
  TMatrixD* LocRes = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc = new TMatrixD(NumObs,1);
  //-- End
  
  // Arrays to hold results for final print out
  Double_t RatVal[NumObs] = {0.};
  Double_t RatSta[NumObs] = {0.};
  Double_t RatSys[NumObs] = {0.};

  // Define formats for figures and latex file
  const TString FilBas = "B_JHEP_03_176";
  TString FilNam;
  char Buffer[100];
  TString ForVal = "%7.3f";
  TString ForUnc = "%6.3f";
  const TString ForWei = "%4.2f";
  const TString ForRho = ForWei;
  const TString ForPul = ForRho;
  const TString ForUni = "GeV";
  
  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);
  myBlue->PrintStatus();
  
  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);

  // Fill all estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }
  
  // Fill correlation
  for(Int_t k = 0; k<NumUnc; k++)myBlue->FillCor(k, RhoVal[k]);
  
  // Deactivate all estimates
  for(Int_t i = 0; i<NumEst; i++)myBlue->SetInActiveEst(i);

  // Loop over combinations for the four ratios
  for(Int_t m = 0; m<NumObs; m++){
    
    // Set formats per observable
    if(m == 0){
      ForVal = "%8.4f";
      ForUnc = "%7.4f";
    }else if(m == 3){
      ForVal = "%7.2f";
      ForUnc = "%6.2f";
    }else{
      ForVal = "%7.3f";
      ForUnc = "%6.3f";
    }
    myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);

    // Select the estimates for observable m
    myBlue->SetActiveEst(2*m);
    myBlue->SetActiveEst(2*m+1);

    // Fix, and solve
    myBlue->FixInp();
    myBlue->Solve();

    // Get result into local structures
    myBlue->GetResult(LocRes);
    RatVal[m] = LocRes->operator()(0,0);
    RatSta[m] = LocRes->operator()(0,1);
    myBlue->GetUncert(LocUnc);
    RatSys[m] = LocUnc->operator()(0,0);
    RatSys[m] = TMath::Sqrt(RatSys[m]*RatSys[m] - RatSta[m]*RatSta[m]);

    // Set filename
    sprintf(Buffer,"_Rat_%1i", m);
    FilNam = FilBas + &Buffer[0];
    myBlue->DisplayPair(m, m+1, FilNam);
    myBlue->LatexResult(FilNam);
    myBlue->ReleaseInp();

    // De-Select the estimates for observable m
    myBlue->SetInActiveEst(2*m);
    myBlue->SetInActiveEst(2*m+1);
  }

  // Print out what we got
  printf("... B_JHEP_03_176: \n");
  printf("... B_JHEP_03_176: The combined results using the units they");
  printf("are presented in the paper are: \n");
  for(Int_t m = 0; m<NumObs; m++){
    printf("... B_JHEP_03_176:");
    if(m == 0){ 
      printf(" Ratio 0 = (%5.3f +- %5.3f +- %5.3f) x 10^{-4} \n",
	     10*RatVal[m], 10.*RatSta[m], 10.*RatSys[m]);
    }else if(m == 1){ 
      printf(" Ratio 1 = (%5.3f +- %5.3f +- %5.3f) x 10^{-3} \n",
	     RatVal[m], RatSta[m], RatSys[m]);
    }else if(m == 2){ 
      printf(" Ratio 2 = (%5.3f +- %5.3f +- %5.3f) x 10^{-3} \n",
	     RatVal[m], RatSta[m], RatSys[m]);
    }else if(m == 3){ 
      printf(" Ratio 1 = (%5.3f +- %5.3f +- %5.3f) x 10^{-2} \n",
	     RatVal[m]/10., RatSta[m]/10., RatSys[m]/10.);
    }
  }
  printf("... B_JHEP_03_176: \n");
  
  // Delete Object and TMatrices
  delete myBlue;
  myBlue = NULL;
  LocRes->Delete(); LocRes = NULL;
  LocUnc->Delete(); LocUnc = NULL;

  // Return
  return;
}
