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

void B_JHEP_10_072(Int_t Flag = 0){
  //----------------------------------------------------------------------------
  //  Flag == 0: CMS combination full breakdown
  //  Flag == 1: CMS combination with scan of correlation for total sytematics
  //----------------------------------------------------------------------------
  
  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  4;
  static const Int_t NumUnc = 11;
  static const Int_t NumObs =  2;
  
  // Array of names of estimates
  const TString NamEst[NumEst] = {
    "Wgg_{el}", "Wgg_{mu}", "Zgg_{el}", "Zgg_{mu}"
  };
  
  // Array of names of uncertainties
  TString NamUnc[NumUnc] = {
    //    0         1         2         3         4         5
    "  Stat", "SimSta", " Trigg", "Lepton", "PTmiss", "Pileup", 
    //    6         7         8         9        10
    "PDF+Sc", "MisIdJ", "MisIdE", "GPromp", "  Lumi"
  };
  
  // Array of names of observables
  const TString NamObs[NumObs] = {
    "    Wgg", "    Zgg"
  };
  
  // Index for which estimates determines which observable
  const Int_t IWhichObs[NumEst] = {0, 0, 1, 1};
  
  // The steering flag
  if(Flag == 0){    
    printf("... B_JHEP_10_072: ----------------------------------------- \n");
    printf("... B_JHEP_10_072: CMS combination full breakdown, Flag = %2i \n",
	   Flag);
  }else if(Flag == 1){
    printf("... B_JHEP_10_072: --------------------------------------------");
    printf("---------------------------- \n");
    printf("... B_JHEP_10_072: CMS combination with scan of correlation");
    printf(" for total sytematics, Flag = %2i \n",Flag);
    NamUnc[1] = " SysTot";
  }else{
    printf("... B_JHEP_10_072: ------------------------------------- \n");
    printf("... B_JHEP_10_072:            Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  // Define array for SysTot
  Double_t SysTot[NumEst] = {0.};
  
  // Fill estimates and percentage uncertainties
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    4.2,  47.8, 2.8, 0.5, 4.1, 1.5, 0.5, 1.5, 36.6, 6.9, 6.7, 2.6,
    6.0,  29.6, 2.4, 0.3, 3.0, 1.4, 0.2, 1.6, 37.2, 0.0, 5.8, 2.6,
    12.5, 16.6, 3.3, 1.3, 5.3, 0.0, 1.3, 1.2, 15.1, 0.0, 0.2, 2.6,
    12.8, 13.7, 2.9, 1.2, 4.3, 0.0, 0.4, 1.3, 12.5, 0.0, 0.3, 2.6
  };
  
  // Multiply percentage uncertainties with measured value, and calculate SysTot
  for(Int_t i = 0; i<NumEst; i++){
    for(Int_t k = 0; k<NumUnc; k++){
      XEst[i*(NumUnc+1)+k+1] = XEst[i*(NumUnc+1)] * XEst[i*(NumUnc+1)+k+1]/100.;
      if(k>0 && k<NumUnc){
	SysTot[i] = SysTot[i] +  XEst[i*(NumUnc+1)+k+1]* XEst[i*(NumUnc+1)+k+1];
      }
    }
    SysTot[i] = TMath::Sqrt(SysTot[i]);
  }
  
  // Define assumed correlations:
  Double_t RhoVal[NumUnc] = {
    //0    1    2    3    4    5    6    7    8    9   10
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0
  };
  const Double_t RhoFut = 0.4;
  
  // Define formats for figures and latex file
  const TString FilBas = "B_JHEP_10_072";
  TString FilNam;
  char Buffer[100];
  const TString ForVal = "%5.1f";
  const TString ForUnc = "%4.1f";
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

  // Diplay the assumed correlations
  printf("... B_JHEP_10_072: --------------------------------------------");
  printf("--------------------------------------- \n");
  printf("... B_JHEP_10_072: Since the paper  is not explicit about the");
  printf(" estimator correlations for the various\n");
  printf("... B_JHEP_10_072: sources of uncertainty listed in Table 3,");
  printf(" the following assumptions have been made. \n");
  printf("... B_JHEP_10_072: for Flag==0.");
  printf("... B_JHEP_10_072: \n");
  printf("... B_JHEP_10_072:     Simulation stat uncertainty: Rho = %3.1f \n",
	 RhoVal[1]);
  printf("... B_JHEP_10_072:                         Trigger: Rho = %3.1f \n",
	 RhoVal[2]);
  printf("... B_JHEP_10_072: Lepton and photon ID and energy: Rho = %3.1f \n",
	 RhoVal[3]);
  printf("... B_JHEP_10_072:                  p_T_miss scale: Rho = %3.1f \n",
	 RhoVal[4]);
  printf("... B_JHEP_10_072:                          Pileup: Rho = %3.1f \n",
	 RhoVal[5]);
  printf("... B_JHEP_10_072:  PDFs, renorm. and fact. scales: Rho = %3.1f \n",
	 RhoVal[6]);
  printf("... B_JHEP_10_072:               Misidentified jet: Rho = %3.1f \n",
	 RhoVal[7]);
  printf("... B_JHEP_10_072:          Misidentified electron: Rho = %3.1f \n",
	 RhoVal[8]);
  printf("... B_JHEP_10_072:                 Prompt diphoton: Rho = %3.1f \n",
	 RhoVal[9]);
  printf("... B_JHEP_10_072:               Total statistical: Rho = %3.1f \n",
	 RhoVal[0]);
  printf("... B_JHEP_10_072:           Integrated luminosity: Rho = %3.1f \n",
	 RhoVal[10]);
  printf("... B_JHEP_10_072: \n");
  printf("... B_JHEP_10_072: Flag==0:\n");
  printf("... B_JHEP_10_072: With the above the Zgg combination comes out on");
  printf(" spot, while for the Wgg combination\n");
  printf("... B_JHEP_10_072: a slightly larger systematic");
  printf(" uncertainty is seen.\n");
  printf("... B_JHEP_10_072: Flag==1:\n");
  printf("... B_JHEP_10_072: Using an estimator correlation for the total");
  printf(" systematic uncertainty of %3.1f, the\n", RhoFut);
  printf("... B_JHEP_10_072: paper result for Wgg");
  printf(" is recovered (while now Zgg is slightly off).\n");
  printf("... B_JHEP_10_072: --------------------------------------------");
  printf("--------------------------------------- \n");

  // Prepare for only using the total uncertainty
  // 1) Keep stat uncertainty k==0
  // 2) Put SysTot in k==1
  // 3) Set correlation for k==1 to RhoFut
  // 4) Before combining disable uncertainties k>1
  if(Flag == 1){
    for(Int_t i = 0; i<NumEst; i++)XEst[i*(NumUnc+1)+1+1] = SysTot[i];
    RhoVal[1] = RhoFut;
  }
  
  // Fill all estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }
  
  // Fill correlation
  for(Int_t k = 0; k<NumUnc; k++)myBlue->FillCor(k, RhoVal[k]);
  
  // Prepare for combination with stat + tot_sys + lumi
  if(Flag == 1){
    for(Int_t k = 2; k<NumUnc-1; k++)myBlue->SetInActiveUnc(k);
  }
  
  // Combination for the W
  myBlue->SetInActiveEst(2);
  myBlue->SetInActiveEst(3);
  sprintf(Buffer,"_Wgg_%1i", Flag);
  FilNam = FilBas + &Buffer[0];
  myBlue->FixInp();
  myBlue->PrintEst();
  myBlue->Solve();
  if(Flag == 1){
    printf("... B_JHEP_10_072:\n");
    printf("... B_JHEP_10_072: Using an estimator correlation of %3.1f, the\n",
	   RhoFut);
    printf("... B_JHEP_10_072: CMS Wgg result is retained, while now the \n");
    printf("... B_JHEP_10_072: Zgg combination is slighly off.\n");
    printf("... B_JHEP_10_072:\n");
  }  
  myBlue->PrintResult();
  printf("... B_JHEP_10_072:\n");
  printf("... B_JHEP_10_072: The CMS result for Wgg is: \n");
  printf("... B_JHEP_10_072: Wgg Value +-  Full (stat +- syst +- lumi) fb\n");
  printf("... B_JHEP_10_072: Wgg   4.9 +- %5.1f ( 1.4 +-  1.6 +-  0.1) fb\n",
	 TMath::Sqrt(1.4*1.4+1.6*1.6+0.1*0.1));
  printf("... B_JHEP_10_072:\n");
  myBlue->DisplayPair(0, 1, FilNam, 2, 10, 0, 5);
  myBlue->LatexResult(FilNam);
  myBlue->ReleaseInp();

  // Combination for the Z
  myBlue->SetActiveEst(2);
  myBlue->SetActiveEst(3);
  myBlue->SetInActiveEst(0);
  myBlue->SetInActiveEst(1);
  sprintf(Buffer,"_Zgg_%1i", Flag);
  FilNam = FilBas + &Buffer[0];
  myBlue->FixInp();
  myBlue->PrintEst();
  myBlue->Solve();
  myBlue->PrintResult();
  printf("... B_JHEP_10_072:\n");
  printf("... B_JHEP_10_072: The CMS result for Zgg is: \n");
  printf("... B_JHEP_10_072: Zgg Value +-  Full (stat +- syst +- lumi) fb\n");
  printf("... B_JHEP_10_072: Zgg  12.7 +- %5.1f ( 1.4 +-  1.8 +-  0.3) fb\n",
	 TMath::Sqrt(1.4*1.4+1.8*1.8+0.3*0.3));
  myBlue->DisplayPair(2, 3, FilNam, 10, 19, 0, 5);
  myBlue->LatexResult(FilNam);

  // Delete Object and TMatrices
  delete myBlue;
  myBlue = NULL;

  // Return
  return;
}
