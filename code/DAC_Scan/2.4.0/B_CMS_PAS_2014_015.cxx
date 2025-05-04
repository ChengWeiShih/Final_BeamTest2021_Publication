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

//---------------------------------------------------------------------------
// Function prototype for utility
//---------------------------------------------------------------------------
void FilMat(const Double_t RD, const Double_t RS, TMatrixD *const M);

void B_CMS_PAS_2014_015(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: CMS combination with --reduced-- correlations
  //  1: CMS combination with --nominal-- correlations
  //  2: CMS combination with changed uncorrelated JES (uncJES=4)
  //  3: CMS combination with changed flavour dependent had. unc. (flaHad=16)
  //  4: CMS combination with changed jet energy resolution (JER=8)
  //  5: CMS combination per decay channel
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  8;
  static const Int_t NumUnc = 25;
  Int_t NumObs =  1;
  static const Int_t MaxObs =  3;

  const TString NamEst[NumEst] = {" 10-dil", " 10-l+j", 
				  " 11-dil", " 11-l+j", " 11-had",
				  " 12-dil", " 12-l+j", " 12-had"};
  const TString NamUnc[NumUnc] = {
    //     0          1          2          3          4          5          6 
    "   Stat", "   iJES", "JEC-Int", "JEC-Ins", " uncJES", " remJES", "   Lept",
    //     7          8          9         10         11         12         13 
    " ETmiss", "    JER", "  b-tag", "   Trig", " Pileup", "   BGDT", "   BGMC",
    //    14         15         16         17         18         19         20 
    " FitCal", " flaJES", " flaHad", " b-frag", "    PDF", "  muR/F", "   MEPS",
    //    21         22         23         24  
    "     MC", " top-pt", "     UE", "    CR",};

  TString NamObs[MaxObs] = {"mas-l+j", "mas-dil", "mas-had"};

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {NumEst*0};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_CMS_PAS_2014_015: ----------------------------------- \n");
    printf("... B_CMS_PAS_2014_015: The 2014 CMS combination, Flag = %2i \n",
	   Flag);
    NamObs[0] = "   mtop";
  }else if(Flag == 1){
    printf("... B_CMS_PAS_2014_015: ------------------------------");
    printf("------------------------------- \n");
    printf("... B_CMS_PAS_2014_015: The 2014 CMS combination, with");
    printf(" normal correlations, Flag = %2i \n", Flag);
    NamObs[0] = "   mtop";
  }else if(Flag == 2){
    printf("... B_CMS_PAS_2014_015: --------------------------------\n");
    printf("... B_CMS_PAS_2014_015: Test uncorrelated JES, Flag = %2i\n", Flag);
    NamObs[0] = "   mtop";
  }else if(Flag == 3){
    printf("... B_CMS_PAS_2014_015: ------------------------------------- \n");
    printf("... B_CMS_PAS_2014_015: Test jet energy resolution,");
    printf(" Flag = %2i \n", Flag);
    NamObs[0] = "   mtop";
  }else if(Flag == 4){
    printf("... B_CMS_PAS_2014_015: ------------------------------------- \n");
    printf("... B_CMS_PAS_2014_015: Test flavour JES component,");
    printf(" Flag = %2i \n", Flag);
    NamObs[0] = "   mtop";
  }else if(Flag == 5){
    printf("... B_CMS_PAS_2014_015: ---------------------------------");
    printf("---------------------------- \n");
    printf("... B_CMS_PAS_2014_015: The combination for l+jets,");
    printf(" di-lepton and all-jets, Flag = %2i \n",Flag);
    NumObs = 3;
    IWhichObs[0] = 1;
    IWhichObs[1] = 0;
    IWhichObs[2] = 1;
    IWhichObs[3] = 0;
    IWhichObs[4] = 2;
    IWhichObs[5] = 1;
    IWhichObs[6] = 0;
    IWhichObs[7] = 2;
    NamObs[0] = " mt-l+j";
    NamObs[1] = " mt-dil";
    NamObs[2] = " mt-had";
  }else{
    printf("... B_CMS_PAS_2014_015: ------------------------ \n");
    printf("... CMS_PAS_2014_015: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Estimates 0-7
  // 0 ==   CMS_2010   di-lepton
  // 1 ==   CMS_2010 lepton+jets
  // 2 ==   CMS_2011   di-lepton
  // 3 ==   CMS_2011 lepton+jets
  // 4 ==   CMS_2011    all-jets
  // 5 ==   CMS_2012   di-lepton
  // 6 ==   CMS_2012 lepton+jets
  // 7 ==   CMS_2012    all-jets
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  const Double_t XEst[LenXEst] = {
    //         0      1      2      3      4      5      6      7      8
    //      Stat   iJES intJES MPFJES uncJES remJES   Lept ETmiss    JER
    //         9     10     11     12     13     14     15     16     17
    //     b-tag   Trig Pileup   BGDT   BGMC FitCal flaJES flaHad b-frag 
    //        18     19     20     21     22     23     24
    //       PDF mu_R/F   MEPS     MC top-pt     UE     CR
    175.50, 4.60,  0.00,  0.17,  0.76,  1.48,  3.28,  0.30,  0.10,  0.50,
    0.40,   0.00,  1.00,  0.00,  0.10,  0.20,  1.21,  0.90,  0.00,  0.50,
    0.63,   0.70,  0.50,  0.00,  1.30,  0.00,
    173.10, 2.10,  0.00,  0.08,  0.16,  1.90,  0.00,  0.00,  0.40,  0.10,
    0.10,   0.00,  0.10,  0.40,  0.20,  0.10,  0.87,  0.90,  0.00,  0.10,
    1.12,   0.40,  0.00,  0.00,  0.20,  0.00,
    172.50, 0.43,  0.00,  0.08,  0.35,  0.69,  0.00,  0.14,  0.12,  0.14,
    0.09,   0.00,  0.11,  0.00,  0.05,  0.40,  0.58,  0.76,  0.00,  0.09,
    0.55,   0.19,  0.04,  0.00,  0.05,  0.13,

    173.49, 0.27,  0.33,  0.01,  0.02,  0.24,  0.00,  0.02,  0.06,  0.23,
    0.12,   0.00,  0.07,  0.00,  0.13,  0.06,  0.11,  0.61,  0.00,  0.07,
    0.24,   0.18,  0.02,  0.00,  0.15,  0.54,
    173.49, 0.69,  0.00,  0.08,  0.35,  0.69,  0.00,  0.00,  0.00,  0.15,
    0.06,   0.24,  0.06,  0.13,  0.00,  0.13,  0.58,  0.49,  0.00,  0.06,
    0.22,   0.24,  0.19,  0.00,  0.20,  0.15,
    172.47, 0.17,  0.00,  0.03,  0.31,  0.53,  0.00,  0.12,  0.07,  0.09,
    0.04,   0.00,  0.15,  0.02,  0.01,  0.01,  0.00,  0.28,  0.69,  0.18,
    0.87,   0.13,  0.24,  0.27,  0.04,  0.16,
    172.04, 0.11,  0.15,  0.01,  0.01,  0.12,  0.00,  0.03,  0.09,  0.26,
    0.02,   0.00,  0.27,  0.00,  0.11,  0.10,  0.00,  0.39,  0.17,  0.09,
    0.13,   0.15,  0.14,  0.20,  0.17,  0.15,
    172.08, 0.27,  0.24,  0.01,  0.01,  0.25,  0.00,  0.00,  0.00,  0.10,
    0.02,   0.18,  0.31,  0.22,  0.00,  0.06,  0.00,  0.31,  0.14,  0.02,
    0.19,   0.20,  0.21,  0.08,  0.28,  0.25
   };

  // Correlations:
  //--------------
  // rho_year <==> pairs in -different- decay channels but same year -> Rho_Dif
  // rho_chan <==> pairs in - the same- decay channels but diff year -> Rho_Sam
  // Rho_Dif = 0, Rho_Sam = 0 for: 0, 1, 5, 10, 12, 14 
  // Rho_Dif = 1, Rho_Sam = 1 rho: 2, 3, 6-9, 13, 15-24
  // Rho_Dif = 1, Rho_Sam = 1 for: 4, 11

  // rho_year 
  Double_t RhoDif[NumUnc] = {
    //0    1    2    3    4    5    6    7    8    9   10   11   12
    0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    //13  14   15   16   17   18   19   20   21   22   23   24
    1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  // rho_chan
  Double_t RhoSam[NumUnc] = {
    //0    1    2    3    4    5    6    7    8    9   10   11   12
    0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,
    //13  14   15   16   17   18   19   20   21   22   23   24    
    1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  // Change correlations according to Flag, for systematics on page 10
  if(Flag == 2){RhoDif[4] = 1.0;
  }else if(Flag == 3){RhoDif[16] = 0.0;
  }else if(Flag == 4){RhoDif[ 8] = 0.0;
  }

  //-- Local Structures for Blue input
  // TMatrices
  TMatrixD *const RhoSou = new TMatrixD(NumEst,NumEst);
  //-- End

  //-- Local Structures for Blue output
  // TMatrices
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRhoRes = new TMatrixD(NumObs,NumObs);
  //-- End

  const TString FilBas = "B_CMS_PAS_2014_015";
  TString FilNam;
  char Buffer[100];

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%6.4f";
  const TString ForRho = "%5.2f";
  const TString ForPul = "%4.2f";
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";

  // Construct Object
  Blue* myBlue;
  if(Flag == 5){myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  }else{myBlue = new Blue(NumEst, NumUnc);
  }
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);

  // Fill estimates and correlations
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }
  for(Int_t k = 0; k<NumUnc; k++){
    FilMat(RhoDif[k], RhoSam[k], RhoSou);
    //printf("... CMS_PAS_2014_015: Source k = %2i \n", k);
    //myBlue->PrintMatrix(RhoSou);
    myBlue->FillCor(k,RhoSou);
  }

  // Set reduced correlations but for Flag == 1
  if(Flag == 1){
    printf("... B_CMS_PAS_2014_015: Nominal correlations of the estimates\n");
  }else{
    printf("... B_CMS_PAS_2014_015: Reduced correlations of the estimates\n");
    for(Int_t k = 0; k<NumUnc; k++){
      if(RhoDif[k]+RhoSam[k] > 0)myBlue->SetRhoRedUnc(k);
    }
  }

  // Do different things according to Flag
  if(Flag <= 4){
    myBlue->FixInp();

    // The CMS default
    if(Flag == 0){
      myBlue->PrintCompatEst();
      myBlue->PrintEst();
      printf("... B_CMS_PAS_2014_015: For the correlation of the input");
      printf(" see Table 4 \n");
      myBlue->GetRho(LocRho);
      myBlue->PrintMatrix(LocRho, ForRho);
    }
    myBlue->Solve();
    myBlue->PrintResult();

    // Write files for CMS default and nominal correlations
    if(Flag <= 1){
      // Get the output into files
      sprintf(Buffer,"%1i", Flag);
      FilNam = &Buffer[0];       
      FilNam = FilBas + "_" + FilNam;
      myBlue->LatexResult(FilNam);
      myBlue->DisplayResult(0,FilNam);

      // Solve according to importance
      myBlue->ReleaseInp();
      myBlue->FixInp();
      myBlue->SolveAccImp(1.0);
      printf("... B_CMS_PAS_2014_015: Result of the successive combination \n");
      myBlue->DisplayAccImp(0,FilNam);
    }

  }else if(Flag == 5){

    myBlue->FixInp();
    myBlue->PrintStatus();
    myBlue->Solve();
    printf("... B_CMS_PAS_2014_015: The combination for l+jets,");
    printf(" di-lepton and all-jets, Flag = %2i \n",Flag);
    printf("... B_CMS_PAS_2014_015: For the results see Table 3 \n");
    myBlue->PrintResult();
    myBlue->GetRhoRes(LocRhoRes);
    LocRhoRes->operator*=(100.);
    myBlue->PrintMatrix(LocRhoRes, "%4.1f%%");
    myBlue->PrintCompatObs();

    // Get the output into files
    sprintf(Buffer,"%1i", Flag);
    FilNam = &Buffer[0];       
    FilNam = FilBas + "_" + FilNam;
    myBlue->LatexResult(FilNam);
    for(Int_t n = 0; n<NumObs; n++)myBlue->DisplayResult(n,FilNam);

  }
  
  // Delete object and matrices
  delete myBlue;  myBlue = NULL;

  LocRho->Delete(); LocRho = NULL;
  LocRhoRes->Delete(); LocRhoRes = NULL;

  // Return
  return;
};

//------------------------------------------------------------------------
// Utility to fill a correlation matrix for any pair of RhoDif and RhoSam
//------------------------------------------------------------------------

void FilMat(const Double_t RD, const Double_t RS, TMatrixD *const M){

  // Fill upper half of the correlation matrix
  M->operator()(0,0) =  1;
  M->operator()(0,1) = RS;
  M->operator()(0,2) = RD;
  M->operator()(0,3) = RD;
  M->operator()(0,4) = RD;
  M->operator()(0,5) = RD;
  M->operator()(0,6) = RD;
  M->operator()(0,7) = RD;

  M->operator()(1,1) =  1;
  M->operator()(1,2) = RD;
  M->operator()(1,3) = RD;
  M->operator()(1,4) = RD;
  M->operator()(1,5) = RD;
  M->operator()(1,6) = RD;
  M->operator()(1,7) = RD;

  M->operator()(2,2) =  1;
  M->operator()(2,3) = RS;
  M->operator()(2,4) = RS;
  M->operator()(2,5) = RD;
  M->operator()(2,6) = RD;
  M->operator()(2,7) = RD;

  M->operator()(3,3) =  1;
  M->operator()(3,4) = RS;
  M->operator()(3,5) = RD;
  M->operator()(3,6) = RD;
  M->operator()(3,7) = RD;

  M->operator()(4,4) =  1;
  M->operator()(4,5) = RD;
  M->operator()(4,6) = RD;
  M->operator()(4,7) = RD;

  M->operator()(5,5) =  1;
  M->operator()(5,6) = RS;
  M->operator()(5,7) = RS;

  M->operator()(6,6) =  1;
  M->operator()(6,7) = RS;

  M->operator()(7,7) =  1;

  // Now copy to lower half
  const Int_t NN = M->GetNrows();
  for(Int_t i = 0; i<NN; i++){
    for(Int_t j = i; j<NN; j++){
      M->operator()(j,i) = M->operator()(i,j);
    }
  }
  return;
};
