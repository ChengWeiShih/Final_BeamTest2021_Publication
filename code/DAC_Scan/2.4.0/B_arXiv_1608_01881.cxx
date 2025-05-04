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

void B_arXiv_1608_01881(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // [Shows how to use the Logo with SetLogo()]
  // Flag steers which of the results should be calculated
  //  0: The Tevatron mtop combination (August 5, 2016)
  //  1: Suggest a combination according to importance up to the point where
  //     each added measurement adds less than 1% in precision
  //  2: The combination for four observables, namely M(hadronic),
  //     M(lepton+jets), M(di-lepton) and M(Missing Et)
  //  3: The combination for two correlated observables: M(CDF), M(DO)
  //  4: The combination for two correlated observables: M(RUN-I), M(RUN-II)
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 12;
  static const Int_t NumUnc = 16;
  static const Int_t MaxObs =  4;
  Int_t NumObs = 1;

  // The names of estimates, uncertainties and observables
  //                               0          1          2          3
  TString NamEst[NumEst] = {"CI  l+j", "CI  dil", "CI  had", "DI  l+j", 
  //                               4          5          6          7
			    "DI  dil", "CII l+j", "CII Lxy", "CII MEt", 
  //                               8          9         10         11
			    "CII dil", "CII had", "DII l+j", "DII dil"};
  //                               0          1          2          3

  TString NamUnc[NumUnc] = {"   Stat", "   iJES", "   aJES", "   bJES", 
  //                               4          5          6          7
			    "   cJES", "   rJES", "   dJES", "  LepPt", 
  //                               8          9         10         11
			    " Signal", " DeTMod", "  b-tag", "   BGMC", 
  //                              12         13         14         15
			    " BGData", "   Meth", "  UN/MI", "    MHI"};

  TString NamObs[MaxObs] = {"   mtop", "   mtop", "   mtop", "   mtop"};
  
  // Index for which estimate determines which observable
  Int_t IWhichObs[NumEst] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // Flag: steers which of the results should be calculated
  if(Flag == 0 || Flag == 1){
    printf("... B_arXiv_1608_01881: -----------------------------------");
    printf("---------- \n");
    printf("... B_arXiv_1608_01881: The 2016 Tevatron mtop combination,");
    printf(" Flag = %2i \n", Flag);
    NumObs = 1;
  }else if(Flag == 2){
    printf("... B_arXiv_1608_01881: -----------------------------------------");
    printf("----------- \n");
    printf("... B_arXiv_1608_01881: The mtop combination for four");
    printf(" observables, Flag = %2i \n",Flag);
    // 0 = had, 1=l+j, 2=dil, 3=MEt
    NamObs[0] = " mt-had";
    NamObs[1] = " mt-l+j";
    NamObs[2] = " mt-dil";
    NamObs[3] = " mt-MEt";
    NumObs = 4;
    IWhichObs[ 0] = 1;
    IWhichObs[ 1] = 2;
    IWhichObs[ 2] = 0;
    IWhichObs[ 3] = 1;
    IWhichObs[ 4] = 2;
    IWhichObs[ 5] = 1;
    IWhichObs[ 6] = 1;
    IWhichObs[ 7] = 3;
    IWhichObs[ 8] = 2;
    IWhichObs[ 9] = 0;
    IWhichObs[10] = 1;
    IWhichObs[11] = 2;
  }else if(Flag == 3){
    printf("... B_arXiv_1608_01881: ----------------------------------------");
    printf("------------------------- \n");
    printf("... B_arXiv_1608_01881: The mtop combination for two observables");
    printf(" M(CDF), M(DO), Flag = %2i \n", Flag);
    // 0 = CDF, 1=DO
    NamObs[0] = " mt-CDF";
    NamObs[1] = "  mt-D0";
    NumObs = 2;
    IWhichObs[ 0] = 0;
    IWhichObs[ 1] = 0;
    IWhichObs[ 2] = 0;
    IWhichObs[ 3] = 1;
    IWhichObs[ 4] = 1;
    IWhichObs[ 5] = 0;
    IWhichObs[ 6] = 0;
    IWhichObs[ 7] = 0;
    IWhichObs[ 8] = 0;
    IWhichObs[ 9] = 0;
    IWhichObs[10] = 1;
    IWhichObs[11] = 1;
  }else if(Flag == 4){
    printf("... B_arXiv_1608_01881: ----------------------------------------");
    printf("------------------------------- \n");
    printf("... B_arXiv_1608_01881: The mtop combination for two observables");
    printf(" M(RUN-I), M(RUN-II), Flag = %2i \n", Flag);
    // 0 = RUN-I, 1 = RUN-II
    NamObs[0] = " mt-RUN-I";
    NamObs[1] = "mt-RUN-II";
    NumObs = 2;
    IWhichObs[ 0] = 0;
    IWhichObs[ 1] = 0;
    IWhichObs[ 2] = 0;
    IWhichObs[ 3] = 0;
    IWhichObs[ 4] = 0;
    IWhichObs[ 5] = 1;
    IWhichObs[ 6] = 1;
    IWhichObs[ 7] = 1;
    IWhichObs[ 8] = 1;
    IWhichObs[ 9] = 1;
    IWhichObs[10] = 1;
    IWhichObs[11] = 1;
  }else{
    printf("... B_arXiv_1608_01881: Not implemented Flag = %2i \n", Flag);
    return;
  }
  
  // Estimates
  // The order of estimates follows Table 1 of the reference
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  // Num:       0      1      2      3      4      5      6      7      8      9
  //         Stat   iJES   aJES   bJES   cJES   rJES   dJES  LepPt Signal DeTMod
  // Rho:       0  Cor01  Cor02      1      1  Cor05  Cor02  Cor02      1  Cor05
  //      10   11     12     13     14     15
  //   b-tag BGMC BGData   Meth  UN/MI    MHI
  //  Cor02 Cor11  Cor12      0   Cor05  Cor02
  Double_t XEst[LenXEst] = {
    176.10,  5.10,  0.00,  0.00,  0.60,  2.70,  3.35,  0.70,  0.00,  2.62, 0.00,
    0.40,    1.30,  0.00,  0.00,  0.00,  0.00,
    167.40, 10.30,  0.00,  0.00,  0.80,  2.60,  2.65,  0.60,  0.00,  2.86, 0.00,
    0.00,    0.30,  0.00,  0.70,  0.00,  0.00,
    186.00, 10.00,  0.00,  0.00,  0.60,  3.00,  4.00,  0.30,  0.00,  1.97, 0.00,
    0.00,    0.00,  1.70,  0.60,  0.00,  0.00,
    180.10,  3.60,  0.00,  0.00,  0.71,  2.00,  0.00,  2.53,  0.00,  1.10, 0.00,
    0.00,    1.00,  0.00,  0.58,  1.30,  0.00,
    168.40, 12.30,  0.00,  0.00,  0.71,  2.00,  0.00,  1.12,  0.00,  1.80, 0.00,
    0.00,    1.10,  0.00,  1.14,  1.30,  0.00,
    172.85,  0.52,  0.49,  0.09,  0.16,  0.21,  0.48,  0.07,  0.03,  0.61, 0.00,
    0.03,    0.12,  0.16,  0.05,  0.00,  0.07,
    166.90,  9.00,  0.00,  0.00,  0.00,  0.36,  0.24,  0.06,  0.00,  0.90, 0.00,
    0.00,    0.80,  0.20,  2.50,  0.00,  0.00,
    173.93,  1.26,  1.05,  0.10,  0.17,  0.18,  0.40,  0.04,  0.00,  0.63, 0.00,
    0.03,    0.00,  0.15,  0.21,  0.00,  0.18,
    171.50,  1.91,  0.00,  0.16,  0.26,  1.47,  1.56,  0.37,  0.41,  1.01, 0.00,
    0.05,    0.24,  0.31,  0.20,  0.00,  0.27,
    175.07,  1.19,  0.97,  0.01,  0.20,  0.37,  0.42,  0.09,  0.00,  0.53, 0.00,
    0.04,    0.00,  0.15,  0.87,  0.00,  0.22,
    174.98,  0.41,  0.41,  0.16,  0.09,  0.00,  0.00,  0.21,  0.01,  0.35, 0.07,
    0.10,    0.06,  0.09,  0.07,  0.00,  0.06,
    173.50,  1.31,  0.47,  0.28,  0.13,  0.00,  0.00,  0.31,  0.08,  0.43, 0.14,
    0.22,    0.00,  0.08,  0.14,  0.00,  0.07
  };
  
  // Dimension for correlation matrices
  static const Int_t LenCor = NumEst * NumEst;

  // Correlation matrix Cor01, for uncertainty 1
  // C1   C1   C1   D1   D1   C2   C2   C2   C2   C2   D2   D2
  Double_t Cor01[LenCor] = {
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0
  };

  // Correlation matrix Cor02, for uncertainties 2, 6, 7, 10, 15
  // C1   C1   C1   D1   D1   C2   C2   C2   C2   C2   D2   D2
  Double_t Cor02[LenCor] = {
    // RUNI CDF
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    // RUNI D0
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    // RUNII CDF
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    // RUNII D0
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
  };

  // Correlation matrix Cor05, for uncertainties 5, 9, 14
  // C1   C1   C1   D1   D1   C2   C2   C2   C2   C2   D2   D2
  Double_t Cor05[LenCor] = {
    // CDF
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    // D0
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    // CDF
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    // D0
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
  };

  // Correlation matrix Cor11, for uncertainty 11
  // C1   C1   C1   D1   D1   C2   C2   C2   C2   C2   D2   D2
  //  1    2    0    1    2    1    1    3    2    0    1    2
  Double_t Cor11[LenCor] = {
    //L+j
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    //dil
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    //had
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    //L+j
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    //dil
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    //L+j
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    //Met
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    //dil
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    //had
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    //L+j
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    //dil
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0
  };
  
  // Correlation matrix Cor12, for uncertainty 12
  // C1   C1   C1   D1   D1   C2   C2   C2   C2   C2   D2   D2
  //  1    2    0    1    2    1    1    3    2    0    1    2
  Double_t Cor12[LenCor] = {
    // L+j CDF RUN-I
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    // Dil CDF RUN-I
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    // Had CDF RUN-I
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    // L+j D0 RUN-I
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    // Dil D0 RUN-I
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    // L+j CDF RUN-II
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    // Met CDF RUN-II
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    // Dil CDF RUN-II
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    // Had CDF RUN-II
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    // L+j D0 RUN-II
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    // Dil D0 RUN-II
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
  };
  
  // The values from the paper combination
  const Double_t TopVal = 174.30;
  const Double_t TopSta =   0.35;
  const Double_t TopSys =   0.54;
  const Double_t TopFul =   0.65;


  //-- Local Structures for BLUE output
  // TMatrices
  TMatrixD*    LocRho  = new TMatrixD(NumEst,NumEst);
  TMatrixD*    LocWei  = new TMatrixD(NumEst,NumObs);
  TMatrixD* LocRhoRes  = new TMatrixD(NumObs,NumObs);
  //-- End

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%5.2f";
  const TString ForRho = ForVal;
  const TString ForPul = ForVal;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";
  const TString FilBas = "B_arXiv_1608_01881";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->PrintStatus();
  myBlue->SetLogo("BLUE", "2.4.0", kBlue, 0.08, 0.92);
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

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

  // Fill all correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(k ==  0 || k == 13){
      myBlue->FillCor(k, 0.0);
    }else if(k == 1){
      myBlue->FillCor(k, &Cor01[0]);
    }else if(k == 2 || k == 6 || k ==  7 || k == 10 || k == 15){
      myBlue->FillCor(k, &Cor02[0]);
    }else if(k == 5 || k == 9 || k == 14){
      myBlue->FillCor(k, &Cor05[0]);
    }else if(k == 11){
      myBlue->FillCor(k, &Cor11[0]);
    }else if(k == 12){
      myBlue->FillCor(k, &Cor12[0]);
    }else{
      myBlue->FillCor(k, 1.0);
    }
  }
  
  // Fix input, solve depending on Flag, and finally delete
  if(Flag == 0){
    myBlue->FixInp();
    myBlue->PrintCompatEst("B_arXiv_1608_01881_0_All");
    myBlue->Solve();
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - full combination,");
    printf(" for correlations see Table 2 \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->GetRho(LocRho);
    myBlue->PrintMatrix(LocRho,"%+5.2f");

    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - full combination,");
    printf(" for results see Table 3 \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintResult();

    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - full combination,");
    printf(" for pulls see Table 4\n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintPull();

    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - full combination,");
    printf(" for weights see Table 4\n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->GetWeight(LocWei);
    LocWei->operator*=(100.);
    myBlue->PrintMatrix(LocWei," %+4.0f%%");

    myBlue->LatexResult(FilBas,ForVal,ForUnc,ForWei,ForRho,ForPul);
    myBlue->DisplayResult(0,FilBas);

    // Independent CDF Combination
    myBlue->ReleaseInp();
    for(Int_t i =  3; i <=  4; i++)myBlue->SetInActiveEst(i);
    for(Int_t i = 10; i <= 11; i++)myBlue->SetInActiveEst(i);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - independent CDF combination. \n");
    printf("... B_arXiv_1608_01881: Table 6 reports the two observable");
    printf(" combination, see Flag == 3 \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintResult();

    // Independent D0 Combination
    myBlue->ResetInp();
    for(Int_t i = 0; i <= 2; i++)myBlue->SetInActiveEst(i);
    for(Int_t i = 5; i <= 9; i++)myBlue->SetInActiveEst(i);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - independent D0 combination. \n");
    printf("... B_arXiv_1608_01881: Table 6 reports the two observable,");
    printf(" combination, see Flag == 3 \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintResult();

    // Independent RUN-I Combination
    myBlue->ResetInp();
    for(Int_t i = 5; i <= 11; i++)myBlue->SetInActiveEst(i);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - independent RUN-I combination,");
    printf(" not in reference \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintResult();

    // Independent RUN-II Combination
    myBlue->ResetInp();
    for(Int_t i = 0; i <= 4; i++)myBlue->SetInActiveEst(i);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - independent RUN-II combination,");
    printf(" not in reference \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintResult();

  }else if(Flag == 1){

    myBlue->FixInp();
    myBlue->SolveAccImp(0,1.0);
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: The results of the combination");
    printf(" according to importance. \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->DisplayAccImp(0,"B_arXiv_1608_01881_1");

    // Get the indices from SolveAccImp()
    Int_t ImpInd[NumEst] = {0};
    myBlue->GetAccImpIndEst(0,ImpInd);
    Int_t ImpLas = myBlue->GetAccImpLasEst(0);
    Int_t ImpFir = ImpInd[0];
    Int_t ImpNum = 0;
    for(Int_t i = 0; i<NumEst; i++)if(ImpInd[i] == ImpLas)ImpNum = i;
    printf("... B_arXiv_1608_01881: The list of importance is:");
    for(Int_t i = 0; i<NumEst; i++)printf(" %3i", ImpInd[i]);
    printf("\n");
    printf("... B_arXiv_1608_01881: The first and the last of the");
    printf(" %3i estimates are %3i and %3i \n", ImpNum+1, ImpFir, ImpLas);

    // Disable all insignificant estimates
    myBlue->ReleaseInp();
    for(Int_t i = ImpNum+1; i<NumEst; i++)myBlue->SetInActiveEst(ImpInd[i]);
    myBlue->FixInp();

    // Now solve for the selected estimates
    myBlue->PrintCompatEst("B_arXiv_1608_01881_1_Sel");
    myBlue->Solve();
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: The result based");
    printf(" on the selected estimates. \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintResult();

    // Show the paper result
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: The paper combination is:\n");
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881:%s",myBlue->GetNamObs(0).Data());
    printf(" = %5.2f (+- %4.2f +- %4.2f) = +- %4.2f\n",
	   TopVal, TopSta, TopSys, TopFul);
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: The difference is marginal! \n");
    printf("... B_arXiv_1608_01881:\n");

    myBlue->DisplayResult(0,"B_arXiv_1608_01881_1");

    // Get the weights
    printf("... B_arXiv_1608_01881: The Blue weights. \n");
    myBlue->GetWeight(LocWei);
    LocWei->operator*=(100.);
    myBlue->PrintMatrix(LocWei,ImpNum+1,1," %+4.1f%%");

    // Scan the rho parameters for the selected estimates
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: The individual scans");
    printf(" of the correlations of the selected estimates. \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->ReleaseInp();
    myBlue->FixInp();
    myBlue->SolveScaRho(0);
    myBlue->PrintScaRho("B_arXiv_1608_01881_1");

  }else if(Flag == 2){

    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop -");
    printf(" (4-observables had, l+j, dil and MEt) \n");
    printf("... B_arXiv_1608_01881: for the results see Table 5. \n");
    printf("... B_arXiv_1608_01881:\n");
    for(Int_t n = 0; n<NumObs; n++){
      myBlue->DisplayResult(n,"B_arXiv_1608_01881_2");
    }
    myBlue->PrintResult();
    myBlue->PrintCompatObs();
    myBlue->GetRhoRes(LocRhoRes);
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: The correlations of the observables,");
    printf(" see Table 5. \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintMatrix(LocRhoRes,"%+5.2f");

  }else if(Flag == 3){

    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintCompatObs();
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - (two correlated observables");
    printf(" CDF and D0)\n");
    printf("... B_arXiv_1608_01881: For the results, see Table 6. \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintResult();

    myBlue->GetRhoRes(LocRhoRes);
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: The correlation of the observables. \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintMatrix(LocRhoRes,"%+5.2f");

  }else if(Flag == 4){

    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: mtop - (two correlated observables");
    printf(" RUN-I and RUN-II) \n");
    printf("... B_arXiv_1608_01881: Not provided in reference. \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintResult();

    myBlue->GetRhoRes(LocRhoRes);
    printf("... B_arXiv_1608_01881:\n");
    printf("... B_arXiv_1608_01881: The correlation of the observables. \n");
    printf("... B_arXiv_1608_01881:\n");
    myBlue->PrintMatrix(LocRhoRes,"%+5.2f");
  }

  // Delete objects
  delete myBlue; myBlue = NULL; 

  LocRho->Delete(); LocRho = NULL;
  LocWei->Delete(); LocWei = NULL;
  LocRhoRes->Delete(); LocRhoRes = NULL;

  // Return
  return;
}
