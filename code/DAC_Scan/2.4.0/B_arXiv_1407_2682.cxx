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

void B_arXiv_1407_2682(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The Tevatron mtop combination (v2 July 12, 2014)
  //  1: Suggest a combination according to importance up to the point where
  //     each added measurement adds less than 1% in precision
  //  2: The combination for four observables, namely M(all-had),
  //     M(lepton+jest), M(di-lepton) and M(missing Et)
  //  3: The combination for two correlated observables: M(CDF), M(DO)
  //  4: The combination for two correlated observables: M(RUN-I), M(RUN-II)
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 12;
  static const Int_t NumUnc = 16;
  static const Int_t MaxObs =  4;
  Int_t NumObs = 1;

  // The names of estimates, uncertainties and observables
  TString NamEst[NumEst] = {"CI  l+j", "CI  dil", "CI  had", "DI  l+j", 
			    "DI  dil", "CII l+j", "CII Lxy", "CII Met", 
			    "DII l+j", "DII dil", "CII dil", "CII had"};
  TString NamUnc[NumUnc] = {"   Stat", "   iJES", "   aJES", "   bJES", 
			    "   cJES", "   rJES", "   dJES", "  LepPt", 
			    " Signal", " DeTMod", "  b-tag", "   BGMC", 
			    " BGData", "   Meth", "  UN/MI", 
			    "    MHI"};
  TString NamObs[MaxObs] = {"   mtop", "   mtop", "   mtop", "   mtop"};
  
  // Index for which estimate determines which observable
  Int_t IWhichObs[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // Flag: steers which of the results should be calculated
  if(Flag == 0 || Flag == 1){
    printf("... B_arXiv_1407_2682: -----------------------------------");
    printf("---------- \n");
    printf("... B_arXiv_1407_2682: The 2013 Tevatron mtop combination,");
    printf(" Flag = %2i \n", Flag);
    NumObs = 1;
  }else if(Flag == 2){
    printf("... B_arXiv_1407_2682: ------------------------------------------");
    printf("---------- \n");
    printf("... B_arXiv_1407_2682: The mtop combination for four observables,");
    printf(" Flag = %2i \n",Flag);
    // 0 = allh, 1=l+j, 2=di-l, 3=Met
    NamObs[0] = "mt-allh";
    NamObs[1] = " mt-l+j";
    NamObs[2] = "mt-di-l";
    NamObs[3] = " mt-Met";
    NumObs = 4;
    IWhichObs[ 0] = 1;
    IWhichObs[ 1] = 2;
    IWhichObs[ 2] = 0;
    IWhichObs[ 3] = 1;
    IWhichObs[ 4] = 2;
    IWhichObs[ 5] = 1;
    IWhichObs[ 6] = 1;
    IWhichObs[ 7] = 3;
    IWhichObs[ 8] = 1;
    IWhichObs[ 9] = 2;
    IWhichObs[10] = 2;
    IWhichObs[11] = 0;
  }else if(Flag == 3){
    printf("... B_arXiv_1407_2682: ----------------------------------------");
    printf("------------------------- \n");
    printf("... B_arXiv_1407_2682: The mtop combination for two observables");
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
    IWhichObs[ 8] = 1;
    IWhichObs[ 9] = 1;
    IWhichObs[10] = 0;
    IWhichObs[11] = 0;
  }else if(Flag == 4){
    printf("... B_arXiv_1407_2682: ----------------------------------------");
    printf("------------------------------- \n");
    printf("... B_arXiv_1407_2682: The mtop combination for two observables");
    printf(" M(RUN-I), M(RUN-II), Flag = %2i \n", Flag);
    // 0 = RUN-I, 1=RUN-II
    NamObs[0] = " mt-RUI";
    NamObs[1] = "mt-RUII";
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
    printf("... B_arXiv_1407_2682: Not implemented Flag = %2i \n", Flag);
    return;
  }
  
  // Estimates
  // The order of estimates follows Table 2 of the reference
  // What: 'Stat' 'iJES' 'aJES' 'bJES' 'cJES' 'rJES' 'dJES' 'Lept' 'Sign' 'DTMO' 'btag' 'BGMC' 'BGDT' 'Meth' 'UN/MI'  'MHI'
  // Num:      0      1      2      3      4      5      6      7      8      9     10     11     12     13      14     15
  // Rho:      0  Cor01  Cor02      1      1  Cor05  Cor02  Cor02      1  Cor05  Cor02  Cor11  Cor12      0   Cor05  Cor02
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    176.10,  5.10,  0.00,  0.00,  0.60,  2.70,  3.40,  0.70,  0.00,  2.60,  0.00,  0.40,  1.30,  0.00,  0.00,  0.00,  0.00,
    167.40, 10.30,  0.00,  0.00,  0.80,  2.60,  2.70,  0.60,  0.00,  2.90,  0.00,  0.00,  0.30,  0.00,  0.70,  0.00,  0.00,
    186.00, 10.00,  0.00,  0.00,  0.60,  3.00,  4.00,  0.30,  0.00,  2.00,  0.00,  0.00,  0.00,  1.70,  0.60,  0.00,  0.00,
    180.10,  3.60,  0.00,  0.00,  0.70,  2.00,  0.00,  2.50,  0.00,  1.10,  0.00,  0.00,  1.00,  0.00,  0.60,  1.30,  0.00,
    168.40, 12.30,  0.00,  0.00,  0.70,  2.00,  0.00,  1.10,  0.00,  1.80,  0.00,  0.00,  1.10,  0.00,  1.10,  1.30,  0.00,
    172.85,  0.52,  0.49,  0.09,  0.16,  0.21,  0.48,  0.07,  0.03,  0.61,  0.00,  0.03,  0.12,  0.16,  0.05,  0.00,  0.07,
    166.90,  9.00,  0.00,  0.00,  0.00,  0.36,  0.24,  0.06,  0.00,  0.90,  0.00,  0.00,  0.80,  0.20,  2.50,  0.00,  0.00,
    173.93,  1.26,  1.05,  0.10,  0.17,  0.18,  0.40,  0.04,  0.00,  0.63,  0.00,  0.03,  0.00,  0.15,  0.21,  0.00,  0.18,
    174.98,  0.41,  0.41,  0.16,  0.09,  0.00,  0.00,  0.21,  0.01,  0.35,  0.07,  0.10,  0.06,  0.09,  0.07,  0.00,  0.06,
    174.00,  2.36,  0.55,  0.40,  0.20,  0.00,  0.00,  0.56,  0.35,  0.86,  0.50,  0.00,  0.00,  0.20,  0.51,  0.00,  0.00,
    170.80,  1.83,  0.00,  0.18,  0.28,  1.65,  1.72,  0.46,  0.36,  0.96,  0.00,  0.05,  0.30,  0.33,  0.19,  0.00,  0.30,
    175.07,  1.19,  0.97,  0.02,  0.20,  0.37,  0.42,  0.09,  0.00,  0.53,  0.00,  0.04,  0.00,  0.15,  0.87,  0.00,  0.22
  };
  
  // Dimension for correlation matrices
  static const Int_t LenCor = NumEst * NumEst;

  // Correlation matrix Cor01, for uncertainties 1
  // C1   C1   C1   D1   D1   C2   C2   C2   D2   D2   C2   C2
  Double_t Cor01[LenCor] = {
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
  };

  // Correlation matrix Cor02, for uncertainties 2, 6, 7, 10, 15
  // C1   C1   C1   D1   D1   C2   C2   C2   D2   D2   C2   C2
  Double_t Cor02[LenCor] = {
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0
  };

  // Correlation matrix Cor05, for uncertainties 5, 9, 14
  // C1   C1   C1   D1   D1   C2   C2   C2   D2   D2   C2   C2
  Double_t Cor05[LenCor] = {
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0
  };

  // Correlation matrix Cor11, for uncertainty 11
  // C1   C1   C1   D1   D1   C2   C2   C2   D2   D2   C2   C2
  //  1    2    0    1    2    1    1    3    1    2    2    0
  Double_t Cor11[LenCor] = {
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
  };
  
  // Correlation matrix Cor12, for uncertainty 12
  // C1   C1   C1   D1   D1   C2   C2   C2   D2   D2   C2   C2
  //  1    2    0    1    2    1    1    3    1    2    2    0
  Double_t Cor12[LenCor] = {
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
  };
  
  //-- Local Structures for BLUE output
  // TMatrices
  TMatrixD*    LocRho  = new TMatrixD(NumEst,NumEst);
  TMatrixD*    LocWei  = new TMatrixD(NumEst,NumObs);
  TMatrixD* LocRhoRes  = new TMatrixD(NumObs,NumObs);
  //-- End

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%5.3f";
  const TString ForRho = ForVal;
  const TString ForPul = ForVal;
  const TString ForChi = ForWei;
  const TString ForUni = "GeV";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->PrintStatus();
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
      myBlue->FillCor(k,0.0);
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
      myBlue->FillCor(k,1.0);
    }
  }
  
  // Fix input, solve depending on Flag, and finally delete
  if(Flag == 0){
    myBlue->FixInp();
    myBlue->PrintCompatEst();
    myBlue->Solve();
    printf("... B_arXiv_1407_2682: mtop full combination,");
    printf(" for correlations see Table 2 \n");
    myBlue->GetRho(LocRho);
    myBlue->PrintMatrix(LocRho,"%+5.2f");

    printf("... B_arXiv_1407_2682: mtop full combination,");
    printf(" for results see Table 3 \n");
    myBlue->PrintResult();

    printf("... B_arXiv_1407_2682: mtop full combination pulls, Table 4\n");
    myBlue->PrintPull();

    printf("... B_arXiv_1407_2682: mtop full combination weights, Table 4\n");
    myBlue->GetWeight(LocWei);
    LocWei->operator*=(100.);
    myBlue->PrintMatrix(LocWei," %+4.1f%%");

    myBlue->LatexResult("B_arXiv_1407_2682");
    myBlue->DisplayResult(0,"B_arXiv_1407_2682");

    // Independent CDF Combination
    myBlue->ReleaseInp();
    myBlue->SetInActiveEst(3);
    myBlue->SetInActiveEst(4);
    myBlue->SetInActiveEst(8);
    myBlue->SetInActiveEst(9);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1407_2682: mtop independent CDF combination. \n");
    printf("... B_arXiv_1407_2682: Table 6 reports the correlated");
    printf(" combination, see Flag == 3 \n");
    myBlue->PrintResult();

    // Independent D0 Combination
    myBlue->ResetInp();
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(1);
    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(5);
    myBlue->SetInActiveEst(6);
    myBlue->SetInActiveEst(7);
    myBlue->SetInActiveEst(10);
    myBlue->SetInActiveEst(11);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1407_2682: mtop independent D0 combination. \n");
    printf("... B_arXiv_1407_2682: Table 6 reports the correlated");
    printf(" combination, see Flag == 3 \n");
    myBlue->PrintResult();

    // Independent RUN-I Combination
    myBlue->ResetInp();
    myBlue->SetInActiveEst(5);
    myBlue->SetInActiveEst(6);
    myBlue->SetInActiveEst(7);
    myBlue->SetInActiveEst(8);
    myBlue->SetInActiveEst(9);
    myBlue->SetInActiveEst(10);
    myBlue->SetInActiveEst(11);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1407_2682: mtop RUN-I combination");
    printf(", not in reference\n");
    myBlue->PrintResult();

    // Independent RUN-II Combination
    myBlue->ResetInp();
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(1);
    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(3);
    myBlue->SetInActiveEst(4);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1407_2682: mtop RUN-II combination");
    printf(", not in reference\n");
    myBlue->PrintResult();

  }else if(Flag == 1){

    myBlue->FixInp();
    myBlue->PrintCompatEst("B_arXiv_1407_2682_1_All");
    myBlue->CorrelPair(5,10,"B_arXiv_1407_2682");
    myBlue->SolveAccImp(0,1.0);
    printf("... B_arXiv_1407_2682: The results of the successive");
    printf(" combination \n");
    myBlue->DisplayAccImp(0,"B_arXiv_1407_2682_1");

    // Get the indices from SolveAccImp()
    Int_t ImpInd[NumEst] = {0};
    myBlue->GetAccImpIndEst(0,ImpInd);
    Int_t ImpLas = myBlue->GetAccImpLasEst(0);
    Int_t ImpFir = ImpInd[0];
    Int_t ImpNum = 0;
    for(Int_t i = 0; i<NumEst; i++)if(ImpInd[i] == ImpLas)ImpNum = i;
    printf("... B_arXiv_1407_2682: The list of importance is:");
    for(Int_t i = 0; i<NumEst; i++)printf(" %3i", ImpInd[i]);
    printf("\n");
    printf("... B_arXiv_1407_2682: The first and the last of the");
    printf(" %3i estimates are %3i and %3i \n", ImpNum+1, ImpFir, ImpLas);

    // Disable all insignificant estimates
    myBlue->ReleaseInp();
    for(Int_t i = ImpNum+1; i<NumEst; i++)myBlue->SetInActiveEst(ImpInd[i]);
    myBlue->FixInp();

    // Now solve for the selected estimates
    myBlue->PrintCompatEst("B_arXiv_1407_2682_1_Sel");
    myBlue->Solve();
    printf("... B_arXiv_1407_2682: The result based on the selected");
    printf(" estimates\n");
    myBlue->PrintResult();
    myBlue->DisplayResult(0,"B_arXiv_1407_2682_1");

    // Get the weights
    printf("... B_arXiv_1407_2682: The Blue weights \n");
    myBlue->GetWeight(LocWei);
    LocWei->operator*=(100.);
    myBlue->PrintMatrix(LocWei,ImpNum+1,1," %+4.1f%%");

    // Scan the rho parameters for the selected estimates
    printf("... B_arXiv_1407_2682: The individual scans");
    printf(" of the correlations of the selected estimates \n");
    myBlue->ReleaseInp();
    myBlue->FixInp();
    myBlue->SolveScaRho(0);
    myBlue->PrintScaRho("B_arXiv_1407_2682_1");

  }else if(Flag == 2){

    myBlue->FixInp();
    myBlue->PrintCompatEst("B_arXiv_1407_2682_2");
    myBlue->Solve();
    printf("... B_arXiv_1407_2682: mtop");
    printf(" (4-observables all-had, l+j, di-l and met) \n");
    printf("... B_arXiv_1407_2682: for the results see Table 5 \n");
    for(Int_t n = 0; n<NumObs; n++){
      myBlue->DisplayResult(n,"B_arXiv_1407_2682_2");
    }
    myBlue->PrintResult();
    myBlue->PrintCompatObs();
    myBlue->GetRhoRes(LocRhoRes);
    printf("... B_arXiv_1407_2682: The correlations of the observables,");
    printf(" see Table 5 \n");
    myBlue->PrintMatrix(LocRhoRes,"%+5.2f");

  }else if(Flag == 3){

    myBlue->FixInp();
    myBlue->PrintCompatEst();
    myBlue->Solve();
    myBlue->PrintCompatObs();
    printf("... B_arXiv_1407_2682: mtop (2 correlated observables CDF + D0)\n");
    printf("... B_arXiv_1407_2682: For the results, see Table 6 \n");
    myBlue->PrintResult();
    myBlue->GetRhoRes(LocRhoRes);
    printf("... B_arXiv_1407_2682: The correlation of the observables \n");
    myBlue->PrintMatrix(LocRhoRes,"%+5.2f");

  }else if(Flag == 4){

    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1407_2682: mtop (2 correlated observables");
    printf(" RUN-I and RUN-II) \n");
    printf("... B_arXiv_1407_2682: Not provided in reference \n");
    myBlue->PrintResult();
    myBlue->GetRhoRes(LocRhoRes);
    printf("... B_arXiv_1407_2682: The correlation of the observables \n");
    myBlue->PrintMatrix(LocRhoRes,"%+5.2f");
  }
  // delet objects
  delete myBlue; myBlue = NULL; 
  LocRho->Delete(); LocRho = NULL;
  LocWei->Delete(); LocWei = NULL;
  LocRhoRes->Delete(); LocRhoRes = NULL;
  return;
}
