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

void B_arXiv_1305_3929(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The Tevatron mtop combination (v2 December 4, 2013)
  //     In addition the D0 and CDF Combinations as well as Run-I and Run-II
  //     combinations are quoted
  //  1: Suggest a combination according to importance up to the point where
  //     each added measurement adds less than 1% in precision
  //     [Shows how to use SolveAccImp()]
  //  2: The combination for four observables, namely M(all-had),
  //     M(lepton+jest), M(di-lepton) and M(missing Et)
  //     [Shows how to retrieve results into TMatrices]
  //  3: The combination for two correlated observables: M(CDF), M(DO)
  //  4: The combination for two correlated observables: M(RUN-I), M(RUN-II)
  //  5: The combination for with relative uncertainties (numerical
  //     exercise only)
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 12;
  static const Int_t NumUnc = 15;
  static const Int_t MaxObs =  4;
  Int_t NumObs = 1;

  // The names of estimates, uncertainties and observables
  TString NamEst[NumEst] = {
    "CDF-I lepton+jets", "    CDF-I dilepton", 
    "    CDF-I alljets",     "  D0-I lepton+jets", 
    "    D0-I dilepton",     "CDF-II lepton+jets",
    "  CDF-II dilepton",   "    CDF-II alljets", 
    "     CDF-II track",      " D0-II lepton+jets",
    "  D0(II) dilepton",   "   CDF-II MET+Jets"};  
  TString NamUnc[NumUnc] = {"   Stat", "   iJES", "   aJES", "   bJES", 
			    "   cJES", "   dJES", "   rJES", "   Lept", 
			    "   Sign", "   DTMO", "  UN/MI", "   BGMC", 
			    "   BGDT", "   Meth", "    MHI"};
  TString NamObs[MaxObs] = {"   mtop", "   mtop", "   mtop", "   mtop"};
 
  // A char to fill fill names
  TString ArcNam = "EPJC_74_3004";
 
  // Index for which estimate determines which observable
  Int_t IWhichObs[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // Flag: steers which of the results should be calculated
  if(Flag == 0 || Flag == 1 || Flag == 5){
    printf("... B_arXiv_1305_3929: -----------------------");
    printf("---------------------- \n");
    printf("... B_arXiv_1305_3929: The 2013 Tevatron mtop combination,");
    printf(" Flag = %2i \n", Flag);
    NumObs = 1;
    if(Flag == 1){
      printf("... B_arXiv_1305_3929: Solve according to importance");
    } else if(Flag == 5){
      printf("... B_arXiv_1305_3929: Numerical exercise mtop combination,");
      printf(" using relative syst uncertainties \n");
    }
  }else if(Flag == 2){
    printf("... B_arXiv_1305_3929: -----------------------");
    printf("----------------------------- \n");
    printf("... B_arXiv_1305_3929: The mtop combination for four observables,");
    printf(" Flag = %2i \n", Flag);
    printf("... B_arXiv_1305_3929: In addition, solve according to importance");
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
    IWhichObs[ 6] = 2;
    IWhichObs[ 7] = 0;
    IWhichObs[ 8] = 1;
    IWhichObs[ 9] = 1;
    IWhichObs[10] = 2;
    IWhichObs[11] = 3;
  }else if(Flag == 3){
    printf("... B_arXiv_1305_3929: -----------------------");
    printf("------------------------------------------ \n");
    printf("... B_arXiv_1305_3929: The mtop combination for");
    printf(" two observables M(CDF), M(DO), Flag = %2i \n", Flag);
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
    IWhichObs[ 9] = 1;
    IWhichObs[10] = 1;
    IWhichObs[11] = 0;
  }else if(Flag == 4){
    printf("... B_arXiv_1305_3929: -----------------------");
    printf("------------------------------------------------ \n");
    printf("... B_arXiv_1305_3929: The mtop combination for");
    printf(" two observables M(RUN-I), M(RUN-II), Flag = %2i \n",Flag);
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
  }else if(Flag == 5){
    printf("... B_arXiv_1305_3929: ------------------------------------ \n");
    printf("... B_arXiv_1305_3929: The 2013 Tevatron mtop combination,  \n");
    printf("... B_arXiv_1305_3929: but with relative uncertainties, a   \n");
    printf("... B_arXiv_1305_3929: purely numerical exercise, Flag = %2i \n", 
	   Flag);
  }else{
    printf("... B_arXiv_1305_3929: Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  // Estimates
  // The order of estimates follows Table 2 of the reference
  // What: 'Stat' 'iJES' 'aJES' 'bJES' 'cJES' 'dJES' 'rJES' 'Lept' 'Sign' 'DTMO' 'UN/MI' 'BGMC' 'BGDT' 'Meth'  'MHI'
  // Num:      0      1      2      3      4      5      6      7      8      9      10     11     12     13     14
  // Rho:      0  Cor01  Cor02      1      1  Cor02  Cor06  Cor02      1  Cor06   Cor06  Cor11   Cor12     0  Cor02
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    176.1,   5.1,  0.0,   0.0,   0.6,   2.7,   0.7,   3.4,   0.0,   2.6,   0.0,    0.0,   1.3,   0.0,   0.0,   0.0,
    167.4,  10.3,  0.0,   0.0,   0.8,   2.6,   0.6,   2.7,   0.0,   2.9,   0.0,    0.0,   0.3,   0.0,   0.7,   0.0,
    186.0,  10.0,  0.0,   0.0,   0.6,   3.0,   0.3,   4.0,   0.0,   2.0,   0.0,    0.0,   1.7,   0.0,   0.6,   0.0,
    180.1,   3.6,  0.0,   0.0,   0.7,   2.0,   2.5,   0.0,   0.0,   1.1,   0.0,    1.3,   1.0,   0.0,   0.6,   0.0,
    168.4,  12.3,  0.0,   0.0,   0.7,   2.0,   1.1,   0.0,   0.0,   1.8,   0.0,    1.3,   1.1,   0.0,   1.1,   0.0,
    172.85,  0.52, 0.49,  0.09,  0.16,  0.21,  0.07,  0.48,  0.03,  0.61,  0.0,    0.00,  0.12,  0.16,  0.00,  0.07,
    170.28,  1.95, 0.00,  0.14,  0.33,  2.13,  0.58,  2.01,  0.27,  0.73,  0.0,    0.00,  0.24,  0.14,  0.12,  0.23,
    172.47,  1.43, 0.95,  0.03,  0.15,  0.24,  0.04,  0.38,  0.00,  0.62,  0.0,    0.00,  0.00,  0.56,  0.38,  0.08,
    166.90,  9.00, 0.00,  0.00,  0.00,  0.36,  0.06,  0.24,  0.00,  0.90,  0.0,    0.00,  0.80,  0.20,  2.50,  0.00,
    174.94,  0.83, 0.53,  0.0,   0.07,  0.0,   0.63,  0.00,  0.17,  0.77,  0.36,   0.00,  0.18,  0.23,  0.16,  0.05, 
    174.00,  2.36, 0.55,  0.40,  0.20,  0.0,   0.56,  0.0,   0.35,  0.86,  0.50,   0.0,   0.0,   0.20,  0.51,  0.00,
    173.95,  1.26, 1.05,  0.10,  0.17,  0.18,  0.04,  0.40,  0.00,  0.64,  0.0,    0.00,  0.00,  0.12,  0.31,  0.18
  };

  // Dimension for correlation matrices
  static const Int_t LenCor = NumEst * NumEst;

  // Correlation matrix Cor01, for uncertainties 1
  // C1   C1   C1   D1   D1   C2   C2   C2   C2   D2   D2   C2
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
  };

  // Correlation matrix Cor02, for uncertainties 2, 5, 7, 14
  // C1   C1   C1   D1   D1   C2   C2   C2   C2   D2   D2   C2
  Double_t Cor02[LenCor] = {
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
  };

  // Correlation matrix Cor06, for uncertainties 6, 9, 10
  // C1   C1   C1   D1   D1   C2   C2   C2   C2   D2   D2   C2
  Double_t Cor06[LenCor] = {
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0
  };

  // Correlation matrix Cor11, for uncertainty 11
  //  1    2    0    1    2    1    2    0    1    1    2    3
  Double_t Cor11[LenCor] = {
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
  };
  
  // Correlation matrix Cor12, for uncertainty 12
  //  1    2    0    1    2    1    2    0    1    1    2    3
  Double_t Cor12[LenCor] = {
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
  };
  
  //-- Local Structures for BLUE output
  // TMatrices
  TMatrixD* LocCov  = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocCovI = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRho  = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocWei  = new TMatrixD(NumEst,NumObs);
  TMatrixD* LocInsLik = new TMatrixD(NumObs,7);
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
    }else if(k == 2 || k == 5 || k ==  7|| k == 14){
      myBlue->FillCor(k, &Cor02[0]);
    }else if(k == 6 || k == 9 || k == 10){
      myBlue->FillCor(k, &Cor06[0]);
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
    myBlue->PrintStatus();
    myBlue->SolveInfWei();
    myBlue->PrintInfWei();
    printf("... B_arXiv_1305_3929: mtop full combination, for correlations");
    printf(" see Table 2 \n");

    myBlue->PrintRho();
    printf("... B_arXiv_1305_3929: mtop full combination \n");
    printf("... B_arXiv_1305_3929: for weight matrix and results");
    printf(" see Table 4+3 \n");
    myBlue->PrintResult();
    myBlue->PrintPull();

    myBlue->LatexResult("B_arXiv_1305_3929");
    myBlue->DisplayResult(0,"B_arXiv_1305_3929");

    // Independent CDF Combination
    myBlue->ReleaseInp();
    myBlue->SetInActiveEst(3);
    myBlue->SetInActiveEst(4);
    myBlue->SetInActiveEst(9);
    myBlue->SetInActiveEst(10);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1305_3929: mtop CDF combination. The text after");
    printf(" Table 5 has the correlated combination, see Flag == 3 \n");
    myBlue->PrintResult();

    // Independent D0 Combination
    myBlue->ResetInp();
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(1);
    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(5);
    myBlue->SetInActiveEst(6);
    myBlue->SetInActiveEst(7);
    myBlue->SetInActiveEst(8);
    myBlue->SetInActiveEst(11);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1305_3929: mtop  D0 combination. The text after");
    printf(" Table 5 has the correlated combination, see Flag == 3 \n");
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
    printf("... B_arXiv_1305_3929: mtop RUN-I combination, not in reference \n");
    myBlue->PrintResult();

    // Independent RUN-II Combination
    myBlue->ResetInp();
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(1);
    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(3);
    myBlue->SetInActiveEst(4);
    myBlue->FixInp();
    myBlue->PrintCompatEst();
    myBlue->Solve();
    printf("... B_arXiv_1305_3929: mtop RUN-II combination, not in reference \n");
    myBlue->PrintResult();
    myBlue->PrintCompatObs();

  }else if(Flag == 1){

    myBlue->FixInp();
    myBlue->PrintCompatEst("B_arXiv_1305_3929_1_All");
    myBlue->SolveAccImp(0,1.0);
    printf("... B_arXiv_1305_3929: The results of the successive combination");
    printf(" Figure 6 of %s \n",ArcNam.Data());
    myBlue->PrintAccImp();
    myBlue->DisplayAccImp(0,"B_arXiv_1305_3929_1");

    // Get the indices from SolveAccImp()
    Int_t ImpInd[NumEst] = {0};
    myBlue->GetAccImpIndEst(0,ImpInd);
    Int_t ImpLas = myBlue->GetAccImpLasEst(0);
    Int_t ImpFir = ImpInd[0];
    Int_t ImpNum = 0;
    for(Int_t i = 0; i<NumEst; i++)if(ImpInd[i] == ImpLas)ImpNum = i;
    printf("... B_arXiv_1305_3929: The list of importance is:");
    for(Int_t i = 0; i<NumEst; i++)printf(" %3i", ImpInd[i]);
    printf("\n");
    printf("... B_arXiv_1305_3929: The first and the last of the");
    printf(" %3i estimates are %3i and %3i \n", ImpNum+1, ImpFir, ImpLas);

    // Print the parameters
    printf("... B_arXiv_1305_3929: The parameters of the pairwise");
    printf("  combination, Table 3 of %s \n",ArcNam.Data());
    printf("... B_arXiv_1305_3929: The most precise estimate is %s \n",
	   NamEst[ImpFir].Data());
    myBlue->PrintParams();

    // Get the weights
    printf("... B_arXiv_1305_3929: The Blue weights \n");
    myBlue->GetWeight(LocWei);
    LocWei->operator*=(100.);
    myBlue->PrintMatrix(LocWei," %+4.1f");

    // Look at pairs wrt the most precise according to the list from
    // SolveAccImp(). Produce sub-figures for Figure 7 from EPJC_74_3004
    Int_t IndFig = 0;
    printf("... B_arXiv_1305_3929: Inspect all pairs with");
    printf(" the most precise estimate \n");
    for(Int_t i = 1; i<NumEst; i++){
      IndFig = 0;
      if(i <= ImpNum && LocWei->operator()(ImpInd[i],0) < 0)IndFig = 1;
      myBlue->InspectPair(ImpFir, ImpInd[i],"B_arXiv_1305_3929",IndFig);
    }

    // Disable all, include one at a time according to SolveAccImp, Print Result
    myBlue->ReleaseInp();
    for(Int_t i = 1; i<NumEst; i++)myBlue->SetInActiveEst(ImpInd[i]);
    printf("... B_arXiv_1305_3929: Solve successively according to");
    printf(" importance for all estimates \n");
    for(Int_t i = 1; i<NumEst; i++){
      myBlue->SetActiveEst(ImpInd[i]);
      myBlue->FixInp();
      myBlue->Solve();
      myBlue->PrintResult();
      myBlue->ReleaseInp();
    }

    // Disable all insignificant estimates, solve for the selected estimates
    myBlue->ReleaseInp();
    for(Int_t i = ImpNum+1; i<NumEst; i++)myBlue->SetInActiveEst(ImpInd[i]);
    myBlue->FixInp();
    myBlue->PrintCompatEst("B_arXiv_1305_3929_1_Sel");
    myBlue->Solve();
    printf("... B_arXiv_1305_3929: Result based on the selected estimates \n");
    myBlue->PrintResult();
    myBlue->ReleaseInp();

    // Scan the rho parameters for the selected estimates
    printf("... B_arXiv_1305_3929: The individual scans of the correlations");
    printf(" of the selected estimates\n");
    for(Int_t l=0; l<2; l++){
      myBlue->FixInp();
      myBlue->SolveScaRho(l);
      myBlue->PrintScaRho("B_arXiv_1305_3929_1");
      myBlue->ReleaseInp();
    } 

  }else if(Flag == 2){

    myBlue->FixInp();
    myBlue->PrintCompatEst("B_arXiv_1305_3929_2");
    myBlue->PrintListEst();
    myBlue->PrintListUnc();
    myBlue->PrintListObs();
    for(Int_t n = 0; n<NumObs; n++){
      printf("... B_arXiv_1305_3929: The most prescise estimate of");
      printf(" observable %2i is %2i \n", n, myBlue->GetPreEst(n));
    }
    myBlue->Solve();
    printf("... B_arXiv_1305_3929: mtop (4-observables all-had, l+j,");
    printf(" di-l and  met) \n");
    printf("... B_arXiv_1305_3929: Results and Correlations, see Table 5 \n");
    myBlue->PrintResult();
    myBlue->PrintCompatObs();
    myBlue->PrintRhoRes();

    //-- Examples of how to extract results into local structures
    // Covariance matrix of estimates
    printf("... B_arXiv_1305_3929: Get the Covariance matrix as TMatrixD \n");
    myBlue->GetCov(LocCov);
    myBlue->PrintMatrix(LocCov,"%+7.2f");

    // Inverse covariance matrix of estimates
    printf("... B_arXiv_1305_3929: The inverse covariance matrix \n");
    myBlue->GetCovInvert(LocCovI);
    myBlue->PrintMatrix(LocCovI,"%+7.4f");

    // Correlation of estimates
    printf("... B_arXiv_1305_3929: The correlation matrix of the estimates \n");
    myBlue->GetRho(LocRho);
    myBlue->PrintMatrix(LocRho);
    //-- End

    // Also solve according to importance
    myBlue->ReleaseInp();
    myBlue->FixInp();
    myBlue->SolveAccImp(0,1.0);
    myBlue->PrintAccImp();
    Int_t IndEst[NumEst] = {0};
    for(Int_t n = 0; n<NumObs; n++){
      printf("... B_arXiv_1305_3929: Number of estimates for observable");
      printf(" = %2i is %2i \n", n,myBlue->GetActEst(n));
      myBlue->GetAccImpIndEst(n,IndEst);

      printf("... B_arXiv_1305_3929: List of estimates for observable");
      printf(" = %2i is: %2i",n, IndEst[0]);
      for(Int_t l = 1; l<myBlue->GetActEst(n); l++)printf(", %2i",IndEst[l]);
      printf("\n");

      printf("... B_arXiv_1305_3929: The last to be used estimates for");
      printf(" observable = %2i is %2i \n", n,myBlue->GetAccImpLasEst(n));

      myBlue->DisplayResult(n,"B_arXiv_1305_3929_2");
      myBlue->DisplayAccImp(n,"B_arXiv_1305_3929_2");
    }

  }else if(Flag == 3){

    myBlue->FixInp();
    myBlue->PrintCompatEst();
    myBlue->Solve();
    myBlue->PrintCompatObs();
    printf("... B_arXiv_1305_3929: mtop (2 correlated observables CDF + D0) \n");
    printf("... B_arXiv_1305_3929: For results, see text after Table 5 \n");
    myBlue->PrintResult();
    myBlue->PrintRhoRes();

  }else if(Flag == 4){

    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1305_3929: mtop (2 correlated observables RUN-I");
    printf(" and RUN-II) \n");
    printf("... B_arXiv_1305_3929: Not provided in reference \n");
    myBlue->PrintResult();
    myBlue->PrintRhoRes();

  }else if(Flag == 5){

    myBlue->SetRelUnc();
    myBlue->SetNotRelUnc(0);
    myBlue->FixInp();
    //myBlue->PrintStatus();
    myBlue->SolveRelUnc(1.0);
    printf("... B_arXiv_1305_3929: Numerical exercise mtop combination,");
    printf(" using relative syst uncertainties \n");
    myBlue->PrintResult();
    myBlue->InspectLike(0,"B_arXiv_1305_3929");
    myBlue->GetInspectLike(LocInsLik);
    myBlue->PrintMatrix(LocInsLik,myBlue->GetActObs(),7);

  }
  // delet objects
  delete myBlue; myBlue = NULL; 
  LocCov->Delete(); LocCov = NULL;
  LocCovI->Delete(); LocCovI = NULL;
  LocRho->Delete(); LocRho = NULL;
  LocWei->Delete(); LocWei = NULL;
  LocInsLik->Delete(); LocInsLik = NULL;
  return;
}
