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

void B_arXiv_1107_5255(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The Tevatron mtop combination (v3 September 8, 2011)
  //     In addition the D0 and CDF Combinations as well as RunI and RunII
  //     combinations are quoted
  //  1: Suggest a combination by importance up to the point were each added
  //     measurement improves the result in precision by less than 1%.
  //     [Shows how to use SolveAccImp()]
  //  2: The Tevatron combination for four observables, namely M(all-had),
  //     M(lepton+jest), M(di-lepton) and M(missing Et)
  //     [Shows how to retrieve results into TMatrices]
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 12;
  static const Int_t NumUnc = 15;
  Int_t NumObs = 1;

  // Index for which estimate determines which observable
  Int_t IWhichObs[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // Flag: steers which of the results should be calculated
  if(Flag == 0 || Flag == 1){
    NumObs = 1;
  }else if(Flag == 2){
    // 0 = allh, 1=l+j, 2=di-l, 3=Met
    NumObs = 4;
    IWhichObs[ 0] = 1;
    IWhichObs[ 1] = 2;
    IWhichObs[ 2] = 0;
    IWhichObs[ 3] = 1;
    IWhichObs[ 4] = 2;
    IWhichObs[ 5] = 1;
    IWhichObs[ 6] = 2;
    IWhichObs[ 7] = 1;
    IWhichObs[ 8] = 1;
    IWhichObs[ 9] = 2;
    IWhichObs[10] = 0;
    IWhichObs[11] = 3;
  }else{
    printf("... B_arXiv_1107_5255: Not implemente Flag = %2i \n",Flag);
    return;
  }
  
  // Estimates
  // The order of estimates follows Table 2 of the reference
  // What: 'Stat' 'iJES' 'aJES' 'bJES' 'cJES' 'dJES' 'rJES' 'Lept' 'Sign' 'DTMO' 'UN/MI' 'BGMC' 'BGDT' 'Meth'  'MHI'
  // Num:      0      1      2      3      4      5      6      7      8      9      10     11     12     13     14
  // Rho:      0      0  Cor02      1      1  Cor02  Cor06  Cor02      1  Cor06   Cor06  Cor11      0      0  Cor02
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    176.1,   5.1,  0.0,   0.0,   0.6,   2.7,   0.7,   3.4,   0.0,   2.6,   0.0,    0.0,   1.3,   0.0,   0.0,   0.0,
    167.4,  10.3,  0.0,   0.0,   0.8,   2.6,   0.6,   2.7,   0.0,   2.9,   0.0,    0.0,   0.3,   0.0,   0.7,   0.0,
    186.0,  10.0,  0.0,   0.0,   0.6,   3.0,   0.3,   4.0,   0.0,   2.0,   0.0,    0.0,   1.7,   0.0,   0.6,   0.0,
    180.1,   3.6,  0.0,   0.0,   0.7,   2.0,   0.0,   2.5,   0.0,   1.1,   0.0,    1.3,   1.0,   0.0,   0.6,   0.0,
    168.4,  12.3,  0.0,   0.0,   0.7,   2.0,   0.0,   1.1,   0.0,   1.8,   0.0,    1.3,   1.1,   0.0,   1.1,   0.0,
    173.00,  0.65, 0.58,  0.13,  0.23,  0.27,  0.01,  0.41,  0.14,  0.56,  0.0,    0.00,  0.27,  0.06,  0.10,  0.10,
    170.28,  1.95, 0.00,  0.14,  0.33,  2.13,  0.58,  2.01,  0.27,  0.73,  0.0,    0.00,  0.24,  0.14,  0.12,  0.23,
    166.90,  9.00, 0.00,  0.00,  0.00,  0.36,  0.06,  0.24,  0.00,  0.90,  0.0,    0.00,  0.80,  0.20,  2.50,  0.00,
    174.94,  0.83, 0.53,  0.0,   0.07,  0.0,   0.63,  0.00,  0.18,  0.77,  0.36,   0.00,  0.18,  0.23,  0.16,  0.05, 
    173.97,  1.83, 0.0,   1.57,  0.40,  0.0,   1.50,  0.0,   0.49,  0.74,  0.33,   0.0,   0.0,   0.47,  0.10,  0.00,
    172.47,  1.43, 0.95,  0.03,  0.15,  0.24,  0.04,  0.38,  0.00,  0.62,  0.0,    0.00,  0.00,  0.56,  0.38,  0.08,
    172.32,  1.80, 1.54,  0.12,  0.26,  0.20,  0.05,  0.45,  0.00,  0.74,  0.0,    0.00,  0.00,  0.12,  0.14,  0.16
  };

  // Dimension for correlation matrices
  static const Int_t LenCor = NumEst * NumEst;

  // Correlation matrix Cor02, for uncertainties 2, 5, 7, 14
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

  // Correlation matrix Cor06, for uncertainties 6, 9, 10
  // C1   C1   C1   D1   D1   C2   C2   C2   D2   D2   C2   C2
  Double_t Cor06[LenCor] = {
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
  //  1    2    0    1    2    1    2    1    1    2    0    3
  Double_t Cor11[LenCor] = {
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
  };
  
  //-- Local Structures for BLUE output
  // TMatrices
  Int_t iok = 0;
  TMatrixD* LocCov = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocCovI = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);

  // Double_t Arrays
  Double_t LocCovArr[LenCor] = {LenCor*0.0};
  Double_t LocCovIArr[LenCor]= {LenCor*0.0};
  Double_t LocRhoArr[LenCor] = {LenCor*0.0};  
  //-- End

  // Define formats for Figures and Latex file
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

  // Fill all estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }

  // Fill all correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(k ==  0 || k ==  1|| k == 12 || k == 13){
      myBlue->FillCor(k,0.0);
    }else if(k == 2 || k == 5 || k ==  7|| k == 14){
      myBlue->FillCor(k, &Cor02[0]);
    }else if(k == 6 || k == 9 || k == 10){
      myBlue->FillCor(k, &Cor06[0]);
    }else if(k == 11){
      myBlue->FillCor(k, &Cor11[0]);
    }else{
      myBlue->FillCor(k,1.0);
    }
  }

  // Fix input, solve several times, and finally delete
  if(Flag == 0){
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1107_5255: mtop full combination,");
    printf(" for the correlation see Table 2 \n");
    myBlue->PrintRho();
    printf("... B_arXiv_1107_5255: mtop full combination, for the weight matrix");
    printf(" and result see Table 4+3 \n");
    myBlue->PrintResult();

    //  CDF Combination
    myBlue->ReleaseInp();
    myBlue->SetInActiveEst(3);
    myBlue->SetInActiveEst(4);
    myBlue->SetInActiveEst(8);
    myBlue->SetInActiveEst(9);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1107_5255: mtop CDF combination,");
    printf(" see text after Table 5 \n");
    myBlue->PrintResult();

    //   D0 Combination
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
    printf("... B_arXiv_1107_5255: mtop D0 combination,");
    printf(" see text after Table 5 \n");
    myBlue->PrintResult();

    // RUN-I Combination
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
    printf("... B_arXiv_1107_5255: mtop RUN-I combination,");
    printf(" see text after Table 5 \n");
    myBlue->PrintResult();

    // RUN-II Combination
    myBlue->ResetInp();
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(1);
    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(3);
    myBlue->SetInActiveEst(4);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1107_5255: mtop RUN-II combination,");
    printf(" see text after Table 5 \n");
    myBlue->PrintResult();
  }else if(Flag == 1){
    myBlue->FixInp();
    myBlue->SolveAccImp(1.0);
    myBlue->PrintAccImp();
  }else if(Flag == 2){
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1107_5255: mtop");
    printf(" (4-observables all-had, l+j, di-l, met)\n");
    printf("... B_arXiv_1107_5255: Results and correlations, see Table 5 \n");
    myBlue->PrintResult();
    myBlue->PrintRhoRes();
    myBlue->PrintCompatObs();

    // Examples of how to use the display feature
    for(Int_t n = 0; n<NumObs; n++){
      myBlue->DisplayResult(n,"B_arXiv_1107_5255");
    }
    myBlue->LatexResult("B_arXiv_1107_5255");

    // Examples of how to extract results into local structures
    // Covariance matrix of estimates
    printf("... B_arXiv_1107_5255: The covariance matrix \n");
    iok = myBlue->GetCov(LocCov);
    myBlue->PrintMatrix(LocCov,"%6.2f");
    printf("... B_arXiv_1107_5255: The same as array of doubles \n");
    iok = myBlue->GetCov(LocCovArr);
    if(iok == 1)myBlue->PrintDouble(LocCovArr,myBlue->GetActEst(),
				    myBlue->GetActEst(),"%6.2f");

    // Inverse covariance matrix of estimates
    printf("... B_arXiv_1107_5255: The inverted covariance matrix \n");
    iok = myBlue->GetCovInvert(LocCovI);
    if(iok == 1)myBlue->PrintMatrix(LocCovI," %7.3f");
    printf("... B_arXiv_1107_5255: The same as array of doubles \n");
    iok = myBlue->GetCovInvert(LocCovIArr);
    if(iok == 1)myBlue->PrintDouble(LocCovIArr,myBlue->GetActEst(),
				    myBlue->GetActEst(),"%7.3f");

    // Correlation of estimates
    printf("... B_arXiv_1107_5255: The correlations of the estimates \n");
    myBlue->GetRho(LocRho);
    LocRho->operator*=(100.);
    myBlue->PrintMatrix(LocRho,"%6.1f%%");

    printf("... B_arXiv_1107_5255: The same as array of doubles \n");
    iok = myBlue->GetRho(LocRhoArr);
    if(iok == 1)myBlue->PrintDouble(LocRhoArr,myBlue->GetActEst(),
				    myBlue->GetActEst()," %6.3f");
  }
  // Delete objects
  delete myBlue; myBlue = NULL; 
  LocCov->Delete(); LocCov = NULL;
  LocCovI->Delete(); LocCovI = NULL;
  LocRho->Delete(); LocRho = NULL;
  return;
}
