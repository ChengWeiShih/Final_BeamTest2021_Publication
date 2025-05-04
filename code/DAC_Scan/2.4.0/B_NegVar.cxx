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

void B_NegVar(Int_t Flag = 0){

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 3;
  static const Int_t NumUnc = 4;
  static const Int_t NumObs = 1;

  // The names of estimates, uncertainties and observables
  TString NamEst[NumEst] = {"     x0", "     x1", "     x2"};
  TString NamUnc[NumUnc] = {"   Stat", "  Sys_1", "  Sys_2", "  Sys_3"};
  TString NamObs[NumObs] = {"      x"};

  // Fill three estimates with four uncertainties
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //       k=0     1     2     3
    172.50, 0.06, 0.15, 0.69, 0.00,
    172.35, 0.05, 0.12, 0.10, 0.04,
    172.80, 0.09, 0.14, 0.28, 0.12
  };

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_NegVar: --------------------------------------------");
    printf("------------------------------------\n");
    printf("... B_NegVar: Get a -negative- uncertainty for a specific");
    printf(" source from -negative-      weights, \n");
    printf("... B_NegVar: also resulting in a -negative- variance");
    printf(" for the observable, Flag = %2i \n", Flag);
  }else if(Flag == 1){
  // Clean up and return
    printf("... B_NegVar: --------------------------------------------");
    printf("----------------------------------------------\n");
    printf("... B_NegVar: Get a -negative- uncertainty for a specific");
    printf(" source from -negative- correlations, Flag = %2i \n", Flag);
  }else{
    printf("... B_NegVar: Not implemented Flag = %2i \n",Flag);
    return;
  }
  printf("... B_NegVar:\n");

  // Set the correlation for k=2,3
  Double_t Rho = 1.0 - Flag * 2.0;

  // Fill the correlation matrix to be used for k=2,3
  static const Int_t LenCor = NumEst * NumEst;
  Double_t RhoMat[LenCor] = {
    1.0, Rho, 0.0,
    Rho, 1.0, Rho,
    0.0, Rho, 1.0
  };

  // Local Structures for BLUE output
  TMatrixD* LocCov = new TMatrixD(NumEst,NumEst);
  TMatrixD* TotCov = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocWei = new TMatrixD(NumEst,NumObs);
  // End

  // Matrices to calculate the uncertainty of the combined 
  // result for source Sys_2
  TMatrixD* LocWeiT = new TMatrixD(NumObs,NumEst);
  TMatrixD* H = new TMatrixD(NumEst,NumObs);
  TMatrixD* F = new TMatrixD(NumObs,NumObs);
  // End

  // Define file name
  const TString FilBas = "B_NegVar";
  TString FilNam = "To be filled later";
  char Buffer[100];
  sprintf(Buffer,"%1i", Flag);
  FilNam = &Buffer[0];
  FilNam = FilBas + "_" + FilNam;

  // Define formats
  const TString ForVal = "%6.3f";
  const TString ForUnc = ForVal;
  const TString ForWei = ForUnc;
  const TString ForRho = ForUnc;
  const TString ForPul = ForUnc;
  const TString ForUni = "None";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc);
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);
  myBlue->PrintStatus();
  myBlue->PrintStatus();

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
  
  // Fill correlation
  for(Int_t k = 0; k<NumUnc; k++){
    if(k == 0){myBlue->FillCor(k, 0.0);
    }else if(k == 1){myBlue->FillCor(k, 1.0);
    }else{myBlue->FillCor(k, &RhoMat[0]);
    }
  }
  
  // Get covariance matrix of uncertainty Sys_2 == Mij[k=2] for later use,
  // by disabling k = 0, 1, 3 and retrieving the total covariance matrix
  myBlue->SetInActiveUnc(0);
  myBlue->SetInActiveUnc(1);
  myBlue->SetInActiveUnc(3);
  myBlue->FixInp(); 
  myBlue->GetCov(LocCov);
  
  // Reset and fix
  myBlue->ResetInp(); 
  myBlue->FixInp();

  // The estimates
  printf("... B_NegVar:\n");
  printf("... B_NegVar: The estimates \n");
  printf("... B_NegVar:\n");
  myBlue->PrintEst();
  
  // Their correlation
  printf("... B_NegVar:\n");
  printf("... B_NegVar: The correlations of the estimates \n");
  printf("... B_NegVar:\n");
  myBlue->GetRho(LocRho);
  myBlue->PrintMatrix(LocRho,ForRho);
  
  // The total covariance matrix
  printf("... B_NegVar:\n");
  printf("... B_NegVar: The total covariance matrix \n");
  printf("... B_NegVar:\n");
  myBlue->GetCov(TotCov);
  myBlue->PrintMatrix(TotCov,ForUnc);
  
  // Solve it  
  myBlue->Solve();
  
  // The result    
  printf("... B_NegVar:\n");
  printf("... B_NegVar: The result of the combination \n");
  printf("... B_NegVar:\n");
  myBlue->PrintResult();

  myBlue->SetPrintLevel(1);
  printf("... B_NegVar: InspectResult = %i \n", myBlue->InspectResult());
  
  // Write out the result
  myBlue->LatexResult(FilNam);
  
  // Write what we intend to do
  printf("... B_NegVar:\n");
  printf("... B_NegVar: Calculate the -variance- for the source Sys_2 \n");
  printf("... B_NegVar: using Eq.18 of NIMA500(2003)391 \n");
  printf("... B_NegVar: Lambda_i^T * Mij[k=2] * Lambda_j =");
  printf(" LocWeiT * LocCov * LocWei\n");
  
  // Print LocCov
  printf("... B_NegVar: The contribution to the covariance matrix");
  printf(" from source %s. This equals Mij[k=2] \n",NamUnc[2].Data());
  printf("... B_NegVar:\n");
  myBlue->PrintMatrix(LocCov,ForUnc);
      
  // Get the Weight matrix and its transpose
  printf("... B_NegVar:\n");
  printf("... B_NegVar: The matrix of weights \n");
  printf("... B_NegVar:\n");
  myBlue->GetWeight(LocWei);
  myBlue->PrintMatrix(LocWei,ForWei);
  LocWeiT->Transpose(*LocWei);
      
  // Calculate the "variance" for Sys_2 (k=2)
  H->Mult(*LocCov, *LocWei);
  F->Mult(*LocWeiT,*H);
      
  // Print an intermediate result
  printf("... B_NegVar:\n");
  printf("... B_NegVar: Mij[k=2] * Lambda_j =  \n");
  printf("... B_NegVar:\n");
  myBlue->PrintMatrix(H,ForUnc);
  
  // Print the "variance" and the uncertainty
  printf("... B_NegVar:\n");
  printf("... B_NegVar: The -variance- for k=2 \n");
  printf("... B_NegVar:\n");
  myBlue->PrintMatrix(F,ForUnc);
  Double_t Sys_2 = -1.0*TMath::Sqrt(-1.0*F->operator()(0,0));
  printf("... B_NegVar: Thus the -negative uncertainty- for the source");
  printf(" %s evaluates to %5.2f \n",NamUnc[2].Data(), Sys_2);
  printf("... B_NegVar:\n");
  if(Flag == 0){
    printf("... B_NegVar: Since the absolute value of this negative" );
    printf(" variance is larger than \n" );
    printf("... B_NegVar: the statistical uncertainty,");
    printf(" the total variance is negative. \n");
    printf("... B_NegVar: This constitutes a non-solvable problem for");
    printf(" the BLUE combination. \n");
    printf("... B_NegVar: It is caused by: \n");
    printf("... B_NegVar: 1) the partially large estimator correlations \n");
    printf("... B_NegVar: 2) the imbalanced size of this systematic");
    printf(" uncertainty across the estimators \n");
    printf("... B_NegVar: and \n");
    printf("... B_NegVar: 3) the large ratio of systematic");
    printf(" and statistical uncertainties. \n");
  }else{
    printf("... B_NegVar: Since the absolute value of this negative" );
    printf(" variance is smaller than \n" );
    printf("... B_NegVar: the statistical uncertainty,");
    printf(" the total variance is positive. \n");
    printf("... B_NegVar: This still constitutes a solvable problem for");
    printf(" the BLUE combination. \n");
    printf("... B_NegVar: However, the interpretation of source k=2 as an");
    printf(" uncertainty in the \n");
    printf("... B_NegVar: combined value from this");
    printf(" systematic variation is not possible.\n");
  }
  printf("... B_NegVar:\n");

  // Delete object
  delete myBlue;
  myBlue = NULL;

  // Delete matrices
  delete LocCov; LocCov = NULL;
  delete TotCov; TotCov = NULL;
  delete LocRho; LocRho = NULL;
  delete LocWei; LocWei = NULL;
  delete LocWeiT; LocWeiT = NULL;
  delete H; H = NULL;
  delete F; F = NULL;
  
  // Return
  return;
}
