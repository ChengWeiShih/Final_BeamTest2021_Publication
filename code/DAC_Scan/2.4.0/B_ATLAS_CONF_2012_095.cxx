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

void B_ATLAS_CONF_2012_095(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The   LHC Combination
  //     [Shows how to retrieve quantities with GetCovRes(), GetRhoRes()]
  //     [Shows how to retrieve quantities with GetWeight(), GetResult()]
  //  1: The ATLAS Combination
  //  2: The   CMS Combination
  //  3: The Combination for 3 Observables = l+jets, di-lepton and all-jets
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  7;
  static const Int_t NumUnc = 16;
  Int_t NumObs =  1;
  static const Int_t LocObs = 1;

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {7*0};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_ATLAS_CONF_2012_095: The 2011   LHC Combination = %2i\n",
	   Flag);
  }else if(Flag == 1){
    printf("... B_ATLAS_CONF_2012_095: The 2011 ATLAS Combination = %2i\n",
	   Flag);
  } else if(Flag == 2){
    printf("... B_ATLAS_CONF_2012_095: The 2011   CMS Combination = %2i\n",
	   Flag);
  } else if(Flag == 3){
    printf("... B_ATLAS_CONF_2012_095: The combination for l+jets,");
    printf(" di-lepton and all-jets = %2i \n",Flag);
    NumObs =  3;
    IWhichObs[0] = 0;
    IWhichObs[1] = 0;
    IWhichObs[2] = 2;
    IWhichObs[3] = 1;
    IWhichObs[4] = 0;
    IWhichObs[5] = 1;
    IWhichObs[6] = 0;
  }else{
    printf("... ATLAS_CONF_2012_095: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Estimates 0-6 
  // 0 == ATLAS_2010 lepton+jets
  // 1 == ATLAS_2011 lepton+jets
  // 2 == ATLAS_2011    all jets
  // 3 ==   CMS_2010   di-lepton
  // 4 ==   CMS_2010 lepton+jets
  // 5 ==   CMS_2011   di-lepton
  // 6 ==   CMS_2011   muon+jets

  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //       0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
    169.3, 4.0, 0.0, 2.5, 2.1, 0.0, 0.0, 1.0, 2.5, 0.6, 0.5, 1.2, 0.6, 1.8, 0.6, 0.4, 0.7,
    174.5, 0.6, 0.4, 1.6, 0.7, 0.0, 0.0, 0.4, 1.0, 0.6, 0.1, 0.3, 0.6, 0.1, 0.5, 0.1, 0.0,
    174.9, 2.1, 0.0, 1.4, 2.1, 0.0, 0.0, 0.5, 1.7, 0.6, 0.6, 0.5, 0.6, 0.0, 1.9, 1.0, 0.0,
    175.5, 4.6, 0.0, 0.9, 2.1, 3.3, 0.3, 0.4, 0.9, 0.5, 0.5, 0.7, 1.4, 0.1, 0.0, 0.3, 1.0,
    173.1, 2.1, 0.0, 0.9, 2.1, 0.0, 0.0, 0.0, 1.2, 0.5, 0.1, 0.4, 0.2, 0.2, 0.4, 0.1, 0.1,
    173.3, 1.2, 0.0, 1.1, 2.0, 0.0, 0.2, 0.1, 0.8, 0.5, 0.4, 0.7, 0.6, 0.0, 0.4, 0.4, 0.2,
    172.6, 0.4, 0.4, 0.7, 0.2, 0.0, 0.0, 0.0, 0.8, 0.5, 0.1, 0.3, 0.6, 0.1, 0.0, 0.2, 0.4,
  };
 
  static const Int_t LenCor = NumEst * NumEst;
  //  0,  1,  4, 13, 14 <==}  rho_exp=0, rho_LHC=  0 <==> rho == 0 
  //          2,  6,  7 <==}  rho_exp=1, rho_LHC=0.5 <==> Cor02 
  //      3,  5, 10, 11 <==}  rho_exp=1, rho_LHC=  0 <==> Cor03
  //      8,  9, 12, 15 <==}  rho_exp=1, rho_LHC=  1 <==> rho == 1

  Double_t Cor02[LenCor] = {
    1.00, 1.00, 1.00, 0.50, 0.50, 0.50, 0.50,
    1.00, 1.00, 1.00, 0.50, 0.50, 0.50, 0.50,
    1.00, 1.00, 1.00, 0.50, 0.50, 0.50, 0.50,
    0.50, 0.50, 0.50, 1.00, 1.00, 1.00, 1.00,
    0.50, 0.50, 0.50, 1.00, 1.00, 1.00, 1.00,
    0.50, 0.50, 0.50, 1.00, 1.00, 1.00, 1.00,
    0.50, 0.50, 0.50, 1.00, 1.00, 1.00, 1.00,
  };

  Double_t Cor03[LenCor] = {
    1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00,
    1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00,
    1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00,
    0.00, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00,
    0.00, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00,
    0.00, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00,
  };

  //-- Local Structures for Blue output
  // TMatrices
  Int_t iok = 0;
  TMatrixD* LocCovRes = new TMatrixD(NumObs,NumObs);
  TMatrixD* LocRhoRes = new TMatrixD(NumObs,NumObs);
  TMatrixD* LocWeight = new TMatrixD(NumEst,NumObs);
  TMatrixD* LocResult = new TMatrixD(NumObs,NumUnc+1);

  // Double_t Arrays
  static const Int_t LenRes = LocObs * LocObs;
  Double_t LocCovResArr[LenRes] = {LenRes*0};
  Double_t LocRhoResArr[LenRes] = {LenRes*0.0};
  
  static const Int_t LenWei = NumEst * LocObs;
  Double_t LocWeightArr[LenWei] = {LenWei*0.0};
  
  static const Int_t LenResult = LocObs * NumUnc+1;
  Double_t LocResultArr[LenResult] = {LenResult*0.0};
  //-- End


  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%5.3f";
  const TString ForRho = ForVal;
  const TString ForPul = ForVal;
  const TString ForUni = "GeV";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);

  // Fill estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }

  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(k == 0 || k == 1 || k == 4 || k == 13 || k == 14){
      myBlue->FillCor(k,0.0);
    }else if(k == 2 || k == 6 || k == 7){
      myBlue->FillCor(k,&Cor02[0]);
    }else if(k == 3 || k == 5 || k == 10 || k == 11){
      myBlue->FillCor(k,&Cor03[0]);
    }else{
      myBlue->FillCor(k,1.0);
    }
  }

  if(Flag == 0){
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2012_095: The 2011   LHC Combination = %2i \n",
	   Flag);
    printf("... B_ATLAS_CONF_2012_095: For the results see Table 4 \n");
    myBlue->PrintResult();
    printf("... B_ATLAS_CONF_2012_095: For the correlations");
    printf(" see Matrix on page 9\n");
    myBlue->PrintRho();

    // Examples of how to extract results into local structures
    // Covariance matrix of observables
    printf("... B_ATLAS_CONF_2012_095: Covariance matrix of observables \n");
    iok = myBlue->GetCovRes(LocCovRes);
    if(iok == 1)myBlue->PrintMatrix(LocCovRes);
    iok = myBlue->GetCovRes(LocCovResArr);
    if(iok == 1)myBlue->PrintDouble(LocCovResArr,myBlue->GetActObs(),
				    myBlue->GetActObs());
    
    // Correlation matrix of observables
    printf("... B_ATLAS_CONF_2012_095: Correlation matrix of observables \n");
    myBlue->PrintRhoRes();
    iok = myBlue->GetRhoRes(LocRhoRes);
    if(iok == 1)LocRhoRes->Print();
    iok = myBlue->GetRhoRes(LocRhoResArr);
    if(iok == 1)myBlue->PrintDouble(LocRhoResArr,myBlue->GetActObs(),
				    myBlue->GetActObs());

    // Weight matrix of estimates for all observables
    printf("... B_ATLAS_CONF_2012_095: Weight matrix \n");
    printf("... B_ATLAS_CONF_2012_095: For the weights see Fig 1b \n");
    iok = myBlue->GetWeight(LocWeight);
    if(iok == 1)myBlue->PrintMatrix(LocWeight," %+5.2f");
    iok = myBlue->GetWeight(LocWeightArr);
    if(iok == 1)myBlue->PrintDouble(LocWeightArr,myBlue->GetActEst(),
				    myBlue->GetActObs()," %+5.2f");
    
    // The result
    printf("... B_ATLAS_CONF_2012_095: Result \n");
    iok = myBlue->GetResult(LocResult);
    if(iok == 1)myBlue->PrintMatrix(LocResult);
    iok = myBlue->GetResult(LocResultArr);
    if(iok == 1)myBlue->PrintDouble(LocResultArr,myBlue->GetActObs(),
				    myBlue->GetActUnc()+1);

  }else if(Flag == 1){
    myBlue->SetInActiveEst(3);
    myBlue->SetInActiveEst(4);    
    myBlue->SetInActiveEst(5);
    myBlue->SetInActiveEst(6);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2012_095: The 2011 ATLAS Combination = %2i\n",
	   Flag);
    myBlue->PrintResult();
  }else if(Flag == 2){
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(1);    
    myBlue->SetInActiveEst(2);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2012_095: The 2011   CMS Combination = %2i\n",
	   Flag);
    myBlue->PrintResult();
  } else if(Flag == 3){
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2012_095: The combination for l+jets,");
    printf(" di-lepton and all-jets = %2i \n",Flag);
    myBlue->PrintResult();
  }
  // Delete Object
  delete myBlue; myBlue = NULL;
  LocCovRes->Delete(); LocCovRes = NULL;
  LocRhoRes->Delete(); LocRhoRes = NULL;
  LocWeight->Delete(); LocWeight = NULL;
  LocResult->Delete(); LocResult = NULL;
  return;
  }
