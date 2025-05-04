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

void B_NIMA_500_391(Int_t Flag = 0){
  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //
  // 0-3 == Lepton NON-universality and SINGLE uncertainty
  //  0: uncorrelated uncertainties Eq(38)
  //  1: 15% correlated uncertainties within Exp.A Eq(42)
  //     [Shows how to use InspectPair()]
  //  2: 99.5% correlated uncertainties within Exp.B Eq(46)
  //  3: -99.5% correlated uncertainties Eq(49)
  //
  // 10-13 == as 0-3 but imposing Lepton universality
  // 10: uncorrelated uncertainties Eq(40)
  // 11: 15% correlated uncertainties within Exp.A Eq(43)
  // 12: 99.5% correlated uncertainties within Exp.B Eq(48)
  // 13: -99.5% correlated uncertainties Eq(50)
  //
  // 4-7 == Lepton NON-universality and TWO uncertainties stat. + syst.
  //  4: 100% correlated syst. uncertainties Eq(53)
  //  5: un-correlated syst. uncertainties Eq(54)
  //  6: -100% correlated syst. uncertainties Eq(55)
  //  7: neglegted syst. uncertainties Eq(56)
  //
  // 14-17 == same as 4-7, but imposing Lepton universality
  // 14: 100% correlated syst. uncertainties Eq(57)
  // 15: un-correlated syst. uncertainties Eq(58)
  // 16: -100% correlated syst. uncertainties Eq(59)
  // 17: neglegted syst. uncertainties Eq(60)
  //----------------------------------------------------------------------------

  // Field housing the equations together with its index
  const Int_t NumEq[16] = {38, 42, 46, 49, 53, 54, 55, 56,
			   40, 43, 48, 50, 57, 58, 59, 60};
  Int_t IndEq;
  
  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 4;
  static const Int_t OneUnc = 1;
  static const Int_t TwoUnc = 2;
  Int_t NumUnc;
  Int_t NumObs;

  // Index for which estimate determines which observable
  Int_t IWhichObs[NumEst] = {0, 0, 0, 0};

  // Steer the input depending on Flag
  if(Flag ==  0 || Flag ==  1 || Flag ==  2 || Flag ==  3){
    IWhichObs[2] = 1;
    IWhichObs[3] = 1;
    NumObs = 2;
    NumUnc = OneUnc;
    IndEq  = Flag;
  }else if(Flag == 10 || Flag == 11 || Flag == 12 || Flag == 13){
    NumObs = 1;
    NumUnc = OneUnc;
    IndEq  = Flag-2;
  }else if(Flag ==  4 || Flag ==  5 || Flag ==  6 || Flag ==  7){
    IWhichObs[2] = 1;
    IWhichObs[3] = 1;
    NumObs = 2;
    NumUnc = TwoUnc;
    IndEq  = Flag;
  }else if(Flag == 14 || Flag == 15 || Flag == 16 || Flag == 17){
    NumObs = 1;
    NumUnc = TwoUnc;
    IndEq  = Flag-2;
  }else{
    IWhichObs[0] = 0;
    printf("B_NIMA_500_391: Not implemented Flag IGNORED %2i \n", Flag);
    return;
  }

  // Fill estimates for single uncertainty
  // Estimates 0=Bele(A), 0=Bele(B), 0=Btau(A), 0=Btau(b)
  static const Int_t LenXEstOne = NumEst * (OneUnc+1);
  Double_t XEstOne[LenXEstOne] = {
    10.5, 1.,
    13.5, 3.,
     9.5, 3.,
    14.0, 3.
  };

  // Fill estimates for two uncertainties
  static const Int_t LenXEstTwo = NumEst * (TwoUnc+1);
  Double_t XEstTwo[LenXEstTwo] = {
    10.5, 1.00, 0.00,
    13.5, 0.21, 2.99,
     9.5, 3.00, 0.00,
    14.0, 0.21, 2.99
  };

  // Fill estimates for two uncertainties with syst == 0
  Double_t XEstNull[LenXEstTwo]= {
    10.5, 1.00, 0.00,
    13.5, 0.21, 0.00,
     9.5, 3.00, 0.00,
    14.0, 0.21, 0.00
  };

  // Fill correlation for Flag == 1, 11
  static const Int_t LenCor = NumEst * NumEst;
  Double_t Cor1[LenCor] = {
    1.00, 0.15, 0.00, 0.00,
    0.15, 1.00, 0.00, 0.00,
    0.00, 0.00, 1.00, 0.00,
    0.00, 0.00, 0.00, 1.00
  };

  // Fill correlation for Flag == 2, 12
  Double_t Cor2[LenCor] = {
    1.00,  0.00, 0.00,  0.00,
    0.00,  1.00, 0.00, 0.995,
    0.00,  0.00, 1.00,  0.00,
    0.00, 0.995, 0.00,  1.00
  };

  // Fill correlation for Flag == 3, 13
  Double_t Cor3[LenCor] = {
    1.00,  0.00, 0.00,  0.00,
    0.00,  1.00, 0.00,-0.995,
    0.00,  0.00, 1.00,  0.00,
    0.00,-0.995, 0.00,  1.00
  };
  
  // Define formats for figures and latex file
  const TString FilBas = "B_NIMA_500_391";
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%4.3f";
  const TString ForRho = "%5.3f";
  const TString ForPul = ForRho;
  const TString ForUni = "None";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);

  // Fill all estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    if(Flag <= 3 || (Flag >= 10 && Flag-10 <= 3)){
      myBlue->FillEst(i,&XEstOne[ind]);
      ind = ind + OneUnc + 1;
    }else if(Flag == 7 || Flag == 17){
      myBlue->FillEst(i,&XEstNull[ind]);
      ind = ind + TwoUnc + 1;
    }else{
      myBlue->FillEst(i,&XEstTwo[ind]);
      ind = ind + TwoUnc + 1;
    }
  }

  // Fill all Correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(Flag == 1 || Flag == 11){
      myBlue->FillCor(k, &Cor1[0]);
    }else if(Flag == 2 || Flag == 12){
      myBlue->FillCor(k, &Cor2[0]);
    }else if(Flag == 3 || Flag == 13){
      myBlue->FillCor(k, &Cor3[0]);
    }else if((Flag == 4 || Flag == 14) && k == 1){
      myBlue->FillCor(k, 1.0);
    }else if((Flag == 6 || Flag == 16) && k == 1){
      myBlue->FillCor(k, -1.0);
    }else{
      myBlue->FillCor(k, 0.0);
    }
  }

  // Fix input and inspect input
  myBlue->FixInp();
  myBlue->PrintEst();
  myBlue->PrintCor(0);
  myBlue->PrintCov();
  myBlue->PrintRho();

  // Solve and inspect result
  myBlue->Solve();
  printf("... B_NIMA_500_391: Calculate Equation %2i \n", NumEq[IndEq]);
  myBlue->PrintResult();
  if(Flag < 10){
    myBlue->PrintRhoRes();
  }

  // Example of inspecting a Pair for Flag == 1
  if(Flag == 1){
    myBlue->ResetInp();
    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(3);
    myBlue->FixInp();
    myBlue->PrintEst(0);
    myBlue->PrintEst(1);   
    myBlue->InspectPair(0,1,FilBas);
  }

  // Delete Object
  delete myBlue;
  myBlue = NULL;
  return;
}
