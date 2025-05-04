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

void B_ATLAS_CONF_2012_134(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The LHC Combination
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  2;
  static const Int_t NumUnc = 13;
  static const Int_t NumObs =  1;

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {0, 0};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_ATLAS_CONF_2012_134: The 2011 LHC Combination = %2i \n",Flag);
  }else{
    printf("... ATLAS_CONF_2012_134: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Estimates, 0 == ATLAS, 1 == CMS
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //       0    1    2    3    4    5    6    7    8    9   10   11   12
    177.0, 3.2, 2.7, 5.3, 4.2, 1.3, 0.8, 1.9, 1.5, 1.6, 2.4, 1.0, 5.3, 4.3,
    165.8, 2.2, 3.5, 8.8, 1.1, 2.2, 4.1, 4.1, 3.4, 1.6, 0.0, 1.0, 5.1, 5.9
  };
 
  // Correlations of uncertainty sources
  //      0,  1,  2,  7,  9, 12 <==}  rho=0,
  //  3,  4,  5,  6,  8, 10, 11 <==}  rho=1

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.1f";
  const TString ForUnc = "%4.1f";
  const TString ForWei = "%4.3f";
  const TString ForRho = "%4.2f";
  const TString ForPul = ForRho;
  const TString ForUni = "pb";

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
    if(k == 0 || k == 1 || k == 2 || k == 7 || k == 9|| k == 12){
      myBlue->FillCor(k,0.0);
    }else{
      myBlue->FillCor(k,1.0);
    }
  }

  // Calculate accroding to Flag
  if(Flag == 0){
    myBlue->FixInp();
    myBlue->PrintEst(0);
    myBlue->PrintEst(1);
    myBlue->PrintRho();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2012_134: The 2011 LHC Combination = %2i \n",Flag);
    printf("... B_ATLAS_CONF_2012_134: For the results see Table 1 \n");
    myBlue->PrintResult();
  }
  delete myBlue;
}


