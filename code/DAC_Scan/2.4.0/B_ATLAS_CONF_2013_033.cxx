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

void B_ATLAS_CONF_2013_033(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The   LHC combination
  //  1: The ATLAS combination
  //  2: The ATLAS and CMS 2011 only combination
  //  3: The single lepton only combination
  //  4: The 2010 vs 2011 compatibility combination
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  8;
  static const Int_t NumUnc = 13;
  static const Int_t MaxObs =  4;
  Int_t NumObs =  2;
  
  // The names
  TString NamEst[NumEst] = {"F0A10LJ", "F0A11LJ", "F0A10DL", "F0C11LJ",
			    "FLA10LJ", "FLA11LJ", "FLA10DL", "FLC11LJ"};
  TString NamUnc[NumUnc] = {"   Stat", " DetMod", "    JES", "   Lumi", 
			    "     MC", "    Rad", "   mtop", "    PDF", 
			    " BMCQCD", " BMCWpJ", " BMCOth", "  BData", 
			    "   Meth"};
  TString NamObs[MaxObs] = {"     F0", "     FL", " FL(10)", " FL(11)"};

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {0, 0, 0, 0, 1, 1, 1, 1};
  
  // Preset according to Flag
  if(Flag == 0){
    printf("... B_ATLAS_CONF_2013_033: The 2011 LHC combination,");
    printf(" Flag = %2i \n", Flag);
  }else if(Flag == 1){
    printf("... B_ATLAS_CONF_2013_033: The ATLAS combination,");
    printf(" Flag = %2i \n", Flag);
  }else if(Flag == 2){
    printf("... B_ATLAS_CONF_2013_033: The ATLAS and CMS 2011 only");
    printf("combination, Flag = %2i \n", Flag);
  }else if(Flag == 3){
    printf("... B_ATLAS_CONF_2013_033: The single lepton only combination,");
    printf(" Flag = %2i \n", Flag);
  }else if(Flag == 4){
    printf("... B_ATLAS_CONF_2013_033: The 2010 vs. 2011 compatibility");
    printf(" combination, Flag = %2i \n", Flag);
    NumObs = 4;
    NamObs[0] = " F0(10)";
    NamObs[1] = " F0(11)";
    IWhichObs[0] = 0;
    IWhichObs[1] = 1;
    IWhichObs[2] = 0;
    IWhichObs[3] = 1;
    IWhichObs[4] = 2;
    IWhichObs[5] = 3;
    IWhichObs[6] = 2;
    IWhichObs[7] = 3;
  }else{
    printf("... ATLAS_CONF_2013_033: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Character to set the file names
  char Buffer[100];

  // Estimates 0-6 
  // First the F0
  // 0 == ATLAS_2010 lepton+jets
  // 1 == ATLAS_2011 lepton+jets
  // 2 == ATLAS_2011   di-lepton
  // 3 ==   CMS_2011 lepton+jets
  // 4-7 == the FL using the same order

  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //         0      1      2      3      4      5      6      7      8      9     10     11     12
    //      Stat DetMod    JES   Lumi     MC    Rad   mtop    PDF BMCQCD BMCWpJ BMCOth  BData   Meth
    0.652, 0.134, 0.047, 0.043, 0.013, 0.020, 0.033, 0.018, 0.004, 0.000, 0.023, 0.005, 0.035, 0.026,
    0.642, 0.030, 0.032, 0.027, 0.012, 0.019, 0.030, 0.027, 0.009, 0.000, 0.000, 0.008, 0.027, 0.015,
    0.744, 0.050, 0.012, 0.056, 0.002, 0.023, 0.028, 0.028, 0.028, 0.000, 0.000, 0.006, 0.018, 0.032,
    0.567, 0.074, 0.020, 0.018, 0.000, 0.000, 0.026, 0.009, 0.001, 0.007, 0.020, 0.019, 0.000, 0.000,
    0.359, 0.088, 0.029, 0.027, 0.006, 0.012, 0.016, 0.012, 0.002, 0.000, 0.013, 0.003, 0.023, 0.017,
    0.344, 0.020, 0.019, 0.014, 0.005, 0.014, 0.019, 0.014, 0.005, 0.000, 0.000, 0.005, 0.017, 0.011,
    0.276, 0.031, 0.005, 0.036, 0.001, 0.015, 0.014, 0.016, 0.015, 0.000, 0.000, 0.004, 0.011, 0.016,
    0.393, 0.045, 0.015, 0.011, 0.000, 0.000, 0.008, 0.010, 0.001, 0.002, 0.006, 0.007, 0.000, 0.000
  };
 
  // Due to the calculation of correlation from uncertainties within the note,
  // all correlation matrices for all uncertainty sources are different
  static const Int_t LenCor = NumEst * NumEst;
  Double_t Cor00[LenCor] = {
    +1.00,  0.00,  0.00,  0.00, -0.94,  0.00,  0.00,  0.00,
    +0.00,  1.00,  0.00,  0.00,  0.00, -0.91,  0.00,  0.00,
    +0.00,  0.00,  1.00,  0.00,  0.00,  0.00, -0.91,  0.00,
    +0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00, -0.94,
    -0.94,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00, 
    +0.00, -0.91,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00, 
    +0.00,  0.00, -0.91,  0.00,  0.00,  0.00,  1.00,  0.00, 
    +0.00,  0.00,  0.00, -0.94,  0.00,  0.00,  0.00,  1.00
  };
  Double_t Cor01[LenCor] = {
    +1.00,  1.00,  1.00,  0.00, -0.91, -1.00, -1.00,  0.00,
    +1.00,  1.00,  1.00,  0.00, -1.00, -0.78, -1.00,  0.00,
    +1.00,  1.00,  1.00,  0.00, -1.00, -1.00, -0.89,  0.00,
    +0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00, -0.95,
    -0.91, -1.00, -1.00,  0.00,  1.00,  1.00,  1.00,  0.00, 
    -1.00, -0.78, -1.00,  0.00,  1.00,  1.00,  1.00,  0.00, 
    -1.00, -1.00, -0.89,  0.00,  1.00,  1.00,  1.00,  0.00, 
    +0.00,  0.00,  0.00, -0.95,  0.00,  0.00,  0.00,  1.00
  };
  Double_t Cor02[LenCor] = {
    +1.00,  1.00,  1.00,  0.00, -0.92, -1.00, -1.00,  0.00,
    +1.00,  1.00,  1.00,  0.00, -1.00, -0.31, -1.00,  0.00,
    +1.00,  1.00,  1.00,  0.00, -1.00, -1.00, -0.49,  0.00,
    +0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00, -0.99,
    -0.92, -1.00, -1.00,  0.00,  1.00,  1.00,  1.00,  0.00, 
    -1.00, -0.31, -1.00,  0.00,  1.00,  1.00,  1.00,  0.00, 
    -1.00, -1.00, -0.49,  0.00,  1.00,  1.00,  1.00,  0.00, 
    +0.00,  0.00,  0.00, -0.99,  0.00,  0.00,  0.00,  1.00
  };
  Double_t Cor03[LenCor] = {
    +1.00,  1.00,  1.00,  1.00, -0.86, -1.00, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -0.86, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -0.94, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -1.00,  0.00,
    -0.86, -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -0.86, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -0.94, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -1.00,  0.00,  1.00,  1.00,  1.00,  1.00
  };
  Double_t Cor04[LenCor] = {
    +1.00,  1.00,  1.00,  1.00, -0.89, -1.00, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -0.92, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -0.92, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -1.00,  0.00,
    -0.89, -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -0.92, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -0.92, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -1.00,  0.00,  1.00,  1.00,  1.00,  1.00
  };
  Double_t Cor05[LenCor] = {
    +1.00,  1.00,  1.00,  0.50, -0.88, -1.00, -1.00, -0.50,
    +1.00,  1.00,  1.00,  0.50, -1.00, -0.58, -1.00, -0.50,
    +1.00,  1.00,  1.00,  0.50, -1.00, -1.00, -0.85, -0.50,
    +0.50,  0.50,  0.50,  1.00, -0.50, -0.50, -0.50,  0.21,
    -0.88, -1.00, -1.00, -0.50,  1.00,  1.00,  1.00,  0.50, 
    -1.00, -0.58, -1.00, -0.50,  1.00,  1.00,  1.00,  0.50, 
    -1.00, -1.00, -0.85, -0.50,  1.00,  1.00,  1.00,  0.50, 
    -0.50, -0.50, -0.50,  0.21,  0.50,  0.50,  0.50,  1.00
  };
  Double_t Cor06[LenCor] = {
    +1.00,  1.00,  1.00,  1.00, -0.92, -1.00, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -0.09, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -0.44, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -1.00, -0.87,
    -0.92, -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -0.09, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -0.44, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -1.00, -0.87,  1.00,  1.00,  1.00,  1.00
  };
  Double_t Cor07[LenCor] = {
    +1.00,  1.00,  1.00,  1.00, -0.92, -1.00, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -0.88, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -0.88, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -1.00, -1.00,
    -0.92, -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -0.88, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -0.88, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00
  };
  Double_t Cor08[LenCor] = {
    +1.00,  1.00,  1.00,  1.00,  0.00, -1.00, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00,  0.00, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00,  0.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -1.00, -1.00,
     0.00, -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00,  0.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00,  0.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00
  };
  Double_t Cor09[LenCor] = {
    +1.00,  1.00,  1.00,  1.00, -0.91, -1.00, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00,  0.00, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00,  0.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -1.00,  1.00,
    -0.91, -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00,  0.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00,  0.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00,  1.00
  };
  Double_t Cor10[LenCor] = {
    +1.00,  1.00,  1.00,  1.00, -0.95, -1.00, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -0.89, -1.00, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -0.91, -1.00,
    +1.00,  1.00,  1.00,  1.00, -1.00, -1.00, -1.00, -0.59,
    -0.95, -1.00, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -0.89, -1.00, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -0.91, -1.00,  1.00,  1.00,  1.00,  1.00, 
    -1.00, -1.00, -1.00, -0.59,  1.00,  1.00,  1.00,  1.00
  };
  Double_t Cor11[LenCor] = {
    +1.00,  0.00,  0.00,  0.00, -0.93,  0.00,  0.00,  0.00,
    +0.00,  1.00,  0.00,  0.00,  0.00, -0.93,  0.00,  0.00,
    +0.00,  0.00,  1.00,  0.00,  0.00,  0.00, -1.00,  0.00,
    +0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
    -0.93,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00, 
    +0.00, -0.93,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00, 
    +0.00,  0.00, -1.00,  0.00,  0.00,  0.00,  1.00,  0.00, 
    +0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00
  };
  Double_t Cor12[LenCor] = {
    +1.00,  0.00,  0.00,  0.00, -0.96,  0.00,  0.00,  0.00,
    +0.00,  1.00,  0.00,  0.00,  0.00, -0.78,  0.00,  0.00,
    +0.00,  0.00,  1.00,  0.00,  0.00,  0.00, -0.95,  0.00,
    +0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
    -0.96,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00, 
    +0.00, -0.78,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00, 
    +0.00,  0.00, -0.95,  0.00,  0.00,  0.00,  1.00,  0.00, 
    +0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00
  };
    
  //-- Local Structures for Blue output
  // TMatrices
  TMatrixD* LocRho    = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRhoRes = new TMatrixD(NumObs,NumObs);
  TMatrixD* LocWeight = new TMatrixD(NumEst,NumObs);
  //-- End

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.3f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%5.3f";
  const TString ForRho = ForWei;
  const TString ForPul = ForVal;
  const TString ForChi = "%5.3f";
  const TString ForUni = "";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

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
  
  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(k == 0){myBlue->FillCor(k,&Cor00[0]);
    }else if(k ==  1){myBlue->FillCor(k,&Cor01[0]);
    }else if(k ==  2){myBlue->FillCor(k,&Cor02[0]);
    }else if(k ==  3){myBlue->FillCor(k,&Cor03[0]);
    }else if(k ==  4){myBlue->FillCor(k,&Cor04[0]);
    }else if(k ==  5){myBlue->FillCor(k,&Cor05[0]);
    }else if(k ==  6){myBlue->FillCor(k,&Cor06[0]);
    }else if(k ==  7){myBlue->FillCor(k,&Cor07[0]);
    }else if(k ==  8){myBlue->FillCor(k,&Cor08[0]);
    }else if(k ==  9){myBlue->FillCor(k,&Cor09[0]);
    }else if(k == 10){myBlue->FillCor(k,&Cor10[0]);
    }else if(k == 11){myBlue->FillCor(k,&Cor11[0]);
    }else if(k == 12){myBlue->FillCor(k,&Cor12[0]);
    }
  }
  
  if(Flag == 0){

    myBlue->FixInp();
    printf("... B_ATLAS_CONF_2013_033: The input in the following order \n");
    printf("... B_ATLAS_CONF_2013_033: 0 == ATLAS_2010 lepton+jets \n");
    printf("... B_ATLAS_CONF_2013_033: 1 == ATLAS_2011 lepton+jets \n");
    printf("... B_ATLAS_CONF_2013_033: 2 == ATLAS_2011   di-lepton \n");
    printf("... B_ATLAS_CONF_2013_033: 3 ==   CMS_2011 lepton+jets \n");
    printf("... B_ATLAS_CONF_2013_033: 4-7 == FL same order \n");
    for(Int_t i = 0; i<NumEst; i++){
      myBlue->PrintEst(i);
    }
    sprintf(Buffer,"B_ATLAS_CONF_2013_033_%i",Flag);
    myBlue->PrintCompatEst(Buffer);
    printf("... B_ATLAS_CONF_2013_033: The correlations");
    printf(" of the estimates in %%\n");
    myBlue->GetRho(LocRho);
    LocRho->operator*=(100);
    myBlue->PrintMatrix(LocRho,"%+4.0f");

    myBlue->Solve();
    printf("... B_ATLAS_CONF_2013_033: The 2011 LHC Combination");
    printf(" Flag = %2i.\n",Flag);
    printf("... B_ATLAS_CONF_2013_033: For the weight matrix and result");
    printf(" and the uncertainties see Table 7 and 6 \n");
    myBlue->PrintResult();
    myBlue->LatexResult("B_ATLAS_CONF_2013_033_0");

    printf("... B_ATLAS_CONF_2013_033: The correlations");
    printf(" of the observables in %%\n");
    myBlue->GetRhoRes(LocRhoRes);
    LocRhoRes->operator*=(100);
    myBlue->PrintMatrix(LocRhoRes,"%+4.0f");

    printf("... B_ATLAS_CONF_2013_033: A difference is observed in the weights");
    printf(" that is under discussion with the authors \n");
    printf("... B_ATLAS_CONF_2013_033: The observed weights\n");
    myBlue->GetWeight(LocWeight);
    myBlue->PrintMatrix(LocWeight," %+5.3f");
    myBlue->PrintCompatObs();

  }else if(Flag == 1){

    myBlue->SetInActiveEst(3);
    myBlue->SetInActiveEst(7);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2013_033: The ATLAS Combination");
    printf(" Flag = %2i \n",Flag);
    myBlue->PrintResult();

  }else if(Flag == 2){

    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(4);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2013_033: The ATLAS and CMS 2011 only combination,");
    printf(" Flag = %2i \n",Flag);
    myBlue->PrintResult();

  }else if(Flag == 3){

    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(6);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2013_033: The single lepton only combination,");
    printf(" Flag = %2i \n",Flag);
    myBlue->PrintResult();

  }else if(Flag == 4){

    myBlue->FixInp();
    myBlue->PrintCompatEst();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2013_033: The 2010 vs 2011 compatibility");
    printf(" combination, Flag = %2i \n",Flag);
    myBlue->PrintResult();
    myBlue->LatexResult("B_ATLAS_CONF_2013_033_4");
    printf("... B_ATLAS_CONF_2013_033: Only the compatibility pairs");
    printf(" F0(2010,2011)=( F0(10), F0(11)) and");
    printf(" FL(2010,2011)=(FL(10), FL(11)) make sense \n");
    myBlue->PrintCompatObs();
  }

  // Delete Object
  delete myBlue; myBlue = NULL;
  LocRho->Delete(); LocRho = NULL;
  LocRhoRes->Delete(); LocRhoRes = NULL;
  LocWeight->Delete(); LocWeight = NULL;
  return;
}


