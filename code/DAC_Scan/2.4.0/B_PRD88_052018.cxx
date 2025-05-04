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

void B_PRD88_052018(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The 2013 Tevatron Combination with reduced correlations
  //     [Shows how to use SetRhoRedUnc() and retrieve results]
  //     [Shows how to retrieve results with GetResult() and GetUncert()]
  //  1: The 2013 Tevatron Combination with all but "Rest" as fully correlated 
  //     [Shows how to use SolveScaRho()]
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 7;
  static const Int_t NumUnc = 5;
  static const Int_t NumObs = 1;

  // Set the names
  TString NamEst[NumEst] = {"    C[8]", "    C[9]", "   C[10]", 
			    "D[12-15]", "   D[16]", "   C[17]", "   D[18]"};
  TString NamUnc[NumUnc] = {"   Rest", "    PDF", "RadCorC", 
			    "RadCorU", " GammaW"};
  TString NamObs[NumObs] = {"  WMass"};

  //-- Local Structures for BLUE output
  // TMatrices
  Int_t iok = 0;
  TMatrixD* LocWei = new TMatrixD(NumEst,NumObs);
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);

  TMatrixD* LocRes = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc = new TMatrixD(NumObs,1);
  //-- End

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_PRD88_052018: -----------------------------------------");
    printf("-------------------\n");
    printf("... B_PRD88_052018: The 2013 Tevatron combination with");
    printf(" reduced correlations = %2i \n",Flag);
  }else if(Flag == 1){
    printf("... B_PRD88_052018: --------------------------------------------");
    printf("----------------------------------------------------------- \n");
    printf("... B_PRD88_052018: The 2013 Tevatron combination with");
    printf(" all but the uncertainty Rest=Total-Specific");
    printf(" as fully correlated = %2i \n",Flag);
  }else{
    printf("... B_PRD88_052018: ------------------------------- \n");
    printf("... B_PRD88_052018: Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  // The original input numbers
  Double_t ValMas[NumEst] = {79927.7, 80377.3, 80470.5, 80478.5, 80401.8, 80387.3, 80368.6};
  Double_t PDFUnc[NumEst] = {   60.0,    50.0,    15.0,     8.0,    10.0,    10.0,    11.0};
  Double_t RadUnc[NumEst] = {   10.0,    20.0,     5.0,    12.0,     7.0,     4.0,     7.0};
  Double_t GaWUnc[NumEst] = {    0.5,     1.4,     0.3,     1.5,     0.4,     0.2,     0.5};
  Double_t TotUnc[NumEst] = {  390.0,   181.0,    89.0,    84.0,    43.0,    19.0,    26.0};
  
  // The fully correlated constant value for Radiation is 3.5 MeV
  Double_t RadCor = 3.5;
 
  // Now calculate the separation into the uncertainty sources
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {0};
  Int_t idn = 0;
  for(Int_t i = 0; i<NumEst; i++){
    // Value
    XEst[idn] = ValMas[i];
    idn = idn + 1;
    // Rest 
    XEst[idn] = TMath::Sqrt(TotUnc[i]*TotUnc[i] - PDFUnc[i]*PDFUnc[i] -
			    RadUnc[i]*RadUnc[i] - GaWUnc[i]*GaWUnc[i]);
    idn = idn + 1;
    //PDF
    XEst[idn] = PDFUnc[i];
    idn = idn + 1;
    //RadCorc
    XEst[idn] = RadCor;
    idn = idn + 1;
    //RadCoru
    XEst[idn] = TMath::Sqrt(RadUnc[i]*RadUnc[i] - RadCor*RadCor);
    idn = idn + 1;
    // GammaW
    XEst[idn] = GaWUnc[i];
    idn = idn + 1;
  }
  
  // The default correlations: 0    1    2    3    4
  Double_t RhoVal[NumUnc] = {0.0, 1.0, 1.0, 0.0, 1.0};

  // Make 3 fully correlated for Flag == 1
  if(Flag == 1) RhoVal[3] = 1.0;

  // For Flag == 0
  // The correlations will be changed to reduced correlation for 
  // the PDF uncertainty: k = 1
  // The split for k = 2, 3 is -- the same as reduced correlations -- for all
  // estimates but for 4, 6
  
  // Matrix for the 'uncorrelated' Radiation part k=3 that is still 
  // correlated for Est 4, 6
  Double_t RhoRad[NumEst*NumEst] = {
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0
  };

  // Define formats for figures and latex file
  const TString FilBas = "B_PRD88_052018";
  const TString ForVal = "%7.1f";
  const TString ForUnc = "%4.1f";
  const TString ForWei = "%5.3f";
  const TString ForRho = "%5.3f";
  const TString ForPul = "%3.3f";
  const TString ForUni = "MeV";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);
  
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
    if(Flag == 0 && k == 3){myBlue->FillCor(k,&RhoRad[0]);
    }else{myBlue->FillCor(k,RhoVal[k]);
    }
  }
    
  // Solve according to Flag
  if(Flag == 0){

    // The central paper result with reduced correlations
    myBlue->SetRhoRedUnc(1);
    myBlue->FixInp();

    // The correlation matrices per source
    for(Int_t k = 0; k<NumUnc; k++){
      printf("... B_PRD88_052018: The correlations, for source %s \n",
	     NamUnc[k].Data());
      myBlue->GetCor(k, LocRho);
      myBlue->PrintMatrix(LocRho, ForRho);
    }
    
    // The estimates
    printf("... B_PRD88_052018: The estimates, see Table 1 \n");
    myBlue->PrintEst();
    myBlue->Solve();

    // The weights
    printf("... B_PRD88_052018: The weights, see Table 2 \n");
    iok = myBlue->GetWeight(LocWei);
    LocWei->operator*=(100);
    if(iok == 1)myBlue->PrintMatrix(LocWei,"%5.1f%%");

    printf("... B_PRD88_052018: The correlations, see Table 3 \n");
    printf("... B_PRD88_052018: Small differences for the least precise");
    printf(" measurements [8-10] are observed. \n");
    printf("... B_PRD88_052018:                C[8]    C[9]    C[10] D[12-15]");
    printf("    D[16]    C[17]    D[18] \n");
    iok = myBlue->GetRho(LocRho);
    myBlue->PrintMatrix(LocRho,"   %5.3f");

    printf("... B_PRD88_052018: The chi-squared of the fit,");
    printf(" see page 9 right column \n");
    myBlue->PrintChiPro();

    printf("... B_PRD88_052018: The result, see page 9 right column \n");
    myBlue->PrintResult();

    // Print out result
    iok = myBlue->GetResult(LocRes);
    iok = myBlue->GetUncert(LocUnc);
    printf("... B_PRD88_052018: ------------------------------------");
    printf("-------------------------------- \n");
    printf("... B_PRD88_052018: The default == REDUCED CORRELATIONS result");
    printf(" is: %s = %5.0f +- %2.0f\n", 
	   NamObs[0].Data(), LocRes->operator()(0,0), LocUnc->operator()(0,0));
    printf("... B_PRD88_052018: ------------------------------------");
    printf("-------------------------------- \n");
    myBlue->LatexResult(FilBas,ForVal,ForUnc,ForWei,ForRho,ForPul);
    myBlue->DisplayResult(0,FilBas);

  }else if(Flag == 1){

    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

    // Print out result
    iok = myBlue->GetResult(LocRes);
    iok = myBlue->GetUncert(LocUnc);
    printf("... B_PRD88_052018: ------------------------------------");
    printf("-------------------------------- \n");
    printf("... B_PRD88_052018: The                full correlation result");
    printf(" is: %s = %5.0f +- %2.0f\n", 
	   NamObs[0].Data(), LocRes->operator()(0,0), LocUnc->operator()(0,0));
    printf("... B_PRD88_052018: The difference to the default result is");
    printf(" just one MeV \n");
    printf("... B_PRD88_052018: ------------------------------------");
    printf("-------------------------------- \n");

    // The rho variations
    myBlue->ReleaseInp();
    myBlue->FixInp();
    myBlue->SolveScaRho(0);
    myBlue->PrintScaRho();
  }

  // Delete Object
  delete myBlue;
  myBlue = NULL;  

  LocWei->Delete(); LocWei = NULL;
  LocRho->Delete(); LocRho = NULL;
  LocRes->Delete(); LocRes = NULL;
  LocUnc->Delete(); LocUnc = NULL;

  // Return
  return;
}
