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

void B_PRD79_092005(Int_t Flag = 0){
  
  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The paper combination,     vary rho(k>0)
  //  1: Combine all systematics,   vary rho(k=1)
  //  2: Combine all uncertainties, vary rho(k=0)
  //
  // In the Gaussian approximation the combined result can not be reached, since
  // the total correlation is bound to be smaller than about 0.4 if estimates
  // are uncorrelated for their statistical uncertainty.
  //
  // If the total uncertainty is used, the correlation would be about 0.7.
  //
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  2;
  static const Int_t NumUnc = 11;
  static const Int_t NumObs =  1;

  // The names
  TString NamEst[NumEst] = {"   LpJ", "   DIL"};
  //                              0         1         2         3
  TString NamUnc[NumUnc] = {"  Stat", "  RJES", "   Gen", "   PDF",
  //                              4         5         6         7
			    "  bJES", "  BSha", " gFrac", " I/FSR",
  //                              8         9        10   
			    "MCStat", "LepSca", "  Pile"};
  TString NamObs[1]      = {"  MTop"};
   
  // Preset according to Flag
  if(Flag == 0){
    printf("... B_PRD79_092005: ------------------------------------------ \n");
    printf("... B_PRD79_092005: Paper combination vary rho(k>0), Flag = %2i \n",
	   Flag);
  }else if(Flag == 1){
    printf("... B_PRD79_092005: ---------------------------------------- \n");
    printf("... B_PRD79_092005: Add systematics vary rho(k=1), Flag = %2i \n",
	   Flag);
    NamUnc[1] = "  Syst";
  }else if(Flag == 2){
    printf("... B_PRD79_092005: ----------------------------- \n");
    printf("... B_PRD79_092005: Add total vary rho, Flag = %2i \n", Flag);
    NamUnc[0] = "  Full";
  }else{
    printf("... B_PRD79_092005: ------------------------- \n");
    printf("... B_PRD79_092005: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Estimates 0=LpJ 1=Dil
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    171.8, 1.9, 0.7, 0.8, 0.3, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1,   
    171.2, 3.5, 3.5, 1.3, 0.5, 0.2, 0.3, 0.2, 0.2, 0.5, 0.3, 0.1
  };

  // Sum of systematics for adapting k=1 for Flag=1
  const Double_t SysLpJ = 1.1;
  const Double_t SysDil = 3.8;
  if(Flag == 1){
    XEst[0*(NumUnc+1) + 2] = SysLpJ;
    XEst[1*(NumUnc+1) + 2] = SysDil;
  }
  // Total uncertainty for adapting k=0 for Flag=2
  const Double_t TotLpJ = 2.2;
  const Double_t TotDil = 0.5*(5.3 + 5.1);
  if(Flag == 2){
    XEst[0*(NumUnc+1) + 1]= TotLpJ;
    XEst[1*(NumUnc+1) + 1]= TotDil;
  }
  
  // The CDF combination from the paper
  const Double_t TopVal = 171.9;
  const Double_t TopSta =   1.7;
  const Double_t TopSys =   1.1;
  const Double_t TopFul =   2.0;
  
  // Correlation for the systematic uncertainty
  const Double_t RhoSys = 1.0;

  // Correlation for best fit for Flag == 2
  const Double_t RhoTwo = 0.72;

  //-- Local Structures for Blue output
  // TMatrices
  TMatrixD* LocRes = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc = new TMatrixD(NumObs,1);
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);

  // Define file names
  const TString FilBas = "B_PRD79_092005";
  TString FilNam = "To be filled later";
  char Buffer[100];
  sprintf(Buffer,"%1i", Flag);
  FilNam = &Buffer[0];
  FilNam = FilBas + "_" + FilNam;

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.1f";
  const TString ForUnc = "%4.1f";
  const TString ForWei = ForUnc;
  const TString ForRho = "%4.2f";
  const TString ForPul = ForRho;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";
 
  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc);
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
    if(k == 0){
      myBlue->FillCor(k,0.);
    }else {
      myBlue->FillCor(k,RhoSys);
    }
  }

  // Adjust active uncertainties of the estimates according to Flag
  Int_t kl = 1;
  if(Flag != 0){
    for(Int_t k = 2; k<NumUnc; k++)myBlue->SetInActiveUnc(k);
    if(Flag == 2){
      myBlue->SetInActiveUnc(1);
      kl = 0;
    }
  }

  // Fix and print input
  myBlue->FixInp();
  myBlue->PrintEst();
  myBlue->SetQuiet();

  // Show the paper result
  printf("... B_PRD79_092005: The paper combination is: %s",
	 myBlue->GetNamObs(0).Data());
  printf(" = %5.1f (+- %4.1f +- %4.1f) = +- %4.1f\n",
	 TopVal, TopSta, TopSys, TopFul);

  // Solve for various rho assumptions from -1 to 0.96 in steps of 0.08
  printf("... B_PRD79_092005: The combination as a function of Rho is\n");
  for(Double_t r = 0.00; r < 1.00; r=r+0.08){
    myBlue->ReleaseInp(); 

    // Set Rho values for all k or k=kl
    if(Flag == 0){for(Int_t k = 1; k<NumUnc; k++)myBlue->SetRhoValUnc(k, r);
    }else{myBlue->SetRhoValUnc(kl, r);}

    // Fix and solve
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->GetResult(LocRes);
    myBlue->GetUncert(LocUnc);
    
    // Display the result as a function of rho
    printf("... B_PRD79_092005: Rho = %4.2f", r);
    if(Flag <= 1){
      printf(" M = %5.2f (+- %4.1f +- %4.1f) = +- %4.1f\n",
	     LocRes->operator()(0,0), LocRes->operator()(0,1),
	     LocRes->operator()(0,2), LocUnc->operator()(0,0));
    }else if(Flag == 2){
      printf(" M = %5.2f +- %4.1f\n",
	     LocRes->operator()(0,0), LocUnc->operator()(0,0));
      // For the best fit: rho=0.72
      if(r == RhoTwo){
	printf("... B_PRD79_092005: This is the best fit to the paper");
	printf("... combination\n");
	myBlue->PrintPull();
	myBlue->PrintCompatEst();
	myBlue->DisplayPair(0,1,FilNam);
	myBlue->LatexResult(FilNam);
	myBlue->GetRho(LocRho);
	printf("... B_PRD79_092005:\n");
      }
    }
  }
  
  // For Flag<=1, print compatibility for the best fit = full correlation
  if(Flag <= 1){
    myBlue->PrintPull();
    myBlue->PrintCompatEst();
    myBlue->LatexResult(FilNam);
    myBlue->GetRho(LocRho);
  }
  
  // See what we could be achieved depending on Flag
  printf("... B_PRD79_092005:\n");
  if(Flag == 0){
    printf("... B_PRD79_092005: The value of M = %5.1f cannot", TopVal);
    printf(" be reached since the correlation is bound to rho < %4.2f.\n",
	   LocRho->operator()(0,1));
  }else if(Flag == 1){
    printf("... B_PRD79_092005: If one first combines all systematic");
    printf(" uncertainties and assigns rho=1 to this sum, the\n");
    printf("... B_PRD79_092005: total estimator correlation is larger by");
    printf(" construction. Still, the value of M = %5.1f\n", TopVal);
    printf("... B_PRD79_092005: cannot be reached since the correlation");
    printf(" is still bound to rho < %4.2f.\n",
	   LocRho->operator()(0,1));
  }else if(Flag == 2){
    printf("... B_PRD79_092005: The value of M %5.1f",TopVal);
    printf(" is now reached for a total correlation of about rho = %4.2f.\n",
	   LocRho->operator()(0,1));
  }
  printf("... B_PRD79_092005:\n");

  // Delete Object and Matrices
  delete myBlue;  myBlue = NULL;

  LocRho->Delete(); LocRho = NULL;
  LocRes->Delete(); LocRes = NULL;
  LocUnc->Delete(); LocUnc = NULL;

  // Return
  return;
}
