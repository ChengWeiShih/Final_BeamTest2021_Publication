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

void B_PRD41_982(Int_t Flag = 0){
  
  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The paper combination
  //     [Shows how to use names and InspectLike()]
  //  1: The combination with the default choice for relative uncertainties
  //  2: The Blue combination with absolute uncertainties
  //----------------------------------------------------------------------------
  
  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  8;
  static const Int_t NumUnc =  2;
  static const Int_t NumObs = 1;

  // The names
  TString NamEst[NumEst] = {"  DELCO", "    HRS", "JADE(e)", "JADE(m)", 
			    "JADE(d)", "   MAC ", "  MarkJ", "  TASSO"};
  TString NamUnc[NumUnc] = {"   Stat", "   Syst"};
  TString NamObs[1]      = {" B-Life"};
   
  // Preset according to Flag
  if(Flag == 0){
    printf("... B_PRD41_982: -------------------------------- \n");
    printf("... B_PRD41_982: The paper combination, Flag = %2i \n", Flag);
  }else if(Flag == 1){
    printf("... B_PRD41_982: ---------------------------------");
    printf("----------------------------------------------- \n");
    printf("... B_PRD41_982: The combination with my default");
    printf(" choice for the relative uncertainties, Flag = %2i \n", Flag);
  }else if(Flag == 2){
    printf("... B_PRD41_982: -------------------------------------- \n");
    printf("... B_PRD41_982: The normal Blue combination, Flag = %2i \n", Flag);
  }else{
    printf("... B_PRD41_982: -------------------------- \n");
    printf("... B_PRD41_982: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // The input from Table II
  Double_t Values[NumEst] = {1.17, 1.02, 1.27, 1.36, 1.46, 1.29, 0.98, 1.35};
  Double_t StaUpw[NumEst] = {0.27, 0.41, 0.35, 0.32, 0.22, 0.20, 0.12, 0.10};
  Double_t StaDnw[NumEst] = {0.22, 0.37, 0.29, 0.27, 0.21, 0.20, 0.12, 0.10};
  Double_t SysUpw[NumEst] = {0.17, 0.00, 0.24, 0.21, 0.34, 0.07, 0.13, 0.24};
  Double_t SysDnw[NumEst] = {0.16, 0.00, 0.24, 0.21, 0.34, 0.07, 0.13, 0.24};
  Double_t StaCo2[NumEst] = {1.30, 1.30, 1.40, 1.40, 0.50, 1.30, 1.30, 1.00};
  Double_t StaCo0[NumEst] = {1.60, 1.70, 1.80, 1.80, 0.70, 1.70, 0.90, 1.10};
  Double_t SysCo2[NumEst] = {0.09, 0.00, 0.17, 0.15, 0.20, 0.15, 0.10, 0.16};
  Double_t SysCo0[NumEst] = {0.14, 0.00, 0.10, 0.07, 0.18, 0.07, 0.07, 0.13};
  Double_t RhoSys = 0.5, RhoUnc = 0.3;

  // Calculate estimates with symmeterised uncertainties
  // Calculate the Ni = NumEff from Eq. 9. Rescale the coefficients.
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {0};
  Int_t ind = 0;
  Double_t NumEff = 0;
  for(Int_t i = 0; i<NumEst; i++){
    // Fill estimates
    XEst[  ind] = Values[i];
    XEst[ind+1] = (StaUpw[i]+StaDnw[i])/2;
    XEst[ind+2] = (SysUpw[i]+SysDnw[i])/2;

    // Rescale stat coefficients for relative uncertainties
    NumEff = (StaCo0[i] * StaCo0[i] + StaCo2[i] * XEst[ind] * XEst[ind]) /
      (XEst[ind+1] * XEst[ind+1]);
    StaCo0[i] = StaCo0[i] * StaCo0[i] / NumEff;
    StaCo2[i] = StaCo2[i] / NumEff;

    // Rescale syst coefficients for relative uncertainties
    // Not for estimate 1 which has no syst
    if(i != 1){
      SysCo0[i] = SysCo0[i] * SysCo0[i];
      SysCo2[i] = SysCo2[i] * SysCo2[i];
    }

    // Increment the pointer
    ind = ind + 3;
  }

  // Local structure for Blue output and systematic uncertainties Table 3
  // Arrays:
  Int_t iok = 0;
  static const Int_t LenRes = NumObs*(NumUnc+1);
  Double_t ValDef[LenRes] = {0.};
  Double_t ValAct[LenRes] = {0.};
  Double_t UncDef[LenRes] = {0l};
  static const Int_t NumSys = 2;
  Double_t ValDif[NumSys] = {0.};
  TString LocNamEst = " ";
  // Matrices
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);
  // End

  // Define formats for figures and latex file
  const TString FilBas = "B_PRD41_982";
  TString FilNam = "To be filled later";
  char Buffer[100];
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%4.3f";
  const TString ForRho = "%5.3f";
  const TString ForPul = "%4.2f";
  const TString ForUni = "None";

  // Construct object
  Blue *myBlue = new Blue(NumEst, NumUnc);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);

  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);

  // Fill estimates and correlations
  ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }
  for(Int_t k = 0; k<NumUnc; k++){
    if(k == 0){
      myBlue->FillCor(k,0.0);
    }else{
      myBlue->FillCor(k,RhoSys);
    }
  }
  
  // Set relative uncertainties if needed
  // Flag == 0: The functions are both a0 + 0*x^1 + a2*x^2
  // Flag == 1: Use default case
  // Flag == 2: Do nothing
  static const Int_t LenCof = 3;
  Double_t ActCof[LenCof] = {0};
  if(Flag == 0){
    for(Int_t i = 0; i<NumEst; i++){
      for(Int_t k = 0; k<NumUnc; k++){
	if(k == 0){
	  ActCof[0] = StaCo0[i];
	  ActCof[2] = StaCo2[i];
	}else{
	  ActCof[0] = SysCo0[i];
	  ActCof[2] = SysCo2[i];
	}
	myBlue->SetRelUnc(i,k,&ActCof[0]);
      }
    }
  }else if(Flag == 1){
    myBlue->SetRelUnc();
  }

  // Do the relative uncertainty (Flag <= 1) or the normal Blue (Flag == 2) 
  if(Flag <= 1){
    // Do the combinations
    // Loop over rho variations first save the difference
    printf("... B_PRD41_982: Perform the rho variations first:");
    printf("... B_PRD41_982: Rho = %3.1f +- %3.1f \n", RhoSys, RhoUnc);
    for(Int_t l = 0; l<2; l++){
      if(l == 0){
	myBlue->SetRhoValUnc(1,RhoSys-RhoUnc);
      }else{
	myBlue->SetRhoValUnc(1,RhoSys+RhoUnc);
      }
      myBlue->FixInp();

      // Retrieve one name
      LocNamEst = myBlue->GetNamEst(2);
      printf("... B_PRD41_982: Estimator 2 is called \"%s\" \n",
	     LocNamEst.Data());

      // Solve
      myBlue->SolveRelUnc(0.1);

      // Save to calculate difference to default
      iok = myBlue->GetResult(ValAct);      
      if(iok == 1)ValDif[l] = ValAct[0];
      myBlue->ReleaseInp();
      myBlue->SetNotRhoValUnc(1);
    }
    
    // Produce the central value, print the result
    printf("... B_PRD41_982: Now Produce the central result \n");
    myBlue->FixInp();
    myBlue->PrintNamEst();
    myBlue->PrintNamUnc();
    myBlue->PrintNamObs();
    myBlue->SolveRelUnc(0.1);
    myBlue->PrintCofRelUnc();
    myBlue->PrintEst();
    myBlue->PrintResult();

    // The correlation matrix
    printf("\n");
    printf("... B_PRD41_982: The correlation matrix \n");
    myBlue->GetRho(LocRho);
    myBlue->PrintMatrix(LocRho,ForRho);    

    // Compare to likelihood fit
    if(Flag ==0){
      sprintf(Buffer,"%s_%i", FilBas.Data(),Flag);
      FilNam = &Buffer[0];
      myBlue->InspectLike(0,FilNam);
      myBlue->PrintInspectLike();
    }
    
    // We are done
    iok = myBlue->GetResult(ValDef);
    iok = myBlue->GetUncert(UncDef);
    printf("... B_PRD41_982: --------------------------------------- \n");
    if(Flag == 0){      
      printf("... B_PRD41_982:  The final paper result combination is: \n");
    }else{
      printf("... B_PRD41_982: The result using the default RelUnc is: \n");
    }
    printf("... B_PRD41_982: %s = %5.3f +- %5.3f +- %5.3f ps \n",
	   myBlue->GetNamObs(0).Data(),ValDef[0],UncDef[0],
	   TMath::Abs(ValDif[0] - ValDif[1])/2);
    printf("... B_PRD41_982: --------------------------------------- \n");
    printf(" \n");   

    // Get the Latex and Display output
    if(Flag == 0){
      myBlue->LatexResult(FilBas);
      myBlue->DisplayResult(0,FilBas);
    }   
  }else{
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintEst();
    myBlue->PrintResult();

    // The correlation matrix
    printf("\n");
    printf("... B_PRD41_982: The correlation matrix \n");
    myBlue->GetRho(LocRho);
    myBlue->PrintMatrix(LocRho,ForRho);    

    // We are done
    iok = myBlue->GetResult(ValDef);
    iok = myBlue->GetUncert(UncDef);
    printf("... B_PRD41_982: ------------------------------------------ \n");
    printf("... B_PRD41_982: The BLUE result without rho variations is: \n"); 
    printf("... B_PRD41_982: %s = %5.3f +- %5.3f +- --- ps \n",
	    myBlue->GetNamObs(0).Data(),ValDef[0], UncDef[0]);
    printf("... B_PRD41_982: ------------------------------------------ \n");
    printf(" \n");   
  }

  // Delete Object and TMatrices
  delete myBlue;
  myBlue = NULL;

  delete LocRho; LocRho = NULL;
  
  // Return
  return;
}
