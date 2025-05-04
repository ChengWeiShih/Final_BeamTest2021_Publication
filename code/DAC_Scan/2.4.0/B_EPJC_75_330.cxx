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

void B_EPJC_75_330(Int_t Flag = 0){
  
  //--------------------------------------------------------------------------
  // [Shows how to fill a TMatrix using FillSta()]
  // Flag steers which of the results should be calculated
  // 0: The paper combination of the l+j and dil results,
  //    including the calculation of the correlations per group 
  // 1: The combination using traditional correlation assignments
  //--------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =   2;
  static const Int_t NumUnc =  21;
  static const Int_t NumObs =   1;
  
  // The names
  const TString NamEst[NumEst] = {"    l+j", "    dil"};
  const TString NamUnc[NumUnc] = {"   Stat", "   Meth", "  MCgen", "    Had",
			    "   IFSR", "     UE", "     CR", "    PDF",
                            "  Wnorm", "   Wsha", "  Qnorm", "   Qsha",
			    "    JES", "   bJES", "    JER", "   JEFF",
			    "    JVF", "   btag", "    MET", "   Lept",
			    "   Pile"};
  const TString NamObs[NumObs] = {"   mtop"};
  
  // Preset according to Flag
  if(Flag == 0){
    printf("... B_EPJC_75_330: -------------------------------- \n");
    printf("... B_EPJC_75_330: The ATLAS combination, Flag = %2i \n",
	   Flag);
  }else if(Flag == 1){
    printf("... B_EPJC_75_330: -------------------------------");
    printf("---------------------------------\n");
    printf("... B_EPJC_75_330: Combination using traditional correlation");
    printf(" assignments, Flag = %2i \n", Flag);
  }else{
    printf("... B_EPJC_75_330: ------------------------ \n");
    printf("... EPJC_75_330: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Estimates 0-1 == l+j, dil
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  const Double_t XEst[LenXEst] = {
    //         0     1     2     3     4     5     6     7     8     9    10
    //      Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Wnorm  Wsha Qnorm
    //        11    12    13    14    15    16    17    18    19    20
    //      Qsha   JES  bJES   JER  JEFF   JVF  btag   MET  Lept  Pile
    172.33, 0.75, 0.11, 0.22, 0.18, 0.32, 0.15, 0.11, 0.25, 0.02, 0.29, 0.10, 
    0.05,         0.58, 0.06, 0.22, 0.12, 0.01, 0.50, 0.15, 0.04, 0.02,
    173.79, 0.54, 0.09, 0.26, 0.53, 0.47, 0.05, 0.14, 0.11, 0.01, 0.00, 0.04,
    0.01,         0.75, 0.68, 0.19, 0.07, 0.00, 0.07, 0.04, 0.13, 0.01
  };

  // Fill statistical unc of systematic uncertainties 
  static const Int_t LenSEst = NumEst * NumUnc;
  const Double_t StaUnc[LenSEst] = {
    //   0     1     2     3     4     5     6     7     8     9    10
    //Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Wnorm  Wsha Qnorm
    //  11    12    13    14    15    16    17    18    19    20
    //Qsha   JES  bJES   JER  JEFF   JVF  btag   MET  Lept  Pile
    0.00,   0.10, 0.21, 0.12, 0.06, 0.07, 0.07, 0.00, 0.00, 0.00, 0.00,
    0.00,   0.11, 0.03, 0.11, 0.00, 0.00, 0.00, 0.04, 0.00, 0.01,
    0.00,   0.07, 0.16, 0.09, 0.05, 0.05, 0.05, 0.00, 0.00, 0.00, 0.00,
    0.00,   0.08, 0.02, 0.04, 0.00, 0.00, 0.00, 0.03, 0.00, 0.00
  };

  // The correlations
  // 1) From the paper
  const Double_t RhoPap[NumUnc] = {
    //   0     1     2     3     4     5     6     7     8     9    10
    //Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Wnorm  Wsha Qnorm
    //  11    12    13    14    15    16    17    18    19    20
    //Qsha   JES  bJES   JER  JEFF   JVF  btag   MET  Lept  Pile
    0.00,   0.00, 1.00, 1.00,-1.00,-1.00,-1.00, 0.57, 1.00, 0.00, 1.00,
    0.23,  -0.23, 1.00,-1.00, 1.00,-1.00,-0.77,-0.15,-0.34, 0.00
  };

  const Double_t RhoTra[NumUnc] = {
    // 2) The traditionally used ones
    //   0     1     2     3     4     5     6     7     8     9    10
    //Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Wnorm  Wsha Qnorm
    //  11    12    13    14    15    16    17    18    19    20
    //Qsha   JES  bJES   JER  JEFF   JVF  btag   MET  Lept  Pile
    0.00,   0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.00,
    0.00,   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
  };

  // The ones used in the combinations below according to Flag
  Double_t RhoVal[NumUnc] = {0.0};
  if(Flag == 0){
    for(Int_t k = 0; k<NumUnc; k++){RhoVal[k] = RhoPap[k];}
  }else if(Flag == 1){
    for(Int_t k = 0; k<NumUnc; k++){RhoVal[k] = RhoTra[k];}
  }

  // The 21 individual JES components that make up:
  // 1) The combined uncertainty and correlation in Table 1
  // 2) The group uncertainties and correlations in Table 4
  static const Int_t NumSou =  21;
  const Double_t UncLpj[NumSou] = {
    //   0     1      2      3      4      5      6
    -0.17, +0.02, -0.01, -0.07, -0.30, +0.03, -0.01,
    //  7      8      9     10     11     12     13
    -0.01, +0.07, -0.01, -0.05, -0.02, +0.00, +0.00, 
    // 14     15     16     17     18     19     20
    +0.00, -0.11, -0.10, -0.24, -0.28, -0.22, +0.06};
  const Double_t UncDil[NumSou] = {
    //  0      1      2      3      4      5      6
    +0.01, +0.05, +0.12, +0.10, +0.22, +0.14, -0.15,
    //  7      8      9     10     11     12     13
    +0.02, +0.43, +0.45, +0.03, +0.02, +0.02, +0.00,
    // 14     15     16     17     18     19     20
    +0.03, -0.02, +0.03, -0.02, +0.03, +0.25, +0.68};
  const Int_t IndUnc[NumSou] = {
    0,         0,     0,     0,     1,     1,     1,      
    1,         1,     2,     2,     3,     3,     4,
    5,         6,     6,     7,     7,     8,     9};

  // The 10 groups of JES components from Table 4
  static const Int_t NumGrp =  10;
  Double_t SigLpj[NumGrp+1] = {0.0};
  Double_t SigDil[NumGrp+1] = {0.0};
  Double_t RhoUnc[NumGrp+1] = {0.0};
  Int_t    ActUnc = 0;

  // Declare names for output files
  const TString FilBas = "B_EPJC_75_330";
  TString FilNam;
  char Buffer[100];

  //-- Local Structures for BLUE 
  // input TMatrices
  TMatrixD* InpSta = new TMatrixD(NumEst,NumUnc,&StaUnc[0]);
  // output TMatrices
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);
  //-- End

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%6.3f";
  const TString ForRho = "%5.2f";
  const TString ForPul = ForVal;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";

  // Local structure for Blue output
  Double_t UncAct[NumObs] = {0};
  Double_t EstUnc[NumEst] = {NumEst*0};

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

  // Fill all stat uncertainties
  myBlue->FillSta(InpSta);

  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++){
      myBlue->FillCor(k, RhoVal[k]);
  }

  // Now calculate
  if(Flag == 0){

    // First show how to evaluate the correlations in Table 4
    Int_t ka = 0;
    for(Int_t k = 0; k<NumSou; k++){
      // Print components
      if(k == 0){
	printf("... B_EPJC_75_330:  The individual JES uncertainty");
	printf(" components from Table 4. \n");
	printf("... B_EPJC_75_330:  Group  Compo SigLpj SigDil RhoCom\n");
	printf("... B_EPJC_75_330: ----------------------------------\n");
      }
      if(IndUnc[k] != ka)
	printf("... B_EPJC_75_330: ----------------------------------\n");
      ka =  IndUnc[k];
      printf("... B_EPJC_75_330: %6i %6i %+6.2f %+6.2f %+6.2f \n",
	     IndUnc[k], k, UncLpj[k], UncDil[k], 
	     TMath::Sign(1.0,UncLpj[k]*UncDil[k]));
      if(k == NumSou-1){
	printf("... B_EPJC_75_330: ----------------------------------\n");
	printf("\n");
      }
      
      // Sum individual
      ActUnc = IndUnc[k];
      SigLpj[ActUnc] = SigLpj[ActUnc] + UncLpj[k]*UncLpj[k];
      SigDil[ActUnc] = SigDil[ActUnc] + UncDil[k]*UncDil[k];
      RhoUnc[ActUnc] = RhoUnc[ActUnc] + 
	TMath::Sign(UncLpj[k]*UncDil[k],UncLpj[k]*UncDil[k]);
      
      // Sum for total but bjes == IndUnc[k] == 9
      if(ActUnc != 9){
	SigLpj[NumGrp] = SigLpj[NumGrp] + UncLpj[k]*UncLpj[k];
	SigDil[NumGrp] = SigDil[NumGrp] + UncDil[k]*UncDil[k];
	RhoUnc[NumGrp] = RhoUnc[NumGrp] + 
	  TMath::Sign(UncLpj[k]*UncDil[k],UncLpj[k]*UncDil[k]);
      }
    }

    // Print results per source
    printf("... B_EPJC_75_330:  The groups of JES uncertainty");
    printf(" components from Table 4 \n");
    printf("... B_EPJC_75_330:  For (SigLpj*SigDil)== 0  RhoMea == -2.00\n");
    printf("... B_EPJC_75_330:  The RhoMea values suffer most from");
    printf(" rounding of the inputs.\n");
    printf("... B_EPJC_75_330:  Group SigLpj  SigDil  RhoMea\n");
    printf("... B_EPJC_75_330: -----------------------------\n");
    for(Int_t k = 0; k<NumGrp; k++){
      SigLpj[k] = TMath::Sqrt(SigLpj[k]);
      SigDil[k] = TMath::Sqrt(SigDil[k]);
      if(SigLpj[k]*SigDil[k]>0.0){RhoUnc[k] = RhoUnc[k] / (SigLpj[k]*SigDil[k]);
      }else{RhoUnc[k] = -2.;
      }
      printf("... B_EPJC_75_330: %6i %6.2f  %6.2f  %+6.2f \n", k,
	     SigLpj[k], SigDil[k], RhoUnc[k]);
    }
    
    // Get final result for groups and print it
    SigLpj[NumGrp] = TMath::Sqrt(SigLpj[NumGrp]);
    SigDil[NumGrp] = TMath::Sqrt(SigDil[NumGrp]);
    if(SigLpj[NumGrp]*SigDil[NumGrp]>0.0){
      RhoUnc[NumGrp] = RhoUnc[NumGrp] / (SigLpj[NumGrp]*SigDil[NumGrp]);
    }else{RhoUnc[NumGrp] = -2.;
    }
    printf("... B_EPJC_75_330: -----------------------------\n");
    printf("... B_EPJC_75_330:  Total %6.2f  %6.2f  %+6.2f (without bJES)\n",
	   SigLpj[NumGrp], SigDil[NumGrp], RhoUnc[NumGrp]);
  }

  // Perform the combination
  myBlue->FixInp();
  myBlue->Solve();
  
  // Do some print out
  if(Flag == 0){
    printf("... B_EPJC_75_330: \n");
    printf("... B_EPJC_75_330: For the BLUE weights, the compatibility");
    printf(" and the correlation see Section 8.1. \n");     
    printf("... B_EPJC_75_330: \n");
  }
  myBlue->PrintEst();
  myBlue->PrintCompatEst();

  printf("... B_EPJC_75_330: The estimator correlation \n");
  myBlue->GetRho(LocRho);
  myBlue->PrintMatrix(LocRho, ForRho);
  myBlue->PrintPull();
  if(Flag == 0){
    printf("... B_EPJC_75_330: \n");
    printf("... B_EPJC_75_330: For the breakdown of uncertainties of");
    printf(" the result see Table 3. \n");
    printf("... B_EPJC_75_330: \n");
  }
  myBlue->PrintResult();
  
  // Do the numerical difference for the combination scenarios
  myBlue->GetEstUnc(EstUnc);
  myBlue->GetUncert(UncAct);
  printf("... B_EPJC_75_330: \n");
  if(Flag == 0){printf("... B_EPJC_75_330: Evaluated correlations:");
  }else{printf("... B_EPJC_75_330: Assigned correlations:");
  }
  printf(" -- the gain in precision is \n");
  printf("... B_EPJC_75_330: from sigma_1 = %4.2f", EstUnc[0]);
  printf(" to sigma_x = %4.2f,\n", UncAct[0]);
  printf("... B_EPJC_75_330: which corresponds to %2.0f%%\n",
	 100*(EstUnc[0]-UncAct[0])/ EstUnc[0]);

  // Get the output into tex files and plots
  printf("... B_EPJC_75_330: \n");
  sprintf(Buffer,"%1i", Flag);
  FilNam = &Buffer[0];       
  FilNam = FilBas + "_" + FilNam;
  myBlue->LatexResult(FilNam);

  // Display result
  myBlue->DisplayResult(0,FilNam);

  // Show the dependence on rho
  myBlue->DisplayPair(0,1,FilNam);

  // Exercise SolveScaSta()
  if(Flag == 0){
    myBlue->ReleaseInp();
    myBlue->FixInp();
    myBlue->SolveScaSta(2);
    myBlue->PrintScaSta(FilNam);
  }    
  
  // Delete object and matrices
  delete myBlue;  myBlue = NULL;  

  InpSta->Delete(); InpSta = NULL;
  LocRho->Delete(); LocRho = NULL;

  // Return
  return;
}
