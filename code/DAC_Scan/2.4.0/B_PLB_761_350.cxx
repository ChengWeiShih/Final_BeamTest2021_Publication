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

//---------------------------------------------------------------------------
// Function prototype for utility
//---------------------------------------------------------------------------
void FilMat(Double_t R01, Double_t R02, Double_t R12, Double_t *M);

void B_PLB_761_350(Int_t Flag = 0){
 
  //---------------------------------------------------------------------------
  //  [Shows how to use SetFormat()]
  //  [Shows how to fill names with FillNamXXX() with XXX = Est, Unc, Obs]
  //  [Shows how to fill stat. precions of systematics with FillSta()]
  // Flag == 0 ATLAS combination - one  observable -   evaluated correlations
  //           [Shows how to use CorrelPair() and DiplayPair()]
  // Flag == 1 ATLAS combination - one  observable - traditional correlations
  // Flag == 2 ATLAS combination - two observables -   evaluated correlations
  //           [Shows how to use SolveScaSta(), GetSta() and GetStaRes()]
  // Flag == 3 ATLAS combination - two observables - traditional correlations
  //---------------------------------------------------------------------------
 
  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  3;
  static const Int_t NumUnc = 20;
  static const Int_t NumObs =  2;
  Int_t IndCor = 0;
  Int_t ActObs = 1;
  Int_t IWhichObs[NumEst] = {NumEst*0};

  // The names
  const TString NamEst[NumEst] = {"l+jets (7 TeV)", "dil (7 TeV)", 
				  "dil (8 TeV)"};
  const TString NamUnc[NumUnc] = {
    //     0          1          2          3          4          5
    "   Stat", "   Meth", "  MCgen", "    Had", "   IFSR", "     UE", 
    //     6          7          8          9         10         11
    "     CR", "    PDF", "  Bnorm", "   Wsha", "   Qsha", "    JES", 
    //    12         13         14         15         16         17
    "   bJES", "    JER", "   JEFF", "    JVF", "   btag", "   Lept", 
    //    18         19
    "    MET", "   Pile"};
  TString NamObs[NumObs] = {"m_{top}", "m_{top} (dil)"};

  // Preset according to Flag
  if(Flag == 0 || Flag == 1){
    printf("... B_PLB_761_350: --------------------------------- \n");
    printf("... B_PLB_761_350: The ATLAS combinations, Flag = %2i \n",
	   Flag);
    if(Flag == 1)IndCor = 1;
  }else if(Flag == 2 || Flag == 3){
    printf("... B_PLB_761_350: -------------------------------");
    printf("--------------------\n");
    printf("... B_PLB_761_350: The ATLAS combination for two observables");
    printf(" Flag = %2i \n", Flag);
    ActObs = 2;
    NamObs[0] = "m_{top} (l+j)";
    IWhichObs[0] = 0;
    IWhichObs[1] = 1;
    IWhichObs[2] = IWhichObs[1];
    if(Flag == 3)IndCor = 1;
  }else{
    printf("... B_PLB_761_350: ------------------------ \n");
    printf("... B_PLB_761_350: Not implemented Flag = %2i \n",Flag);
    return;
  }
  if(Flag == 1 || Flag == 3){
    printf("... B_PLB_761_350: The traditional correlation assumptions\n");
  }

  // Fill estimates (0, 1, 2 == l+jets7, dil7, dil8) and uncertainties 
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  const Double_t XEst[LenXEst] = {
    //         0     1     2     3     4     5     6     7     8     9    10
    //      Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Bnorm  Wsha  Qsha
    // 11     12    13    14    15    16    17    18    19
    //JES   bJES   JER  JEFF   JVF  btag  Lept   MET  Pile
    172.33, 0.75, 0.11, 0.22, 0.18, 0.32, 0.15, 0.11, 0.25, 0.10, 0.29, 0.05,
    0.58,   0.06, 0.22, 0.12, 0.01, 0.50, 0.04, 0.15, 0.02,
    173.79, 0.54, 0.09, 0.26, 0.53, 0.47, 0.05, 0.14, 0.11, 0.04, 0.00, 0.01,
    0.75,   0.68, 0.19, 0.07, 0.00, 0.07, 0.13, 0.04, 0.01,
    172.99, 0.41, 0.05, 0.09, 0.22, 0.23, 0.10, 0.03, 0.05, 0.03, 0.00, 0.08,
    0.54,   0.30, 0.09, 0.01, 0.02, 0.03, 0.14, 0.01, 0.05
  };

  // Fill statistical uncertainty of systematic uncertainties 
  static const Int_t LenSEst = NumEst * NumUnc;
  const Double_t StaUnc[LenSEst] = {
    //  0     1     2     3     4     5     6     7     8     9    10
    //Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Bnorm  Wsha  Qsha
    //  11    12    13    14    15    16    17    18    19
    // JES  bJES   JER  JEFF   JVF  btag  Lept   MET  Pile
    0.00,   0.10, 0.21, 0.12, 0.06, 0.07, 0.07, 0.00,  0.00, 0.00, 0.00,
    0.11,   0.03, 0.11, 0.00, 0.00, 0.00, 0.00, 0.04,  0.01,
    0.00,   0.07, 0.16, 0.09, 0.05, 0.05, 0.05, 0.00,  0.00, 0.00, 0.00,
    0.08,   0.02, 0.04, 0.00, 0.00, 0.00, 0.00, 0.03,  0.00,
    0.00,   0.07, 0.15, 0.09, 0.07, 0.14, 0.14, 0.00,  0.00, 0.00, 0.00,
    0.04,   0.01, 0.05, 0.00, 0.00, 0.02, 0.01, 0.01,  0.01
  };

  // 1) Fill actual correlations:
  //  (01, 02, 12 == l+jets7-dil7,l+jets7-dil8, dil7-dil8)
  static const Int_t LenRho = NumEst * NumUnc;  
  const Double_t RhoMij[LenRho] = {
    //  0     1     2     3     4     5     6     7     8     9    10
    //Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Bnorm  Wsha  Qsha
    //  11    12    13    14    15    16    17    18    19
    // JES  bJES   JER  JEFF   JVF  btag  Lept   MET  Pile
    +0.00,  0.00, 1.00, 1.00,-1.00,-1.00,-1.00, 0.57, 1.00, 0.00, 0.23,
    -0.23,  1.00,-1.00, 1.00,-1.00,-0.77,-0.34,-0.15, 0.00,
    +0.00,  0.00, 1.00, 1.00,-1.00,-1.00,-1.00,-0.29, 0.23, 0.00, 0.20, 
    +0.06,  1.00, 0.00, 1.00, 1.00, 0.00,-0.52, 0.25, 0.00,
    +0.00,  0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.03, 0.23, 0.00,-0.08,
    +0.35,  1.00, 0.00, 1.00,-1.00, 0.00, 0.96,-0.24, 0.00
  };

  // 2) fill the traditionally used ones
  const Double_t RhoTra[LenRho] = {
    //  0     1     2     3     4     5     6     7     8     9    10
    //Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Bnorm  Wsha  Qsha
    //  11    12    13    14    15    16    17    18    19
    // JES  bJES   JER  JEFF   JVF  btag  Lept   MET  Pile
      0.00, 0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.50, 0.00, 0.00, 
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      0.00, 0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.50, 0.00, 0.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      0.00, 0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.50, 0.00, 0.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
  };

  // Copy the ones used in the combinations below according to Flag
  Double_t RhoVal[LenRho] = {0.0};
  if(Flag == 0 || Flag == 2){
    for(Int_t k = 0; k<LenRho; k++)RhoVal[k] = RhoMij[k];
  }else{
    for(Int_t k = 0; k<LenRho; k++)RhoVal[k] = RhoTra[k];
  }

  // The dummy correlation matrix to be filled dynamically
  static const Int_t LenCor = NumEst * NumEst;
  Double_t CorSou[LenCor] = {0};

  // The default file name
  const TString FilBas = "B_PLB_761_350";
  TString FilNam;
  char Buffer[100];
  
  //-- Local Structures for BLUE
  // Input TMatrices
  TMatrixD* InpEst = new TMatrixD(NumEst,NumUnc+1,&XEst[0]);
  TMatrixD* InpSta = new TMatrixD(NumEst,NumUnc,&StaUnc[0]);
  //-- End

  //-- Local Structures for BLUE output
  // TMatrices
  TMatrixD* LocWei = new TMatrixD(NumEst,NumObs);
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocSta    = new TMatrixD(NumEst,NumUnc);
  TMatrixD* LocStaRes = new TMatrixD(NumObs,6);
  // Arrays
  Double_t  ArrSta[NumEst*NumUnc] = {0};
  Double_t  ArrStaRes[NumObs*6] = {0};
  //-- End

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = ForUnc;
  const TString ForRho = ForUnc;
  const TString ForPul = ForRho;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";
 
  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, ActObs, &IWhichObs[0]);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

  // Retrieve some parameters
  printf("... B_PLB_761_350: Object constructed for NumEst = %i",
	 myBlue->GetNumEst());
  printf(" NumUnc = %i and NumObs = %i \n",
	 myBlue->GetNumUnc(), myBlue->GetNumObs());

  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);

  // Fill all estimates
  myBlue->FillEst(InpEst);

  // Fill all statistical uncertainties of the systematic uncertainties
  myBlue->FillSta(InpSta);

  // Fill correlations, fill dummy and pass on to Blue
  for(Int_t k = 0; k<NumUnc; k++){
    FilMat(RhoVal[k], RhoVal[k+NumUnc], RhoVal[k+2*NumUnc], &CorSou[0]);
    myBlue->FillCor(k,&CorSou[0]);
  }

  // Solve and inspect result
  const Double_t xmin = 166.;
  const Double_t xmax = 180.;
  const Double_t smin = 0.0;
  const Double_t smax = 2.5;
  Int_t nmax = 3;
  if(Flag >= 2)nmax = 1;
  // Loop over sets of selected estimates
  // n = 0/1/2 = all/dil/old
  for(Int_t n = 0; n < nmax; n++){
    // Set file name and select estimates
    if(n == 0){
      sprintf(Buffer,"all_NObs_%1i_Cor_%1i", ActObs, IndCor);
    }else if(n == 1){
      sprintf(Buffer,"dil_Cor_%1i", IndCor);
      myBlue->SetInActiveEst(0);
    }else if(n == 2){
      sprintf(Buffer,"old_Cor_%1i", IndCor);
      myBlue->SetInActiveEst(2);      
    }
    FilNam = &Buffer[0];

    // Fix and print input
    myBlue->FixInp();
    myBlue->PrintEst();

    // Look at pairs compatibility and correlatioons, when all are enabled
    if(n == 0 && Flag == 0){
      myBlue->DisplayPair(0,1,FilBas,xmin,xmax,smin,smax,ForVal,ForUnc,ForRho);
      myBlue->DisplayPair(0,2,FilBas,xmin,xmax,smin,smax,ForVal,ForUnc,ForRho);
      myBlue->DisplayPair(1,2,FilBas,xmin,xmax,smin,smax,ForVal,ForUnc,ForRho);

      printf("... B_PLB_761_350: Compared to the original publications,");
      printf(" the figures produced by CorrelPair\n");
      printf("... B_PLB_761_350: only show subset of all uncertainties,");
      printf(" and has no information on the quadrant\n");
      printf("... B_PLB_761_350: the point was located in.");
      printf(" As a result, quadrants (Q1,Q3) are shown \n");
      printf("... B_PLB_761_350: in Q1, and similarly (Q2,Q4) in Q2.");
      printf(" In addition, only combined sources for \n");
      printf("... B_PLB_761_350: which the estimators have rho=+-1");
      printf(" can be shown, see Table 2. \n");
      myBlue->CorrelPair(0,1,FilBas,-0.85,0.85,-0.15,0.85,ForUnc);
      myBlue->CorrelPair(0,2,FilBas,-0.85,0.85,-0.15,0.85,ForUnc);
      myBlue->CorrelPair(1,2,FilBas,-0.85,0.85,-0.15,0.85,ForUnc);

      // The compatibility
      myBlue->PrintCompatEst();

      // The correlations
      printf("... B_PLB_761_350: The correlations, see Table 2 \n");
      printf("... B_PLB_761_350:");
      for(Int_t n = 0; n < NumEst; n++)printf(" %s ", NamEst[n].Data());
      printf("\n");
      myBlue->GetRho(LocRho);
      myBlue->PrintMatrix(LocRho,"   %+4.2f");
    }

    // Solve and print result
    myBlue->Solve();
    myBlue->PrintResult();

    // The weights
    printf("... B_PLB_761_350: The weights, see Section 7 \n");
    myBlue->GetWeight(LocWei);
    LocWei->operator*=(100);
    myBlue->PrintMatrix(LocWei,myBlue->GetActEst(),
			myBlue->GetActObs(),"%+6.1f%%");
    
    // Get the output into files
    FilNam = FilBas + "_" + FilNam; 
    myBlue->LatexResult(FilNam);
    for(Int_t n = 0; n < myBlue->GetActObs(); n++){
      myBlue->DisplayResult(n,FilNam);
    }
    
    // Show how to use SolveScaSta()
    if(n == 0 &&  Flag == 2){ 
      for(Int_t IScSta = 0; IScSta<3; IScSta++){
	myBlue->ReleaseInp();
	sprintf(Buffer,"_%1i", Flag);
	FilNam = &Buffer[0];
	FilNam = FilBas + FilNam; 
	myBlue->FixInp();
	myBlue->SolveScaSta(IScSta);
	myBlue->PrintScaSta(FilNam, 172.00, 174.00, 0.70, 1.90);

	// Show how to use GetSta and GetStaRes
	if(IScSta == 0){
	  myBlue->GetSta(LocSta);
	  printf("... B_PLB_761_350: LocSta \n");
	  myBlue->PrintMatrix(LocSta, "%4.2f");

	  myBlue->GetSta(ArrSta);
	  printf("... B_PLB_761_350: ArrSta \n");
	  myBlue->PrintDouble(ArrSta,NumEst,NumUnc,"%4.2f");
	  
	  myBlue->GetStaRes(LocStaRes);
	  printf("... B_PLB_761_350: LocStaRes \n");
	  myBlue->PrintMatrix(LocStaRes,ActObs,6, "%4.2f");

	  myBlue->GetStaRes(ArrStaRes);
	  printf("... B_PLB_761_350: ArrStaRes \n");
	  myBlue->PrintDouble(ArrStaRes,ActObs,6,"%4.2f");
	}
      }
    }
    // Reset for next n
    myBlue->ResetInp();
  }

  // Delete Object
  delete myBlue; myBlue = NULL;
  InpEst->Delete(); InpEst = NULL;
  InpSta->Delete(); InpSta = NULL;

  LocWei->Delete(); LocWei = NULL;
  LocRho->Delete(); LocRho = NULL;
  LocSta->Delete(); LocSta = NULL;
  LocStaRes->Delete(); LocStaRes = NULL;

  return;
}

//------------------------------------------------------------------------
// Utility to fill a correlation matrix for any source of uncertainty
//------------------------------------------------------------------------

void FilMat(Double_t R01, Double_t R02, Double_t R12, Double_t *M){
  M[0] =   1; M[1] = R01; M[2] = R02;
  M[3] = R01; M[4] =   1; M[5] = R12;
  M[6] = R02; M[7] = R12; M[8] =   1;
  return;
};
