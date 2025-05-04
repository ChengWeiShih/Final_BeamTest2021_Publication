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
// The combination
//---------------------------------------------------------------------------
void B_EPJC_79_290(Int_t Flag = 0){
  
  //---------------------------------------------------------------------------
  // [Shows how to de-activate estimates with SetActiveEst()]
  // [Shows how to use SolveAccImp() and DisplayAccImp()]
  // Flag == 0 ATLAS combination - 7 TeV data
  // Flag == 1 ATLAS combination - 8 TeV data
  // Flag == 2 ATLAS combination - three  observables
  //         [Shows how to use LatexResult() with different formats]
  // Flag == 3 ATLAS combination - three significant estimates
  //           [Shows how to use SolveScaSta() and PrintScaSta()]
  // Flag == 4 ATLAS combination - full combination
  //           [Shows how to use CompatEst()]
  //---------------------------------------------------------------------------
 
  // The full list of estimates, uncertainties and observables
  static const Int_t NumEst =   6;
  static const Int_t NumUnc =  23;
  static const Int_t MaxObs =   3;
  Int_t NumObs = 1;
  
  // Array of names of estimates
  const TString NamEst[NumEst] = {
    "m_{top}^{dilepton}(7TeV)", 
    "m_{top}^{l+jets}(7TeV)", 
    "m_{top}^{all jets}(7TeV)",
    "m_{top}^{dilepton}(8TeV)", 
    "m_{top}^{l+jets}(8TeV)", 
    "m_{top}^{all jets}(8TeV)"};

  // Array of names of uncertainties
  const TString NamUnc[NumUnc] = {
    //     0          1          2          3          4          5
    "   Stat", "   Meth", "  MCgen", "    Had", "   IFSR", "     UE", 
    //     6          7          8          9         10         11
    "     CR", "    PDF", "  Bnorm", "   Wsha", "   FSha", "   BAJE",
    //    12         13         14         15         16         17
    "    JES", "   bJES", "    JER", "   JEFF", "    JVF", "   btag",
    //    18         19         20         21         22
    "   Lept", "    MET", "   Pile", " AJTrig", "   FaFu"};

  // Array of names of observables
  TString NamObs[MaxObs] = {"m_{top}",
			    "m_{top}^{l+jets}", 
			    "m_{top}^{all jets}"};

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {0};
  
  // The steering flag
  if(Flag == 0){    
    printf("... B_EPJC_79_290: ------------------------------------- \n");
    printf("... B_EPJC_79_290:     The 7 TeV combination, Flag = %2i \n",Flag);
    NamObs[0] = "m_{top}^{7TeV}";
  }else if(Flag == 1){
    printf("... B_EPJC_79_290: ------------------------------------- \n");
    printf("... B_EPJC_79_290:     The 8 TeV combination, Flag = %2i \n",Flag);
    NamObs[0] = "m_{top}^{8TeV}";
  }else if(Flag == 2){
    printf("... B_EPJC_79_290: ------------------------------------- \n");
    printf("... B_EPJC_79_290: The 3-channel combination, Flag = %2i \n",Flag);
    NamObs[0] = "m_{top}^{dilepton}";
    NumObs = 3;
    IWhichObs[0] = 0;
    IWhichObs[1] = 1;
    IWhichObs[2] = 2;
    IWhichObs[3] = 0;
    IWhichObs[4] = 1;
    IWhichObs[5] = 2;
  }else if(Flag == 3){
    printf("... B_EPJC_79_290: --------------------------------------- \n");
    printf("... B_EPJC_79_290:   The mtop(3) combination, Flag = %2i \n",Flag);
    NamObs[0] = "m_{top}^{(3)}";
  }else if(Flag == 4){
    printf("... B_EPJC_79_290: ------------------------------------- \n");
    printf("... B_EPJC_79_290:      The full combination, Flag = %2i \n",Flag);
  }else{
    printf("... B_EPJC_79_290: ------------------------------------- \n");
    printf("... B_EPJC_79_290:            Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  // Fill estimates and uncertainties 
  // 0-5 == dil7, l+j7, aje7, dil8, l+j8, aje8
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  const Double_t XEst[LenXEst] = {
    //         0     1     2     3     4     5     6     7     8     9    10
    //      Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Bnorm  Wsha  FSha
    // 11     12    13    14    15    16    17    18    19   20     21    22
    //BAJE   JES  bJES   JER  JEFF   JVF  btag  Lept   MET Pile AJTrig  FaFu
    173.79, 0.54, 0.09, 0.26, 0.53, 0.47, 0.05, 0.14, 0.10, 0.04, 0.00, 0.01,
    0.00,   0.76, 0.68, 0.19, 0.07, 0.00, 0.07, 0.13, 0.04, 0.01, 0.00, 0.00,
    172.33, 0.75, 0.11, 0.22, 0.18, 0.32, 0.15, 0.11, 0.25, 0.10, 0.29, 0.05,
    0.00,   0.58, 0.06, 0.22, 0.12, 0.01, 0.50, 0.04, 0.15, 0.02, 0.00, 0.00,
    175.06, 1.35, 0.42, 0.30, 0.50, 0.22, 0.08, 0.22, 0.09, 0.00, 0.00, 0.00,
    0.35,   0.50, 0.62, 0.01, 0.01, 0.01, 0.16, 0.00, 0.02, 0.02, 0.01, 0.24,
    172.99, 0.41, 0.05, 0.09, 0.22, 0.23, 0.10, 0.03, 0.05, 0.03, 0.00, 0.07,
    0.00,   0.54, 0.30, 0.09, 0.01, 0.02, 0.04, 0.14, 0.01, 0.05, 0.00, 0.00,
    172.08, 0.39, 0.13, 0.16, 0.15, 0.08, 0.08, 0.19, 0.09, 0.08, 0.11, 0.00,
    0.00,   0.54, 0.03, 0.20, 0.02, 0.09, 0.38, 0.16, 0.05, 0.15, 0.00, 0.00,
    173.72, 0.55, 0.11, 0.18, 0.64, 0.10, 0.12, 0.12, 0.09, 0.00, 0.00, 0.00,
    0.17,   0.60, 0.34, 0.10, 0.00, 0.03, 0.10, 0.01, 0.01, 0.01, 0.08, 0.00
  };

  // Fill statistical uncertainty of systematic uncertainties 
  static const Int_t LenSEst = NumEst * NumUnc;
  const Double_t SUnc[LenSEst] = {
    //   0     1     2     3     4     5     6     7     8     9    10
    //Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Bnorm  Wsha  FSha
    // 11     12    13    14    15    16    17    18    19   20     21    22
    //BAJE   JES  bJES   JER  JEFF   JVF  btag  Lept   MET Pile AJTrig  FaFu
    0.00,   0.07, 0.16, 0.09, 0.05, 0.05, 0.05, 0.00, 0.00, 0.00, 0.00,
    0.00,   0.09, 0.02, 0.04, 0.00, 0.00, 0.00, 0.00, 0.03, 0.00, 0.00, 0.00,
    0.00,   0.10, 0.21, 0.12, 0.06, 0.07, 0.07, 0.00, 0.00, 0.00, 0.00,
    0.00,   0.11, 0.03, 0.11, 0.00, 0.00, 0.00, 0.00, 0.04, 0.01, 0.00, 0.00,
    0.00,   0.01, 0.30, 0.15, 0.11, 0.10, 0.10, 0.00, 0.00, 0.00, 0.00,
    0.21,   0.05, 0.05, 0.08, 0.01, 0.01, 0.00, 0.00, 0.05, 0.00, 0.01, 0.18,
    0.00,   0.07, 0.15, 0.09, 0.07, 0.14, 0.14, 0.00, 0.00, 0.00, 0.00,
    0.00,   0.04, 0.01, 0.05, 0.00, 0.00, 0.02, 0.01, 0.01, 0.01, 0.00, 0.00,
    0.00,   0.11, 0.17, 0.10, 0.11, 0.15, 0.15, 0.00, 0.00, 0.00, 0.00,
    0.00,   0.02, 0.01, 0.04, 0.01, 0.01, 0.00, 0.01, 0.01, 0.01, 0.00, 0.00, 
    0.00,   0.00, 0.21, 0.15, 0.28, 0.16, 0.16, 0.00, 0.00, 0.00, 0.00,
    0.00,   0.03, 0.02, 0.04, 0.00, 0.01, 0.00, 0.00, 0.01, 0.00, 0.01, 0.00
  };

  // Fill actual correlations:
  // 01, 02, 03, 04, 05 // 12, 13, 14, 15 // 23, 24, 25 // 34, 35 // 45
  static const Int_t NumPai = 0.5 * NumEst * (NumEst-1);
  static const Int_t LenRho = NumPai * NumUnc;
  const Double_t RhoVal[LenRho] = {
    //   0     1     2     3     4     5     6     7     8     9    10
    //Stat  Meth MCgen   Had  IFSR    UE    CR   PDF Bnorm  Wsha  FSha
    // 11     12    13    14    15    16    17    18    19   20     21    22
    //BAJE   JES  bJES   JER  JEFF   JVF  btag  Lept   MET Pile AJTrig  FaFu
    //-01
    +0.00,  0.00, 1.00, 1.00,-1.00,-1.00,-1.00, 0.53, 1.00, 0.00, 0.20, 
    +0.00, -0.24, 1.00,-1.00, 1.00,-1.00,-0.80,-0.35, 0.00, 0.00, 0.00, 0.00,
    //-02
    +0.00,  0.00,-1.00, 1.00, 1.00, 1.00, 1.00, 0.22, 0.00, 0.00, 0.00, 
    +0.00,  0.86, 1.00, 1.00, 1.00, 1.00,-0.03, 0.00,-0.26, 0.00, 0.00, 0.00,
    //-03
    +0.00,  0.00, 1.00, 1.00, 1.00, 1.00, 1.00,-0.02, 0.31, 0.00, 0.0,
    +0.00,  0.36, 1.00, 0.00, 1.00,-1.00, 0.00, 0.93,-0.26, 0.00, 0.00, 0.00,
    //-04
    +0.00,  0.00, 1.00,-1.00, 1.00, 1.00,-1.00, 0.72,-0.77, 0.00, 0.00, 
    +0.00,  0.18, 1.00, 0.00, 1.00,-1.00, 0.00,-0.08,-0.12, 0.00, 0.00, 0.00,
    //-05
    +0.00,  0.00, 1.00,-1.00, 1.00,-1.00,-1.00,-0.61, 0.00, 0.00, 0.00, 
    +0.00,  0.36, 1.00, 0.00, 0.00,-1.00, 0.00, 0.42, 0.04, 0.00, 0.00, 0.00,
    //-12
    +0.00,  0.00,-1.00, 1.00,-1.00,-1.00,-1.00,-0.36, 0.00, 0.00, 0.00, 
    +0.00,  0.10, 1.00,-1.00, 1.00,-1.00, 0.00, 0.00, 0.84, 0.00, 0.00, 0.00,
    //-13
    +0.00,  0.00, 1.00, 1.00,-1.00,-1.00,-1.00,-0.32, 0.31, 0.00, 0.27, 
    +0.00,  0.04, 1.00, 0.00, 1.00, 1.00, 0.00,-0.51, 0.26, 0.00, 0.00, 0.00,
    //-14
    +0.00,  0.00, 1.00,-1.00,-1.00,-1.00, 1.00, 0.72,-0.74, 0.00, 0.00, 
    +0.00, -0.29, 1.00, 0.00, 1.00, 1.00, 0.00,-0.17, 0.22, 0.00, 0.00, 0.00,
    //-15
    +0.00,  0.00, 1.00,-1.00,-1.00, 1.00, 1.00,-0.81, 0.00, 0.00, 0.00, 
    +0.00,  0.13, 1.00, 0.00, 0.00, 1.00, 0.00, 0.02, 0.16, 0.00, 0.00, 0.00,
    //-23
    +0.00,  0.00,-1.00, 1.00, 1.00, 1.00, 1.00, 0.41, 0.00, 0.00, 0.00, 
    +0.00,  0.41, 1.00, 0.00, 1.00,-1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    //-24
    +0.00,  0.00,-1.00,-1.00, 1.00, 1.00,-1.00,-0.05, 0.00, 0.00, 0.00, 
    +0.00,  0.09, 1.00, 0.00, 1.00,-1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    //-25
    +0.00,  0.00,-1.00,-1.00, 1.00,-1.00,-1.00, 0.27, 0.00, 0.00, 0.00, 
    +1.00,  0.42, 1.00, 0.00, 0.00,-1.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00,
    //-34
    +0.00,  0.00, 1.00,-1.00, 1.00, 1.00,-1.00,-0.48,-0.06, 0.00, 0.00, 
    +0.00, -0.54, 1.00, 0.22, 1.00, 1.00,-0.23, 0.11, 0.97, 0.00, 0.00, 0.00,
    //-35
    +0.00,  0.00, 1.00,-1.00, 1.00,-1.00,-1.00, 0.40, 0.00, 0.00, 0.00, 
    +0.00,  0.98, 1.00,-0.07, 0.00, 1.00, 1.00, 0.28, 0.86, 0.00, 0.00, 0.00,
    //-45
    +0.00,  0.00, 1.00, 1.00, 1.00,-1.00, 1.00,-0.76, 0.00, 0.00, 0.00, 
    +0.00, -0.57, 1.00,-0.17, 0.00, 1.00, 1.00,-0.36, 0.96, 0.00, 0.00, 0.00
  };

  //-- Local Structures for BLUE
  // input TMatrices
  TMatrixD* InpEst = new TMatrixD(NumEst,NumUnc+1,&XEst[0]);
  TMatrixD* InpSta = new TMatrixD(NumEst,NumUnc,  &SUnc[0]);
  TMatrixD* InpRho = new TMatrixD(NumPai,NumUnc,&RhoVal[0]);
  // output TMatrices
  TMatrixD* LocRhoRes  = new TMatrixD(NumObs,NumObs);
  //-- End

  // The dummy correlation array to be filled dynamically
  static const Int_t LenCor = 0.5 * NumEst * (NumEst-1);
  Double_t CorSou[LenCor] = {0};

  // The default file names
  char Buffer[100]; 
  sprintf(Buffer,"_%1i", Flag);
  const TString FilBas = "B_EPJC_79_290";
  TString FilNam = FilBas + &Buffer[0];
  const TString WebNam =
    "//atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2017-03/";

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = ForUnc;
  const TString ForRho = "%+4.2f";
  const TString ForPul = ForUnc;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";

  // Construct object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);
  myBlue->PrintStatus();
  
  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);
  
  // Fill estimates
  myBlue->FillEst(InpEst);
       
  // Fill statistical precision
  myBlue->FillSta(InpSta);

  // Fill correlations. First fill dummy array, then pass it on to Blue
  Int_t ind = 0;
  for(Int_t k = 0; k<NumUnc; k++){
    if(k == 0){myBlue->FillCor(k, 0.0);
    }else{
      ind = 0;
      for(Int_t i = 0; i<NumPai; i++){
	CorSou[ind] = InpRho->operator()(i,k);
	ind = ind + 1;
      }
      myBlue->FillCor(-k, &CorSou[0]);
    }
  }

  // Deactivate all estimates then re-activate according to Flag
  for(Int_t i = 0; i<NumEst; i++)myBlue->SetInActiveEst(i);
  if(Flag == 0){
    // 7 TeV
    myBlue->SetActiveEst(0);
    myBlue->SetActiveEst(1);
    myBlue->SetActiveEst(2);
  }else if(Flag == 1){
    // 8 TeV
    myBlue->SetActiveEst(3);
    myBlue->SetActiveEst(4);
    myBlue->SetActiveEst(5);
  }else if(Flag == 2){
    // 3 Obs
    for(Int_t i = 0; i<NumEst; i++)myBlue->SetActiveEst(i);
  }else if(Flag == 3){
    // mtop(3)
    myBlue->SetActiveEst(1);
    myBlue->SetActiveEst(3);
    myBlue->SetActiveEst(4);
  }else if(Flag == 4){
    // mtop
    for(Int_t i = 0; i<NumEst; i++)myBlue->SetActiveEst(i);
  }

  // Fix and solve depending on Flag
  myBlue->FixInp();
  if(Flag == 0){
    // Solve According to Importance == Fig13a
    printf("... B_EPJC_79_290: This is Figure 13a \n");
    myBlue->SolveAccImp(0.001);
    myBlue->DisplayAccImp(0,FilNam);
  }else if(Flag == 1){
    // Solve According to Importance == Fig13b
    printf("... B_EPJC_79_290: This is Figure 13b \n");
    myBlue->SolveAccImp(0.001);
    myBlue->DisplayAccImp(0,FilNam);
  }else if(Flag == 2){
    // Solve According to Importance == Fig13c-e
    printf("... B_EPJC_79_290: These are Figures 13c-e \n");
    myBlue->SolveAccImp(0.001);
    for(Int_t n = 0; n<NumObs; n++)myBlue->DisplayAccImp(n,FilNam);
  }else if(Flag == 3){
    // Compatibity of estimates
    myBlue->PrintCompatEst();

    printf("... B_EPJC_79_290: These are the auxiliary figures 3a,b from");
    printf(" %s \n", WebNam.Data());
    myBlue->SolveScaSta();
    myBlue->PrintScaSta(FilNam, 172.38, 172.98, 0.34, 0.78);
  }else if(Flag == 4){
    // Compatibity of estimates
    myBlue->PrintCompatEst();

    // Fig7 a/ b / c / d
    printf("... B_EPJC_79_290: This is Figure 7, but only for the single \n");
    printf("... B_EPJC_79_290: component uncertainty sources with rho=+-1,\n");
    printf("... B_EPJC_79_290: and using two quadrants only.\n");

    printf("... B_EPJC_79_290: The information about the sub-components\n");
    printf("... B_EPJC_79_290: and the signs of the shifts is lost.\n");
    myBlue->CorrelPair(3,4,FilBas,-0.85,0.85,-0.10,0.85,ForUnc);
    myBlue->CorrelPair(3,5,FilBas,-0.85,0.85,-0.10,0.85,ForUnc);
    myBlue->CorrelPair(4,5,FilBas,-0.85,0.85,-0.10,0.85,ForUnc);
    myBlue->CorrelPair(4,1,FilBas,-0.85,0.85,-0.10,0.85,ForUnc);
    
    // Fig8 a,b / c,d  / e,f
    printf("... B_EPJC_79_290: This is Figure 8 \n");
    myBlue->DisplayPair(3,4,FilBas);
    myBlue->DisplayPair(3,5,FilBas);
    myBlue->DisplayPair(4,5,FilBas);

    // Solve According to Importance == Fig10b
    printf("... B_EPJC_79_290: This is Figure 10b \n");
    myBlue->SolveAccImp(0.001);
    myBlue->DisplayAccImp(0,FilNam);

    // Solve scanning the statistical uncertainties
    printf("... B_EPJC_79_290: These are the auxiliary figures 3c,d from");
    printf(" %s \n", WebNam.Data());
    myBlue->ReleaseInp();
    myBlue->FixInp();
    myBlue->SolveScaSta();
    myBlue->PrintScaSta(FilNam, 172.38, 172.98, 0.34, 0.78);
  }

  // A final solve to obtain the Latex files
  myBlue->ReleaseInp();
  myBlue->FixInp();
  myBlue->Solve();
  if(Flag == 0){    
    printf("... B_EPJC_79_290: The 7 TeV combination. See Table 6\n");
  }else if(Flag == 1){
    printf("... B_EPJC_79_290: The 8 TeV combination. See Table 6\n");
  }else if(Flag == 2){
    printf("... B_EPJC_79_290: The 3-channel combination. See Table 6\n");
  }else if(Flag == 3){
    printf("... B_EPJC_79_290: The mtop(3) combination. See Table 6\n");
    printf("... B_EPJC_79_290: There is a typy introduced by the journal");
    printf(" in the paper.\n");
    printf("... B_EPJC_79_290: The combined value is mtop(3)=173.51 GeV\n");
    printf("... B_EPJC_79_290: instead of mtop(3)=173.13 GeV, see third line");
    printf(" of Fig 10b.\n");
    printf("... B_EPJC_79_290: Table 6 in arXiv:1810.01772 is correct.\n");
  }else if(Flag == 4){
    printf("... B_EPJC_79_290: The full combination. See Table 6\n");
    printf("... B_EPJC_79_290: There is a typy introduced by the journal");
    printf(" in the paper.\n");
    printf("... B_EPJC_79_290: The combined value is mtop=172.69 GeV\n");
    printf("... B_EPJC_79_290: instead of mtop(3)=174.08 GeV, see last line");
    printf(" of Fig 10b.\n");
    printf("... B_EPJC_79_290: Table 6 in arXiv:1810.01772 is correct.\n");
  }  
  myBlue->PrintResult();
  if(Flag == 2){
    // The compatibility of the three observables
    myBlue->PrintCompatObs();

    // The correlation of the three observables
    myBlue->GetRhoRes(LocRhoRes);
    printf("... B_EPJC_79_290:\n");
    printf("... B_EPJC_79_290: The correlations of the observables,");
    printf(" see Section 10.3. \n");
    printf("... B_EPJC_79_290:\n");
    myBlue->PrintMatrix(LocRhoRes, ForRho);
  }
  myBlue->LatexResult(FilNam);

  // LatexResult with different formats
  if(Flag == 2){
    const TString DifVal = "%7.3f";
    const TString DifUnc = "%4.1f";
    const TString DifWei = DifUnc;
    const TString DifRho = "%6.3f";
    const TString DifPul = "%+4.1f";
    const TString DifChi = "%+5.1f";
    const TString DifUni = "NewUnit";
    FilNam = FilNam + "_diff_format";
    myBlue->SetFormat(DifUni);
    myBlue->LatexResult(FilNam, DifVal, DifUnc, DifWei, DifRho, DifPul, DifChi);
  }

  // Delete object and matrices
  delete myBlue; myBlue = NULL;
  
  InpEst->Delete(); InpEst = NULL; 
  InpSta->Delete(); InpSta = NULL; 
  InpRho->Delete(); InpRho = NULL; 

  LocRhoRes->Delete(); LocRhoRes = NULL;

  // Return
  return;
};
