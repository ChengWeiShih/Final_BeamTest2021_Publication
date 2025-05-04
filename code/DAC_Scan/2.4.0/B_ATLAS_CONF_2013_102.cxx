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

void FilMat(Double_t RE, Double_t RL, Double_t *M);

void B_ATLAS_CONF_2013_102(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The   LHC combination
  //     [Shows how to retrieve the Chiq information of the result]
  //     [Shows how to use DisplayResult() and LatexResult()]
  //  1: The ATLAS combination
  //  2: The   CMS combination
  //  3: The combination for three observables = l+jets, di-lepton and all-jets
  //     [Shows how to retrieve the Pull information of the result]
  //     [Shows how to use DisplayResult() with user formats]
  //  4: The combination for   two observables = ATLAS and CMS 
  //  5: The systematic variations performed for Figure 2a,b
  //  6: The systematic variations performed for Figure 2c,d
  //  7: Numerical exercise Flag = 0 but with unpysical ad-hoc 
  //     relative systematic uncertainties
  //  8: The part of 5 scanning Rho_exp and rho_LHC simultaneously, 
  //     but using SolveScaRho
  //  9: The part of 5 scanning Rho_exp and rho_LHC  independently, 
  //     but using SolveScaRho
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  5;
  static const Int_t NumUnc = 19;
  Int_t NumObs =  1;
  static const Int_t MaxObs =  3;

  TString NamEst[NumEst] = {"ATL l+j", "ATL dil", "CMS l+j", 
			    "CMS dil", "CMS had"};
  TString NamUnc[NumUnc] = {"   Stat", "   iJES", " uncJES", " insJES", 
			    " intJES", " flaJES", "   bJES", "     MC", 
			    "    Rad", "     CR", "     UE", "    PDF", 
			    "   DTMO", "  b-tag", "   Lept", "   BGMC", 
			    "   BGDT", "   Meth", " Pileup"};
  TString NamObs[MaxObs] = {"   mtop","   mtop","   mtop"};

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {NumEst*0};

  // Index for the scale factors
  Int_t IWhichFac[NumEst*NumEst] = {
    0, 0, 1, 1, 1,
    0, 0, 1, 1, 1,
    1, 1, 0, 0, 0,
    1, 1, 0, 0, 0,
    1, 1, 0, 0, 0};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_ATLAS_CONF_2013_102: -----------------------------------\n");
    printf("... B_ATLAS_CONF_2013_102: The 2013 LHC combination, Flag = %2i\n",
	   Flag);
  }else if(Flag == 1){
    printf("... B_ATLAS_CONF_2013_102: ----------------------------------");
    printf("---\n");
    printf("... B_ATLAS_CONF_2013_102: The 2013 ATLAS combination");
    printf(", Flag = %2i\n", Flag);
    NamObs[0] = "ATLmtop";
  }else if(Flag == 2){
    printf("... B_ATLAS_CONF_2013_102: ----------------------------------- \n");
    printf("... B_ATLAS_CONF_2013_102: The 2013 CMS combination, Flag = %2i \n",
	   Flag);
    NamObs[0] = "CMSmtop";
  }else if(Flag == 3){
    printf("... B_ATLAS_CONF_2013_102: ---------------------------------");
    printf("--------------------------------- \n");
    printf("... B_ATLAS_CONF_2013_102: The 2013 combination for l+jets,");
    printf(" di-lepton and all-jets, Flag = %2i \n",Flag);
    NumObs = 3;
    IWhichObs[0] = 0;
    IWhichObs[1] = 1;
    IWhichObs[2] = 0;
    IWhichObs[3] = 1;
    IWhichObs[4] = 2;
    NamObs[0] = " mt-l+j";
    NamObs[1] = " mt-dil";
    NamObs[2] = " mt-had";
  }else if(Flag == 4){
    printf("... B_ATLAS_CONF_2013_102: --------------------------------");
    printf("----------------------------------------- \n");
    printf("... B_ATLAS_CONF_2013_102: The 2013 combination for 2 correlated");
    printf(" observables (ATLAS, CMS), Flag = %2i \n",Flag);
    NumObs = 2;
    IWhichObs[0] = 0;
    IWhichObs[1] = 0;
    IWhichObs[2] = 1;
    IWhichObs[3] = 1;
    IWhichObs[4] = 1;
    NamObs[0] = "ATLmtop";
    NamObs[1] = "CMSmtop";
  }else if(Flag == 5){
    printf("... B_ATLAS_CONF_2013_102: --------------------------------");
    printf("------------------ \n");
    printf("... B_ATLAS_CONF_2013_102: The systematic variations for Fig");
    printf(" 2a,b, Flag = %2i \n",Flag);
  }else if(Flag == 6){
    printf("... B_ATLAS_CONF_2013_102: --------------------------------");
    printf("----------------- \n");
    printf("... B_ATLAS_CONF_2013_102: The systematic variations for Fig");
    printf(" 2c,d, Flag = %2i \n",Flag);
  }else if(Flag == 7){
    printf("... B_ATLAS_CONF_2013_102: --------------------------------");
    printf("-------------------------------------------------- \n");
    printf("... B_ATLAS_CONF_2013_102: Numerical exercise the 2013 LHC");
    printf(" combination with relative uncertainties, Flag = %2i \n",Flag);
  }else if(Flag == 8){
    printf("... B_ATLAS_CONF_2013_102: --------------------------------");
    printf("------------------ \n");
    printf("... B_ATLAS_CONF_2013_102: The systematic variations for Fig");
    printf(" 2a,b, Flag = %2i \n",Flag);
    printf("... B_ATLAS_CONF_2013_102: The part of 5 scanning Rho_exp and");
    printf(" rho_LHC simultaneously using SolveScaRho() \n");
  }else if(Flag == 9){
    printf("... B_ATLAS_CONF_2013_102: --------------------------------");
    printf("------------------ \n");
    printf("... B_ATLAS_CONF_2013_102: The systematic variations for Fig");
    printf(" 2a,b, Flag = %2i \n",Flag);
    printf("... B_ATLAS_CONF_2013_102: The part of 5 scanning Rho_exp and");
    printf(" rho_LHC independently using SolveScaRho() \n");
  }else{
    printf("... B_ATLAS_CONF_2013_102: ------------------------ \n");
    printf("... ATLAS_CONF_2013_102: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Estimates 0-4
  // 0 == ATLAS_2011 lepton+jets
  // 1 == ATLAS_2011   di-lepton
  // 2 ==   CMS_2011 lepton+jets
  // 3 ==   CMS_2011   di-lepton
  // 4 ==   CMS_2011    all-jets

  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //             0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18
    //          Stat   iJES uncJES insJES intJES flaJES   bJES     MC    Rad     CR     UE    PDF   DTMO  b-tag   Lept   BGMC   BGDT   Meth Pileup
    // rho =       0      0  Cor02  Cor02  Cor04  Cor02  Cor04      1      1      1      1      1  Cor02  Cor04  Cor02      1      0      0      1
    172.31,     0.23,  0.72,  0.61,  0.29,  0.19,  0.36,  0.08,  0.33,  0.45,  0.32,  0.12,  0.17,  0.23,  0.81,  0.04,  0.00,  0.10,  0.13,  0.03,
    173.09,     0.64,  0.00,  0.73,  0.31,  0.39,  0.02,  0.71,  0.48,  0.37,  0.29,  0.42,  0.12,  0.22,  0.46,  0.12,  0.14,  0.00,  0.07,  0.01,
    173.49,     0.27,  0.33,  0.24,  0.02,  0.01,  0.11,  0.61,  0.02,  0.30,  0.54,  0.15,  0.07,  0.24,  0.12,  0.02,  0.13,  0.00,  0.06,  0.07,
    172.50,     0.43,  0.00,  0.69,  0.35,  0.08,  0.58,  0.76,  0.04,  0.58,  0.13,  0.05,  0.09,  0.18,  0.09,  0.14,  0.05,  0.00,  0.40,  0.11,
    173.49,     0.69,  0.00,  0.69,  0.35,  0.08,  0.58,  0.49,  0.19,  0.33,  0.15,  0.20,  0.06,  0.28,  0.06,  0.00,  0.00,  0.13,  0.13,  0.06
  };

  // Correlations 
  //              0,  1, 16, 17 <==}  rho_exp=0, rho_LHC=   0 <==> rho == 0 
  //          2,  3,  5, 12, 14 <==}  rho_exp=1, rho_LHC=   0 <==> Cor02 
  //                  4,  6, 13 <==}  rho_exp=1, rho_LHC= 0.5 <==> Cor04
  //  7,  8,  9, 10, 11, 15, 18 <==}  rho_exp=1, rho_LHC=   1 <==> rho == 1

  // The default values 
  //                             0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
  Double_t   RhoExp[NumUnc] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0};
  Double_t   RhoLHC[NumUnc] = {0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.5, 0.0, 1.0, 0.0, 0.0, 1.0};

  // The systematic variations see Fig 2a,b
  // We have 10 tests for f*RhoExp, f*RhoLHC and f*both
  static const Int_t NumRho = 10;
  //                    j =    0    1    2    3    4    5    6    7    8    9
  Double_t RhoFac[NumRho] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};
  Double_t ValExp[NumRho] = {0};
  Double_t UncExp[NumRho] = {0};
  Double_t ValLHC[NumRho] = {0};
  Double_t UncLHC[NumRho] = {0};
  Double_t ValAll[NumRho] = {0};
  Double_t UncAll[NumRho] = {0};
  Double_t ActExp = 0, ActLHC = 0;

  // The systematic variations see Fig 2c,d
  // We have 13 tests, 10 correlations and three changed sys components
  // The correlations
  static const Int_t NumCor = 10;
  static const Int_t NumSys = 13;
  //                    j =    0    1    2    3    4    5    6    7    8    9
  Int_t    IndSys[NumCor] = {  3,   4,   5,   5,   6,  12,  13,  13,  13,  11};
  Double_t SysExp[NumCor] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 1.0};
  Double_t SysLHC[NumCor] = {0.5, 1.0, 0.5, 1.0, 1.0, 0.5, 0.0, 1.0, 0.5, 0.0};
  Double_t ValSys[NumSys] = {0};
  Double_t UncSys[NumSys] = {0};

  // The different systematic assumptions for Sys = 10-12
  // 10)   CMS: add quadratically    0.58, 0.76, 0.93   to Unc=7 for Est=2,3,4
  // 11) ATLAS: remove quadratically 0.27, 0.44       from Unc=7 for Est=0,1
  // 12)   CMS: use 10) and add quadratically 0.18      to Unc=6 for Est=2,3,4
  // The values to be used in the systematic variations
  Double_t CMSHad[3] = {0.58, 0.76, 0.93};
  Double_t   CMSBJes = 0.18;
  Double_t ATLHad[2] = {0.27, 0.44};
  // The values to do the reset after the previous systematic variations
  Double_t       CMSMC[3] = {XEst[2 * (NumUnc+1)+7+1],XEst[3 * (NumUnc+1)+7+1],
			     XEst[4 * (NumUnc+1)+7+1]};
  Double_t       ATLMC[2] = {XEst[0 * (NumUnc+1)+7+1],XEst[1 * (NumUnc+1)+7+1]};

  static const Int_t LenCor = NumEst * NumEst;
  Double_t Corxx[LenCor] = {0};

  //-- Local Structures for Blue output
  Int_t iok = 0;
  static const Int_t LocObs = 1;
  static const Int_t LenResult = LocObs * NumUnc+1;

  // Number of groups in IWhichFac
  static const Int_t NumSca =  2;

  // Number of rho values scanned, the default is 10
  static const Int_t NumVal = 10;
  static const Int_t LenScaRho = NumUnc * NumSca * NumVal;

  // TMatrices
  TMatrixD* LocScaVal = new TMatrixD(NumUnc*NumSca,NumVal);
  TMatrixD* LocScaUnc = new TMatrixD(NumUnc*NumSca,NumVal);

  // Double_t Arrays
  Double_t LocResultArr[LenResult] = {0};
  Double_t LocUncertArr[LocObs] = {0};
  Double_t LocScaValArr[LenScaRho] = {0};
  Double_t LocScaUncArr[LenScaRho] = {0};
  //-- End
  
  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%5.3f";
  const TString ForRho = ForWei;
  const TString ForPul = ForVal;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";

  // Construct Object
  Blue* myBlue;
  if(Flag == 9){
    myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0], &IWhichFac[0]);
  }else{
    myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  }
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);

  // Fill estimates and correlations
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }
  for(Int_t k = 0; k<NumUnc; k++){
    FilMat(RhoExp[k], RhoLHC[k], &Corxx[0]);
    myBlue->FillCor(k,&Corxx[0]);
  }
  Double_t Value = 0, Sigma = 0;
  if(Flag == 0){    
    
    myBlue->FixInp();
    myBlue->PrintCompatEst();
    myBlue->PrintEst();

    // The parameters of the combination and the pairs of estimates
    myBlue->PrintParams();
    TString FilBas = "B_ATLAS_CONF_2013_102";
    TString FilNam;
    char Buffer[100];
    for(Int_t i = 0; i<NumEst; i++){
      for(Int_t j = i+1; j<NumEst; j++){
	sprintf(Buffer,"%1i%1i",i,j);
	FilNam = &Buffer[0];       
	FilNam = FilBas + "_" + FilNam;
	myBlue->InspectPair(i,j,FilNam);
      }
    }
    printf("... B_ATLAS_CONF_2013_102: The 2013 LHC Combination Flag = %2i. \n",
	   Flag);
    printf("... B_ATLAS_CONF_2013_102: For the correlation of the input");
    printf(" see matrix on page 9 \n");
    myBlue->PrintRho();
    myBlue->Solve();
    myBlue->PrintResult();
    myBlue->LatexResult("B_ATLAS_CONF_2013_102");
    myBlue->DisplayResult(0,FilBas);

    // The information weights
    myBlue->ReleaseInp();
    myBlue->FixInp();
    myBlue->SolveInfWei();
    printf("... B_ATLAS_CONF_2013_102: For the results see Table 1+7 and");
    printf(" for the weight matrix see Table 2 \n");
    myBlue->PrintResult();
    printf("... B_ATLAS_CONF_2013_102: For the pulls see Fig1c \n");
    myBlue->PrintPull();
    printf("... B_ATLAS_CONF_2013_102: ChiQua = %5.3f for", myBlue->GetChiq());
    printf(" NDof = %2i, Probability = %5.2f %% \n", myBlue->GetNdof(),
	   100*myBlue->GetProb());
    printf("... B_ATLAS_CONF_2013_102: The information weights, Table 2 \n");
    myBlue->PrintInfWei();

    printf("... B_ATLAS_CONF_2013_102: The maximaization of the variance");
    printf(" see Table 5 \n");
    for(Int_t l = 0; l<3; l++){
      myBlue->ResetInp();
      myBlue->FixInp();
      myBlue->SolveMaxVar(l);
      myBlue->PrintMaxVar();
      myBlue->PrintResult();
    }
  }else if(Flag == 1){
    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(3);    
    myBlue->SetInActiveEst(4);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2013_102: The 2013 ATLAS combination Flag = %2i.",
	   Flag);
    printf(" For the results see Table 7 \n");
    myBlue->PrintResult();
    myBlue->PrintPull();
  }else if(Flag == 2){
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(1);    
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2013_102: The 2013 CMS combination Flag = %2i.",
	   Flag);
    printf(" For the results see Table 7 \n");
    myBlue->PrintResult();
    myBlue->PrintPull();
  }else if(Flag == 3){
    myBlue->FixInp();
    myBlue->PrintNamEst();
    myBlue->PrintNamUnc();
    myBlue->PrintNamObs();
    myBlue->PrintCompatEst();
    myBlue->PrintParams();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2013_102: The 2013 combination for l+jets,");
    printf(" di-lepton and all-jets Flag = %2i \n",Flag);
    printf("... B_ATLAS_CONF_2013_102: For the results see Table 3 \n");
    myBlue->PrintResult();
    myBlue->PrintPull();
    myBlue->PrintRhoRes();
    myBlue->PrintCompatObs();
    for(Int_t i = 0; i<NumEst; i++){
      printf("... B_ATLAS_CONF_2013_102: Pull(%2i) = %5.3f \n", 
	     i, myBlue->GetPull(i));
    }
    myBlue->LatexResult("B_ATLAS_CONF_2013_102_3Obs");
    for(Int_t n = 0; n<NumObs; n++){
      myBlue->DisplayResult(n,"B_ATLAS_CONF_2013_102_3Obs",ForVal,ForUnc);
    }
  }else if(Flag == 4){
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_ATLAS_CONF_2013_102: The 2013 combination for 2 correlated");
    printf(" observables (ATLAS, CMS) Flag = %2i \n",Flag);
    printf("... B_ATLAS_CONF_2013_102: For the results see Table 4 \n");
    printf("... B_ATLAS_CONF_2013_102: The Chiq Probabaility differs.");
    printf(" This is under investigation with the authors.\n");
    myBlue->PrintResult();
    myBlue->PrintRhoRes();
    myBlue->PrintPull();
    myBlue->PrintCompatObs();
  }else if(Flag == 5){
    printf("... B_ATLAS_CONF_2013_102: The systematic variations for Fig");
    printf(" 2a,b. Flag = %2i \n",Flag);
    myBlue->FixInp();
    myBlue->Solve();

    // First get result and total uncertainty
    iok = myBlue->GetResult(LocResultArr);
    if(iok == 1){Value = LocResultArr[0];
    }else{printf("... B_ATLAS_CONF_2013_102: Error in GetResult \n");
    }
    iok = myBlue->GetUncert(LocUncertArr);
    if(iok == 1){Sigma = LocUncertArr[0];
    }else{printf("... B_ATLAS_CONF_2013_102: Error in GetUncert \n");
    }
    printf("... B_ATLAS_CONF_2013_102: The Default result:");
    printf(" Value +- Sigma  = %5.3F +- %5.3F \n", Value, Sigma);
    
    // Fig 2a, b
    // Loop over the three scenarios and the ten fractions
    for(Int_t l = 0; l<3; l++){
      for(Int_t j = 0; j<NumRho; j++){
	// Construct Object
	Blue *myBlueI = new Blue(NumEst, NumUnc);
	myBlueI->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);
	
	// Fill estimates
	ind = 0;
	for(Int_t i = 0; i<NumEst; i++){
	  myBlueI->FillEst(i,&XEst[ind]);
	  ind = ind + NumUnc + 1;
	}

	// Fill changed correlations
	for(Int_t k = 0; k<NumUnc; k++){
	  ActExp = RhoExp[k];
	  ActLHC = RhoLHC[k];
	  if(l == 0 || l == 2)ActExp = RhoFac[j] * RhoExp[k];
	  if(l == 1 || l == 2)ActLHC = RhoFac[j] * RhoLHC[k];
	  FilMat(ActExp, ActLHC, &Corxx[0]);
	  myBlueI->FillCor(k,&Corxx[0]);
	}
	myBlueI->FixInp();
	myBlueI->Solve();
	// Get individual result and uncertainty
	iok = myBlueI->GetResult(LocResultArr);
	iok = myBlueI->GetUncert(LocUncertArr);
	if(l == 0){
	  ValExp[j] = LocResultArr[0];
	  UncExp[j] = LocUncertArr[0];
	}else if(l == 1){
	  ValLHC[j] = LocResultArr[0];
	  UncLHC[j] = LocUncertArr[0];
	}else if(l == 2){
	  ValAll[j] = LocResultArr[0];
	  UncAll[j] = LocUncertArr[0];
	}
	delete myBlueI;
	myBlueI = NULL;
      }
    }
    // Report the findings
    printf("... B_ATLAS_CONF_2013_102:         ");
    printf("The systematic variations from Fig 2a,b -- \n");
    printf("... B_ATLAS_CONF_2013_102:         ");
    printf("-- The differences to the default result -- \n");
    printf("... B_ATLAS_CONF_2013_102:         ");
    printf("----- Values -------   --- Uncertainties -- \n");
    printf("... B_ATLAS_CONF_2013_102: Factor  ");
    printf("RhoExp RhoLHC   Both   RhoExp RhoLHC   Both \n");
    for(Int_t j = 0; j<NumRho; j++){
      printf("... B_ATLAS_CONF_2013_102:    ");
    printf("%3.1F  %+5.3F %+5.3F %+5.3F   %+5.3F %+5.3F %+5.3F \n", RhoFac[j],
	     ValExp[j] - Value, ValLHC[j] - Value, ValAll[j] - Value,
	     UncExp[j] - Sigma, UncLHC[j] - Sigma, UncAll[j] - Sigma);
    }
  }else if(Flag == 6){
    printf("... B_ATLAS_CONF_2013_102: The systematic variations for Fig");
    printf(" 2c,d. Flag = %2i \n",Flag);
    myBlue->FixInp();
    myBlue->Solve();

    // Get result and total uncertainty
    iok = myBlue->GetResult(LocResultArr);
    if(iok == 1){Value = LocResultArr[0];
    }else{printf("... B_ATLAS_CONF_2013_102: Error in GetResult \n");
    }
    iok = myBlue->GetUncert(LocUncertArr);
    if(iok == 1){Sigma = LocUncertArr[0];
    }else{printf("... B_ATLAS_CONF_2013_102: Error in GetUncert \n");
    }
    printf("... B_ATLAS_CONF_2013_102: The Default result:");
    printf(" Value +- Sigma = %5.3F +- %5.3F \n", Value, Sigma);
    
    // Fig 2c, d
    for(Int_t j = 0; j<NumSys; j++){
      // Construct Object
      Blue *myBlueI = new Blue(NumEst, NumUnc);
      myBlueI->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);
      
      // Change estimates for the last three cases
      Int_t Ioff = 0;
      // If needed revert to original, then change to new
      if(j == NumCor){
	// CMS add had
	Ioff = 2 * (NumUnc+1)+7+1;
	for(Int_t n = 0; n<3; n++){
	  XEst[Ioff] = TMath::Sqrt(XEst[Ioff]*XEst[Ioff] + CMSHad[n]*CMSHad[n]);
	  Ioff = Ioff + (NumUnc+1);
	}
      }else if(j == NumCor+1){
	// First revert CMS
	Ioff = 2 * (NumUnc+1)+7+1;
	for(Int_t n = 0; n<3; n++){
	  XEst[Ioff] = CMSMC[n];
	  Ioff = Ioff + (NumUnc+1);
	}
	// Now remove had for ATLAS
	Ioff = 0 * (NumUnc+1)+7+1;
	for(Int_t n = 0; n<2; n++){
	  XEst[Ioff] = TMath::Sqrt(XEst[Ioff]*XEst[Ioff] - ATLHad[n]*ATLHad[n]);
	  Ioff = Ioff + (NumUnc+1);
	}
      }else if(j == NumCor+2){
	// First revert ATLAS
	Ioff = 0 * (NumUnc+1)+7+1;
	for(Int_t n = 0; n<2; n++){
	  XEst[Ioff] = ATLMC[n];
	  Ioff = Ioff + (NumUnc+1);
	}
	// CMS use the had from above and change the bJes
	//had
	Ioff = 2 * (NumUnc+1)+7+1;
	for(Int_t n = 0; n<3; n++){
	  XEst[Ioff] = TMath::Sqrt(XEst[Ioff]*XEst[Ioff] + CMSHad[n]*CMSHad[n]);
	  Ioff = Ioff + (NumUnc+1);
	}
	//bjes
	Ioff = 2 * (NumUnc+1)+6+1;
	for(Int_t n = 0; n<3; n++){
	  XEst[Ioff] = CMSBJes;
	  Ioff = Ioff + (NumUnc+1);
	}
      }

      // Now fill the estimates
      ind = 0;
      for(Int_t i = 0; i<NumEst; i++){
	myBlueI->FillEst(i,&XEst[ind]);
	ind = ind + NumUnc + 1;
      }
      
      // Fill correlations change for the first NumCor cases
      for(Int_t k = 0; k<NumUnc; k++){
	ActExp = RhoExp[k];
	ActLHC = RhoLHC[k];
	if(j < NumCor && k == IndSys[j]){
	  ActExp = SysExp[j];
	  ActLHC = SysLHC[j];
	  printf("... B_ATLAS_CONF_2013_102: %2i  %2i %5.3F %5.3F \n", 
		 j,IndSys[j],ActExp,ActLHC);
	}
	FilMat(ActExp, ActLHC, &Corxx[0]);
	myBlueI->FillCor(k,&Corxx[0]);
      }
      myBlueI->FixInp();
      if(j >= NumCor){
	if(j == NumCor){
	  printf("... B_ATLAS_CONF_2013_102: Changed   CMS estimates k=7\n");
	}else if(j == NumCor+1){
	  printf("... B_ATLAS_CONF_2013_102: Changed ATLAS estimates k=7\n");
	}else if(j == NumCor+2){
	  printf("... B_ATLAS_CONF_2013_102: Changed   CMS estimates k=6,7\n");
	}
	myBlueI->PrintEst();
      }
      myBlueI->Solve();
      // Get individual result and uncertainty
      iok = myBlueI->GetResult(LocResultArr);
      iok = myBlueI->GetUncert(LocUncertArr);
      ValSys[j] = LocResultArr[0];
      UncSys[j] = LocUncertArr[0];
      delete myBlueI;
      myBlueI = NULL;
    }

    // Report the findings
    printf("... B_ATLAS_CONF_2013_102:  ");
    printf("The systematic variations from Fig 2c,d-- \n");
    printf("... B_ATLAS_CONF_2013_102:  ");
    printf("-- The differences to the default result -- \n");
    printf("... B_ATLAS_CONF_2013_102:  ");
    printf("j =  0 -  9 <==> F*Rho_i \n");
    printf("... B_ATLAS_CONF_2013_102:  ");
    printf("j = 10 - 12 <==> Changed systematics \n");
    printf("... B_ATLAS_CONF_2013_102:  ");
    printf("------------------------------------------- \n");
    printf("... B_ATLAS_CONF_2013_102:  ");
    printf("   Variation   Value  Uncertainty \n");
    for(Int_t j = 0; j<NumSys; j++){
      printf("... B_ATLAS_CONF_2013_102:");
      printf("            %2i  %+5.3F       %+5.3F \n", j,
	     ValSys[j] - Value, UncSys[j] - Sigma);
    }
  }else if(Flag == 7){
    myBlue->SetRelUnc();
    myBlue->SetNotRelUnc(0);
    myBlue->FixInp();
    myBlue->SolveRelUnc(0.1);
    printf("... B_ATLAS_CONF_2013_102: Numerical exercise mtop full");
    printf(" combination, using relative syst uncertainties \n");
    myBlue->PrintResult();
    myBlue->LatexResult("B_ATLAS_CONF_2013_102_Rel");
    myBlue->InspectLike(0,"B_ATLAS_CONF_2013_102_7");
    myBlue->PrintInspectLike();  
  }else if(Flag == 8){
    // Both at the same time
    printf("... B_ATLAS_CONF_2013_102: Simultaneous variations of");
    printf(" rho_exp and rho_LHC \n");
    myBlue->FixInp();
    myBlue->SolveScaRho(1);
    myBlue->PrintScaRho("B_ATLAS_CONF_2013_102");
  }else if(Flag == 9){
    // Independent
    printf("... B_ATLAS_CONF_2013_102: Independent variations of");
    printf(" rho_exp and rho_LHC \n");
    myBlue->FixInp();
    myBlue->SolveScaRho(0);
    myBlue->PrintScaRho("B_ATLAS_CONF_2013_102");

    // Get the differences in values
    myBlue->GetScaVal(0, LocScaVal); 
    //LocScaVal->Print();
    myBlue->GetScaVal(0, LocScaValArr);
    ind = 0;
    for(Int_t l = 0; l<myBlue->GetNumScaFac()+1; l++){
      printf("... B_ATLAS_CONF_2013_102:");
      printf(" Next Group of Differences in Values l = %2i \n",l);
      for(Int_t k = 0; k<myBlue->GetActUnc(); k++){
	printf("... B_ATLAS_CONF_2013_102:");
	for(Int_t ll = 0; ll<myBlue->GetNumScaRho(); ll++){
	  printf(" %+5.3f",LocScaValArr[ind]);
	  if(ll == NumVal-1)printf("\n");
	  ind = ind + 1;
	}
      }
    }
    
    // Get the differences in uncertainties
    myBlue->GetScaUnc(0, LocScaUnc); 
    //LocScaUnc->Print();
    myBlue->GetScaUnc(0, LocScaUncArr);
    ind = 0;
    for(Int_t l = 0; l<myBlue->GetNumScaFac(); l++){
      printf("... B_ATLAS_CONF_2013_102:");
      printf(" Next Group of Differences in Uncertainties l = %2i \n",l);
      for(Int_t k = 0; k<myBlue->GetActUnc(); k++){
	printf("... B_ATLAS_CONF_2013_102:");
	for(Int_t ll = 0; ll<myBlue->GetNumScaRho(); ll++){
	  printf(" %+5.3f",LocScaUncArr[ind]);
	  if(ll == NumVal-1)printf("\n");
	  ind = ind + 1;
	}
      }
    }
  }
  //Delete Object
  delete myBlue;  myBlue = NULL;
  LocScaVal->Delete(); LocScaVal = NULL;
  LocScaUnc->Delete(); LocScaUnc = NULL;
};

//------------------------------------------------------------------------
// Utility to fill a correlation matrix for any pair of RhoExp and RhoLHC
//------------------------------------------------------------------------

void FilMat(Double_t RE, Double_t RL, Double_t *M){
  M[ 0] =  1; M[ 1] = RE; M[ 2] = RL; M[ 3] = RL; M[ 4] = RL;
  M[ 5] = RE; M[ 6] =  1; M[ 7] = RL; M[ 8] = RL; M[ 9] = RL;
  M[10] = RL; M[11] = RL; M[12] =  1; M[13] = RE; M[14] = RE;
  M[15] = RL; M[16] = RL; M[17] = RE; M[18] =  1; M[19] = RE;
  M[20] = RL; M[21] = RL; M[22] = RE; M[23] = RE; M[24] =  1;
  return;
};
