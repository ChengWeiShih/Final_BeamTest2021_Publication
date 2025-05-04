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

void FilMat(Double_t CDF, Double_t DZE, Double_t ATL, Double_t CMS,
	    Double_t LHC, Double_t TEV, Double_t ATE, Double_t CTE,
	    Double_t *M);

void B_arXiv_1403_4427(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  // 0: The 2014 World mtop combination
  // 1: The combined rho variations
  // 2: The individual rho variations, Group = 0/1/2/3 = EXP/LHC/TEV/COL
  // 3: The systematic variations
  // 4: The combination per decay channel
  // 5: The combination per experiment
  // 6: The combination per collider
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 11;
  static const Int_t NumUnc = 16;
  static const Int_t MaxObs =  4;
  Int_t NumObs = 1;
  static const Int_t NumFac =  4;

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {NumEst*0};

  TString NamObs[MaxObs] = {"   mtop","   mtop","   mtop","   mtop"};

  TString NamEst[NumEst] = {"CDF l+j", "CDF dil", "CDF had", "CDF Met", 
			    "DZE l+j", "DZE dil",
			    "ATL l+j", "ATL dil", 
			    "CMS l+j", "CMS dil", "CMS had"};
  TString NamUnc[NumUnc] = {"   Stat", "   iJES", " stdJES", "flavJES", 
			    "   bJES", "     MC", "    Rad", "     CR", 
			    "    PDF", " DeTMod", "  b-tag", "  LepPt", 
			    "   BGMC", " BGData", "   Meth", "    MHI"};
  // The filemanes
 TString NamFil = "B_arXiv_1403_4427";
 TString NamCha = NamFil + "_Channel";
 TString NamExp = NamFil + "_Experiment";
 TString NamCol = NamFil + "_Collider";

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_arXiv_1403_4427: -----------------------------------------\n");
    printf("... B_arXiv_1403_4427: The 2014 World mtop combination Flag = %2i\n",
	   Flag);
  }else if(Flag == 1){
    printf("... B_arXiv_1404427: -------------------------------------------\n");
    printf("... B_arXiv_1403_4427: The rho variation Group_0= ALL, Flag = %2i\n",
	   Flag);
  }else if(Flag == 2){
    printf("... B_arXiv_1403_4427: -----------------------------------\n");
    printf("... B_arXiv_1403_4427: The individual rho scans Flag = %2i\n",Flag);
    printf("... B_arXiv_1403_4427: Group = 0/1/2/3 = EXP/LHC/TEV/COL \n");
  }else if(Flag == 3){
    printf("... B_arXiv_1403_4427: ------------------------- \n");
    printf("... B_arXiv_1403_4427: The systematics Flag = %2i \n",Flag);
  }else if(Flag == 4){
    printf("... B_arXiv_1403_4427: ----------------------------------");
    printf("----------\n");
    printf("... B_arXiv_1403_4427: The combination per decay channel,");
    printf(" Flag = %2i\n", Flag);

    NumObs = 4;
    IWhichObs[ 0] = 0;
    IWhichObs[ 1] = 1;
    IWhichObs[ 2] = 2;
    IWhichObs[ 3] = 3;
    IWhichObs[ 4] = 0;
    IWhichObs[ 5] = 1;
    IWhichObs[ 6] = 0;
    IWhichObs[ 7] = 1;
    IWhichObs[ 8] = 0;
    IWhichObs[ 9] = 1;
    IWhichObs[10] = 2;

    NamObs[0] = " mt-l+j";  
    NamObs[1] = " mt-dil";
    NamObs[2] = " mt-had";
    NamObs[3] = " mt-met";
  }else if(Flag == 5){
    printf("... B_arXiv_1403_4427: --------------------------------");
    printf(" ---------- \n");
    printf("... B_arXiv_1403_4427: The combination per experiment,");
    printf(" Flag = %2i \n",Flag);

    NumObs = 4;
    IWhichObs[ 0] = 0;
    IWhichObs[ 1] = 0;
    IWhichObs[ 2] = 0;
    IWhichObs[ 3] = 0;
    IWhichObs[ 4] = 1;
    IWhichObs[ 5] = 1;
    IWhichObs[ 6] = 2;
    IWhichObs[ 7] = 2;
    IWhichObs[ 8] = 3;
    IWhichObs[ 9] = 3;
    IWhichObs[10] = 3;

    NamObs[0] = " mt-CDF";
    NamObs[1] = " mt-DZE";
    NamObs[2] = " mt-ATL";
    NamObs[3] = " mt-CMS";
  }else if(Flag == 6){
    printf("... B_arXiv_1403_4427: ------------------------");
    printf("-------------------\n");
    printf("... B_arXiv_1403_4427: The combination per collider,");
    printf(" Flag = %2i \n",Flag);

    NumObs = 2;
    IWhichObs[ 0] = 0;
    IWhichObs[ 1] = 0;
    IWhichObs[ 2] = 0;
    IWhichObs[ 3] = 0;
    IWhichObs[ 4] = 0;
    IWhichObs[ 5] = 0;
    IWhichObs[ 6] = 1;
    IWhichObs[ 7] = 1;
    IWhichObs[ 8] = 1;
    IWhichObs[ 9] = 1;
    IWhichObs[10] = 1;

    NamObs[0] = " mt-TEV";
    NamObs[1] = " mt-LHC";
  }else{
    printf("... arXiv_1403_4427: Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  // Store the estimates Table 3
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //         0     1     2     3     4     5     6     7     8     9    10    11    12    13    14     15 
    //      Stat  iJES sdJES fvJES  bJES    MC   Rad    CR   PDF  DTMO  btag  Lept  BGMC  BGDT  Meth Pileup   
    172.85, 0.52, 0.49, 0.53, 0.09, 0.16, 0.56, 0.06, 0.21, 0.08, 0.00, 0.03, 0.03, 0.12, 0.16, 0.05,  0.07, 
    170.28, 1.95, 0.00, 2.99, 0.14, 0.33, 0.36, 0.22, 0.51, 0.31, 0.00, 0.00, 0.27, 0.24, 0.14, 0.12,  0.23, 
    172.47, 1.43, 0.95, 0.45, 0.03, 0.15, 0.49, 0.10, 0.32, 0.19, 0.00, 0.10, 0.00, 0.00, 0.56, 0.38,  0.08, 
    173.93, 1.26, 1.05, 0.44, 0.10, 0.17, 0.48, 0.28, 0.28, 0.16, 0.00, 0.00, 0.00, 0.00, 0.15, 0.21,  0.18, 
    174.94, 0.83, 0.47, 0.63, 0.26, 0.07, 0.63, 0.26, 0.28, 0.21, 0.36, 0.10, 0.18, 0.18, 0.21, 0.16,  0.05, 
    174.00, 2.36, 0.55, 0.56, 0.40, 0.20, 0.50, 0.30, 0.55, 0.30, 0.50, 0.00, 0.35, 0.00, 0.20, 0.51,  0.00, 
    172.31, 0.23, 0.72, 0.70, 0.36, 0.08, 0.35, 0.45, 0.32, 0.17, 0.23, 0.81, 0.04, 0.00, 0.10, 0.13,  0.03, 
    173.09, 0.64, 0.00, 0.89, 0.02, 0.71, 0.64, 0.37, 0.29, 0.12, 0.22, 0.46, 0.12, 0.14, 0.00, 0.07,  0.01, 
    173.49, 0.27, 0.33, 0.24, 0.11, 0.61, 0.15, 0.30, 0.54, 0.07, 0.24, 0.12, 0.02, 0.13, 0.00, 0.06,  0.07, 
    172.50, 0.43, 0.00, 0.78, 0.58, 0.76, 0.06, 0.58, 0.13, 0.09, 0.18, 0.09, 0.14, 0.05, 0.00, 0.40,  0.11, 
    173.49, 0.69, 0.00, 0.78, 0.58, 0.49, 0.28, 0.33, 0.15, 0.06, 0.28, 0.06, 0.00, 0.00, 0.13, 0.13,  0.06 
  };

  // Store the indices per experiment
  static const Int_t NumCDF = 4; Int_t IndCDF[NumCDF] = {0,  1, 2, 3};
  static const Int_t NumDZE = 2; Int_t IndDZE[NumDZE] = {4,  5};
  static const Int_t NumATL = 2; Int_t IndATL[NumATL] = {6,  7};
  static const Int_t NumCMS = 3; Int_t IndCMS[NumCMS] = {8,  9, 10};

  // Store the indices per CHANNEL
  static const Int_t NumLPJ = 4; Int_t IndLPJ[NumLPJ] = {0,  4, 6, 8};
  static const Int_t NumDIL = 4; Int_t IndDIL[NumDIL] = {1,  5, 7, 9};
  static const Int_t NumHAD = 2; Int_t IndHAD[NumHAD] = {2, 10};
  static const Int_t NumMET = 1; Int_t IndMET[NumMET] = {3};

  static const Int_t LenRho = NumEst * NumEst;
  Double_t RhoMat[LenRho] = {0};
  
  // Index for the scale factors
  // 0/1/2/3 = EXP/LHC/TEV/COL
  Int_t IWhichFac[NumEst*NumEst] = {
    0, 0, 0, 0, 2, 2, 3, 3, 3, 3, 3,
    0, 0, 0, 0, 2, 2, 3, 3, 3, 3, 3,
    0, 0, 0, 0, 2, 2, 3, 3, 3, 3, 3,
    0, 0, 0, 0, 2, 2, 3, 3, 3, 3, 3,
    2, 2, 2, 2, 0, 0, 3, 3, 3, 3, 3,
    2, 2, 2, 2, 0, 0, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 0, 0, 1, 1, 1,
    3, 3, 3, 3, 3, 3, 0, 0, 1, 1, 1,
    3, 3, 3, 3, 3, 3, 1, 1, 0, 0, 0,
    3, 3, 3, 3, 3, 3, 1, 1, 0, 0, 0,
    3, 3, 3, 3, 3, 3, 1, 1, 0, 0, 0
  };

  //-----------------------------------------------------
  // Store the correlations Table 4
  //-----------------------------------------------------
  // The structure of the correlation matrix
  // ATE = ATL-TEV, CTE = CMS-TEV 
  //   0   1   2   3   4   5   6   7   8   9  10
  // 1.0
  // CDF 1.0
  // CDF CDF 1.0
  // CDF CDF CDF 1.0 
  // TEV TEV TEV TEV 1.0
  // TEV TEV TEV TEV DZE 1.0
  // ATE ATE ATE ATE ATE ATE 1.0
  // ATE ATE ATE ATE ATE ATE ATL 1.0
  // CTE CTE CTE CTE CTE CTE LHC LHC 1.0
  // CTE CTE CTE CTE CTE CTE LHC LHC CMS 1.0
  // CTE CTE CTE CTE CTE CTE LHC LHC CMS CMS 1.0
  //   0   1   2   3   4   5   6   7   8   9  10
  // ---------------------------=-----------------
  // The different matrices per uncertainty source
  // ---------------------------------------------
  //         0, 13, 14 <==} CDF=0.0, DZE=0.0, ATL=0.0, CMS=0.0, LHC=0.0, TEV=0.0, ATE=0.0, CTE=0.0 <==> 0
  //                 1 <==} CDF=0.0, DZE=1.0, ATL=0.0, CMS=0.0, LHC=0.0, TEV=0.0, ATE=0.0, CTE=0.0 <==> 1
  // 2,  3,  9, 10, 11 <==} CDF=1.0, DZE=1.0, ATL=1.0, CMS=1.0, LHC=0.0, TEV=0.0, ATE=0.0, CTE=0.0 <==> 2
  //                 4 <==} CDF=1.0, DZE=1.0, ATL=1.0, CMS=1.0, LHC=0.5, TEV=1.0, ATE=1.0, CTE=0.5 <==> 3
  //             5,  7 <==} CDF=1.0, DZE=1.0, ATL=1.0, CMS=1.0, LHC=1.0, TEV=1.0, ATE=1.0, CTE=1.0 <==> 4
  //             6,  8 <==} CDF=1.0, DZE=1.0, ATL=1.0, CMS=1.0, LHC=1.0, TEV=1.0, ATE=0.5, CTE=0.5 <==> 5
  //                15 <==} CDF=1.0, DZE=1.0, ATL=1.0, CMS=1.0, LHC=1.0, TEV=0.0, ATE=0.0, CTE=0.0 <==> 6
  //

  // The index to the correct rho values per source k
  // Uncertainty =        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 
  Int_t RhoInd[NumUnc] = {0, 1, 2, 2, 3, 4, 5, 4, 5, 2, 2, 2,-1, 0, 0, 6};

  // The rho values per matrix
  static const Int_t IndMat = 7;
  //                           0    1    2    3    4    5    6
  Double_t RhoCDF[IndMat] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  Double_t RhoDZE[IndMat] = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  Double_t RhoATL[IndMat] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  Double_t RhoCMS[IndMat] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  Double_t RhoLHC[IndMat] = {0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
  Double_t RhoTEV[IndMat] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0};
  Double_t RhoATE[IndMat] = {0.0, 0.0, 0.0, 1.0, 1.0, 0.5, 0.0};
  Double_t RhoCTE[IndMat] = {0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 0.0};

  // Source k=12 is special it is like k=5 but rho=1 only for
  // the same decay channel ==> Make it a special matrix
  Double_t CoBGMC[LenRho] = {
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
    };

  //-- Local Structures for Blue output
  // TMatrices
  Int_t iok = 0;
  TMatrixD* LocRho    = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRhoRes = new TMatrixD(MaxObs,MaxObs);
  TMatrixD* LocRes    = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc    = new TMatrixD(NumObs,1);
  // -- End

  // Define formats for Figures and Latex file
  const TString ForVal = "%4.2f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%5.3f";
  const TString ForRho = ForVal;
  const TString ForPul = ForVal;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";

  // Construct Object according to Flag
  Blue *myBlue;
  if(Flag == 2 || Flag == 3){
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

  // Fill estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }
  
  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(k == 12){
      myBlue->FillCor(k,&CoBGMC[0]);
    }else{
      ind = RhoInd[k];    
      FilMat(RhoCDF[ind], RhoDZE[ind], RhoATL[ind], RhoCMS[ind],
	     RhoLHC[ind], RhoTEV[ind], RhoATE[ind], RhoCTE[ind], &RhoMat[0]);
      myBlue->FillCor(k,&RhoMat[0]);
    }
  }

  // Solve according to Flag
  if(Flag == 0){
    printf("... B_arXiv_1403_4427: The 2014 World mtop combination \n");

    // The central result
    myBlue->FixInp();
    myBlue->PrintCompatEst(NamFil);
    myBlue->SolveInfWei();
    printf("... B_arXiv_1403_4427: For the estimates and results see Table 3 \n");
    myBlue->PrintEst();
    myBlue->PrintResult();

    printf("... B_arXiv_1403_4427: For the correlation see Table 5 \n");
    iok = myBlue->GetRho(LocRho);
    if(iok == 1)myBlue->PrintMatrix(LocRho,"%5.2f");

    printf("... B_arXiv_1403_4427: For the information weights see Table 6 \n");
    myBlue->PrintInfWei();

    myBlue->LatexResult(NamFil);
    for(Int_t n = 0; n<NumObs; n++)myBlue->DisplayResult(n,NamFil);

  }else if(Flag == 1){

    // Rho_ALL
    // The correlated scan
    myBlue->FixInp();
    printf("... B_arXiv_1403_4427: The correlated scan of RHO_ALL see");
    printf(" Figure 4a,b \n");
    printf("... B_arXiv_1403_4427: The rightmost end of the histogram should be");
    printf(" be compared to the \n");
    printf("... B_arXiv_1403_4427: values at 40%%, 70%% and 100%% of the");
    printf(" respective Figure of the note \n");
    myBlue->SolveScaRho(1);
    myBlue->PrintScaRho(NamFil);

    // The independent scan
    myBlue->ResetInp();
    myBlue->FixInp();
    printf("... B_arXiv_1403_4427: The independent scan of RHO_ALL,");
    printf(" see Figure 7c,d \n");
    myBlue->SolveScaRho(0);

    myBlue->PrintScaRho(NamFil);

  }else if(Flag == 2){

    // Rho_Gro
    // The correlated scan
    myBlue->FixInp();
    printf("... B_arXiv_1403_4427: The correlated scan of RHO_Group,");
    printf(" see Figure 4a,b \n");
    printf("... B_arXiv_1403_4427: Group = 0/1/2/3 = EXP/LHC/TEV/COL\n");
    printf("... B_arXiv_1403_4427: The rightmost end of the histogram should be");
    printf(" be compared to the \n");
    printf("... B_arXiv_1403_4427: values at 40%%, 70%% and 100%% of the");
    printf(" respective Figure of the note \n");
    myBlue->SolveScaRho(1);
    myBlue->PrintScaRho(NamFil);

    // The independent scan
    myBlue->ResetInp();
    myBlue->FixInp();
    printf("... B_arXiv_1403_4427: The independent scan of RHO_Group,");
    printf(" see Figure 6 and 7 \n");
    printf("... B_arXiv_1403_4427: Group = 0/1/2/3 = EXP/LHC/TEV/COL =");
    printf(" 6a,b / 6e,f/ 6c,d /7a,b \n");
    myBlue->SolveScaRho(0);
    myBlue->PrintScaRho(NamFil);

  }else if(Flag == 3){

    // Arrays for systematic variations
    static const Int_t NumSys = 14;
    static const Int_t NumCor = 11;
    static const Int_t NumFla = NumFac*NumCor;
    Double_t ValSys[NumSys] = {0}, ValDef = 0;
    Double_t UncSys[NumSys] = {0}, UncDef = 0;
    Int_t    IndSys[NumCor] = {  2,   3,   3,   4,   5,   6,   6,   7,   8,   8,   9};
    Double_t RhoSys[NumCor] = {0.5, 0.5, 1.0, 1.0, 0.5, 0.0, 1.0, 0.5, 0.0, 1.0, 0.5};
    Double_t FlaSys[NumFla] = {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
				 1,   1,   1,   0,   1,   0,   0,   0,   0,   0,   1,
				 1,   1,   1,   0,   1,   0,   0,   0,   0,   0,   1,
				 1,   1,   1,   0,   1,   1,   1,   1,   1,   1,   1};
    //myBlue->PrintDouble(FlaSys,NumFac,NumCor,"  %3.1f");
    
    // The default result
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->GetResult(LocRes);
    myBlue->GetUncert(LocUnc);
    ValDef = LocRes->operator()(0,0);
    UncDef = LocUnc->operator()(0,0);

    // Now loop over correlation changes
    for(Int_t l = 0; l<NumCor; l++){
      myBlue->ResetInp();
      for(Int_t ll = 0; ll<NumFac; ll++){
	if(FlaSys[l + ll*NumCor] > 0.5)myBlue->SetRhoValUnc(IndSys[l], ll, RhoSys[l]);
      }
      if(l == 3)myBlue->SetRhoValUnc(IndSys[l], RhoSys[l]);
      myBlue->FixInp();
      myBlue->Solve();
      myBlue->GetResult(LocRes);
      myBlue->GetUncert(LocUnc);
      ValSys[l] = LocRes->operator()(0,0) - ValDef;
      UncSys[l] = LocUnc->operator()(0,0) - UncDef;
    }

    //-- Add or remove systematics -- This part needs a new Blue object
    // Save the MC systematics
    Double_t HadSys[NumEst] = {0};
    for(Int_t i = 0; i<NumEst; i++)HadSys[i] = XEst[i*(NumUnc+1) + 5 + 1];

    // Loop
    // l == 0 <==> Add CMS hadronisation k == 5 amounting to 0.56/0.76/0.93 
    //             for the l+j/dil/had channels i=8/9/10
    Double_t AddHad[3] = {0.58, 0.76, 0.93};

    // l == 1 Remove CDF, D0 and ATLAS hadronisation
    // Tevatron: Subtract HAD, numbers are from G. Cortiana priv comm  i=0-5
    //    ATLAS: Subtract HAD, numbers are from the original notes Refs.[14,15],
    //            0.27/0.44 for i=6/7
    //      CMS: Keep original numbers 
    Double_t TEVHad[6] = {0.56, 0.32, 0.48, 0.36, 0.58, 0.50};
    Double_t ATLHad[2] = {0.27, 0.44};
	  
    // l == 2 <==> CMS alternative categorisation, change bjes k == 4 to 0.18
    //             for all channels i=8/9/10 keep hadronisation from l==0
    Double_t SetBJe = 0.18;

    Int_t idx = NumCor;
    for(Int_t l = 0; l<3; l++){
      Blue *myBlueI = new Blue(NumEst, NumUnc);
      if(l == 0){
	// CMS Add
	for(Int_t i = 8; i<11; i++){
	  XEst[i*(NumUnc+1) + 5 + 1] = 
	    TMath::Sqrt(HadSys[i]*HadSys[i] + AddHad[i-8]*AddHad[i-8]);
	}
      }
      if(l == 1){
	// Tevatron subtract
	for(Int_t i = 0; i< 6; i++){
	  XEst[i*(NumUnc+1) + 5 + 1] =
	    TMath::Sqrt(HadSys[i]*HadSys[i] - TEVHad[i]*TEVHad[i]);
	}
	// ATLAS subtract
	for(Int_t i = 6; i< 8; i++){
	  XEst[i*(NumUnc+1) + 5 + 1] =
	    TMath::Sqrt(HadSys[i]*HadSys[i] - ATLHad[i-6]*ATLHad[i-6]);
	}
	// CMS set to original
	for(Int_t i = 8; i<11; i++)XEst[i*(NumUnc+1) + 5 + 1] = HadSys[i];
      }
      if(l == 2){
	// Tevatron + ATLAS set to original
	for(Int_t i = 0; i< 8; i++)XEst[i*(NumUnc+1) + 5 + 1] = HadSys[i];
	
	// CMS add Had and set BJes
	for(Int_t i = 8; i<11; i++){
	  XEst[i*(NumUnc+1) + 5 + 1] = 
	    TMath::Sqrt(HadSys[i]*HadSys[i] + AddHad[i-8]*AddHad[i-8]);
	}
	for(Int_t i = 8; i<11; i++)XEst[i*(NumUnc+1) + 4 + 1] = SetBJe;
      }

      // Fill estimates
      ind = 0;
      for(Int_t i = 0; i<NumEst; i++){
	myBlueI->FillEst(i,&XEst[ind]);
	ind = ind + NumUnc + 1;
      }
      
      // Fill correlations
      for(Int_t k = 0; k<NumUnc; k++){
	if(k == 12){
	  myBlueI->FillCor(k,&CoBGMC[0]);
	}else{
	  ind = RhoInd[k];    
	  FilMat(RhoCDF[ind], RhoDZE[ind], RhoATL[ind], RhoCMS[ind],
		 RhoLHC[ind], RhoTEV[ind], RhoATE[ind], RhoCTE[ind], &RhoMat[0]);
	  myBlueI->FillCor(k,&RhoMat[0]);
	}
      }
      // Now Solve and get the result
      myBlueI->FixInp();
      myBlueI->Solve();
      myBlueI->GetResult(LocRes);
      myBlueI->GetUncert(LocUnc);
      ValSys[idx] = LocRes->operator()(0,0) - ValDef;
      UncSys[idx] = LocUnc->operator()(0,0) - UncDef;
      idx = idx + 1;
      delete myBlueI;
      myBlueI = NULL;
    }
    //-- End

    // Print what we got for the systematics
    printf("... B_arXiv_1403_4427: The systematics variations,");
    printf(" see Figures 3c,d \n");
    for(Int_t l = 0; l<NumSys; l++){
      printf("... B_arXiv_1403_4427: Syst %2i:",l);
      printf(" Delta(mtop) = %+5.3f Delta(Sigma_mtop) = %+5.3f \n", 
	     ValSys[l], UncSys[l]);
    }

  }else if(Flag == 4){

    // The correlated  result
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1403_4427: For the results per channel see Table 7");
    myBlue->PrintResult();
    iok = myBlue->GetRhoRes(LocRhoRes);
    if(iok == 1)myBlue->PrintMatrix(LocRhoRes,"%5.2f");

    myBlue->PrintCompatObs();
    myBlue->LatexResult(NamCha);
    for(Int_t n = 0; n<NumObs; n++)myBlue->DisplayResult(n,NamCha);
    for(Int_t n = 0; n<NumObs; n++)myBlue->InspectLike(n,NamCha);

    // The individual Results
    // LPJ
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumDIL; l++)myBlue->SetInActiveEst(IndDIL[l]);
    for(Int_t l = 0; l<NumHAD; l++)myBlue->SetInActiveEst(IndHAD[l]);
    for(Int_t l = 0; l<NumMET; l++)myBlue->SetInActiveEst(IndMET[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

    // DIL
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumLPJ; l++)myBlue->SetInActiveEst(IndLPJ[l]);
    for(Int_t l = 0; l<NumDIL; l++)myBlue->SetActiveEst(IndDIL[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

    // HAD
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumDIL; l++)myBlue->SetInActiveEst(IndDIL[l]);
    for(Int_t l = 0; l<NumHAD; l++)myBlue->SetActiveEst(IndHAD[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

    // MET
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumHAD; l++)myBlue->SetInActiveEst(IndHAD[l]);
    for(Int_t l = 0; l<NumMET; l++)myBlue->SetActiveEst(IndMET[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

  }else if(Flag == 5){

    // The correlated  result
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1403_4427: For the results per experiment see Table 8");
    myBlue->PrintResult();
    iok = myBlue->GetRhoRes(LocRhoRes);
    if(iok == 1)myBlue->PrintMatrix(LocRhoRes,"%5.2f");
    myBlue->PrintCompatObs();
    myBlue->LatexResult(NamExp);
    for(Int_t n = 0; n<NumObs; n++)myBlue->DisplayResult(n,NamExp);

    // The individual Results
    // CDF
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumDZE; l++)myBlue->SetInActiveEst(IndDZE[l]);
    for(Int_t l = 0; l<NumATL; l++)myBlue->SetInActiveEst(IndATL[l]);
    for(Int_t l = 0; l<NumCMS; l++)myBlue->SetInActiveEst(IndCMS[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

    // DZE
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumCDF; l++)myBlue->SetInActiveEst(IndCDF[l]);
    for(Int_t l = 0; l<NumDZE; l++)myBlue->SetActiveEst(IndDZE[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

    // ATL
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumDZE; l++)myBlue->SetInActiveEst(IndDZE[l]);
    for(Int_t l = 0; l<NumATL; l++)myBlue->SetActiveEst(IndATL[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

    // CMS
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumATL; l++)myBlue->SetInActiveEst(IndATL[l]);
    for(Int_t l = 0; l<NumCMS; l++)myBlue->SetActiveEst(IndCMS[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

  }else if(Flag == 6){

    // The correlated  result
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1403_4427: For the results per collider see Table 9");
    myBlue->PrintResult();
    iok = myBlue->GetRhoRes(LocRhoRes);
    if(iok == 1)myBlue->PrintMatrix(LocRhoRes,"%5.2f");
    myBlue->PrintCompatObs();
    myBlue->LatexResult(NamCol);
    for(Int_t n = 0; n<NumObs; n++)myBlue->DisplayResult(n,NamCol);

    // The individual Results
    // TEV
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumATL; l++)myBlue->SetInActiveEst(IndATL[l]);
    for(Int_t l = 0; l<NumCMS; l++)myBlue->SetInActiveEst(IndCMS[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

    // LHC
    myBlue->ReleaseInp();
    for(Int_t l = 0; l<NumCDF; l++)myBlue->SetInActiveEst(IndCDF[l]);
    for(Int_t l = 0; l<NumDZE; l++)myBlue->SetInActiveEst(IndDZE[l]);
    for(Int_t l = 0; l<NumATL; l++)myBlue->SetActiveEst(IndATL[l]);
    for(Int_t l = 0; l<NumCMS; l++)myBlue->SetActiveEst(IndCMS[l]);
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();

  }
  // Delete objects
  delete myBlue;  myBlue = NULL;
  LocRho->Delete(); LocRho = NULL;
  LocRhoRes->Delete(); LocRhoRes = NULL;
  LocRes->Delete(); LocRes = NULL;
  LocUnc->Delete(); LocUnc = NULL;
  return;
}

//------------------------------------------------------------------------
// Utility to fill a correlation matrix for any pair of RhoExp and RhoLHC
//------------------------------------------------------------------------

void FilMat(Double_t CDF, Double_t DZE, Double_t ATL, Double_t CMS,
	    Double_t LHC, Double_t TEV, Double_t ATE, Double_t CTE,
	    Double_t *M){
  //       0           1           2           3           4           5           6           7           8           9          10
  M[  0]=  1; M[  1]=CDF; M[  2]=CDF; M[  3]=CDF; M[  4]=TEV; M[  5]=TEV; M[  6]=ATE; M[  7]=ATE; M[  8]=CTE; M[  9]=CTE; M[ 10]=CTE;
  M[ 11]=CDF; M[ 12]=  1; M[ 13]=CDF; M[ 14]=CDF; M[ 15]=TEV; M[ 16]=TEV; M[ 17]=ATE; M[ 18]=ATE; M[ 19]=CTE; M[ 20]=CTE; M[ 21]=CTE;
  M[ 22]=CDF; M[ 23]=CDF; M[ 24]=  1; M[ 25]=CDF; M[ 26]=TEV; M[ 27]=TEV; M[ 28]=ATE; M[ 29]=ATE; M[ 30]=CTE; M[ 31]=CTE; M[ 32]=CTE;
  M[ 33]=CDF; M[ 34]=CDF; M[ 35]=CDF; M[ 36]=  1; M[ 37]=TEV; M[ 38]=TEV; M[ 39]=ATE; M[ 40]=ATE; M[ 41]=CTE; M[ 42]=CTE; M[ 43]=CTE;
  M[ 44]=TEV; M[ 45]=TEV; M[ 46]=TEV; M[ 47]=TEV; M[ 48]=  1; M[ 49]=DZE; M[ 50]=ATE; M[ 51]=ATE; M[ 52]=CTE; M[ 53]=CTE; M[ 54]=CTE;
  M[ 55]=TEV; M[ 56]=TEV; M[ 57]=TEV; M[ 58]=TEV; M[ 59]=DZE; M[ 60]=  1; M[ 61]=ATE; M[ 62]=ATE; M[ 63]=CTE; M[ 64]=CTE; M[ 65]=CTE;
  M[ 66]=ATE; M[ 67]=ATE; M[ 68]=ATE; M[ 69]=ATE; M[ 70]=ATE; M[ 71]=ATE; M[ 72]=  1; M[ 73]=ATL; M[ 74]=LHC; M[ 75]=LHC; M[ 76]=LHC;
  M[ 77]=ATE; M[ 78]=ATE; M[ 79]=ATE; M[ 80]=ATE; M[ 81]=ATE; M[ 82]=ATE; M[ 83]=ATL; M[ 84]=  1; M[ 85]=LHC; M[ 86]=LHC; M[ 87]=LHC;
  M[ 88]=CTE; M[ 89]=CTE; M[ 90]=CTE; M[ 91]=CTE; M[ 92]=CTE; M[ 93]=CTE; M[ 94]=LHC; M[ 95]=LHC; M[ 96]=  1; M[ 97]=CMS; M[ 98]=CMS;
  M[ 99]=CTE; M[100]=CTE; M[101]=CTE; M[102]=CTE; M[103]=CTE; M[104]=CTE; M[105]=LHC; M[106]=LHC; M[107]=CMS; M[108]=  1; M[109]=CMS;
  M[110]=CTE; M[111]=CTE; M[112]=CTE; M[113]=CTE; M[114]=CTE; M[115]=CTE; M[116]=LHC; M[117]=LHC; M[118]=CMS; M[119]=CMS; M[120]=  1;
  //       0           1           2           3           4           5           6           7           8           9          10
  return;
};

