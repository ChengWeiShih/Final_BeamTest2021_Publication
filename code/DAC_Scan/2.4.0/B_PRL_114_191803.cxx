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
//------------------------------------------------------------------------------
// Utility prototypes
//------------------------------------------------------------------------------
    void FilMat(const Double_t RE, const Double_t RL, TMatrixD *const M);
Double_t DifUnc(const Double_t ful, const Double_t par);
Double_t TotUnc(const Double_t sta, const Double_t sys);
Double_t BetVal(const Double_t xva, const Double_t xv1, Double_t xv2);
Double_t RhoVal(const Double_t bet, const Double_t zva);

void B_PRL_114_191803(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The mass using all four inputs
  //  1: The mass in the GG channel 
  //  2: The mass in the LL channel 
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  4;
  static const Int_t NumUnc = 14;
  static const Int_t NumObs =  1;

  // Set the names  
  const TString NamEst[NumEst] = {"M_{ATL,GG}", "M_{ATL,4L}", 
				  "M_{CMS,GG}", "M_{CMS,4L}"};
  const TString NamUnc[NumUnc] = {
    //     0          1          2          3          4          5         6
    "   Stat", " NonLin", " Matter", "LongRes", " LatRes", " GamRes", "Vertex",
    //     7          8          9         10         11         12        13
    " ZeeCal", " CMSElE", " MuoMom", " GGBack ", "  Lumi", " AddExp ", " Theo"
  };
  TString NamObs[NumObs] = {"M_{LHC}"};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_PRL_114_191803: ----------------------------\n");
    printf("... B_PRL_114_191803: The combined mass, Flag = %2i \n",Flag);
  }else if(Flag == 1){
    printf("... B_PRL_114_191803: -------------------------------------\n");
    printf("... B_PRL_114_191803: The mass in the GG channel, Flag = %2i \n",
	   Flag);
    NamObs[0] = "M_{LHC,GG}";
  }else if(Flag == 2){
    printf("... B_PRL_114_191803: -------------------------------------\n");
    printf("... B_PRL_114_191803: The mass in the 4L channel, Flag = %2i \n",
	   Flag);
    NamObs[0] = "M_{LHC,4L}";
  }else{
    printf("... PRL_114_191803: Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  // Estimates 0=ATL,GG 1=ATL,4L, 2=CMS,GG, 3=CMS,4L
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  const Double_t XEst[LenXEst] = {
    // Observed uncertainties
    //     0      1      2       3      4      5     6 
    //  Stat NonLin Matter LongRes LatRes GamRes Vertex
    //     7      8      9      10     11     12     13
    //ZeeCal CMSElE MuoMom  GGBack   Lumi AddExp   Theo
    126.02,
    0.43,      0.14,  0.15,  0.12,   0.09,  0.03,  0.05,
    0.05,      0.00,  0.00,  0.04,   0.01,  0.03,  0.01,
    124.51,
    0.52,      0.00,  0.00,  0.00,   0.00,  0.00,  0.00,
    0.03,      0.00,  0.03,  0.00,   0.01,  0.01,  0.01,
    124.70, 
    0.31,      0.10,  0.07,  0.02,   0.06,  0.01,  0.00, 
    0.05,      0.00,  0.00,  0.00,   0.01,  0.02,  0.02,
    125.59, 
    0.42,      0.00,  0.00,  0.00,   0.00,  0.00,  0.00,
    0.00,      0.12,  0.11,  0.00,   0.01,  0.01,  0.01
  };
  const Double_t XExp[LenXEst] = {
    // Expected uncertainties
    //     0      1      2       3      4      5     6 
    //  Stat NonLin Matter LongRes LatRes GamRes Vertex
    //     7      8      9      10     11     12     13
    //ZeeCal CMSElE MuoMom  GGBack   Lumi AddExp   Theo
    126.02,
    0.45,      0.16,  0.13,   0.13,  0.08,  0.01,  0.05,
    0.04,      0.00,  0.00,   0.03,  0.01,  0.01,  0.01,
    124.51,
    0.66,      0.00,  0.00,   0.00,  0.00,  0.00,  0.00,
    0.02,      0.00,  0.04,   0.00,  0.01,  0.01,  0.01,
    124.70,
    0.32,      0.13,  0.07,   0.01,  0.06,  0.01,  0.00,
    0.05,      0.00,  0.00,   0.00,  0.01,  0.01,  0.01,
    125.59,
    0.57,      0.00,  0.00,   0.00,  0.00,  0.00,  0.00,
    0.00,      0.09,  0.10,   0.00,  0.01,  0.01,  0.01
  };

  // Estimates with stat and full syst uncertainty only
  static const Int_t RedUnc = 2;
  static const Int_t LenYEst = NumEst * (RedUnc+1);
  const Double_t YEst[LenYEst] = {
    //         0     1
    //      Stat  Syst
    126.02, 0.43, 0.27,
    124.51, 0.52, 0.04,
    124.70, 0.31, 0.15,
    125.59, 0.42, 0.16
  };

  // The combined values from the paper
  // 0=all, 1=GG, 2=4L
  static const Int_t NumPap = 3;
  const Double_t XPap[LenYEst] = {
    //         0     1
    //      Stat  Syst
    125.09, 0.21, 0.11,
    125.07, 0.25, 0.14,
    125.15, 0.37, 0.15
  };

  // The correlations within experiments
  const Double_t RhoExp[NumUnc] = {
    //     0      1      2       3      4      5     6 
    //  Stat NonLin Matter LongRes LatRes GamRes Vertex
    //     7      8      9      10     11     12     13
    //ZeeCal CMSElE MuoMom  GGBack   Lumi AddExp   Theo
    0.00,      0.00,  0.00,   0.00,  0.00,  0.00,  0.00,
    1.00,      0.00,  0.00,   0.00,  1.00,  1.00,  1.00};

  // The correlations across experiments
  const Double_t RhoLHC[NumUnc] = {
    //     0      1      2       3      4      5     6 
    //  Stat NonLin Matter LongRes LatRes GamRes Vertex
    //     7      8      9      10     11     12     13
    //ZeeCal CMSElE MuoMom  GGBack   Lumi AddExp   Theo
    0.00,      0.00,  0.00,   0.00,  0.00,  0.00,  0.00,
    0.00,      0.00,  0.00,   0.00,  1.00,  0.00,  1.00};
  
  //-- Local Structures for Blue input
  //- Estimates
  TMatrixD* InpEst = new TMatrixD(NumEst,NumUnc+1,&XEst[0]);
  //- Estimates for stat + syst only
  TMatrixD* RedEst = new TMatrixD(NumEst,RedUnc+1,&YEst[0]);
  //- Results of the paper
  TMatrixD* PapRes = new TMatrixD(NumPap,RedUnc+1,&XPap[0]);
  //- Correlations from the paper combination of pairs
  TMatrixD* PapRho = new TMatrixD(NumEst,NumEst);

  //- Correlation per source filled dynamically
  TMatrixD* RhoSou = new TMatrixD(NumEst,NumEst);
  //-- End

  //-- Local Structures for Blue output
  //- TMatrices
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRes = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc = new TMatrixD(NumObs,1);
  //-- End

  // Define formats for figures and latex file
  const TString FilBas = "B_PRL_114_191803";
  TString FilNam = "To be filled later";
  char Buffer[100];
  sprintf(Buffer,"%1i", Flag);
  FilNam = &Buffer[0];
  FilNam = FilBas + "_" + FilNam;

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%+5.3f";
  const TString ForRho = "%7.4f";
  const TString ForPul = ForUnc;
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
  myBlue->FillEst(InpEst);

  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++){
    FilMat(RhoExp[k], RhoLHC[k], RhoSou);
    myBlue->FillCor(k, RhoSou);
  }

  // Fill correlation matrix for same channel pairs from paper info
  // RedEst(ii) = more precise <==> x1
  Double_t ti = 0., tj = 0., zva = 0.;
  Int_t    ii = 0,  jj = 0,  ir = 1;
  PapRho->UnitMatrix();
  for(Int_t i = 0; i<2; i++){
    for(Int_t j = i+2; j<NumEst; j=j+7){
      ii = i;
      jj = j;
      ti = TotUnc(RedEst->operator()(i,1), RedEst->operator()(i,2));
      tj = TotUnc(RedEst->operator()(j,1), RedEst->operator()(j,2));
      zva = tj/ti;
      if(ti > tj){
	ii=j;
	jj=i;
	zva = 1/zva;
      }
      PapRho->operator()(ii,jj) = RhoVal(BetVal(PapRes->operator()(ir,0),
						RedEst->operator()(ii,0),
						RedEst->operator()(jj,0)),
					 zva);
      PapRho->operator()(jj,ii) = PapRho->operator()(ii,jj);
    }
    ir = 2;
  }

  // Solve according to Flag
  Int_t MinRow = 0;
  Int_t MaxRow = NumEst;
  Int_t MinCol = MinRow;
  Int_t MaxCol = MaxRow;
  const Int_t n = Flag;

  // Selection of inputs
  if(Flag == 0){
    // The 4 channel result
  }else if(Flag == 1){
    // The GG result
    myBlue->SetInActiveEst(1);
    myBlue->SetInActiveEst(3);
  }else if(Flag == 2){
    // The 4L result
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(2);
  }

  // Perform the combination
  myBlue->FixInp();
  myBlue->Solve();
  printf("... B_PRL_114_191803: For the estimates and the result");
  printf(" see Table 1 \n");
  myBlue->PrintEst();
  myBlue->PrintResult();

  // Plot compatibility
  if(Flag == 0)myBlue->PrintCompatEst(FilNam);

  // Calculate differences Paper - Blue
  myBlue->GetResult(LocRes);
  myBlue->GetUncert(LocUnc);
  printf("... B_PRL_114_191803: \n");
  printf("... B_PRL_114_191803: ----------   The difference of");
  printf(" the two approaches is: -------- \n");

  //- Blue
  printf("... B_PRL_114_191803:    The Blue Combination: ");
  printf("%s = %7.2f (+-%6.2f +-%6.2f) = +-%6.2f\n",
	 NamObs[0].Data(), LocRes->operator()(0,0), LocRes->operator()(0,1), 
	 DifUnc(LocUnc->operator()(0,0), LocRes->operator()(0,1)),
	 LocUnc->operator()(0,0));
  
  //- Paper
  printf("... B_PRL_114_191803:        The paper result: ");
  printf("%s = %7.2f (+-%6.2f +-%6.2f) = +-%6.2f\n",
	 NamObs[0].Data(), PapRes->operator()(n,0), 
	 PapRes->operator()(n,1), PapRes->operator()(n,2), 
	 TotUnc(PapRes->operator()(n,1), PapRes->operator()(n,2)));
  
  //- Difference
  printf("... B_PRL_114_191803: Difference = Paper-Blue: ");
  printf("%s = %7.2f (  %+6.2f +-%+6.2f) =   %+6.2f\n",
	 NamObs[0].Data(), 
	 PapRes->operator()(n,0) - LocRes->operator()(0,0), 
	 PapRes->operator()(n,1) - LocRes->operator()(0,1), 
	 PapRes->operator()(n,2) -
	 DifUnc(LocUnc->operator()(0,0), LocRes->operator()(0,1)),
	 TotUnc(PapRes->operator()(n,1), PapRes->operator()(n,2))-
	 LocUnc->operator()(0,0));
  printf("... B_PRL_114_191803: \n");
  
  // Put results into latex files
  myBlue->LatexResult(FilNam);
  
  // Investigate the estimator correlations
  if(Flag == 0){
    printf("... B_PRL_114_191803: The Blue  correlation of");
    printf(" the estimates is:\n");
    myBlue->GetRho(LocRho);
    LocRho->operator*=(100.);
    myBlue->PrintMatrix(LocRho,myBlue->GetActEst(),
			myBlue->GetActEst(),"%5.1f%%");

    // Display the hypothetical combination of pairs
    myBlue->DisplayPair(0, 2, FilBas, 121., 128., 0.0, 0.8);    
    myBlue->DisplayPair(1, 3, FilBas);

    // The partly filled correlation matrix from the paper
    PapRho->operator*=(100.);
    printf("... B_PRL_114_191803: The estimator correlations of the paper:");
    printf(" Only the pairs 0,2\n");
    printf("... B_PRL_114_191803: and 1,3 can be obtained from");
    printf(" the pairwise results in the paper.\n");
    myBlue->PrintMatrix(PapRho, MinRow, MaxRow, MinCol, MaxCol,"%5.1f%%");

    // Summary
    printf("... B_PRL_114_191803: Given the description in the paper, the");
    printf(" apparent estimator correlations calculated\n");
    printf("... B_PRL_114_191803: from the the pairwise combinations in");
    printf(" the paper are surprisingly large. \n");
    printf("... B_PRL_114_191803: \n");
  }

  // Delete Object, clean up and return 
  delete myBlue; myBlue = NULL;

  InpEst->Delete(); InpEst = NULL;
  RedEst->Delete(); RedEst = NULL;
  PapRes->Delete(); PapRes = NULL;

  PapRho->Delete(); PapRho = NULL;
  RhoSou->Delete(); RhoSou = NULL;

  LocRho->Delete(); LocRho = NULL;
  LocRes->Delete(); LocRes = NULL;
  LocUnc->Delete(); LocUnc = NULL;

  return;
};

//------------------------------------------------------------------------------
// Fill a correlation matrix for any pair of RhoExp and RhoLHC
//------------------------------------------------------------------------------
void FilMat(const Double_t RE, const Double_t RL, TMatrixD *const M){

   // Check input
   const Int_t ir = M->GetNrows();
   const Int_t ic = M->GetNcols();
   if(ir != ic){
     printf("... FilMat(): Should be a square matrix");
     printf(" NRow = %i, NCol = %i \n",ir, ic);
     return;
   }

   // Fill the diagonal
   for(Int_t i = 0; i<ir; i++)M->operator()(i,i) = 1.0;

   // Fill the upper half
   M->operator()(0,1) = RE;
   M->operator()(0,2) = RL;
   M->operator()(0,3) = RL;
   M->operator()(1,2) = RL;
   M->operator()(1,3) = RL;
   M->operator()(2,3) = RE;

   // Copy to the lower half
   for(Int_t i = 0; i<ir; i++){
     for(Int_t j = i+1; j<ir; j++){
       M->operator()(j,i) = M->operator()(i,j);
     }
   }
   return;
};

//------------------------------------------------------------------------------
// Calculate total uncertainty from stat and syst
//------------------------------------------------------------------------------
Double_t TotUnc(const Double_t sta, const Double_t sys){
  // Calculate sum in quadrature
  return TMath::Sqrt(sta*sta + sys*sys);
};

//------------------------------------------------------------------------------
// Calculate difference from full and par=stat/syst
//------------------------------------------------------------------------------
Double_t DifUnc(const Double_t ful, const Double_t par){
  // Calculate difference in quadrature
  return TMath::Sqrt(ful*ful - par*par);
};

//------------------------------------------------------------------------------
// Calculate beta from x, x1 and x2
//------------------------------------------------------------------------------
Double_t BetVal(const Double_t xva, const Double_t xv1, Double_t xv2){
  // Calculate beta from manual Eq.2 first vs. second term
  return (xva - xv1) / (xv2 - xv1);
};

//------------------------------------------------------------------------
// Calculate rho from beta and z
//------------------------------------------------------------------------
Double_t RhoVal(const Double_t bet, const Double_t zva){
  // Calculate rho from manual Eq.2 first vs. third term
  return (1. - bet*(1 + zva*zva)) / (zva*(1. - 2.*bet));
};
