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
    void FilMat(const Double_t RDif, const Double_t RSsam, TMatrixD *const M);
Double_t SysUnc(const Double_t ful, const Double_t sta);

void B_PRD93_072004(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: CMS combination with --reduced-- correlations
  //  1: CMS combination with --nominal-- correlations
  //  2: CMS combination per decay channel with --reduced-- correlations
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 11;
  static const Int_t NumUnc = 27;
  static const Int_t MaxObs =  3;
  Int_t NumObs =  1;

  const TString NamEst[NumEst] = {
    //      0
    "10-dil-1",
    //      1           2           3
    "11-l+j-2", "11-had-1", "11-dil-1",
    //      4           5           6
    "12-l+j-2", "12-l+j-1", "12-l+j-h",
    //      7           8           9
    "12-had-2", "12-had-1", "12-had-h",
    //     10
    "12-dil-1"};

  const TString NamUnc[NumUnc] = {
    //     0          1          2          3          4          5    
    "   Stat", " FitCal", "JEC-Int", "JEC-Ins", "JEC-Unp", "JEC-Uwp", 
    //     6          7          8          9         10         11
    "   Lept", " ETmiss", "    JER", "  b-tag", "   Trig", " Pileup",
    //    18         19         20         21         22         23
    "   BGDT", "   BGMC", "JEC-lig", "JEC-cha", "JEC-bot", "JEC-glu",
    //    12         13         14         15         16         17
    " b-frag", " b-slep", "    PDF", "  muR/F","   MEPS", "     MC", 
    //    24         25         26
    " top-pt", "     UE", "     CR"};
  
  TString NamObs[MaxObs] = {"mas-top", "mas-had", "mas-dil"};

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {NumEst*0};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_PRD93_072004: ------------------------------");
    printf("------------------------------ \n");
    printf("... B_PRD93_072004: The various scenarios of the 2014");
    printf(" CMS combination, Flag = %2i \n", Flag);
  }else if(Flag == 1){
    printf("... B_PRD93_072004: ----------------------------------");
    printf("-------------------------------- \n");
    printf("... B_PRD93_072004: The 2014 CMS combination, but with");
    printf(" nominal correlations, Flag = %2i \n", Flag);
  }else if(Flag == 2){
    printf("... B_PRD93_072004: ---------------------------------");
    printf("--------------------------- \n");
    printf("... B_PRD93_072004: The combination for l+jets,");
    printf(" dilepton and all-jets, Flag = %2i \n", Flag);
    NumObs = 3;
    IWhichObs[ 0] = 2;
    IWhichObs[ 1] = 0;
    IWhichObs[ 2] = 1;
    IWhichObs[ 3] = 2;
    IWhichObs[ 4] = 0;
    IWhichObs[ 5] = 0;
    IWhichObs[ 6] = 0;
    IWhichObs[ 7] = 1;
    IWhichObs[ 8] = 1;
    IWhichObs[ 9] = 1;
    IWhichObs[10] = 2;
    NamObs[0] = "mas-l+j";
  }else{
    printf("... B_PRD93_072004: ------------------------ \n");
    printf("... PRD93_072004: Not implemented Flag = %2i \n", Flag);
    return;
  }

  // Estimates 0-10
  //  0 ==   CMS_2010   di-lepton
  //  1 ==   CMS_2011 lepton+jets
  //  2 ==   CMS_2011    all-jets
  //  3 ==   CMS_2011   di-lepton
  //  4 ==   CMS_2012 lepton+jets 2D
  //  5 ==   CMS_2012 lepton+jets 1D
  //  6 ==   CMS_2012 lepton+jets Hybrid
  //  7 ==   CMS_2012    all-jets 2D
  //  8 ==   CMS_2012    all-jets 1D
  //  9 ==   CMS_2012    all-jets Hybrid
  // 10 ==   CMS_2012   di-lepton

  // Estimates i=0-3 have unclear input. There is no table in this paper. This
  // input had to be constructed from the Refs given, and from PAS-TOP-14-015.
  // Choices were made as follows:
  //
  // i=0: Paper says [54] AMWT, but quotes 4.6 instead of 4.9 as stats
  // i=0: and 175.5 instead of 175.8 as value.
  // i=0: In contrast PAS-TOP-14-015 uses the combination Klnb and AMWT
  // i=0: Assume PAS-TOP-14-015 is correct.
  //
  // i=1: Paper says [55] 2D.
  // i=1: k=2-5+14: [55] says JES = 0.28
  // i=1: PAS-TOP-14-015 = 2: 0.01, 3: 0.02, 4: 0.24, 5:n/a, 14: 0.11 = 0.27
  // i=1: take split from PAS-TOP-14-015
  // i=1: k=23: [55] does not have this ==> use PAS-TOP-14-015
  //
  // i=2: Paper says [48] 1D.
  // i=2: k=2-5+14: [48] says JES = 0.97
  // i=2: PAS-TOP-14-015 = 2: 0.08, 3: 0.35, 4: 0.69, 5:n/a, 14: 0.58 = 0.97
  // i=2: take split from PAS-TOP-14-015
  // i=2: k=23: [48] does not have this ==> use PAS-TOP-14-015
  //
  // i=3: Paper says [43] AMWT
  // i=3: k=2-5+14: [43] says =90/-0.97 
  // i=3: PAS-TOP-14-015 has the same split as for i=2, take this.
  // i=3: k=16: [43] says +76/-66, PAS-TOP-14-015 uses max of both, keep this.
  // i=3: k=23: [43] does not have this ==> use PAS-TOP-14-015

  // Add uncorrelated JES for Est=0
  Double_t UncU = TMath::Sqrt(1.48*1.48+3.28*3.28);

  // The estimates
  //           0     1     2     3     4     5     6     7     8     9    10
  //    11    12    13    14    15    16    17    18    19    20    21    22
  //    23    24    25    26
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  const Double_t XEst[LenXEst] = {
    // 2010 dil == 0
    175.50, 4.60, 0.20, 0.17, 0.76, UncU, 0.00, 0.30, 0.10, 0.50, 0.40, 0.00,
    1.00,   0.00, 0.10, 1.21, 0.00, 0.90, 0.00, 0.00, 0.00, 0.50, 0.63, 0.70,
    0.50,   0.00, 1.30, 0.00,
    // 2011 l+j, had, dil == 1, 2, 3
    173.49, 0.43, 0.06, 0.01, 0.02, 0.24, 0.00, 0.02, 0.06, 0.23, 0.12, 0.00,
    0.07,   0.00, 0.13, 0.11, 0.00, 0.61, 0.00, 0.00, 0.00, 0.07, 0.24, 0.18,
    0.02,   0.00, 0.15, 0.54, 
    173.49, 0.69, 0.13, 0.08, 0.35, 0.69, 0.00, 0.00, 0.00, 0.15, 0.06, 0.24,
    0.06,   0.13, 0.00, 0.58, 0.00, 0.49, 0.00, 0.00, 0.00, 0.06, 0.22, 0.24,
    0.19  , 0.00, 0.20, 0.15, 
    172.50, 0.43, 0.40, 0.08, 0.35, 0.69, 0.00, 0.14, 0.12, 0.14, 0.09, 0.00,
    0.11,   0.00, 0.05, 0.58, 0.00, 0.76, 0.00, 0.00, 0.00, 0.09, 0.55, 0.19,
    0.04,   0.00, 0.05, 0.13,
    // 2012 l+j 2, 1, h == 4, 5, 6
    172.14, 0.20, 0.04, 0.01, 0.01, 0.09, 0.06, 0.01, 0.04, 0.11, 0.06, 0.00,
    0.12,   0.00, 0.05, 0.11, 0.03, 0.32, 0.22, 0.06, 0.16, 0.09, 0.17, 0.11,
    0.07,   0.16, 0.15, 0.11, 
    172.56, 0.12, 0.04, 0.02, 0.24, 0.26, 0.11, 0.01, 0.03, 0.05, 0.04, 0.00,
    0.05,   0.00, 0.01, 0.02, 0.01, 0.31, 0.05, 0.06, 0.15, 0.06, 0.24, 0.07,
    0.16,   0.11, 0.07, 0.09, 
    172.35, 0.16, 0.04, 0.01, 0.12, 0.10, 0.04, 0.01, 0.04, 0.03, 0.06, 0.00,
    0.04,   0.00, 0.03, 0.05, 0.01, 0.32, 0.08, 0.01, 0.16, 0.04, 0.09, 0.03,
    0.12,   0.02, 0.08, 0.01, 
    // 2012 had 2, 1, h == 7, 8, 9
    171.64, 0.32, 0.06, 0.01, 0.01, 0.06, 0.04, 0.00, 0.00, 0.10, 0.02, 0.04,
    0.09,   0.61, 0.00, 0.10, 0.03, 0.30, 0.17, 0.08, 0.14, 0.06, 0.29, 0.18,
    0.04,   0.04, 0.27, 0.35,
    172.46, 0.23, 0.06, 0.02, 0.23, 0.19, 0.08, 0.00, 0.00, 0.03, 0.01, 0.01,
    0.02,   0.14, 0.00, 0.02, 0.01, 0.29, 0.02, 0.03, 0.13, 0.03, 0.19, 0.12,
    0.18,   0.08, 0.13, 0.14,
    172.32, 0.25, 0.06, 0.02, 0.19, 0.16, 0.06, 0.00, 0.00, 0.02, 0.02, 0.01,
    0.01,   0.20, 0.00, 0.00, 0.01, 0.29, 0.02, 0.04, 0.13, 0.03, 0.12, 0.13,
    0.16,   0.06, 0.14, 0.16,
    // 2012 dil == 10
    172.82, 0.19, 0.03, 0.03, 0.24, 0.28, 0.12, 0.12, 0.06, 0.06, 0.04, 0.00,
    0.04,   0.00, 0.02, 0.02, 0.02, 0.34, 0.06, 0.69, 0.17, 0.16, 0.75, 0.12,
    0.24, 0.25, 0.04, 0.11
   };

  // Correlations:
  //--------------
  // rho_year <==> pairs in -different- decay channels but same year -> Rho_Dif
  // rho_chan <==> pairs in - the same- decay channels but diff year -> Rho_Sam
  // Rho_Dif = 0, Rho_Sam = 0 for: 0-1, 10, 12
  // Rho_Dif = 1, Rho_Sam = 1 rho: 2-3, 6-9, 13, 14-26
  // Rho_Dif = 1, Rho_Sam = 1 for: 4, 11

  // rho_year 
  Double_t RhoDif[NumUnc] = {
  //  0    1    2    3    4    5?   6    7    8    9   10
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
  // 11   12   13   14   15   16   17   18   19   20   21
    1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
  // 22   23   24   25   26
    1.0, 1.0, 1.0, 1.0, 1.0};

  // rho_chan
  Double_t RhoSam[NumUnc] = {
  //  0    1    2    3    4    5?   6    7    8    9   10
    0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0,
  // 11   12   13   14   15   16   17   18   19   20   21
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
  // 22   23   24   25   26
    1.0, 1.0, 1.0, 1.0, 1.0};

  //-- Local structures for BLUE input
  // TMatrices
  TMatrixD* RhoSou = new TMatrixD(NumEst,NumEst);
  //-- End

  //-- Local structures for BLUE output
  // TMatrices
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRes = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc = new TMatrixD(NumObs,1);

  // Arrays for the Numijk = 7 results of Table 7
  // The index Defijk = 5 is the CMS default
  // Comijk has the pairs i,j = l+j, had results to be included in combination
  // k==Est(10) is used in all cases
  static const Int_t Numijk = 7;
  static const Int_t Defijk = 5;
  const Int_t    Comijk[Numijk*2] = {4,8,  5,7,  4,7,  5,8,  6,8,  6,9,  4,9};
  const TString  Namijk[Numijk]   = {"211", "121", "221", "111", 
				     "h11", "hh1", "2h1"};
  // Arrays for combinations
  Double_t Valijk[Numijk] = {0};
  Double_t Staijk[Numijk] = {0};
  Double_t Sysijk[Numijk] = {0};
  Double_t Fulijk[Numijk] = {0};

  // Original indices for the order of Table 8
  //- The names
  const Int_t IndPap[Numijk] = {0, 3, 1, 2, 10, 6, 9};
  //- Their indices in the array of active estimates
  const Int_t IndRho[Numijk] = {0, 3, 1, 2,  6, 4, 5};
  //- Indices for the negative variances occuring for nominal correlations
  const Int_t NegVar[2] = {5, 11};

  // Arrays for combination per channel
  Double_t ValCha[MaxObs] = {0};
  Double_t StaCha[MaxObs] = {0};
  Double_t SysCha[MaxObs] = {0};
  Double_t FulCha[MaxObs] = {0};
  //-- End

  // Define formats for figures and latex file
  const TString FilBas = "B_PRD93_072004";
  TString FilNam;
  char Buffer[100];
  const TString ForVal = "%5.2f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%4.3f";
  const TString ForRho = "%4.2f";
  const TString ForPul = ForRho;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";

  // Construct Object
  Blue* myBlue;
  if(Flag == 2){myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  }else{myBlue = new Blue(NumEst, NumUnc);
  }
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);
  myBlue->SetQuiet();

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
    FilMat(RhoDif[k], RhoSam[k], RhoSou);
    myBlue->FillCor(k,RhoSou);
  }

  // Reduced or nominal correlations
  if(Flag == 1){
    printf("... B_PRD93_072004: -- Nominal -- correlations \n");
  }else{
    printf("... B_PRD93_072004: -- Reduced -- correlations \n");
    for(Int_t k = 1; k<NumUnc; k++){
      if(RhoDif[k] + RhoSam[k] > 0)myBlue->SetRhoRedUnc(k);
    }
  }

  // For Flag==0 store a dummy Latex file that has all inputs
  if(Flag == 0){
    myBlue->FixInp();
    myBlue->Solve();
    sprintf(Buffer,"AllInp");
    FilNam = &Buffer[0];       
    FilNam = FilBas + "_" + FilNam;
    myBlue->LatexResult(FilNam);
    myBlue->ReleaseInp();
  }
  
  // Solve according to Flag
  if(Flag < 2){
    // Run over scenarios
    for(Int_t l = 0; l<Numijk; l++){
      // De-activate all 2012 l+j and had results
      for(Int_t i = 4; i<NumEst-1; i++)myBlue->SetInActiveEst(i);
      // Activate the l+j - had pair to use
      for(Int_t i = l*2; i<(l+1)*2; i++)myBlue->SetActiveEst(Comijk[i]);
      myBlue->FixInp();
      myBlue->Solve();
      printf("... B_PRD93_072004: The %s combination", Namijk[l].Data());
      if(Flag == 0 && l == Defijk){
	printf(" == the final 2014 CMS combination, see Table 9");
      }
      printf("\n");
      myBlue->PrintResult();
      if(Flag == 0 && l == Defijk){
	printf("... B_PRD93_072004: The weights and pulls see Fig.13\n");
	myBlue->PrintWeight();
	myBlue->PrintPull();
      }

      // The default == hh1 combination
      if(l == Defijk){
	myBlue->PrintCompatEst();
	myBlue->PrintEst();
	printf("... B_PRD93_072004: \n");
	printf("... B_PRD93_072004: The correlations of the inputs");
	printf(" see Table 8");
	if(Flag == 1)printf(", but for nominal correlations");
	printf("\n");
	printf("... B_PRD93_072004: \n");
	myBlue->GetRho(LocRho);
	printf("... B_PRD93_072004:");
	for(Int_t i = 0; i<myBlue->GetActEst(); i++){
	  printf("%s  ", NamEst[IndPap[i]].Data());
	}	
	printf("\n");
	printf("... B_PRD93_072004: --------------------------------------");
	printf("-----------------------------\n");
	for(Int_t i = 0; i<myBlue->GetActEst(); i++){
	  printf("... B_PRD93_072004:");
	  for(Int_t j = 0; j<=i; j++){	    	 
	    printf(" %7.2f  ", LocRho->operator()(IndRho[i], IndRho[j]));
	  }
	  printf(" \n");
	}
	printf("... B_PRD93_072004: --------------------------------------");
	printf("-----------------------------\n");
    
	// Get output into file
	sprintf(Buffer,"%1i", Flag);
	FilNam = &Buffer[0];       
	FilNam = FilBas + "_" + FilNam;
	myBlue->LatexResult(FilNam);
	myBlue->DisplayResult(0,FilNam);

	// Show some interesting pairs
	myBlue->DisplayPair(6,9,FilNam);
	myBlue->DisplayPair(3,6,FilNam);
      }
      
      // Store result for later print out
      myBlue->GetResult(LocRes);
      myBlue->GetUncert(LocUnc);
      Valijk[l] = LocRes->operator()(0,0);
      Staijk[l] = LocRes->operator()(0,1);
      Sysijk[l] = SysUnc(LocUnc->operator()(0,0), LocRes->operator()(0,1));
      Fulijk[l] = LocUnc->operator()(0,0);
      myBlue->ReleaseInp();

      // For the default combination, also solve according to importance
      if(l == Defijk){
	sprintf(Buffer,"%1i", Flag);
	FilNam = &Buffer[0];       
	FilNam = FilBas + "_" + FilNam;
	myBlue->FixInp();
	myBlue->SolveAccImp(1.0);
	printf("... B_PRD93_072004: Result of the successive combination \n");
	myBlue->DisplayAccImp(0,FilNam);
	myBlue->ReleaseInp();
      }

      // For the default combination, also solve with positive weights only
      if(l == Defijk){
	sprintf(Buffer,"%1i_PosWei", Flag);
	FilNam = &Buffer[0];       
	FilNam = FilBas + "_" + FilNam;
	myBlue->FixInp();
	myBlue->SolvePosWei();
	printf("... B_PRD93_072004: Result with positive weights only \n");
	myBlue->PrintResult();
	myBlue->LatexResult(FilNam);
	myBlue->DisplayResult(0,FilNam);
	myBlue->ReleaseInp();
      }
    }
    
    // After all combinations are performed, write out the findings
    printf("... B_PRD93_072004: \n");
    printf("... B_PRD93_072004: The results of Table 7");
    if(Flag == 1)printf(", but for nominal correlations");
    printf("\n");
    printf("... B_PRD93_072004: \n");
    printf("... B_PRD93_072004: Combination:   Value  Stat  Syst  Full\n");
    printf("... B_PRD93_072004: --------------------------------------\n");
    for(Int_t l = 0; l<Numijk; l++){
      sprintf(Buffer,"%s=(%1i,%1i,%2i):", 
	      Namijk[l].Data(),Comijk[l*2], Comijk[l*2+1], 10);
      FilNam = &Buffer[0];
      printf("... B_PRD93_072004: %s %5.2f %5.2f %5.2f %5.2f \n",
	     FilNam.Data(), Valijk[l], Staijk[l], Sysijk[l], Fulijk[l]);
      if(l == 3 || l == 5 || l == 6)
	printf("... B_PRD93_072004: --------------------------------------\n");

    }
    if(Flag == 1){
      printf("... B_PRD93_072004: \n");
      printf("... B_PRD93_072004: When using -nominal- correlations");
      printf(" for some combinations the uncertainties are\n");
      printf("... B_PRD93_072004: much smaller then when using the concept");
      printf(" of -reduced- correlations.\n");
      printf("... B_PRD93_072004: This is caused by some pairs that now");
      printf(" are located to the right of the peak in the figure of\n");
      printf("... B_PRD93_072004: the combined uncertainty as a function of");
      printf(" the correlation, see *DisPai_Unc.pdf.\n");
      printf("... B_PRD93_072004: In addition, negative variances occur");
      printf(" for the sources %s and %s.\n",
	     NamUnc[NegVar[0]].Data(), NamUnc[NegVar[1]].Data());
      printf("... B_PRD93_072004: See the example B_NegVar.cxx");
      printf(" for an explanation of how this comes about.\n");
      printf("... B_PRD93_072004: \n");
      printf("... B_PRD93_072004: Although this feature is 'removed' by");
      printf(" construction when using --reduced-- correlations,\n");
      printf("... B_PRD93_072004: there is no scientific reason, why");
      printf(" this difference in sensitivity of the various estimators\n");
      printf("... B_PRD93_072004: to this source can be discussed away, by");
      printf(" an ad-hoc postulation of an additional systematic\n");
      printf("... B_PRD93_072004: source per pair of estimates, for which");
      printf(" the estimators are artifically assumed to\n");
      printf("... B_PRD93_072004: be uncorrelated. See EPJC74(2014)3004 for");
      printf(" a detailed discussion.\n");
      printf("... B_PRD93_072004: \n");     
      printf("... B_PRD93_072004: It rather seems advisable to investigate");
      printf(" the correlations in detail.\n");
      printf("... B_PRD93_072004: \n");
    }

  }else if(Flag == 2){

    // First the independent combinations per channel as in the paper
    // This is strictly speaking wrong since the estimates are correlated
    printf("... B_PRD93_072004: The independent combinations per channel,");
    printf(" see Figure 14\n");
    for(Int_t l = 0; l<NumObs; l++){
      for(Int_t i = 0; i<NumEst; i++)myBlue->SetInActiveEst(i);
      if(l ==0){
	//lepton+jets
	myBlue->SetActiveEst(1);
	myBlue->SetActiveEst(6);
	printf("... B_PRD93_072004: The independent lepton+jets result\n");
      }else if(l == 1){
	//all-jets
	myBlue->SetActiveEst(2);
	myBlue->SetActiveEst(9);
	printf("... B_PRD93_072004: The independent    hadronic result\n");
      }else if(l == 2){
	// dilepton
	myBlue->SetActiveEst( 0);
	myBlue->SetActiveEst( 3);
	myBlue->SetActiveEst(10);
	printf("... B_PRD93_072004: The independent    dilepton result\n");
      }
      myBlue->FixInp();
      myBlue->Solve();
      myBlue->PrintResult();
      myBlue->GetResult(LocRes);
      myBlue->GetUncert(LocUnc);
      ValCha[l] = LocRes->operator()(0,0);
      StaCha[l] = LocRes->operator()(0,1);
      SysCha[l] = SysUnc(LocUnc->operator()(0,0), LocRes->operator()(0,1));
      FulCha[l] = LocUnc->operator()(0,0);
      myBlue->ReleaseInp();
    }

    // The three observable combination for the hh1 combination
    myBlue->ResetInp();
    myBlue->SetInActiveEst(4);
    myBlue->SetInActiveEst(5);
    myBlue->SetInActiveEst(7);
    myBlue->SetInActiveEst(8);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_PRD93_072004: The combination for l+jets,");
    printf(" di-lepton and all-jets as correlated observables \n");
    myBlue->PrintResult();
    myBlue->PrintRhoRes();
    myBlue->PrintCompatObs();

    // Get the output into files
    sprintf(Buffer,"%1i", Flag);
    FilNam = &Buffer[0];       
    FilNam = FilBas + "_" + FilNam;
    myBlue->LatexResult(FilNam);
    for(Int_t n = 0; n<NumObs; n++)myBlue->DisplayResult(n,FilNam);

    // Store combined result and compare to individual ones
    myBlue->GetResult(LocRes);
    myBlue->GetUncert(LocUnc);

    printf("... B_PRD93_072004:");
    printf(" Observable Combination  Value  Stat  Syst Full \n");
    printf("... B_PRD93_072004:");
    printf("--------------------------------------------- \n");
    for(Int_t l = 0; l<NumObs; l++){
      printf("... B_PRD93_072004:");
      printf(" %s Individual: %5.2f %5.2f %5.2f %5.2f \n",
	     NamObs[l].Data(), ValCha[l], StaCha[l], SysCha[l], FulCha[l]);
      printf("... B_PRD93_072004:");
      printf("           Combined: %5.2f %5.2f %5.2f %5.2f \n",
	     LocRes->operator()(l,0), LocRes->operator()(l,1),
	     SysUnc(LocUnc->operator()(l,0), LocRes->operator()(l,1)),
	     LocUnc->operator()(l,0));
      printf("... B_PRD93_072004:");
      printf("--------------------------------------------- \n");
    }
  }

  // Delete object and matrices
  delete myBlue;  myBlue = NULL;

  RhoSou->Delete(); RhoSou = NULL;
  LocRho->Delete(); LocRho = NULL;
  LocRes->Delete(); LocRes = NULL;
  LocUnc->Delete(); LocUnc = NULL;
  
  // Return
  return;
};

//------------------------------------------------------------------------------
// Utility to fill a correlation matrix for any pair of RhoDif and RhoSam
//------------------------------------------------------------------------------

void FilMat(const Double_t RDif, const Double_t RSam, TMatrixD *const M){

  // Fill upper half of the correlation matrix
  M->operator()( 0, 1) = RDif;
  M->operator()( 0, 2) = RDif;
  M->operator()( 0, 3) = RSam;
  M->operator()( 0, 4) = RDif;
  M->operator()( 0, 5) = RDif;
  M->operator()( 0, 6) = RDif;
  M->operator()( 0, 7) = RDif;
  M->operator()( 0, 8) = RDif;
  M->operator()( 0, 9) = RDif;
  M->operator()( 0,10) = RSam;

  M->operator()( 1, 2) = RDif;
  M->operator()( 1, 3) = RDif;
  M->operator()( 1, 4) = RSam;
  M->operator()( 1, 5) = RSam;
  M->operator()( 1, 6) = RSam;
  M->operator()( 1, 7) = RDif;
  M->operator()( 1, 8) = RDif;
  M->operator()( 1, 9) = RDif;
  M->operator()( 1,10) = RDif;

  M->operator()( 2, 3) = RDif;
  M->operator()( 2, 4) = RDif;
  M->operator()( 2, 5) = RDif;
  M->operator()( 2, 6) = RDif;
  M->operator()( 2, 7) = RSam;
  M->operator()( 2, 8) = RSam;
  M->operator()( 2, 9) = RSam;
  M->operator()( 2,10) = RDif;

  M->operator()( 3, 4) = RDif;
  M->operator()( 3, 5) = RDif;
  M->operator()( 3, 6) = RDif;
  M->operator()( 3, 7) = RDif;
  M->operator()( 3, 8) = RDif;
  M->operator()( 3, 9) = RDif;
  M->operator()( 3,10) = RSam;

  M->operator()( 4, 5) = RSam;
  M->operator()( 4, 6) = RSam;
  M->operator()( 4, 7) = RDif;
  M->operator()( 4, 8) = RDif;
  M->operator()( 4, 9) = RDif;
  M->operator()( 4,10) = RDif;

  M->operator()( 5, 6) = RSam;
  M->operator()( 5, 7) = RDif;
  M->operator()( 5, 8) = RDif;
  M->operator()( 5, 9) = RDif;
  M->operator()( 5,10) = RDif;

  M->operator()( 6, 7) = RDif;
  M->operator()( 6, 8) = RDif;
  M->operator()( 6, 9) = RDif;
  M->operator()( 6,10) = RDif;

  M->operator()( 7, 8) = RSam;
  M->operator()( 7, 9) = RSam;
  M->operator()( 7,10) = RDif;

  M->operator()( 8, 9) = RSam;
  M->operator()( 8,10) = RDif;

  M->operator()( 9,10) = RDif;


  // Copy to lower half, add diagonal
  const Int_t NN = M->GetNrows();
  for(Int_t i = 0; i<NN; i++){
    M->operator()(i,i) = 1.;
    for(Int_t j = i; j<NN; j++){
      M->operator()(j,i) = M->operator()(i,j);
    }
  }
  return;
};

//------------------------------------------------------------------------------
// Calculate systematic from full and stat
//------------------------------------------------------------------------------
Double_t SysUnc(const Double_t ful, const Double_t sta){
  // Calculate sum in quadrature
  return TMath::Sqrt(ful*ful - sta*sta);
};
