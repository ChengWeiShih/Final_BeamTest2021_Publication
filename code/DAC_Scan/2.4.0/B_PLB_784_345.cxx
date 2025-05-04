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

void B_PLB_784_345(Int_t Flag = 0){

  //---------------------------------------------------------------------------
  // [Shows how to fill a TMatrix using FillEst()]
  // Flag == 0 The mass in 4L channel
  // Flag == 1 The mass in GG channel
  // Flag == 2 The mass for Run 1
  // Flag == 3 The mass for Run 2
  // Flag == 4 The mass from all 4 results
  // Flag == 5 The mass without GGR1
  //---------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  4;
  static const Int_t NumUnc =  2;
  static const Int_t NumObs =  1;

  // Set the names  
  const TString NamEst[NumEst] = {"M_{4L, R1}", "M_{GG, R1}", 
				  "M_{4L, R2}", "M_{GG, R2}"};
  const TString NamUnc[NumUnc] = {"   Stat", "  Syst"};
  TString NamObs[NumObs] = {"  M_{AT}"};

  // Steer according to Flag
  printf("... B_PLB_784_345: \n");
  if(Flag == 0){
    printf("... B_PLB_784_345: ---------------------------------\n");
    printf("... B_PLB_784_345: The mass in 4L channel, Flag = %2i \n", Flag);
    NamObs[0] = "  M_{4L}";
  }else if(Flag == 1){
    printf("... B_PLB_784_345: ---------------------------------\n");
    printf("... B_PLB_784_345: The mass in GG channel, Flag = %2i \n", Flag);
    NamObs[0] = "  M_{GG}";
  }else if(Flag == 2){
    printf("... B_PLB_784_345: -----------------------------\n");
    printf("... B_PLB_784_345: The mass for Run 1, Flag = %2i \n",  Flag);
    NamObs[0] = "  M_{R1}";
  }else if(Flag == 3){
    printf("... B_PLB_784_345: -----------------------------\n");
    printf("... B_PLB_784_345: The mass for Run 2, Flag = %2i \n", Flag);
    NamObs[0] = "  M_{R2}";
  }else if(Flag == 4){
    printf("... B_PLB_784_345: ----------------------------\n");
    printf("... B_PLB_784_345: The combined mass, Flag = %2i \n",Flag);
  }else if(Flag == 5){
    printf("... B_PLB_784_345: --------------------------------\n");
    printf("... B_PLB_784_345: The mass without GGR1, Flag = %2i \n",Flag);
    NamObs[0] = "M_{No GGR1}";
  }else{
    printf("... B_PLB_784_345: Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  // Estimates 0=4L(R1), 1=GG(R1), 2=4L(R2), 3=GG(R1)
  // From Fig4 simply guess that syst 4L(R1) = syst 4L(R2)
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  const Double_t XEst[LenXEst] = {
    // Observed uncertainties
    //         0     1
    //      Stat  Syst
    124.51, 0.52, DifUnc(0.37, 0.36),
    126.02, 0.43, DifUnc(0.51, 0.43),
    124.79, 0.36, DifUnc(0.37, 0.36),
    124.93, 0.21, DifUnc(0.40, 0.21)
  };

  // The combined values from the paper
  // 0=4L, 1=GG, 2=Run1, 3=Run2, 4,5=all
  static const Int_t NumPap = 6;
  static const Int_t LenYEst = NumPap * (NumUnc+1);
  const Double_t XPap[LenYEst] = {
    //         0     1
    //      Stat  Syst
    124.71, 0.30, DifUnc(0.30, 0.30),
    125.32, 0.19, DifUnc(0.35, 0.19),
    125.38, 0.37, DifUnc(0.41, 0.37),
    124.86, 0.18, DifUnc(0.27, 0.18),
    124.97, 0.16, DifUnc(0.24, 0.16),
    124.97, 0.16, DifUnc(0.24, 0.16)
  };

  //-- Local Structures for Blue input
  //- Estimates
  TMatrixD* InpEst = new TMatrixD(NumEst,NumUnc+1,&XEst[0]);
  //- Correlation per source filled dynamically
  TMatrixD* RhoSou = new TMatrixD(NumEst,NumEst);

  //- Results of the paper
  TMatrixD* PapRes = new TMatrixD(NumPap,NumUnc+1,&XPap[0]);
  //- Correlations from the paper combination of pairs
  TMatrixD* PapRho = new TMatrixD(NumEst,NumEst);
  //-- End

  //-- Local Structures for Blue output
  //- TMatrices
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRes = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc = new TMatrixD(NumObs,1);
  //-- End
 
  // Fill correlation matrix for the four channels from the information of the
  // pairwise combination in the paper
  Double_t ti = 0., tj = 0., zva = 0.;
  Int_t    ii = 0,  jj = 0,  ir = 1;
  PapRho->UnitMatrix();
  for(Int_t i = 0; i<NumEst; i++){
    for(Int_t j = i; j<NumEst; j++){
      ii = i;
      jj = j;
      // Set the combination we use
      ir = -1.;
      if(i == 0 && j == 1){ir = 2;
      }else if(i == 0 && j == 2){ir = 0;
      }else if(i == 1 && j == 3){ir = 1;
      }else if(i == 2 && j == 3){ir = 3;
      }
      if(ir >= 0){
	ti = TotUnc(InpEst->operator()(i,1), InpEst->operator()(i,2));
	tj = TotUnc(InpEst->operator()(j,1), InpEst->operator()(j,2));
	zva = tj/ti;
	if(ti > tj){
	  ii=j;
	  jj=i;
	  zva = 1. / zva;
	}
	PapRho->operator()(ii,jj) = RhoVal(BetVal(PapRes->operator()(ir,0),
						  InpEst->operator()(ii,0),
						  InpEst->operator()(jj,0)),
					   zva);
	PapRho->operator()(jj,ii) = PapRho->operator()(ii,jj);
      }
    }
  }

  // Define file names
  const TString FilBas = "B_PLB_784_345";
  TString FilNam = "To be filled later";
  char Buffer[100];
  sprintf(Buffer,"%1i", Flag);
  FilNam = &Buffer[0];
  FilNam = FilBas + "_" + FilNam;

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%+5.3f";
  const TString ForRho = ForUnc;
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
    if(k == 0){myBlue->FillCor(k, 0.);
    }else{myBlue->FillCor(k, PapRho);
    }
  }

  // Process according to Flag
  if(Flag == 0){
    // The mass in the 4L channel
    myBlue->SetInActiveEst(1);
    myBlue->SetInActiveEst(3);
  }else if(Flag == 1){
   // The mass in the GG channel
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(2);
  }else if(Flag == 2){
    // The mass in Run 1
    myBlue->SetInActiveEst(2);
    myBlue->SetInActiveEst(3);
  }else if(Flag == 3){
    // The mass in Run 2
    myBlue->SetInActiveEst(0);
    myBlue->SetInActiveEst(1);
  }else if(Flag == 4){
    // The combined mass
  }else if(Flag == 5){
    // The mass without GGR1
    myBlue->SetInActiveEst(1);
  }

  // Fix input
  myBlue->FixInp();

  // Look at input when all estimates are used i.e. for Flag == 4
  if(Flag == 4){
    myBlue->PrintMatrix(PapRho);
    myBlue->PrintEst();
    myBlue->PrintCompatEst();
  }

  // Solve and print result
  myBlue->Solve();
  printf("... B_PLB_784_345: Combination from calculated Rho matrix \n");
  myBlue->PrintResult();

  // Compare to paper result
  Double_t PapUnc =
    TMath::Sqrt(PapRes->operator()(Flag,1)*PapRes->operator()(Flag,1) +
		PapRes->operator()(Flag,2)*PapRes->operator()(Flag,2));
  printf("... B_PLB_784_345: Combination from ");
  if(Flag == 5){printf("Flag=4");
  }else{printf(" paper");
  }
  printf(": %s", myBlue->GetNamObs(0).Data());
  printf(" = %5.2f (+- %4.2f +- %4.2f) = +- %4.2f\n",
	 PapRes->operator()(Flag,0), PapRes->operator()(Flag,1),
	 PapRes->operator()(Flag,2), PapUnc);

  if(Flag != 5){
    // Loop over correlation for Syst k = 1
    Int_t m = 0;
    Int_t IndRho = 0;
    Double_t DifPap = 50.;
    Double_t DifAct =  0.;
    printf("... B_PLB_784_345: \n");
    printf("... B_PLB_784_345: Loop over correlation for k = 1 for %s \n",
	   myBlue->GetNamObs(0).Data());
    myBlue->SetQuiet();
    for(Double_t Rho = -1.0; Rho<1.02; Rho = Rho+0.1){
      myBlue->ReleaseInp();
      myBlue->SetRhoValUnc(1, Rho);
      myBlue->FixInp();
      myBlue->Solve();
      myBlue->GetResult(LocRes);
      myBlue->GetUncert(LocUnc);
      
      // Use sum of differences of full and statistical uncertainties
      DifAct =
	TMath::Abs(LocRes->operator()(0,0) - PapRes->operator()(Flag,0)) +
	TMath::Abs(LocUnc->operator()(0,0) - PapUnc);
      if(DifAct < DifPap){
	DifPap = DifAct;
	IndRho = m;
      }
      printf("... B_PLB_784_345: %2i: Rho = %+5.2f", m, Rho);
      printf(" M = %5.2f (+- %4.2f +- %4.2f) = +- %4.2f\n",
	     LocRes->operator()(0,0), LocRes->operator()(0,1),
	     LocRes->operator()(0,2), LocUnc->operator()(0,0));
      m = m + 1;
    }
    printf("... B_PLB_784_345: Best match to paper result for m = %2i \n",
	   IndRho);
  }
  
  // Delete Object, clean up matrices and return 
  delete myBlue; myBlue = NULL;

  InpEst->Delete(); InpEst = NULL;
  RhoSou->Delete(); RhoSou = NULL;
  PapRes->Delete(); PapRes = NULL;
  PapRho->Delete(); PapRho = NULL;

  LocRho->Delete(); LocRho = NULL;
  LocRes->Delete(); LocRes = NULL;
  LocUnc->Delete(); LocUnc = NULL;

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
  Double_t Bet = (xva - xv1) / (xv2 - xv1);

  // Check whether within [-infty, 0.5] set to 0.5 if needed
  if(Bet > 0.5){
    printf("... BetVal xva = %+6.3f, xv1 = %6.3f, xv2 = %6.3f",
	   xva, xv1, xv2);
    printf(" bet = %6.3f -->  bet = 0.5 \n", Bet);
    Bet = 0.5;
  }

  // Return
  return Bet;
};

//------------------------------------------------------------------------
// Calculate rho from beta and z
//------------------------------------------------------------------------
Double_t RhoVal(const Double_t bet, const Double_t zva){

  // Calculate rho from manual Eq.2 first vs. third term
  // Check whether bet = 0.5 (if yes check wether zva = 1.0
  Double_t Rho;
  if(TMath::Abs(bet - 0.5) <= 0.0001){
    Rho = -1.;
    if(TMath::Abs(zva - 1.0) <= 0.0001){
      printf("... RhoVal bet = %+6.3f, zva = %6.3f -> rho = undefined \n", 
	     bet, zva);
    }
  }else{
    Rho = (1. - bet*(1 + zva*zva)) / (zva*(1. - 2.*bet));
  }

  // See what we got
  //printf("... RhoVal bet =%+6.3f, zva =%6.3f, rho =%6.3f \n", bet,zva, Rho);

  // Check whether within [-1.0, 1.0] set to boundary if needed
  if(Rho < -1.0 || Rho > 1.0){
    printf("... RhoVal bet = %+6.3f, zva = %6.3f, rho = %6.3f \n", 
	   bet,zva, Rho);    
    Rho = TMath::Max(-1.0, Rho);
    Rho = TMath::Min( Rho, 1.0);
  }

  return Rho;
};
