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

void B_arXiv_1307_4003(Int_t Flag = 0){
  
  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The paper combination
  //  1: The combination of   two estimates from Table 1
  //  2: The combination of three estimates from Table 1
  //     [Shows how to use SolvePosWei()]
  //  3: The calculation of the covariance
  //----------------------------------------------------------------------------
  
  // The number of estimates, uncertainties and observables
  static const Int_t NF0Est = 2;
  static const Int_t NF0Unc = 1;
  static const Int_t NF1Est = 3;
  static const Int_t NF1Unc = 1;
  static const Int_t NF2Est = 4;
  static const Int_t NF2Unc = 3;

  Int_t NumEst = 0, NumUnc = 0;

  TString NamF0Est[NF0Est] = {"      A", "      B"};
  TString NamF0Unc[NF0Unc] = {"   Full"};
  TString NamF1Est[NF1Est] = {"      A", "     B1", "     B2"};
  TString NamF1Unc[NF1Unc] = {"   Full"};
  TString NamF2Est[NF2Est] = {"     YA", "     YB", "     YC", "     YD"};
  TString NamF2Unc[NF2Unc] = {"    Unc", "   BKGD", "   Lumi"};

  static const Int_t MaxObs =  1;
  TString NamObs[MaxObs] = {" Combin"};

  // Preset according to Flag
  if(Flag == 0){
    NumEst = NF0Est;
    NumUnc = NF0Unc;
    printf("... B_arXiv_1307_4003: ------------------------------------------");
    printf("-----------\n");
    printf("... B_arXiv_1307_4003: The two  estimate combination from");
    printf(" Table 1, Flag = %2i \n", Flag);
  }else if(Flag == 1){
    NumEst = NF1Est;
    NumUnc = NF1Unc;
    printf("... B_arXiv_1307_4003: ------------------------------------------");
    printf("------------\n");
    printf("... B_arXiv_1307_4003: The three estimate combination from");
    printf(" Table 1, Flag = %2i \n", Flag);
  }else if(Flag == 2){
    NumEst = NF2Est;
    NumUnc = NF2Unc;
    printf("... B_arXiv_1307_4003: ----------------------------------------\n");
    printf("... B_arXiv_1307_4003: The combination from Table 4,");
    printf(" Flag = %2i \n", Flag);
  }else if(Flag == 3){
    NumEst = NF2Est;
    NumUnc = NF2Unc;
    printf("... B_arXiv_1307_4003: ------------------------------------------");
    printf("----------\n");
    printf("... B_arXiv_1307_4003: The covariance for f_ijk = f_ij for all k,");
    printf(" Flag = %2i \n", Flag);
  }else{
    printf("... B_arXiv_1307_4003: -------------------------- \n");
    printf("... B_arXiv_1307_4003: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // The input from Table 1
  static const Int_t LenXF0Est = NF0Est * (NF0Unc+1);
  Double_t XF0Est[LenXF0Est] = {
    103.00, TMath::Sqrt(15),
     98.00, TMath::Sqrt(10)
  };
  static const Int_t LenXF1Est = NF1Est * (NF1Unc+1);
  Double_t XF1Est[LenXF1Est] = {
    103.00, TMath::Sqrt(15),
     99.00, 4.00,
    101.00, 8.00
  };
  static const Int_t LenF1Cor = NF1Est * NF1Est;
  Double_t CorF1[LenF1Cor] = {
    1.0,   0.0,   0.0,
    0.0,   1.0, 0.875,
    0.0, 0.875,   1.0
  };

  // The input from Table 4
  static const Int_t LenXF2Est = NF2Est * (NF2Unc+1);
  Double_t XF2Est[LenXF2Est] = {
     95.00, 10.00, 10.00, 11.00,
    144.00, 14.00, 40.00, 14.00,
    115.00, 18.00,  3.00, 10.00,
    122.00, 25.00,  0.00,  0.00
  };

  static const Int_t LenF2Cor = NF2Est * NF2Est;
  Double_t CorF2[LenF2Cor] = {
          1.0, 272./554., 140./140., 0.0,
    272./554.,       1.0, 219./260., 0.0,
    140./140., 219./260.,       1.0, 0.0, 
          0.0,       0.0,       0.0, 1.0
  };

  // Local structure for Blue output and systematic uncertainties Table 4
  Int_t iok = 0;
  static const Int_t NumObs = 1;
  Double_t UncDef[NumObs] = {0};
  static const Int_t LenRes = NumObs*(NF2Unc+1);
  Double_t ValDef[LenRes] = {0};
  Double_t Weight[NF2Est] = {0};
  Double_t Chiq = 0;
  Int_t    Ndof = 0;
  static const Int_t NumSys = 7;
  Double_t Values[NumSys] = {0};
  Double_t Uncert[NumSys] = {0};
  Double_t UncUnc[NumSys] = {0};
  Double_t UncBKD[NumSys] = {0};
  Double_t UncLUM[NumSys] = {0};
  Double_t UncChi[NumSys] = {0};
  Int_t    Numdof[NumSys] = {0};
  Double_t UncWgA[NumSys] = {0};
  Double_t UncWgB[NumSys] = {0};
  Double_t UncWgC[NumSys] = {0};
  Double_t UncWgD[NumSys] = {0};
  // End

  // Define formats for figures and latex file
  const TString FilBas = "B_arXiv_1307_4003";
  const TString ForVal = "%5.3f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%5.3f";
  const TString ForRho = "%5.3f";
  const TString ForPul = ForRho;
  const TString ForUni = "None";

  // Construct object
  Blue *myBlue = new Blue(NumEst, NumUnc);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);
 
  // Fill the stuff
  Int_t ind = 0;
  myBlue->FillNamObs(&NamObs[0]);
  if(Flag == 0){
    // Fill names
    myBlue->FillNamEst(&NamF0Est[0]);
    myBlue->FillNamUnc(&NamF0Unc[0]);
    // Fill estimates and correlations
    ind = 0;
    for(Int_t i = 0; i<NumEst; i++){
      myBlue->FillEst(i,&XF0Est[ind]);
      ind = ind + NumUnc + 1;
    }
    // Fill correlations
    myBlue->FillCor(0,0.0);
  }else if(Flag == 1){
    // Fill names
    myBlue->FillNamEst(&NamF1Est[0]);
    myBlue->FillNamUnc(&NamF1Unc[0]);
    // Fill estimates and correlations
    ind = 0;
    for(Int_t i = 0; i<NumEst; i++){
      myBlue->FillEst(i,&XF1Est[ind]);
      ind = ind + NumUnc + 1;
    }
    // Fill correlations
    myBlue->FillCor(0,&CorF1[0]);
  }else if(Flag >= 2){
    printf("... B_arXiv_1307_4003: The combination from Table 4 = %2i \n",Flag);
    // Fill names
    myBlue->FillNamEst(&NamF2Est[0]);
    myBlue->FillNamUnc(&NamF2Unc[0]);
    // Fill estimates
    ind = 0;
    for(Int_t i = 0; i<NumEst; i++){
      myBlue->FillEst(i,&XF2Est[ind]);
      ind = ind + NumUnc + 1;
    }
    // Fill correlations
    for(Int_t k = 0; k<NumUnc; k++){
      if(k == 0){myBlue->FillCor(k,0.0);
      }else{
	if(Flag == 2){myBlue->FillCor(k,1.0);
	}else{myBlue->FillCor(k,&CorF2[0]);
	}
      }
    }
  }
  
  // Do the calculations
  if(Flag < 2){
    printf("... B_arXiv_1307_4003: Now Produce the central result \n");
    myBlue->FixInp();
    myBlue->SolveInfWei();
    if(Flag == 0){
      printf("... B_arXiv_1307_4003: The combination of A, B from Table 1,");
      printf(" Flag = %2i \n", Flag);
    }else if(Flag == 1){
      printf("... B_arXiv_1307_4003: The combination of A, B1, B2");
      printf(" from Table 1, Flag = %2i \n", Flag);
    }
    // Produce the central value, print the result
    myBlue->PrintInfWei();
    myBlue->PrintEst();
    myBlue->PrintResult();
  }else if(Flag == 2){
    printf("... B_arXiv_1307_4003: The combination from Table 4, %2i \n", Flag);
    // 1) The nominal correlation
    ind = -1;
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1307_4003: The nominal correlations\n");
    printf("... B_arXiv_1307_4003: For the covariance matrizes of the");
    printf(" uncertainties see Table 5 \n");
    printf("... B_arXiv_1307_4003: BKGD \n");
    myBlue->PrintCov(1);
    printf("... B_arXiv_1307_4003: LUMI \n");
    myBlue->PrintCov(2);
    printf("... B_arXiv_1307_4003: For the total covariance matrix");
    printf(" see Table 6 \n");
    myBlue->PrintCov();

    // 1) Get results store locally
    ind = ind + 1;
    iok = myBlue->GetResult(ValDef); 
    iok = myBlue->GetUncert(UncDef); 
    iok = myBlue->GetWeight(Weight);
    if(iok == 0)printf("... B_arXiv_1307_4003: Error retrieving result \n");
    Chiq = myBlue->GetChiq(); Ndof = myBlue->GetNdof();
    Values[ind] = ValDef[0]; Uncert[ind] = UncDef[0];
    UncUnc[ind] = ValDef[1]; UncBKD[ind] = ValDef[2]; 
    UncLUM[ind] = ValDef[3];
    UncChi[ind] = Chiq; Numdof[ind] = Ndof;
    UncWgA[ind] = Weight[0]; UncWgB[ind] = Weight[1]; 
    UncWgC[ind] = Weight[2]; UncWgD[ind] = Weight[3];

    // 2) Fijk = F    for all ijk
    // 3) Fijk = Fk   for all ij
    // 4) Fijk = Fij  for all k
    for(Int_t l = 0; l<3; l++){
      myBlue->ResetInp();
      myBlue->FixInp();
      myBlue->SolveMaxVar(l);
      if(l == 0){printf("... B_arXiv_1307_4003: Fijk = F   for all ijk \n");
      }else if(l == 1){printf("... B_arXiv_1307_4003: Fijk = Fk  for all ij\n");
      }else if(l == 2){printf("... B_arXiv_1307_4003: Fijk = Fij for all k\n");
      }
      myBlue->PrintMaxVar();
      printf("... B_arXiv_1307_4003: For the total covariance matrix");
      printf(" see Table 6 \n");
      myBlue->PrintCov();
      // Get results store locally
      ind = ind + 1;
      iok = myBlue->GetResult(ValDef); iok = myBlue->GetUncert(UncDef); 
      iok = myBlue->GetWeight(Weight);
      if(iok == 0)printf("... B_arXiv_1307_4003: Error retrieving result \n");
      Chiq = myBlue->GetChiq(); Ndof = myBlue->GetNdof();
      Values[ind] = ValDef[0]; Uncert[ind] = UncDef[0];
      UncUnc[ind] = ValDef[1]; UncBKD[ind] = ValDef[2]; UncLUM[ind] = ValDef[3];
      UncChi[ind] = Chiq; Numdof[ind] = Ndof;
      UncWgA[ind] = Weight[0]; UncWgB[ind] = Weight[1]; 
      UncWgC[ind] = Weight[2]; UncWgD[ind] = Weight[3];
      if(l == 2){
	printf("... B_arXiv_1307_4003: Due to a different minimum found");
	printf(" the covariance matrix does not coincide \n");
	printf("... B_arXiv_1307_4003: with the preprint. However, the result");
	printf(" is identical, see Table 4. \n");
	printf("... B_arXiv_1307_4003: The difference in the variance of");
	printf(" the result %8.5f is minimal, see Flag 3 \n", 
	       Uncert[ind]*Uncert[ind]);
      }
    }

    // 5) Do not use estimates with negative weights
    myBlue->ResetInp();
    myBlue->FixInp();
    myBlue->PrintListEst();
    myBlue->SolvePosWei();

    printf("... B_arXiv_1307_4003: No estimates with negative weights \n");
    myBlue->PrintCov();
    // Get results store locally
    ind = ind + 1;
    iok = myBlue->GetResult(ValDef); 
    iok = myBlue->GetUncert(UncDef); 
    iok = myBlue->GetWeight(Weight);
    if(iok == 0)printf("... B_arXiv_1307_4003: Error retrieving result \n");
    Chiq = myBlue->GetChiq(); Ndof = myBlue->GetNdof();
    Values[ind] = ValDef[0]; Uncert[ind] = UncDef[0];
    UncUnc[ind] = ValDef[1]; UncBKD[ind] = ValDef[2]; UncLUM[ind] = ValDef[3];
    UncChi[ind] = Chiq; Numdof[ind] = Ndof;

    // Get the indices of the active estimates
    Int_t IndEst[NF2Est] = {0};
    iok = myBlue->GetIndEst(IndEst);
    printf("... B_arXiv_1307_4003: Active     estimates:");
    for(Int_t i = 0; i<myBlue->GetActEst(); i++)printf("%2i",IndEst[i]);
    printf("\n");

    // First reset, then put weights at the right place
    UncWgA[ind] = 0, UncWgB[ind] = 0, UncWgC[ind] = 0, UncWgD[ind] = 0;
    Int_t j = -1;
    for(Int_t i = 0; i<myBlue->GetActEst(); i++){
      j = IndEst[i];
      if(j == 0){UncWgA[ind] = Weight[i];
      }else if(j == 1){UncWgB[ind] = Weight[i];
      }else if(j == 2){UncWgC[ind] = Weight[i];
      }else if(j == 3){UncWgD[ind] = Weight[i];
      }
    }

    // 6) Reduced correlations
    myBlue->ResetInp();
    myBlue->SetRhoRedUnc();
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1307_4003: Reduced correlations \n");
    printf("... B_arXiv_1307_4003: For the covariance matrizes of the");
    printf(" uncertainties see Table 5 \n");
    printf("... B_arXiv_1307_4003: BKGD \n");
    myBlue->PrintCov(1);
    printf("... B_arXiv_1307_4003: LUMI \n");
    myBlue->PrintCov(2);
    printf("... B_arXiv_1307_4003: For the total covariance matrix");
    printf(" see Table 6 \n");
    myBlue->PrintCov();
    // Get results store locally
    ind = ind + 1;
    iok = myBlue->GetResult(ValDef); 
    iok = myBlue->GetUncert(UncDef); 
    iok = myBlue->GetWeight(Weight);
    if(iok == 0)printf("... B_arXiv_1307_4003: Error retrieving result \n");
    Chiq = myBlue->GetChiq(); Ndof = myBlue->GetNdof();
    Values[ind] = ValDef[0]; Uncert[ind] = UncDef[0];
    UncUnc[ind] = ValDef[1]; UncBKD[ind] = ValDef[2]; 
    UncLUM[ind] = ValDef[3];
    UncChi[ind] = Chiq; Numdof[ind] = Ndof;
    UncWgA[ind] = Weight[0]; UncWgB[ind] = Weight[1]; 
    UncWgC[ind] = Weight[2]; UncWgD[ind] = Weight[3];

    // 7) No correlations
    myBlue->ResetInp();
    myBlue->SetRhoValUnc(1,0);
    myBlue->SetRhoValUnc(2,0);
    myBlue->FixInp();
    myBlue->Solve();
    printf("... B_arXiv_1307_4003: No correlations \n");
    // Get results store locally
    ind = ind + 1;
    iok = myBlue->GetResult(ValDef); 
    iok = myBlue->GetUncert(UncDef); 
    iok = myBlue->GetWeight(Weight);
    if(iok == 0)printf("... B_arXiv_1307_4003: Error retrieving result \n");
    Chiq = myBlue->GetChiq(); Ndof = myBlue->GetNdof();
    Values[ind] = ValDef[0]; Uncert[ind] = UncDef[0];
    UncUnc[ind] = ValDef[1]; UncBKD[ind] = ValDef[2]; UncLUM[ind] = ValDef[3];
    UncChi[ind] = Chiq; Numdof[ind] = Ndof;
    UncWgA[ind] = Weight[0]; UncWgB[ind] = Weight[1]; 
    UncWgC[ind] = Weight[2]; UncWgD[ind] = Weight[3];

    // Now report the Table
    printf("... B_arXiv_1307_4003: The final result for Table 4 are : \n");
    printf("... B_arXiv_1307_4003: ------------------------------------------");
    printf("------------------------------------");
    printf("---------------------------------\n");
    for(Int_t i = 0; i<NumSys; i++){
      printf("... B_arXiv_1307_4003: %1i",i);
      if(i ==0){printf(       "                 Nominal: ");
      }else if(i == 1){printf("  Fijk = F   for all ijk: ");
      }else if(i == 2){printf("   Fijk = Fk  for all ij: ");
      }else if(i == 3){printf("    Fijk = Fij for all k: ");
      }else if(i == 4){printf(" No estimate with Wgt< 0: ");
      }else if(i == 5){printf("     Reduced Correlatons: ");
      }else if(i == 6){printf("            Uncorrelated: ");
      }
      printf("%5.2f +- %5.2f | +- %5.2f +-  %4.2f +-  %4.2f | %3.1f / %1i", 
	     Values[i], Uncert[i], UncUnc[i], UncBKD[i], UncLUM[i], 
	     UncChi[i], Numdof[i]);
      printf(" | %+5.1f%% %+5.1f%% %+5.1f%% %+5.1f%% \n", 
	     100*UncWgA[i], 100*UncWgB[i], 100*UncWgC[i], 100*UncWgD[i]);
    }
    printf("... B_arXiv_1307_4003: ------------------------------------------");
    printf("------------------------------------");
    printf("---------------------------------\n");
  }else if(Flag == 3){
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->PrintResult();
    myBlue->PrintWeight();

    ind = 3;
    iok = myBlue->GetResult(ValDef);
    iok = myBlue->GetUncert(UncDef); 
    iok = myBlue->GetWeight(Weight);
    if(iok == 0)printf("... B_arXiv_1307_4003: Error retrieving result \n");
    Chiq = myBlue->GetChiq(); Ndof = myBlue->GetNdof();
    Values[ind] = ValDef[0]; Uncert[ind] = UncDef[0];
    UncUnc[ind] = ValDef[1]; UncBKD[ind] = ValDef[2]; UncLUM[ind] = ValDef[3];
    UncChi[ind] = Chiq; Numdof[ind] = Ndof;
    UncWgA[ind] = Weight[0]; UncWgB[ind] = Weight[1]; 
    UncWgC[ind] = Weight[2]; UncWgD[ind] = Weight[3];
    printf("... B_arXiv_1307_4003: The covariance matrix from the preprint");
    printf(" for Fijk + Fij \n");
    myBlue->PrintCov();
    printf("... B_arXiv_1307_4003: The difference in the variance of the");
    printf(" result %8.5f is minimal, see Flag 2 \n", Uncert[ind]*Uncert[ind]);

    // Now report the Table
    printf("... B_arXiv_1307_4003: The final result for Table 4 are : \n");
    printf("... B_arXiv_1307_4003: ----------------------------------------");
    printf("---------------------------------------------------------------");
    printf("-----\n");
    for(Int_t i = 3; i<4; i++){
      if(i == 0){printf("... B_arXiv_1307_4003:                 Nominal");
      }else if(i == 1){printf("... B_arXiv_1307_4003:  Fijk = F   for all ijk");
      }else if(i == 2){printf("... B_arXiv_1307_4003:   Fijk = Fk  for all ij");
      }else if(i == 3){printf("... B_arXiv_1307_4003:    Fijk = Fij for all k");
      }else if(i == 4){printf("... B_arXiv_1307_4003: No estimate with Wgt< 0");
      }else if(i == 5){printf("... B_arXiv_1307_4003:     Reduced Correlatons");
      }else if(i == 6){printf("... B_arXiv_1307_4003:            Uncorrelated");
      }
      printf(": %5.2f +- %5.2f | +- %5.2f +-  %4.2f +-  %4.2f | %3.1f / %1i",
	     Values[i], Uncert[i], 
	     UncUnc[i], UncBKD[i], UncLUM[i], 
	     UncChi[i], Numdof[i]);
      printf(" | %+5.1f%% %+5.1f%% %+5.1f%% %+5.1f%% \n", 
	     100*UncWgA[i], 100*UncWgB[i], 100*UncWgC[i], 100*UncWgD[i]);
    }
    printf("... B_arXiv_1307_4003: ----------------------------------------");
    printf("---------------------------------------------------------------");
    printf("-----\n");
  }
  // Delete object
  delete myBlue; myBlue = NULL;
  return;
}
