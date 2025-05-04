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

void B_EPJC_74_3004(Int_t Flag = 0){
  
  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  [Shows how to use SolveMaxVar()]
  //  0: Table 2 Scenario A, Figure 4
  //  1: Table 2 Scenario B 
  //  2: Table 2 Scenario C
  //  3: Table 2 Scenario D
  //----------------------------------------------------------------------------
  
  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 2;
  static const Int_t NumUnc = 3;
  static const Int_t NumObs = 1;

  // The scenarios
  static const Int_t NumSce = 4;
  TString SceNam[NumSce] = {"A", "B", "C", "D"};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_EPJC_74_3004: ----------------------------- \n");
    printf("... B_EPJC_74_3004: Table 1 Scenario A Flag = %2i \n",Flag);
  }else if(Flag == 1){
    printf("... B_EPJC_74_3004: ----------------------------- \n");
    printf("... B_EPJC_74_3004: Table 1 Scenario B Flag = %2i \n",Flag);
  }else if(Flag == 2){
    printf("... B_EPJC_74_3004: ----------------------------- \n");
    printf("... B_EPJC_74_3004: Table 1 Scenario C Flag = %2i \n",Flag);
  }else if(Flag == 3){
    printf("... B_EPJC_74_3004: ----------------------------- \n");
    printf("... B_EPJC_74_3004: Table 1 Scenario D Flag = %2i \n",Flag);
  }else{
    printf("... B_EPJC_74_3004: -------------------------- \n");
    printf("... B_EPJC_74_3004: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // names of estimates uncertaintries and observables
  const TString NamEst[NumEst] = {"     x1", "     x2"};
  const TString NamUnc[NumUnc] = {"   Stat", "   Sys1", "   Sys2"};
  const TString NamObs[NumObs] = {"      x"};

  // Define default estimates and correlations
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    172.10, 0.60, 0.50, 0.70,
    173.10, 0.40, 0.70, 1.40
  };
  Double_t RhoVal[NumUnc] = {0.0, 1.0, 1.0};

  // Change assumptions according to Flag
  if(Flag == 1){RhoVal[1] = 0.0;
  }else if(Flag == 2){XEst[7] = 0.7;
  }else if(Flag == 3){XEst[3] = 1.4;
  }

  // Define formats for figures and latex file
  const TString FilBas = "B_EPJC_74_3004";
  TString FilNam;
  char Buffer[100];
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%4.3f";
  const TString ForRho = "%4.2f";
  const TString ForPul = ForRho;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";

  // Output file names for InspectPair
  sprintf(Buffer,"_%s", SceNam[Flag].Data());
  FilNam = &Buffer[0];
  FilNam = FilBas + FilNam;
  const TString InsDef = FilNam + "_Defa";
  const TString InsRed = FilNam + "_Redc";
  const TString InsMax = FilNam + "_Maxv";
  const TString TexDef = FilNam + "_Defa_t";
  const TString TexRed = FilNam + "_Redc_t";
  const TString TexMax = FilNam + "_Maxv_t";
  const TString DisDef = FilNam;

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

  // Fill names 
  myBlue->FillNamEst(&NamEst[0]); 
  myBlue->FillNamUnc(&NamUnc[0]); 
  myBlue->FillNamObs(&NamObs[0]); 

  // Fill all estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }

  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++)myBlue->FillCor(k,RhoVal[k]);

  // Solve 
  myBlue->FixInp();
  myBlue->Solve();
  printf("... B_EPJC_74_3004: BLUE result for Flag = %2i \n", Flag);
  if(Flag == 0){myBlue->InspectPair(0,1,InsDef,1);
  }else{myBlue->InspectPair(0,1,InsDef);
  }
  myBlue->PrintEst();
  myBlue->PrintRho();
  myBlue->PrintCompatEst();
  myBlue->PrintParams();
  myBlue->PrintResult();
  myBlue->DisplayResult(0,DisDef);
  myBlue->LatexResult(TexDef);

  // Solve with reduced correlations
  myBlue->ResetInp();
  myBlue->SetRhoRedUnc(1);
  myBlue->SetRhoRedUnc(2);
  myBlue->FixInp();
  printf("... B_EPJC_74_3004: Reduced correlations for Flag = %2i \n", Flag);
  myBlue->InspectPair(0,1,InsRed);
  myBlue->Solve();
  myBlue->PrintCompatEst();
  myBlue->PrintRho();
  myBlue->PrintResult();
  myBlue->LatexResult(TexRed);

  // Solve by maximising the variance
  Int_t IntMax[3] ={0, -1, 2};
  for(Int_t l=0; l<3 ;l++){
    myBlue->ResetInp();
    myBlue->FixInp();
    myBlue->SolveMaxVar(IntMax[l]);
    myBlue->PrintCompatEst();
    myBlue->PrintMaxVar();
    printf("... B_EPJC_74_3004: MaxVar(%2i) results for Flag = %2i \n",
	   IntMax[l], Flag);
    if(Flag == 0 && l == 0){
      myBlue->InspectPair(0,1,InsMax);
      myBlue->LatexResult(TexMax);
    }
    myBlue->PrintResult();
  }
  
  // Delete object
  delete myBlue; myBlue = NULL;
  return;
}
