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

void B_Peelles(Int_t Flag = 0){
  //---------------------------------------------------------------------------
  // Flag steers which of the results should be calculated.
  // For each scenario the results with both absolute and relative 
  // uncertainties are calculated.
  // 0: The original version of the Puzzle
  // 1: = 0 but with twice the uncertainties
  // 2: = 0 but with a more consistent estimate x_2 but unchanged uncertainties
  // 3: = 0 but with a more consistent estimate x_2 and    scaled uncertainties
  // 4: = 0 but with more consistent input due to an adapted correlation
  // 5: Another numerical version of the Puzzle
  //---------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 2;
  static const Int_t NumUnc = 2;
  static const Int_t NumObs = 1;

  // The names
  const TString NamEst[NumEst] = {" Est_{1}"," Est_{2}"};
  const TString NamUnc[NumEst] = {" Unc_{1}"," Unc_{2}"};
  const TString NamObs[NumEst] = {" Obs_{1}"};

  // The scenario name
  TString SceNam = "to be filled later";
  TString ModNam = "to be filled later";
  TString FilNam = "to be filled later";
  TString PubNam = "Peelles";

  // Steer the input depending on Flag
  if(Flag == 0){
    printf("... B_Peelles: ------------------------------ \n");
    printf("... B_Peelles: The original Puzzle Flag = %2i \n", Flag);
    SceNam = "A";
  }else if(Flag == 1){
    printf("... B_Peelles: --------------------------");
    printf("------------------------ \n");
    printf("... B_Peelles: The Puzzle with twice");
    printf(" the uncertainties Flag = %2i \n", Flag);
    SceNam = "B";
  }else if(Flag == 2){
    printf("... B_Peelles: --------------------------");
    printf("----------------------------------------------------------- \n");
    printf("... B_Peelles: The Puzzle with a more consistent estimate x_2 and");
    printf(" unchanged uncertainties Flag = %2i \n", Flag);
    SceNam = "C";
  }else if(Flag == 3){
    printf("... B_Peelles: --------------------------");
    printf("----------------------------------------------------------- \n");
    printf("... B_Peelles: The Puzzle with a more consistent estimate x_2 and");
    printf("    scaled uncertainties Flag = %2i \n", Flag);
    SceNam = "D";
  }else if(Flag == 4){
    printf("... B_Peelles: --------------------------");
    printf("---------------------------------------------------- \n");
    printf("... B_Peelles: The Puzzle with more consistent input due to an\n");
    printf(" adapted correlation Flag = %2i \n", Flag);
    SceNam = "E";
  }else if(Flag == 5){
    printf("... B_Peelles: --------------------------");
    printf("------------------------------------------ \n");
    printf("... B_Peelles: Another version of the Puzzle with absolute \n");
    printf(" uncertainties Flag = %2i \n", Flag);
    SceNam = "AltVers";
  }else{
    printf("... B_Peelles: -------------------------- \n");
    printf("... B_Peelles: Not implemented Flag = %2i \n", Flag);
    return;
  }

  //  Store all possible estimates
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t X0Est[LenXEst] = {1.00, 0.10, 0.20, 1.50, 0.150, 0.30};
  Double_t X1Est[LenXEst] = {1.00, 0.20, 0.40, 1.50, 0.300, 0.60};
  Double_t X2Est[LenXEst] = {1.00, 0.10, 0.20, 1.25, 0.150, 0.30};
  Double_t X3Est[LenXEst] = {1.00, 0.10, 0.20, 1.25, 0.125, 0.25};
  Double_t X4Est[LenXEst] = {1.00, 0.10, 0.20, 1.50, 0.150, 0.30};
  Double_t X5Est[LenXEst] = {8.00, 0.16, 0.80, 8.50, 0.170, 0.85};

  // Correlation for systematic uncertainties
  Double_t RhoVal = 1.0;
  if(Flag == 4)RhoVal = 0.05;

  // Array for relative uncertainties
  Double_t ActCof[3] = {0., 0., 0.};

  //-- Local Structures for BLUE output
  // TMatrices
  TMatrixD* LocEst    = new TMatrixD(NumEst,NumUnc+1);
  TMatrixD* LocEstVal = new TMatrixD(NumEst,1);
  TMatrixD* LocEstUnc = new TMatrixD(NumEst,1);
  //-- End

  // Define formats for figures and latex file
  const TString FilBas = "B_Peelles";
  const TString ForVal = "%5.2f";
  const TString ForUnc = "%4.2f";
  const TString ForWei = "%4.3f";
  const TString ForRho = "%5.3f";
  const TString ForPul = ForRho;
  const TString ForUni = "mb";

  //
  // ----------------------------- End of Input ------------------------------
  //

  // Select input according to Flag
  Double_t XEst[LenXEst] = {0};
  for(Int_t l = 0; l<LenXEst; l++){
    if(Flag == 0){XEst[l] = X0Est[l];
    }else if(Flag == 1){XEst[l] = X1Est[l];
    }else if(Flag == 2){XEst[l] = X2Est[l];
    }else if(Flag == 3){XEst[l] = X3Est[l];
    }else if(Flag == 4){XEst[l] = X4Est[l];
    }else if(Flag == 5){XEst[l] = X5Est[l];
    }
  }

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);

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

  // Fill correlation
  myBlue->FillCor(0,0.0);
  myBlue->FillCor(1,RhoVal);
 
  // The header
  if(Flag < 5){
    printf("... B_Peelles: ---------- Scenario %s: Table 1 of %s ----------\n",
	   SceNam.Data(), PubNam.Data());
  }else{
    printf("... B_Peelles: ---------- Scenario %s: NOT discussed in %s ----\n",
	   SceNam.Data(), PubNam.Data());
  }

  // Quiet mode
  myBlue->SetQuiet();

  // Loop over absolute or relative uncertainties
  for(Int_t l = 0; l < 2; l++){    
    if(l == 0){
      //-- The absolute uncertainties
      ModNam = "absolute";
      myBlue->FixInp();
      myBlue->Solve();
      printf("... B_Peelles: Scenario %s(%s): The estimates \n",
	     SceNam.Data(), ModNam.Data());
      myBlue->PrintEst();
    }else{
      //-- The relative uncertainties
      ModNam = "relative";
      for(Int_t i = 0; i < NumEst; i++){
	for(Int_t k = 0; k < NumUnc; k++){
	  ActCof[2] = XEst[i*(NumUnc+1)+k+1]/XEst[i*(NumUnc+1)] *
	    XEst[i*(NumUnc+1)+k+1]/XEst[i*(NumUnc+1)];
	  myBlue->SetRelUnc(i, k, &ActCof[0]);
	}
      }
      myBlue->FixInp();
      myBlue->SolveRelUnc(0.1);
      printf("... B_Peelles: Scenario %s(%s): The --modified-- estimates \n",
	     SceNam.Data(), ModNam.Data());
      myBlue->PrintEst();
      printf("... B_Peelles: Their value \n");
      myBlue->GetEstVal(LocEstVal); 
      myBlue->PrintMatrix(LocEstVal);
      printf("... B_Peelles: Their uncertainties \n");
      myBlue->GetEstUnc(LocEstUnc); 
      myBlue->PrintMatrix(LocEstUnc);
    }
    printf("\n");
    
    printf("... B_Peelles: Scenario %s(%s): The result \n",
	   SceNam.Data(), ModNam.Data());
    myBlue->PrintResult();
    printf("\n");
    
    printf("... B_Peelles: Scenario %s(%s): The chi-squared of the estimates\n",
	   SceNam.Data(), ModNam.Data());
    myBlue->PrintChiPro();
    printf("\n");
    
    printf("... B_Peelles: Scenario %s(%s): The parameters \n",
	   SceNam.Data(), ModNam.Data());
    myBlue->PrintParams();
    printf("\n");
    
    printf("... B_Peelles: Scenario %s(%s): The likelihood fit \n",
	   SceNam.Data(), ModNam.Data());
    FilNam = FilBas + "_" + SceNam + "_" + ModNam;
    myBlue->InspectLike(0, FilNam);
    printf("\n");
    
    myBlue->ReleaseInp();
  }

  // Delete Object
  delete myBlue; myBlue = NULL;
  LocEst->Delete()   ; LocEst    = NULL;
  LocEstVal->Delete(); LocEstVal = NULL;
  LocEstUnc->Delete(); LocEstUnc = NULL;
  return;
}
