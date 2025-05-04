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

void B_ATLAS_CONF_2013_098(Int_t Flag = 0){
  
  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  // Combination of LHC single top cross-sections
  //  0: The LHC Combination with  COMBINED categories 
  //     [Shows how to use SetRelUnc()]
  //  1: The LHC Combination with INDIVIDUAL categories and the sys variations
  //     [Shows how to use SetRhoValUnc()]
  //  2: The LHC Combination with INDIVIDUAL categories but changes assumptions
  //     on what is a relative uncertainty
  //  3: The calculation of the combined categories in Table 1
  //     [Shows how to Activate and De-activate uncertainty sources]
  //----------------------------------------------------------------------------
  
  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  2;
  static const Int_t NCoUnc =  6;
  static const Int_t NFuUnc = 22;
  Int_t NumUnc = NFuUnc;

  // Preset according to Flag
  if(Flag == 0){
    NumUnc = NCoUnc;
    printf("... B_ATLAS_CONF_2013_098: -------------------------------");
    printf("------------------------------ \n");
    printf("... B_ATLAS_CONF_2013_098: The 2011 LHC Combination with");
    printf("   COMBINED categories Flag = %2i \n", Flag);
  }else if(Flag == 1){
    printf("... B_ATLAS_CONF_2013_098: -------------------------------");
    printf("------------------------------ \n");
    printf("... B_ATLAS_CONF_2013_098: The 2011 LHC Combination with");
    printf(" INDIVIDUAL categories Flag = %2i \n", Flag);
  }else if(Flag == 2){
    printf("... B_ATLAS_CONF_2013_098: -------------------------------");
    printf("-------------------------------------------------------------- \n");
    printf("... B_ATLAS_CONF_2013_098: The 2011 LHC Combination with");
    printf(" INDIVIDUAL categories and only relative uncertainties");
    printf(" Flag = %2i \n", Flag);
  }else if(Flag == 3){
    printf("... B_ATLAS_CONF_2013_098: -------------------------------");
    printf("------------------------------------------------ \n");
    printf("... B_ATLAS_CONF_2013_098: Calcultate input for the COMBINED");
    printf(" categories using the full breakdown Flag = %2i \n", Flag);
  }else{
    printf("... ATLAS_CONF_2013_098: Not implemented Flag = %2i \n", Flag);
    return;
  }

  // Character to set the file names
  char Buffer[100];

  // Estimates, 0 == ATLAS, 
  //            1 == CMS
  // Fill percentage uncertainties as negative entries.

  // Full (Fu) breakdown of individual categories, see Table 1
  static const Int_t LenXFuEst = NumEst * (NFuUnc+1);
  Double_t XFuEst[LenXFuEst] = {
    //       0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21
    95.1, -2.4, -2.9, -3.0, -2.0, -9.1, -2.8, -7.1, -3.3, -0.8, -7.7, -3.0, -1.6, -3.1,  0.0, -8.5, -2.3, -1.6,  0.0, -4.1,  0.0, -2.2, -2.1,
    80.1, -7.1, -2.2, -4.1, -1.6, -3.1, -4.6, -5.5,  0.0,  0.0, -6.8, -0.7, -2.1, -0.9, -4.5, -4.6, -1.0,  0.0, -0.5,  0.0, -5.1,  0.0,  0.0
  };
  Double_t XFuCor[NFuUnc] = {
             0,    0,    1,    0,    1,    1,    1,    0,    0,    0,    0,    1,    0,    0,   0.5,    0,    0,    0,   0,    0,    0,    0};

  // Combined (Co) categories, bold numbers in Table 1. 
  // The values and correlations have been calculated from full info see Flag == 1
  static const Int_t LenXCoEst = NumEst * (NCoUnc+1);
  Double_t XCoEst[LenXCoEst] = { 
    //             0             1            2            3            4            5
    // The values in the Table
    //95.1,        -3.8,        -3.6,       -12.3,        -8.3,        -3.5,       -10.3,
    //80.1,        -7.5,        -4.4,        -7.8,        -6.8,        -5.0,        -6.9
    // The calculated values from the full breakdown
    95.1, TMath::Sqrt(12.82), TMath::Sqrt(11.76), TMath::Sqrt(138.0), TMath::Sqrt(61.76), TMath::Sqrt(11.01), TMath::Sqrt(96.01), 
    80.1, TMath::Sqrt(35.45), TMath::Sqrt(12.43), TMath::Sqrt(39.15), TMath::Sqrt(29.98), TMath::Sqrt(16.34), TMath::Sqrt(31.07)
  };
  Double_t XCoCor[NCoUnc] = {
    // The values in the Table
    //              0,        0.78,        0.83,           0,        0.19,       0.27
    // The calculated values from the full breakdown
                    0,      0.7751,      0.8305,           0,      0.1908,     0.2727
  };

  // Local structure for Blue output and systematic uncertainties Table 3
  Int_t iok = 0;
  static const Int_t NumObs = 1;
  static const Int_t LenRes = NumObs*(NFuUnc+1);
  Double_t ValDef[LenRes] = {0}, UncDef[NumObs] = {0};
  Double_t ValAct[LenRes] = {0}, UncAct[NumObs] = {0};

  static const Int_t NumSys = 8;
  Double_t ValDif[NumSys] = {0}, UncDif[NumSys] = {0};
  
  // Define formats for Figures and Latex file
  const TString ForVal = "%5.1f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%5.3f";
  const TString ForRho = ForVal;
  const TString ForPul = ForVal;
  const TString ForChi = "%5.3f";
  const TString ForUni = "GeV";

  // Create object
  Blue *myBlue = new Blue(NumEst, NumUnc);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

  // Fill estimates and correlations
  Int_t ind = 0;
  if(Flag == 0){
    for(Int_t i = 0; i<NumEst; i++){
      myBlue->FillEst(i,&XCoEst[ind]);
      ind = ind + NCoUnc + 1;
    }
    for(Int_t k = 0; k<NCoUnc; k++){
      myBlue->FillCor(k,XCoCor[k]);
    }
    myBlue->SetRelUnc();
    myBlue->SetNotRelUnc(0);
  }else{
    for(Int_t i = 0; i<NumEst; i++){
      myBlue->FillEst(i,&XFuEst[ind]);
      ind = ind + NFuUnc + 1;
    }
    for(Int_t k = 0; k<NFuUnc; k++){
      myBlue->FillCor(k,XFuCor[k]);
    }
  }  
  if(Flag == 1){
    // Get result with full breakdown
    myBlue->SetRelUnc();
    myBlue->SetNotRelUnc(0);
    myBlue->SetNotRelUnc(1);
  }else if(Flag == 2){
    // Get result with full breakdown but some changed sources
    myBlue->SetRelUnc();
  }else if(Flag == 3){
    // Show how to get the combined values  
    // Do not solve in this case
    for(Int_t k = 0; k<NFuUnc; k++){myBlue->SetInActiveUnc(k);};
    myBlue->SetActiveUnc(0);
    myBlue->SetActiveUnc(1);
    myBlue->FixInp();
    myBlue->PrintCov();
    myBlue->PrintRho();
    myBlue->ResetInp();
    
    for(Int_t k = 0; k<NFuUnc; k++){myBlue->SetInActiveUnc(k);};
    myBlue->SetActiveUnc(2);
    myBlue->SetActiveUnc(3);
    myBlue->FixInp();
    myBlue->PrintCov();
    myBlue->PrintRho();
    myBlue->ResetInp();

    for(Int_t k = 0; k<NFuUnc; k++){myBlue->SetInActiveUnc(k);};
    myBlue->SetActiveUnc(4);
    myBlue->SetActiveUnc(5);
    myBlue->SetActiveUnc(6);
    myBlue->SetActiveUnc(7);
    myBlue->SetActiveUnc(8);
    myBlue->FixInp();
    myBlue->PrintCov();
    myBlue->PrintRho();
    myBlue->ResetInp();
    
    for(Int_t k = 0; k<NFuUnc; k++){myBlue->SetInActiveUnc(k);};
    myBlue->SetActiveUnc(9);
    myBlue->SetActiveUnc(10);
    myBlue->FixInp();
    myBlue->PrintCov();
    myBlue->PrintRho();
    myBlue->ResetInp();
    
    for(Int_t k = 0; k<NFuUnc; k++){myBlue->SetInActiveUnc(k);};
    myBlue->SetActiveUnc(11);
    myBlue->SetActiveUnc(12);
    myBlue->SetActiveUnc(13);
    myBlue->FixInp();
    myBlue->PrintCov();
    myBlue->PrintRho();
    myBlue->ResetInp();
    
    for(Int_t k = 0; k<NFuUnc; k++){myBlue->SetInActiveUnc(k);};
    myBlue->SetActiveUnc(14);
    myBlue->SetActiveUnc(15);
    myBlue->SetActiveUnc(16);
    myBlue->SetActiveUnc(17);
    myBlue->SetActiveUnc(18);
    myBlue->SetActiveUnc(19);
    myBlue->SetActiveUnc(20);
    myBlue->SetActiveUnc(21);
    myBlue->FixInp();
    myBlue->PrintCov();
    myBlue->PrintRho();
  }

  // Now solve if needed
  if(Flag != 3){
    myBlue->FixInp();
    printf("... B_ATLAS_CONF_2013_098: Estimators after last user call");
    printf(" to FixInp() before SolveRelUnc()\n");
    myBlue->PrintEst();
    myBlue->SolveRelUnc(0.1);
    printf("... B_ATLAS_CONF_2013_098: Estimators after last call to");
    printf(" FixInp() from SolveRelUnc()\n");
    myBlue->PrintEst();

    // Digest results
    if(Flag == 0){
      printf("----------------------------------------------------");
      printf("------------------- \n");
      printf("... B_ATLAS_CONF_2013_098: The results see Table 1+2");
      printf(" and text on page 7 \n");
      printf("----------------------------------------------------");
      printf("------------------- \n");
      myBlue->PrintCofRelUnc();
    }
    myBlue->PrintCompatEst();
    myBlue->PrintRho();
    myBlue->PrintCov();
    myBlue->PrintResult();
    myBlue->PrintPull();
    if(Flag == 0){
      // Inspect the likelihood
      sprintf(Buffer,"B_ATLAS_CONF_2013_098_%i",Flag);
      myBlue->InspectLike(0,Buffer); 
      myBlue->PrintInspectLike();  
      printf("----------------------------------------------------");
      printf("--------------------------- \n");
      printf("... B_ATLAS_CONF_2013_098: For comparison Table 1+2");
      printf(" with Absolute Uncertainties \n");
      printf("----------------------------------------------------");
      printf("--------------------------- \n");
      myBlue->ReleaseInp();
      myBlue->SetNotRelUnc();
      myBlue->FixInp();
      myBlue->Solve();
      myBlue->PrintResult();
    }
	
    // Systematic variations. There are four pairs of tests see Table 3
    // Use ReleaseInp() to not break SetRelUnc() from above
    if(Flag == 1){
      // Safe the Default values
      iok = myBlue->GetResult(ValDef);
      iok = myBlue->GetUncert(UncDef);
      if(iok == 1){
	printf("... B_ATLAS_CONF_2013_098: The default result is:");
	printf(" %5.3f +- %5.3f \n", ValDef[0], UncDef[0]);
      }else{
	printf("... B_ATLAS_CONF_2013_098: Error when retrieving result %2i\n", 
	       iok);
      }
      myBlue->ReleaseInp();

      // Now the eight cases
      for(Int_t ll = 0; ll<8; ll++){
	if(ll == 0){
	  // Lumi Calibration k = 2 (rho = 0.5, 0)
	  printf("... B_ATLAS_CONF_2013_098: Syst for the luminosity \n");
	  myBlue->SetRhoValUnc( 2,0.5);
	}else if(ll == 1){
	  myBlue->SetRhoValUnc( 2,0.0);
	}else if(ll == 2){
	  printf("... B_ATLAS_CONF_2013_098: Syst for simulation");
	  printf(" and modelling \n");
	  // Simulation and modelling k = 4, 5, 6, 11 (rho = 0.5, 0)
	  myBlue->SetRhoValUnc( 4,0.5);
	  myBlue->SetRhoValUnc( 5,0.5);
	  myBlue->SetRhoValUnc( 6,0.5);
	  myBlue->SetRhoValUnc(11,0.5);
	}else if(ll == 3){
	  myBlue->SetRhoValUnc( 4,0.0);
	  myBlue->SetRhoValUnc( 5,0.0);
	  myBlue->SetRhoValUnc( 6,0.0);
	  myBlue->SetRhoValUnc(11,0.0);
	}else if(ll == 4){
	  printf("... B_ATLAS_CONF_2013_098: Syst for the JES \n");
	  // JES k = 9 (rho = 0.5, 1)
	  myBlue->SetRhoValUnc( 9,0.5);
	}else if(ll == 5){
	  myBlue->SetRhoValUnc( 9,1.0);
	}else if(ll == 6){
	  // b-tagging k = 14 (rho = 0, 1)
	  printf("... B_ATLAS_CONF_2013_098: Syst for b-tagging \n");
	  myBlue->SetRhoValUnc(14,0.0);
	}else if(ll == 7){
	  myBlue->SetRhoValUnc(14,1.0);
	}
	myBlue->FixInp();
	myBlue->SolveRelUnc(0.1);
	iok = myBlue->GetResult(ValAct); iok = myBlue->GetUncert(UncAct);
	ValDif[ll] = ValAct[0]-ValDef[0];
	UncDif[ll] = UncAct[0]-UncDef[0];
	myBlue->ReleaseInp();
	myBlue->SetNotRhoValUnc();
      }

      // Report the results, see table 3
      printf("------------------------------------------------------------ \n");
      printf("... B_ATLAS_CONF_2013_098: The systematic variations Table 3 \n");
      printf("------------------------------------------------------------ \n");
      for(Int_t ll = 0; ll<8; ll++){
	if(ll == 0){
	  printf("... B_ATLAS_CONF_2013_098:");
	  printf("    Luminosity calibration -- rho = 0.5 / 0.0,");
	}else if(ll == 2){
	  printf("... B_ATLAS_CONF_2013_098:");
	  printf(" Simulation and modelling  -- rho = 0.5 / 0.0,");
	}else if(ll == 4){
	  printf("... B_ATLAS_CONF_2013_098:");
	  printf("                      JES  -- rho = 0.5 / 1.0,");
	}else if(ll == 6){
	  printf("... B_ATLAS_CONF_2013_098:");
	  printf("                b-tagging  -- rho = 0.0 / 1.0,");
	}
	if(ll == 0 || ll == 2 || ll == 4 || ll == 6){
	  printf(" Shift: %4.1f / %4.1f,", ValDif[ll], ValDif[ll+1]);
	  printf(" Uncertainty: %4.1f / %4.1f \n",UncDif[ll], UncDif[ll+1]);
	}
      }
    }
  }
  // Delete Object
  delete myBlue; myBlue = NULL;
  return;
}
