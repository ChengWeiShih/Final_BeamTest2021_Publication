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

void B_ATLAS_CONF_2014_052(Int_t Flag = 0){
  
  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The LHC Combination with COMBINED   uncertainties
  //  1: The LHC Combination with INDIVIDUAL uncertainties
  //  2: The LHC Combination with INDIVIDUAL and absolute uncertainties
  //----------------------------------------------------------------------------
  
  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  2;
  static const Int_t NumUnc = 20;

  static const Int_t NumObs =  1;
  static const Int_t NCoUnc =  6;
  static const Int_t NVaSys =  6;
  static const Int_t NVaRho = 12;

 // The names of estimates, uncertainties and observables
  TString NamEst[NumEst] = {"  ATLAS", "    CMS"};
  TString NamUnc[NumUnc] = {
    " DaStat", " SiStat", "   Lumi", 
    "  I/FSR", " tW Gen", " tt Gen", "    PDF", "tW/ttOv", " Top pt", 
    "  BGMod", "    Z+j", 
    "   JESc", "   JESf", "   JESi", "   JESr", 
    " LepMod", "   METs", "   METr", "   btag", "   Pile"};
  TString NamObs[NumObs] = {"  Sigma"};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_ATLAS_CONF_2014_052: -------------------------------");
    printf("--------------------------------- \n");
    printf("... B_ATLAS_CONF_2014_052: The 2011 LHC Combination with");
    printf(" COMBINED   uncertainties Flag = %2i \n", Flag);
    NamUnc[0] = " DaStat";
    NamUnc[1] = " SiStat";
    NamUnc[2] = "   Lumi";
    NamUnc[3] = " Theory";
    NamUnc[4] = "Backgrd";
    NamUnc[5] = "   Jets";
    NamUnc[6] = " DetMod";
  }else if(Flag == 1){
    printf("... B_ATLAS_CONF_2014_052: -------------------------------");
    printf("-------------------------------------------------------------- \n");
    printf("... B_ATLAS_CONF_2014_052: The 2011 LHC Combination with");
    printf(" INDIVIDUAL uncertainties");
    printf(" Flag = %2i \n", Flag);
  }else if(Flag == 2){
    printf("... B_ATLAS_CONF_2014_052: -------------------------------");
    printf("-------------------------------------------------------------- \n");
    printf("... B_ATLAS_CONF_2014_052: The 2011 LHC Combination with");
    printf(" INDIVIDUAL and absolute uncertainties");
    printf(" Flag = %2i \n", Flag);
  }else{
    printf("... ATLAS_CONF_2014_052: Not implemented Flag = %2i \n", Flag);
    return;
  }

  // Character to set the file names
  char Buffer[100];
  sprintf(Buffer,"B_ATLAS_CONF_2014_052_%i",Flag);

  // Estimates, 0 == ATLAS, 
  //            1 == CMS
  // Fill percentage uncertainties as negative entries.
  // Full breakdown of individual categories, see Table 1
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //       0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19
    27.2, -7.1, -2.8, -3.7, -5.9,-11.0, -7.5, -2.5, -1.4,  0.0, -3.6,  0.0,-10.0, -5.0,  0.2, -0.7, -2.4, -4.1, -4.5, -8.4,  0.0,
    23.4, -8.1, -2.4, -3.0,-12.4,  0.0,-14.1, -1.7, -2.1, -0.4, -1.7, -2.6, -3.8,  0.0,  0.0, -0.9, -1.8, -0.4,  0.0, -0.9, -0.4
  };
  // Save default values
  Double_t XEstSav[LenXEst] = {0};
  for(Int_t i = 0; i<LenXEst; i++)XEstSav[i] = XEst[i];

  // The correlation assumptions for individual categories, see Table 1
  Double_t RhoVal[NumUnc] = {
    //       0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19
             0,    0, 0.31,    1,    0,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,  0.5,    0};
  // Save default values
  Double_t RhoSav[NumUnc] = {0};
  for(Int_t k = 0; k<NumUnc; k++)RhoSav[k] = RhoVal[k];

  // RhoVar = The correlation variations for various categories, see Table 3
  // RhoInd = The corresponding indices k>0 / k<0 == original/ combined sources
  //                         0    1  2    3   4    5    6   7   8   9   10  11
  Double_t RhoVar[NVaRho] = {0, 0.5, 0, 0.5,  0, 0.3, 0.5,  1,  0,  1, 0.5,  1};
  Int_t    RhoInd[NVaRho] = {2,   2, 5,   5, -3,  -3,  -4, -4, 18, 18,  -5, -5};

  // The default correlations for Table 3. Some can be taken from Table 1. The
  // ones preset to -1 are calculated below from the COMBINED uncertainty
  // combination, namely from the groups Theory modelling and Jets.
  Double_t RhoDef[NVaSys] = {RhoVal[2], RhoVal[5], -1, RhoVal[9], RhoVal[18], -1};

  //--Local structure for Blue output and systematic uncertainties Table 2
  TMatrixD* LocRes = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc = new TMatrixD(NumObs,1);
  TMatrixD* Values = new TMatrixD(NumObs,NVaRho+1);
  TMatrixD* Uncert = new TMatrixD(NumObs,NVaRho+1);
  //--End
  
  // Define formats for Figures and Latex file
  const TString ForVal = "%5.1f";
  const TString ForUnc = "%4.1f";
  const TString ForWei = "%5.3f";
  const TString ForRho = "%5.2f";
  const TString ForPul = ForWei;
  const TString ForChi = "%5.3f";
  const TString ForUni = "None";

  // Now do the flag dependent stuff
  if(Flag == 0){

    // The default kk == 0, plus the twelve cases of Rho variations
    for(Int_t kk = 0; kk<NVaRho+1; kk++){

      // Change the correlations for the --original-- sources if needed
      if(kk >0 && RhoInd[kk-1] > 0)RhoVal[RhoInd[kk-1]] = RhoVar[kk-1];

      // Calculate the six combined categories and their correlation overwrite
      // the original array
      Double_t ATLSum = 0;
      Double_t CMSSum = 0;
      Double_t RhoSum = 0;
      Double_t ATLVal = XEst[0];
      Double_t CMSVal = XEst[NumUnc+1];
      // Data stats    0 = 0
      // Simul stats   1 = 1
      // Luminosity    2 = 2
      // Theory Model  3 = 3-8
      // Background    4 = 9-10
      // Jets          5 = 11-14
      // Det modelling 6 = 15-19
      Int_t kmin = 0, kmax = 0;
      if(kk == 0){
	printf("... B_ATLAS_CONF_2014_052: Uncertainties and correlation");
	printf(" per category \n");
	printf("... B_ATLAS_CONF_2014_052: see bold face values in Table 1 \n");
      }
      for(Int_t ll = 0; ll<NCoUnc+1; ll++){
	if(ll == 0){kmin=0;kmax=0;
	}else if(ll == 1){kmin=1;kmax=1;
	}else if(ll == 2){kmin=2;kmax=2;
	}else if(ll == 3){kmin=3;kmax=8;
	}else if(ll == 4){kmin=9;kmax=10;
	}else if(ll == 5){kmin=11;kmax=14;
	}else if(ll == 6){kmin=15;kmax=19;
	}
	ATLSum = 0;
	CMSSum = 0;
	RhoSum = 0;
	for(Int_t k = kmin; k<kmax+1; k++){
	  ATLSum = ATLSum + XEst[k+1]*XEst[k+1];
	  CMSSum = CMSSum + XEst[k+1+NumUnc+1]*XEst[k+1+NumUnc+1];
	  RhoSum = RhoSum + RhoVal[k] * (ATLVal*XEst[k+1]/100) * (CMSVal*XEst[k+1+NumUnc+1]/100);
	}
	ATLSum = ATLVal * TMath::Sqrt(ATLSum)/100;
	CMSSum = CMSVal * TMath::Sqrt(CMSSum)/100;
	RhoSum = RhoSum / (ATLSum*CMSSum);
	
	// Update Arrays for Categories with more than one uncertainty
	if(ll > 2){
	  XEst[ll+1] = ATLSum;
	  XEst[ll+1+NumUnc+1] = CMSSum;
	  RhoVal[ll] = RhoSum;
	}

	// Change the correlations for the COMBINED sources if needed
	if(kk >0 && RhoInd[kk-1] == -ll)RhoVal[ll] = RhoVar[kk-1];

	// Print Table 1
	if(kk == 0){
	  printf("... B_ATLAS_CONF_2014_052: Category = %2i, ATLAS = %4.1f%%,", 
		 ll, 100*ATLSum/ATLVal);
	  printf(" CMS = %4.1f%%, Correlation = %4.2f \n", 
		 100*CMSSum/CMSVal, RhoVal[ll]);
	}
      }

      // Update default correlations for display of Table 3
      if(kk == 0){
	RhoDef[2] = RhoVal[3];
	RhoDef[5] = RhoVal[5];
      }

      // Create object
      Blue *myBlue = new Blue(NumEst, NumUnc);
      myBlue->PrintStatus();
      myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);
      if(kk > 0){
	myBlue->SetQuiet();
	if(kk == 1){
	  printf("... B_ATLAS_CONF_2014_052: Now do the systematic variations");
	  printf(" for Table 3 in quiet mode \n");
	}
      }
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
	myBlue->FillCor(k,RhoVal[k]);
      }
      
      // Solve and digest
      myBlue->SetRelUnc();
      for(Int_t k = NCoUnc+1; k<NumUnc; k++)myBlue->SetInActiveUnc(k);
      myBlue->FixInp();
      myBlue->SolveRelUnc(0.1);
      if(kk == 0){
	myBlue->PrintEst();
	printf("... B_ATLAS_CONF_2014_052: For the result see Table 2.");
	printf(" For the weights, the chi-squared,\n");
	printf("... B_ATLAS_CONF_2014_052: the total correlation and the");
	printf(" pulls, see text in Section 5. \n");
	myBlue->PrintResult();
	myBlue->PrintCompatEst();
	myBlue->PrintRho();
	myBlue->PrintPull();
	myBlue->LatexResult(Buffer);
	myBlue->DisplayResult(0,Buffer,ForVal,ForUnc);
      }

      // Store results
      myBlue->GetResult(LocRes);
      Values->operator()(0,kk) = LocRes->operator()(0,0);
      myBlue->GetUncert(LocUnc);
      Uncert->operator()(0,kk) = LocUnc->operator()(0,0);
      
      // Delete Object
      delete myBlue; myBlue = NULL;

      // Reset the input to the original values for a fresh start of next
      // systematic variation
      for(Int_t k = 0; k<NumUnc; k++)RhoVal[k] = RhoSav[k];
      for(Int_t i = 0; i<LenXEst; i++)XEst[i] = XEstSav[i];
    }

    // Now all is calculated look at the results of Table 3
    printf("... B_ATLAS_CONF_2014_052: -------------- Results in Table 3");
    printf("-------------------\n");
    printf("... B_ATLAS_CONF_2014_052: Category Default    Test-Rho");
    printf("   Shift-Val   Shift-Unc\n");
    Int_t ind = 0;
    for(Int_t kk = 0; kk<NVaRho/2; kk++){
      ind=2*kk;
      printf("... B_ATLAS_CONF_2014_052:       %2i", kk);
      printf("    %+3.1f", RhoDef[kk]);
      printf("   %+3.1f/%+3.1f", RhoVar[ind],RhoVar[ind+1]);
      printf("   %+3.1f/%+3.1f",
	     Values->operator()(0,ind+1)-Values->operator()(0,0),
	     Values->operator()(0,ind+2)-Values->operator()(0,0));
      printf("   %+3.1f/%+3.1f\n",
	     Uncert->operator()(0,ind+1)-Uncert->operator()(0,0),
	     Uncert->operator()(0,ind+2)-Uncert->operator()(0,0));
    }
    printf("... B_ATLAS_CONF_2014_052: ---------------------------------");
    printf("-------------------\n");
    
  }else{

    // Create object
    Blue *myBlue = new Blue(NumEst, NumUnc);
    myBlue->PrintStatus();
    myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);

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
      myBlue->FillCor(k,RhoVal[k]);
    }
    
    if(Flag == 1){

      // Get result with INDIVIDUAL uncertainties but relative uncertainties
      myBlue->SetRelUnc();
      myBlue->FixInp();
      myBlue->PrintEst();
      myBlue->SolveRelUnc(0.1);
      myBlue->PrintResult();
      myBlue->InspectLike(0,Buffer); 
      myBlue->PrintInspectLike();  

      myBlue->DisplayResult(0,Buffer,ForVal,ForUnc);

    }else if(Flag == 2){

      // Get result with INDIVIDUAL uncertainties and absolute uncertainties
      myBlue->FixInp();
      myBlue->PrintEst();
      myBlue->Solve();
      myBlue->PrintResult();

      myBlue->LatexResult(Buffer);
      myBlue->DisplayResult(0,Buffer,ForVal,ForUnc);

      // Scan the correlations
      myBlue->ReleaseInp();
      myBlue->FixInp();
      myBlue->SolveScaRho(0);
      myBlue->PrintScaRho(Buffer);
     }
    
    // Delete Object
    delete myBlue; myBlue = NULL;
  }

  return;
}
