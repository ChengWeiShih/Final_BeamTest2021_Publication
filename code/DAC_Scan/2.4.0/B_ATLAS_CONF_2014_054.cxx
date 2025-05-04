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

void B_ATLAS_CONF_2014_054(Int_t Flag = 0){
  
  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  // Combination of LHC top-quark pair cross section
   //  0: The LHC Combination
  //----------------------------------------------------------------------------
  
  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  2;
  static const Int_t NumUnc = 23;
  static const Int_t NumObs =  1;
  static const Int_t NumSys = 19;

 // The names of estimates, uncertainties and observables
  TString NamEst[NumEst] = {"  ATLAS", "    CMS"};
  TString NamUnc[NumUnc] = {
    "   Stat", "   Trig", "   LepS", "   LepI", "   JetR", "   JetI", 
    "   btag", "   Pile", "   JESu", "   JESi", "   JESc", "   JESf", 
    "   JESb", "  Scale", "    Rad", "  MC+PS", "    PDF", "    Z+j", 
    "   LepM", "  DiBos", "   STop", "    Vdm", "   Lumi"};
  TString NamObs[NumObs] = {"  Sigma"};

 // The names of the systematic variations
  TString NamSys[NumSys] = {"    Vdm 0.6", "    Vdm 0.3", "    Vdm 0.0",
			    "G/S 0.0 0.0", "G/S 0.0 0.5", "G/S 0.0 1.0",
			    "G/S 0.5 0.0", "G/S 0.5 0.5", "G/S 0.5 1.0",
			    "G/S 1.0 0.0", "G/S 1.0 0.5", "G/S 1.0 1.0",
			    "    PDF 0.5", "    PDF 0.0",
			    "   JESi 0.5", "   JESc 1.0", "   JESf 0.0",
			    "   btag 1.0", "DiB/STo 0.0"};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_ATLAS_CONF_2014_054: -------------------------------");
    printf("--------------------------------- \n");
    printf("... B_ATLAS_CONF_2014_054: The 2011 LHC Combination with");
    printf(" COMBINED   uncertainties Flag = %2i \n", Flag);
  }else{
    printf("... ATLAS_CONF_2014_054: Not implemented Flag = %2i \n", Flag);
    return;
  }
  
  // Character to set the file names
  char Buffer[100];
  sprintf(Buffer,"B_ATLAS_CONF_2014_054_%i",Flag);

  // Estimates, 0 == ATLAS, 
  //            1 == CMS
  // Full breakdown of individual categories, see Table 1
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //       0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22
    242.4, 1.7, 0.4, 1.2, 1.7, 1.2, 0.1, 1.0, 0.0, 0.6, 0.6, 0.3, 0.9, 0.1, 0.7, 0.0, 3.0, 2.7, 0.1, 0.8, 0.3, 2.0, 2.9, 6.9,
    239.0, 2.6, 3.6, 0.2, 4.0, 3.0, 0.0, 1.7, 2.0, 4.3, 0.6, 0.1, 2.9, 0.0, 5.6, 3.8, 3.3, 0.5, 1.5, 1.9, 0.5, 2.3, 5.0, 3.6
  };
  // The correlation assumptions for individual categories, see Table 1
  Double_t RhoVal[NumUnc] = {
    //0  1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0.5,   1,   0, 0.5,   0, 0.5,   1,   0,   0,   1,   1,   1,   0};

  //--Local structure for Blue output and systematic uncertainties Table 2
  TMatrixD* LocRes = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc = new TMatrixD(NumObs,1);
  TMatrixD* Values = new TMatrixD(NumObs,NumSys+1);
  TMatrixD* Uncert = new TMatrixD(NumObs,NumSys+1);
  //--End

  // Define formats for Figures and Latex file
  const TString ForVal = "%6.1f";
  const TString ForUnc = "%4.1f";
  const TString ForWei = "%4.3f";
  const TString ForRho = ForWei;
  const TString ForPul = ForWei;
  const TString ForChi = "%5.3f";
  const TString ForUni = "pb";

  // Create object
  Blue *myBlue = new Blue(NumEst, NumUnc);
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
    myBlue->FillCor(k,RhoVal[k]);
  }
  

  if(Flag == 0){

    // Solve and save result
    myBlue->FixInp();
    myBlue->Solve();
    myBlue->GetResult(LocRes);
    Values->operator()(0,0) = LocRes->operator()(0,0);
    myBlue->GetUncert(LocUnc);
    Uncert->operator()(0,0) = LocUnc->operator()(0,0);

    // Digest
    myBlue->PrintEst();
    printf("... B_ATLAS_CONF_2014_054: For the result see Table 2.");
    printf(" For the weights, the chi-squared,\n");
    printf("... B_ATLAS_CONF_2014_054: the total correlation and the pulls");
    printf(" see text in Section 5. \n");
    myBlue->PrintResult();
    myBlue->PrintCompatEst();
    myBlue->PrintRho();
    myBlue->PrintPull();
    myBlue->LatexResult(Buffer);
    myBlue->DisplayResult(0,Buffer,ForVal,ForUnc);

    // Scan the correlations
    myBlue->ReleaseInp();
    myBlue->FixInp();
    myBlue->SolveScaRho(0);
    myBlue->PrintScaRho(Buffer);
    
    // Do systematic variations store all results
    myBlue->ResetInp();
    myBlue->SetQuiet();
    for(Int_t kk = 0; kk<NumSys; kk++){

      if(kk == 0){
	// Vdm: k = 21, rho = 0.6
	myBlue->SetRhoValUnc(21,0.6);
      }else if(kk ==  1){
	// Vdm: k = 21, rho = 0.3
	myBlue->SetRhoValUnc(21,0.3);
      }else if(kk ==  2){
	// Vdm: k = 21, rho = 0.0
	myBlue->SetRhoValUnc(21,0.0);
      }else if(kk ==  3){
	// Generator: k = 15, rho = 0.0
	//     Scale: k = 13, rho = 0.0
	myBlue->SetRhoValUnc(15,0.0);
	myBlue->SetRhoValUnc(13,0.0);
      }else if(kk ==  4){
	// Generator: k = 15, rho = 0.0
	//     Scale: k = 13, rho = 0.5
	myBlue->SetRhoValUnc(15,0.0);
	myBlue->SetRhoValUnc(13,0.5);
      }else if(kk ==  5){
	// Generator: k = 15, rho = 0.0
	//     Scale: k = 13, rho = 1.0
	myBlue->SetRhoValUnc(15,0.0);
	myBlue->SetRhoValUnc(13,1.0);
      }else if(kk ==  6){
	// Generator: k = 15, rho = 0.5
	//     Scale: k = 13, rho = 0.0
	myBlue->SetRhoValUnc(15,0.5);
	myBlue->SetRhoValUnc(13,0.0);
      }else if(kk ==  7){
	// Generator: k = 15, rho = 0.5
	//     Scale: k = 13, rho = 0.5
	myBlue->SetRhoValUnc(15,0.5);
	myBlue->SetRhoValUnc(13,0.5);
      }else if(kk ==  8){
	// Generator: k = 15, rho = 0.5
	//     Scale: k = 13, rho = 1.0
	myBlue->SetRhoValUnc(15,0.5);
	myBlue->SetRhoValUnc(13,1.0);
      }else if(kk ==  9){
	// Generator: k = 15, rho = 1.0
	//     Scale: k = 13, rho = 0.0
	myBlue->SetRhoValUnc(15,1.0);
	myBlue->SetRhoValUnc(13,0.0);
      }else if(kk == 10){
	// Generator: k = 15, rho = 1.0
	//     Scale: k = 13, rho = 0.5
	myBlue->SetRhoValUnc(15,1.0);
	myBlue->SetRhoValUnc(13,0.5);
      }else if(kk == 11){
	// Generator: k = 15, rho = 1.0
	//     Scale: k = 13, rho = 1.0
	myBlue->SetRhoValUnc(15,1.0);
	myBlue->SetRhoValUnc(13,1.0);
      }else if(kk == 12){
	// PDF: k = 16, rho = 0.5
	myBlue->SetRhoValUnc(16,0.5);
      }else if(kk == 13){
	// PDF: k = 16, rho = 0.0
	myBlue->SetRhoValUnc(16,0.0);
      }else if(kk == 14){
	// JESi: k =  9, rho = 0.5
	myBlue->SetRhoValUnc( 9,0.5);
      }else if(kk == 15){
	// JESc: k = 10, rho = 1.0
	myBlue->SetRhoValUnc(10,1.0);
      }else if(kk == 16){
	// JESf: k = 11, rho = 0.0
	myBlue->SetRhoValUnc(11,0.0);
      }else if(kk == 17){
	// btag: k = 6, rho = 1.0
	myBlue->SetRhoValUnc( 6,1.0);
      }else if(kk == 18){
	// DiBos: k = 19, rho = 0.0
	//  STop: k = 20, rho = 0.0
	myBlue->SetRhoValUnc(19,0.0);
	myBlue->SetRhoValUnc(20,0.0);
      }

      // Solve
      myBlue->FixInp();
      myBlue->Solve();

      // Store results
      myBlue->GetResult(LocRes);
      Values->operator()(0,kk+1) = LocRes->operator()(0,0);
      myBlue->GetUncert(LocUnc);
      Uncert->operator()(0,kk+1) = LocUnc->operator()(0,0);
      
      // Reset
      myBlue->ResetInp();
    }

    // Print the findings
    printf("... B_ATLAS_CONF_2014_054: -------------------------------------");
    printf("----------------------------------\n");
    printf("... B_ATLAS_CONF_2014_054: ------------------------");
    printf(" Results of systematics -----------------------\n");
    printf("... B_ATLAS_CONF_2014_054: -------------------------------------");
    printf("----------------------------------\n");
    printf("... B_ATLAS_CONF_2014_054: N =  0 -  2: Vdm Scan variation");
    printf(" see Table 2\n");
    printf("... B_ATLAS_CONF_2014_054: N =  3 - 11: Generator and");
    printf(" parton shower vs Scale variation, see Table 3\n");
    printf("... B_ATLAS_CONF_2014_054: N = 12 - 13: PDF variation,");
    printf(" see Section 6.2\n");
    printf("... B_ATLAS_CONF_2014_054: N = 14 - 16: JESx variation,");
    printf(" see Table 3\n");
    printf("... B_ATLAS_CONF_2014_054: N =      17: b-tagging variation");
    printf(" see Section 6.4\n");
    printf("... B_ATLAS_CONF_2014_054: N =      18: Background variation,");
    printf(" see Section 6.5\n");
    printf("... B_ATLAS_CONF_2014_054: -------------------------------------");
    printf("----------------------------------\n");
    printf("\n");
    printf("... B_ATLAS_CONF_2014_054:               Scenario");
    printf("  Value / Uncert\n");
    printf("... B_ATLAS_CONF_2014_054: -----------------------------------");
    printf("---\n");
    for(Int_t kk = 0; kk<NumSys; kk++){
      printf("... B_ATLAS_CONF_2014_054: N=%2i:      %s", kk,
	     NamSys[kk].Data());
      printf("   %+3.1f / %+3.1f\n", 
	     Values->operator()(0,kk+1)-Values->operator()(0,0),
	     Uncert->operator()(0,kk+1)-Uncert->operator()(0,0));
    }
  }  

  // Delete Object
  delete myBlue; myBlue = NULL;
  return;
}
