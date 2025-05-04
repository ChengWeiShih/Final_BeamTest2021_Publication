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

void B_EPJC_74_3109(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The combination of the top quark pole mass
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =   2;
  static const Int_t NumUnc =   6;
  static const Int_t NumObs =   1;

  // The names
  TString NamEst[NumEst] = {"7TeV", "8TeV"};
  TString NamUnc[NumUnc] = {"   Stat", "   Syst", 
			    "   Lumi", "   Beam",
			    "   PDFs", "   QCDS"};
  TString NamObs[NumObs] = {"   mtop"};

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_EPJC_74_3109: -------------------------------- \n");
    printf("... B_EPJC_74_3109: The ATLAS combination, Flag = %2i \n",
	   Flag);
  }else{
    printf("... B_EPJC_74_3109: ------------------------ \n");
    printf("... EPJC_74_3109: Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Estimates 0-1 == 7TeV, 8TeV
  // Values from CT10 in text
  // Uncertainties from Table 7
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //       0     1     2     3     4            5
    //    Stat, Syst, Lumi, Beam, PDFs,        QCDS
    171.4, 0.6,  0.8,  0.7,  0.7,  1.8, (0.9+1.2)/2.,
    174.1, 0.3,  0.9,  1.2,  0.6,  1.7, (0.9+1.3)/2.
  };

  // The ATLAS result
  const Double_t MasVal = 172.9;
  const Double_t UncHig =   2.5;
  const Double_t UncLow =   2.6;

  // The correlations          0     1     2     3     4     5
  //                        Stat, Syst, Lumi, Beam, PDFs, QCDS
  Double_t RhoVal[NumUnc] = {0.0,  1.0,  1.0,  1.0,  1.0,  1.0};

  // Assumptions on correlation for "Analysis systematics" k=1
  const Int_t NumRho = 6;
  Double_t    RhoSys[NumRho] = {-1.0, -0.6, -0.2, 0.2, 0.6, 1.0};

  //--Local structure for Blue output and systematic uncertainties Table 2
  TMatrixD* LocRes = new TMatrixD(NumObs,NumUnc+1);
  TMatrixD* LocUnc = new TMatrixD(NumObs,1);
  TMatrixD* Values = new TMatrixD(NumRho,1);
  TMatrixD* Uncert = new TMatrixD(NumRho,1);
  //--End

  // Declare names for output files
  TString FilBas = "B_EPJC_74_3109";
  TString FilNam;
  char Buffer[100];

  // Define formats for Figures and Latex file
  const TString ForVal = "%6.1f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%6.2f";
  const TString ForRho = "%5.2f";
  const TString ForPul = ForVal;
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
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }

  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++){
      myBlue->FillCor(k, RhoVal[k]);
  }

  // Fix and print estimates
  myBlue->FixInp();
  myBlue->PrintEst();

  // Loop over correlations for k=1 from 0 to 1 in steps of 0.2
  myBlue->SetQuiet();
  for(Int_t i = 0; i < NumRho; i++){
    printf("... EPJC_74_3109: Next correlation for source k=1 = %4.2f \n",
	   RhoSys[i]);
    myBlue->ReleaseInp();
    myBlue->SetRhoValUnc(1, RhoSys[i]);
    myBlue->FixInp();
    myBlue->Solve();

    // Digest results
    myBlue->PrintChiPro();

    // Get the output into tex files and plots
    if(i == 0){
      sprintf(Buffer,"%1i", Flag);
      FilNam = &Buffer[0];       
      FilNam = FilBas + "_" + FilNam;
      myBlue->LatexResult(FilNam);
      myBlue->DisplayResult(0,FilNam);
      myBlue->PrintResult();
    }

    // Store results in local structures
    myBlue->GetResult(LocRes);
    Values->operator()(i,0) = LocRes->operator()(0,0);
    myBlue->GetUncert(LocUnc);
    Uncert->operator()(i,0) = LocUnc->operator()(0,0);
  }

  // Write out what we got
  printf("... EPJC_74_3109:\n");
  printf("... EPJC_74_3109: Combined result as a function of rho(k=1) \n");
  printf("... EPJC_74_3109: See Table 6. Input is from Table 7 \n");
  for(Int_t i = 0; i < NumRho; i++){
    printf("... EPJC_74_3109: rho(k=1) = %+4.2f => mtop +- dmtop =", RhoSys[i]);
    printf(" %5.1f +- %5.1f \n", Values->operator()(i,0),
	   Uncert->operator()(i,0));
  }

  // The final findings
  printf("... EPJC_74_3109: \n");
  printf("... EPJC_74_3109: To achieve a combined value close to the quoted");
  printf(" result of %5.1f + %3.1f - %3.1f,\n", MasVal, UncHig, UncLow);
  printf("... EPJC_74_3109: the indvidual measurements need to be anti-");
  printf("correlated for the -- Analysis systematics--. \n");
  printf("... EPJC_74_3109:\n");

  //Delete Object
  delete myBlue;  myBlue = NULL;  
}
