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
#include "Riostream.h"
#include "TString.h"
#include "TMatrixD.h"

void ForttoBlue(TString FilNam = "EPJC_72_2046Fort", 
		TString ForVal = "%5.2f",
		TString ForUnc = "%4.2f",
		TString ForRho = "%4.2f"){

  // Dummy dimension extend if this is too small
  const Int_t DumDim = 500;

  // Fill TString variables
  TString FilInp = FilNam +".in";
  TString FilOut = "B_" + FilNam + ".cxx";
  TString FilSte = "B_" + FilNam + ".inp";
  TString FilCod = "B_" + FilNam;

  printf("... ForttoBlue(): I will convert the Fortran code input file '%s' \n", FilInp.Data());
  printf("... ForttoBlue(): into a function '%s' to be used with the Blue software \n", FilOut.Data());

  // Print warning
  printf("... ForttoBlue(): ---------------- WARNING !!! ---------------- \n");
  printf("... ForttoBlue(): This software has been tested on a number of input files \n");
  printf("... ForttoBlue(): However, is not very robust against additional content the file may have. \n");
  printf("... ForttoBlue(): It will work, provided the input file follows the following rules: \n");
  printf("... ForttoBlue(): 1) The file starts with three irrelevant lines \n");
  printf("... ForttoBlue(): 2) Line 4 starts with the numbers of \n");
  printf("... ForttoBlue():    observables, estimates und uncertainties \n");
  printf("... ForttoBlue():    and has an additnal string at the end \n");
  printf("... ForttoBlue(): 3) Line 5 starts with the names of the observables \n");
  printf("... ForttoBlue():    and has an additnal string at the end \n");
  printf("... ForttoBlue(): 4) Line 6 only contains the names of the uncertainties \n");
  printf("... ForttoBlue(): 5) Line 7ff contain the names of the estimate and observable, \n");
  printf("... ForttoBlue():    the value and the uncertainties of the estimate \n");
  printf("... ForttoBlue(): 6) The following lines contain the correlation matrices \n");
  printf("... ForttoBlue():    where the first row ends with an additional string \n");
  printf("... ForttoBlue(): 7) The last line of the file only contains a ! \n");
  printf("... ForttoBlue(): ---------------- WARNING !!! ---------------- \n");

  printf("... ForttoBlue(): ---------------- Example  !!! ---------------- \n");
  printf("... ForttoBlue(): This software works e.g with EPJC_72_2046Fort.in \n");
  printf("... ForttoBlue(): The produced function B_EPJC_72_2046Fort.cxx \n");
  printf("... ForttoBlue(): reproduces the result of B_EPJC_72_2046(0). \n");

  // Open input file
  FILE *fpinp = fopen(FilInp.Data(),"r");   
  if(fpinp == NULL){
    printf("... ForttoBlue(): File %s does not exist \n", FilInp.Data());
    return;
  }else{
    printf("... ForttoBlue(): Input file %s has been openend \n", FilInp.Data());
  }

  // Dummy return variable for sprintf and fscanf
  Int_t IRet = 0;
  if(IRet == 1)return;
  
  // Dummy return variable for fgets
  char* ret;

  // Set variables needed for read 
  char buf[200];
  TString NamBuf = " ";
  Float_t valdum = 0;

  // Read from input 
  ret = fgets(buf,200,fpinp);
  if(ret == NULL)return;
  ret = fgets(buf,200,fpinp);
  ret = fgets(buf,200,fpinp);
  Int_t NumEst = 1, NumUnc = 1, NumObs = 1;

  IRet = fscanf(fpinp,"%d %d %d",&NumObs,&NumEst,&NumUnc);
  printf("... ForttoBlue(): On Line 4 I found NumObs NumEst NumUnc: %i %i %i \n",NumObs,NumEst,NumUnc);
  ret = fgets(buf,200,fpinp);

  // Get the names of the observables
  TString NamObs[DumDim];
  printf("... ForttoBlue(): NamObs =");
  for(Int_t n = 0; n < NumObs; n++){
    IRet = fscanf(fpinp,"%s",&buf[0]);
    NamObs[n] = &buf[0];
    NamObs[n] = NamObs[n](1, NamObs[n].Length()-2);
    printf(" %s", NamObs[n].Data());
  }
  printf("\n");
  ret = fgets(buf,200,fpinp);
  
  // Get the names of the uncertanties
  //TString NamUnc[NumUnc];
  TString NamUnc[DumDim];
  printf("... ForttoBlue(): NamUnc =");
  for(Int_t k = 0; k < NumUnc; k++){
    IRet = fscanf(fpinp,"%s",&buf[0]);
    NamUnc[k] = &buf[0];
    NamUnc[k] = NamUnc[k](1, NamUnc[k].Length()-2);
    printf(" %s", NamUnc[k].Data());
  }
  printf("\n");
  
  // Book IWhichObs 
  Int_t IWhichObs[NumEst];
  for(Int_t i = 0; i < NumEst; i++){IWhichObs[i] = 0;};
  
  // Get names of the estimates the related observables and finally 
  // the values of uncertainties
  //  TString NamEst[NumEst];
  TString NamEst[DumDim];
  Int_t LenXEst = NumEst *(NumUnc+1);
  Double_t XEst[LenXEst];

  Int_t ind = 0, IObs = -1;
  for(Int_t i = 0; i < NumEst; i++){
    // The estimate
    IRet = fscanf(fpinp,"%s",&buf[0]);
    NamEst[i] = &buf[0];
    NamEst[i] = NamEst[i](1, NamEst[i].Length()-2);
    
    // The corresponding observable
    IRet = fscanf(fpinp,"%s",&buf[0]);     
    NamBuf = &buf[0];
    NamBuf = NamBuf(1, NamBuf.Length()-2);
    IObs = -1;
    if(NumObs == 1){
      IWhichObs[i] = 0;
      IObs = 0;
    }else{
      for(Int_t n = 0; n < NumObs; n++){	
	if(NamObs[n] == &buf[0])IWhichObs[i] = n;
	IObs = n;
      }
    }    
    printf("... ForttoBlue(): NamEst: %2i=%s %2i=%s \n",i,NamEst[i].Data(),IObs,NamBuf.Data());

    // Loop over values
    IRet = fscanf(fpinp,"%f",&valdum);
    XEst[ind] = static_cast<Double_t>(valdum);
    ind = ind + 1;
    for(Int_t k = 0; k < NumUnc; k++){	
      IRet = fscanf(fpinp,"%f",&valdum);
      XEst[ind] = static_cast<Double_t>(valdum);
      ind = ind + 1;
    }
  }

  // The correlation (Rho) matricees
  TMatrixD* Rho[NumUnc];
  for(Int_t k = 0; k < NumUnc; k++){  
    Rho[k] = new TMatrixD(NumEst,NumEst);
    Float_t FWhat;
    for(Int_t i = 0; i < NumEst; i++){
      for(Int_t j = 0; j < NumEst; j++){
	IRet = fscanf(fpinp,"%f ",&FWhat);
	Rho[k]->operator()(i,j) = static_cast<Double_t>(FWhat);
      }
      if(i == 0)ret = fgets(buf,200,fpinp);
    }
  }

  //
  // Try to see which martices are identical
  //
  TMatrixD* NCor = new TMatrixD(NumEst,NumEst);
  TMatrixD* FCor = new TMatrixD(NumEst,NumEst);
  TMatrixD* Null = new TMatrixD(NumEst,NumEst);
  TMatrixD* Comp = new TMatrixD(NumEst,NumEst);
  NCor->UnitMatrix();
  for(Int_t i = 0; i < NumEst; i++){
    for(Int_t j = 0; j < NumEst; j++){
	FCor->operator()(i,j) = 1;
	Null->operator()(i,j) = 0;
    }
  }
  // printf("... ForttoBlue: NCor \n"); NCor->Print();
  // printf("... ForttoBlue: FCor \n"); FCor->Print();
  // printf("... ForttoBlue: Null \n"); Null->Print();

  // Loop over matricees compare to FCor, NCor or previously found matrices
  Int_t IsWhichRho[NumUnc];
  Int_t ITest = 0;
  for(Int_t k = 0; k < NumUnc; k++){  
    IsWhichRho[k] = -5;
    ITest = 0;
    // Fully-correlated
    Comp->operator=(*Null);
    Comp->operator+=(*FCor);
    Comp->operator-=(*Rho[k]);
    for(Int_t i = 0; i < NumEst; i++){
      for(Int_t j = i; j < NumEst; j++){
	if(Comp->operator()(j,i) != 0)ITest = 1;
      }     
    }
    if(ITest == 0){
      printf("... ForttoBlue: %2i = %s is Fully-correlated \n",k,NamUnc[k].Data());
      IsWhichRho[k] = -1;
    }else{
      // Un-correlated
      Comp->operator=(*Null);
      Comp->operator+=(*NCor);
      Comp->operator-=(*Rho[k]);
      ITest = 0;
      for(Int_t i = 0; i < NumEst; i++){
	for(Int_t j = i; j < NumEst; j++){
	  if(Comp->operator()(j,i) != 0)ITest = 1;
	}     
      }
      if(ITest == 0){
	printf("... ForttoBlue: %2i = %s is Un-correlated \n",k,NamUnc[k].Data());
	IsWhichRho[k] = -2;
      }else{
	for(Int_t kk = 0; kk < k; kk++){  
	  if(IsWhichRho[k] != -2 && IsWhichRho[k] != -1 && ITest > -1){
	    Comp->operator=(*Null);
	    Comp->operator+=(*Rho[kk]);
	    Comp->operator-=(*Rho[k]);
	    for(Int_t i = 0; i < NumEst; i++){
	      for(Int_t j = i; j < NumEst; j++){
		if(Comp->operator()(j,i) != 0)ITest = 1;
	      }     
	    }
	    if(ITest == 0){
	      printf("... ForttoBlue: %2i = %s equals Matrix %2i \n",k,NamUnc[k].Data(),kk);
	      IsWhichRho[k] = kk;
	      ITest = -2;
	    }else{
	      ITest = 0;
	    }
	  }
	}
	if(IsWhichRho[k] == -5){
	  IsWhichRho[k] = k;
	  printf("... ForttoBlue: %2i = %s is a new Matrix \n",k,NamUnc[k].Data());
	}
      }
    }
  }
  Comp->Delete();
  
  //
  // ------------ We are done. That is what we got 
  //
  printf("... ForttoBlue: That is what we got: \n");
  printf("... ForttoBlue: IWhichObs =");
  for(Int_t i = 0; i < NumEst; i++){printf(" %2i",IWhichObs[i]);};
  printf("\n");

  printf("... ForttoBlue: NamEst =");
  for(Int_t i = 0; i < NumEst; i++){printf(" %s",NamEst[i].Data());};
  printf("\n");

  printf("... ForttoBlue: NamUnc =");
  for(Int_t k = 0; k < NumUnc; k++){printf(" %s",NamUnc[k].Data());};
  printf("\n");

  printf("... ForttoBlue: NamObs =");
  for(Int_t n = 0; n < NumObs; n++){printf(" %s",NamObs[n].Data());};
  printf("\n");

  ind = 0;
  for(Int_t i = 0; i < NumEst; i++){
    printf("%2i %5.2f ",i,XEst[ind]);
    ind = ind + 1;
    for(Int_t k = 0; k < NumUnc; k++){
      printf(" %5.2f ",XEst[ind]);
      ind = ind + 1;
    }
    printf("\n");
  }

  for(Int_t k = 0; k < NumUnc; k++){  
    printf("... ForttoBlue: Correlation matrix for %s\n",NamUnc[k].Data());
    for(Int_t i = 0; i < NumEst; i++){
      for(Int_t j = 0; j < NumEst; j++){
	printf(" %5.2f ",Rho[k]->operator()(i,j));
      }
      printf("\n");
    }        
  }

  printf("... ForttoBlue: IsWhichRho =");
  for(Int_t k = 0; k < NumUnc; k++){printf(" %2i",IsWhichRho[k]);};
  printf("\n");


  // Close input file
  fclose(fpinp);

  // ---------------------------------------------------------------

  // Open output file
  std::ofstream ofs (FilOut, std::ofstream::out);
  
  if(ofs.good() == false){
    printf("... ForttoBlue(): File %s cannot be opened \n", FilOut.Data());
    return;
  }else{
    printf("... ForttoBlue(): Output file %s has been openend \n", FilOut.Data());
  }

  // Declare variables
  char c[100];
  TString Format = "To be filled later";
  TString Fordum = "To be filled later";

  // Write header
  IRet = sprintf(c,"#include \"Riostream.h\" \n"); ofs<<c;
  IRet = sprintf(c,"#include \"Blue.h\" \n"); ofs<<c;
  IRet = sprintf(c," \n"); ofs<<c;
  IRet = sprintf(c,"void %s (Int_t Flag = 0){ \n",FilCod.Data()); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  // Declare variables
  IRet = sprintf(c,"\t // The number of estimates, uncertainties and observables\n"); ofs<<c;
  IRet = sprintf(c,"\t static const Int_t NumEst = %2i;\n",NumEst); ofs<<c;
  IRet = sprintf(c,"\t static const Int_t NumUnc = %2i;\n",NumUnc); ofs<<c;
  IRet = sprintf(c,"\t static const Int_t NumObs = %2i;\n",NumObs); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  Int_t ISumCo = 0;
  
  // Declare and fill names
  // Estimates
  ISumCo = 0;
  IRet = sprintf(c,"\t TString NamEst[NumEst] = {"); ofs<<c;
  IRet = sprintf(c,"\"%s\"",NamEst[0].Data()); ofs<<c;
  for(Int_t i = 1; i < NumEst; i++){
    IRet = sprintf(c,", \"%s\"",NamEst[i].Data()); ofs<<c;
    ISumCo = ISumCo + 1;
    if(ISumCo == 8){ISumCo = 0;IRet = sprintf(c,"\n\t\t\t"); ofs<<c;};
  }
  IRet = sprintf(c,"};\n"); ofs<<c;

  // Uncertainties
  ISumCo = 0;
  IRet = sprintf(c,"\t TString NamUnc[NumUnc] = {"); ofs<<c;
  IRet = sprintf(c,"\"%s\"",NamUnc[0].Data()); ofs<<c;
  for(Int_t k = 1; k < NumUnc; k++){
    IRet = sprintf(c,", \"%s\"",NamUnc[k].Data()); ofs<<c;
    ISumCo = ISumCo + 1;
    if(ISumCo == 8){ISumCo = 0;IRet = sprintf(c,"\n\t\t\t"); ofs<<c;};
  }
  IRet = sprintf(c,"};\n"); ofs<<c;

  // Observables
  ISumCo = 0;
  IRet = sprintf(c,"\t TString NamObs[NumObs] = {"); ofs<<c;
  IRet = sprintf(c,"\"%s\"",NamObs[0].Data()); ofs<<c;
  for(Int_t n = 1; n < NumObs; n++){
    IRet = sprintf(c,", \"%s\"",NamObs[n].Data()); ofs<<c;
    ISumCo = ISumCo + 1;
    if(ISumCo == 8){ISumCo = 0;IRet = sprintf(c,"\n\t\t\t"); ofs<<c;};
  }
  IRet = sprintf(c,"};\n\n"); ofs<<c;

  // Set IWhichObs
  IRet = sprintf(c,"\t // Index for which estimates determines which observable \n"); ofs<<c;
  IRet = sprintf(c,"\t Int_t IWhichObs[NumEst] = {"); ofs<<c;
  Format = "%2i,";
  for(Int_t i = 0; i < NumEst; i++){
    if(i==NumEst-1)Format = "%2i};\n";
    IRet = sprintf(c,Format,IWhichObs[i]); ofs<<c;
  }
  IRet = sprintf(c,"\n"); ofs<<c;

  // Set the Flag 
  IRet = sprintf(c,"\t // Preset according to Flag \n"); ofs<<c;
  IRet = sprintf(c,"\t if(Flag == 0){ \n"); ofs<<c;
  Format = "\t\t printf(\"... %s: I use Flag = %%2i \\n\",Flag);\n";
  IRet = sprintf(c,Format,FilCod.Data()); ofs<<c;
  IRet = sprintf(c,"\t }else if(Flag == 1){\n"); ofs<<c;
  IRet = sprintf(c,Format,FilCod.Data()); ofs<<c;
  IRet = sprintf(c,"\t }else{\n"); ofs<<c;
  Format = "\t\t printf(\"... %s: Not implemented Flag = %%2i \\n\",Flag);\n";
  IRet = sprintf(c,Format,FilCod.Data()); ofs<<c;
  IRet = sprintf(c,"\t\t return;\n"); ofs<<c;
  IRet = sprintf(c,"\t}\n"); ofs<<c;

  // The estimates 
  IRet = sprintf(c,"\n"); ofs<<c;
  IRet = sprintf(c,"\t //The estimates \n"); ofs<<c;
  IRet = sprintf(c,"\t static const Int_t LenXEst = NumEst * (NumUnc+1);\n"); ofs<<c;
  IRet = sprintf(c,"\t Double_t XEst[LenXEst] = {\n"); ofs<<c;
  ind = 0;
  for(Int_t i = 0; i < NumEst; i++){
    Format = "\t" + ForVal + ",";
    IRet = sprintf(c,Format,XEst[ind]); ofs<<c;
    Format = " " + ForUnc + ",";
    ind = ind + 1;
    for(Int_t k = 0; k < NumUnc; k++){
      if(i == NumEst-1 && k == NumUnc-1)Format = " " + ForUnc;
      IRet = sprintf(c,Format,XEst[ind]); ofs<<c;
      ind = ind + 1;
    }
    IRet = sprintf(c,"\n"); ofs<<c;
  }
  IRet = sprintf(c,"\t}; \n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  // The correlation matrices
  Int_t IFound = 0;
  Int_t IntMat = 0;
  IRet = sprintf(c,"\t static const Int_t LenCor = NumEst * NumEst; \n"); ofs<<c;
  for(Int_t k = 0; k < NumUnc; k++){
    if(IsWhichRho[k] >= 0){
      IFound = 0;
      for(Int_t kk = 0; kk < k; kk++){
	if(IsWhichRho[k] == IsWhichRho[kk])IFound = 1;
      }

      if(IFound == 0){
	// Write header
	if(k < 10){
	  IRet = sprintf(c,"\t Double_t Cor0%1i[LenCor] = {",k); ofs<<c;
	}else{
	  IRet = sprintf(c,"\t Double_t Cor%2i[LenCor] = {",k); ofs<<c;
	}
	// Write matrix
	IntMat = IntMat + 1;	
	for(Int_t i = 0; i < NumEst; i++){
	  Format = "\n\t\t" + ForRho + ",";
	  IRet = sprintf(c,Format,Rho[k]->operator()(i,0)); ofs<<c;
	  Format = " " + ForRho + ",";
	  for(Int_t j = 1; j < NumEst; j++){
	    if(j == NumEst-1 && i == NumEst-1)Format = " " + ForRho;
	    IRet = sprintf(c,Format,Rho[k]->operator()(i,j)); ofs<<c;
	  }
	}
	IRet = sprintf(c,"\n\t\t};\n"); ofs<<c;
      }
    }
  }
  if(IntMat == 0){IRet = sprintf(c,"\t if(LenCor==-1)return;\n"); ofs<<c;};


  // Construct the object 
  IRet = sprintf(c,"\n"); ofs<<c;
  IRet = sprintf(c,"\t // Construct Object \n"); ofs<<c;
  IRet = sprintf(c,"\t Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]); \n"); ofs<<c;
  IRet = sprintf(c,"\t myBlue->PrintStatus(); \n\n"); ofs<<c;

  // Fill names
  IRet = sprintf(c,"\t // Fill names \n"); ofs<<c;
  IRet = sprintf(c,"\t myBlue->FillNamEst(&NamEst[0]); \n"); ofs<<c;
  IRet = sprintf(c,"\t myBlue->FillNamUnc(&NamUnc[0]); \n"); ofs<<c;
  IRet = sprintf(c,"\t myBlue->FillNamObs(&NamObs[0]); \n"); ofs<<c;

  // Fill estimates
  IRet = sprintf(c,"\n"); ofs<<c;
  IRet = sprintf(c,"\t // Fill estimates \n"); ofs<<c;
  IRet = sprintf(c,"\t Int_t ind = 0; \n"); ofs<<c;
  IRet = sprintf(c,"\t for(Int_t i = 0; i<NumEst; i++){ \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->FillEst(i,&XEst[ind]);\n"); ofs<<c;
  IRet = sprintf(c,"\t\t ind = ind + NumUnc + 1; \n"); ofs<<c;
  IRet = sprintf(c,"\t}\n"); ofs<<c;

  // Fill correlations
  IRet = sprintf(c,"\n"); ofs<<c;
  IRet = sprintf(c,"\t // Fill correlations \n"); ofs<<c;
  IRet = sprintf(c,"\t for(Int_t k = 0; k<NumUnc; k++){\n"); ofs<<c;

  // First the uncorrelated
  ISumCo = 0;
  IFound = 0;
  for(Int_t k = 0; k < NumUnc; k++){
    if(IsWhichRho[k] == -2){
      if(IFound == 0){
	IRet = sprintf(c,"\t\t if(k == %2i",k); ofs<<c;
	IFound = 1;
      }else{
	ISumCo = ISumCo + 1;
	IRet = sprintf(c,"|| k==%2i ",k); ofs<<c;
	if(ISumCo == 8){
	  ISumCo = 0;
	  IRet = sprintf(c,"\n\t\t"); ofs<<c;
	}
      }
    }
  }
  if(IFound == 1){
    IRet = sprintf(c,"){\n"); ofs<<c;
    IRet = sprintf(c,"\t\t\t myBlue->FillCor(k,0.0); \n"); ofs<<c;
  }

  // Now the fully correlated
  IFound = 0;
  ISumCo = 0;
  for(Int_t k = 0; k < NumUnc; k++){
    if(IsWhichRho[k] == -1){
      if(IFound == 0){
	IRet = sprintf(c,"\t\t }else if(k == %2i",k); ofs<<c;
	IFound = 1;
      }else{
	ISumCo = ISumCo + 1;
	IRet = sprintf(c,"|| k==%2i ",k); ofs<<c;
	if(ISumCo == 8){
	  ISumCo = 0;
	  IRet = sprintf(c,"\n\t\t"); ofs<<c;
	}
      }
    }
  }
  if(IFound == 1){
    IRet = sprintf(c,"){\n"); ofs<<c;
    IRet = sprintf(c,"\t\t\t myBlue->FillCor(k,1.0); \n"); ofs<<c;
  }

  // Now the individual matrices
  for(Int_t k = 0; k < NumUnc; k++){
    if(IsWhichRho[k] > -1){
      IRet = sprintf(c,"\t\t }else if(k == %2i){\n",k); ofs<<c;
      if(IsWhichRho[k] < 10 ){
	IRet = sprintf(c,"\t\t\t myBlue->FillCor(k,&Cor0%1i[0]); \n",IsWhichRho[k]); ofs<<c;
      }else{
	IRet = sprintf(c,"\t\t\t myBlue->FillCor(k,&Cor%2i[0]); \n",IsWhichRho[k]); ofs<<c;
      }
    }
  }
  IRet = sprintf(c,"\t\t}\n"); ofs<<c;
  IRet = sprintf(c,"\t}\n"); ofs<<c;

  // Solve depending of Flag
  // The default
  IRet = sprintf(c,"\t // Fix input, solve according to Flag and finally delete \n"); ofs<<c;
  IRet = sprintf(c,"\t if(Flag == 0){\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->FixInp(); \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->PrintEst();\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->Solve();  \n"); ofs<<c;
  IRet = sprintf(c,"\t\t printf(\"... B_EPJC_72_2046:"); ofs<<c;
  IRet = sprintf(c,"  mtop(1d + 2d full combination)\\n\");\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->PrintResult();\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->PrintCompatEst();\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->PrintPull();\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->LatexResult(\"%s\");\n",FilCod.Data()); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->DisplayResult(0,\"%s\");\n",FilCod.Data()); ofs<<c;
  // the alternative
  IRet = sprintf(c,"\t }else if(Flag == 1){\n"); ofs<<c;
  IRet = sprintf(c,"\t\t // Set more detailed print level \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->SetPrintLevel(1); \n\n"); ofs<<c;
  IRet = sprintf(c,"\t\t // Change, fix and solve \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->SetInActiveEst(2); \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->SetInActiveEst(3); \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->FixInp(); \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->Solve();  \n"); ofs<<c;
  IRet = sprintf(c,"\t\t printf(\"... B_EPJC_72_2046:"); ofs<<c;
  IRet = sprintf(c," mtop(1d independent combination)\\n\");\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->PrintResult();\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->PrintPull();\n\n"); ofs<<c;
  IRet = sprintf(c,"\t\t // Reset, change, fix and solve \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->ResetInp(); \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->SetInActiveEst(0); \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->SetInActiveEst(1); \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->FixInp(); \n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->Solve();  \n"); ofs<<c;
  IRet = sprintf(c,"\t\t printf(\"... B_EPJC_72_2046:"); ofs<<c;
  IRet = sprintf(c," mtop(2d independent combination)\\n\");\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->PrintResult();\n"); ofs<<c;
  IRet = sprintf(c,"\t\t myBlue->PrintPull();\n"); ofs<<c;
  IRet = sprintf(c,"\t }\n"); ofs<<c;
  IRet = sprintf(c,"\t delete myBlue;\n"); ofs<<c;

  // Write footer
  IRet = sprintf(c,"} \n"); ofs<<c;

  // Close output file
  ofs.close();

  // ---------------------------------------------------------------

  // Open steering file
  std::ofstream ofs2 (FilSte, std::ofstream::out);
  
  IRet = sprintf(c,"gSystem->Load(\"libBlue.so\"); \n"); ofs2<<c;
  IRet = sprintf(c,"\n"); ofs2<<c;

  Format = ".L B_" + FilNam + ".cxx++\n";
  IRet = sprintf(c,Format,FilNam.Data()); ofs2<<c;

  Format = "B_" + FilNam + "(0)\n";
  IRet = sprintf(c,Format,FilNam.Data()); ofs2<<c;

  Format = "B_" + FilNam + "(1)\n";
  IRet = sprintf(c,Format,FilNam.Data()); ofs2<<c;

  IRet = sprintf(c,"\n"); ofs2<<c;
  Format = ".L B_" + FilNam + "_DisRes_Obs_0.cxx++\n";
  IRet = sprintf(c,Format,FilNam.Data()); ofs2<<c;

  Format = "B_" + FilNam + "_DisRes_Obs_0()\n";
  IRet = sprintf(c,Format,FilNam.Data()); ofs2<<c;

  IRet = sprintf(c,"\n"); ofs2<<c;
  IRet = sprintf(c,".q \n"); ofs2<<c;

  // Close steering file
  ofs2.close();

  return;
};
 
