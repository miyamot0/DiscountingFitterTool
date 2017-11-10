/**
   Copyright 2017 Shawn Gilroy

   This file is part of Discounting Fitting Tool.

   Discounting Model Selector is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, version 3.

   Discounting Model Selector is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Discounting Model Selector.  If not, see http://www.gnu.org/licenses/.

   The Discounting Fitting Tool is a tool to assist researchers in behavior economics.

   Email: shawn(dot)gilroy(at)temple.edu

   The Discounting Fitting Tool uses ALGLIB to provide access to mathematical methods.

   ====================================================================================

   ALGLIB 3.11.0 (source code generated 2017-05-11)
   Copyright (c) Sergey Bochkanov (ALGLIB project).

   >>> SOURCE LICENSE >>>
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation (www.fsf.org); either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is available at
   http://www.fsf.org/licensing/licenses
   >>> END OF LICENSE >>>

  */

#ifndef MODELSELECTION_H
#define MODELSELECTION_H

#include <iostream>
#include "interpolation.h"

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"

using namespace alglib;

class ModelSelection
{
public:
    ModelSelection();

    void SetX(const char *mString);
    void SetY(const char *mString);
    void SetStarts(const char *mString);
    void SetLowerUpperBounds(const char *mUpperString, const char *mLowerString);

    real_1d_array GetParams();
    lsfitstate GetState();
    ae_int_t GetInfo();
    lsfitreport GetReport();

    double GetNoiseMean();
    double getED50ep();
    double getED50crdi();
    double getED50rodriguez();

    double rodriguez_logue_plotting(double k, double s, double x);
    double ebert_prelec_plotting(double k, double s, double x);
    double bleichrodt_plotting(double k, double s, double beta, double x);

    std::string formatStringResult(int value);
    std::string getED50BestModel();
    std::string getAUCBestModel();
    std::string getLog10AUCBestModel();

    double getErrorExponential(double lnK);
    double getErrorHyperbolic(double lnK);
    double getErrorQuasiHyperbolic(double beta, double delta);
    double getErrorGreenMyerson(double lnK, double s);
    double getErrorRachlin(double lnK, double s);
    double getErrorRodriguezLogue(double lnK, double beta);
    double getErrorEbertPrelec(double lnK, double s);
    double getErrorBleichrodt(double lnK, double s, double beta);

    void FitNoise();
    void FitExponential(const char *mStarts);
    void FitHyperbolic(const char *mStarts);
    void FitQuasiHyperbolic(const char *mStarts);
    void FitMyerson(const char *mStarts);
    void FitRachlin(const char *mStarts);
    void FitRodriguezLogue(const char *mStarts);
    void FitEbertPrelec(const char *mStarts);
    void FitBleichrodt(const char *mStarts);

    void PrepareProbabilities();

    void InitializeDefaults();

    double NoiseBIC;

    std::vector<std::pair<std::string, double> > mBicList;
    std::vector<std::pair<std::string, double> > mProbList;

    /** AICs
     *
     */
    double aicNoise;
    double aicHyperbolic;
    double aicExponential;
    double aicQuasiHyperbolic;
    double aicMyerson;
    double aicRachlin;
    double aicRodriguezLogue;
    double aicEbertPrelec;
    double aicBleichrodt;

    /** BICs
     *
     */
    double bicNoise;
    double bicHyperbolic;
    double bicExponential;
    double bicQuasiHyperbolic;
    double bicMyerson;
    double bicRachlin;
    double bicRodriguezLogue;
    double bicEbertPrelec;
    double bicBleichrodt;

    /** Bayes Factors
      *
      */
    double bfNoise;
    double bfHyperbolic;
    double bfExponential;
    double bfQuasiHyperbolic;
    double bfMyerson;
    double bfRachlin;
    double bfRodriguezLogue;
    double bfEbertPrelec;
    double bfBleichrodt;

    /** Probs
      *
      */
    double probsNoise;
    double probsHyperbolic;
    double probsExponential;
    double probsQuasiHyperbolic;
    double probsMyerson;
    double probsRachlin;
    double probsRodriguezLogue;
    double probsEbertPrelec;
    double probsBleichrodt;

    /** Fits
      *
      */
    double fitHyperbolicK;
    double fitExponentialK;
    double fitQuasiHyperbolicBeta;
    double fitQuasiHyperbolicDelta;
    double fitMyersonK;
    double fitMyersonS;
    double fitRachlinK;
    double fitRachlinS;
    double fitRodriguezLogueK;
    double fitRodriguezLogueBeta;
    double fitEbertPrelecK;
    double fitEbertPrelecS;
    double fitBleichrodtBeta;
    double fitBleichrodtK;
    double fitBleichrodtS;

    double AVE;
    double SSR;
    double sumErr;
    double N;
    double S2;

    /** BIC Stuff
      *
      */
    double L;
    double DF;
    double PROJ;
    double holder;

    double sumBayesFactors;

private:
    real_2d_array x;
    real_1d_array y;
    real_1d_array c;

    real_1d_array bndl;
    real_1d_array bndu;

    ae_int_t maxits;
    ae_int_t info;
    lsfitstate state;
    lsfitreport rep;

    double epsx;
    double diffstep;

    double leastSquaresError;
};

#endif // MODELSELECTION_H
