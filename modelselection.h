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
    std::string formatStringResult(int value);
    std::string getED50BestModel(std::string model);
    std::string getAUCBestModel(std::string model);

    void FitNoise();
    void FitExponential(const char *mStarts);
    void FitHyperbolic(const char *mStarts);
    void FitQuasiHyperbolic(const char *mStarts);
    void FitMyerson(const char *mStarts);
    void FitRachlin(const char *mStarts);

    void PrepareProbabilities();

    void InitializeDefaults();

    double NoiseBIC;

    std::vector<std::pair<std::string, double> > mBicList;
    std::vector<std::pair<std::string, double> > mProbList;

    //QList<QPair<QString, double>> mBicList;
    //QList<QPair<QString, double>> mProbList;

    /** AICs
     *
     */
    double aicNoise;
    double aicHyperbolic;
    double aicExponential;
    double aicQuasiHyperbolic;
    double aicMyerson;
    double aicRachlin;

    /** BICs
     *
     */
    double bicNoise;
    double bicHyperbolic;
    double bicExponential;
    double bicQuasiHyperbolic;
    double bicMyerson;
    double bicRachlin;

    /** Bayes Factors
      *
      */
    double bfNoise;
    double bfHyperbolic;
    double bfExponential;
    double bfQuasiHyperbolic;
    double bfMyerson;
    double bfRachlin;

    /** Probs
      *
      */
    double probsNoise;
    double probsHyperbolic;
    double probsExponential;
    double probsQuasiHyperbolic;
    double probsMyerson;
    double probsRachlin;

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

    double AVE;
    double SSR;
    double sumErr;
    double N;
    double S2;

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

    /** BIC Stuff
      *
      */
    double L;
    double DF;
    double PROJ;
    double holder;

    double sumBayesFactors;
};

#endif // MODELSELECTION_H
