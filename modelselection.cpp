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

#include "modelselection.h"

#include <iostream>
#include "interpolation.h"

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include "optimization.h"
#include "integration.h"

using namespace alglib;
using namespace std;

void exponential_discounting(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = exp(-exp(c[0])*x[0]);
}

void exponential_integration(double x, double xminusa, double bminusx, double &y, void *ptr)
{
    //QList<double> *param = (QList<double> *) ptr;
    double *param = (double *) ptr;
    //double k = param->at(0);
    double k = param[0];
    y = exp(-exp(k)*x);
}

void hyperbolic(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = pow((1+exp(c[0])*x[0]), -1);
}

void hyperbolic_integration(double x, double xminusa, double bminusx, double &y, void *ptr)
{
    //QList<double> *param = (QList<double> *) ptr;
    double *param = (double *) ptr;
    //double k = param->at(0);
    double k = param[0];
    y = pow((1+exp(k)*x), -1);
}

void hyperboloid_myerson(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = pow((1+exp(c[0])*x[0]), -c[1]);
}

void hyperboloid_myerson_integration(double x, double xminusa, double bminusx, double &y, void *ptr)
{
    //QList<double> *param = (QList<double> *) ptr;
    double *param = (double *) ptr;
    //double k = param->at(0);
    //double s = param->at(1);
    double k = param[0];
    double s = param[1];

    y = pow((1+exp(k)*x), -s);
}

void hyperboloid_rachlin(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = pow((1+exp(c[0])*pow(x[0], c[1])), -1);
}

void hyperboloid_rachlin_integration(double x, double xminusa, double bminusx, double &y, void *ptr)
{
    //QList<double> *param = (QList<double> *) ptr;
    double *param = (double *) ptr;
    //double k = param->at(0);
    //double s = param->at(1);
    double k = param[0];
    double s = param[1];
    y = pow((1+exp(k)*pow(x, s)), -1);
}

void quasi_hyperboloid(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = c[0] * pow(c[1], x[0]);
}

void quasi_hyperboloid_integration(double x, double xminusa, double bminusx, double &y, void *ptr)
{
    //QList<double> *param = (QList<double> *) ptr;
    double *param = (double *) ptr;
    //double b = param->at(0);
    //double d = param->at(1);
    double b = param[0];
    double d = param[1];
    y = b * pow(d, x);
}

void ModelSelection::SetX(const char *mString)
{
    x = mString;
}

void ModelSelection::SetY(const char *mString)
{
    y = mString;
}

void ModelSelection::SetStarts(const char *mString)
{
    c = mString;
}

void ModelSelection::SetLowerUpperBounds(const char *mUpperString, const char *mLowerString)
{
    bndu = mUpperString;
    bndl = mLowerString;
}

real_1d_array ModelSelection::GetParams()
{
    return c;
}

lsfitstate ModelSelection::GetState()
{
    return state;
}

ae_int_t ModelSelection::GetInfo()
{
    return info;
}

lsfitreport ModelSelection::GetReport()
{
    return rep;
}

double ModelSelection::GetNoiseMean()
{
    return AVE;
}

/** Calculate Noise BIC
 * @brief ModelSelection::FitNoise
 */
void ModelSelection::FitNoise()
{
    N = y.length();
    AVE = 0;

    for (int i = 0; i < N; i++)
    {
        AVE += (double) y[i];
    }

    AVE = AVE / N;
    SSR = 0;

    for (int i = 0; i < N; i++)
    {
        SSR += pow(((double) y[i] - AVE), 2);
    }

    S2 = SSR / N;

    L = pow((1.0 / sqrt(2 * M_PI * S2)), N) * exp(-SSR / (S2 * 2.0));

    DF = 2;

    aicNoise = (-2 * log(L)) + (2 * DF);
    bicNoise = -2 * log(L) + log(N) * DF;
}

/** Exponentials
 * @brief
 */
void ModelSelection::FitExponential(const char *mStarts)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponential_discounting);

    lsfitresults(state, info, c, rep);

    N = y.length();
    SSR = 0;

    for (int i = 0; i < N; i++)
    {
        holder = (exp(-exp( (double) c[0])* (double) x[i][0]));

        SSR += pow(((double) y[i] - holder), 2);
    }

    S2 = SSR / N;

    L = pow((1.0 / sqrt(2 * M_PI * S2)), N) * exp(-SSR / (S2 * 2.0));

    DF = 2;

    aicExponential = (-2 * log(L)) + (2 * DF);
    bicExponential = -2 * log(L) + log(N) * DF;
    fitExponentialK = (double) c[0];
}

/** Hyperbolics
 * @brief
 */
void ModelSelection::FitHyperbolic(const char *mStarts)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, hyperbolic);

    lsfitresults(state, info, c, rep);

    N = y.length();
    SSR = 0;

    for (int i = 0; i < N; i++)
    {
        holder = pow((1+exp( (double) c[0])* (double) x[i][0]), -1);
        SSR += pow(((double) y[i] - holder), 2);
    }

    S2 = SSR / N;

    L = pow((1.0 / sqrt(2 * M_PI * S2)), N) * exp(-SSR / (S2 * 2.0));

    DF = 2;

    aicHyperbolic = (-2 * log(L)) + (2 * DF);
    bicHyperbolic = -2 * log(L) + log(N) * DF;
    fitHyperbolicK = (double) c[0];
}

/** Hyperboloid Myerson
 * @brief
 */
void ModelSelection::FitMyerson(const char *mStarts)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, hyperboloid_myerson);

    lsfitresults(state, info, c, rep);

    N = y.length();
    SSR = 0;

    for (int i = 0; i < N; i++)
    {
        holder = pow((1+exp( (double) c[0])* (double) x[i][0]),  (double) -c[1]);
        SSR += pow(((double) y[i] - holder), 2);
    }

    S2 = SSR / N;

    L = pow((1.0 / sqrt(2 * M_PI * S2)), N) * exp(-SSR / (S2 * 2.0));

    DF = 3;

    aicMyerson = (-2 * log(L)) + (2 * DF);
    bicMyerson = -2 * log(L) + log(N) * DF;
    fitMyersonK = (double) c[0];
    fitMyersonS = (double) c[1];
}

/** Beta delta
  *  @brief
  */
void ModelSelection::FitQuasiHyperbolic(const char *mStarts)
{
    SetStarts(mStarts);
    SetLowerUpperBounds("[1.0, 1.0]", "[0.0, 0.0]");

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);
    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, quasi_hyperboloid);

    lsfitresults(state, info, c, rep);

    N = y.length();
    SSR = 0;

    for (int i = 0; i < N; i++)
    {
        holder = (double) c[0] * pow( (double) c[1], (double) x[i][0]);
        SSR += pow(((double) y[i] - holder), 2);
    }

    S2 = SSR / N;

    L = pow((1.0 / sqrt(2 * M_PI * S2)), N) * exp(-SSR / (S2 * 2.0));

    DF = 3;

    aicQuasiHyperbolic = (-2 * log(L)) + (2 * DF);
    bicQuasiHyperbolic = -2 * log(L) + log(N) * DF;
    fitQuasiHyperbolicBeta = (double) c[0];
    fitQuasiHyperbolicDelta = (double) c[1];
}

/** Hyperboloid Rachlin
  *  @brief
  */
void ModelSelection::FitRachlin(const char *mStarts)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, hyperboloid_rachlin);

    lsfitresults(state, info, c, rep);

    N = y.length();
    SSR = 0;

    for (int i = 0; i < N; i++)
    {
        holder = pow((1+exp( (double) c[0])*pow( (double) x[i][0], (double) c[1])), -1);
        SSR += pow(((double) y[i] - holder), 2);
    }

    S2 = SSR / N;

    L = pow((1.0 / sqrt(2 * M_PI * S2)), N) * exp(-SSR / (S2 * 2.0));

    DF = 3;

    aicRachlin = (-2 * log(L)) + (2 * DF);
    bicRachlin = -2 * log(L) + log(N) * DF;
    fitRachlinK = (double) c[0];
    fitRachlinS = (double) c[1];
}

double ScaleFactor(double modelBic, double noiseBic)
{
    return exp(-0.5 * (modelBic - noiseBic));
}

std::string ModelSelection::formatStringResult(int value)
{
    if (value == -7)
    {
        return std::string("gradient verification failed");
    }
    else if (value == 2)
    {
        return std::string("Success: relative step is no more than EpsX");
    }
    else if (value == 5)
    {
        return std::string("MaxIts steps was taken");
    }
    else if (value == 7)
    {
        return std::string("stopping conditions are too stringent, further improvement is impossible");
    }
    else
    {
        return std::string("No notes supplied");
    }
}

std::string ModelSelection::getED50BestModel(std::string model)
{
    std::ostringstream strs;

    if (model.find("Exponential") != std::string::npos)
    {
        double result = log(1/(exp(fitExponentialK)));
        strs << result;

        return strs.str();
    }
    else if (model.find("Hyperbolic") != std::string::npos)
    {
        double result = log(1/(exp(fitHyperbolicK)));
        strs << result;

        return strs.str();
    }
    else if (model.find("Beta") != std::string::npos)
    {
        double result = log(log((1/(2*fitQuasiHyperbolicBeta)))/log(fitQuasiHyperbolicDelta));
        strs << result;

        return strs.str();
    }
    else if (model.find("Myerson") != std::string::npos)
    {
        double result = log((pow(2, (1/fitMyersonS))-1)/exp(fitMyersonK));
        strs << result;

        return strs.str();
    }
    else if (model.find("Rachlin") != std::string::npos)
    {
        double result = log( pow( (1/ (exp(fitRachlinK))), (1/fitRachlinS)));
        strs << result;

        return strs.str();
    }
    else
    {
        strs << "NA";

        return strs.str();
    }
}

std::string ModelSelection::getAUCBestModel(std::string model)
{
    std::ostringstream strs;

    double result = -1;
    double a = x[0][0];
    double b = x[x.rows() - 1][0];

    double mParams[2];

    //QList<double> mParams;
    autogkstate s;
    double v;
    autogkreport rep;

    if (model.find("Exponential") != std::string::npos)
    {
        mParams[0] = fitExponentialK;
        //mParams << fitExponentialK;

        autogksmooth(a, b, s);
        alglib::autogkintegrate(s, exponential_integration, &mParams);
        autogkresults(s, v, rep);

        result = double(v) / (b - a);

        strs << result;

        return strs.str();
    }
    else if (model.find("Hyperbolic") != std::string::npos)
    {
        mParams[0] = fitHyperbolicK;
        //mParams << fitHyperbolicK;

        autogksmooth(a, b, s);
        alglib::autogkintegrate(s, hyperbolic_integration, &mParams);
        autogkresults(s, v, rep);

        result = double(v) / (b - a);

        strs << result;

        return strs.str();
    }
    else if (model.find("Beta") != std::string::npos)
    {
        mParams[1] = fitQuasiHyperbolicBeta;
        mParams[2] = fitQuasiHyperbolicDelta;

        //mParams << fitQuasiHyperbolicBeta;
        //mParams << fitQuasiHyperbolicDelta;

        autogksmooth(a, b, s);
        alglib::autogkintegrate(s, quasi_hyperboloid_integration, &mParams);
        autogkresults(s, v, rep);

        result = double(v) / (b - a);

        strs << result;

        return strs.str();
    }
    else if (model.find("Myerson") != std::string::npos)
    {
        mParams[1] = fitMyersonK;
        mParams[2] = fitMyersonS;

        //mParams << fitMyersonK;
        //mParams << fitMyersonS;

        autogksmooth(a, b, s);
        alglib::autogkintegrate(s, hyperboloid_myerson_integration, &mParams);
        autogkresults(s, v, rep);

        result = double(v) / (b - a);

        strs << result;

        return strs.str();
    }
    else if (model.find("Rachlin") != std::string::npos)
    {
        mParams[1] = fitRachlinK;
        mParams[2] = fitRachlinS;

        //mParams << fitRachlinK;
        //mParams << fitRachlinS;

        autogksmooth(a, b, s);
        alglib::autogkintegrate(s, hyperboloid_rachlin_integration, &mParams);
        autogkresults(s, v, rep);

        result = double(v) / (b - a);

        strs << result;

        return strs.str();
    }
    else
    {
        strs << "NA";

        return strs.str();
    }

    //autogksmooth(a, b, s);
    //alglib::autogkintegrate(s, hyperbolic_integration, &mParams);
    //autogkresults(s, v, rep);
    //return QString::number(double(v));
}

void ModelSelection::InitializeDefaults()
{
	maxits = 100;

    epsx = 0.000001;
    diffstep = 0.00001;

    SSR = 0;
    N = 0;
    S2 = 0;
    L = 0;
    DF = 0;
    PROJ = 0;
    holder = 0;

    sumBayesFactors = 0;

    AVE = 0;
}

void ModelSelection::PrepareProbabilities()
{
    bfNoise = ScaleFactor(NoiseBIC, NoiseBIC);
    bfExponential = -1;
    bfHyperbolic = -1;
    bfQuasiHyperbolic = -1;
    bfMyerson = -1;
    bfRachlin = -1;

    sumBayesFactors = 0;

    /** Set up BF's
      */
    for (int i=0; i<mBicList.size(); i++)
    {
        std::string mModel = mBicList[i].first;

        if (mModel.find("Exponential") != std::string::npos)
        {
            bfExponential = ScaleFactor(mBicList[i].second, NoiseBIC);
            sumBayesFactors = sumBayesFactors + bfExponential;
        }
        else if (mModel.find("Hyperbolic") != std::string::npos)
        {
            bfHyperbolic = ScaleFactor(mBicList[i].second, NoiseBIC);
            sumBayesFactors = sumBayesFactors + bfHyperbolic;
        }
        else if (mModel.find("Beta") != std::string::npos)
        {
            bfQuasiHyperbolic = ScaleFactor(mBicList[i].second, NoiseBIC);
            sumBayesFactors = sumBayesFactors + bfQuasiHyperbolic;
        }
        else if (mModel.find("Myerson") != std::string::npos)
        {
            bfMyerson = ScaleFactor(mBicList[i].second, NoiseBIC);
            sumBayesFactors = sumBayesFactors + bfMyerson;
        }
        else if (mModel.find("Rachlin") != std::string::npos)
        {
            bfRachlin = ScaleFactor(mBicList[i].second, NoiseBIC);
            sumBayesFactors = sumBayesFactors + bfRachlin;
        }
    }

    probsNoise = bfNoise/sumBayesFactors;
    probsHyperbolic = -1;
    probsExponential = -1;
    probsQuasiHyperbolic = -1;
    probsMyerson = -1;
    probsRachlin = -1;

    mProbList.clear();
    mProbList.push_back(std::pair<std::string, double>("Noise Model", probsNoise));

    for (int i=0; i<mBicList.size(); i++)
    {
        std::string mModel = mBicList[i].first;

        if (mModel.find("Exponential") != std::string::npos)
        {
            probsExponential = bfExponential/sumBayesFactors;
            mProbList.push_back(std::pair<std::string, double>("Exponential Model", probsExponential));
        }
        else if (mModel.find("Hyperbolic") != std::string::npos)
        {
            probsHyperbolic = bfHyperbolic/sumBayesFactors;
            mProbList.push_back(std::pair<std::string, double>("Hyperbolic Model", probsHyperbolic));
        }
        else if (mModel.find("Beta") != std::string::npos)
        {
            probsQuasiHyperbolic = bfQuasiHyperbolic/sumBayesFactors;
            mProbList.push_back(std::pair<std::string, double>("Beta Delta Model", probsQuasiHyperbolic));
        }
        else if (mModel.find("Myerson") != std::string::npos)
        {
            probsMyerson = bfMyerson/sumBayesFactors;
            mProbList.push_back(std::pair<std::string, double>("Myerson Model", probsMyerson));
        }
        else if (mModel.find("Rachlin") != std::string::npos)
        {
            probsRachlin = bfRachlin/sumBayesFactors;
            mProbList.push_back(std::pair<std::string, double>("Rachlin Model", probsRachlin));
        }
    }
}

ModelSelection::ModelSelection()
{

}