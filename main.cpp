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

#include <iostream>
#include <sstream>
#include <algorithm>
#include "modelselection.h"

using namespace std;

struct BruteForce {
  double p1;
  double p2;
  double p3;

  double err;

  bool operator < (const BruteForce& r1) const {
      return (err < r1.err);
  }
};

struct BruteForceValues {
    BruteForce oneParamStartingValueArray[100];
    BruteForce twoParamStartingValueArray[1000];
    BruteForce threeParamStartingValueArray[1000];
};

BruteForceValues provisionalValues;

bool BruteSorter(BruteForce const& lhs, BruteForce const& rhs) {
    return lhs.err < rhs.err;
}

double p1Span, p1Step;
double p2Span, p2Step;
double p3Span, p3Step;

int grandLoop;

int main(int argc, char *argv[])
{
	ModelSelection mFitter;
	mFitter.InitializeDefaults();

	//std::string mX = "[[1],[30],[180],[540],[1080],[2160]]";
	//std::string mY = "[1.0,0.9,0.8,0.7,0.6,0.5]";

	//std::string mX = "[[1],[7],[14],[30],[183],[365],[1825],[9125]]";
	//std::string mY = "[1.0,0.95,0.91,0.83,0.81,0.54,0.09,0.04]";

	mFitter.SetX(argv[1]);
    mFitter.SetY(argv[2]);
	//mFitter.SetX(mX.c_str());
	//mFitter.SetY(mY.c_str());

    mFitter.mBicList.clear();

    std::ostringstream startWriter;

    std::ostringstream out;
    out << "{";

    /*
     * Hyperbolic
     *
     * Work on start points here
     *
     */

    startWriter.clear();
    startWriter.str("");

    p1Span = abs(-12) + abs(12); // -12 to 12
    p1Step = p1Span / 100;

    grandLoop = 0;

    for (int kLoop = 0; kLoop < 100; kLoop++)
    {
        provisionalValues.oneParamStartingValueArray[grandLoop].p1 = 12 - ((kLoop + 1) * p1Step);

        grandLoop++;
    }

    for (int mStep = 0; mStep < 100; mStep++)
    {
    	provisionalValues.oneParamStartingValueArray[mStep].err = mFitter.getErrorExponential(provisionalValues.oneParamStartingValueArray[mStep].p1);
    }

    std::sort(provisionalValues.oneParamStartingValueArray, provisionalValues.oneParamStartingValueArray + 100);

    startWriter << "[" << provisionalValues.oneParamStartingValueArray[0].p1 << "]";

    /*
     * Perform fitting with start values
     */

    mFitter.FitHyperbolic(startWriter.str().c_str());

    /*
     * Report back results
     */

    if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
    {
        out << "\"Hyperbolic\":" << mFitter.fitHyperbolicK << ",";
        out << "\"HyperbolicRMS\":" << mFitter.GetReport().rmserror << ",";
        out << "\"Hyperbolicavgerr\":" << mFitter.GetReport().avgerror << ",";
        out << "\"HyperbolicBIC\":" << mFitter.bicHyperbolic << ",";
        out << "\"HyperbolicAIC\":" << mFitter.aicHyperbolic << ",";
        out << "\"HyperbolicCode\":" << (int) mFitter.GetInfo() << ",";

        mFitter.mBicList.push_back(std::pair<std::string, double>("Hyperbolic", mFitter.bicHyperbolic));
    }
    else
    {
        out << "\"Hyperbolic\":" << "\"\"" << ",";
        out << "\"HyperbolicRMS\":" << "\"\"" << ",";
        out << "\"Hyperbolicavgerr\":" << "\"\"" << ",";
        out << "\"HyperbolicBIC\":" << "\"\"" << ",";
        out << "\"HyperbolicAIC\":" << "\"\"" << ",";
        out << "\"HyperbolicCode\":" << (int) mFitter.GetInfo() << ",";
    }

    /*
     * Exponential
     *
     * Work on start points here
     */

    startWriter.clear();
    startWriter.str("");

    p1Span = abs(-12) + abs(12); // -12 to 12
    p1Step = p1Span / 100;

    grandLoop = 0;

    for (int kLoop = 0; kLoop < 100; kLoop++)
    {
        provisionalValues.oneParamStartingValueArray[grandLoop].p1 = 12 - ((kLoop + 1) * p1Step);

        grandLoop++;
    }

    for (int mStep = 0; mStep < 100; mStep++)
    {
    	provisionalValues.oneParamStartingValueArray[mStep].err = mFitter.getErrorHyperbolic(provisionalValues.oneParamStartingValueArray[mStep].p1);
    }

    std::sort(provisionalValues.oneParamStartingValueArray, provisionalValues.oneParamStartingValueArray + 100);

    startWriter << "[" << provisionalValues.oneParamStartingValueArray[0].p1 << "]";

    /*
     * Perform fitting with start values
     */

    mFitter.FitExponential(startWriter.str().c_str());

    /*
     * Report back results
     */

    if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
    {
        out << "\"Exponential\":" << mFitter.fitExponentialK << ",";
        out << "\"ExponentialRMS\":" << mFitter.GetReport().rmserror << ",";
        out << "\"Exponentialavgerr\":" << mFitter.GetReport().avgerror << ",";
        out << "\"ExponentialBIC\":" << mFitter.bicExponential << ",";
        out << "\"ExponentialAIC\":" << mFitter.aicExponential << ",";
        out << "\"ExponentialCode\":" << (int) mFitter.GetInfo() << ",";

        mFitter.mBicList.push_back(std::pair<std::string, double>("Exponential", mFitter.bicExponential));
    }
    else
    {
        out << "\"Exponential\":" << "\"\"" << ",";
        out << "\"ExponentialRMS\":" << "\"\"" << ",";
        out << "\"Exponentialavgerr\":" << "\"\"" << ",";
        out << "\"ExponentialBIC\":" << "\"\"" << ",";
        out << "\"ExponentialAIC\":" << "\"\"" << ",";
        out << "\"ExponentialCode\":" << (int) mFitter.GetInfo() << ",";
    }

    /*
     * Beta Delta
     *
     * Work on start points here
     *
     */

    startWriter.clear();
    startWriter.str("");

    p1Span = 1; // 0 to 1
    p1Step = p1Span / 10;

    p2Span = 1;
    p2Step = p2Span / 100;

    grandLoop = 0;

    for (int bLoop = 0; bLoop < 10; bLoop++)
    {
        for (int dLoop = 0; dLoop < 100; dLoop++)
        {
            provisionalValues.twoParamStartingValueArray[grandLoop].p1 = ((bLoop + 1) * p1Step);
            provisionalValues.twoParamStartingValueArray[grandLoop].p2 = ((dLoop + 1) * p2Step);

            grandLoop++;
        }
    }

    for (int mStep = 0; mStep < 100; mStep++)
    {
    	provisionalValues.oneParamStartingValueArray[mStep].err = mFitter.getErrorQuasiHyperbolic(provisionalValues.oneParamStartingValueArray[mStep].p1, provisionalValues.oneParamStartingValueArray[mStep].p2);
    }

    std::sort(provisionalValues.twoParamStartingValueArray, provisionalValues.twoParamStartingValueArray + 1000);

    /*
     * Perform fitting with start values
     */

    startWriter << "[" << provisionalValues.twoParamStartingValueArray[0].p1 << "," << provisionalValues.twoParamStartingValueArray[0].p2 << "]";

    mFitter.FitQuasiHyperbolic(startWriter.str().c_str());

    /*
     * Report back results
     */

    if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
    {
        out << "\"QuasiHyperbolicBeta\":" << mFitter.fitQuasiHyperbolicBeta << ",";
        out << "\"QuasiHyperbolicDelta\":" << mFitter.fitQuasiHyperbolicDelta << ",";
        out << "\"QuasiHyperbolicRMS\":" << mFitter.GetReport().rmserror << ",";
        out << "\"QuasiHyperbolicavgerr\":" << mFitter.GetReport().avgerror << ",";
        out << "\"QuasiHyperbolicBIC\":" << mFitter.bicQuasiHyperbolic << ",";
        out << "\"QuasiHyperbolicAIC\":" << mFitter.aicQuasiHyperbolic << ",";
        out << "\"QuasiHyperbolicCode\":" << (int) mFitter.GetInfo() << ",";

        mFitter.mBicList.push_back(std::pair<std::string, double>("Beta Delta", mFitter.bicQuasiHyperbolic));
    }
    else
    {
        out << "\"QuasiHyperbolicBeta\":" << "\"\"" << ",";
        out << "\"QuasiHyperbolicDelta\":" << "\"\"" << ",";
        out << "\"QuasiHyperbolicRMS\":" << "\"\"" << ",";
        out << "\"QuasiHyperbolicavgerr\":" << "\"\"" << ",";
        out << "\"QuasiHyperbolicBIC\":" << "\"\"" << ",";
        out << "\"QuasiHyperbolicAIC\":" << "\"\"" << ",";
        out << "\"QuasiHyperbolicCode\":" << (int) mFitter.GetInfo() << ",";
    }

    /*
     * Myerson Green
     *
     * Work on start points here
     */

    startWriter.clear();
    startWriter.str("");

    p1Span = abs(-12) + abs(12); // -12 to 12
    p1Step = p1Span / 10;

    p2Span = 1;
    p2Step = p2Span / 100;

    grandLoop = 0;

    for (int kLoop = 0; kLoop < 10; kLoop++)
    {
        for (int sLoop = 0; sLoop < 100; sLoop++)
        {
            provisionalValues.twoParamStartingValueArray[grandLoop].p1 = 12 - ((kLoop + 1) * p1Step);
            provisionalValues.twoParamStartingValueArray[grandLoop].p2 = ((sLoop + 1) * p2Step);

            grandLoop++;
        }
    }

    for (int mStep = 0; mStep < 100; mStep++)
    {
    	provisionalValues.oneParamStartingValueArray[mStep].err = mFitter.getErrorGreenMyerson(provisionalValues.oneParamStartingValueArray[mStep].p1, provisionalValues.oneParamStartingValueArray[mStep].p2);
    }

    std::sort(provisionalValues.twoParamStartingValueArray, provisionalValues.twoParamStartingValueArray + 1000);

    startWriter << "[" << provisionalValues.twoParamStartingValueArray[0].p1 << "," << provisionalValues.twoParamStartingValueArray[0].p2 << "]";

    /*
     * Perform fitting with start values
     */

    mFitter.FitMyerson(startWriter.str().c_str());

    /*
     * Report back results
     */

    if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
    {
        out << "\"MyersonK\":" << mFitter.fitMyersonK << ",";
        out << "\"MyersonS\":" << mFitter.fitMyersonS << ",";
        out << "\"MyersonRMS\":" << mFitter.GetReport().rmserror << ",";
        out << "\"Myersonavgerr\":" << mFitter.GetReport().avgerror << ",";
        out << "\"MyersonBIC\":" << mFitter.bicMyerson << ",";
        out << "\"MyersonAIC\":" << mFitter.aicMyerson << ",";
        out << "\"MyersonCode\":" << (int) mFitter.GetInfo() << ",";

        mFitter.mBicList.push_back(std::pair<std::string, double>("Myerson", mFitter.bicMyerson));
    }
    else
    {
        out << "\"MyersonK\":" << "\"\"" << ",";
        out << "\"MyersonS\":" << "\"\"" << ",";
        out << "\"MyersonRMS\":" << "\"\"" << ",";
        out << "\"Myersonavgerr\":" << "\"\"" << ",";
        out << "\"MyersonBIC\":" << "\"\"" << ",";
        out << "\"MyersonAIC\":" << "\"\"" << ",";
        out << "\"MyersonCode\":" << (int) mFitter.GetInfo() << ",";
    }

    /*
     * Rachlin
     *
     * Work on start points here
     */

    startWriter.clear();
    startWriter.str("");

    p1Span = abs(-12) + abs(12); // -12 to 12
    p1Step = p1Span / 10;

    p2Span = 1;
    p2Step = p2Span / 100;

    grandLoop = 0;

    for (int kLoop = 0; kLoop < 10; kLoop++)
    {
        for (int sLoop = 0; sLoop < 100; sLoop++)
        {
            provisionalValues.twoParamStartingValueArray[grandLoop].p1 = 12 - ((kLoop + 1) * p1Step);
            provisionalValues.twoParamStartingValueArray[grandLoop].p2 = ((sLoop + 1) * p2Step);

            grandLoop++;
        }
    }

    for (int mStep = 0; mStep < 100; mStep++)
    {
    	provisionalValues.oneParamStartingValueArray[mStep].err = mFitter.getErrorRachlin(provisionalValues.oneParamStartingValueArray[mStep].p1, provisionalValues.oneParamStartingValueArray[mStep].p2);
    }

    std::sort(provisionalValues.twoParamStartingValueArray, provisionalValues.twoParamStartingValueArray + 1000);

    startWriter << "[" << provisionalValues.twoParamStartingValueArray[0].p1 << "," << provisionalValues.twoParamStartingValueArray[0].p2 << "]";

    /*
     * Perform fitting with start values
     */

    mFitter.FitRachlin(startWriter.str().c_str());

    /*
     * Report back results
     */

    if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
    {
        out << "\"RachlinK\":" << mFitter.fitRachlinK << ",";
        out << "\"RachlinS\":" << mFitter.fitRachlinS << ",";
        out << "\"RachlinRMS\":" << mFitter.GetReport().rmserror << ",";
        out << "\"Rachlinavgerr\":" << mFitter.GetReport().avgerror << ",";
        out << "\"RachlinBIC\":" << mFitter.bicRachlin << ",";
        out << "\"RachlinAIC\":" << mFitter.aicRachlin << ",";
        out << "\"RachlinCode\":" << (int) mFitter.GetInfo() << ",";

        mFitter.mBicList.push_back(std::pair<std::string, double>("Rachlin", mFitter.bicRachlin));
    }
    else
    {
        out << "\"RachlinK\":" << "\"\"" << ",";
        out << "\"RachlinS\":" << "\"\"" << ",";
        out << "\"RachlinRMS\":" << "\"\"" << ",";
        out << "\"Rachlinavgerr\":" << "\"\"" << ",";
        out << "\"RachlinBIC\":" << "\"\"" << ",";
        out << "\"RachlinAIC\":" << "\"\"" << ",";
        out << "\"RachlinCode\":" << (int) mFitter.GetInfo() << ",";
    }

    /*
     * RodriguezLogue
     *
     * Work on start points here
     */

    startWriter.clear();
    startWriter.str("");

    p1Span = abs(-12) + abs(12); // -12 to 12
    p1Step = p1Span / 10;

    p2Span = 1;
    p2Step = p2Span / 100;

    grandLoop = 0;

    for (int kLoop = 0; kLoop < 10; kLoop++)
    {
        for (int sLoop = 0; sLoop < 100; sLoop++)
        {
            provisionalValues.twoParamStartingValueArray[grandLoop].p1 = 12 - ((kLoop + 1) * p1Step);
            provisionalValues.twoParamStartingValueArray[grandLoop].p2 = ((sLoop + 1) * p2Step);

            grandLoop++;
        }
    }

    for (int mStep = 0; mStep < 100; mStep++)
    {
    	provisionalValues.oneParamStartingValueArray[mStep].err = mFitter.getErrorEbertPrelec(provisionalValues.oneParamStartingValueArray[mStep].p1,
    			provisionalValues.oneParamStartingValueArray[mStep].p2);
    }

    std::sort(provisionalValues.twoParamStartingValueArray, provisionalValues.twoParamStartingValueArray + 1000);

    startWriter << "[" << provisionalValues.twoParamStartingValueArray[0].p1 << "," << provisionalValues.twoParamStartingValueArray[0].p2 << "]";

    /*
     * Perform fitting with start values
     */

    mFitter.FitRodriguezLogue(startWriter.str().c_str());

    /*
     * Report back results
     */

    if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
    {
        out << "\"RodriguezLogueK\":" << mFitter.fitRodriguezLogueK << ",";
        out << "\"RodriguezLogueBeta\":" << mFitter.fitRodriguezLogueBeta << ",";
        out << "\"RodriguezLogueRMS\":" << mFitter.GetReport().rmserror << ",";
        out << "\"RodriguezLogueavgerr\":" << mFitter.GetReport().avgerror << ",";
        out << "\"RodriguezLogueBIC\":" << mFitter.bicRodriguezLogue << ",";
        out << "\"RodriguezLogueAIC\":" << mFitter.aicRodriguezLogue << ",";
        out << "\"RodriguezLogueCode\":" << (int) mFitter.GetInfo() << ",";

        mFitter.mBicList.push_back(std::pair<std::string, double>("RodriguezLogue", mFitter.bicRodriguezLogue));
    }
    else
    {
        out << "\"RodriguezLogueK\":" << "\"\"" << ",";
        out << "\"RodriguezLogueBeta\":" << "\"\"" << ",";
        out << "\"RodriguezLogueRMS\":" << "\"\"" << ",";
        out << "\"RodriguezLogueavgerr\":" << "\"\"" << ",";
        out << "\"RodriguezLogueBIC\":" << "\"\"" << ",";
        out << "\"RodriguezLogueAIC\":" << "\"\"" << ",";
        out << "\"RodriguezLogueCode\":" << (int) mFitter.GetInfo() << ",";
    }

    /*
     * Ebert Prelec
     *
     * Work on start points here
     */

    startWriter.clear();
    startWriter.str("");

    p1Span = abs(-12) + abs(12); // -12 to 12
    p1Step = p1Span / 10;

    p2Span = 1; // -12 to 12
    p2Step = p2Span / 100;

    grandLoop = 0;

    for (int kLoop = 0; kLoop < 10; kLoop++)
    {
        for (int sLoop = 0; sLoop < 100; sLoop++)
        {
            provisionalValues.twoParamStartingValueArray[grandLoop].p1 = 12 - ((kLoop + 1) * p1Step);
            provisionalValues.twoParamStartingValueArray[grandLoop].p2 = ((sLoop + 1) * p2Step);

            grandLoop++;
        }
    }

    for (int mStep = 0; mStep < 100; mStep++)
    {
    	provisionalValues.oneParamStartingValueArray[mStep].err = mFitter.getErrorRodriguezLogue(provisionalValues.oneParamStartingValueArray[mStep].p1,
    			provisionalValues.oneParamStartingValueArray[mStep].p2);
    }

    std::sort(provisionalValues.twoParamStartingValueArray, provisionalValues.twoParamStartingValueArray + 1000);

    startWriter << "[" << provisionalValues.twoParamStartingValueArray[0].p1 << "," << provisionalValues.twoParamStartingValueArray[0].p2 << "]";

    /*
     * Perform fitting with start values
     */

    mFitter.FitEbertPrelec(startWriter.str().c_str());

    /*
     * Report back results
     */

    if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
    {
        out << "\"EbertPrelecK\":" << mFitter.fitEbertPrelecK << ",";
        out << "\"EbertPrelecS\":" << mFitter.fitEbertPrelecS << ",";
        out << "\"EbertPrelecRMS\":" << mFitter.GetReport().rmserror << ",";
        out << "\"EbertPrelecavgerr\":" << mFitter.GetReport().avgerror << ",";
        out << "\"EbertPrelecBIC\":" << mFitter.bicEbertPrelec << ",";
        out << "\"EbertPrelecAIC\":" << mFitter.aicEbertPrelec << ",";
        out << "\"EbertPrelecCode\":" << (int) mFitter.GetInfo() << ",";

        mFitter.mBicList.push_back(std::pair<std::string, double>("EbertPrelec", mFitter.bicEbertPrelec));
    }
    else
    {
        out << "\"EbertPrelecK\":" << "\"\"" << ",";
        out << "\"EbertPrelecS\":" << "\"\"" << ",";
        out << "\"EbertPrelecRMS\":" << "\"\"" << ",";
        out << "\"EbertPrelecavgerr\":" << "\"\"" << ",";
        out << "\"EbertPrelecBIC\":" << "\"\"" << ",";
        out << "\"EbertPrelecAIC\":" << "\"\"" << ",";
        out << "\"EbertPrelecCode\":" << (int) mFitter.GetInfo() << ",";
    }

    /*
     * Beleichrodt
     *
     * Work on start points here
     */

    startWriter.clear();
    startWriter.str("");

    p1Span = abs(-12) + abs(12); // -12 to 12
    p1Step = p1Span / 10;

    p2Span = 1; //
    p2Step = p2Span / 10;

    p3Span = 1; //
    p3Step = p3Span / 10;

    grandLoop = 0;

    for (int kLoop = 0; kLoop < 10; kLoop++)
    {
        for (int sLoop = 0; sLoop < 10; sLoop++)
        {
        	for (int bLoop = 0; bLoop < 10; bLoop++)
        	{
                provisionalValues.twoParamStartingValueArray[grandLoop].p1 = 12 - ((kLoop + 1) * p1Step);
                provisionalValues.twoParamStartingValueArray[grandLoop].p2 = ((sLoop + 1) * p2Step);
                provisionalValues.twoParamStartingValueArray[grandLoop].p3 = ((bLoop + 1) * p3Step);

                grandLoop++;
        	}
        }
    }

    for (int mStep = 0; mStep < 100; mStep++)
    {
    	provisionalValues.oneParamStartingValueArray[mStep].err = mFitter.getErrorBleichrodt(provisionalValues.oneParamStartingValueArray[mStep].p1,
    			provisionalValues.oneParamStartingValueArray[mStep].p2,
    			provisionalValues.oneParamStartingValueArray[mStep].p3);
    }

    std::sort(provisionalValues.twoParamStartingValueArray, provisionalValues.twoParamStartingValueArray + 1000);

    startWriter << "[" << provisionalValues.twoParamStartingValueArray[0].p1 << "," << provisionalValues.twoParamStartingValueArray[0].p2 << "]";

    /*
     * Perform fitting with start values
     */

    mFitter.FitBleichrodt(startWriter.str().c_str());

    /*
     * Report back results
     */

    if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
    {
    	out << "\"BleichrodtBeta\":" << mFitter.fitBleichrodtBeta << ",";
    	out << "\"BleichrodtK\":" << mFitter.fitBleichrodtK << ",";
        out << "\"BleichrodtS\":" << mFitter.fitBleichrodtS << ",";
        out << "\"BleichrodtRMS\":" << mFitter.GetReport().rmserror << ",";
        out << "\"Bleichrodtavgerr\":" << mFitter.GetReport().avgerror << ",";
        out << "\"BleichrodtBIC\":" << mFitter.bicBleichrodt << ",";
        out << "\"BleichrodtAIC\":" << mFitter.aicBleichrodt << ",";
        out << "\"BleichrodtCode\":" << (int) mFitter.GetInfo() << ",";

        mFitter.mBicList.push_back(std::pair<std::string, double>("Bleichrodt", mFitter.bicBleichrodt));
    }
    else
    {
    	out << "\"BleichrodtBeta\":" << "\"\"" << ",";
    	out << "\"BleichrodtK\":" << "\"\"" << ",";
        out << "\"BleichrodtS\":" << "\"\"" << ",";
        out << "\"BleichrodtRMS\":" << "\"\"" << ",";
        out << "\"Bleichrodtavgerr\":" << "\"\"" << ",";
        out << "\"BleichrodtBIC\":" << "\"\"" << ",";
        out << "\"BleichrodtAIC\":" << "\"\"" << ",";
        out << "\"BleichrodtCode\":" << (int) mFitter.GetInfo() << ",";
    }

    mFitter.FitNoise();

    out << "\"NoiseMean\":" << mFitter.AVE << ",";
    out << "\"NoiseRMS\":" << sqrt(mFitter.S2) << ",";
    out << "\"Noiseavgerr\":" << mFitter.sumErr << ",";
    out << "\"NoiseBIC\":" << mFitter.bicNoise << ",";
    out << "\"NoiseAIC\":" << mFitter.aicNoise << ",";
    out << "\"NoiseCode\":" << (int) mFitter.GetInfo() << ",";

    mFitter.mBicList.push_back(std::pair<std::string, double>("Noise", mFitter.bicNoise));

    mFitter.PrepareProbabilities();

    out << "\"HyperbolicBF\":" << mFitter.bfHyperbolic << ",";
    out << "\"ExponentialBF\":" << mFitter.bfExponential << ",";
    out << "\"QuasiHyperbolicBF\":" << mFitter.bfQuasiHyperbolic << ",";
    out << "\"MyersonBF\":" << mFitter.bfMyerson << ",";
    out << "\"RachlinBF\":" << mFitter.bfRachlin << ",";
    out << "\"RodriguezLogueBF\":" << mFitter.bfRodriguezLogue << ",";
    out << "\"EbertPrelecBF\":" << mFitter.bfEbertPrelec << ",";
    out << "\"BleichrodtBF\":" << mFitter.bfBleichrodt << ",";
    out << "\"NoiseBF\":" << mFitter.bfNoise << ",";

    out << "\"HyperbolicProb\":" << mFitter.probsHyperbolic << ",";
    out << "\"ExponentialProb\":" << mFitter.probsExponential << ",";
    out << "\"QuasiHyperbolicProb\":" << mFitter.probsQuasiHyperbolic << ",";
    out << "\"MyersonProb\":" << mFitter.probsMyerson << ",";
    out << "\"RachlinProb\":" << mFitter.probsRachlin << ",";
    out << "\"RodriguezLogueProb\":" << mFitter.probsRodriguezLogue << ",";
    out << "\"EbertPrelecProb\":" << mFitter.probsEbertPrelec << ",";
    out << "\"BleichrodtProb\":" << mFitter.probsBleichrodt << ",";
    out << "\"NoiseProb\":" << mFitter.probsNoise << ",";

    out << "\"ProbableModel\":" << "\"" << mFitter.mProbList[0].first << "\",";
    out << "\"ProbableED50\":" << mFitter.getED50BestModel() << ",";
    out << "\"ProbableArea\":" << mFitter.getAUCBestModel() << ",";
    out << "\"ProbableAreaLog10\":" << mFitter.getLog10AUCBestModel();

    out << "}";

	cout << out.str() << endl;

	return 0;
}
