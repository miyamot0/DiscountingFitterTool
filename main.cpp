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
#include "modelselection.h"

using namespace std;

int main(int argc, char *argv[])
{
	ModelSelection mFitter;
	mFitter.InitializeDefaults();

	mFitter.SetX(argv[1]);
    mFitter.SetY(argv[2]);
    mFitter.mBicList.clear();

    std::ostringstream out;
    out << "{";

    mFitter.FitHyperbolic("[0.3]");

    out << "\"Hyperbolic\":" << mFitter.fitHyperbolicK << ",";
    out << "\"HyperbolicRMS\":" << mFitter.GetReport().rmserror << ",";
    out << "\"HyperbolicBIC\":" << mFitter.bicHyperbolic << ",";
    out << "\"HyperbolicAIC\":" << mFitter.aicHyperbolic << ",";
    out << "\"HyperbolicCode\":" << (int) mFitter.GetInfo() << ",";

    mFitter.FitExponential("[0.3]");

    out << "\"Exponential\":" << mFitter.fitExponentialK << ",";
    out << "\"ExponentialRMS\":" << mFitter.GetReport().rmserror << ",";
    out << "\"ExponentialBIC\":" << mFitter.bicExponential << ",";
    out << "\"ExponentialAIC\":" << mFitter.aicExponential << ",";
    out << "\"ExponentialCode\":" << (int) mFitter.GetInfo() << ",";

    mFitter.FitQuasiHyperbolic("[0.3, 0.3]");

    out << "\"QuasiHyperbolicBeta\":" << mFitter.fitQuasiHyperbolicBeta << ",";
    out << "\"QuasiHyperbolicDelta\":" << mFitter.fitQuasiHyperbolicDelta << ",";
    out << "\"QuasiHyperbolicRMS\":" << mFitter.GetReport().rmserror << ",";
    out << "\"QuasiHyperbolicBIC\":" << mFitter.bicQuasiHyperbolic << ",";
    out << "\"QuasiHyperbolicAIC\":" << mFitter.aicQuasiHyperbolic << ",";
    out << "\"QuasiHyperbolicCode\":" << (int) mFitter.GetInfo() << ",";

    mFitter.FitMyerson("[0.3, 0.3]");

    out << "\"MyersonK\":" << mFitter.fitMyersonK << ",";
    out << "\"MyersonS\":" << mFitter.fitMyersonS << ",";
    out << "\"MyersonRMS\":" << mFitter.GetReport().rmserror << ",";
    out << "\"MyersonBIC\":" << mFitter.bicMyerson << ",";
    out << "\"MyersonAIC\":" << mFitter.aicMyerson << ",";
    out << "\"MyersonCode\":" << (int) mFitter.GetInfo() << ",";

    mFitter.FitRachlin("[0.3, 0.3]");

    out << "\"RachlinK\":" << mFitter.fitRachlinK << ",";
    out << "\"RachlinS\":" << mFitter.fitRachlinS << ",";
    out << "\"RachlinRMS\":" << mFitter.GetReport().rmserror << ",";
    out << "\"RachlinBIC\":" << mFitter.bicRachlin << ",";
    out << "\"RachlinAIC\":" << mFitter.aicRachlin << ",";
    out << "\"RachlinCode\":" << (int) mFitter.GetInfo() << ",";

    mFitter.FitNoise();

    out << "\"NoiseMean\":" << mFitter.AVE << ",";
    out << "\"NoiseRMS\":" << mFitter.GetReport().rmserror << ",";
    out << "\"NoiseBIC\":" << mFitter.bicNoise << ",";
    out << "\"NoiseAIC\":" << mFitter.aicNoise << ",";
    out << "\"NoiseCode\":" << (int) mFitter.GetInfo() << ",";

    mFitter.mBicList.push_back(std::pair<std::string, double>("Noise", mFitter.bicNoise));
    mFitter.mBicList.push_back(std::pair<std::string, double>("Exponential", mFitter.bicExponential));
    mFitter.mBicList.push_back(std::pair<std::string, double>("Hyperbolic", mFitter.bicHyperbolic));
    mFitter.mBicList.push_back(std::pair<std::string, double>("Beta Delta", mFitter.bicQuasiHyperbolic));
    mFitter.mBicList.push_back(std::pair<std::string, double>("Myerson", mFitter.bicMyerson));
    mFitter.mBicList.push_back(std::pair<std::string, double>("Rachlin", mFitter.bicRachlin));

    mFitter.PrepareProbabilities();

    out << "\"HyperbolicBF\":" << mFitter.bfHyperbolic << ",";
    out << "\"ExponentialBF\":" << mFitter.bfExponential << ",";
    out << "\"QuasiHyperbolicBF\":" << mFitter.bfQuasiHyperbolic << ",";
    out << "\"MyersonBF\":" << mFitter.bfMyerson << ",";
    out << "\"RachlinBF\":" << mFitter.bfRachlin << ",";
    out << "\"NoiseBF\":" << mFitter.bfNoise << ",";

    out << "\"HyperbolicProb\":" << mFitter.probsHyperbolic << ",";
    out << "\"ExponentialProb\":" << mFitter.probsExponential << ",";
    out << "\"QuasiHyperbolicProb\":" << mFitter.probsQuasiHyperbolic << ",";
    out << "\"MyersonProb\":" << mFitter.probsMyerson << ",";
    out << "\"RachlinProb\":" << mFitter.probsRachlin << ",";
    out << "\"NoiseProb\":" << mFitter.probsNoise;

    out << "}";

	cout << out.str() << endl;

	return 0;
}
