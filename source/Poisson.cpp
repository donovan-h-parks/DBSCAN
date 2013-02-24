//=======================================================================
// Author: Donovan Parks
//
// Copyright 2013 Donovan Parks
//
// This file is part of DBSCAN.
//
// DBSCAN is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// DBSCAN is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public 
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with DBSCAN. If not, see 
// <http://www.gnu.org/licenses/>.
//=======================================================================

#include "Precompiled.hpp"

#include "Poisson.hpp"

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>

double Poisson::PDF(uint k)
{
	boost::math::poisson_distribution<> poissonDist(m_lambda);
	return pdf(poissonDist, k);
}

double Poisson::CDF(uint k)
{
	boost::math::poisson_distribution<> poissonDist(m_lambda);
	return cdf(poissonDist, k);
}

uint Poisson::PPF(double p)
{
	boost::math::poisson_distribution<> poissonDist(m_lambda);
	return uint(quantile(poissonDist, p));
}

bool Poisson::Test(double t1, double t2, double success1, double success2, double criticalValue)
{
	double c = t2/t1;
	double numerator = c*success1 - success2;
	double denomenator = sqrt(c*c*success1 + success2);

	boost::math::normal_distribution<> normal(0, 1);
	double Z = -boost::math::quantile(normal, criticalValue);

	double diff = success1/t1 - success2/t2;
	double lowerCI = (numerator - Z*denomenator)/t2;
	double upperCI = (numerator + Z*denomenator)/t2;

	return diff >= lowerCI && diff <= upperCI;
}