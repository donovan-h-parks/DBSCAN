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

#ifndef _POISSION_
#define _POISSION_

#include "DataTypes.hpp"

class Poisson
{
public:
	Poisson(double lambda): m_lambda(lambda) {}

	static bool Test(double t1, double t2, double success1, double success2, double criticalValue);

	double Mean() { return m_lambda; }
	double Mode() { return floor(m_lambda); }
	double Variance() { return m_lambda; }

	double PDF(uint k);
	double CDF(uint k);

	// Finds the smallest k such that P(X <= k) is greater than or equal to p.
	// The PPF function rounded outwards. That is to say lower quantiles (where 
	// the probability is less than 0.5) are rounded downward, and upper 
	// quantiles (where the probability is greater than 0.5) are rounded 
	// upwards.
	uint PPF(double p); // aka, inverse of CDF

private:
	double Factorial(uint n);

	double m_lambda;
};

#endif