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

#ifndef _POINT_
#define _POINT_

#include "DataTypes.hpp"

class Point
{
public:
	static const int NO_INFO = -1;

public:
	static bool SortPointsByLength(const Point* pt1, const Point* pt2)
	{
		return pt1->GetLength() > pt2->GetLength();
	}

	static bool SortPointsByCoverage(const Point* pt1, const Point* pt2)
	{
		return pt1->GetCoverage() > pt2->GetCoverage();
	}

public:
	Point(): m_gc(float(NO_INFO)), m_length(NO_INFO), m_coverage(float(NO_INFO)) {}

	void SetName(const std::string& name) { m_name = name; }
	std::string GetName() const { return m_name; }

	void SetKmerFreqs(const std::vector<float>& freqs) { m_freqs = freqs; }
	std::vector<float> GetKmerFreqs() const { return m_freqs; }

	void SetGC(float gc) { m_gc = gc; }
	float GetGC() const { return m_gc; }

	void SetLength(uint length) { m_length = length; }
	uint GetLength() const { return m_length; }

	void SetCoverage(float coverage) { m_coverage = coverage; }
	float GetCoverage() const { return m_coverage; }

private:
	std::string m_name;

	std::vector<float> m_freqs;
	float m_gc;
	uint m_length;
	float m_coverage;
};

#endif