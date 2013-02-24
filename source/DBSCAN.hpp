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

#ifndef _DBSCAN_
#define _DBSCAN_

#include "DataTypes.hpp"

#include "Point.hpp"

struct ClusterStats
{
	float meanLength;
	float stdLength;

	float meanGC;
	float stdGC;

	float meanCoverage;
	float stdCoverage;

	static bool SortClustersByCoverage(const ClusterStats& c1, const ClusterStats& c2)
	{
		return c1.meanCoverage > c2.meanCoverage;
	}
		
};

class DBSCAN
{
public:
	static const int NOISE;

public:
	DBSCAN() {}

	void DBSCAN::AddPoint(const Point& pt) { m_pts.push_back(pt); }

	void Run(float eps, uint minPts, float coverageBound);

	uint GetNumPoints() const { return m_pts.size(); }
	std::vector< std::vector<Point*> > GetClusters() const { return m_clusters; }
	std::vector<Point*> GetNoisePts() const { return m_noise; }

private:
	std::vector<Point*> ExpandCluster(const uint ptIndex, std::set<uint>& neighbours);
	std::set<uint> FindNeighbours(const uint ptIndex);

	inline float CalculateDistance(const std::vector<float>& freq1, const std::vector<float>& freq2);
	inline bool CoverageCheck(float avgCoverage1, float avgCoverage2, uint len1, uint len2);

private:
	float m_eps;
	float m_coverageBound;
	uint m_minPts;

	std::vector<Point> m_pts;

	std::vector< std::vector<Point*> > m_clusters;
	std::vector<Point*> m_noise;

	std::vector<bool> m_visited;
	std::vector<bool> m_clustered;
};

#endif
