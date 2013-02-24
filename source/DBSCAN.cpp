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

#include "DBSCAN.hpp"

#include "Poisson.hpp"

const int DBSCAN::NOISE = std::numeric_limits<int>::min();

void DBSCAN::Run(float eps, uint minPts, float coverageBound)
{
	m_eps = eps;
	m_minPts = minPts;
	m_coverageBound = coverageBound;

	// cluster points
	m_visited.resize(m_pts.size(), false);
	m_clustered.resize(m_pts.size(), false);

	for(uint i = 0; i < GetNumPoints(); ++i)
	{
		if(m_visited[i])
			continue;

		m_visited[i] = true;

		std::set<uint> neighbours = FindNeighbours(i);
		if(neighbours.size() >= m_minPts)
		{
			std::vector<Point*> cluster = ExpandCluster(i, neighbours);
			m_clusters.push_back(cluster);
		}
	}

	// identify unclustered points and label them as noise
	for(uint i = 0; i < m_clustered.size(); ++i)
	{
		if(!m_clustered[i])
			m_noise.push_back(&m_pts[i]);
	}
}

std::vector<Point*> DBSCAN::ExpandCluster(const uint ptIndex, std::set<uint>& neighbours)
{
	std::vector<Point*> cluster;
	
	cluster.push_back(&m_pts[ptIndex]);
	m_clustered[ptIndex] = true;

	while(!neighbours.empty())
	{
		std::set<uint>::iterator it = neighbours.begin();
		uint neighbourIndex = *it;
		neighbours.erase(it);

		Point& neighbourPt = m_pts[neighbourIndex];
		
		if(!m_visited[neighbourIndex])
		{
			m_visited[neighbourIndex] = true;

			std::set<uint> curNeighbours = FindNeighbours(neighbourIndex);
			if(curNeighbours.size() >= m_minPts) // check if core point
				neighbours.insert(curNeighbours.begin(), curNeighbours.end());
		}

		if(!m_clustered[neighbourIndex])
		{
			cluster.push_back(&neighbourPt);
			m_clustered[neighbourIndex] = true;
		}
	}

	return cluster;
}

/*
Find all points within P's eps-neighborhood
*/
std::set<uint> DBSCAN::FindNeighbours(const uint ptIndex)
{
	std::set<uint> neighbours;

	std::vector<float> freq1 = m_pts[ptIndex].GetKmerFreqs();

	for(uint i = 0; i < m_pts.size(); ++i)
	{
		if(ptIndex == i)
			continue;

		std::vector<float> freq2 = m_pts[i].GetKmerFreqs();

		if(CoverageCheck(m_pts[ptIndex].GetCoverage(), m_pts[i].GetCoverage(), m_pts[ptIndex].GetLength(), m_pts[i].GetLength()))
		{
			float dist = CalculateDistance(freq1, freq2);
			if(dist < m_eps)
				neighbours.insert(i);
		}
	}

	return neighbours;
}

bool DBSCAN::CoverageCheck(float avgCoverage1, float avgCoverage2, uint len1, uint len2)
{
	float avgCoverage = avgCoverage1;
	float testCoverage = avgCoverage2;
	if(avgCoverage2 > avgCoverage1)
	{
		avgCoverage = avgCoverage2;
		testCoverage = avgCoverage1;
	}

	Poisson poisson(avgCoverage);
	int lowerBound = poisson.PPF(m_coverageBound);
	int upperBound = poisson.PPF(1.0 - m_coverageBound);

	if(testCoverage <= lowerBound || testCoverage > upperBound)
	{
		// Note: the inequalities above are correct. The CDF (and by proxy the PPF)
		// is inclusive of the random variable. For more details, see:
		// http://www.boost.org/doc/libs/1_35_0/libs/math/doc/sf_and_dist/html/math_toolkit/policy/pol_tutorial/understand_dis_quant.html
		return false;
	}

	return true;
}

float DBSCAN::CalculateDistance(const std::vector<float>& freq1, const std::vector<float>& freq2)
{
	/*
	float dist = 0.0f;
	for(uint i = 0; i < freq1.size(); ++i)
	{
		float diff = std::abs(freq1[i] - freq2[i]);
		dist += diff;
	}

	return dist/freq1.size();
	*/

	double logProb = 0.0f;
	for(uint i = 0; i < freq1.size(); ++i)
		logProb += freq1[i]*log(freq2[i]) + freq2[i]*log(freq1[i]);

	double nullProb = 2 * log(1.0/freq1.size());

	if (logProb > log(1.4) + nullProb)
		return 0;
	
	return 100000;
}
	
	

