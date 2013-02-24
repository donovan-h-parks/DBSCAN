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

#include "Point.hpp"
#include "DBSCAN.hpp"

#include "getopt_pp.hpp"

#include "Poisson.hpp"

int main(int argc, char *argv[])
{
	double totalElapsedTime = 0;

	// parse command line
	bool bShowHelp, bVerbose;
	std::string inputFileStr, epsStr, minPtsStr, coverageBoundStr, outputFileStr;
	GetOpt::GetOpt_pp opts(argc, argv);
	opts >> GetOpt::OptionPresent('h', "help", bShowHelp);
	opts >> GetOpt::OptionPresent('v', "verbose", bVerbose);
	opts >> GetOpt::Option('i', "input-file", inputFileStr);
	opts >> GetOpt::Option('e', "eps", epsStr, "0");
	opts >> GetOpt::Option('m', "minPts", minPtsStr, "4");
	opts >> GetOpt::Option('c', "coverage", coverageBoundStr, "0.025");
	opts >> GetOpt::Option('o', "output-file ", outputFileStr, "clusters.txt");

	float eps = float(atof(epsStr.c_str()));
	int minPts = atoi(minPtsStr.c_str());
	float coverageBound = float(atof(coverageBoundStr.c_str()));

	if(bShowHelp || argc <= 1) 
	{		
		std::cout << std::endl;
		std::cout << "DBSCAN v1.0.0 (Feb. 5, 2013)" << std::endl;
		std::cout << "  by Donovan Parks (donovan.parks@gmail.com)" << std::endl;
		std::cout << std::endl;
		std::cout << " Usage: " << opts.app_name() << " -i <input file> -e <eps> -m <minimum pts> -o <output file>" << std::endl;
		std::cout << "  -h, --help        Produce help message." << std::endl;
		std::cout << "  -i, --input-file  Input file specifying k-mer frequencies." << std::endl;
		std::cout << "  -e, --eps         Eps threshold used for defining clusters." << std::endl;
		std::cout << "  -m, --minPts      Minimum points threshold used for defining clusters (default = 4)." << std::endl;
		std::cout << "  -c, --coverage    Constrain clusters to sequences with similar coverage based on the Poisson distribution (default = 0.025)." << std::endl;
		std::cout << "  -o, --output-file Output file indicating clustering." << std::endl;
		std::cout << std::endl;
		std::cout << "  -v, --verbose     Provide additional information on program execution." << std::endl;
	}

	if(epsStr == "0")
	{
		std::cout << "The eps threshold (-e) must be specified. Use -h for help." << std::endl;
		return 0;
	}

	if(bVerbose)
	{
		std::cout << "Parameters:" << std::endl;
		std::cout << "  Eps: " << eps << std::endl;
		std::cout << "  Min. points: " << minPts << std::endl;
		std::cout << "  Coverage bound: " << coverageBound << std::endl << std::endl;
	}

	// read kmer frequencies and create points to cluster
	if(bVerbose)
		std::cout << "Reading in k-mer frequencies." << std::endl;

	std::clock_t start = std::clock();
	std::ifstream fin(inputFileStr.c_str());
	if(!fin.is_open())
	{
		std::cerr << "Failed to open input file: " << inputFileStr << std::endl;
		return -1;
	}

	DBSCAN dbscan;

	bool bHeader = true;
	std::string line;
	while(std::getline(fin, line, '\n'))
	{
		if(bHeader)
		{
			bHeader = false;
			continue;
		}

		Point pt;
		std::vector<float> freqs;
		std::string token;
		std::istringstream data(line);
		while(std::getline(data, token, '\t'))
		{
			if(pt.GetName().empty())
				pt.SetName(token);
			else if(pt.GetGC() == Point::NO_INFO)
				pt.SetGC(float(atof(token.c_str())));
			else if(pt.GetLength() == Point::NO_INFO)
				pt.SetLength(atoi(token.c_str()));
			else if(pt.GetCoverage() == Point::NO_INFO)
				pt.SetCoverage(float(atof(token.c_str())));
			else
				freqs.push_back(float(atof(token.c_str())));
		}

		pt.SetKmerFreqs(freqs);

		dbscan.AddPoint(pt);
	}

	fin.close();

	std::clock_t end = std::clock();
	double elapsedTime = ( end - start ) / (double)CLOCKS_PER_SEC;
	totalElapsedTime += elapsedTime;
	if(bVerbose)
		std::cout << "  Elapsed time (s): " << elapsedTime << std::endl;

	// perform clustering
	if(bVerbose)
		std::cout << "Performing clustering." << std::endl;

	start = std::clock();

	dbscan.Run(eps, minPts, coverageBound);
	std::vector< std::vector<Point*> > clusters = dbscan.GetClusters();
	std::vector<Point*> noisePts = dbscan.GetNoisePts();

	end = std::clock();
	elapsedTime = ( end - start ) / (double)CLOCKS_PER_SEC;
	totalElapsedTime += elapsedTime;
	if(bVerbose)
		std::cout << "  Elapsed time (s): " << elapsedTime << std::endl;

	// get coverage information for each cluster
	uint ptsInClusters = 0;
	std::vector<ClusterStats> clusterStats(clusters.size());

	for(uint i = 0; i < clusters.size(); ++i)
	{
		ptsInClusters += clusters[i].size();

		float uGC = 0.0f;
		float uLength = 0.0f;
		float uCoverage = 0.0f;
		for(uint j = 0; j < clusters[i].size(); ++j)
		{
			uGC += clusters[i][j]->GetGC();
			uLength += clusters[i][j]->GetLength();
			uCoverage += clusters[i][j]->GetCoverage();
		}
		uGC /= clusters[i].size();
		uLength /= clusters[i].size();
		uCoverage /= clusters[i].size();

		float sGC = 0.0f;
		float sLength = 0.0f;
		float sCoverage = 0.0f;
		for(uint j = 0; j < clusters[i].size(); ++j)
		{
			float diff = clusters[i][j]->GetGC() - uGC;
			sGC += diff*diff;

			diff = clusters[i][j]->GetLength() - uLength;
			sLength += diff*diff;

			diff = clusters[i][j]->GetCoverage() - uCoverage;
			sCoverage += diff*diff;
		}
		sGC = sqrt(sGC/clusters[i].size());
		sLength = sqrt(sLength/clusters[i].size());
		sCoverage = sqrt(sCoverage/clusters[i].size());

		clusterStats[i].meanLength = uLength;
		clusterStats[i].stdLength = sLength;
		clusterStats[i].meanGC = uGC;
		clusterStats[i].stdGC = sGC;
		clusterStats[i].meanCoverage = uCoverage;
		clusterStats[i].stdCoverage = sCoverage;
	}

	// write out clustering results
	std::ofstream fout(outputFileStr.c_str());
	fout.precision(3);
	fout.flags(std::ios::fixed);
	fout << "Input parameters:" << std::endl;
	fout << "  Input file: " << inputFileStr << std::endl;
	fout << "  Eps: " << eps << std::endl;
	fout << "  Min. Pts: " << minPts << std::endl;
	fout << "  Coverage bound: " << coverageBound << std::endl << std::endl;

	fout << "Clustering info:" << std::endl;
	fout << "  Number of clusters: " << clusters.size() << std::endl;
	fout << "  Points to cluster: " << dbscan.GetNumPoints() << std::endl;
	fout << "  Points contained in clusters: " << ptsInClusters << std::endl;
	fout << "  Points considered noise: " << noisePts.size() << std::endl;
	fout << "  Total time to cluster (minutes): " << totalElapsedTime/60.0 << std::endl << std::endl;

	fout << "Cluster summary:" << std::endl;
	fout << "Cluster Id\tPoint in Cluster\tAverage length +/- std\tAverage GC +/- std\tAverage coverage +/- std" << std::endl;
	std::sort(clusterStats.begin(), clusterStats.end(), ClusterStats::SortClustersByCoverage);
	for(uint i = 0; i < clusters.size(); ++i)
	{
		fout << "Cluster " << i << "\t" << clusters[i].size() << "\t" << clusterStats[i].meanLength << " +/- " << clusterStats[i].stdLength;
		fout << "\t" << clusterStats[i].meanGC << " +/- " << clusterStats[i].stdGC;
		fout << "\t" << clusterStats[i].meanCoverage << " +/- " << clusterStats[i].stdCoverage <<std::endl;
	}
	fout << std::endl;

	for(uint i = 0; i < clusters.size(); ++i)
	{
		fout << "Cluster " << i << " (" << clusters[i].size() << " points, Seq. Length: " << clusterStats[i].meanLength << " +/- " << clusterStats[i].stdLength;
		fout << ", GC: " << clusterStats[i].meanGC << " +/- " << clusterStats[i].stdGC;
		fout << ", Coverage: " << clusterStats[i].meanCoverage << " +/- " << clusterStats[i].stdCoverage << "):" << std::endl;

		std::sort(clusters[i].begin(), clusters[i].end(), Point::SortPointsByCoverage);
		for(uint j = 0; j < clusters[i].size(); ++j)
			fout << clusters[i][j]->GetName() << "\t" << clusters[i][j]->GetLength() << "\t" << clusters[i][j]->GetGC() << "\t" << clusters[i][j]->GetCoverage() << std::endl;
		fout << std::endl;
	}

	fout << "Noise (" << noisePts.size() << " points, Seq. Length, GC, Coverage):" << std::endl;
	std::sort(noisePts.begin(), noisePts.end(), Point::SortPointsByLength);
	for(uint i = 0; i < noisePts.size(); ++i)
		fout << noisePts[i]->GetName() << "\t" << noisePts[i]->GetLength() << "\t" << noisePts[i]->GetGC() << "\t" << noisePts[i]->GetCoverage() << std::endl;
	fout << std::endl;

	fout.close();

	return 0;
}
