import random

numClusters = 100
numPtsPerCluster = 250

fout = open('random_' + str(numClusters) + 'by' + str(numPtsPerCluster) + '.tsv', 'w')
fout.write('Id\tAverage Coverage\tx\ty\n')
for c in xrange(0, numClusters):
	xCentre = random.randint(0, numClusters)
	yCentre = random.randint(0, numClusters)
	
	xStd = random.uniform(0.05, 0.15)
	yStd = random.uniform(0.05, 0.15)
	
	for i in xrange(0, numPtsPerCluster):
		x = random.gauss(xCentre, xStd)
		y = random.gauss(yCentre, yStd)
		fout.write(str(c) + '_' + str(i) + '\t1\t' + str(x) + '\t' + str(y) + '\n')

fout.close()
