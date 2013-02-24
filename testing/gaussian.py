import random

numPtsPerCluster = 10
centers = [(0,0), (0,1), (1,0), (1,1)]
std = [(0.15, 0.1), (0.1, 0.15), (0.15, 0.1), (0.1, 0.15)]

fout = open('gaussian4.txt', 'w')
fout.write('Id\tAverage Coverage\tx\ty\n')
for c in xrange(0, len(centers)):
	for i in xrange(0, numPtsPerCluster):
		x = random.gauss(centers[c][0], std[c][0])
		y = random.gauss(centers[c][1], std[c][1])
		fout.write(str(c) + '_' + str(i) + '\t1\t' + str(x) + '\t' + str(y) + '\n')

fout.close()
