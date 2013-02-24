# This is a simple test to see if the mean of a
# large number of independently and identically
# distributed Poisson random variables is well
# approximated by a Normal distrubtion with
# mean u and variance v/n.

import numpy as np
import matplotlib.pyplot as plt

trials = 10000
numRV = 1000
u = 5.0

poissonDist = []
for i in xrange(0, trials):
	sumRV = 0.0
	for j in xrange(0, numRV):
		sumRV += np.random.poisson(u)
	poissonDist.append(sumRV/numRV)
	
normalDist = []
for i in xrange(0, trials):
	normalDist.append(np.random.normal(u, np.sqrt(u/numRV)))

#testDist = []
#for i in xrange(0, trials):
#	testDist.append(np.random.poisson(u))
	
count, bins, ignored = plt.hist(poissonDist, bins=30, normed=True, align='mid', alpha=0.5, color='red', label='Mean Poisson RVs')
count, bins, ignored = plt.hist(testDist, bins=30, normed=True, align='mid', alpha=0.5, color='blue', label='Normal RV')
plt.legend(loc='upper right')
plt.show()