import numpy as np
import matplotlib.pyplot as plt

G = 5000
L = 250
N = 1000
c = float(L*N)/G
print c

x = [0]*(G+L)
for i in xrange(0, N):
	t = np.random.randint(0, G)
	for j in xrange(t, t+L):
		x[j] += 1
		
avgCoverage = 0.0
for i in xrange(0, G):
	avgCoverage += x[i]
	
print str(avgCoverage/G)

p = list(np.random.poisson(c, G))

for i in xrange(int(c)-10, int(c)+10):
	print str(i) + ': ' + str(p.count(i)) + '   ' + str(x.count(i))
	
print float(sum(p))/len(p)

#bins = max(x) #max(max(x), max(p))
#count, bins, ignored = plt.hist(x[0:G], bins=30, normed=True, align='mid', alpha=0.5, color='red', label='Simulated')
#count, bins, ignored = plt.hist(p, bins=18, normed=True, alpha=0.5, color='blue', label='Poisson RV')
#plt.legend(loc='upper right')
#plt.show()
