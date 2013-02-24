# get contigs in each bin
fin = open('ph28_illumina_contigs_5k.clusters.txt')
data = fin.readlines()
fin.close()

bins = {}
binId = ''
scaffoldIds = set([])
multiHit = 0
for line in data:
	if line.strip() == '':
		continue
		
	if 'Cluster' in line and 'Seq. Length' in line:
		binId = line.split('(')[0].strip()
	elif 'Noise' in line:
		break
	elif binId != '':
		contigId = line.split('_')[1]
		
		temp = bins.get(binId, set([]))
		
		if contigId in temp:
			multiHit += 1
		
		temp.add(contigId)
		bins[binId] = temp
		
		scaffoldIds.add(contigId)
		
# check for contigs from same scaffold that span multiple bins
clusterCounter = []
for scaffoldId in scaffoldIds:
	clusterCount = 0
	for b in bins:
		if scaffoldId in bins[b]:
			clusterCount += 1
			
	clusterCounter.append(clusterCount)

# report results
for i in xrange(0, len(bins)):
	print 'Occuring in ' + str(i) + ' bins: ' + str(clusterCounter.count(i))
		
print 'Times contigs from same scaffold was assigned to the same bin: ' + str(multiHit)
print 'Total clustered contigs: ' + str(multiHit + sum(clusterCounter))