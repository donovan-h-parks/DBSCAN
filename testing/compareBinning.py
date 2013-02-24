#get ESOM bins
fin = open('ph28_esom_bin_contig_list.txt')
data = fin.readlines()
fin.close()

esomBins = {}
for line in data:
	lineSplit = line.split(':>')
	
	binId = lineSplit[0]
	binId = binId[0:binId.rfind('.')]
	
	contigId = lineSplit[1].split()[0]
	contigId = contigId[contigId.rfind('_')+1:]
	
	temp = esomBins.get(binId, set([]))
	temp.add(contigId)
	esomBins[binId] = temp
	
# get DBSCAN bins
fin = open('ph28_illumina_contigs_5k.clusters.txt')
data = fin.readlines()
fin.close()

bins = {}
binId = ''
bNoise = False
noiseBin = set([])
for line in data:
	if line.strip() == '':
		continue
		
	if 'Cluster' in line and 'Seq. Length' in line:
		binId = line.split('(')[0].strip()
	elif 'Noise' in line:
		bNoise = True
	elif binId != '':
		contigId = line.split('_')[1]
		
		if bNoise:
			noiseBin.add(contigId)
		else:
			temp = bins.get(binId, set([]))
			temp.add(contigId)
			bins[binId] = temp
		
# create table of contig mappings
fout = open('bin_mapping.tsv', 'w')
for esomId in esomBins.keys():
	fout.write('\t' + esomId)
fout.write('\tUnassigned\n')

for binId in bins.keys():
	fout.write(binId)
	unassignedContigs = bins[binId]
	for esomId in esomBins.keys():
		fout.write('\t' + str(len(bins[binId].intersection(esomBins[esomId]))))
		unassignedContigs = unassignedContigs - esomBins[esomId]
	fout.write('\t' + str(len(unassignedContigs)))
	fout.write('\n')
	
inNoise = []
fout.write('Unassigned')
for esomId in esomBins.keys():
	unassignedContigs = esomBins[esomId]
	for binId in bins.keys():
		unassignedContigs = unassignedContigs - bins[binId]
	fout.write('\t' + str(len(unassignedContigs)))
	inNoise.append(len(unassignedContigs.intersection(noiseBin)))
fout.write('\n')

fout.write('In Noise:')
for c in inNoise:
	fout.write('\t' + str(c))
fout.write('\n')
	
fout.close()