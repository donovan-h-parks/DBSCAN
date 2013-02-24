# get DBSCAN bins
fin = open('ph28_illumina_contigs_5k.clusters.txt')
data = fin.readlines()
fin.close()

bins = {}
binId = ''
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
		temp.add(contigId)
		bins[binId] = temp
			
# get scaffolds
fin = open('../ph28_illumina_scaffolds_5k.fna')
data = fin.readlines()
fin.close()

seqId = ''
seqs = {}
for line in data:
	if line[0] == '>':
		if seqId != '':
			seqs[seqId] = seq
		seqId = line[1:].strip()
		seqId = seqId[seqId.find('_')+1:]
		seq = ''
	else:
		seq += line
seqs[seqId] = seq
		
for key in bins.keys():
	bin = key.split()[1]
	fout = open('./bins/Bin.' + bin + '.fna', 'w')
	for id in bins[key]:
		fout.write('>' + id + '\n')
		fout.write(seqs[id])
	fout.close()