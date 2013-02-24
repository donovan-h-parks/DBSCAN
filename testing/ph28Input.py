fin = open('ph28_nb_kmer_table.contigs_5k.tsv')
data = fin.readlines()
fin.close()

header = data[0].split('\t')
header = [header[0]] + ['GC\tSeq. Length\tAverage Coverage'] + header[1:]
header = '\t'.join(header)
freqs = {}
for i in xrange(1, len(data)):
	lineSplit = data[i].split('\t')
	contigId = lineSplit[0].split()[0]
	freq = lineSplit[1:]
	
	freqs[contigId] = freq
	
fin = open('ph28_illumina_contigs_5k.avgCoverage.tsv')
data = fin.readlines()
fin.close()

contigInfo = {}
for line in data:
	lineSplit = line.split(',')
	contigId = lineSplit[0]
	seqLen = lineSplit[1]
	gc = lineSplit[2]
	coverage = lineSplit[3].strip()
	
	contigInfo[contigId] = [gc, seqLen, coverage]
	
fout = open('ph28_illumina_contigs_nb_5k.tsv', 'w')
fout.write(header)
for key in freqs:
	freqStr = '\t'.join(freqs[key])
	if key in contigInfo:
		fout.write(key + '\t' + contigInfo[key][0] + '\t' + contigInfo[key][1] + '\t' + contigInfo[key][2] + '\t' + freqStr)
fout.close()