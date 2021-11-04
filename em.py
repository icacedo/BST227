import math
import random
import sys

def make_matrix(seqs):

	def shape(fp):

		shape = []
		n = 0 
		N = 0
		count = 0
		for line in fp.readlines():
			line = line.rstrip()
			N += 1
			if count == 0:
				for ACTG in line:
					n += 1
			count += 1
		shape.append(n)
		shape.append(N)
		return shape 
		
	fp = open(seqs, 'r')
	shape = shape(fp) #shape = [n, N]
	fp.close()

	fp = open(seqs, 'r') # needs to be opened twice
	matrix = [[None for i in range(shape[0])] for j in range(shape[1])]
	#Note: matrix[row][column][A/C/G/T]

	i = 0
	for line in fp.readlines():
		line = line.rstrip()
		j = 0
		for n in line:
			if n.upper() == 'A':
				matrix[i][j] = [1,0,0,0]
			if n.upper() == 'C':
				matrix[i][j] = [0,1,0,0]
			if n.upper() == 'G':
				matrix[i][j] = [0,0,1,0]
			if n.upper() == 'T':
				matrix[i][j] = [0,0,0,1]
			j += 1
		i += 1
	return matrix

def freq_dist(length, decimal):

	vector = [None] * length
	isone = 0
	while (isone != 1):
		total = 0 
		isone = 0
		for i in range(len(vector)):
			freq = float(random.uniform(0, 1))
			vector[i] = freq
			total += freq
		for i in range(len(vector)):
			vector[i] = round(vector[i]/total, decimal)
		for i in range(len(vector)):
			isone += vector[i]
		if isone == 1: 
			return vector
			break
			
def psi(length, decimal):
	psi = []
	for i in range(4):
		psi.append(freq_dist(length, decimal))
	return psi
	
def shape(matrix):
	rows = 0
	for i in seq_mx:
		rows += 1
		columns = 0
		for j in i:
			columns += 1
	return (rows, columns)
	
# initialize psis
# 1 is foreground, 0 is background
# matrix of parameters size 4 x P
# P = 6, length of motif
# dec means number of decimal places
dec = 2
psi_0mx = psi(6, dec)
psi_1mx = psi(6, dec)

# read in data to matrix and get matrix shape
seq_mx = make_matrix(sys.argv[1])
seq_mx_shape = shape(seq_mx)
# rows aka sequences
rows = seq_mx_shape[0]
# columns aka positions
columns = seq_mx_shape[1]

# specify size of motif
P = 6

# initialize lambdas
lambdajs = []
for i in range(columns):
	lambdajs.append(round(random.uniform(0, 1), dec))

# A=0, C=1, G=2, T=3
def unencode(encode):
	if encode == [1, 0, 0, 0]:
		return 0
	if encode == [0, 1, 0, 0]:
		return 1
	if encode == [0, 0, 1, 0]:
		return 2
	if encode == [0, 0, 0, 1]:
		return 3
		
# position j is where the motif starts
# i for each sequence
numerators = []
for i in range(rows):
	# at each position
	# motif can start at up to L-P+1 positions
	Cijs_list = []
	for j in range(columns-(P-1)):
		# probability that Cij is the start of a motif
		# that probability is lambdaj
		# times the sum of probabilities for psi
		Cij = lambdajs[j] 
		for p in range(P):
			Xijp = unencode(seq_mx[i][j+p])
			Cij *= psi_1mx[Xijp][p] 
		for jp in range(columns-(P-1)):
			# jp is j'
			# to get background frequencies
			# they cannot come from the same j as the foreground
			if jp == j: continue
			for p in range(P):	
				Xijp = unencode(seq_mx[i][j+p])
				Cij *= psi_0mx[Xijp][p]
		Cijs_list.append(Cij)
	numerators.append(Cijs_list)

denominators = []
for freqs in numerators:
	total = 0
	for Cij in freqs:
		total += Cij
	denominators.append(total)

posteriors = []
for i in range(len(numerators)):
	one_row = []
	for j in range(len(numerators[i])):
		# get a domain error only sometimes
		# what is meant by subtract by the log of the smallest number?
		# this doesn't work
		# smallest = min(numerators[i][j], denominators[i])
		# idk how to implement the logsumexp trick
		one_row.append(numerators[i][j] \
			/denominators[i])
	posteriors.append(one_row)
	


				


























	
