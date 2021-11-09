import math
import random
import sys

# i know i am not using one hot encoding correctly
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
		Cij = math.log(lambdajs[j]) 
		for p in range(P):
			Xijp = unencode(seq_mx[i][j+p])
			# logsumexp trick
			# change from *= to +=
			Cij += math.log(psi_1mx[Xijp][p]) 
		for jp in range(columns-(P-1)):
			# jp is j'
			# to get background frequencies
			# they cannot come from the same j as the foreground
			if jp == j: continue
			for p in range(P):	
				Xijp = unencode(seq_mx[i][j+p])
				# logsumexp
				Cij += math.log(psi_0mx[Xijp][p])
		Cijs_list.append(math.exp(Cij))
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
		one_row.append(numerators[i][j] \
			/denominators[i])
	posteriors.append(one_row)

#######################################
# posteriors[0] length of 33 is correct

# adding back in zeroes 
# getting indexing error in m step
# i think there is an indexing error somewhere
# code gets sequence length of 39, not 38
# might be ignorable for now
for i in range(len(posteriors)):
	# first position in motif length should not be 0
	for p in range(P-1):
		posteriors[i].append(0.0001)

# according to a1sol.pdf, posteriors = entropy?
# is lambaj_hat the sum of all posteriors at one position?

column_sums = [0 for i in range(len(posteriors[0])-P+1)]
for i in range(len(posteriors)):
	for j in range(len(posteriors[i])-P+1):
		column_sums[j] += posteriors[i][j]
		
lambdaj_hats = []
for i in range(len(column_sums)):
	N = len(posteriors)
	lambdaj_hats.append(column_sums[i]/N)

position_sumsA = [0 for i in range(P)]
position_sumsC = [0 for i in range(P)]
position_sumsG = [0 for i in range(P)]
position_sumsT = [0 for i in range(P)]
for i in range(len(posteriors)):
	for j in range(len(posteriors[i])-P+1):
		for p in range(P):
			base = unencode(seq_mx[i][j+p])
			if base == 0:
				position_sumsA[p] += posteriors[i][j+p]
			if base == 1:
				position_sumsC[p] += posteriors[i][j+p]
			if base == 2:
				position_sumsG[p] += posteriors[i][j+p]
			if base == 3:
				position_sumsT[p] += posteriors[i][j+p]

psi1_hats = [[] for i in range(4)]
Alist1 = []
for i in range(len(position_sumsA)):
	Alist1.append(position_sumsA[i]/len(posteriors))
Clist1 = []
for i in range(len(position_sumsC)):
	Clist1.append(position_sumsC[i]/len(posteriors))
Glist1 = []
for i in range(len(position_sumsG)):
	Glist1.append(position_sumsG[i]/len(posteriors))	
Tlist1 = []
for i in range(len(position_sumsT)):
	Tlist1.append(position_sumsT[i]/len(posteriors))
psi1_hats[0] = Alist1
psi1_hats[1] = Clist1
psi1_hats[2] = Glist1
psi1_hats[3] = Tlist1

#####################################################
psi0_hats = [[] for i in range(4)]
Alist0 = []
for i in range(len(position_sumsA)):
	Alist0.append((1-position_sumsA[i])/((len(posteriors[0])- \
		P+1-1)*len(posteriors)))
Clist0 = []
for i in range(len(position_sumsC)):
	Clist0.append((1-position_sumsC[i])/((len(posteriors[0])- \
		P+1-1)*len(posteriors)))
Glist0 = []
for i in range(len(position_sumsG)):
	Glist0.append((1-position_sumsG[i])/((len(posteriors[0])- \
		P+1-1)*len(posteriors)))	
Tlist0 = []
for i in range(len(position_sumsT)):
	Tlist0.append((1-position_sumsT[i])/((len(posteriors[0])- \
		P+1-1)*len(posteriors)))
psi0_hats[0] = Alist0
psi0_hats[1] = Clist0
psi0_hats[2] = Glist0
psi0_hats[3] = Tlist0

##############################################
# log likelihood

term1 = 0
for i in range(len(posteriors)):
	for j in range(len(posteriors[i])-P+1):
		term1 += posteriors[i][j] * math.log(lambdaj_hats[j])

term2 = 0 
for i in range(len(posteriors)):
	for j in range(len(posteriors[i])-P+1):
		for p in range(P):
			base = unencode(seq_mx[i][j+p])
			this = (posteriors[i][j])*(math.log(psi1_hats[base][p]))
			print(psi0_hats[base][p])
			that = (1-posteriors[i][j])*(math.log(psi0_hats[base][p]))
			term2 += this + that
	









				


























	
