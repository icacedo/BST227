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

def lambdas(decimal):
	lambda0 = round(random.uniform(0, 1), decimal)
	lambda1 = round(1 - lambda0, decimal)
	return lambda0, lambda1
	
def shape(matrix):
	rows = 0
	for i in seq_mx:
		rows += 1
		columns = 0
		for j in i:
			columns += 1
	return (rows, columns)
	
# initialize psis
# matrix of parameters size 4 x P
# P = 6, length of motif
psi_0 = psi(6, 2)
psi_1 = psi(6, 2)

# initialize lambdas
# 0 = fg, 1 = bg
tup = lambdas(2)
lambda_0 = tup[0]
lambda_1 = tup[1]




# A = 0, C = 1, G = 2, T = 3
# positions 0 to 7


# read in data to matrix and get matrix shape
seq_mx = make_matrix(sys.argv[1])
seq_mx_shape = shape(seq_mx)
# rows aka sequences
rows = seq_mx_shape[0]
# columns aka positions
columns = seq_mx_shape[1]

# position j is where the motif starts
# for each sequence
for i in range(rows):
	# at each position
	for j in range(columns):
		print(j)

























	
