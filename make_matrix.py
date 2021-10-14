'''
Code that takes in a file of sequences and produces the matrix X_{ijk} that was
described in class. Assume the text file input is a set of N sequences (one per
line, all the same length, with characters A, C, G and T). 
'''

import math
import random
import sys

def random_dna(N, length, A, C, G, T):
	
	assert(math.isclose(1, A+C+G+T, abs_tol = 0.001))	
	seqs = []
	for i in range(N):
		seq = []
		for i in range(length):
			r = random.random()
			if r < A:		seq.append('A')
			elif r < A + C:	seq.append('C')
			elif r < A + C +G:	seq.append('G')
			else:			seq.append('T')
		seqs.append(''.join(seq))
	return '\n'.join(seqs)

'''
file = open('sum_seqs.txt', 'w')
file.write(random_dna(10, 10, 0.25, 0.25, 0.25, 0.25))
file.close
'''

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

#print(make_matrix(sys.argv[1]))




























