'''
Code that takes in a file of sequences and produces the matrix X_{ijk} that was
described in class. Assume the text file input is a set of N sequences (one per
line, all the same length, with characters A, C, G and T). 
'''

import math
import random

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

fp = open('sum_seqs.txt', 'r')

n = 10
N = 10
'''
matrix = [[None for i in range(n)] for j in range(N)]
#Note: matrix[row][column]

i = 0
for line in fp.readlines():
	line = line.rstrip()
	j = 0
	for n in line:
		if n.upper() == 'A':
			matrix[i][j] = 1
		if n.upper() == 'C':
			matrix[i][j] = 2
		if n.upper() == 'G':
			matrix[i][j] = 3
		if n.upper() == 'T':
			matrix[i][j] = 4
		j += 1
	i += 1
print(matrix)
'''






























