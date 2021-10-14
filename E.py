import random
import make_matrix as mm
import sys

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
			
def psis(length, decimal):
	psi = []
	for i in range(2):
		psi.append(freq_dist(length, decimal))
	return psi[0], psi[1]

def lambdas(decimal):
	lambda0 = round(random.uniform(0, 1), decimal)
	lambda1 = round(1 - lambda0, decimal)
	return lambda0, lambda1

# do these need to be in a loop?
# to calculate posteriors for different values of lambda and psi?
# only need to calculate posterior 4 times, but for different lambda and psi values?
tup_l = lambdas(3)
lambda0 = tup_l[0]
lambda1 = tup_l[1]
print(lambda0, lambda1)

tup_p = psis(4, 3)
psi0 = tup_p[0]
psi1 = tup_p[1]
print(psi0, psi1)

# L is the total number of a particular nucleotide that was seen?
# if L = 0, nothing was seen, what does L = -1 mean? L = infinity?
# but L can only be 0 or 1?

'''	
for sequence in mm.make_matrix(sys.argv[1]):
	for position in seq:
		print(position)
'''
	
	
	
	
	
	
