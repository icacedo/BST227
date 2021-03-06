import sys
import random
import math
import numpy as np

##### read in sequences #######################################################
fp = open(sys.argv[1], 'r')

seqs = []
for line in fp.readlines():
	line = line.rstrip()
	seqs.append(line)
	
##### functions used to initialize psi ########################################
def freq_dist(length, decimal):

	vector = [None] * length
	isone = 0
	while (isone != 1):
		total = 0 
		isone = 0
		for i in range(len(vector)):
			# don't want 0s in the initializations
			freq = float(random.uniform(0.01, 1))
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

##### initialize parameters ###################################################
P = 6
psi_0s = psi(P, 2)
psi_1s = psi(P, 2)
lambdajs = []
for i in range(len(seqs[0])):
	# don't want 0s in the initializations
	lambdajs.append(round(random.uniform(0.01, 1), 2))
	
##### E step ##################################################################
def E_step(seqs, P, psi_0s, psi_1s, lambdajs):

	letters = {'A':0, 'C':0, 'G':0, 'T':0}
	numerators = []
	# position j is where the motif starts
	# i for each sequence
	for i in range(len(seqs)):
		# at each position, motif can start at up to L-P+1 positions
		Cijs_list = []
		for j in range(len(seqs[0])-(P-1)):
			# probability that Cij is the start of a motif is lambdaj
			# times the sum of probabilities for psi
			Cij = math.log(lambdajs[j]) 
			for p in range(P):
				Xijp = letters[seqs[i][j+p]]
				# logsumexp trick
				# change from *= to +=
				Cij += math.log(psi_1s[Xijp][p]) 
			for jp in range(len(seqs[0])-(P-1)):
				# jp is j'
				# to get background frequencies
				# they cannot come from the same j as the foreground
				if jp == j: continue
				for p in range(P):	
					Xijp = letters[seqs[i][j+p]]
					# logsumexp
					Cij += math.log(psi_0s[Xijp][p])
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
		
	# add back in zeros to keep L at 38
	# can't actually be zero
	for i in range(len(posteriors)):
	# first position in motif length should not be 0
		for p in range(P-1):
			posteriors[i].append(0.0001)
		
	return posteriors
	
##### M step ##################################################################
posteriors = E_step(seqs, P, psi_0s, psi_1s, lambdajs)
theta = {'lmbda': lambdajs, 'psi_0': psi_0s, 'psi_1': psi_1s}

def M_step(seqs, posteriors, P):

	theta_hat = {'lmbda_hat': None, 'psi_1_hat': None, 'psi_0_hat': None}

	arr_lmbda = np.array([np.array(row) for row in posteriors])
	arrsum_lmbda = np.sum(arr_lmbda, axis=0)
	theta_hat['lmbda_hat'] = np.divide(arrsum_lmbda, len(posteriors))

	letters = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
	psi1_arr = np.zeros(shape=(4,P))
	for i in range(len(posteriors)):
		for j in range(len(posteriors[0])-P+1):
			for p in range(P):
				Xijp = letters[seqs[i][j+p]]
				k = letters[seqs[i][j]]
				if Xijp == k:
					Xijpk = 1
				else:
					Xijpk = 0
				psi1_arr[Xijp][p] += Xijpk * posteriors[i][j]
				
	theta_hat['psi_1_hat'] = np.divide(psi1_arr, len(posteriors))

	psi0_arr = np.zeros(shape=(4,P))		
	for i in range(len(posteriors)):
		for j in range(len(posteriors[0])-P+1):
			for p in range(P):
				Xijp = letters[seqs[i][j+p]]
				k = letters[seqs[i][j]]
				if Xijp == k:
					Xijpk = 1
				else:
					Xijpk = 0
				psi0_arr[Xijp][p] += Xijpk * (1- \
					posteriors[i][j])

	theta_hat['psi_0_hat'] = np.divide(psi0_arr, ((len(posteriors[0]) \
		-P+1-1)*(len(posteriors))))
	
	return theta_hat


##### test ####################################################################
# wanted to see if our e steps get the same-ish output
# they are the same-ish
np.random.seed(10)

def init_EM(seq_length, motif_length):
    lmbda = np.random.uniform(0,1,size=(seq_length,))
    lmbda = lmbda/np.sum(lmbda)  # normalization
    psi_0 = np.random.uniform(0,1,size=(4,motif_length))
    psi_0 = psi_0/psi_0.sum(axis=0)
    psi_1 = np.random.uniform(0,1,size=(4,motif_length))
    psi_1 = psi_1/psi_1.sum(axis=0)
    theta = {'lmbda': lmbda, 'psi_0': psi_0, 'psi_1': psi_1}
    return theta
    
# katherine's e step
def E_step_kat(data, theta, P):
    dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    C = []
    for i in range(len(data)):
        C_i = []
        for j in range(len(data[0])-P+1):   # 0 to 38-6+1
            C_ij = np.log(theta['lmbda'][j])
            # Iterate through all positions of the motif
            for p in range(P):
                base = data[i][j+p]
                k = dict[base]
                C_ij += np.log(theta['psi_0'][k][p])
            # Iterate through all positions of the non-motif
            for jpr in range(len(data[0])-P+1): # j' is the start position of a non-motif sequence
                if jpr == j: # if j:j+p includes a base that is non motif, score it as background
                    continue
                for p in range(P):
                    base = data[i][jpr+p]
                    k = dict[base]
                    C_ij += np.log(theta['psi_0'][k][p])
            C_i.append(np.exp(C_ij))  # move cij back to probability space
        sm = sum(C_i) # denominator
        C_i = [item/sm for item in C_i]  # normalization
        C.append(C_i)
    return C

theta = init_EM(len(seqs[0]), P)

#C = E_step_kat(seqs, theta, P)

posteriors = E_step(seqs, P, theta['psi_0'], theta['psi_1'], theta['lmbda'])

#print(posteriors)

# testing what np.sum and np.divide does
'''
arr = np.array([[1,2,3],[4,5,6]])
print(arr)
arrsum = np.sum(arr, axis=0)
print(arrsum)
arrdiv = np.divide(arrsum, 2)
print(arrdiv)
'''
# katherine's m step
def M_step_kat(data, C, P):
    dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    lmbda = np.array(C).sum(axis=0)  # sum by column of matrix
    lmbda = lmbda / 357  # divide all elements in list by N (normalization)

    # Initialize new psi matrices
    psi_1 = np.zeros((4, P))
    psi_0 = np.zeros((4, P))
    for p in range(P):
        for i in range(len(data)):
            for j in range(0, len(data[0])-P+1):
                base = data[i][j+p]
                k = dict[base]
                psi_1[k, p] += C[i][j]
                psi_0[k, p] += 1 - (C[i][j])
    psi_1 /= len(data)  # normalization
    psi_0 /= len(data)*(len(data[0])-P)  # normalization
    theta = {'lmbda': lmbda, 'psi_1': psi_1, 'psi_0': psi_0}
    return theta

np.random.seed(10)
theta_hat_kat = M_step_kat(seqs, posteriors, P)
theta_hat = M_step(seqs, posteriors, P)

print(theta_hat)




























