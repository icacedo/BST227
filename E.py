import random

psi0 = [None] * 4
# returns a vector of a frequency distribution
def freq_dist(vector):

	isone = 0
	while (isone != 1):
		total = 0 
		isone = 0
		for i in range(len(vector)):
			freq = float(random.uniform(0, 1))
			vector[i] = freq
			total += freq
		for i in range(len(vector)):
			vector[i] = round(vector[i]/total, 3)
		for i in range(len(vector)):
			isone += vector[i]
		if isone == 1: 
			return vector
			break
print(freq_dist(psi0))
	
