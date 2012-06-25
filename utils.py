
def EOF(f):
	pos = f.tell()					#f.tell() returns the current position of the file read/write pointer within the file (Integer Bytes).
	if(f.readline()):
		f.seek(pos)					#moves #pos bytes ahead
		return False

	else:
		return True	

def decideFSMorNot(mat):
	for i in range(len(mat)):
		#print sum(mat[i,:])
		if sum(mat[i,:]) > 1.5: return True
	return False
	
def readSeq(fileIn):									#seqs-size LIMITATION to one-line
	f = open(fileIn)
	#return f.readline().upper().replace('\n','')
	myseqs = []
	for line in f:
		if (line[0] !=	'>'):
			myseqs.append(line.upper().replace('\n',''))
	return myseqs

def computeMaxMin(A,row):
	_max = -1
	_min = 2
	maxV =  max(A[1][row])
	
	for i in range(16): 
		mu = A[0][row][i]
		v = maxV - A[1][row][i]
		curSim = mu * v
		if curSim > _max: _max = curSim
		if curSim < _min: _min = curSim
	
	return _max,_min
