import utils
import types
import sys
from numpy import *


class Motif:
	
	def __init__(self, fileIn):
		
		self.name = None
		self.ID = None
		self.FSM = None
		self.n = None
		self.numSamples = None
		self.matFSM = None
		self.matPSSM = None
		self.seq = []
		
		if(fileIn != None):
			self.readFromFile(fileIn)
		else:
			print "Wrong usage for Motif.__init__()"

		self.intuitM = zeros((16,(self.n * (self.n -1))/2))	# Matrix for the membership degree "\mu"
		self.intuitV = zeros((16,(self.n * (self.n -1))/2))	# Matrix for the non-membership degree "\nu"
	

	
	def calculatePSSMFromFSM(self):
		self.matPSSM = array(self.matFSM,float)
		self.matPSSM /= self.numSamples

	def calculateFSMFromPSSM(self,numSamples=None):
		if numSamples == None:
			numSamples = self.numSamples

		self.numSamples = numSamples
		self.matFSM = array(self.matPSSM*self.numSamples).astype(int)
		for i in range(self.n):
			s = sum(self.matFSM[i,:])
			ind = self.matFSM[i,:].tolist().index(max(self.matFSM[i,:]))
			self.matFSM[i,ind] += self.numSamples - s
		


	def pseudoCountFSM(self):
		if not all(self.matFSM):
			self.matFSM += 1
			self.numSamples += 4
			self.calculatePSSMFromFSM()

	def pseudoCountPSSM(self):
		if not all(self.matPSSM):
			
			self.matPSSM += 0.000001
			self.calculateFSMFromPSSM()

						 

	def ICmatrix(self):
	
		self.pseudoCountPSSM()
		matAux = self.matPSSM * (log(self.matPSSM)/log(2))		
		return matAux
	
	

				
	def readFromFile(self,fileIn):

		if type(fileIn) == types.StringType: fileIn = open(fileIn)
		
		header = fileIn.readline().split()
		AC = header[0].replace('>','')
		name = AC
		fileIn.readline()
		mat = []
		col = fileIn.readline().split()
		pos = fileIn.tell()
		while(col and col[0] != 'SEQUENCES'):
			mat.append([float(s) for s in col])
			pos = fileIn.tell()
			col = fileIn.readline().split()
				
		col = fileIn.readline().upper().replace('\n','')
		
		while col and not utils.EOF(fileIn) and col[0] != '>':
			self.seq.append(col)
			col = fileIn.readline().upper().replace('\n','')
			
				
		mat = array(mat,float)
	
		
		if self.FSM == None:
			self.FSM = utils.decideFSMorNot(mat)
		
		if self.FSM == True:
			#self.numSamples = sum(mat[0,:]) # 26/8/08 CAMBIADO xq el numSamples seria la maxima suma de las posiciones del motivo
			self.numSamples = self.computeNumSamples(mat) # 26/8/08 Esto es lo nuevo!
			self.matFSM = mat.astype(int)
			self.n = mat.shape[0]
			self.calculatePSSMFromFSM()
			self.name = name
			self.ID = AC
		elif self.FSM == False:
			self.matPSSM = array(mat,float)
			self.n = mat.shape[0]
			self.calculateFSMFromPSSM()
			self.name = name
			self.ID = AC
		
			

	def computeNumSamples(self,mat):
		numSamples = 0
		for i in range(len(mat)):
			if numSamples < sum(mat[i,:]):
				numSamples = sum(mat[i,:])
		return numSamples

	
	def computeIntuitionisticMotif(self):
		t = float(len(self.seq))
		F = zeros((16,(self.n * (self.n -1))/2))
		a = 0.000001 * 0.000001
	

		dposmatFSM = {}
		dposmatFSM['A'] = 0
		dposmatFSM['C'] = 1
		dposmatFSM['G'] = 2
		dposmatFSM['T'] = 3

		drow = {}
		drevrow = {}
	
		conti = 0
		for i in ['A','C','G','T']:
			for j in ['A','C','G','T']:
				drow[str(i)+str(j)] = conti
				drevrow[conti] = str(i)+str(j)
				conti += 1

		dcol = {}
		drevcol = {}
		conti = 0
		for i in range(self.n):
			for j in range(i+1,self.n):
				dcol[(i,j)] = conti
				drevcol[conti] = (i,j)
				conti += 1
	
		
		for s in self.seq:
			for indexi,posi in enumerate(s[0:-1]):
				for indexj,posj in enumerate(s[indexi+1:]):
					F[drow[str(posi)+str(posj)],dcol[(indexi,indexj+indexi+1)]] += 1
			
		self.intuitM = (F /t) + a
		
		for indexi in range(self.n):
			for indexj in range(indexi+1,self.n):
				for posi in ['A','C','G','T']:
					for posj in ['A','C','G','T']:
						izq = self.matPSSM[indexi][dposmatFSM[str(posi)]]
						dcha = self.matPSSM[indexj][dposmatFSM[str(posj)]]
						tot1 = (izq + dcha) / 2.
						oldM = self.intuitM[drow[str(posi)+str(posj)],dcol[(indexi,indexj)]]
						self.intuitM[drow[str(posi)+str(posj)],dcol[(indexi,indexj)]] += (1-oldM)*tot1						
				
		ic = self.ICmatrix()
		
		for indexi in range(self.n):
			for indexj in range(indexi+1,self.n):
				for posi in ['A','C','G','T']:
					for posj in ['A','C','G','T']:
						izq = (2 + ic[indexi,dposmatFSM[str(posi)]]) / 2
						dcha = (2 + ic[indexj, dposmatFSM[str(posj)]]) / 2
						tot1 = (izq + dcha) / 2
						tot2 = tot1 * (1- self.intuitM[drow[str(posi)+str(posj)],dcol[(indexi,indexj)]])
						self.intuitV[drow[str(posi)+str(posj)],dcol[(indexi,indexj)]] = tot2
		
	def SC_intuit(self,seq):

		self.computeIntuitionisticMotif()
		
		drow = {}
		drevrow = {}
		
		conti = 0
		for i in ['A','C','G','T']:
			for j in ['A','C','G','T']:
				drow[str(i)+str(j)] = conti
				drevrow[conti] = str(i)+str(j)
				conti += 1
	
		dcol = {}
		drevcol = {}
		conti = 0
		for i in range(self.n):
			for j in range(i+1,self.n):
				dcol[(i,j)] = conti
				drevcol[conti] = (i,j)
				conti += 1
	
		
		A = [self.intuitM.transpose(),self.intuitV.transpose()]
		
		
		simsNorm = []
		for indexi,posi in enumerate(seq[0:-1]):
			for indexj,posj in enumerate(seq[indexi+1:]):
				maxSim,minSim = utils.computeMaxMin(A,dcol[(indexi,indexj+ indexi +1)])
				mu = A[0][dcol[(indexi,indexj+ indexi +1)]][drow[posi+posj]]
				maxV =  max(A[1][dcol[(indexi,indexj+ indexi +1)]])
				v = maxV - A[1][dcol[(indexi,indexj+ indexi +1)]][drow[posi+posj]]
				curSim = mu * v
				curSimNorm = (curSim - minSim) / (maxSim -minSim)
				simsNorm.append(curSimNorm)
		
		return average(simsNorm)
