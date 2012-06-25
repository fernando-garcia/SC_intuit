from Motif import *

class WebMotif(Motif):

	def SC_intuit_Web(self,seqs):
		print seqs
		num_seqs = len(seqs)
		result = []
		resume = ""
		cont = 0

		for i in range(num_seqs):
			seq = seqs[cont]
			num_operations = len(seq) - self.n + 1
			dic = {}
			inicont = 0
			endcont = self.n

			for i in range(num_operations):
				subseq = seq[inicont:endcont]
				#print subseq, len(subseq), inicont, endcont
				dic[subseq] = self.SC_intuit(subseq)
			
				inicont += 1
				endcont += 1

			cont = cont+1
			best_Scored = sorted(dic, key=lambda key: dic[key], reverse = True)[0]
			resume = resume + "sequence: "+seq+'\t'+"best subsequence: "+ best_Scored+'\t'+'score: '+str(dic[best_Scored])+'\n'
			result.append( [seq, best_Scored, dic[best_Scored]] )

		return resume
