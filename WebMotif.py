from Motif import *

class WebMotif(Motif):

	def SC_intuit_Web(self,seq):
		num_operations = len(seq) - self.n + 1
		results = {}
		inicont = 0
		endcont = self.n
		for i in range(num_operations):
			subseq = seq[inicont:endcont]
			print subseq, len(subseq), inicont, endcont
			results[self.SC_intuit(subseq)] = subseq
			
			inicont += 1
			endcont += 1

		keys = results.keys()
		keys.sort()
		return [dict[key] for key in keys]
		#return results.items().sort()
