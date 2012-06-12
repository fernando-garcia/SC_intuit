import sys
import os
from Motif import *
import utils

def helpProgram():
	print "\nCorrect usage:\n\n\tpython intuit.py -fMotif <fMotif> -fSeq <FileSeq>"
	print "\npython intuit -h prints this help\n"

fMotif = None
fSeq = None
out = None
verbose = False

i = 1
while i < len(sys.argv):
	if sys.argv[i] == "-fMotif":
		i += 1
		fMotif = sys.argv[i]
	elif sys.argv[i] == "-fSeq":
		i += 1
		fSeq = sys.argv[i]
	elif sys.argv[i] == "-h":
		helpProgram()
		sys.exit(-1)
	else:
		print "Unknown flag " + sys.argv[i]
		helpProgram()
		sys.exit(-1)
	i+=1
	
if (fSeq != None and fMotif == None):
	print("-fMotif not found")
	helpProgram()
	sys.exit(-1)

if (fSeq == None and fMotif != None):
	print("-fSeq not found")
	helpProgram()
	sys.exit(-1)


if (fSeq == None and fMotif == None):
	print("Input files not found")
	helpProgram()
	sys.exit(-1)

if fMotif != None and not os.path.isfile(fMotif):
	print "Input file fMotif " + fMotif + " not found"
	sys.exit(-1)

if fSeq != None and not os.path.isfile(fSeq):
	print "Input file fSeq " + fSeq + " not found"
	sys.exit(-1)



if fMotif != None and fSeq != None:
	m = Motif(fileIn = fMotif)
	s = utils.readSeq(fileIn = fSeq)
	
	print m.SC_intuit(s)

else:
	helpProgram()
	sys.exit(-1)

