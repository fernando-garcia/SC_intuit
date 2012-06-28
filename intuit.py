import sys
import os
from Motif import *
from WebMotif import *
import utils

def helpProgram():
	print "\nCorrect usage:\n\n\tpython intuit.py -fMotif <fMotif> -fSeq <FileSeq>"
	print "\npython intuit -h prints this help\n"

fMotif = None
fSeq = None
out = None
verbose = False
TRANSFACpath = '/home/pulido/Documents/TRANSFAC/matrix.dat'

def ReadTransfacFile(pathfile)
	
	ID_BEGINING = 'AC'
	ID_END = '//'
	OUTPUT_FILENAME = '.transafc.txt'
	tempath = '/home/pulido/Documents/TRANSFAC/matrix_temp/'

	source_f = open(pathfile, 'r')
	transfile = ""
	canprint = False

	for line in source_f:

		if line[0:2] == ID_BEGINING:
			ofm = tempath+line.replace('\n','').split('  ')[1]+OUTPUT_FILENAME
			fh = open(ofm, 'w')
			canprint = True

		if canprint:
			fh.write(line)
			#transfile = transfile + line
		if (line[0:2] == ID_END and canprint and fh):
			fh.close()
			canprint = False



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
	m = WebMotif(fileIn = fMotif)
	seqs = utils.readSeq(fileIn = fSeq)
	
	print m.SC_intuit_Web(seqs)
else:
	helpProgram()
	sys.exit(-1)

