import os, sys, glob

ID_BEGINING = 'AC'
ID_END = '//'
OUTPUT_FILENAME = '.transafc.txt'
tempath = '/home/pulido/Documents/TRANSFAC/matrix_temp/'

source_f = open('/home/pulido/Documents/TRANSFAC/matrix.dat', 'r')
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

insideFiles = os.listdir(tempath)
insideFiles.sort()
