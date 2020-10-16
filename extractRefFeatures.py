#!/usr/bin/env python3
'''
ExtractRefFeatures.py generates a features table csv for ISseeker.
Written by Mauricio J. Lozano
UNLP - CONICET - Instituto de Biotecnología y Biología Molecular (IBBM)
'''
VERSION="1"
REF="None yet"
GITHUB="https://github.com/maurijlozano/ISMapper"

#modules
from Bio import SeqIO
import argparse, os, re

#argument parsing
def parseArgs():
	'''
	Argument parsing is done.
	'''
	parser = argparse.ArgumentParser(description='ExtractRefFeatures is a program designed to extract genome features to use with ISseeker..')
	parser.add_argument("-i", "--input",help="Input file in genbank format...", dest="file", action='store', required=True)
	args = parser.parse_args()
	return args

#Print presentation

def printSoftName():
	print("\n\n\n")
	print("   ********************************************")
	print("   *****   ExtractRefFeatures")
	print("   *****   Version: "+str(VERSION))
	print("   *****   Developed by Mauricio J. Lozano")
	print("   *****   github.com/maurijlozano")
	print("   ********************************************")
	print("   Please cite: "+REF)
	print("   Downloaded from: "+GITHUB)
	print("\n\n\n")

#********************************************************************************************************
# FUNCTIONS

def extractFeaturesFromGB(file):
	bin = 0
	for record in SeqIO.parse(file, "genbank"):
		fileName = record.id
		f = open(fileName+'.csv', "w+")
		for feature in record.features:
			if feature.type == "CDS":
				if feature.location.strand == 1:
					strand = "+"
				elif feature.location.strand == -1: 
					strand = "-"
				else:
					strand = "unknown"
				f.write(str(feature.location.start) + "," + str(feature.location.end) + "," + str(bin) + ","  + strand + "," + str(feature.qualifiers.get('locus_tag')[0]) + ","  + str(re.sub(',',' ',feature.qualifiers.get('product')[0])) +"\n")
		f.close()


#********************************************************************************************************
#********************************************************************************************************
#********************************************************************************************************
# Main
if __name__ == "__main__":
	#Presentation
	printSoftName()
	#Argument parsing
	args = parseArgs()
	file = args.file
	extractFeaturesFromGB(file)
	print('Done!')

