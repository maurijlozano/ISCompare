#!/usr/bin/env python3
'''
ExtractISbyName.py extracts IS sequences in fasta format from IS.fna file.
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
	parser = argparse.ArgumentParser(description='ExtractISbyName.py extracts IS sequences in fasta format from IS.fna file...')
	parser.add_argument("-i", "--input",help="Input the IS name from ISCompare Final results table...", dest="ISname", action='store', required=True)
	args = parser.parse_args()
	return args

#Print presentation

def printSoftName():
	print("\n\n\n")
	print("   ********************************************")
	print("   *****   ExtractISbyName")
	print("   *****   Version: "+str(VERSION))
	print("   *****   Developed by Mauricio J. Lozano")
	print("   *****   github.com/maurijlozano")
	print("   ********************************************")
	print("   Please cite: "+REF)
	print("   Downloaded from: "+GITHUB)
	print("\n\n\n")

#********************************************************************************************************
# FUNCTIONS

def extractISbyName(file,output,ISname):
	f = open(output, "w+")
	bin = 0
	for record in SeqIO.parse(file, "fasta"):
		if record.id == ISname:
			f.write(">"+str(record.id) + "\n"+ str(record.seq) +"\n")
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
	ISname = args.ISname
	cwd = os.getcwd()
	output = cwd+"/"+ISname+".fasta"
	file='IS.fna'
	extractISbyName(file,output,ISname)
	print('Done!')

