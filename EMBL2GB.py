#!/usr/bin/env python3
'''
EMBL to GB format conversion script.
Written by Mauricio J. Lozano
UNLP - CONICET - Instituto de Biotecnología y Biología Molecular (IBBM)
'''
VERSION="1.0"
REF="\n\n   Easy identification of insertion sequence mobilization events\n   in related bacterial strains with ISCompare. \n   E.G. Mogro, N. Ambrosis, M.J. Lozano\n   doi: https://doi.org/10.1101/2020.10.16.342287\n   Instituto de Biotecnología y Biología Molecular\n   CONICET - CCT La Plata - UNLP - FCE\n"
GITHUB="https://github.com/maurijlozano/ISCompare"

import os, argparse
#Third party modules
from Bio import SeqIO

#argument parsing
def parseArgs():
	'''
	Argument parsing is done.
	'''
	parser = argparse.ArgumentParser(description='EMBL2GB -> Convert genome files from EMBL to Genbank format.')
	parser.add_argument("-f", "--Files",help="EMBL files to convert", dest="files", action='append', nargs='+', required=True)
	#
	args = parser.parse_args()
	return args

#Print presentation

def printSoftName():
	print("\n\n\n")
	print("   ********************************************")
	print("   *****   EMBL2GB")
	print("   *****   Version: "+str(VERSION))
	print("   *****   Developed by Mauricio J. Lozano")
	print("   *****   github.com/maurijlozano")
	print("   ********************************************")
	print("   Please cite: "+REF)
	print("   Downloaded from: "+GITHUB)
	print("\n\n\n")

def checkFiles(filePath):
	if os.path.exists(filePath) and os.path.getsize(filePath) > 0:
		pass
	else:
		print('Error: '+ str(filePath) + ' file missing or empty.')
		exit()


# Main
if __name__ == "__main__":
	#Presentation
	printSoftName()
	#Argument parsing
	args = parseArgs()
	files=args.files[0]
	for filePath in files:
		checkFiles(filePath)
	for EMBLfile in files:
		basename = os.path.basename(EMBLfile)
		basename = os.path.splitext(basename)[0]
		SeqIO.convert(EMBLfile,'embl',basename+'.gb','gb')
		print(EMBLfile + ' converted to Genbank format.')
#end
