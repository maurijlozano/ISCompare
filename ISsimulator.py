#!/usr/bin/env python3
'''
ISsimulator is a program designed to generate random insertions of a desired IS into a provided genome.
Version:0.1
Written by Mauricio J. Lozano
UNLP - CONICET - Instituto de Biotecnología y Biología Molecular (IBBM)
'''
VERSION="0.1"
REF="None yet"
GITHUB="https://github.com/maurijlozano/ISMapper"

#modules
import os
import argparse
import random

#third party
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

#argument parsing
def parseArgs():
	'''
	Argument parsing is done.
	'''
	parser = argparse.ArgumentParser(description='ISsimulator is a program designed to randomly insert a desired IS into a provided genome.')
	parser.add_argument("-i", "--input",help="Genome in genbank format...", dest="input", action='store', required=True)
	parser.add_argument("-s", "--IS",help="Insertion sequence.", dest="IS", action='store', required=True)
	parser.add_argument("-o", "--Output",help="Output Genbank file name (e.g. name.gb).",action="store", dest="output", required=False)
	parser.add_argument("-dr", "--directRepeat",help="Direct repeat length [Default 0].",action="store", dest="directRepeats", required=False)
	parser.add_argument("-n", "--nisertions",help="Number of Is insertions [Default 30].", action='store', dest="ninsertions", required=False)
	args = parser.parse_args()
	if not args.ninsertions:
		args.ninsertions = 30
	if not args.output:
		args.output = str(os.path.splitext(args.input)[0]) + "_" +str(args.ninsertions) + "_IS.gb"
	return args

#Print presentation

def printSoftName():
	print("\n\n\n")
	print("   ********************************************")
	print("   *****   ISsimulator")
	print("   *****   Version: "+str(VERSION))
	print("   *****   Developed by Mauricio J. Lozano")
	print("   *****   github.com/maurijlozano")
	print("   ********************************************")
	print("   Please cite: "+REF)
	print("   Downloaded from: "+GITHUB)
	print("\n\n\n")

#********************************************************************************************************
# FUNCTIONS
def annotateIS(IS):
	featstart = SeqFeature.ExactPosition(0)
	featend = SeqFeature.ExactPosition(len(IS.seq))
	featLoc = FeatureLocation(featstart,featend,strand=+1)
	featType = SeqFeature.SeqFeature(featLoc, "mobile_element", qualifiers={"mobile_element_type":"artificially inserted IS"})
	IS.features.append(featType)
	pass

def generateRandomInsertion(record,IS,insertionCoord):
	start = 0
	end = len(record.seq)
	recordLeft = record[start:insertionCoord]
	dr = record[insertionCoord-directRepeats:insertionCoord]
	recordRight = record[insertionCoord:end]
	strand = random.choice(["-1","+1"])
	if strand == "-1":
		ISid=IS.id
		ISname = IS.name
		ISdescription = IS.description
		IS = IS.reverse_complement()
		IS.id=ISid
		IS.name = ISname
		IS.description = ISdescription
	newrecord = recordLeft+IS+dr+recordRight
	return newrecord

def getRecordInfo(record):
	record.id = str(record.id)+"_art_ISs"
	record.name =str(record.name)+"_art_ISs"
	record.description = str(record.description)+" with artificially inserted ISs"
	return [record.id,record.name,record.description]

def setRecordInfo(record,recordInfo):
	record.id = recordInfo[0]
	record.name = recordInfo[1]
	record.description =recordInfo[2]
	pass

#********************************************************************************************************
#********************************************************************************************************
#********************************************************************************************************
# Main
if __name__ == "__main__":
	#Presentation
	printSoftName()
	#Argument parsing
	args = parseArgs()
	input = args.input
	output = args.output
	outputTable = str(os.path.splitext(output)[0])+".csv"
	outputfasta = str(os.path.splitext(output)[0])+".fa"
	ninsertions = int(args.ninsertions)
	if not args.directRepeats:
		directRepeats = 0
	isSeq = args.IS
	#Load IS and Genome
	records = list(SeqIO.parse(input, "gb"))
	IS = list(SeqIO.parse(isSeq, "fasta"))[0]
	annotateIS(IS)
	#insert record into GB file!! record slices can be added with all their features...!!!! record[start:end] + record...
	f = open(outputTable,'w+')
	# if len(records) == 1:
	# 	record = records[0]
	# 	recInfo = getRecordInfo(record)
	# 	insertionCoords = []
	# 	for i in range(0,ninsertions):
	# 		insertionCoords.append(random.randint(0,len(record.seq)))
	# 	insertionCoords.sort()
	# 	for insertionCoord in insertionCoords:
	# 		record = generateRandomInsertion(record,IS,insertionCoord)
	# 		f.write(str(i)+','+str(recInfo[0])+','+str(insertionCoord)+"\n")
	# 	record.seq.alphabet = IUPAC.IUPACUnambiguousDNA()
	# 	setRecordInfo(record,recInfo)
	# 	SeqIO.write(record,output,'gb')
	# 	SeqIO.write(record,outputfasta,'fasta')
	# # The insertions must be randomly distributed in all records... derive a probability depending on record length.
	# else:
	recordIDXs = list(range(0,len(records)))
	recordLengths = [len(record.seq) for record in records]
	totalLength = sum(recordLengths)
	weights = [recordLength/totalLength for recordLength in recordLengths]
	insertionsperrecord = {}
	for rec in recordIDXs:
		insertionsperrecord[rec] = 0
	for i in range(0,ninsertions):
		randomRecord = random.choices(recordIDXs,weights=weights,k=1)[0]
		insertionsperrecord[randomRecord] += 1
	for rec in range(0,len(records)):
		recInfo = getRecordInfo(records[rec])
		insertionCoords = []
		for i in range(0,insertionsperrecord[rec]):
			insertionCoords.append(random.randint(0,len(records[rec].seq)))
		insertionCoords.sort()
		for insertionCoord in insertionCoords:
			records[rec] = generateRandomInsertion(records[rec],IS,insertionCoord)
			f.write(str(recInfo[0])+','+str(insertionCoord)+"\n")
		records[rec].seq.alphabet = IUPAC.IUPACUnambiguousDNA()
		setRecordInfo(records[rec],recInfo)
	SeqIO.write(records,output,'gb')
	SeqIO.write(records,outputfasta,'fasta')
	f.close()
	print('ISs inserted!!!')
