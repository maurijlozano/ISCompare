#!/usr/bin/env python3
'''
ISCompare is a program designed to look for and compare insertion sequence position between a query and a eference genome (or WGS Assembly).
Written by Mauricio J. Lozano
UNLP - CONICET - Instituto de Biotecnología y Biología Molecular (IBBM)
'''
VERSION="1.0.5"
REF="\n\n   Easy identification of insertion sequence mobilization events\n   in related bacterial strains with ISCompare. \n   E.G. Mogro, N. Ambrosis, M.J. Lozano\n   doi: https://doi.org/10.1101/2020.10.16.342287\n   Instituto de Biotecnología y Biología Molecular\n   CONICET - CCT La Plata - UNLP - FCE\n"
GITHUB="https://github.com/maurijlozano/ISCompare"

#Local modules
import ISFinderBlast
#Third party modules
from Bio import Entrez
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature as seqFeature
from os import sys
from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import re, os, subprocess, multiprocessing, shutil
import pandas as pd
import numpy as np
import errno

#argument parsing
def parseArgs():
	'''
	Argument parsing is done.
	'''
	parser = argparse.ArgumentParser(description='ISCompare is a program designed to look for and compare insertion sequence position between a query and a reference genomes (or WGS Assembly)..')
	parser.add_argument("-q", "--QueryFiles",help="Query genome in Genbank format...", dest="qfiles", action='append', nargs='+', required=False)
	parser.add_argument("-r", "--RefFiles",help="Reference genome in Genbank format...", dest="rfiles", action='append', nargs='+', required=False)
	parser.add_argument("-i", "--ISfile",help="IS database file [Multifasta nucleotide file]...", dest="isfile", action='store', required=False)
	parser.add_argument("-I", "--ISscan",help="Scan IS on query and reference genomes using ISFinder [Generates IS.fna file in results folder]...", dest="isscan", action='store_true', required=False)
	parser.add_argument("-Q", "--QueryACC",help="Accession numbers for the query genome to download from NCBI.", dest="qacc", action='append', nargs='+', required=False)
	parser.add_argument("-R", "--RefACC",help="Accession numbers for the reference genome to download from NCBI.", dest="racc", action='append', nargs='+', required=False)
	parser.add_argument("-e", "--email",help="User email. Required for accession number download mode.",action="store", dest="yourEmail", required=False)
	parser.add_argument("-E", "--evalue",help="E-value cutoff for blastn search.",action="store", dest="evalue", required=False)
	parser.add_argument("-o", "--OutputDir",help="Output folder.",action="store", dest="output", required=False)
	parser.add_argument("-x", "--IDScaffolds",help="OMIT the remove identical scaffolds steps.",action="store_false", dest="remove", required=False)
	parser.add_argument("-l", "--minAlnLength",help="Minimal required length of the alignment of the QIFs to the reference genome.", action="store", dest="minAlnLength", required=False)
	parser.add_argument("-m", "--minLength",help="Minimum QIFS length to be considered in the analysis. For complete genomes it should be < 2*surroundingLen. For genome assemblies with scaffolds/contigs it should be < surroundingLen.", action="store", dest="minLength", required=False)
	parser.add_argument("-s", "--surroundingLen",help="Nucleotides to extract from upstream and downstream IS blast hits [QIFs].",action="store", dest="surroundingLen", required=False)
	parser.add_argument("-s2", "--surroundingLen2",help="Nucleotides to extract from QIFs blast hits [RAFs].",action="store", dest="surroundingLen2", required=False)
	parser.add_argument("-S", "--shift",help="Shifts the start and end of the QIFs an specified number of nucleotides from the IS.",action="store", dest="shift", required=False)
	parser.add_argument("-d", "--ISdiff",help="ISdiff: minimal difference in nucloetids between qlen and IS alingnment length (discard IS scaffolds).",action="store", dest="ISdiff", required=False)
	parser.add_argument("-f", "--scaffoldDiff",help="scaffoldDiff --> maximal difference in nucloetides for two scaffolds to be considered identical.",action="store", dest="scaffoldDiff", required=False)
	parser.add_argument("-p", "--plot",help="Plot IS surroundings with genomic features in both query and reference genomes...", dest="plot", action='store_true', required=False)
	parser.add_argument("-c", "--clean",help="Clean files...", dest="clean", action='store_true', required=False)
	parser.add_argument("-rs", "--SLIS",help="Report IS with the same location (SLIS)...", dest="reportSLIS", action='store_true', required=False)
	#
	args = parser.parse_args()
	return args

#Print presentation

def printSoftName():
	print("\n\n\n")
	print("   ********************************************")
	print("   *****   ISCompare")
	print("   *****   Version: "+str(VERSION))
	print("   *****   Developed by Mauricio J. Lozano")
	print("   *****   github.com/maurijlozano")
	print("   ********************************************")
	print("   Please cite: "+REF)
	print("   Downloaded from: "+GITHUB)
	print("\n\n\n")

#********************************************************************************************************
# FUNCTIONS
def getRefSeqFromNCBIbyACC(acc,name):
	#Retrieve GB file for the reference genome using accession numbers and NCBItools.
	sequences = acc
	if outputDirProvided:
		fileName= os.path.join(outputDir,name+".gb")
		fileNameFasta= os.path.join(outputDir,name+".fasta")
	else:
		fileName="results/"+name+".gb"
		fileNameFasta="results/"+name+".fasta"
		# Create target Directory if don't exist
	if not os.path.isfile(fileName):
		f=open(fileName, "w+")
		for sequence in sequences:
			print("Downloading "+sequence+" reference sequence in Genbank format...")
			log.write("Downloading "+sequence+" reference sequence in Genbank format...\n")
			try:
				handle = Entrez.efetch(db="nucleotide", id=sequence, rettype="gbwithparts", retmode="text")
			except:
				print("There was an error retrieving the sequence from Genbank. Verify the accession number.. ")
				raise
			f.write(handle.read())
			handle.close()
		f.close()
	if not os.path.isfile(fileNameFasta):
		f=open(fileNameFasta, "w+")
		for sequence in sequences:
			print("Downloading "+sequence+" reference sequence in fasta format...")
			log.write("Downloading "+sequence+" reference sequence in fasta format...\n")
			try:
				handle = Entrez.efetch(db="nucleotide", id=sequence, rettype="fasta", retmode="text")
			except:
				print("There was an error retrieving the sequence from Genbank. Verify the accession number.. ")
				raise
			f.write(handle.read())
			handle.close()
		f.close()
	return name

def extractFastaFromGB(seq):
	f = open(os.path.join(outputDir,seq+".fasta"), "w+")
	fileBaseName = os.path.join(outputDir,seq)
	SeqIO.convert(fileBaseName+'.gb', "gb", f, "fasta")
	f.close()

def getNumberOfScaffoldsFromGB(qseq):
	#Get the number of scaffolds
	fileBaseName = os.path.join(outputDir , qseq)
	records = list(SeqIO.parse(fileBaseName + ".gb", "genbank"))
	maxlen = 0
	for record in records:
		if len(record) > maxlen:
			maxlen = len(record)
	return(len(records),maxlen)

def blastQuery2IS(query,ref,DBNAME,resFile):
	#Blastn function to search for ISs
	ncpu = multiprocessing.cpu_count()
	dbIn = ref
	dbOut = os.path.join(outputDir,"blastDB/"+DBNAME)
	query=query
	blastres= os.path.join(outputDir, resFile)
	subprocess.call(['makeblastdb', '-in',dbIn, '-out', dbOut, '-dbtype', 'nucl', '-hash_index' ], stdout=log, stderr=log, shell=False)
	subprocess.call(['blastn','-query', query, '-db', dbOut, '-out', blastres, '-evalue', evalue, '-word_size', '11','-num_threads', str(int(ncpu/2)), '-outfmt', '6 qseqid sseqid qlen slen length qstart qend sstart send evalue bitscore qcovs','-culling_limit','1'], stdout=log, stderr=log, shell=False)
	log.flush()
	pass

def blastQuery2IS2(query,ref,DBNAME,resFile):
	#Blastn function to search for ISs in RAFS
	ncpu=multiprocessing.cpu_count()
	dbIn=ref
	dbOut= os.path.join(outputDir,"blastDB/"+DBNAME)
	query=query
	blastres= os.path.join(outputDir, resFile)
	subprocess.call(['makeblastdb', '-in',dbIn, '-out', dbOut, '-dbtype', 'nucl', '-hash_index' ], stdout=log, stderr=log, shell=False)
	subprocess.call(['blastn','-query', query, '-db', dbOut, '-out', blastres, '-evalue', evalue,  '-word_size', '11', '-num_threads', str(int(ncpu/2)), '-outfmt', '6 qseqid sseqid qlen slen length qstart qend sstart send evalue bitscore qcovs'], stdout=log, stderr=log, shell=False)
	log.flush()
	pass

def blastQuery2Ref(query,ref,DBNAME,resFile):
	#Blastn function to search for QIFS in reference genome
	ncpu=multiprocessing.cpu_count()
	dbIn=ref
	dbOut= os.path.join(outputDir,"blastDB/"+DBNAME)
	query=query
	blastres= os.path.join(outputDir,resFile)
	subprocess.call(['makeblastdb', '-in',dbIn, '-out', dbOut, '-dbtype', 'nucl', '-hash_index' ], stdout=log, stderr=log, shell=False)
	subprocess.call(['blastn','-query', query, '-db', dbOut, '-out', blastres, '-evalue', evalue, '-word_size', '11','-num_threads', str(int(ncpu/2)), '-outfmt', '6 qseqid sseqid qlen slen length qstart qend sstart send evalue bitscore qcovs', '-max_target_seqs', '10','-max_hsps','5','-culling_limit','3'], stdout=log, stderr=log, shell=False)
	log.flush()
	pass

def blastQuery2Ref2(query,ref,DBNAME,resFile):
	#Blastn function to search for identical scaffolds
	ncpu=multiprocessing.cpu_count()
	dbIn=ref
	dbOut=os.path.join(outputDir,"blastDB/"+DBNAME)
	query=query
	blastres=os.path.join(outputDir,resFile)
	subprocess.call(['makeblastdb', '-in',dbIn, '-out', dbOut, '-dbtype', 'nucl', '-hash_index' ], stdout=log, stderr=log, shell=False)
	subprocess.call(['blastn','-query', query, '-db', dbOut, '-out', blastres, '-evalue', evalue, '-num_threads', str(int(ncpu/2)), '-outfmt', '6 qseqid sseqid qlen slen length qstart qend sstart send evalue bitscore qcovs qcovhsp', '-max_target_seqs', '20','-max_hsps','20', '-qcov_hsp_perc','10'], stdout=log, stderr=log, shell=False)
	log.flush()
	pass

def removeIdenticalScaffolds(fileName,identicalToRef):
	#Function to remove identical scaffolds from the analysis.
	records = list(SeqIO.parse(fileName, "fasta"))
	fileName2 = os.path.splitext(fileName)[0]+"_reduced.fasta"
	f = open(fileName, "w+")
	f2 = open(fileName2, "w+")
	for record in records:
		f.write(">"+record.id + "\n" + str(record.seq) + "\n")
		if record.id in list(identicalToRef["qseqid"]):
			continue
		else:
			f2.write(">"+record.id + "\n" + str(record.seq) + "\n")
	f.close()
	f2.close()
	pass

def findpairs(keep):
	#Function to find QIFS pairs, by location on the refgenome
	uniqueID=0
	modedBlastResTab = pd.DataFrame(columns=list(keep.columns))
	scaffoldNames = list(keep["qseqid"].unique())
	if len(scaffoldNames) > 0:
		for Scaffold in scaffoldNames:
			howManyScaffolds = len(keep[keep["qseqid"] == Scaffold])
			scaffoldEntries = pd.DataFrame(keep[keep["qseqid"] == Scaffold])
			scaffoldEntries['uniqueID'] = uniqueID
			scaffoldEntries['Reversed'] = ""
			scaffoldEntries['Extract'] = ""
			for idx,scaffold in scaffoldEntries.iterrows():
				scaffoldEntries.loc[idx,'sstart'] = min(scaffold['sstart'],scaffold['send'])
				scaffoldEntries.loc[idx,'send'] = max(scaffold['sstart'],scaffold['send'])
				#reversed sstart and send
				if scaffold['sstart'] < scaffold['send']:
					scaffoldEntries.loc[idx,'Reversed'] = False
				else:
					scaffoldEntries.loc[idx,'Reversed'] = True
				#only left
				if ((scaffold['Flank'] == "L") and not scaffoldEntries.loc[idx,'Reversed']):
					scaffoldEntries.loc[idx,'Extract'] = "R"
				elif ((scaffold['Flank'] == "L") and scaffoldEntries.loc[idx,'Reversed']):
					scaffoldEntries.loc[idx,'Extract'] = "L"
				elif ((scaffold['Flank'] == "R") and not scaffoldEntries.loc[idx,'Reversed']):
					scaffoldEntries.loc[idx,'Extract'] = "L"
				else:
					scaffoldEntries.loc[idx,'Extract'] = "R"
			if howManyScaffolds >= 2:
				newScaffoldEntries = pd.DataFrame(columns=scaffoldEntries.columns)
				newScaffoldEntries[['qlen', 'slen', 'length', 'qstart', 'qend','sstart', 'send', 'bitscore', 'start','end', 'Lend', 'Rstart', 'Llen', 'Rlen']] = newScaffoldEntries[['qlen', 'slen', 'length', 'qstart', 'qend','sstart', 'send', 'bitscore', 'start','end', 'Lend', 'Rstart', 'Llen', 'Rlen']].astype("Int32")
				for idx,scaffold in scaffoldEntries.iterrows():
					newScaffoldEntries['Reversed'] = newScaffoldEntries['Reversed'].astype('bool')
					newScaffoldEntries.loc[idx] = scaffold
				newScaffoldEntries = newScaffoldEntries.sort_values(["evalue"])
				ISMatchLen = int(re.sub(r'^.*\|ISMatchLen:([0-9]*)\*?\*?$','\\1',newScaffoldEntries.iloc[0,0]))
				Llen = newScaffoldEntries["Llen"].iloc[0] #.iloc[0,17]
				Rlen = newScaffoldEntries["Rlen"].iloc[0] #.iloc[0,18]
				#Remove multiple hits with only left or right flaks.
				if Llen > minLength and Rlen > minLength:
					left = 	newScaffoldEntries[newScaffoldEntries['Flank'] == "L"]
					right = newScaffoldEntries[newScaffoldEntries['Flank'] == "R"]
					if len(left) == 1:
						leftend = float(left['send'])
						leftstart = float(left['sstart'])
						right = right.sort_values(["evalue"])
						notAsigned=True
						for idx,entry in right.iterrows():
							rightstart = float(entry['sstart'])
							rightend = float(entry['send'])
							if  ((rightstart < leftend + ISMatchLen + 100 and rightstart > leftend - surroundingLen/2) or (rightend > leftstart - ISMatchLen - 100 and rightend < leftstart + surroundingLen/2)) and entry['sseqid'] == left.iloc[0,1]:
								modedBlastResTab = concatPdTables(modedBlastResTab,left)
								modedBlastResTab = concatPdTables(modedBlastResTab,entry)
								notAsigned=False
								break
						if notAsigned:
							for idx,entry in right.iterrows():
								rightstart = float(entry['sstart'])
								rightend = float(entry['send'])
								if  ((rightstart < leftend + ISMatchLen + 10000 and rightstart > leftend) or (rightend > leftstart - ISMatchLen - 10000 and rightend < leftstart)) and entry['sseqid'] == left.iloc[0,1]:
									modedBlastResTab = concatPdTables(modedBlastResTab,left)
									modedBlastResTab = concatPdTables(modedBlastResTab,entry)
									notAsigned=False
									break
						if notAsigned:
							for idx,entry in right.iterrows():
								rightstart = float(entry['sstart'])
								rightend = float(entry['send'])
								contigLen = entry['slen']
								if (entry['length'] > Rlen*.9) and ( rightstart < 1000  or rightend > contigLen - 1000 ):
									modedBlastResTab = concatPdTables(modedBlastResTab,left)
									modedBlastResTab = concatPdTables(modedBlastResTab,entry)
									notAsigned=False
									break
					elif len(right) == 1:
						rightstart = float(right['sstart'])
						rightend = float(right['send'])
						left = left.sort_values(["evalue"])
						notAsigned=True
						for idx,entry in left.iterrows():
							leftend = float(entry['send'])
							leftstart = float(entry['sstart'])
							if ((leftend > rightstart - ISMatchLen - 100 and leftend < rightstart + surroundingLen/2) or (leftstart < rightend + ISMatchLen + 100 and leftstart > rightend - surroundingLen/2) ) and entry['sseqid'] == right.iloc[0,1]:
								modedBlastResTab = concatPdTables(modedBlastResTab,right)
								modedBlastResTab = concatPdTables(modedBlastResTab,entry)
								notAsigned=False
								break
						if notAsigned:
							for idx,entry in left.iterrows():
								leftend = float(entry['send'])
								leftstart = float(entry['sstart'])
								if ((leftend > rightstart - ISMatchLen - 10000 and leftend < rightstart) or (leftstart < rightend + ISMatchLen + 10000 and leftstart > rightend)) and entry['sseqid'] == right.iloc[0,1]:
									modedBlastResTab = concatPdTables(modedBlastResTab,right)
									modedBlastResTab = concatPdTables(modedBlastResTab,entry)
									notAsigned=False
									break
						if notAsigned:
							for idx,entry in left.iterrows():
								leftend = float(entry['send'])
								leftstart = float(entry['sstart'])
								contigLen = entry['slen']
								if (entry['length'] > Llen*.9) and ( leftstart < 1000  or leftend > contigLen - 1000 ):
									modedBlastResTab = concatPdTables(modedBlastResTab,right)
									modedBlastResTab = concatPdTables(modedBlastResTab,entry)
									notAsigned=False
									break
					elif (len(left) > 1) and (len(right) > 1):
						left = left.sort_values(["evalue"])
						right = right.sort_values(["evalue"])
						notAsigned=True
						if (left.iloc[0,9] < left.iloc[1,9]):
							entryL = left.iloc[0,:]
							leftend = float(entryL['send'])
							leftstart = float(entryL['sstart'])
							right = right.sort_values(["evalue"])
							for idx,entry in right.iterrows():
								rightstart = float(entry['sstart'])
								rightend = float(entry['send'])
								if  ((rightstart < leftend + ISMatchLen + 100 and rightstart > leftend - surroundingLen/2) or (rightend > leftstart - ISMatchLen - 100 and rightend < leftstart + 50)) and entry['sseqid'] == left.iloc[0,1]:
									modedBlastResTab = concatPdTables(modedBlastResTab,entryL)
									modedBlastResTab = concatPdTables(modedBlastResTab,entry)
									notAsigned=False
									break
							if notAsigned:
								for idx,entry in right.iterrows():
									rightstart = float(entry['sstart'])
									rightend = float(entry['send'])
									if  ((rightstart < leftend + ISMatchLen + 10000 and rightstart > leftend) or (rightend > leftstart - ISMatchLen - 10000 and rightend < leftstart)) and entry['sseqid'] == left.iloc[0,1]:
										modedBlastResTab = concatPdTables(modedBlastResTab,entryL)
										modedBlastResTab = concatPdTables(modedBlastResTab,entry)
										notAsigned=False
										break
							if notAsigned:
								for idx,entry in right.iterrows():
									rightstart = float(entry['sstart'])
									rightend = float(entry['send'])
									contigLen = entry['slen']
									if (entry['length'] > Rlen*.9) and ( rightstart < 1000  or rightend > contigLen - 1000 ):
										modedBlastResTab = concatPdTables(modedBlastResTab,entryL)
										modedBlastResTab = concatPdTables(modedBlastResTab,entry)
										notAsigned=False
										break
						if notAsigned:
							if (right.iloc[0,9] < right.iloc[1,9]):
								entryR = right.iloc[0,:]
								rightstart = float(entryR['sstart'])
								rightend = float(entryR['send'])
								left = left.sort_values(["evalue"])
								for idx,entry in left.iterrows():
									leftend = float(entry['send'])
									leftstart = float(entry['sstart'])
									if ((leftend > rightstart - ISMatchLen - 100 and leftend < rightstart + surroundingLen/2) or (leftstart < rightend + ISMatchLen + 100 and leftstart > rightend - 50) ) and entry['sseqid'] == right.iloc[0,1]:
										modedBlastResTab = concatPdTables(modedBlastResTab,entryR)
										modedBlastResTab = concatPdTables(modedBlastResTab,entry)
										notAsigned=False
										break
								if notAsigned:
									for idx,entry in left.iterrows():
										leftend = float(entry['send'])
										leftstart = float(entry['sstart'])
										if ((leftend > rightstart - ISMatchLen - 10000 and leftend < rightstart) or (leftstart < rightend + ISMatchLen + 10000 and leftstart > rightend)) and entry['sseqid'] == right.iloc[0,1]:
											modedBlastResTab = concatPdTables(modedBlastResTab,entryR)
											modedBlastResTab = concatPdTables(modedBlastResTab,entry)
											notAsigned=False
											break
								if notAsigned:
									for idx,entry in left.iterrows():
										leftend = float(entry['send'])
										leftstart = float(entry['sstart'])
										contigLen = entry['slen']
										if (entry['length'] > Llen*.9) and ( leftstart < 1000  or leftend > contigLen - 1000 ):
											modedBlastResTab = concatPdTables(modedBlastResTab,entryR)
											modedBlastResTab = concatPdTables(modedBlastResTab,entry)
											notAsigned=False
											break
				else:
					newScaffoldEntries = newScaffoldEntries.sort_values(["evalue"])
					if (newScaffoldEntries.iloc[0,4] > newScaffoldEntries.iloc[0,2]*.9) and (newScaffoldEntries.iloc[0,9] < newScaffoldEntries.iloc[1,9]):
						entry = newScaffoldEntries.iloc[0,:]
						modedBlastResTab = concatPdTables(modedBlastResTab,entry)
			else:
				ISMatchLen = int(re.sub(r'^.*\|ISMatchLen:([0-9]*)\*?\*?$','\\1',scaffoldEntries.iloc[0,0]))
				Llen = scaffoldEntries["Llen"].iloc[0] #.iloc[0,17]
				Rlen = scaffoldEntries["Rlen"].iloc[0]
				#Discriminate cases with 2 QIFS and only one hit...
				if Llen > minLength and Rlen > minLength:
					pass
				else:
					modedBlastResTab = concatPdTables(modedBlastResTab,scaffoldEntries)
			uniqueID+=1
		if len(modedBlastResTab) == 0:
			modedBlastResTab['uniqueID'] = np.nan
			modedBlastResTab['Reversed'] = np.nan
			modedBlastResTab['Extract'] = ""
		return 	modedBlastResTab
	else:
		modedBlastResTab['uniqueID'] = np.nan
		modedBlastResTab['Reversed'] = np.nan
		modedBlastResTab['Extract'] = ""
		return 	modedBlastResTab

def getScaffoldsWithGoodQIFs(table):
	modedBlastResTab = pd.DataFrame(columns=list(table.columns))
	scaffoldNames = list(table["qseqid"].unique())
	if len(scaffoldNames) > 0:
		for Scaffold in scaffoldNames:
			howManyScaffolds = len(table[table["qseqid"] == Scaffold])
			scaffoldEntries = pd.DataFrame(table[table["qseqid"] == Scaffold])
			if howManyScaffolds >= 2:
				newScaffoldEntries = pd.DataFrame(columns=scaffoldEntries.columns)
				newScaffoldEntries[['qlen', 'slen', 'length', 'qstart', 'qend','sstart', 'send', 'bitscore', 'start','end', 'Lend', 'Rstart', 'Llen', 'Rlen']] = newScaffoldEntries[['qlen', 'slen', 'length', 'qstart', 'qend','sstart', 'send', 'bitscore', 'start','end', 'Lend', 'Rstart', 'Llen', 'Rlen']].astype("Int32")
				for idx,scaffold in scaffoldEntries.iterrows():
					newScaffoldEntries['onlyRorL'] = newScaffoldEntries['onlyRorL'].astype('bool')
					newScaffoldEntries.loc[idx] = scaffold
				newScaffoldEntries = newScaffoldEntries.sort_values(["evalue"])
				Llen = newScaffoldEntries['Llen'].iloc[0] #.iloc[0,17]
				Rlen = newScaffoldEntries['Rlen'].iloc[0] #.iloc[0,18]
				left = 	newScaffoldEntries[newScaffoldEntries['Flank'] == 'L']
				right = newScaffoldEntries[newScaffoldEntries['Flank'] == 'R']
				if (len(right) > 1) and (right.iloc[0,9] < right.iloc[1,9]):
					entryR = right.iloc[0,:]
					modedBlastResTab = concatPdTables(modedBlastResTab,entryR)
				else:
					modedBlastResTab = concatPdTables(modedBlastResTab,right)
				if (len(left) > 1) and (left.iloc[0,9] < left.iloc[1,9]):
					entryL = left.iloc[0,:]
					modedBlastResTab = concatPdTables(modedBlastResTab,entryL)
				else:
					modedBlastResTab = concatPdTables(modedBlastResTab,left)
			else:
				modedBlastResTab = concatPdTables(modedBlastResTab,scaffoldEntries)
	return modedBlastResTab

def extractSeqFromGB(gbfileID,IStable,IStableName,surroundingLen, shift):
	#Extracts QIFS from GB file.
	f = open(os.path.join(outputDir,IStableName+".fasta"), "w+")
	fl = open(os.path.join(outputDir,IStableName+"_L.fasta"), "w+")
	fr = open(os.path.join(outputDir,IStableName+"_R.fasta"), "w+")
	records = list(SeqIO.parse(gbfileID, "genbank"))
	for i in range(len(IStable)):
		tag=""
		scaffold = IStable.iloc[i,0]
		coordStart = min(IStable.iloc[i,5],IStable.iloc[i,6])
		coordEnd = max(IStable.iloc[i,5],IStable.iloc[i,6])
		IScoordStart = coordStart
		IScoordEnd = coordEnd
		coordStart = coordStart - shift
		coordEnd = coordEnd + shift
		alnLenght = IStable.iloc[i,4]
		scaffoldLen = IStable.iloc[i,2]
		#RC6
		ISID = re.sub(' ','_',str(IStable.iloc[i,1]))
		if IStable.loc[i,'Observations'] == 'Partial IS Hit':
			tag='*'
		elif IStable.loc[i,'Observations'] == 'Partial or unspecific IS hit':
			tag='**'
		#
		if surroundingLen > (scaffoldLen/3):
			mysurroundingLen = int(scaffoldLen/3)
		else:
			mysurroundingLen = surroundingLen
		#left
		if coordStart > 0:
			if (coordStart - mysurroundingLen - 1) < 0:
				start = 0
				lstart = 1
			elif (coordStart - mysurroundingLen - 1) >= 0:
				start = (coordStart - mysurroundingLen - 1) - 1
				lstart = coordStart - 1 - mysurroundingLen
			end =  coordStart - 1
			if coordStart == 1:
				lend = lstart
			else:
				lend =  coordStart - 1
		else:
			start = 0
			lstart = 1
			end = start
			lend = lstart
		#right
		if coordEnd < scaffoldLen:
			start2 = coordEnd
			if coordEnd == scaffoldLen:
				lstart2 = scaffoldLen
			else:
				lstart2 = coordEnd + 1
			if (coordEnd + 1 + mysurroundingLen) > scaffoldLen:
				end2 = scaffoldLen
			elif (coordEnd + 1 + mysurroundingLen) <= scaffoldLen:
				end2 = coordEnd + 1 + mysurroundingLen
			lend2 = end2
		else:
			start2 = scaffoldLen
			lstart2 = scaffoldLen
			end2 = scaffoldLen
			lend2 = end2
		#write file
		for record in records:
			if record.id == scaffold:
				seq = record.seq[start:end]
				seq2 = record.seq[start2:end2]
				if (len(seq) > minLength ) and (len(seq2)  > minLength ):
					f.write(">" + scaffold + "|LEFT_start:" + str(lstart) + "|LEFT_End:" + str(lend) + "|RIGHT_start:" + str(lstart2) + "|RIGHT_End:" + str(lend2) +  "|IS:" + str(ISID) + "|ISstart:" + str(IScoordStart) +  "|ISend:" + str(IScoordEnd) + "|ISMatchLen:" + str(alnLenght) + str(tag) + "\n" + str(seq) +str(seq2) + "\n")
					fl.write(">" + scaffold + "|LEFT_start:" + str(lstart) + "|LEFT_End:" + str(lend) + "|RIGHT_start:" + str(lstart2) + "|RIGHT_End:" + str(lend2) +  "|IS:" + str(ISID) + "|ISstart:" + str(IScoordStart) +  "|ISend:" + str(IScoordEnd) + "|ISMatchLen:" + str(alnLenght) + str(tag) + "\n" + str(seq) + "\n")
					fr.write(">" + scaffold + "|LEFT_start:" + str(lstart) + "|LEFT_End:" + str(lend) + "|RIGHT_start:" + str(lstart2) + "|RIGHT_End:" + str(lend2) +  "|IS:" + str(ISID) + "|ISstart:" + str(IScoordStart) +  "|ISend:" + str(IScoordEnd) + "|ISMatchLen:" + str(alnLenght) + str(tag) + "\n" + str(seq2) + "\n")
				elif (len(seq) > minLength ) and (len(seq2)  < minLength ):
					fl.write(">" + scaffold + "|LEFT_start:" + str(lstart) + "|LEFT_End:" + str(lend) + "|RIGHT_start:" + str(lstart2) + "|RIGHT_End:" + str(lend2) +  "|IS:" + str(ISID) + "|ISstart:" + str(IScoordStart) +  "|ISend:" + str(IScoordEnd) + "|ISMatchLen:" + str(alnLenght) + str(tag) + "\n" + str(seq) + "\n")
				elif (len(seq) < minLength ) and (len(seq2)  > minLength ):
					fr.write(">" + scaffold + "|LEFT_start:" + str(lstart) + "|LEFT_End:" + str(lend) + "|RIGHT_start:" + str(lstart2) + "|RIGHT_End:" + str(lend2) +  "|IS:" + str(ISID) + "|ISstart:" + str(IScoordStart) +  "|ISend:" + str(IScoordEnd) + "|ISMatchLen:" + str(alnLenght) + str(tag) + "\n" + str(seq2) + "\n")
	f.close()
	fl.close()
	fr.close()
	pass

def extractSeqFromREFGB(gbfileID,IStable,IStableName,surroundingLen2, shift):
	#Extracts RAFS from genbank file.
	f = open(os.path.join(outputDir,IStableName+".fasta"), "w+")
	records = list(SeqIO.parse(gbfileID, "genbank"))
	for i in range(len(IStable)):
		surroundingLen2 = int(surroundingLen2)
		scaffold = IStable.iloc[i,1]
		scaffold2 = re.sub(r'(.*)\|LEFT_start.*$','\\1',IStable.iloc[i,0])
		ISIS = re.sub(r'^.*\|IS:([^|]*).*?$','\\1',IStable.iloc[i,0])
		uniqueID = IStable['uniqueID'].iloc[i]
		reversed = IStable['Reversed'].iloc[i]
		extract = IStable['Extract'].iloc[i]
		# ismatchLen = int(re.sub(r'^.*\|ISMatchLen:([0-9]{1,10}).*$','\\1',IStable.iloc[i,0]))
		# if (ismatchLen > surroundingLen2) & (ismatchLen < surroundingLen):
		# 	surroundingLenby2 = ismatchLen
		coordStart = min(IStable.iloc[i,7],IStable.iloc[i,8])
		coordEnd = max(IStable.iloc[i,7],IStable.iloc[i,8])
		coordStart = coordStart - shift
		coordEnd = coordEnd + shift
		alnLenght = IStable.iloc[i,4]
		scaffoldLen = IStable.iloc[i,3]
		if coordStart > 0:
			if (coordStart - surroundingLen2) < 0:
				start = 1
			elif (coordStart - surroundingLen2) >= 0:
				start = coordStart - surroundingLen2
		else:
			start = 1
			end = start
		if coordEnd < scaffoldLen:
			if (coordEnd + surroundingLen2) > scaffoldLen:
				end2 = scaffoldLen
			elif (coordEnd + surroundingLen2) <= scaffoldLen:
				end2 = coordEnd + surroundingLen2
			end = coordStart
			start2 = coordEnd
		else:
			start2 = scaffoldLen
			end2 = scaffoldLen
		for record in records:
			if record.id == scaffold:
				seq=record.seq[start:end]
				seq2=record.seq[start2:end2]
				if extract == 'L':
					f.write(">" + scaffold2 + "|" + scaffold + "|IS:" +str(ISIS) +"|ID:"+ str(uniqueID) + "|LEFT_start:" + str(start) + "|LEFT_End:" + str(end) + "|RIGHT_start:" + str(start2) + "|RIGHT_End:" + str(end2) + "|refstart:" + str(IStable.iloc[i,7]) + "|refEnd:" + str(IStable.iloc[i,8]) + "\n" + str(seq) + "\n")
				elif  extract == 'R':
					f.write(">" + scaffold2 + "|" + scaffold + "|IS:" +str(ISIS) +"|ID:"+ str(uniqueID) + "|LEFT_start:" + str(start) + "|LEFT_End:" + str(end) + "|RIGHT_start:" + str(start2) + "|RIGHT_End:" + str(end2) + "|refstart:" + str(IStable.iloc[i,7]) + "|refEnd:" + str(IStable.iloc[i,8]) + "\n" + str(seq2) + "\n")
	f.close()
	pass

def concatPdTables(tab1,tab2):
	#Concat a row or a table to a dataframe
	tab1t = tab1.copy()
	if 'ConsecutiveIS' in tab1t.columns:
		tab1t['ConsecutiveIS'] = tab1t['ConsecutiveIS'].astype('bool')
	if 'onlyRorL' in tab1t.columns:
		tab1t['onlyRorL'] = tab1t['onlyRorL'].astype('bool')
	if 'Reversed' in tab1t.columns:
		tab1t['Reversed'] = tab1t['Reversed'].astype('bool')
	if not isinstance(tab2, pd.DataFrame):
		tab2t = pd.DataFrame(tab2).transpose().copy()
	else:
		tab2t = tab2.copy()
	if 'ConsecutiveIS' in tab2t.columns:
		tab2t['ConsecutiveIS'] = tab2t['ConsecutiveIS'].astype('bool')
	if 'onlyRorL' in tab2t.columns:
		tab2t['onlyRorL'] = tab2t['onlyRorL'].astype('bool')
	if 'Reversed' in tab2t.columns:
		tab2t['Reversed'] = tab2t['Reversed'].astype('bool')
	return(pd.concat([tab1t,tab2t]))

def filterLowLen(blastResTab,ISdiff):
	#Filter low length IS blast matches from the analysis. Filters IS scaffolds. Adds information.
	for idx, row in blastResTab.iterrows():
		blastResTab.loc[idx,"IScover"] = (int(blastResTab.loc[idx]["length"]) / int(blastResTab.loc[idx]["slen"]) * 100)
		if blastResTab.loc[idx,"IScover"] >= 100:
			blastResTab.loc[idx,"Observations"] = "Complete IS Hit"
		elif blastResTab.loc[idx,"IScover"] < 100 and blastResTab.loc[idx,"IScover"] > 30:
			blastResTab.loc[idx,"Observations"] = "Partial IS Hit"
		else:
			blastResTab.loc[idx,"Observations"] = "Partial or unspecific IS hit"
	#
	blastResTabfilter = pd.DataFrame(columns=list(blastResTab.columns))
	if 'ConsecutiveIS' in blastResTabfilter.columns:
		blastResTabfilter['ConsecutiveIS'] = blastResTabfilter['ConsecutiveIS'].astype('bool')
	for index, row in blastResTab.iterrows():
		if (row["qlen"]-row["length"] < ISdiff):
			pass
		else:
			if int(row["qstart"]) < 50:
				blastResTabfilter = concatPdTables(blastResTabfilter,row)
			elif int(row["qend"]) > (row["qlen"] - 50):
				blastResTabfilter = concatPdTables(blastResTabfilter,row)
			elif int(row["length"]) > row["slen"]*.3:
				blastResTabfilter = concatPdTables(blastResTabfilter,row)
	return(blastResTabfilter)

def identifyIdenticalScaffolds(rseq, qseq, maxLenScaffolds,	scaffoldDiff):
	#Removes identical scaffolds from the analysis.
	#get number of scaffolds and largest scaffold len
	shiftBasepairs = int(maxLenScaffolds/2)
	refFileFasta = os.path.join(outputDir,rseq+".fasta")
	queryFileFasta = os.path.join(outputDir,qseq+".fasta")
	REFDB = rseq + "_DB"
	resfile =  qseq+"Genome2"+rseq+"Genome.res"
	#add correction for circular genomes..
	records = list(SeqIO.parse(refFileFasta, "fasta"))
	if len(records) < 15:
		f = open( os.path.join(outputDir,"shifted"+rseq+"End.fasta"), 'w+')
		for record in records:
			f.write(">"+record.id+"\n"+str(record.seq)+"\n")
			end = len(record)
			shiftedend = len(record) - shiftBasepairs - 10000
			seq=str(record.seq[shiftedend:end])+str(record.seq[0:shiftBasepairs+10000])
			f.write(">"+record.id+"_shiftedEnd"+"\n"+seq+"\n")
		f.close()
		refFileFastaAndShifted = os.path.join(outputDir,"shifted"+rseq+"End.fasta")
		blastQuery2Ref2(queryFileFasta,refFileFastaAndShifted,REFDB,resfile)
	else:
		blastQuery2Ref2(queryFileFasta,refFileFasta,REFDB,resfile)
	#
	resFilePath =   os.path.join(outputDir,qseq+"Genome2"+rseq+"Genome.res")
	headerBlast = ['qseqid','sseqid','qlen','slen','length','qstart','qend','sstart','send','evalue','bitscore','qcovs','qcovhsp']
	blastResTab = pd.read_csv(resFilePath,sep='\t',names=headerBlast)
	blastResTab = blastResTab[(blastResTab['qcovhsp'] >99 ) & ((blastResTab['qlen'] - blastResTab['length'] ) < scaffoldDiff) ]
	blastResTab = blastResTab.groupby("qseqid", sort=False, as_index=False).agg({'sseqid':['first','count'], 'qlen':'max'})
	identicalToRef = blastResTab[blastResTab['sseqid']['count'] < 10]
	return(identicalToRef)

def tagConsecutiveIS(modedBlastResTab):
	#Adds a flag for consecutive ISs.
	if len(modedBlastResTab) > 0:
		tagconsecutive = modedBlastResTab.sort_values(['qseqid','qstart'])
		tagconsecutive = tagconsecutive.reset_index(drop=True)
		tagconsecutiveLen = len(tagconsecutive)
		tagconsecutive.loc[:,'ConsecutiveIS'] = False
		tagconsecutive['ConsecutiveIS'] = tagconsecutive['ConsecutiveIS'].astype('bool')
		for i in range(0,tagconsecutiveLen):
			if i < tagconsecutiveLen - 1 and tagconsecutive.loc[i,'qseqid'] == tagconsecutive.loc[i+1,'qseqid']:
				difWithNext = tagconsecutive.loc[i+1,'qstart'] - tagconsecutive.loc[i,'qend']
				if difWithNext < surroundingLen*2/5:
					tagconsecutive.loc[i,'ConsecutiveIS'] = True
					tagconsecutive.loc[i+1,'ConsecutiveIS'] = True
	else:
		tagconsecutive = modedBlastResTab
		tagconsecutive['ConsecutiveIS'] = np.nan
	return(tagconsecutive)

def modifyBlastResTab(blastResTab):
	modedBlastResTab = pd.DataFrame(columns=list(blastResTab.columns))
	modedBlastResTab['ConsecutiveIS'] = modedBlastResTab['ConsecutiveIS'].astype('bool')
	scaffoldNames = list(blastResTab["qseqid"].unique())
	for scaffold in scaffoldNames:
		howManyScaffolds = len(blastResTab[blastResTab["qseqid"] == scaffold])
		scaffoldEntries = blastResTab[blastResTab["qseqid"] == scaffold]
		if howManyScaffolds > 1:
			newEntry = pd.DataFrame(columns=list(blastResTab.columns))
			newEntry['ConsecutiveIS'] = newEntry['ConsecutiveIS'].astype('bool')
			scaffoldEntries = scaffoldEntries.sort_values("qstart")
			newEntry.loc[0] = scaffoldEntries.iloc[0,:]
			newIndex = 0
			for ISindex in range(1,howManyScaffolds):
				scaffoldLen = scaffoldEntries["qlen"].iloc[ISindex-1]
				ISstart = scaffoldEntries["qstart"].iloc[ISindex-1]
				ISend = scaffoldEntries["qend"].iloc[ISindex-1]
				nextISstart = scaffoldEntries["qstart"].iloc[ISindex]
				nextISend = scaffoldEntries["qend"].iloc[ISindex]
				#merge [consecutive and]removed overlapped hits
				if (abs(ISstart - nextISstart) < 50) & (abs(ISend - nextISend) < 50):
					if scaffoldEntries.iloc[ISindex,10] > scaffoldEntries.iloc[ISindex-1,10]:
						newEntry['ConsecutiveIS'] = newEntry['ConsecutiveIS'].astype('bool')
						newEntry.loc[newIndex] = scaffoldEntries.iloc[ISindex,:]
					pass
				else:
					newEntry['ConsecutiveIS'] = newEntry['ConsecutiveIS'].astype('bool')
					newIndex = newIndex + 1
					newEntry.loc[newIndex] = scaffoldEntries.iloc[ISindex,:]
			newEntry['ConsecutiveIS'] = newEntry['ConsecutiveIS'].astype('bool')
			modedBlastResTab = 	concatPdTables(modedBlastResTab,newEntry)
		else:
			modedBlastResTab = 	concatPdTables(modedBlastResTab,scaffoldEntries)
	return(modedBlastResTab)

def findISsOnQuery(rseq,qseq,ISdiff,surroundingLen,queryFileFasta):
	#Searches for ISs using blastn.
	refFileGB = os.path.join(outputDir,rseq+".gb")
	refFileFasta = os.path.join(outputDir,rseq+".fasta")
	queryFile = os.path.join(outputDir,qseq+".gb")
	ISFile = isfile
	ISDB = "IS_DB"
	REFDB = rseq + "_DB"
	#Step1: blasting Query genomes against IS database...
	resFile = qseq + "IS.txt"
	blastQuery2IS(queryFileFasta,ISFile,ISDB,resFile)
	#table processing
	headerBlast = ['qseqid','sseqid','qlen','slen','length','qstart','qend','sstart','send','evalue','bitscore','qcovs']
	resFilePath = os.path.join(outputDir,resFile)
	blastResTab = pd.read_csv(resFilePath,sep='\t',names=headerBlast)
	blastResTab = blastResTab.sort_values(by=["qseqid","qstart"])
	modedBlastResTab = blastResTab
	totalISFound=len(modedBlastResTab)
	#filter short IS that are not located near and end of the scaffold
	#Also discard scaffolds that have a small length and match all its length with an IS hit.
	modedBlastResTab = tagConsecutiveIS(modedBlastResTab)
	modedBlastResTab = filterLowLen(modedBlastResTab,ISdiff)
	modedBlastResTab = modifyBlastResTab(modedBlastResTab)
	processedIS=len(modedBlastResTab)
	#stats
	aggModedBlastResTab = modedBlastResTab.groupby('qseqid',as_index=False).agg({'sseqid':['first','count']})
	scaffoldsWithISonBothEnds = list( aggModedBlastResTab[aggModedBlastResTab['sseqid']['count'] >1].iloc[:,0])
	NscaffoldsWithISonBothEnds = len(scaffoldsWithISonBothEnds)
	scaffoldsWithISonOneEnd = list( aggModedBlastResTab[aggModedBlastResTab['sseqid']['count'] == 1].iloc[:,0])
	NscaffoldsWithISonOneEnd = len(scaffoldsWithISonOneEnd)
	stats = list([processedIS,NscaffoldsWithISonBothEnds,scaffoldsWithISonBothEnds,NscaffoldsWithISonOneEnd,scaffoldsWithISonOneEnd,totalISFound])
	modedBlastResTab = modedBlastResTab.reset_index(drop=True)
	return(modedBlastResTab,stats)

def testMiddleIS(rseq, qseq, modedBlastResTab):
	#Searches for IS location chages using full length QIFS matches.
	#rseq and qseq contain the basename of the concatenated gb file!
	refFileGB = os.path.join(outputDir,rseq+".gb")
	refFileFasta = os.path.join(outputDir,rseq+".fasta")
	queryFile = os.path.join(outputDir,qseq+".gb")
	ISFile = isfile
	ISDB = "IS_DB"
	REFDB = rseq + "_DB"
	headerBlast = ['qseqid','sseqid','qlen','slen','length','qstart','qend','sstart','send','evalue','bitscore','qcovs']
	querySurroundings = qseq + "_surroundings"
	extractSeqFromGB(queryFile,modedBlastResTab,querySurroundings,surroundingLen,shift)
	qsurroundings = os.path.join(outputDir, querySurroundings + ".fasta")
	testISfile = qseq + "_" + qseq + "_testIS.txt"
	blastQuery2IS2(qsurroundings,ISFile,ISDB,testISfile)
	resFilePath = os.path.join(outputDir,testISfile)
	consecutiveIS = pd.read_csv(resFilePath,sep='\t',names=headerBlast)
	if len(consecutiveIS) > 0:
		consecutiveIS["ISID"]  = consecutiveIS.apply(lambda x: re.sub(r'^.*\|IS:([^|]*)\|.*$','\\1',x["qseqid"]),axis=1)
		consecutiveIS = consecutiveIS[consecutiveIS['qcovs']>10]
		consecutiveIS['SameIS'] = False
		#
		for idx,row in consecutiveIS.iterrows():
			if len(set(re.sub(r'\/','_',row['sseqid']).split('_')) & set(re.sub(r'\/','_',row['ISID']).split('_'))) > 0:
				consecutiveIS.loc[idx,'SameIS'] = True
		#
		consecutiveIS = consecutiveIS[consecutiveIS['SameIS']]
		consecutiveISUQ = list(consecutiveIS['qseqid'].unique())
		consecutiveIS = consecutiveIS[consecutiveIS['qseqid'].isin(consecutiveISUQ)]
		consecutiveIS = consecutiveIS.reset_index(drop=True)
	else:
		consecutiveIS["ISID"]  = np.nan
		consecutiveIS['SameIS'] = False
		consecutiveISUQ = list(consecutiveIS['qseqid'].unique())
	#******
	outfile = qsurroundings
	resFile = qseq + "_surroundingBlastRes.txt"
	blastQuery2Ref(outfile,refFileFasta,REFDB,resFile)
	resFilePath = os.path.join(outputDir,resFile)
	blastResTab = pd.read_csv(resFilePath,sep='\t',names=headerBlast)
	if len(blastResTab) > 0:
		blastResTab["Scaffold"] = blastResTab.apply(lambda x: re.sub(r'^([^|]*)\|.*','\\1',x["qseqid"]), axis=1)
		blastResTab["start"]  = blastResTab.apply(lambda x: int(re.sub(r'^.*\|LEFT_start:([0-9]{1,10})\|.*$','\\1',x["qseqid"])),axis=1)
		blastResTab["end"]  = blastResTab.apply(lambda x: int(re.sub(r'^.*\|RIGHT_End:([0-9]{1,10})\|.*$','\\1',x["qseqid"])),axis=1)
		blastResTab["leftLen"]  = blastResTab.apply(lambda x: int(re.sub(r'^.*\|LEFT_End:([0-9]{1,10})\|.*$','\\1',x["qseqid"])), axis=1) - blastResTab.apply(lambda x: int(re.sub(r'^.*\|LEFT_start:([0-9]{1,10})\|.*$','\\1',x["qseqid"])), axis=1)
		blastResTab["rightLen"]  = blastResTab.apply(lambda x: int(re.sub(r'^.*\|RIGHT_End:([0-9]{1,10})\|.*$','\\1',x["qseqid"])), axis=1) - blastResTab.apply(lambda x: int(re.sub(r'^.*\|RIGHT_start:([0-9]{1,10})\|.*$','\\1',x["qseqid"])), axis=1)
		blastResTab.iloc[:,13] = pd.to_numeric(blastResTab.iloc[:,13])
		blastResTab.iloc[:,14] = pd.to_numeric(blastResTab.iloc[:,14])
		blastResTab["tag"] = blastResTab.apply(lambda x: re.sub(r'^[^*]*(\*?\*?)$','\\1',x["qseqid"]), axis=1)
	else:
		blastResTab["Scaffold"] = np.nan
		blastResTab["start"]  = np.nan
		blastResTab["end"] = np.nan
		blastResTab["end"] = np.nan
		blastResTab["leftLen"] = np.nan
		blastResTab["rightLen"] = np.nan
		blastResTab["tag"] = np.nan
	#Subjects with > 70 qc are either missing iSs or internal ISs
	missingQueryIS = blastResTab
	#surroundingLen/5 because surroundingLen/(surroundingLen/5+surroundingLen) is less than 90%, and in the case where there is alignment only in one side that will be the max cover
	missingQueryIS = missingQueryIS[(missingQueryIS["rightLen"] > surroundingLen/5 ) & (missingQueryIS["leftLen"] > surroundingLen/5 ) & (missingQueryIS["length"]/missingQueryIS["qlen"] > 0.9)]
	missingQueryIS = missingQueryIS.iloc[:,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,17]]	#new
	if len(missingQueryIS) > 0:
		missingQueryIS.loc[missingQueryIS[~missingQueryIS['qseqid'].isin(consecutiveISUQ)].index,'Observations'] = '1. DLIS'
		missingQueryIS.loc[missingQueryIS[missingQueryIS['qseqid'].isin(consecutiveISUQ)].index,'Observations'] = '2. Verify manually, possible false positive caused by consecutive ISs.'
		missingQueryIS['IS_Match_Type'] = 'Complete IS hit.'
		missingQueryIS.loc[missingQueryIS[missingQueryIS['tag'] == '*'].index,'IS_Match_Type'] = 'Partial IS hit.'
		missingQueryIS.loc[missingQueryIS[missingQueryIS['tag'] == '**'].index,'IS_Match_Type'] = 'Partial or unspecific IS hit.'
	else:
		missingQueryIS['Observations'] = np.nan
		missingQueryIS['IS_Match_Type'] = np.nan
	missingQueryIS = missingQueryIS.iloc[:,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17]]
	missingQueryIS = missingQueryIS.reset_index(drop=True)
	for idx,row in missingQueryIS.iterrows():
		missingQueryIS.loc[idx,'sstart'] = min(row['sstart'],row['send'])
		missingQueryIS.loc[idx,'send'] = max(row['sstart'],row['send'])
	return(missingQueryIS,consecutiveIS)

def testIS(rseq, qseq,MiddleMissingISID,consecutiveIS):
	#Searches for IS location changes using half QIFS matches.
	refFileGB = os.path.join(outputDir,rseq+".gb")
	refFileFasta = os.path.join(outputDir,rseq+".fasta")
	querySurroundings_L = qseq + "_surroundings_L"
	querySurroundings_R = qseq + "_surroundings_R"
	ISFile = isfile
	ISDB = "IS_DB"
	REFDB = rseq + "_DB"
	headerBlast = ['qseqid','sseqid','qlen','slen','length','qstart','qend','sstart','send','evalue','bitscore','qcovs']
	qsurroundings_L = os.path.join(outputDir , querySurroundings_L + ".fasta")
	qsurroundings_R = os.path.join(outputDir , querySurroundings_R + ".fasta")
	resFile_L = qseq + "_surroundingBlastRes_L.txt"
	resFile_R = qseq + "_surroundingBlastRes_R.txt"
	# Blast 2 flanks independently and merge tables
	blastQuery2Ref(qsurroundings_L,refFileFasta,REFDB,resFile_L)
	blastQuery2Ref(qsurroundings_R,refFileFasta,REFDB,resFile_R)
	resFilePath = os.path.join(outputDir , resFile_L)
	blastResTab = pd.read_csv(resFilePath,sep='\t',names=headerBlast)
	if len(blastResTab) > 0:
		blastResTab["Flank"] = "L"
	else:
		blastResTab["Flank"] = np.nan
	resFilePath = os.path.join(outputDir , resFile_R)
	blastResTab_R = pd.read_csv(resFilePath,sep='\t',names=headerBlast)
	if len(blastResTab_R) > 0:
		blastResTab_R["Flank"] = "R"
	else:
		blastResTab_R["Flank"] = np.nan
	blastResTab = concatPdTables(blastResTab,blastResTab_R)
	blastResTab = blastResTab.reset_index(drop=True)
	#Process tables
	blastResTab = blastResTab[~blastResTab['qseqid'].isin(MiddleMissingISID)]
	##
	if len(blastResTab) > 0:
		blastResTab["Scaffold"] = blastResTab.apply(lambda x: re.sub(r'^([^|]*)\|.*','\\1',x["qseqid"]), axis=1)
		blastResTab["start"] = blastResTab.apply(lambda x: re.sub(r'^.*\|LEFT_start:([0-9]{1,10})\|.*$','\\1',x["qseqid"]), axis=1)
		blastResTab["end"] = blastResTab.apply(lambda x: re.sub(r'^.*\|RIGHT_End:([0-9]{1,10})\|.*$','\\1',x["qseqid"]), axis=1)
		#eliminate surrounds that are not present on the refgenome...
		blastResTab["Lend"] = blastResTab.apply(lambda x: re.sub(r'^.*\|LEFT_End:([0-9]{1,10})\|.*$','\\1',x["qseqid"]), axis=1)
		blastResTab["Rstart"] = blastResTab.apply(lambda x: re.sub(r'^.*\|RIGHT_start:([0-9]{1,10})\|.*$','\\1',x["qseqid"]), axis=1)
		# establish cutoff values
		for idx, row in blastResTab.iterrows():
			blastResTab.loc[idx,"Llen"] = int(blastResTab.loc[idx]["Lend"]) - int(blastResTab.loc[idx]["start"])
			blastResTab.loc[idx,"Rlen"] = int(blastResTab.loc[idx]["end"]) - int(blastResTab.loc[idx]["Rstart"])
			if (blastResTab.loc[idx,"Llen"] > minLength) & (blastResTab.loc[idx,"Rlen"] > minLength):
				blastResTab.loc[idx,"onlyRorL"] = False
			else:
				blastResTab.loc[idx,"onlyRorL"] = True
			blastResTab.loc[idx,"cutoff"] = 0.90 #qc cutoff
			blastResTab.loc[idx,"qcovh"] = blastResTab.loc[idx]["length"]/blastResTab.loc[idx]["qlen"]
			blastResTab["ISID"] = blastResTab.apply(lambda x: re.sub(r'^.*\|IS:([^|]*).*$','\\1',x["qseqid"]), axis=1)
	else:
		blastResTab["Scaffold"] = np.nan
		blastResTab["start"] = np.nan
		blastResTab["end"] = np.nan
		blastResTab["Lend"] = np.nan
		blastResTab["Rstart"] = np.nan
		blastResTab["Llen"] = np.nan
		blastResTab["Rlen"] = np.nan
		blastResTab["onlyRorL"] = np.nan
		blastResTab["cutoff"] = np.nan
		blastResTab["qcovh"] = np.nan
		blastResTab["ISID"] = np.nan
	#filter only left or right or both sides sequences blast hits by cutoffs. Partial or bad hits.
	keep = blastResTab[(blastResTab['length'] > minAlnLength) & (blastResTab['qcovh'] > blastResTab["cutoff"])]
	#Find pairs of hits, each half of the concatenated sequence
	keepUL = list(keep['qseqid'].unique())
	keep = findpairs(keep)
	#
	#Discard mathces on scaffold ends,
	matchOnScaffoldEnd = keep[((keep['sstart'] < 25) & (keep['Extract'] == 'L'))]
	matchOnScaffoldEnd = concatPdTables(matchOnScaffoldEnd,keep[((keep['send'] > keep['slen'] - 25) & (keep['Extract'] == 'R'))])
	matchOnScaffoldEnd = list(matchOnScaffoldEnd['qseqid'].unique())
	keepUL2 = list(keep['qseqid'].unique())
	#
	#extract sequences near the surroundings on ref genome...Needs sstart, send...
	#refstart and refend are the positions of query-IS match with ref!
	refISSurroundsFromQuery = qseq + "_" + rseq + "ISSurroundsFrom" + qseq
	extractSeqFromREFGB(refFileGB,keep,refISSurroundsFromQuery,surroundingLen,shift)
	refISSurroundsFromQuery = os.path.join(outputDir , refISSurroundsFromQuery + ".fasta")
	#Blast extracted sequences to IS!
	testISfile = qseq + "_" + rseq + "_testIS.txt"
	blastQuery2IS2(refISSurroundsFromQuery,ISFile,ISDB,testISfile)
	#Load table and determine sequences surrounded by IS on reference
	resFilePath = os.path.join(outputDir,testISfile)
	blastResTab2 = pd.read_csv(resFilePath,sep='\t',names=headerBlast)
	#Merge tables to find those that didn't have an IS. Left_only
	# make tables compatible
	#
	if len(blastResTab2) > 0:
		blastResTab2["refstart"] = blastResTab2.apply(lambda x: re.sub(r".*\|refstart:([0-9]{1,10}).*$","\\1",x["qseqid"]), axis = 1)
		blastResTab2["refend"] = blastResTab2.apply(lambda x: re.sub(r".*\|refEnd:([0-9]{1,10}).*$","\\1",x["qseqid"]), axis = 1)
		blastResTab2[blastResTab2.columns[12]] = pd.to_numeric(blastResTab2[blastResTab2.columns[12]])
		blastResTab2[blastResTab2.columns[13]] = pd.to_numeric(blastResTab2[blastResTab2.columns[13]])
		blastResTab2["Query.ID"] = blastResTab2.apply(lambda x: re.sub(r"^([^|]*).*$","\\1",x["qseqid"]), axis = 1)
		blastResTab2["ISID"] = blastResTab2.apply(lambda x: re.sub(r".*\|IS:([^|]*).*$","\\1",x["qseqid"]), axis = 1)
		blastResTab2["uniqueID"] = blastResTab2.apply(lambda x: float(re.sub(r".*\|ID:([^|]*).*$","\\1",x["qseqid"])), axis = 1)
	else:
		blastResTab2["refstart"] = np.nan
		blastResTab2["refend"] = np.nan
		blastResTab2["Query.ID"] = np.nan
		blastResTab2["ISID"] = np.nan
		blastResTab2["uniqueID"] = np.nan
	#
	#eq IS, this because different IS may represent positives
	for idx,row in blastResTab2.iterrows():
		if len(set(re.sub(r'\/','_',row['sseqid']).split('_')) & set(re.sub(r'\/','_',row['ISID']).split('_'))) > 0:
			blastResTab2.loc[idx,'eqIS']=True
		else:
			blastResTab2.loc[idx,'eqIS']=False
	#
	if len(blastResTab2) > 0:
		blastResTab2 = blastResTab2[blastResTab2['eqIS']]
		#
		blastResTab2 = blastResTab2.groupby(['qseqid','refstart'], sort=False, as_index=False)
		blastResTab2 = blastResTab2.first()
		blastResTab2 = pd.DataFrame(blastResTab2)
		blastResTab2 = blastResTab2.reset_index()
		blastResTab2 = blastResTab2.groupby(['uniqueID'], sort=False, as_index=False).agg({'ISID':['first','count']})
		blastResTab2 = pd.DataFrame(blastResTab2)
		blastResTab2.columns = [re.sub('_$','','_'.join(col).strip()) for col in blastResTab2.columns.values]
		blastResTab2 = blastResTab2.reset_index(drop=True)
	else:
		blastResTab2 = pd.DataFrame(columns=['uniqueID', 'ISID_first', 'ISID_count'])
	#
	#Group blastrestab
	if len(keep) > 0:
		#Keep repeated regions..
		keepRep = keep.groupby(['sseqid','sstart'], sort=False, as_index=False)
		# keepRep = keepRep.agg({'Scaffold': ['count']})
		keepRep = keepRep.agg({'Scaffold': ['count'], 'qseqid':lambda x: list(x)})
		keepRep = keepRep[keepRep['Scaffold']['count'] > 1]['qseqid']['<lambda>']
		keepRepList = []
		for row in keepRep:
			for item in row:
				keepRepList.append(item)
		keepRepList = set(keepRepList)
		#
		keepG = keep.groupby(['qseqid','uniqueID'], sort=False, as_index=False)
		keepG = keepG.agg({'Scaffold': ['count'],'sseqid':['first','last'],'ISID':['first'],'start':['first'], 'Lend':['first'], 'Rstart':['first'],'end':['first'], 'sstart':['first','last'], 'send':['first','last'], 'Flank':['first','last']})
		keepG.columns = [re.sub('_$','','_'.join(col).strip()) for col in keepG.columns.values]
		keepG = keepG.reset_index(drop=True)
		for idx,row in keepG.iterrows():
			if row['Flank_first'] == 'R':
				startend = [row['Rstart_first'],row['end_first']]
				startend2 = [row['start_first'],row['Lend_first']]
				keepG.loc[idx, ['start_first','Lend_first']] = startend2
				keepG.loc[idx, ['Rstart_first','end_first']] = startend
		keepG=keepG.iloc[:,0:14]
		#
		outermerge = pd.merge(keepG, blastResTab2, how='outer', left_on=['uniqueID'], right_on=['uniqueID'], indicator=True)
		missingQueryIS2 = outermerge[outermerge["_merge"] == "left_only"].iloc[:,0:14]
		columns = ['qseqid', 'uniqueID', 'Scaffold_count', 'sseqid_first','sseqid_last', 'ISID', 'start', 'Lend', 'Rstart', 'end', 'sstart_first', 'sstart_last', 'send_first', 'send_last']
		missingQueryIS2.columns = columns
		missingQueryIS2['Observations'] = '1. DLIS'
		#Test Consecutive IS (same IS)
		missingQueryIS2.loc[missingQueryIS2['qseqid'].isin(list(consecutiveIS['qseqid'])),'Observations'] = '2. Verify manually, possible false positive caused by consecutive ISs.'
		#Test repeated sequences (different ISs and other)
		missingQueryIS2.loc[missingQueryIS2['qseqid'].isin(keepRepList),'Observations'] = '4. Verify manually, possible false positive produced by repeated sequences [Different QIFs match with same region].'
		#
		missingQueryIS2.loc[missingQueryIS2['qseqid'].isin(matchOnScaffoldEnd),'Observations'] = '5. Verify manually, possible false positive produced by QIF match on scaffold end.'
		#not found by blast
		#SLIS Same location IS
		SLIS = outermerge[outermerge["_merge"] == "both"].iloc[:,0:14]
		columns = ['qseqid', 'uniqueID', 'Scaffold_count', 'sseqid_first','sseqid_last', 'ISID', 'start', 'Lend', 'Rstart', 'end', 'sstart_first', 'sstart_last', 'send_first', 'send_last']
		SLIS.columns = columns
		SLIS['Observations'] = 'SLIS'
		#
		# Hits with 2 half sequences and only one IS match on RAFS
		outermerge21 = outermerge[(outermerge['Scaffold_count'] == 2) & (outermerge['ISID_count'] == 1) ]
		columns = ['qseqid', 'uniqueID', 'Scaffold_count', 'sseqid_first','sseqid_last', 'ISID', 'start', 'Lend', 'Rstart', 'end', 'sstart_first', 'sstart_last', 'send_first', 'send_last', 'ISID_', 'ISID_','_merge']
		outermerge21.columns = columns
		outermerge21 = outermerge21.iloc[:,0:14]
		outermerge21['Observations'] = '3. Verify manually, Possible false positive produced by consecutive IS [2 RAFs found, only one with an IS].'
		missingQueryIS2 = concatPdTables(missingQueryIS2,outermerge21)
	else:
		keepRepList = set()
		missingQueryIS2 = pd.DataFrame(columns = ['qseqid', 'uniqueID', 'Scaffold_count', 'sseqid_first', 'sseqid_last', 'ISID', 'start', 'Lend', 'Rstart', 'end', 'sstart_first', 'sstart_last', 'send_first', 'send_last','Observations','IS_Match_Type','tag','Scaffold'])
		SLIS = pd.DataFrame(columns = ['qseqid', 'uniqueID', 'Scaffold_count', 'sseqid_first', 'sseqid_last', 'ISID', 'start', 'Lend', 'Rstart', 'end', 'sstart_first', 'sstart_last', 'send_first', 'send_last','Observations','IS_Match_Type','tag','Scaffold'])
	#
	# Analysed scaffolds.. we have to test if the surrounds of the IS are present on the ref genome. If not, the surrounding should be removed..
	analyzedScaffolds = []
	records = list(SeqIO.parse(qsurroundings_L,'fasta'))
	for record in records:
		analyzedScaffolds.append(record.id)
	records = list(SeqIO.parse(qsurroundings_R,'fasta'))
	for record in records:
		analyzedScaffolds.append(record.id)
	#
	analyzedScaffolds = list(set(analyzedScaffolds))
	analyzedScaffolds = pd.DataFrame({'qseqid':analyzedScaffolds}, dtype=object)
	analyzedScaffolds = analyzedScaffolds[~analyzedScaffolds['qseqid'].isin(MiddleMissingISID)]
	analyzedScaffolds['Scaffold_count'] = np.nan
	analyzedScaffolds["uniqueID"] = np.nan
	analyzedScaffolds['sseqid_first'] = np.nan
	analyzedScaffolds['sseqid_last'] = np.nan
	if len(analyzedScaffolds) > 0:
		analyzedScaffolds["ISID"] = analyzedScaffolds.apply(lambda x: re.sub(r'^.*\|IS:([^|]*).*$','\\1',x["qseqid"]), axis=1)
		analyzedScaffolds["start"] = analyzedScaffolds.apply(lambda x: re.sub(r'^.*\|LEFT_start:([0-9]{1,10})\|.*$','\\1',x["qseqid"]), axis=1)
		analyzedScaffolds["end"] = analyzedScaffolds.apply(lambda x: re.sub(r'^.*\|RIGHT_End:([0-9]{1,10})\|.*$','\\1',x["qseqid"]), axis=1)
		analyzedScaffolds["Lend"] = analyzedScaffolds.apply(lambda x: re.sub(r'^.*\|LEFT_End:([0-9]{1,10})\|.*$','\\1',x["qseqid"]), axis=1)
		analyzedScaffolds["Rstart"] = analyzedScaffolds.apply(lambda x: re.sub(r'^.*\|RIGHT_start:([0-9]{1,10})\|.*$','\\1',x["qseqid"]), axis=1)
	else:
		analyzedScaffolds["ISID"] = np.nan
		analyzedScaffolds["start"] = np.nan
		analyzedScaffolds["end"] = np.nan
		analyzedScaffolds["Lend"] = np.nan
		analyzedScaffolds["Rstart"] = np.nan
	analyzedScaffolds['sstart_first'] = np.nan
	analyzedScaffolds['sstart_last'] = np.nan
	analyzedScaffolds['send_first'] = np.nan
	analyzedScaffolds['send_last'] = np.nan
	#not found by Blastn
	notfound = pd.DataFrame(analyzedScaffolds[~analyzedScaffolds['qseqid'].isin(list(blastResTab['qseqid'].unique()))])
	notfound['Observations'] = '8. Not found by Blastn.'
	found = pd.DataFrame(analyzedScaffolds[analyzedScaffolds['qseqid'].isin(list(blastResTab['qseqid'].unique()))])
	#
	#Discarded for the analysis because of low query coverage
	lowqc = pd.DataFrame(found[~found['qseqid'].isin(keepUL)])
	lowqc['Observations'] = '7. Discarded for analysis, blastn hits with low query coverage.'
	goodqc = pd.DataFrame(found[found['qseqid'].isin(keepUL)])
	#Discarded because of multiple blast hits
	multiplehits = pd.DataFrame(goodqc[~goodqc['qseqid'].isin(keepUL2)])
	multiplehits['Observations'] = '6. Discarded for analysis. There were multiple blastn hits, and left and right pairs could not be determined.'
	#append into a single table
	missedScaffolds = concatPdTables(notfound,lowqc)
	missedScaffolds = concatPdTables(missedScaffolds,multiplehits)
	#
	#incorporate data from blastResTab
	for idx,scaffold in blastResTab.iterrows():
		blastResTab.loc[idx,'sstart'] = min(scaffold['sstart'],scaffold['send'])
		blastResTab.loc[idx,'send'] = max(scaffold['sstart'],scaffold['send'])
	#Filter by scaffold for first pair o hits where there are differences on e-value.
	group = getScaffoldsWithGoodQIFs(blastResTab)
	group =group.groupby(['qseqid'], sort=False, as_index=False)
	group = group.agg({'Scaffold': ['count'],'sseqid':['first','last'],'ISID':['first'],'start':['first'], 'Lend':['first'], 'Rstart':['first'],'end':['first'], 'sstart':['first','last'], 'send':['first','last']})
	group = group[group['Scaffold']['count'] < 3]
	group.columns = [re.sub('_$','','_'.join(col).strip()) for col in group.columns.values]
	group = group.reset_index(drop=True)
	groupList = list(group['qseqid'].unique())
	#
	for idx, row in missedScaffolds.iterrows():
		if len(group[group['qseqid'] == row['qseqid']]) == 1:
			missedScaffolds.loc[idx,'sseqid_first'] = group[group['qseqid'] == row['qseqid']]['sseqid_first'].iloc[0]
			missedScaffolds.loc[idx,'sseqid_last'] = group[group['qseqid'] == row['qseqid']]['sseqid_last'].iloc[0]
			missedScaffolds.loc[idx,'sstart_first'] = group[group['qseqid'] == row['qseqid']]['sstart_first'].iloc[0]
			missedScaffolds.loc[idx,'sstart_last'] = group[group['qseqid'] == row['qseqid']]['sstart_last'].iloc[0]
			missedScaffolds.loc[idx,'send_first'] = group[group['qseqid'] == row['qseqid']]['send_first'].iloc[0]
			missedScaffolds.loc[idx,'send_last'] = group[group['qseqid'] == row['qseqid']]['send_last'].iloc[0]
	missedScaffolds = missedScaffolds[['qseqid', 'uniqueID','Scaffold_count','sseqid_first','sseqid_last','ISID','start', 'end', 'Lend', 'Rstart','sstart_first','sstart_last', 'send_first', 'send_last', 'Observations']]
	missingQueryIS2 = concatPdTables(missingQueryIS2,missedScaffolds)
	#
	if args.reportSLIS and len(SLIS) > 0:
		missingQueryIS2 = concatPdTables(missingQueryIS2,SLIS)
	#format table
	if len(missingQueryIS2) > 0:
		missingQueryIS2["Scaffold"] = missingQueryIS2.apply(lambda x: re.sub(r'^([^|]*)\|.*','\\1',x["qseqid"]), axis=1)
		missingQueryIS2["tag"] = missingQueryIS2.apply(lambda x: re.sub(r'^[^*]*(\*?\*?)$','\\1',x["qseqid"]), axis=1)
		missingQueryIS2.loc[missingQueryIS2[missingQueryIS2['tag'] == '*'].index,'IS_Match_Type'] = 'Partial IS hit.'
		missingQueryIS2.loc[missingQueryIS2[missingQueryIS2['tag'] == '**'].index,'IS_Match_Type'] = 'Partial or unspecific IS hit.'
	else:
		missingQueryIS2["Scaffold"] = np.nan
		missingQueryIS2['IS_Match_Type'] = np.nan
	#
	missingQueryIS2 = missingQueryIS2.reset_index(drop=True)
	missingQueryIS2 = missingQueryIS2.iloc[:,0:20]
	for idx,row in missingQueryIS2.iterrows():
		if row['sstart_first'] == row['sstart_last']:
			missingQueryIS2.loc[idx, ['sseqid_last','sstart_last','send_last']] = np.nan
		if row['sseqid_first'] == row['sseqid_last']:
			missingQueryIS2.loc[idx, ['sseqid_last']] = np.nan
	#
	return missingQueryIS2

def annotateTable(results):
	columns = list(results.columns)
	columns.extend(["QUERY_Flank1.locus_tag","QUERY_Flank1.Product","QUERY_Flank2.locus_tag","QUERY_Flank2.Product","REF_Flank1.locus_tag","REF.Product_Flank1","REF_Flank2.locus_tag","REF.Product_Flank2","Scaffold_Size"])
	FinalResults = pd.DataFrame(columns = columns)
	i=0
	#annotate
	nullAnnotation = pd.DataFrame(columns=['ID', 'locus_tag', 'Product', 'Start', 'End','Size'])
	for index, result in results.iterrows():
		resultList = list(result)
		start1toEnd2=False
		#query
		if result["Description"] == "IS only present in query Strain.":
			qr = queryRecords[result["Query.ID1"]]
			featureList = [feature for feature in qr.features if feature.type == 'CDS']
			querySize = len(qr)
			scaffoldSize = len(qr)
			if len(featureList) >0:
				#left flank
				if not pd.isna(result["Start1"]) and not pd.isna(result["End1"]) and result["End1"] < 50:
					resultList.extend(["-","Contig end."])
				elif not pd.isna(result["Start1"]) and not pd.isna(result["End1"]) and result["End1"] >= 50:
					startend = set(range(int(result["Start1"]),int(result["End1"])))
					product = []
					locus_tag = []
					for  f in featureList:
						fspan = set(range(f.location.start,f.location.end))
						if 'locus_tag' in f.qualifiers.keys():
							ltag = f.qualifiers['locus_tag'][0]
						else:
							ltag = '-'
						if 'product' in f.qualifiers.keys():
							prod = f.qualifiers['product'][0]
						else:
							prod = '-'
						if fspan & startend:
							product.append(prod)
							locus_tag.append(ltag)
					if len(locus_tag) != 0:
						resultList.extend([' | '.join(locus_tag),' | '.join(product)])
					else:
						resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
				#right flank
				if not pd.isna(result["Start2"]) and not pd.isna(result["End2"]) and result["Start2"] > scaffoldSize - 50:
						resultList.extend(["-","Contig end."])
				elif not pd.isna(result["Start2"]) and not pd.isna(result["End2"]) and result["Start2"] <= scaffoldSize - 50:
					startend = set(range(int(result["Start2"]),int(result["End2"])))
					product = []
					locus_tag = []
					for  f in featureList:
						fspan = set(range(f.location.start,f.location.end))
						if 'locus_tag' in f.qualifiers.keys():
							ltag = f.qualifiers['locus_tag'][0]
						else:
							ltag = '-'
						if 'product' in f.qualifiers.keys():
							prod = f.qualifiers['product'][0]
						else:
							prod = '-'
						if fspan & startend:
							product.append(prod)
							locus_tag.append(ltag)
					if len(locus_tag) != 0:
						resultList.extend([' | '.join(locus_tag),' | '.join(product)])
					else:
						resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
			else:
				resultList.extend(["-","-","-","-"])
			#Annotate ref
			#left flank
			if not pd.isna(result["REF.ID1"]):
				rr = refRecords[result["REF.ID1"]]
				featureList = [feature for feature in rr.features if feature.type == 'CDS']
				scaffoldSize = len(rr)
				if len(featureList) > 0:
					if not pd.isna(result["REF.Start1"]) and not pd.isna(result["REF.End1"]) and result["REF.End1"] < 50:
						resultList.extend(["-","Contig end."])
					elif not pd.isna(result["REF.Start1"]) and not pd.isna(result["REF.End1"]) and result["REF.Start1"] > scaffoldSize - 50:
						resultList.extend(["-","Contig end."])
					elif not pd.isna(result["REF.Start1"]) and not pd.isna(result["REF.End1"]):
						startend = set(range(int(result["REF.Start1"]),int(result["REF.End1"])))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
					# start1 to end2
					elif not pd.isna(result["REF.Start1"]) and pd.isna(result["REF.Start2"]) and pd.isna(result["REF.End1"]) and not pd.isna(result["REF.End2"]):
						start1toEnd2=True
						startend = set(range(int(result["REF.Start1"]),int(result["REF.Start1"])+surroundingLen))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
						#
						startend = set(range(int(result["REF.End2"])-surroundingLen,int(result["REF.End2"])))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
				else:
					resultList.extend(["-","-"])
			else:
				resultList.extend(["-","-"])
			#right flank on same contig on only one flank
			if not pd.isna(result["REF.ID1"]) and pd.isna(result["REF.ID2"]) and not start1toEnd2:
				rr = refRecords[result["REF.ID1"]]
				featureList = [feature for feature in rr.features if feature.type == 'CDS']
				scaffoldSize = len(rr)
				if len(featureList) > 0:
					if not pd.isna(result["REF.Start2"]) and not pd.isna(result["REF.End2"]) and result["REF.End2"] < 50:
						resultList.extend(["-","Contig end."])
					if not pd.isna(result["REF.Start2"]) and not pd.isna(result["REF.End2"]) and result["REF.Start2"] > scaffoldSize - 50:
						resultList.extend(["-","Contig end."])
					elif not pd.isna(result["REF.Start2"]) and not pd.isna(result["REF.End2"]) and result["REF.Start2"]:
						startend = set(range(int(result["REF.Start2"]),int(result["REF.End2"])))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
					else:
						resultList.extend(["-","-"])
				else:
					resultList.extend(["-","-"])
			#right flank on different contig
			elif not pd.isna(result["REF.ID2"]) and not start1toEnd2:
				rr = refRecords[result["REF.ID2"]]
				featureList = [feature for feature in rr.features if feature.type == 'CDS']
				scaffoldSize = len(rr)
				if len(featureList) > 0:
					if not pd.isna(result["REF.Start2"]) and not pd.isna(result["REF.End2"]) and result["REF.End2"] < 50:
						resultList.extend(["-","Contig end."])
					if not pd.isna(result["REF.Start2"]) and not pd.isna(result["REF.End2"]) and result["REF.Start2"] > scaffoldSize - 50:
						resultList.extend(["-","Contig end."])
					elif not pd.isna(result["REF.Start2"]) and not pd.isna(result["REF.End2"]) and result["REF.Start2"]:
						startend = set(range(int(result["REF.Start2"]),int(result["REF.End2"])))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
				else:
					resultList.extend(["-","-"])
			elif not start1toEnd2:
				resultList.extend(["-","-"])
			resultList.extend([querySize])
		#******************************************************************
		#** ref ***********************************************************
		elif result["Description"] == "IS only present in ref Strain.":
			# annotate query
			if not pd.isna(result["Query.ID1"]):
				qr = queryRecords[result["Query.ID1"]]
				featureList = [feature for feature in qr.features if feature.type == 'CDS']
				scaffoldSize = len(qr)
				#left flank
				if len(featureList) > 0:
					if not pd.isna(result["Start1"]) and not pd.isna(result["End1"]) and result["End1"] < 50:
						resultList.extend(["-","Contig end."])
					if not pd.isna(result["Start1"]) and not pd.isna(result["End1"]) and result["Start1"] > scaffoldSize - 50:
						resultList.extend(["-","Contig end."])
					elif not pd.isna(result["Start1"]) and not pd.isna(result["End1"]) and result["End1"] >= 50:
						startend = set(range(int(result["Start1"]),int(result["End1"])))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
					# start1 to end2
					elif not pd.isna(result["Start1"]) and pd.isna(result["Start2"]) and pd.isna(result["End1"]) and not pd.isna(result["End2"]):
						start1toEnd2=True
						startend = set(range(int(result["Start1"]),int(result["Start1"])+surroundingLen))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
						#
						startend = set(range(int(result["End2"])-surroundingLen,int(result["End2"])))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
				else:
					resultList.extend(["-","-"])
			else:
				resultList.extend(["-","-"])
			#right flank on the same contig or only one flank
			if not pd.isna(result["Query.ID1"]) and pd.isna(result["Query.ID2"]) and not start1toEnd2:
				qr = queryRecords[result["Query.ID1"]]
				featureList = [feature for feature in qr.features if feature.type == 'CDS']
				scaffoldSize = len(qr)
				if len(featureList) > 0:
					if not pd.isna(result["Start2"]) and not pd.isna(result["End2"]) and result["End2"] < 50:
						resultList.extend(["-","Contig end."])
					if not pd.isna(result["Start2"]) and not pd.isna(result["End2"]) and result["Start2"] >= scaffoldSize - 50:
						resultList.extend(["-","Contig end."])
					elif not pd.isna(result["Start2"]) and not pd.isna(result["End2"]) and result["Start2"] <= scaffoldSize - 50:
						startend = set(range(int(result["Start2"]),int(result["End2"])))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
					else:
						resultList.extend(["-","-"])
				else:
					resultList.extend(["-","-"])
			#right flank on differnet contigs
			elif not pd.isna(result["Query.ID2"]) and not start1toEnd2:
				qr = queryRecords[result["Query.ID2"]]
				featureList = [feature for feature in qr.features if feature.type == 'CDS']
				scaffoldSize = len(qr)
				if len(featureList) > 0:
					if not pd.isna(result["Start2"]) and not pd.isna(result["End2"]) and result["End2"] < 50:
						resultList.extend(["-","Contig end."])
					if not pd.isna(result["Start2"]) and not pd.isna(result["End2"]) and result["Start2"] >= scaffoldSize - 50:
						resultList.extend(["-","Contig end."])
					elif not pd.isna(result["Start2"]) and not pd.isna(result["End2"]) and result["Start2"] <= scaffoldSize - 50:
						startend = set(range(int(result["Start2"]),int(result["End2"])))
						product = []
						locus_tag = []
						for  f in featureList:
							fspan = set(range(f.location.start,f.location.end))
							if 'locus_tag' in f.qualifiers.keys():
								ltag = f.qualifiers['locus_tag'][0]
							else:
								ltag = '-'
							if 'product' in f.qualifiers.keys():
								prod = f.qualifiers['product'][0]
							else:
								prod = '-'
							if fspan & startend:
								product.append(prod)
								locus_tag.append(ltag)
						if len(locus_tag) != 0:
							resultList.extend([' | '.join(locus_tag),' | '.join(product)])
						else:
							resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
				else:
					resultList.extend(["-","-"])
			elif not start1toEnd2:
				resultList.extend(["-","-"])
			#REF
			rr = refRecords[result["REF.ID1"]]
			featureList = [feature for feature in rr.features if feature.type == 'CDS']
			refSize = len(rr)
			scaffoldSize = len(rr)
			if len(featureList) > 0:
				#left flank
				if not pd.isna(result["REF.Start1"]) and not pd.isna(result["REF.End1"]) and result["REF.End1"] < 50:
					resultList.extend(["-","Contig end."])
				elif not pd.isna(result["REF.Start1"]) and not pd.isna(result["REF.End1"]) and result["REF.End1"] >= 50:
					startend = set(range(int(result["REF.Start1"]),int(result["REF.End1"])))
					product = []
					locus_tag = []
					for  f in featureList:
						fspan = set(range(f.location.start,f.location.end))
						if 'locus_tag' in f.qualifiers.keys():
							ltag = f.qualifiers['locus_tag'][0]
						else:
							ltag = '-'
						if 'product' in f.qualifiers.keys():
							prod = f.qualifiers['product'][0]
						else:
							prod = '-'
						if fspan & startend:
							product.append(prod)
							locus_tag.append(ltag)
					if len(locus_tag) != 0:
						resultList.extend([' | '.join(locus_tag),' | '.join(product)])
					else:
						resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
				#right flank
				if not pd.isna(result["REF.Start2"]) and not pd.isna(result["REF.End2"]) and result["REF.Start2"] >= scaffoldSize - 50:
					resultList.extend(["-","Contig end."])
				elif not pd.isna(result["REF.Start2"]) and not pd.isna(result["REF.End2"]) and result["REF.Start2"] <= scaffoldSize - 50:
					startend = set(range(int(result["REF.Start2"]),int(result["REF.End2"])))
					product = []
					locus_tag = []
					for  f in featureList:
						fspan = set(range(f.location.start,f.location.end))
						if 'locus_tag' in f.qualifiers.keys():
							ltag = f.qualifiers['locus_tag'][0]
						else:
							ltag = '-'
						if 'product' in f.qualifiers.keys():
							prod = f.qualifiers['product'][0]
						else:
							prod = '-'
						if fspan & startend:
							product.append(prod)
							locus_tag.append(ltag)
					if len(locus_tag) != 0:
						resultList.extend([' | '.join(locus_tag),' | '.join(product)])
					else:
						resultList.extend(["-","Intergenic/Pseudo/Repeaed Region"])
			else:
				resultList.extend(["-","-","-","-"])
			resultList.extend([refSize])
		FinalResults.loc[i] = resultList
		i = i +1
	return FinalResults

#Plot functions
def niter(iterable, n):
	"""
	Function that returns an n-element iterator, i.e.
	sub-lists of a list that are max. n elements long.
	"""
	pos = 0
	while pos < len(iterable):
		yield iterable[pos:pos+n]
		pos += n

class MyCustomTranslator(BiopythonTranslator):
	"""Custom translator implementing the following theme:
	- Show only CDS.
	    """
	def compute_feature_color(self, feature):
		if feature.type == "CDS":
			return "darkgray"
		elif feature.type == "DetectedIS":
			return "red"
		elif feature.type == "QIF" or feature.type == "RA":
			return "teal"
		elif feature.type == "QIF2" or feature.type == "RA2":
			return "mediumaquamarine"
		elif feature.type == "RA3":
			return "turquoise"
		else:
			return "darkgray"
	# def compute_feature_label(self, feature):
	# 	if feature.type == 'restriction_site':
	# 		return None
	# 	elif feature.type == "DetectedIS":
	# 		return "IS"
	# 	else:
	# 		return BiopythonTranslator.compute_feature_label(self, feature)
	#
	def compute_filtered_features(self, features):
		"""Do not display source."""
		return [
			feature for feature in features
			if (feature.type == "CDS" or feature.type == 'QIF' or feature.type == 'QIF2' or feature.type == 'RA' or feature.type == 'RA2' or feature.type == 'RA3' or feature.type == 'DetectedIS')
		]

def add_new_feature(record,start,end,type,product):
	my_start_pos = SeqFeature.ExactPosition(start)
	my_end_pos = SeqFeature.ExactPosition(end)
	my_feature_location = FeatureLocation(my_start_pos,my_end_pos)
	my_feature_type = type
	my_feature = seqFeature(my_feature_location,type=my_feature_type)
	my_feature.qualifiers={'product':product}
	record.features.append(my_feature)
	return record

def modifyWarpingFeatures(record):
	if record.features:
		for feature in record.features:
			if len(feature.location.parts) <= 1:
				pass
			elif len(feature.location.parts) > 1:
				if feature.type == 'CDS':
					for part in feature.location.parts:
						feature_start = part.start
						feature_end = part.end
						product = feature.qualifiers['product'][0]
						record = add_new_feature(record,feature_start,feature_end,"CDS",product)
					record.features.remove(feature)
	return record

def recordsDict(filename):
	recordsDictionary={}
	for record in SeqIO.parse(filename, "genbank"):
		recordsDictionary[record.id] = modifyWarpingFeatures(record)
	return 	recordsDictionary

def plot_funcs(ISs, max_row, page):
	"""
	Function that plots all given ISs in multiple figs, to accomodate max_col * max_row
	IS per page.
	"""
	plt.rcParams.update({'figure.max_open_warning': 0})
	##amount of ISs to put in one page
	N = max_row
	##created figures go here
	figs = []
	##looping through ISs N at a time:
	##figure and subplots
	n = len(ISs)
	fig, axes = plt.subplots(N,3,figsize=(8.27*3,11.7*3),dpi=150)
	fig.suptitle('Differential IS insertions graphic report. Page '+ str(page), fontsize=24)
	if  n < N:
		for ax in axes[n:,:].flat:
			ax.remove()
	i=0
	##plotting functions
	for IS in ISs.iterrows():
		#case IS on query or IS on ref...
		if IS[1]['Description'] == 'IS only present in query Strain.':
			ID1 = IS[1]['Query.ID1']
			start1 = int(IS[1]['Start1'] - 3000)
			# end1 = IS[1]['End1']
			# start2 = IS[1]['Start2']
			end2 = int(IS[1]['End2'] + 3000)
			ISstart = int(IS[1]['ISstart'])
			ISend = int(IS[1]['ISend'])
			name = IS[1]['Query.ID1'] + " - " + IS[1]['ISID'] +" Start: " + str(ISstart) + " End: " + str(ISend)
			#
			record = queryRecords[ID1]
			recordLength = len(record)
			if start1 < 1:
				start1 = 1
			if end2 > recordLength:
				end2 = recordLength
			record = record.upper()
			record = add_new_feature(record,ISstart,ISend,"DetectedIS",IS[1]['ISID'])
			record = add_new_feature(record,int(IS[1]['Start1']),int(IS[1]['End1']),"QIF","Flank 1")
			record = add_new_feature(record,int(IS[1]['Start2']),int(IS[1]['End2']),"QIF2","Flank 2")
			graphic_record = MyCustomTranslator().translate_record(record)
			graphic_record.ticks_resolution = 2000
			croped_graphic_record = graphic_record.crop((start1, end2))
			ax1 = axes[i][0]
			ax2 = axes[i][1]
			ax3 = axes[i][2]
			ax, _ = croped_graphic_record.plot(ax=ax1, figure_width=7, strand_in_label_threshold=7)
			ax1.set_title(name,loc='left',y=-0.2)
			#Ref plots REF1/REF2
			Rstart1 = IS[1]['REF.Start1']
			if 	Rstart1 == Rstart1:
				REFID1 = IS[1]['REF.ID1']
				REFID2 = IS[1]['REF.ID2']
				Rstart1 = int(IS[1]['REF.Start1'])
				Rend1 = IS[1]['REF.End1']
				Rstart2 = IS[1]['REF.Start2']
				Rend2 = IS[1]['REF.End2']
				if Rstart2 == Rstart2:
					Rend1 = int(IS[1]['REF.End1'])
					Rstart2 = int(IS[1]['REF.Start2'])
					Rend2 = int(IS[1]['REF.End2'])
					record1 = refRecords[REFID1]
					recordLength1 = len(record1)
					#REF.ID2 != nan
					if REFID2 == REFID2:
						if Rstart1 - 3000 < 1:
							Rstart1 = 1
						else:
							Rstart1 = Rstart1 - 3000
						if Rend1 + 3000 > recordLength1:
							Rend1 = recordLength1
						else:
							Rend1 = Rend1 + 3000
						record2 = refRecords[REFID2]
						recordLength2 = len(record2)
						if Rstart2 - 3000 < 1:
							Rstart2 = 1
						else:
							Rstart2 = Rstart2 - 3000
						if Rend2 + 3000 > recordLength2:
							Rend2 = recordLength2
						else:
							Rend2 = Rend2 + 3000
						#here two plots
						record1 = record1.upper()
						record2 = record2.upper()
						record1 = add_new_feature(record1,int(IS[1]['REF.Start1']),int(IS[1]['REF.End1']),"RA","Match 1")
						record2 = add_new_feature(record2,int(IS[1]['REF.Start2']),int(IS[1]['REF.End2']),"RA2","Match 2")
						graphic_record1 = MyCustomTranslator().translate_record(record1)
						graphic_record1.ticks_resolution = 2000
						croped_graphic_record1 = graphic_record1.crop((Rstart1, Rend1))
						ax, _ = croped_graphic_record1.plot(ax=ax2, figure_width=7, strand_in_label_threshold=7)
						name = IS[1]['REF.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend1)
						ax2.set_title(name,loc='left',y=-0.2)
						graphic_record2 = MyCustomTranslator().translate_record(record2)
						graphic_record2.ticks_resolution = 2000
						croped_graphic_record2 = graphic_record2.crop((Rstart2, Rend2))
						ax, _ = croped_graphic_record2.plot(ax=ax3, figure_width=7, strand_in_label_threshold=7)
						name = IS[1]['REF.ID2'] + " Start: " + str(Rstart2) + " End: " + str(Rend2)
						ax3.set_title(name,loc='left',y=-0.2)
					else:
						if (abs(Rstart1 - Rstart2) < 10000) or (abs(Rstart1 - Rstart2) < 3*shift):
							# here 1 plot instead of 2
							if Rstart1 > Rstart2:
								Rstart = Rstart2
								Rend = Rend1
							else:
								Rstart = Rstart1
								Rend = Rend2
							if Rstart - 3000 < 1:
								Rstart = 1
							else:
								Rstart = Rstart - 3000
							if Rend + 3000 > recordLength1:
								Rend = recordLength1
							else:
								Rend = Rend + 3000
							record1 = record1.upper()
							record1 = add_new_feature(record1,int(IS[1]['REF.Start1']),int(IS[1]['REF.End1']),"RA","Match 1")
							record1 = add_new_feature(record1,int(IS[1]['REF.Start2']),int(IS[1]['REF.End2']),"RA2","Match 2")
							graphic_record1 = MyCustomTranslator().translate_record(record1)
							graphic_record1.ticks_resolution = 2000
							croped_graphic_record1 = graphic_record1.crop((Rstart, Rend))
							ax, _ = croped_graphic_record1.plot(ax=ax2, figure_width=7, strand_in_label_threshold=7)
							name = IS[1]['REF.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend1)
							ax2.set_title(name,loc='left',y=-0.2)
							ax3.remove()
						else:
							#here two plots
							if Rstart1 - 3000 < 1:
								Rstart1 = 1
							else:
								Rstart1 = Rstart1 - 3000
							if Rend1 + 3000 > recordLength1:
								Rend1 = recordLength1
							else:
								Rend1 = Rend1 + 3000
							if Rstart2 - 3000 < 1:
								Rstart2 = 1
							else:
								Rstart2 = Rstart2 - 3000
							if Rend2 + 3000 > recordLength1:
								Rend2 = recordLength1
							else:
								Rend2 = Rend2 + 3000
							record1 = record1.upper()
							record1 = add_new_feature(record1,int(IS[1]['REF.Start1']),int(IS[1]['REF.End1']),"RA","Match 1")
							record1 = add_new_feature(record1,int(IS[1]['REF.Start2']),int(IS[1]['REF.End2']),"RA2","Match 2")
							graphic_record1 = MyCustomTranslator().translate_record(record1)
							graphic_record1.ticks_resolution = 2000
							croped_graphic_record1 = graphic_record1.crop((Rstart1, Rend1))
							ax, _ = croped_graphic_record1.plot(ax=ax2, figure_width=7, strand_in_label_threshold=7)
							name = IS[1]['REF.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend1)
							ax2.set_title(name,loc='left',y=-0.2)
							croped_graphic_record2 = graphic_record1.crop((Rstart2, Rend2))
							ax, _ = croped_graphic_record2.plot(ax=ax3, figure_width=7, strand_in_label_threshold=7)
							name = IS[1]['REF.ID1'] + " Start: " + str(Rstart2) + " End: " + str(Rend2)
							ax3.set_title(name,loc='left',y=-0.2)
				elif Rend1 == Rend1:
					Rend1 = int(IS[1]['REF.End1'])
					record1 = refRecords[REFID1]
					recordLength1 = len(record1)
					if Rstart1 - 3000 < 1:
						Rstart1 = 1
					else:
						Rstart1 = Rstart1 - 3000
					if Rend1 + 3000 > recordLength1:
						Rend1 = recordLength1
					else:
						Rend1 = Rend1 + 3000
					# Here only 1 plot
					record1 = record1.upper()
					record1 = add_new_feature(record1,int(IS[1]['REF.Start1']),int(IS[1]['REF.End1']),"RA3","Match")
					graphic_record = MyCustomTranslator().translate_record(record1)
					graphic_record.ticks_resolution = 2000
					croped_graphic_record = graphic_record.crop((Rstart1, Rend1))
					ax, _ = croped_graphic_record.plot(ax=ax2, figure_width=7, strand_in_label_threshold=7)
					name = IS[1]['REF.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend1)
					ax2.set_title(name,loc='left',y=-0.2)
					ax3.remove()
				elif Rend2 == Rend2 :
					Rend2 = int(IS[1]['REF.End2'])
					record1 = refRecords[REFID1]
					recordLength1 = len(record1)
					if Rstart1 - 3000 < 1:
						Rstart1 = 1
					else:
						Rstart1 = Rstart1 - 3000
					if Rend2 + 3000 > recordLength1:
						Rend2 = recordLength1
					else:
						Rend2 = Rend2 + 3000
					# Here only 1 plot
					record1 = record1.upper()
					record1 = add_new_feature(record1,int(IS[1]['REF.Start1']),int(IS[1]['REF.End2']),"RA3","Match 1-2")
					graphic_record = MyCustomTranslator().translate_record(record1)
					graphic_record.ticks_resolution = 2000
					croped_graphic_record = graphic_record.crop((Rstart1, Rend2))
					ax, _ = croped_graphic_record.plot(ax=ax2, figure_width=7, strand_in_label_threshold=7)
					name = IS[1]['REF.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend2)
					ax2.set_title(name,loc='left',y=-0.2)
					ax3.remove()
			else:
				ax2.remove()
				ax3.remove()
		# Ref
		elif IS[1]['Description'] == 'IS only present in ref Strain.':
			ID1 = IS[1]['REF.ID1']
			start1 = int(IS[1]['REF.Start1'] - 3000)
			# end1 = IS[1]['End1']
			# start2 = IS[1]['Start2']
			end2 = int(IS[1]['REF.End2'] + 3000)
			ISstart = int(IS[1]['ISstart'])
			ISend = int(IS[1]['ISend'])
			name = IS[1]['REF.ID1'] + " - " + IS[1]['ISID'] +" Start: " + str(ISstart) + " End: " + str(ISend)
			#
			record = refRecords[ID1]
			recordLength = len(record)
			if start1 < 1:
				start1 = 1
			if end2 > recordLength:
				end2 = recordLength
			record = record.upper()
			record = add_new_feature(record,ISstart,ISend,"DetectedIS",IS[1]['ISID'])
			record = add_new_feature(record,int(IS[1]['REF.Start1']),int(IS[1]['REF.End1']),"QIF","Flank 1")
			record = add_new_feature(record,int(IS[1]['REF.Start2']),int(IS[1]['REF.End2']),"QIF2","Flank 2")
			graphic_record = MyCustomTranslator().translate_record(record)
			graphic_record.ticks_resolution = 2000
			croped_graphic_record = graphic_record.crop((start1, end2))
			ax1 = axes[i][0]
			ax2 = axes[i][1]
			ax3 = axes[i][2]
			ax, _ = croped_graphic_record.plot(ax=ax1, figure_width=7, strand_in_label_threshold=7)
			ax1.set_title(name,loc='left',y=-0.2)
			#Ref plots REF1/REF2
			Rstart1 = IS[1]['Start1']
			if 	Rstart1 == Rstart1:
				REFID1 = IS[1]['Query.ID1']
				REFID2 = IS[1]['Query.ID2']
				Rstart1 = int(IS[1]['Start1'])
				Rend1 = IS[1]['End1']
				Rstart2 = IS[1]['Start2']
				Rend2 = IS[1]['End2']
				if Rstart2 == Rstart2:
					Rend1 = int(IS[1]['End1'])
					Rstart2 = int(IS[1]['Start2'])
					Rend2 = int(IS[1]['End2'])
					record1 = queryRecords[REFID1]
					recordLength1 = len(record1)
					#REF.ID2 != nan
					if REFID2 == REFID2:
						if Rstart1 - 3000 < 1:
							Rstart1 = 1
						else:
							Rstart1 = Rstart1 - 3000
						if Rend1 + 3000 > recordLength1:
							Rend1 = recordLength1
						else:
							Rend1 = Rend1 + 3000
						record2 = queryRecords[REFID2]
						recordLength2 = len(record2)
						if Rstart2 - 3000 < 1:
							Rstart2 = 1
						else:
							Rstart2 = Rstart2 - 3000
						if Rend2 + 3000 > recordLength2:
							Rend2 = recordLength2
						else:
							Rend2 = Rend2 + 3000
						#here two plots
						record1 = record1.upper()
						record2 = record2.upper()
						record1 = add_new_feature(record1,int(IS[1]['Start1']),int(IS[1]['End1']),"RA","Match 1")
						record2 = add_new_feature(record2,int(IS[1]['Start2']),int(IS[1]['End2']),"RA2","Match 2")
						graphic_record1 = MyCustomTranslator().translate_record(record1)
						graphic_record1.ticks_resolution = 2000
						croped_graphic_record1 = graphic_record1.crop((Rstart1, Rend1))
						ax, _ = croped_graphic_record1.plot(ax=ax2, figure_width=7, strand_in_label_threshold=7)
						name = IS[1]['Query.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend1)
						ax2.set_title(name,loc='left',y=-0.2)
						graphic_record2 = MyCustomTranslator().translate_record(record2)
						graphic_record2.ticks_resolution = 2000
						croped_graphic_record2 = graphic_record2.crop((Rstart2, Rend2))
						ax, _ = croped_graphic_record2.plot(ax=ax3, figure_width=7, strand_in_label_threshold=7)
						name = IS[1]['Query.ID2'] + " Start: " + str(Rstart2) + " End: " + str(Rend2)
						ax3.set_title(name,loc='left',y=-0.2)
					else:
						if (abs(Rstart1 - Rstart2) < 10000) or (abs(Rstart1 - Rstart2) < 3*shift):
							# here 1 plot instead of 2
							if Rstart1 > Rstart2:
								Rstart = Rstart2
								Rend = Rend1
							else:
								Rstart = Rstart1
								Rend = Rend2
							if Rstart - 3000 < 1:
								Rstart = 1
							else:
								Rstart = Rstart - 3000
							if Rend + 3000 > recordLength1:
								Rend = recordLength1
							else:
								Rend = Rend + 3000
							record1 = record1.upper()
							record1 = add_new_feature(record1,int(IS[1]['Start1']),int(IS[1]['End1']),"RA","Match 1")
							record1 = add_new_feature(record1,int(IS[1]['Start2']),int(IS[1]['End2']),"RA2","Match 2")
							graphic_record1 = MyCustomTranslator().translate_record(record1)
							graphic_record1.ticks_resolution = 2000
							croped_graphic_record1 = graphic_record1.crop((Rstart, Rend))
							ax, _ = croped_graphic_record1.plot(ax=ax2, figure_width=7, strand_in_label_threshold=7)
							name = IS[1]['Query.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend1)
							ax2.set_title(name,loc='left',y=-0.2)
							ax3.remove()
						else:
							#here two plots
							if Rstart1 - 3000 < 1:
								Rstart1 = 1
							else:
								Rstart1 = Rstart1 - 3000
							if Rend1 + 3000 > recordLength1:
								Rend1 = recordLength1
							else:
								Rend1 = Rend1 + 3000
							if Rstart2 - 3000 < 1:
								Rstart2 = 1
							else:
								Rstart2 = Rstart2 - 3000
							if Rend2 + 3000 > recordLength1:
								Rend2 = recordLength1
							else:
								Rend2 = Rend2 + 3000
							record1 = record1.upper()
							record1 = add_new_feature(record1,int(IS[1]['Start1']),int(IS[1]['End1']),"RA","Match 1")
							record1 = add_new_feature(record1,int(IS[1]['Start2']),int(IS[1]['End2']),"RA2","Match 2")
							graphic_record1 = MyCustomTranslator().translate_record(record1)
							graphic_record1.ticks_resolution = 2000
							croped_graphic_record1 = graphic_record1.crop((Rstart1, Rend1))
							ax, _ = croped_graphic_record1.plot(ax=ax2,figure_width=7, strand_in_label_threshold=7)
							name = IS[1]['Query.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend1)
							ax2.set_title(name,loc='left',y=-0.2)
							croped_graphic_record2 = graphic_record1.crop((Rstart2, Rend2))
							ax, _ = croped_graphic_record2.plot(ax=ax3,figure_width=7, strand_in_label_threshold=7)
							name = IS[1]['Query.ID1'] + " Start: " + str(Rstart2) + " End: " + str(Rend2)
							ax3.set_title(name,loc='left',y=-0.2)
				elif Rend1 == Rend1:
					Rend1 = int(IS[1]['End1'])
					record1 = queryRecords[REFID1]
					recordLength1 = len(record1)
					if Rstart1 - 3000 < 1:
						Rstart1 = 1
					else:
						Rstart1 = Rstart1 - 3000
					if Rend1 + 3000 > recordLength1:
						Rend1 = recordLength1
					else:
						Rend1 = Rend1 + 3000
					# Here only 1 plot
					record1 = record1.upper()
					record1 = add_new_feature(record1,int(IS[1]['Start1']),int(IS[1]['End1']),"RA3","Match")
					graphic_record = MyCustomTranslator().translate_record(record1)
					graphic_record.ticks_resolution = 2000
					croped_graphic_record = graphic_record.crop((Rstart1, Rend1))
					ax, _ = croped_graphic_record.plot(ax=ax2, figure_width=7, strand_in_label_threshold=7)
					name = IS[1]['Query.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend1)
					ax2.set_title(name,loc='left',y=-0.2)
					ax3.remove()
				elif Rend2 == Rend2 :
					Rend2 = int(IS[1]['End2'])
					record1 = queryRecords[REFID1]
					recordLength1 = len(record1)
					if Rstart1 - 3000 < 1:
						Rstart1 = 1
					else:
						Rstart1 = Rstart1 - 3000
					if Rend2 + 3000 > recordLength1:
						Rend2 = recordLength1
					else:
						Rend2 = Rend2 + 3000
					# Here only 1 plot
					record1 = record1.upper()
					record1 = add_new_feature(record1,int(IS[1]['Start1']),int(IS[1]['End2']),"RA3","Match 1-2")
					graphic_record = MyCustomTranslator().translate_record(record1)
					graphic_record.ticks_resolution = 2000
					croped_graphic_record = graphic_record.crop((Rstart1, Rend2))
					ax, _ = croped_graphic_record.plot(ax=ax2, figure_width=7, strand_in_label_threshold=7)
					name = IS[1]['Query.ID1'] + " Start: " + str(Rstart1) + " End: " + str(Rend2)
					ax2.set_title(name,loc='left',y=-0.2)
					ax3.remove()
			else:
				ax2.remove()
				ax3.remove()
		i+=1
	fig.tight_layout(pad=7.5)#pad=0.5
	figs.append(fig)
	return figs

def checkFiles(filePath):
	if os.path.exists(filePath) and os.path.getsize(filePath) > 0:
		pass
	else:
		print('Error: '+ str(filePath) + ' file missing or empty.')
		exit()

#********************************************************************************************************
#********************************************************************************************************
#********************************************************************************************************
# Main
if __name__ == "__main__":
	#Presentation
	printSoftName()
	#Argument parsing
	args = parseArgs()
	#Check ISFile
	if args.isfile:
		isfile = args.isfile
		checkFiles(isfile)
	if (args.qfiles and args.rfiles and args.qacc and args.racc):
		print("Using local files only.\n")
		mode = 0
		qfiles=args.qfiles[0]
		rfiles=args.rfiles[0]
		for filePath in qfiles:
			checkFiles(filePath)
		for filePath in rfiles:
			checkFiles(filePath)
	elif args.qfiles and args.rfiles:
		mode = 0
		qfiles=args.qfiles[0]
		rfiles=args.rfiles[0]
		for filePath in qfiles:
			checkFiles(filePath)
		for filePath in rfiles:
			checkFiles(filePath)
	elif args.qacc and args.racc:
		if not args.yourEmail:
			sys.stdout.write("You must supply an email address for accession number mode.\n")
			sys.exit()
		elif re.match('[^@]+@[^@]+\..{1,3}$',args.yourEmail):
			Entrez.email = args.yourEmail
			mode = 1
			qacc=args.qacc[0]
			racc=args.racc[0]
		else:
			print("Incorrect email address. Email adress shoud be: xxxxx@xxxx.xxx .\n")
			sys.exit()
	else:
		sys.stdout.write("You must supply either the reference and query genomes, or NCBI accession numbers.\nRun ./ISCompare.py -h for help.\n\n")
		sys.exit()
	#output folder creation
	if args.output:
		outputDir=args.output
		outputDirProvided=True
		if not os.path.exists(outputDir):
			os.mkdir(outputDir)
			print("Directory " , outputDir ,  " created.")
		else:
			print("Directory " , outputDir ,  " already exists.")
	else:
		outputDir='results'
		outputDirProvided=False
		if not os.path.exists(outputDir):
			os.mkdir(outputDir)
			print("Using default folder, 'results'.")
		else:
			print("Directory " , outputDir ,  " already exists.")
	#
	#minlength
	if args.minLength:
		minLength = int(args.minLength)
	else:
		minLength = 50
	#ISdiff
	if args.ISdiff:
		ISdiff = int(args.ISdiff)
	else:
		ISdiff = 50
	#scaffoldDiff
	if args.scaffoldDiff:
		ISdiff = int(args.scaffoldDiff)
	else:
		scaffoldDiff = 20
	#minAlnLength
	if args.minAlnLength:
		minAlnLength = int(args.minAlnLength)
	else:
		minAlnLength = 50
	#surroundingLen
	if args.surroundingLen:
		surroundingLen = int(args.surroundingLen)
	else:
		surroundingLen = 500
	#surroundingLen2
	if args.surroundingLen2:
		surroundingLen2 = int(args.surroundingLen2)
	else:
		surroundingLen2 = 500
	if args.shift:
		shift = int(args.shift)
	else:
		shift = 0
	#evalue
	if args.evalue:
		evalue = args.evalue
	else:
		evalue = '0.0000000001'
	#
	#Open log files
	logFilePath=outputDir+'/ISCompare.log'
	if os.path.exists(outputDir+'/ISCompare.log'):
		os.remove(logFilePath)
	#
	log = open(logFilePath, 'w+')
	log.write('ISCompare log file.\n')
	log.flush()
	#
	#Mode 0, files
	if mode == 0:
		log.write('Moving files to results folder.\n')
		log.flush()
		print("Query: " + str(qfiles) + " Ref: " + str(rfiles))
		print("Copying files to " + outputDir + " folder, and extracting fasta files from genbank.")
		newPath = outputDir+'/query.gb'
		f=open(newPath, "w+")
		for filewithPath in qfiles:
			file = os.path.basename(filewithPath)
			f2=open(filewithPath,"r")
			shutil.copyfileobj(f2,f)
		f.close()
		f2.close()
		qseq = "query"
		extractFastaFromGB(qseq)
		print(qseq + ".fasta generated.")
		newPath = outputDir+'/ref.gb'
		f=open(newPath, "w+")
		for filewithPath in rfiles:
			file = os.path.basename(filewithPath)
			f2=open(filewithPath,"r")
			shutil.copyfileobj(f2,f)
		rseq = "ref"
		f.close()
		f2.close()
		extractFastaFromGB(rseq)
		print(rseq + ".fasta generated.")
	#
	#Mode 1, download files from NCBI
	elif mode == 1:
		log.write('Downloading genomes from NCBI.\n')
		log.flush()
		print(Entrez.email)
		print("Downloading Query: "+str(qacc)+" and Ref: "+str(racc)+" sequences.")
		qseq = getRefSeqFromNCBIbyACC(qacc,"query")
		rseq = getRefSeqFromNCBIbyACC(racc,"ref")
	#
	#ISFile
	if args.isfile:
		pass
	elif args.isscan:
		print('Retrieving IS sequences from ISFinder database. \nPlease cite: Siguier P. et al. (2006)\nISfinder: the reference centre for bacterial insertion sequences. \nNucleic Acids Res. 34: D32-D36 \nISFinder database URL: http://www-is.biotoul.fr.')
		queryFileFasta = os.path.join(outputDir,qseq+".fasta")
		refFileFasta = os.path.join(outputDir,rseq+".fasta")
		isfile = os.path.join(outputDir,"IS.fasta")
		try:
			IStable = ISFinderBlast.BlastISatISFinder(queryFileFasta)
			IStable = concatPdTables(IStable,ISFinderBlast.BlastISatISFinder(refFileFasta))
			IStable = IStable.groupby('Sequences producing significant alignments',as_index=False).first()
		except:
			print('Something went wrong with the ISFinder blast service.\nTry to run the program with -i IS.fasta option.')
			exit()
		try:	
			ISFinderBlast.getISsequencesFromISFinderDB(isfile,IStable)
			print(str(len(IStable))+' IS sequences downloaded from ISFinder database.')
		except:
			print('Something went wrong while downloading the IS sequences.\nTry to run the program with -i IS.fasta option.')
			exit()
	else:
		print('A fasta file with the IS sequences to look for must be provided, \nor you could run with -I option to search for ISs using ISFinder server.')
		exit()
	#
	#************************************************************************
	#**********************    Query to Ref  ********************************
	#************************************************************************
	#************************************************************************
	queryFileFasta = os.path.join(outputDir,qseq+".fasta")
	totalqueryIS, totalstatsQuery = findISsOnQuery(rseq,qseq,ISdiff,surroundingLen,queryFileFasta)
	print(str(len(totalqueryIS)) + " ISs found on the query genome...")
	log.write(str(len(totalqueryIS)) + " ISs found on the query genome...")
	if len(totalqueryIS) == 0:
		print("The selected Insertion sequences were not found on the query genome...")
		log.write("The selected Insertion sequences were not found on the query genome...")
	log.flush()
	nscaffolds, maxLenScaffolds = getNumberOfScaffoldsFromGB(qseq)
	if args.remove:
		print("Step1: Removing identical scaffolds: Query --> Reference ...")
		log.write("Step1: Removing identical scaffolds: Query --> Reference ...")
		log.flush()
		identicalScaffolds = identifyIdenticalScaffolds(rseq, qseq, maxLenScaffolds,scaffoldDiff)
		removeIdenticalScaffolds(queryFileFasta, identicalScaffolds)
		print(str(len(identicalScaffolds)) + " identical scaffolds removed from the analysis...")
		log.write(str(len(identicalScaffolds)) + " identical scaffolds removed from the analysis...")
		queryFileFasta = os.path.join(outputDir,"query_reduced.fasta")
	else:
		print("Step1: Removing identical scaffolds: Query --> Reference // OMITTED by user: -x opt ...")
		log.write("Step1: Removing identical scaffolds: Query --> Reference // OMITTED by user: -x opt ...")
		log.flush()
		queryFileFasta = os.path.join(outputDir,"query.fasta")
	#
	#************************************************************************
	print("Step2: Testing for new IS on the query...")
	log.write("Step2: Testing for new IS on the query...")
	log.flush()

	queryIS, statsQuery = findISsOnQuery(rseq,qseq,ISdiff,surroundingLen,queryFileFasta)
	# Consecutive iSs, get a list of the consecutive ISs on queryIS
	ConsecutiveISQuery = queryIS[queryIS['ConsecutiveIS']]
	MiddleMissingIS, consecutiveIS= testMiddleIS(rseq, qseq, queryIS)
	MiddleMissingISID = MiddleMissingIS['qseqid']
	MiddleMissingIS = MiddleMissingIS[['qseqid','sseqid','start', 'end', 'sstart', 'send', 'Scaffold', 'Observations','IS_Match_Type']]
	MiddleMissingIS['sseqid_first'] = MiddleMissingIS['sseqid']
	MiddleMissingIS['sseqid_last'] = np.nan
	MiddleMissingIS['sstart_first'] = MiddleMissingIS['sstart']
	MiddleMissingIS['sstart_last'] = MiddleMissingIS['send_first'] = np.nan
	MiddleMissingIS['send_last'] = MiddleMissingIS['send']
	MiddleMissingIS = MiddleMissingIS[['qseqid','sseqid_first','sseqid_last', 'start', 'end', 'sstart_first','send_first','sstart_last','send_last', 'Scaffold', 'Observations', 'IS_Match_Type']]
	#
	missingIS = testIS(rseq, qseq,MiddleMissingISID,consecutiveIS)
	missingIS = missingIS[['qseqid','sseqid_first','sseqid_last', 'start', 'end', 'sstart_first','send_first','sstart_last','send_last', 'Scaffold', 'Observations', 'IS_Match_Type']]
	#join tables
	missingIS = concatPdTables(missingIS,MiddleMissingIS)
	missingIS["Description"] = "IS only present in " + qseq + " Strain."
	if len(missingIS)>0:
		missingIS["ISstart"] = missingIS.apply(lambda x: re.sub(r".*\|ISstart:([0-9]{1,10}).*$","\\1",x["qseqid"]), axis = 1)
		missingIS["ISend"] = missingIS.apply(lambda x: re.sub(r".*\|ISend:([0-9]{1,10}).*$","\\1",x["qseqid"]), axis = 1)
		missingIS["ISID"]  = missingIS.apply(lambda x: re.sub(r'^.*\|IS:([^|]*)\|.*$','\\1',x["qseqid"]),axis=1)
		#
		missingIS["Lend"]  = missingIS.apply(lambda x: int(re.sub(r'^.*\|LEFT_End:([0-9]{1,10})\|.*$','\\1',x["qseqid"])),axis=1)
		missingIS["Rstart"]  = missingIS.apply(lambda x: int(re.sub(r'^.*\|RIGHT_start:([0-9]{1,10})\|.*$','\\1',x["qseqid"])),axis=1)
	else:
		missingIS["ISstart"] = np.nan
		missingIS["ISend"] = np.nan
		missingIS["ISID"]  = np.nan
		#
		missingIS["Lend"]  = np.nan
		missingIS["Rstart"]  = np.nan
	#
	missingIS['Scaffold2'] = np.nan
	results = missingIS[["ISID","ISstart","ISend","Scaffold","start","Lend","Scaffold2","Rstart","end","sseqid_first","sstart_first","send_first","sseqid_last","sstart_last","send_last","Description","Observations","IS_Match_Type"]]
	columns = ["ISID","ISstart","ISend","Query.ID1","Start1","End1","Query.ID2","Start2","End2", "REF.ID1","REF.Start1","REF.End1","REF.ID2","REF.Start2","REF.End2","Description","Observations","IS_Match_Type"]
	results.columns = columns
	#************************************************************************
	#**********************    Ref to query   *******************************
	#************************************************************************
	#************************************************************************
	refFileFasta = os.path.join(outputDir,rseq+".fasta")
	totalrefIS, totalstatsRef = findISsOnQuery(qseq,rseq,ISdiff,surroundingLen,refFileFasta)
	print(str(len(totalrefIS)) + " ISs found on the reference genome...")
	log.write(str(len(totalrefIS)) + " ISs found on the reference genome...")
	if len(totalrefIS) == 0:
		print("The selected Insertion sequences were not found on the reference genome...")
		log.write("The selected Insertion sequences were not found on the reference genome...")
	log.flush()
	nscaffoldsRef, maxLenScaffoldsRef = getNumberOfScaffoldsFromGB(rseq)
	if args.remove:
		print("Step3: Removing identical scaffolds: Reference --> Query ...")
		log.write("Step3: Removing identical scaffolds: Reference --> Query ...")
		log.flush()
		identicalScaffoldsRef = identifyIdenticalScaffolds(qseq, rseq, maxLenScaffoldsRef,scaffoldDiff)
		removeIdenticalScaffolds(refFileFasta, identicalScaffoldsRef)
		print(str(len(identicalScaffoldsRef)) + " identical scaffolds removed from the analysis...")
		log.write(str(len(identicalScaffoldsRef)) + " identical scaffolds removed from the analysis...")
		refFileFasta = os.path.join(outputDir,"ref_reduced.fasta")
	else:
		print("Step3: Removing identical scaffolds: Reference --> Query // OMITTED by user: -x opt ...")
		log.write("Step3: Removing identical scaffolds: Reference --> Query // OMITTED by user: -x opt ...")
		log.flush()
		refFileFasta = os.path.join(outputDir,"ref.fasta")
	#************************************************************************
	#************************************************************************
	print("Step4: Testing for IS missing on the query...")
	log.write("Step4: Testing for IS missing on the query...")
	log.flush()
	refIS, statsRef = findISsOnQuery(qseq,rseq,ISdiff,surroundingLen,refFileFasta)
	# Consecutive iSs, get a list of the consecutive ISs on queryIS
	ConsecutiveISRef = refIS[refIS['ConsecutiveIS']]
	#
	MiddleMissingIS, consecutiveIS = testMiddleIS(qseq, rseq, refIS)
	MiddleMissingISID = MiddleMissingIS['qseqid']
	MiddleMissingIS = MiddleMissingIS[['qseqid','sseqid','start', 'end', 'sstart', 'send', 'Scaffold', 'Observations','IS_Match_Type']]
	MiddleMissingIS['sseqid_first'] = MiddleMissingIS['sseqid']
	MiddleMissingIS['sseqid_last'] = np.nan
	MiddleMissingIS['sstart_first'] = MiddleMissingIS['sstart']
	MiddleMissingIS['sstart_last'] = MiddleMissingIS['send_first'] = np.nan
	MiddleMissingIS['send_last'] = MiddleMissingIS['send']
	missingIS = testIS(qseq,rseq,MiddleMissingISID,consecutiveIS)
	missingIS = missingIS[['qseqid','sseqid_first','sseqid_last', 'start', 'end', 'sstart_first','send_first','sstart_last','send_last', 'Scaffold', 'Observations', 'IS_Match_Type']]
	# Join tables
	missingIS = concatPdTables(missingIS,MiddleMissingIS)
	missingIS["Description"] = "IS only present in " + rseq + " Strain."
	if len(missingIS)>0:
		missingIS["ISstart"] = missingIS.apply(lambda x: re.sub(r".*\|ISstart:([0-9]{1,10}).*$","\\1",x["qseqid"]), axis = 1)
		missingIS["ISend"] = missingIS.apply(lambda x: re.sub(r".*\|ISend:([0-9]{1,10}).*$","\\1",x["qseqid"]), axis = 1)
		missingIS["ISID"]  = missingIS.apply(lambda x: re.sub(r'^.*\|IS:([^|]*)\|.*$','\\1',x["qseqid"]),axis=1)
		#
		missingIS["Lend"]  = missingIS.apply(lambda x: int(re.sub(r'^.*\|LEFT_End:([0-9]{1,10})\|.*$','\\1',x["qseqid"])),axis=1)
		missingIS["Rstart"]  = missingIS.apply(lambda x: int(re.sub(r'^.*\|RIGHT_start:([0-9]{1,10})\|.*$','\\1',x["qseqid"])),axis=1)
	else:
		missingIS["ISstart"] = np.nan
		missingIS["ISend"] = np.nan
		missingIS["ISID"]  = np.nan
		missingIS["Lend"]  = np.nan
		missingIS["Rstart"]  = np.nan
	#
	missingIS['Scaffold2'] = np.nan
	results2 = 	missingIS[["ISID","ISstart","ISend","sseqid_first","sstart_first","send_first","sseqid_last","sstart_last","send_last","Scaffold","start","Lend","Scaffold2","Rstart","end","Description","Observations","IS_Match_Type"]]
	columns = ["ISID","ISstart","ISend","Query.ID1","Start1","End1","Query.ID2","Start2","End2", "REF.ID1","REF.Start1","REF.End1","REF.ID2","REF.Start2","REF.End2","Description","Observations","IS_Match_Type"]
	results2.columns = columns
	#************************************************************************
	#***************   Consolidate and annotate   ***************************
	#************************************************************************
	#************************************************************************
	print("Step5: Consolidating results and annotating data...")
	log.write("Step5: Consolidating results and annotating data...")
	log.flush()
	#Now we have to join with the previous results!
	results = concatPdTables(results,results2)
	results = results.sort_values(['Query.ID1','Description','REF.Start1'])
	if args.remove:
		results = results[~results['Query.ID1'].isin(list(identicalScaffolds['qseqid']))]
		results = results[~results['Query.ID2'].isin(list(identicalScaffolds['qseqid']))]
		results = results[~results['REF.ID1'].isin(list(identicalScaffoldsRef['qseqid']))]
		results = results[~results['REF.ID2'].isin(list(identicalScaffoldsRef['qseqid']))]
	results = results.reset_index(drop=True)
	#
	#Search for annotation info
	refFileGB = os.path.join(outputDir,rseq+".gb")
	refFileFasta = os.path.join(outputDir,rseq+".fasta")
	queryFile = os.path.join(outputDir,qseq+".gb")
	queryFileFasta = os.path.join(outputDir,qseq+".fasta")
	# refAnno = genbankToFT(refFileGB)
	# queryAnno = genbankToFT(queryFile)
	queryRecords = recordsDict(queryFile)
	refRecords = recordsDict(refFileGB)
	FinalResults = annotateTable(results)
	#
	#save tables..
	FinalResults[['ISstart', 'ISend', 'Start1', 'End1', 'REF.Start1','REF.End1','Start2', 'End2', 'REF.Start2','REF.End2']] = FinalResults[['ISstart', 'ISend', 'Start1', 'End1', 'REF.Start1','REF.End1','Start2', 'End2', 'REF.Start2','REF.End2']].astype("float")
	FinalResults = FinalResults.sort_values(['Observations','Query.ID1','ISstart'])
	FinalResults = FinalResults.drop_duplicates(["Query.ID1","ISstart","Start1"])
	FinalResults = FinalResults.reset_index(drop=True)
	#replace commas in annotation..
	for idx,row in FinalResults.iterrows():
		FinalResults.loc[idx,'REF.Product_Flank1'] = re.sub(',',';',FinalResults.loc[idx,'REF.Product_Flank1'])
		FinalResults.loc[idx,'QUERY_Flank1.Product'] = re.sub(',',';',FinalResults.loc[idx,'QUERY_Flank1.Product'])
		FinalResults.loc[idx,'REF.Product_Flank2'] = re.sub(',',';',FinalResults.loc[idx,'REF.Product_Flank2'])
		FinalResults.loc[idx,'QUERY_Flank2.Product'] = re.sub(',',';',FinalResults.loc[idx,'QUERY_Flank2.Product'])
	#
	if args.reportSLIS and len(results) > 0:
		SLIS = FinalResults[FinalResults['Observations'] == 'SLIS']
		FinalResults = FinalResults[FinalResults['Observations'] != 'SLIS']
		SLIS.loc[SLIS['Description'] == 'IS only present in query Strain.','Description'] = 'IS location in query Strain.'
		SLIS.loc[SLIS['Description'] == 'IS only present in ref Strain.','Description'] = 'IS location in ref Strain.'
		SLIS = SLIS.reset_index(drop=True)
		FinalResults = FinalResults.reset_index(drop=True)
		SLIS.to_csv(outputDir+'/SLIS.csv')
	#
	FinalResults.to_csv(outputDir+'/FinalResults.csv')
	#
	totalrefIS = totalrefIS.reset_index(drop=True)
	totalrefIS.to_csv(outputDir+'/RefISs.csv')
	totalqueryIS = totalqueryIS.reset_index(drop=True)
	totalqueryIS.to_csv(outputDir+'/QueryISs.csv')
	#New on RC5
	ConsecutiveIS = concatPdTables(ConsecutiveISQuery,ConsecutiveISRef)
	ConsecutiveIS = ConsecutiveIS.reset_index(drop=True)
	ConsecutiveIS.to_csv(outputDir+'/ConsecutiveIS.csv')
	#*****************************************************************************************
	#*****************************************************************************************
	# processedIS,NscaffoldsWithISonBothEnds,scaffoldsWithISonBothEnds,NscaffoldsWithISonOneEnd,scaffoldsWithISonOneEnd
	print("Step6: Printing stats...")
	log.write("Step6: Printing stats...")
	log.flush()
	stats = os.path.join(outputDir,"stats.txt")
	f = open(stats, 'w+')
	#totalqueryIS, totalstatsQuery
	f.write("Query - total Identified ISs: " + str(totalstatsQuery[5]))
	if args.remove:
		f.write("\nQuery - Scaffolds removed from the analysis [identical to reference genome]: " + str(len(identicalScaffolds)))
		f.write("\nQuery - Scaffolds removed from the analysis [identical to reference genome]: " + ', '.join([ str(q) for q in identicalScaffolds['qseqid']]))
	f.write("\nQuery - total Identified ISs in non identical scaffolds: " + str(statsQuery[5]))
	f.write("\nQuery - total processed ISs (ISs with alignment length > 100): " + str(statsQuery[0]))
	f.write("\nQuery - Scaffolds with more than one ISs [Number]: " + str(statsQuery[1]))
	f.write("\nQuery - Scaffolds with more than one ISs [Name]: " + ', '.join([str(q) for q in statsQuery[2]]))
	f.write("\nQuery - Scaffolds with only one IS [Number]: " + str(statsQuery[3]))
	f.write("\nQuery - Scaffolds with only one IS [Name]: " + ', '.join([str(q) for q in statsQuery[4]]))
	if args.remove:
		f.write("\nReference - Scaffolds removed from the analysis [identical to query genome]: " + str(len(identicalScaffoldsRef)))
		f.write("\nReference - Scaffolds removed from the analysis [identical to query genome]: " + ', '.join([str(q) for q in identicalScaffoldsRef['qseqid']]))
	f.write("\nReference - total Identified ISs: " + str(totalstatsRef[5])) #totalrefIS, totalstatsRef
	f.write("\nReference - total processed ISs (ISs with alignment length > 100): " + str(statsRef[0]))
	f.write("\nReference - Scaffolds with more than one ISs [Number]: " + str(statsRef[1]))
	f.write("\nReference - Scaffolds with more than one ISs [Name]: " + ', '.join([str(r) for r in statsRef[2]]))
	f.write("\nReference - Scaffolds with only one IS [Number]: " + str(statsRef[3]))
	f.write("\nReference - Scaffolds with only one IS [Name]: " + ', '.join([str(r) for r in statsRef[4]]))
	f.close()
	#******************************************************************************************
	#******************************************************************************************
	#******************************************************************************************
	# Plot genomic regions using DNA Feature Viewer and FinalResults table. HTML interactive.
	if args.plot:
		print("Step7: Drawing IS insertion sites...")
		log.write("Step7: Drawing IS insertion sites...")
		log.flush()
		#
		plotsPerPage=8
		page = 1
		with PdfPages(outputDir+'/ISGraphicReport.pdf') as pdf:
			for ISs in niter(FinalResults, plotsPerPage):
				figs = plot_funcs(ISs, plotsPerPage,page)
				for fig in figs:
					plt.figure(fig.number)
					pdf.savefig()
				d = pdf.infodict()
				d['Title'] = 'Differential IS insertions graphic report'
				d['Author'] = 'Mauricio J. Lozano | https://github.com/maurijlozano'
				d['Subject'] = 'Graphical comparison of the genomic surroundings of the reported Insertion Sequences.'
				d['Keywords'] = 'Is ISCompare genome'
				d['Creator'] = 'ISCompare'+VERSION
				#
				plt.close('all')
				page += 1
	#Clean if required
	if args.clean:
		trashFiles=[ os.path.join(outputDir, "ref_query_testIS.txt") , os.path.join(outputDir,"ref_queryISSurroundsFromref.fasta"), os.path.join(outputDir,"ref_surroundingBlastRes_R.txt"), os.path.join(outputDir, "ref_surroundingBlastRes_L.txt"), os.path.join(outputDir, "ref_surroundingBlastRes.txt"), os.path.join(outputDir, "ref_surroundings_R.fasta"), os.path.join(outputDir, "ref_surroundings_L.fasta"), os.path.join(outputDir,"ref_surroundings.fasta"), os.path.join(outputDir,"ref_ref_testIS.txt"), os.path.join(outputDir,"refIS.txt"), os.path.join(outputDir,"ref_reduced.fasta"), os.path.join(outputDir,"ref.fasta"), os.path.join(outputDir,"shiftedqueryEnd.fasta"), os.path.join(outputDir,"refGenome2queryGenome.res"), os.path.join(outputDir,"query_ref_testIS.txt"), os.path.join(outputDir,"query_surroundingBlastRes_R.txt"), os.path.join(outputDir,"query_surroundingBlastRes_L.txt"), os.path.join(outputDir,"query_surroundingBlastRes.txt"), os.path.join(outputDir,"query_refISSurroundsFromquery.fasta"), os.path.join(outputDir,"query_query_testIS.txt"), os.path.join(outputDir,"query_surroundings_R.fasta"), os.path.join(outputDir,"query_surroundings_L.fasta"), os.path.join(outputDir,"query_surroundings.fasta"), os.path.join(outputDir,"query_reduced.fasta"), os.path.join(outputDir,"queryIS.txt"), os.path.join(outputDir,"query.fasta"), os.path.join(outputDir,"shiftedrefEnd.fasta"), os.path.join(outputDir,"queryGenome2refGenome.res"), os.path.join(outputDir, "ref.gb"), os.path.join(outputDir, "query.gb")]
		for trash in trashFiles:
			if os.path.exists(trash):
				os.remove(trash)
			else:
				print(trash + 'not found..')
		shutil.rmtree(os.path.join(outputDir,"blastDB"))
	#
	print('Done, the results are on ' + outputDir + '/FinalResults.csv file.\n* Please also check the consecutivIS.csv table for a report\n  of consecutive insertion sequences.\n** QueryIS and RefIS tables contain information of the insertion\n   sequences detected.\n*** Stats.txt contains a report of the IS found.\n**** Running with -p option the ISGraphicReport.pdf file containing \n     a graphic display of the IS genomic location will be \n     generated.')
	log.flush()
	log.close()
#END
