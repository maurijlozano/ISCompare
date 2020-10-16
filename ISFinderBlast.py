#!/usr/bin/env python3
'''
ISFinder blast: this scripts uses ISFinder to search for IS on a fasta file.
Written by Mauricio J. Lozano
UNLP - CONICET - Instituto de Biotecnología y Biología Molecular (IBBM)
'''
from Bio import SeqIO
import mechanize
import ssl
import pandas as pd
import lxml.html
import time
import re
import os



def BlastISatISFinder(genome):
	#Submit form using mechanize
	br = mechanize.Browser()
	br.set_handle_robots(False)
	br.addheaders = [('User-agent', 'Firefox')]
	try:
		br.open("https://www-is.biotoul.fr/blast.php")
	except:
		print('SSL: CERTIFICATE_VERIFY_FAILED. Do you want to continue without SSL certificate verification?')
		yesOrNo = input('Please input Y/N: ')
		if yesOrNo.upper() == 'Y':
			ssl._create_default_https_context = ssl._create_unverified_context
			br.open("https://www-is.biotoul.fr/blast.php")
		else:
			print('Please try using a custom IS database fasta file.')
			exit()
	#
	br.select_form(enctype="multipart/form-data")
	filename = genome
	records = SeqIO.parse(filename,'fasta')
	f=open('temp.fasta','+w')
	f.write(">concat_seqs\n")
	for record in records:
		f.write(str(record.seq)+"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
	f.close()
	br.form.add_file(open('temp.fasta'), 'text/plain', filename, name='seqfile')
	expect = br.form.find_control('expect')
	expect.value = '0.0000000001'
	submitToISfinder = br.submit()
	os.remove('temp.fasta')
	resultsPage = submitToISfinder.geturl()
	#read table with pandas
	loading=True
	while loading:
		loading=False
		try:
			br.open(resultsPage)
			submitToISfinder = br.response()
			blastResTab = pd.read_html(submitToISfinder.get_data())[0]
		except:
			if bool(re.search('Please use the form', str(submitToISfinder.get_data()))):
				print('An error with the ISFinder blast query occurred...')
				exit()
			elif bool(re.search('You have to give a sequence!', str(submitToISfinder.get_data()))):
				print('An error with the ISFinder blast query occurred... File too big??')
				exit()
			time.sleep(10)
			loading=True
	#
	blastResTab = blastResTab[blastResTab['E. value'] < 0.0000000001 ]
	br.close()
	return blastResTab

def getISsequencesFromISFinderDB(ISFILE,blastResTab):
	f = open(ISFILE, '+w')
	br = mechanize.Browser()
	br.set_handle_robots(False)
	br.addheaders = [('User-agent', 'Firefox')]
	try:
		br.open("https://www-is.biotoul.fr/scripts/ficheIS.php?name=IS30")
	except:
		ssl._create_default_https_context = ssl._create_unverified_context
	for idx,IS in blastResTab.iterrows():
		ISname = IS['Sequences producing significant alignments']
		ISgroup = IS['Group']
		ISfamily =  IS['IS Family']
		br.open("https://www-is.biotoul.fr/scripts/ficheIS.php?name="+ISname)
		ISpage = br.response().read()
		ISpage = lxml.html.fromstring(ISpage)
		DNA_seq = ISpage.xpath('//div[@class="seq"]')[0].text_content()
		if ISgroup != ISgroup:
			ISgroup = ""
		else:
			ISgroup = "_"+ISgroup
		if ISfamily != ISfamily:
			ISfamily = ""
		else:
			ISfamily = "_"+ISfamily
		f.write(">"+ISname+ISgroup+ISfamily+"\n"+DNA_seq+"\n")
	f.close()
	br.close()
	pass

