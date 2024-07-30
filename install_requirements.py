#!/usr/bin/env python3
import os
def main():
	if os.popen('blastn -h').read():
		print('blastn already installed...')	
	else:
		print('Installing blast: sudo password will be required')
		os.system('sudo apt-get install ncbi-blast+ -y')
		

	print('Checking pakages...')
	import subprocess, sys


	try:
		import Bio
		print(f"Bio already installed...")
	except:
		print(f"Installing Bio...")
		subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'biopython==1.84'], stdout=subprocess.DEVNULL)


	try:
		import dna_features_viewer
		print(f"dna_features_viewer already installed...")
	except:
		print(f"Installing dna_features_viewer...")
		subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'dna-features-viewer==3.1.3'], stdout=subprocess.DEVNULL)

	try:
		import pandas
		print(f"pandas already installed...")
	except:
		print(f"Installing pandas...")
		subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandas==2.2.2'], stdout=subprocess.DEVNULL)

	try:
		import numpy
		print(f"numpy already installed...")
	except:
		print(f"Installing numpy...")
		subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy==2.0.1'], stdout=subprocess.DEVNULL)

	try:
		import mechanize
		print(f"mechanize already installed...")
	except:
		print(f"Installing mechanize...")
		subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'mechanize==0.4.10'], stdout=subprocess.DEVNULL)

	try:
		import lxml
		print(f"lxml already installed...")
	except:
		print(f"Installing lxml...")
		subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'lxml==5.2.2'], stdout=subprocess.DEVNULL)

	try:
		import beautifulsoup4
		print(f"beautifulsoup4 already installed...")
	except:
		print(f"Installing beautifulsoup4...")
		subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'beautifulsoup4==4.12.3'], stdout=subprocess.DEVNULL)

	print('All the requirements were installed...')

if __name__ == "__main__":
	main()