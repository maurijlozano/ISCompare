# ISCompare: a new tool for the identification of differentially located insertion sequences on bacteria

Bacterial genomes are composed by a core and an accessory genome. The first composed of housekeeping and essential genes, while the second is composed, in its majority, of mobile genetic elements, including transposable elements (Tes). Insertion sequences (ISs) are the smallest TEs formed by imperfect terminal inverted repeats (IRs) and a transposase. ISs are relevant because they have an important role in genome evolution, and contribute to bacterial genome plasticity and adaptability. ISs can spread in a genome, presenting different locations in nearly related strains, and producing phenotypic variations. We developed ISCompare to profile IS mobilization events in related bacterial strains. ISCompare uses blastn to look for ISs on a query genome assembly, extracts the IS flanks and maps them to the reference genome. After filtering and analysis steps, a list of differentially located ISs is reported. ISCompare was validated using artificial genomes with simulated random IS insertions and real sequences from *Escherichia coli*, *Pseudomonas aeruginosa*, *Bordetella pertussis* and *Ensifer meliloti*. In the first case, ISCompare performed very well, achieving high precision (100%) and sensitivity (94%). For real genomes a precision greater than 89% in average was observed, false positive arising mostly from consecutive IS insertions and repeated sequences. Finally we compared ISCompare with ISSeeker, achieving the same or better results, with the advantage that ISCompare can analyse multiple ISs at the same time, and direclty reports list of candidate DLIS. We think that ISCompare provides an easy and straightforward approach to look for differentially located ISs on bacterial genomes.

## Installation and requirements

ISCompare is implemented on python and has the following aditional requirements:

##### Third party python modules:

* [Biopython](https://biopython.org/) (Tested with Biopython v1.76)
* [DNA_features_viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) (Tested with v3.0.3)
* [Pandas](https://pandas.pydata.org/) (Tested with v1.0.3)
* [Numpy](https://numpy.org/) (Tested with v1.18.4)
* [mechanize](https://mechanize.readthedocs.io/en/latest/) (Tested with v0.4.5)

[Blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* For windows users: For blast to run correctly some environment variables need to be created:
 * A modified **path** environment variable to indicate the location of installed blast+ programs
 * A new **BLASTDB** environment variable as pointer to database location, with "e.g. blast_install_dir\db\" as its value
 * A new **BLASTDB_LMDB_MAP_SIZE**, with **1000000** as its value (needed to optimize *makeblastdb* operation when creating new database files)

To install ISCompare simply download ISCompare folder and install all the required third party modules and programs.

### Algorithm overview
![alt text](https://github.com/maurijlozano/ISCompare/blob/master/overview.png "Algorithm overview")

### Running ISCompare

To run ISCompare open a terminal and type: 

```
./ISCompare.py -h
```

```
usage: ISCompare.py [-h] [-q QFILES [QFILES ...]] [-r RFILES [RFILES ...]]
                    [-i ISFILE] [-I] [-Q QACC [QACC ...]] [-R RACC [RACC ...]]
                    [-e YOUREMAIL] [-E EVALUE] [-o OUTPUT] [-x]
                    [-l MINALNLENGTH] [-m MINLENGTH] [-s SURROUNDINGLEN]
                    [-s2 SURROUNDINGLEN2] [-S SHIFT] [-d ISDIFF]
                    [-f SCAFFOLDDIFF] [-p] [-c] [-rs]

ISCompare is a program designed to look for and compare insertion sequence
position between a query and a reference genomes (or WGS Assembly)..

optional arguments:
  -h, --help            show this help message and exit
  -q QFILES [QFILES ...], --QueryFiles QFILES [QFILES ...]
                        Query genome in Genbank format...
  -r RFILES [RFILES ...], --RefFiles RFILES [RFILES ...]
                        Reference genome in Genbank format...
  -i ISFILE, --ISfile ISFILE
                        IS database file [Multifasta nucleotide file]...
  -I, --ISscan          Scan IS on query and reference genomes using ISFinder
                        [Generates IS.fna file in results folder]...
  -Q QACC [QACC ...], --QueryACC QACC [QACC ...]
                        Accession numbers for the query genome to download
                        from NCBI.
  -R RACC [RACC ...], --RefACC RACC [RACC ...]
                        Accession numbers for the reference genome to download
                        from NCBI.
  -e YOUREMAIL, --email YOUREMAIL
                        User email. Required for accession number download
                        mode.
  -E EVALUE, --evalue EVALUE
                        E-value cutoff for blastn search.
  -o OUTPUT, --OutputDir OUTPUT
                        Output folder.
  -x, --IDScaffolds     OMIT the remove identical scaffolds steps.
  -l MINALNLENGTH, --minAlnLength MINALNLENGTH
                        Minimal required length of the alignment of the QIFs
                        to the reference genome.
  -m MINLENGTH, --minLength MINLENGTH
                        Minimum QIFS length to be considered in the analysis.
                        For complete genomes it should be < 2*surroundingLen.
                        For genome assemblies with scaffolds/contigs it should
                        be < surroundingLen.
  -s SURROUNDINGLEN, --surroundingLen SURROUNDINGLEN
                        Nucleotides to extract from upstream and downstream IS
                        blast hits [QIFs].
  -s2 SURROUNDINGLEN2, --surroundingLen2 SURROUNDINGLEN2
                        Nucleotides to extract from QIFs blast hits [RAFs].
  -S SHIFT, --shift SHIFT
                        Shifts the start and end of the QIFs an specified
                        number of nucleotides from the IS.
  -d ISDIFF, --ISdiff ISDIFF
                        ISdiff: minimal difference in nucloetids between qlen
                        and IS alingnment length (discard IS scaffolds).
  -f SCAFFOLDDIFF, --scaffoldDiff SCAFFOLDDIFF
                        scaffoldDiff --> maximal difference in nucloetides for
                        two scaffolds to be considered identical.
  -p, --plot            Plot IS surroundings with genomic features in both
                        query and reference genomes...
  -c, --clean           Clean files...
  -rs, --SLIS           Report IS with the same location (SLIS)...
```

ISCompare requires query and reference genomes in genbank flat format, or its corresponding accession numbers, as input. An optional multifasta DNA file containing all the ISs sequences to be searched in the query and reference genomes can be supplied, otherwise use the -I option which will ISFinderBlast.py script to launch an IS search at ISFinder webpage (http://www-is.biotoul.fr) and download all the found IS sequences.

Optionally, a compilation of IS sequences from ISFinder database can be downloaded from [here](https://github.com/thanhleviet/ISfinder-sequences).
ISCompare uses the following naming convention for ISs: ISname_ISgroup_ISfamily (Required for discrimination between related ISs).

To run with the default parameters type:

```./ISCompare.py -q [query genbank file] -r [reference genbank file] -o [output folder] -c -I -p```



## ISsimulator

ISsimulator.py script inserts a selected number of copies of the insertion sequence at randomly chosen genomic locations, and outputs modified genbank and fasta files, and a table with the location of the inserted ISs. 

```
usage: ISsimulator.py [-h] -i INPUT -s IS [-o OUTPUT] [-dr DIRECTREPEATS]
                      [-n NINSERTIONS]

ISsimulator is a program designed to randomly insert a desired IS into a
provided genome.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Genome in genbank format...
  -s IS, --IS IS        Insertion sequence.
  -o OUTPUT, --Output OUTPUT
                        Output Genbank file name (e.g. name.gb).
  -dr DIRECTREPEATS, --directRepeat DIRECTREPEATS
                        Direct repeat length [Default 0].
  -n NINSERTIONS, --nisertions NINSERTIONS
                        Number of Is insertions [Default 30].
```

To introduce 30 IS at random positions in the query genome query.gb use:

```./ISsimulator.py -i query.gb -s IS.fasta -o artificialGenome.gb -n 30```

IS.fasta must contain only one IS sequence in fasta format.



### Other useful scripts

* extractISbyName.py: This script can be used to extract a single IS from the IS database (multifasta file, IS.fna) by its name.  Type ```./extractISbyName.py -h``` for usage information.
* extractRefFeatures.py: This script can be used to extract all CDS features from a genbank file. Such file is required to run ISSeeker. Type ```./extractRefFeatures.py -h``` for usage information.



