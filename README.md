UniFam
======
**Full sequence protein annotation with UniProt-based Families**

UniFam includes a database of protein families generated based on UniProt, and provides a pipeline to annotate either a genome or a set of proteins. 
UniFam also provides users an option to run other bioinformatics tools to extract genomic information and reconstruct metabolic pathways.
Go to http://unifam.omicsbio.org to submit your annotation jobs with Unifam!
See our [paper](http://www.biomedcentral.com/1471-2148/14/207) in BMC Evolutionary Biology for a phylogenomic study of prokaryotes based on UniFam annotation, and comparison benchmark with other protein annotation programs. 

Envirionment Requirement
------------------------
**Required**:

	1. Python 2.7.2 or above, and "argparse" module to run the Python scripts
	2. HMMER 3.1 for protein annotation
	3. Pipeline was tested on Linux and MacOS systems, not Windows.

**Optional**:

	1. Prodigal v3.0 for prokaryotic gene calling
	2. RNAmmer for predicting RNAs
	3. tRNAscan for predicting tRNAs
	4. Pathway-Tools for pathway reconstruction

Package contents
------------------------
	1. **data**: UniFam databases, including the database in a whole, and two sub-databases, one for prokaryotic proteins and one for eukaryotic proteins specifically. 
	The annotation files for the three databases and their aliases are also included. [Download](http://unifam.omicsbio.org/downloads).
	2. **example**: a test example with input and output
	3. **src**: python scripts for the pipeline

Instructions
------------------------
	1. Download the package, and extract the files from the archive.

	2. Provide the configuration file.

	An easier way to do this is to copy the sample configuration file and edit it according to the provided direction. 
	Start from section "UniFam" at the bottom. Please make sure to change the names of directories, check the input format, and decide if any of the modules should be executed.
	If a particular module is needed, go to the corresponding section in the configuration file and configure the path to its executable, and parameters.

	3. Run the main script in the src directory.
	```
	#!bash
	python UniFam.py -c configFile -i inputfile
	```
	configFile is as described in step 2;
	inputfile should be in fasta format, either of genome sequence (DNA), or proteins in a genome;
	outputfile is the output annotation file for the provided genome (proteins)

Update 
------------------------
	v1.1 Added script to generate GenBank file for NCBI submission for Eukaryotic organisms. 
	     Contigs, protein sequences, genmark .gtf, and Unifam .annot files are required for this function.
