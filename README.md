# Sort-seq analysis of PhoPQ-human AMP interactions

This repository analyzes the results of a high-throughput experiment designed to evaluate peptide sensing by the Salmonella Typhimurium two-component system PhoPQ. We developed a surface-displayed library of 117 human antimicrobial peptides (AMPs) and over 3000 human AMP variants in a laboratory E. coli strain expressing PhoP and PhoQ from S. Typhimurium. We used sort-seq (fluorescence-activated cell sorting followed by next-generation sequencing) to characterize PhoPQ activation by these surface-displayed AMPs.

****

## Packages, programs, and data access

Analysis was performed in Jupyter notebooks using Python 3.7, which was downloaded as part of an Anaconda3 distribution on a computer running Windows 10. Additional code is provided in .py files.

This code was developed with the following package versions installed:
- Python 3.7.4
- pandas 0.25.1
- numpy 1.16.5
- matplotlib 3.1.1
- scipy 1.3.1
- BioPython 1.74: https://biopython.org/
- FlowCal: https://taborlab.github.io/FlowCal/
- regex 2019.11.1: https://pypi.org/project/regex/
- tqdm 4.36.1: https://tqdm.github.io/
- propy3 1.0.0a2: https://pypi.org/project/propy3/

It also requires the following programs:
- FLASH 1.2.11: https://ccb.jhu.edu/software/FLASH/
- usearch 11.0.667: https://www.drive5.com/usearch/download.html

Relevant FACS and flow cytometry data is included in this repository.

Processed data available from GEO under Series accession GSE174191:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174191

Raw NGS data can be downloaded from SRA via the GEO link.

****

## Getting started:

1. Make sure that all packages and supporting programs are installed (see above).
2. Make fixes to the propy3 distribution using the files in the "propy3 fixes" directory (see below).
3. Download raw NGS data from SRA. You may need to change the file names to match those used in the Jupyter notebooks (sort-seq data: sort-seq_library_R1.fastq, sort-seq_library_R2.fastq; pre-sort data: pre-sort_library_R1.fastq, pre-sort_library_R2.fastq).
4. Update the common_dirs_fns.py file to link to the appropriate directories for program, data, and analysis files.


****

## Fixes to propy3 package

Note: I have identified several bugs in propy3 that affect this analysis. For
files to fix these bugs, see the "propy3 fixes" folder included with this
distribution.

1. In PseudoAAC.py, some APAAC features are mislabeled as PAAC features. This is
corrected in propy_functions.py by re-calculating PAAC and APAAC individually and
labeling them manually.

2. The Grantham chemical distance matrix used in propy calculations is incorrect.
The provided matrix in "grantham-chemical-distance-matrix.json" is identical to the
matrix in "schneider-wrede-physicochemical-distance-matrix.json". I have included
an accurate version of the Grantham chemical distance matrix JSON file in this
distribution. This file ("grantham-chemical-distance-matrix.json") needs to be copied
to the "data" subdirectory of the user's propy3 distribution prior to running the
analysis. The Grantham matrix was created using the values from "Grantham.xls"
(located in the propy3 "instruction" subdirectory).

3. The polarity groupings in CTD.py are incorrect, and the mutability value for K is
negative in Autocorrelation.py where it should be positive. The provided files fix
those errors. After installing propy3, you will need to replace the original .py
files in the propy3 directory and then recompile their associated .pyc files:
- Identify your propy3 directory (e.g., C:\ProgramData\Anaconda3\Lib\site-packages\propy)
- Remove the existing copies of Autocorrelation.py and CTD.py from the propy3 directory
- Copy the new copies of Autocorrelation.py and CTD.py to the propy3 directory
- Open a command line interface and navigate to the propy3 directory
- Run the following line to recompile propy3 with the new files:
	python -m py_compile CTD.py Autocorrelation.py

****
