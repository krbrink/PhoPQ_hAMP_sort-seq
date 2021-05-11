from common_dirs_fns import *
import pandas as pd
import numpy as np
import regex
import subprocess
from collections import Counter
from Bio.Seq import Seq
from tqdm.notebook import tqdm
tqdm.pandas()


# COMBINE FORWARD AND REVERSE READS AND COUNT PAIRED-END READS
def run_flash(flash_in1, flash_in2, flash_out):
	'''Combines forward and reverse reads to generate a paired-end FASTQ file using FLASH'''
	# Total length of each amplicon was designed to be 425 nt. Performed 2 x 250 nt sequencing.
	# The overlap between each pair of reads should therefore ideally be 75 nt.
	# Use 10 nt as minimum overlap and 150 nt as maximum overlap settings to include sequences
	# that deviate from this ideal calculation.
	
	options = {'max-overlap': 150,
				'min-overlap': 10,
				'output-directory': flash_out,
				#changes number of mismatches tolerated in overlap region, default=0.25
				'max-mismatch-density': 0.1 
				}
	
	# Generate command
	options_str = ''
	for key in options.keys():
		options_str = options_str + ' --' + key + '=' + str(options[key])

	to_flash = ' '.join(['"' + flash_path + '"',
						options_str,
						'"' + flash_in1 + '"',
						'"' + flash_in2 + '"'])
	print(to_flash)
	
	p = subprocess.Popen(to_flash,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT,
                     shell=True)

	# Print progress output in real time
	for line in iter(p.stdout.readline, b''):
		print(">>> " + line.decode().rstrip())
		

def get_seqs_only(file_path):
	'''Saves just the reads from the fastq file to a new file'''
	# Open fastq file
	with open(file_path,'r') as fastq_file:
		# Open output file
		with open(file_path+'_sequences_only.txt','w+') as out_file:
			for line_counter,line in enumerate(fastq_file):
				# Ignore lines that don't contain sequence data
				if line_counter%4 == 1:
					out_file.write(line.strip() + '\n')
                    
				#Keep track of progress every 100,000 reads
				if (line_counter/4)%100000 == 0:
					print(int(line_counter/4))
					
					
def count_paired_end_reads(flash_in1, flash_in2, flash_out):
	run_flash(flash_in1, flash_in2, flash_out)
	
	# Generate a text file containing only sequence information from FLASH
	# output FASTQ file
	print('Writing sequences file')
	get_seqs_only(flash_out + '/out.extendedFrags.fastq')
	
	# Count the number of times each unique DNA sequence appears in the file
	# Adapted from Stack Overflow: https://stackoverflow.com/a/14260441
	# This takes up to several minutes to run and may not be compatible with
	# very large files
	print('Counting sequences')
	with open(flash_out + '/out.extendedFrags.fastq_sequences_only.txt') as infile:
		counts = Counter(l.strip() for l in infile)
		
	# Reduce the size of the dictionary by combining counts for forward and reverse reads that are
	# reverse complements of one another

	# Change from Counter object to dictionary
	counts_dict = dict(counts)

	# Generate a separate variable containing all of the keys to avoid
	# changing list while iterating
	static_keys = list(counts_dict.keys())

	for i,seq in enumerate(static_keys):
		if seq not in counts_dict.keys(): # sequence has already been removed
			continue
		rc = str(Seq(seq).reverse_complement())
		if rc in counts_dict.keys(): # sequence and its reverse complement both in counts_dict
			counts_dict[seq] += counts_dict.pop(rc) #remove reverse complement and
													#add counts to original sequence
	
	# Generate dataframe of paired-end read counts
	merged_reads_df = pd.DataFrame.from_dict(data=counts, columns=['Count'], orient='index')
	merged_reads_df.index.name = 'Merged'
	merged_reads_df = merged_reads_df.reset_index()
	
	return merged_reads_df




# ASSIGN READS TO BINS BASED ON BARCODE SEQUENCES
# NOTE: This is currently set up to process reads where fwd and rev barcodes are identical
# This code may require small modifications to handle different fwd/reverse barcodes

def assign_barcode(seq, barcode_lookup_table, barcode_regexes):
	number_mismatches = None
	best_matches = []
	best_score = len(barcode_lookup_table.index[0]) # Initially set best mismatch score to barcode length

	for barcode in barcode_lookup_table.index:
		match = regex.search(barcode_regexes[barcode], seq)
        
		if match == None: # No match for barcode in sequence
			continue
		elif sum(match.fuzzy_counts)==0: # Perfect match
			best_matches = [barcode]
			number_mismatches = 0
			break
		elif sum(match.fuzzy_counts) < best_score: # Better match score than previous
			best_matches = [barcode]
			best_score = sum(match.fuzzy_counts)
			number_mismatches = best_score
		elif sum(match.fuzzy_counts) == best_score: # Same match score as previous best
			best_matches.append(barcode)
            
	if len(best_matches) == 1:
		# return assigned_bin, number_mismatches
		return (barcode_lookup_table.loc[best_matches[0],'Bin'], number_mismatches)
	elif (len(best_matches) > 1) or (len(best_matches) == 0):
		return (None, number_mismatches)
    
def consensus_bin(binf, binr):
	if binf == binr:
		return binf
	else:
		# if one bin assignment is ambiguous, return non-ambiguous assignment
		if (binf == None)^(binr == None): # ^ is xor function
			return binf if binf != None else binr
	return None





# IDENTIFY PEPTIDE-ENCODING DNA SEQUENCE
def find_after_rev_primer(input_seq, rev_primer_regex):
	match = regex.search(rev_primer_regex, input_seq)
	if match != None:
		return input_seq[match.end():]
	else:
		return None

def consensus_after_rev_primer(after_fwd, after_rev):
	# The primer binding site should only be present in either the sequence
	# or its reverse complement and not both (^ is xor operator)
	if (after_rev == None)^(after_fwd == None):
		after = after_rev if after_rev is not None else after_fwd
	else:
		after = None
	return after
	




# TRANSLATE PEPTIDE SEQUENCE
'''
There are two translation functions here to distinguish between two different cases:
raw_translation will translate a DNA sequence in its entirety regardless of the presence of a stop codon
(e.g., it will return output that looks like SDFGH*AKLM).
actual_translation takes the output from raw_translation and crops the sequence at the stop codon *.
Separating these two functions enables us to identify any sequences that do not contain stop codons,
which may not allow functional peptide expression and should be excluded.
'''

def raw_translation(input_seq):
	# Note: This function will raise a warning for any sequences whose lengths are not divisible by 3.
	#       Ignoring this warning is acceptable because ultimately we are interested in the translated
	#       sequence that occurs before the stop codon (see actual_translation) and the extra bases are
	#		irrelevant.
	if input_seq != None:
		return str(Seq(input_seq).translate())
	else:
		return None
    
def actual_translation(input_tln):
	# Input of None returns None
	if type(input_tln) != str:
		return None

	# Sequences lacking a stop codon are encoded as an empty string
	if '*' not in input_tln:
		return ''

	tln = input_tln.split('*')[0]
	if tln == '':
		return '*' # Sequences that start with a stop codon are encoded as '*'
	else: 
		# Sequences that do not start with a stop codon return the
		# amino acid translation without '*'
		return tln