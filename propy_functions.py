from tqdm.notebook import tqdm
import pandas as pd
from propy import PyPro # Python3-compatible version of propy
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import regex
import os
from collections import Counter
tqdm.pandas()

# List of all amino acids
aas = ['A','C','D','E','F',
       'G','H','I','K','L',
       'M','N','P','Q','R',
       'S','T','V','W','Y']
aa_charge = {'D': -1, 'E': -1, 'K': 1, 'R': 1}

# Manually specify autocorrelation and sequence order features
# calculated by propy. Irrelevant features (lag parameter >= sequence
# length) will be removed based on this list.
autocorr = ['MoreauBrotoAuto','MoranAuto','GearyAuto']
autocorr_props = ['AvFlexibility',
                'FreeEnergy',
                'Hydrophobicity',
                'Mutability',
                'Polarizability',
                'ResidueASA',
                'ResidueVol',
                'Steric']

seq_order_features = ['tausw','taugrant']

for ac in autocorr:
    for acp in autocorr_props:
        seq_order_features.append(ac+'_'+acp)

#################################################################
def calculate_propy_features(peptide_sequence, propy_path, save_name=None,
								overwrite=False, lamda=10):
	# Exclude peptides that are smaller than the selected
	# lamda value, which will yield non-sensical resullts
	# for Pseudo Amino Acid Composition (PAAC) parameters.
	if len(peptide_sequence) <= lamda: # default lamda=10
		return None
	
	# Default behavior for save_name is to set it as the
	# peptide sequence
	if save_name is None:
		save_name = peptide_sequence
	
	# Save time by not repeating previously-performed calculations
	if not overwrite:
		if save_name+'.txt' in os.listdir(propy_path):
			return None
	
	# Create a new descriptor object
	DesObject = PyPro.GetProDes(peptide_sequence)
	
	# Get all features (result is a dictionary)
	feats = DesObject.GetALL()
	
	# There's a bug in the propy code where APAAC features over 20 are labeled with PAAC
	# instead of APAAC, thus overwriting some calculated PAAC features.
	# Here we remove APAAC and PAAC features from 'feats' so that they can be
	# re-calculated individually and saved under the correct names.
	for feat in list(feats.keys()):
		if 'PAAC' in feat:
			feats.pop(feat)
	
	# Remove QSO features from features dict
	for feat in list(feats.keys()):
		if 'QSOSW' in feat or 'QSOgrant' in feat:
			feats.pop(feat)

	# Re-calculate APAAC and PAAC and re-add them to 'feats' with the proper labeling
	apaac = DesObject.GetAPAAC()
	for apaac_feat in apaac.keys():
		if apaac_feat[0] != 'A':
			feats['A'+apaac_feat] = apaac[apaac_feat]
		else:
			feats[apaac_feat] = apaac[apaac_feat]
			
	paac = DesObject.GetPAAC()
	for paac_feat in paac.keys():
		feats[paac_feat] = paac[paac_feat]
	
	# Get rid of irrelevant autocorrelation and sequence order features,
	# i.e., those with a lag parameters >= the sequence length.
	for sof in seq_order_features:
		for feat in list(feats.keys()):
			if sof in feat:
				feat_num = int(feat[regex.search(sof, feat).end():])
				if feat_num >= len(peptide_sequence):
					feats.pop(feat)
	
	# Save as tab-delimited file
	with open(propy_path+save_name+'.txt','w+') as out_file:
		for feat in feats.keys():
			out_file.write(feat+'\t'+str(feats[feat])+'\n')


def propy_output_from_fasta(input_fasta, propy_path, overwrite=False, lamda=10):
	'''
	Takes a fasta file as input (input_fasta) and saves propy results for each
	peptide to separate text files in the specified	directory (propy_path).
	'''
	# Open fasta file and read lines, strip white space from lines
	with open(input_fasta,'r+') as in_file:
		lines=in_file.readlines()
	for i, line in enumerate(lines):
		lines[i] = line.strip()
	
	# Generate a dictionary of peptide names and sequences
	pepseq_dict = {}
	for i in range(0,len(lines),2):
		pepseq_dict[lines[i][1:]] = lines[i+1]

	for key in tqdm(pepseq_dict.keys()):
		calculate_propy_features(pepseq_dict[key], propy_path, save_name=key,
									overwrite=overwrite, lamda=lamda)

				
def calculate_charge(peptide):
	total_charge = 0
	for aa in peptide:
		try:
			total_charge += aa_charge[aa]
		except:
			continue
	return total_charge

	
def generate_features_dataframe(peptide_names, propy_path):	
	for pep_name in tqdm(peptide_names):
		propy_file = pep_name+'.txt'
		
		# Generate a dataframe from the propy output file
		pep_df = pd.read_csv(propy_path + propy_file,
							 sep='\t', header=None, names=['Feature', pep_name],
							 na_filter=False)
		pep_df = pep_df.set_index('Feature').transpose()
		
		try:
			# Add new columns if the output file has additional columns that have not yet been included
			# in peptide_features
			peptide_features = peptide_features.reindex(peptide_features.columns.union(pep_df.columns), axis=1)
			peptide_features.loc[pep_name,pep_df.columns] = pep_df.loc[pep_name]
		
		# Create peptide_features if it does not yet exist
		except NameError:
			peptide_features = pep_df
	
	# Calculate relative amino acid compositions
	for aai in aas:
		for aaj in aas:
			if aai == aaj:
				continue
			else:
				peptide_features['pc('+aai+','+aaj+')'] = peptide_features[aai]/(peptide_features[aaj]+peptide_features[aai])
				peptide_features.loc[peptide_features['pc('+aai+','+aaj+')'].isna(),'pc('+aai+','+aaj+')'] = 0

	# Charge
	peptide_features['Charge'] = peptide_features.index.map(calculate_charge)

	# Grand average of hydropathicity
	peptide_features['GRAVY'] = peptide_features.index.map(lambda x: ProteinAnalysis(x).gravy())

	# Isoelectric point
	peptide_features['Isoelectric'] = peptide_features.index.map(lambda x: ProteinAnalysis(x).isoelectric_point())

	# Length
	peptide_features['Length'] = peptide_features.index.map(lambda x: len(x))

	# Molecular weight
	peptide_features['MolecularWeight'] = peptide_features.index.map(lambda x: ProteinAnalysis(x).molecular_weight())

	# Aromaticity
	peptide_features['Aromaticity'] = peptide_features.index.map(lambda x: ProteinAnalysis(x).aromaticity())

	# Instability index
	peptide_features['Instability'] = peptide_features.index.map(lambda x: ProteinAnalysis(x).instability_index())
	
	return peptide_features