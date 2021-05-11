import pandas as pd
import numpy as np
import os
import shutil

# Directories
analysis_path = 'Analysis/'
figures_path = 'Figures/'
facs_path = 'FACS/'
ngs_path = 'NGS/'
propy_path = analysis_path+'propy_out/'

# Files
flash_path = 'FLASH-1.2.11-windows-bin/flash'
usearch_path = 'usearch/usearch11.0.667_win32.exe'
peptides_df_path = 'human_AMPs.xlsx'


# *** SHARED FUNCTIONS ***

def rename_unnamed(df):
    '''
    This function handles importing dataframes with MultiIndex columns
    where some columns have no information on one or more column levels.
    Pandas default is to label these unnamed levels with "Unnamed". This
    code replaces "Unnamed" with an empty string.
    Adapted from Stack Overflow: https://stackoverflow.com/a/57224082
	'''
    for i, columns_old in enumerate(df.columns.levels):
        columns_new = np.where(columns_old.str.contains('Unnamed'), '', columns_old)
        df.rename(columns=dict(zip(columns_old, columns_new)), level=i, inplace=True)
    return df

def rewrite_directory(file_path):
    '''
    By default, os.mkdir cannot make a new directory when there is an
    existing directory with the same name. This function removes the old
    directory (and its contents) and replaces it with a new, empty directory.
    '''	
    try:
        os.mkdir(file_path)
    except FileExistsError:
        # Remove trailing slashes if necessary
        if file_path.endswith('/'):
            new_path = file_path[:-1]
        elif file_path.endswith('\\'):
            new_path = file_path[:-2]
        else:
            new_path = file_path
        # Rename the old directory so that it can be removed while the
        # new directory is being created. Otherwise Windows raises an
        # error because it attempts to create the new directory before
        # the old directory has been completely removed.	
        os.rename(file_path, new_path + '_remove')
        shutil.rmtree(new_path + '_remove')
        os.mkdir(file_path)

def to_fasta(df, name_col, seq_col, filename):
    '''
    Generate fasta files based on names and sequences located in two
    dataframe columns
    '''
    with open(filename, 'w') as output:
        for i in df.index:
            output.write('>'+df.loc[i,name_col]+'\n'+df.loc[i,seq_col]+'\n')
			
def capitalize_first(string):
    return string[:1].upper()+string[1:]