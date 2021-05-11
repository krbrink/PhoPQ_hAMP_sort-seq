from common_dirs_fns import *
import pandas as pd
import numpy as np
from scipy.special import erf
from scipy.optimize import minimize
from collections import Counter
from tqdm.notebook import tqdm
tqdm.pandas()

'''
The lognormal MLE fit workflow is based on the methods from:
Peterman, N., Levine, E. Sort-seq under the hood: implications of design choices on large-scale 
characterization of sequence-function relations. BMC Genomics 17, 206 (2016). 
https://doi.org/10.1186/s12864-016-2533-5
'''


# DEFINE VARIABLES THAT CONTAIN INFORMATION ABOUT FACS BINS
# Import mixed_stats
mixed_stats = pd.read_csv(analysis_path+'mixed_stats.csv',index_col=0)
mixed_stats_bins = mixed_stats.drop('FSC-BSC',axis=0)

# Find bin centers (for guessing mu and sigma in lognormal fits)
mixed_stats_bins['BinGeomMean_Log'] = (np.log(mixed_stats_bins['Minimum'])+np.log(mixed_stats_bins['Maximum']))/2
mixed_stats_bins['BinGeomMean_Linear'] = np.exp(mixed_stats_bins['BinGeomMean_Log'])

# Define upper and lower bin boundaries
upper = mixed_stats_bins['Maximum'].values
lower = mixed_stats_bins['Minimum'].values



# BASIC EQUATIONS REQUIRED FOR LOGNORMAL FITTING AND OPTIMIZATION
# Define functions for optimization of lognormal fit
def lognorm_cdf(x, mu, sigma):
	return 0.5*(1+erf((np.log(x)-mu)/(sigma*np.sqrt(2))))

def mle_func(params, r, lower, upper):
	mu = params[0]
	sigma = params[1]
	return -1.*sum((r)*np.log(lognorm_cdf(upper, mu, sigma)-lognorm_cdf(lower, mu, sigma)))

# Define lognormal PDF function and variance
def lognorm_pdf(x, sigma, mu):	
	return 1/(sigma*x*np.sqrt(2*np.pi))*np.exp(-1*(np.log(x)-mu)**2/(2*sigma**2))

def lognorm_variance(mu, sigma):
	return np.exp(2*mu+sigma**2)*(np.exp(sigma**2)-1)
	
def lognorm_mean(mu, sigma):
	return np.exp(mu+(sigma**2)/2)



# PROCESS COUNT AND FREQUENCY INFORMATION INTO PROBABILITY SPACE
def process_tln_count_df(tln_count_df):
	# This rearranges the columns in tln_count_df such that they are ordered by Bin rather than by Count and Frequency
	tln_count_df = tln_count_df.swaplevel(axis=1).sort_index(axis=1)

	# Weight frequency of each peptide in each bin by the percent of events in that bin on the FACS
	wtd_freq = tln_count_df.loc[:,(mixed_stats_bins.index,'Frequency')].droplevel(level=1,axis=1)*mixed_stats_bins['FracTotalCalculated']
	wtd_freq.columns = pd.MultiIndex.from_product([wtd_freq.columns,['Wtd']])
	tln_count_df = tln_count_df.join(wtd_freq).sort_index(axis=1)
	
	wtd_count_sum = tln_count_df.loc[:,(mixed_stats_bins.index,'Wtd')].sum(axis=1)
	
	# Normalize the weighted frequency so that the weighted frequencies for each peptide add up to 1 across bins
	# This ignores probability loss between and on either side of the fluorescence bins, but since the bins are close together
	# and encompass almost all of the data it should be ok:
	# 95.5% of all events were in the bins based on the following calculation:
	# mixed_stats_bins['FracTotalCalculated'].sum()
	norm_wtd_freq = tln_count_df.loc[:,(mixed_stats_bins.index,'Wtd')].droplevel(level=1,axis=1).divide(wtd_count_sum,axis=0)
	norm_wtd_freq.columns = pd.MultiIndex.from_product([norm_wtd_freq.columns,['NormWtd']])
	tln_count_df = tln_count_df.join(norm_wtd_freq).sort_index(axis=1)
	
	# Convert normalized weighted frequencies into probabilities by dividing by bin width
	# When plotted as a histogram, the sum of the areas of the bins adds to 1, putting us in probability space
	bin_norm = tln_count_df.loc[:,(mixed_stats_bins.index,'NormWtd')].droplevel(level=1,axis=1).divide(mixed_stats_bins['BinWidth'])
	bin_norm.columns = pd.MultiIndex.from_product([bin_norm.columns,['BinNorm']])
	tln_count_df = tln_count_df.join(bin_norm).sort_index(axis=1)
	
	# TotalReads represents the total number of NGS reads in Bins 1 through 8
	tln_count_df['TotalReads'] = tln_count_df.loc[:,(mixed_stats_bins.index,'Count')].sum(axis=1)
	tln_count_df = tln_count_df.reset_index().set_index('Translation')
	return tln_count_df



# PERFORM FITTING
def simple_fit(r, upper, lower):
	'''
	Uses bin centers to estimate mean
	'''
	# mu is the mean of the log-transformed distribution.
	# This is a coarse estimate for the mean based on bin centers in log space (geometric mean).
	simple_mu = sum(r*mixed_stats_bins['BinGeomMean_Log'])

	# sigma is the standard deviation of the log-transformed distribution
	simple_var = (r*(mixed_stats_bins['BinGeomMean_Log']-simple_mu)**2).sum()
	simple_sigma = np.sqrt(simple_var)

	simple_mean = lognorm_mean(simple_mu, simple_sigma)
	return [simple_mu, simple_sigma, simple_mean]

def lognorm_fit(r, upper, lower):
	'''
	Implements MLE for lognormal fitting
	'''
	# Run simple_fit to get initial guesses for mu and sigma
	mu_guess, sigma_guess, mean_guess = simple_fit(r, upper, lower)

	# Minimize the MLE function using mu_guess and sigma_guess as initial parameters
	params0 = [mu_guess, sigma_guess]
	res = minimize(mle_func, params0,
					args=(r.values, lower, upper),
					method='nelder-mead')
	
	# Save the results of the optimization
	mu_est = res.x[0]
	sigma_est = res.x[1]
	mean_est = lognorm_mean(mu_est, sigma_est)
	
	# Save the likelihood calculated from the MLE
	log_likelihood = -1*mle_func([mu_est,sigma_est], r.values, lower, upper)
	likelihood = np.exp(log_likelihood)
	
	# Save whether the optimization was successful
	success = res.success
	
	return [mu_guess, sigma_guess, mean_guess,
			mu_est, sigma_est, mean_est,
			success, log_likelihood, likelihood]


def make_stats_table(tln_count_df_processed, save_name, upper, lower, fit_type):
	# Generate a dataframe with peptides as indices and parameters/results as columns
	# to hold the results of the optimization
	if fit_type == 'Simple':
		stats_cols = ['Mu','Sigma','Mean']
		fit_fn = simple_fit
	elif fit_type == 'LogNorm':
		stats_cols = ['Mu_Simple','Sigma_Simple','Mean_Simple',
					  'Mu_LogNorm','Sigma_LogNorm','Mean_LogNorm',
					  'Success','LogLikelihood','Likelihood']
		fit_fn = lognorm_fit
	else:
		raise ValueError('fit_fn argument is not set to "Simple" or "LogNorm"')

	stats_table = pd.DataFrame(index=tln_count_df_processed.index, columns=stats_cols)
	
	# Reset index such that peptide sequences become a column, for use with progress_apply
	stats_table = stats_table.reset_index()
	
	def apply_fit_function(peptide, fit_fn):
		r = tln_count_df_processed.loc[peptide,(mixed_stats_bins.index,'NormWtd')].droplevel(1)
		return fit_fn(r, upper, lower)
			
	# Fit lognormal distributions to each peptide
	stats_table.loc[:,stats_cols] = stats_table.loc[:].progress_apply(lambda x: apply_fit_function(x['Translation'],fit_fn),
																		  axis=1, result_type='expand').values
		
	# Restore peptide sequences as index
	stats_table = stats_table.set_index('Translation')

	for col in stats_cols:
		if col != 'Success':
		   stats_table[col] = stats_table[col].astype(float)

	# Calculate additional distribution statistics and fold change
	if fit_fn == simple_fit:
		stats_table['Variance'] = lognorm_variance(stats_table['Mu'],stats_table['Sigma'])
		stats_table['StdDev'] = np.sqrt(stats_table['Variance'])
		stats_table['CV'] = stats_table['StdDev']/stats_table['Mean']	
		stats_table['Fold Change'] = stats_table['Mean']/stats_table.loc[stats_table.index=='*','Mean'].values[0]

	elif fit_fn == lognorm_fit:
		for fit_type in ['Simple', 'LogNorm']:
			stats_table['Variance_'+fit_type] = lognorm_variance(stats_table['Mu_'+fit_type],stats_table['Sigma_'+fit_type])
			stats_table['StdDev_'+fit_type] = np.sqrt(stats_table['Variance_'+fit_type])
			stats_table['CV_'+fit_type] = stats_table['StdDev_'+fit_type]/stats_table['Mean_'+fit_type]	
			stats_table['Fold Change_'+fit_type] = stats_table['Mean_'+fit_type]/stats_table.loc[stats_table.index=='*','Mean_'+fit_type].values[0]

	# Save stats_table as a csv file
	stats_table.to_csv(save_name)
	return stats_table


# BOOTSTRAPPING FUNCTIONS
def make_random_bins(tln_count_df, mixed_stats_bins):
	# Create a new dataframe with peptide sequences as an index
	random_bins = pd.DataFrame(index=tln_count_df.index)

	for bin_name in mixed_stats_bins.index:
		# Subset only includes reads in a particular bin
		subset = tln_count_df.loc[:,('Frequency',bin_name)]
		
		# Choose randomly with replacement from tln_count_df based on
		# the frequency of each peptide in that bin in tln_count_df,
		# up to the same number of counts that were originally present
		# in that bin in tln_count_df
		random_choice = np.random.choice(subset.index.values,
						size=tln_count_df.loc[:,('Count',bin_name)].sum(),
						replace=True,
						p=subset.values)
		random_bins[bin_name] = pd.Series(Counter(random_choice))
		
	random_bins.columns = pd.MultiIndex.from_product([['Count'],random_bins.columns])
	random_bins = random_bins.fillna(0)
	random_bins = random_bins.astype(int)

	# Generate frequency data from count data
	random_bins = random_bins.join(
					random_bins[['Count']].div(random_bins[['Count']].sum(axis=0)).rename(columns={'Count':'Frequency'}))
	return random_bins