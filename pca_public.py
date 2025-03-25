#!/bin/python

# Author: Michael L. Parker

from pca_functions import *

### Settings - see readme file, some of these will need changing.
n_errors=10 # Set to zero to not calculate errors
bkgcorr=True 
custom_energies=[] # Example: [0.5,1.,2.,5.,10.]
dir_stem='products_pca*' # Assumed to contain directories of form '10000s', where 10000 is the timestep.
file_stem='src*.fits'
output_stem='pca_output'
rmf_path = '/scratch/mlparker/nustar/agn/new_mcg6/xmm/PN.rmf'
emin=0.4
emax=10.0
###


print_header()

spectrum_errors=None
eigenval_errors=None

# Check custom_energies and read in pca user parameters
if len(custom_energies)> 2:
	energies=custom_energies
	nbins=len(energies)-1
	tstep, n_spectra = get_pca_params(False)
elif len(custom_energies)==0:
	tstep, n_bins, n_spectra = get_pca_params()
	energies=calc_bins(emin,emax,n_bins)
else:
	print 'Something is very wrong with custom_energies - it should have either >2 or 0 values'

# Find relevant directories
dir_list=get_dirs(tstep,dir_stem)

# Find all files within those directories
spec_files = get_files(dir_list,file_stem)

# Get channels corresponding to energy bins from response matrix
channel_bins=read_rmf(rmf_path, energies)

# Read in all spectra, calculate perturbed spectra for errors
spectrum_array, perturbed_spectra =read_data(spec_files,n_errors,channel_bins,bkgcorr,tstep)

# Normalize by subtracting and dividing by the mean spectrum
# This removes the instrumental response and restricts the analysis to variable components
normalized_spectrum_array, mean_spectrum=normalize_spectra(spectrum_array)

# Calculate principal components
pc_array,eigenvals=decompose(normalized_spectrum_array)

# Calculate errors
if n_errors>0:
	spectrum_errors, eigenval_errors = calc_errors(n_errors, perturbed_spectra, mean_spectrum, eigenvals)

# Print out eigenvalues (and corresponding errors)
print_eigenvalues(eigenvals,n_spectra,eigenval_errors)

# Plot spectra
plot_results(n_spectra, pc_array, energies, eigenvals, spectrum_errors)

# Write PCA spectra and eigenvalues to files
write_spectra(output_stem,n_spectra,pc_array,energies,spectrum_errors)
write_eigenvalues(output_stem,eigenvals,eigenval_errors)


