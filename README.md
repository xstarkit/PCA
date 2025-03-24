## PCA Code

This package contains three python scripts. `pca_public` is the main PCA script and `pca_functions.py` contains all the functions used by it. `pcatest.py` is a minimal working example of the PCA code, if you are interested in using PCA for other applications than X-ray spectra, I suggest you start there.

### Required Python Modules
 - `numpy`
 - `scipy`
 - `matplotlib`
 - `pylab`
 - `glob`
 - `pyfits`

Most of these are standard modules, generally included with python distributions. `PyFITS` is used for reading data from FITS files, i.e. spectra. The actual decomposition is done using the `numpy.linalg.svd` function.


### Important

You will need to change some of the settings before the code will run. These are at the start of `pca_public.py`. The code assumes a certain directory structure and file naming system, you can change these by modifying the dir_stem and file_stem variables. The spectral files are assumed to be in directories of the form 'Xs', where $X$ is the timestep in seconds.

The `rmf_path` variable will need to be changed - this should point to a response matrix for the detector being used, and is used to convert between spectral channel and energy bin.


### Errors

Errors are calculated by perturbing the input spectra with random noise, up to an amplitude of $N^{0.5}$, where $N$ is the number of counts, then re-running the analysis on the perturbed spectra and looking at the resulting scatter in the output. Set the `n_errors` variable to 0 if you do note wish for errors to be calculated (note that the error calculations significantly increase the run time for large data sets).


### Background Subtraction

The `bkgcorr` variable toggles background subtraction on and off. I recommend testing both, to make sure that spurious components are not being introduced due to background variability or incorrect subtraction. The code currently assumes that background files follow the same naming convention as source files, but with '`src`' replaced with '`bkg`' in the file stem. If you wish to change this convention, you will need to modify line 165 in `pca_functions.py`:
`backfile = specfile.replace('src' ,'bkg')`


### Different Instruments

The code is currently set up to use XMM-Newton spectra. Generalizing it to other X-ray telescopes is not difficult, however. If you are lucky, simply changing the rmf path variable and pointing it at other data should work. If not, you can probably figure it out. If you're still stuck after a week, email me.


### Output

The code produces a set of `.csv` files, containing the average energies, component strengths, and errors, suitable for plotting. The eigenvalues are also written to a `.csv` file.


### Example Script

`pcatest.py` is a simple example of the PCA technique. It takes three functions (a sine wave, sawtooth and stright line) and adds them together in random amounts with additional noise. It then attempts to reconstruct the origial functions from the messy data generated.


### Help and Support

There isn't any. Good luck!


### Distribution and Credit

Feel free to modify, deface, share, distribute, dismantle, and otherwise disperse any and all of the code. If you publish work based on it, I'd appreciate an acknowledgement and some citations!


### Contact

If you have comments or questions, contact me at mparker@sciops.esa.int (expired!)


### Publications
	
* Parker M. L., Fabian A. C., Matt G., Koljonen K. I. I., Kara E., Alston W., Walton D. J., et al., 2015, [MNRAS, 447, 72](https://ui.adsabs.harvard.edu/abs/2015MNRAS.447...72P/abstract). doi:[10.1093/mnras/stu2424](https://doi.org/10.1093/mnras/stu2424)

* Parker M. L., Reeves J. N., Matzeu G. A., Buisson D. J. K., Fabian A. C., 2018, [MNRAS, 474, 108](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474..108P/abstract). doi:[10.1093/mnras/stx2803](https://doi.org/10.1093/mnras/stx2803)
