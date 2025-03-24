#!/bin/python
from numpy.linalg import svd
import numpy as n
# import nimfa
# import mdp
import glob
import pyfits
import os
import sys
import pylab as p
from matplotlib.ticker import *
import matplotlib.animation as animation
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from time import sleep

def area_correct(spectrum,arf_path,energy_bins):
	print '\nFluxing spectrum...'
	arf=pyfits.open(arf_path)['SPECRESP'].data
	energies=[]
	areas=[]

	# read response stuff
	for row in arf:
		energies.append((row[0]+row[1])/2.)
		areas.append(row[2])

	# rebin response
	area_spline=UnivariateSpline(energies,areas,s=0,k=3)

	binned_areas=[]
	for elow,ehigh in zip(energy_bins[:-1],energy_bins[1:]):
		binned_areas.append(quad(area_spline,elow,ehigh)[0]/(ehigh-elow))

	fluxed_spectrum=[]
	for i,flux in enumerate(spectrum):
		fluxed_spectrum.append(flux/binned_areas[i])

	return fluxed_spectrum

def animate_components(unfolded_mean,pc_array,n_spectra,energies,savedir=None):
	print '\nAnimating...'
	amplitude=4
	colours=['k','r','b','g','m','c','y']
	p.ion()
	p.figure()#(figsize=(15,10))
	if savedir != None:
		if not os.path.exists(savedir):
			os.mkdir(savedir)
		if savedir[-1] != '/':
			savedir+='/'
	for i in range(0,n_spectra):
		p.clf()
		ax=p.subplot(111)
		ax.set_yscale('log')
		ax.set_xscale('log')
		ax.set_xlabel('Energy (keV)')
		ax.set_ylabel("Flux (I've lost all track of units)")
		ax.set_xlim(min(energies),max(energies))
		ax.xaxis.set_major_locator(FixedLocator([0.1,0.2,0.5,1,2,5,10,20,50,100]))
		ax.xaxis.set_minor_locator(FixedLocator([0.3,0.4,0.6,0.7,0.8,0.9,3,4,6,7,8,9,15,30,40,60,70,80,90]))
		ax.xaxis.set_major_formatter(ScalarFormatter())
		mean,=ax.step(energies[:-1],unfolded_mean,color='gray',where='mid')
		line,=ax.step(energies[:-1],unfolded_mean,color=colours[i],linewidth=2,where='mid')
		p.legend(['Mean spectrum','Principal component '+str(i+1)])
		j=0
		for phase in n.linspace(0,2*n.pi,200):
			j+=1
			line.set_ydata([(1+x*amplitude/(i+1)*n.sin(phase))*y for x,y in zip(pc_array[i],unfolded_mean)])
			p.draw()
			if savedir != None:
				p.savefig(savedir+'fig_comp%s_%s.png' % (i+1,j))
	return 0.

def read_int(prompt,default):
	'''Prompts the user for an integer value, or uses default'''
	userval=raw_input(prompt)
	if len(userval)>0:
		try:
			return int(userval)
		except:
			print 'Not a valid integer!'
			exit()
	else:
		return default


def get_pca_params(t_step_default=10000,readbins=True):
	'''Gets the relevant user parameters for each run'''
	tstep=read_int('\tEnter timestep (%s): ' % str(t_step_default),t_step_default)
	n_spectra=read_int('\tEnter number of component spectra to plot (3): ', 3)
	if readbins:
		n_bins=read_int('\tEnter desired number of bins (50): ',50)
		return tstep, n_bins, n_spectra
	else:
		return tstep, n_spectra


def get_dirs(tstep,dir_stem):
	'''Find directories'''
	currentdir=str(os.getcwd())+'/'
	return [currentdir+f+'/%ss/' % str(tstep) for f in glob.glob(dir_stem)]


def print_header(telescope='xmm',fctr_type='pca'):
	'''Prints some header text'''
	if fctr_type=='pca':
		fctr_str='Principle Component Analysis'
	elif fctr_type=='nmf':
		fctr_str='Non-Negative Matrix Factorization'
	elif fctr_type=='ica':
		fctr_str='Independent Component Analysis'
	if telescope=='xmm':
		print '\nXMM %s, M. Parker 26/10/14\n' % fctr_str
	elif telescope=='suzaku_xis':
		print '\nSuzaku XIS %s, M. Parker 05/11/14\n' % fctr_str
	elif telescope=='multi':
		print '\n%s, M. Parker 05/11/14\n' % fctr_str
	else:
		print '\n%s, M. Parker 05/11/14\n' % fctr_str


def get_files_starts(folderlist, file_stem):
	'''Find all files from within folderlist, matching file_stem'''
	filenumbers=[]
	specstarts=[]
	specfiles=[]

	for folder in folderlist:
		count=0
		#Locate and sort all files in given folder
		for specfile in glob.glob(folder+'*'+file_stem):
			count+=1
			spec=pyfits.open(specfile)
			specdata=spec['SPECTRUM'].data
			specheader=spec['SPECTRUM'].header
			gtidata=spec[3].data
			if 'EXPOSURE' not in specheader:
				print 'Warning! Keyword EXPOSURE not in file %s header!'
				print 'Skipping file.'
				count -=1
			else:
				specstarts.append(gtidata[0][0])
				specfiles.append(specfile)
		filenumbers.append(count)

		if count == 0:
			print 'Warning! No files in folder %s\nIs filename correct?' % folder

	specfiles=[s for t,s in sorted(zip(specstarts,specfiles))]
	specstarts=sorted(specstarts)

	nit_time=specstarts[0]
	specstarts=[t-nit_time for t in specstarts]

	print '\nTotal intervals: ',sum(filenumbers)
	return specfiles, specstarts

def get_files(folderlist, file_stem):
	'''Find all files from within folderlist, matching file_stem'''
	filenumbers=[]
	specstarts=[]
	specfiles=[]

	for folder in folderlist:
		count=0
		#Locate and sort all files in given folder
		for specfile in glob.glob(folder+'*'+file_stem):
			count+=1
			spec=pyfits.open(specfile)
			specdata=spec['SPECTRUM'].data
			specheader=spec['SPECTRUM'].header
			if 'EXPOSURE' not in specheader:
				print 'Warning! Keyword EXPOSURE not in file %s header!'
				print 'Skipping file.'
				count -=1
			else:
				specfiles.append(specfile)
		filenumbers.append(count)

		if count == 0:
			print 'Warning! No files in folder %s\nIs filename correct?' % folder

	print '\nTotal intervals: ',sum(filenumbers)
	return specfiles


def calc_bins(emin,emax,n_bins):
	'''calculate logarithmic energy bins'''
	energies = [10.**i for i in n.linspace(n.log10(emin),n.log10(emax),n_bins+1)]
	return energies


def read_rmf(rmf_path, energies):
	'''Finds channels corresponding to energy bins using response matrix'''
	rmf_file=pyfits.open(rmf_path)
	ebounds = rmf_file['EBOUNDS'].data
	emins=[]
	emaxs=[]
	channels=[]
	# Find channels corresponding to energy bins
	for row in ebounds:
		channels.append(row[0])
		emins.append(row[1])
		emaxs.append(row[2])
	channel_bins=[]
	for e in energies:
		for i in range(1,len(channels)-1):
			if emins[i]<=e<emins[i+1]:
				channel_bins.append(channels[i])

	return channel_bins


def toolbar_update(fraction,toolbar_width):
	sys.stdout.write('[')
	sys.stdout.write('%s' % ('=' * max([fraction-1,0])))
	sys.stdout.write('>')
	sys.stdout.write('%s' % (' ' * (toolbar_width-max([fraction,1]))))
	sys.stdout.write(']')
	sys.stdout.write("\b" * (toolbar_width+2))
	sys.stdout.flush()


def read_spec(specfile,first,bkgcorr):
	'''Read spectrum from file'''
	channels=[]
	spec=pyfits.open(specfile)
	specdata=spec['SPECTRUM'].data
	specheader=spec['SPECTRUM'].header

	if bkgcorr:
		backfile = specfile.replace('src' ,'bkg')
		back=pyfits.open(backfile)
		backdata=back['SPECTRUM'].data
		backheader=back['SPECTRUM'].header
		src_backscal=specheader['BACKSCAL']
		backscal=backheader['BACKSCAL']
		back_exptime=backheader['EXPOSURE']

	exptime=specheader['EXPOSURE']
	count_rates=[]
	counts_list=[]

	for i in range(0,len(specdata)):
		row=specdata[i]

		if first:
			channels.append(row[0])
		if bkgcorr:
			backrow=backdata[i]
			counts=row[1]-backrow[1]*(exptime/back_exptime)*(src_backscal/backscal)
		else:
			counts=row[1]
		counts_list.append(counts)
		count_rates.append(counts/exptime)

	if first:
		return count_rates,counts_list,exptime,channels
	else:
		return count_rates,counts_list,exptime


def perturb(spectrum):
	'''Add/subtract a value between N^1/2 and 0 from each bin in a spectrum'''
	newspec=[max([c,0])+n.random.randn()*max([c,0])**0.5 for c in spectrum]
	return newspec


def rebin(chans,countlist,bins):
	'''Rebinning function for spectra'''
	newcounts=[]
	for bmin,bmax in zip(bins[:-1],bins[1:]):
		newcounts.append(sum(countlist[bmin:bmax])/(bmax-bmin))
	return newcounts


def read_data(specfiles, n_errors, channel_bins, bkgcorr, tstep):
	'''Read spectral data from all files in specfiles, correct for background and calculate errors'''

	toolbar_width=60
	rand_spectra=[]
	datalist=[]
	first=True
	exposures=[]

	foldernum=0
	totalspectra=0
	nfiles=len(specfiles)
	filecount=0
	exposure=0.

	if n_errors>0:
		print '\nReading data and calculating perturbed spectra...\n'
		errorcalc=True
	else:
		print '\nReading data...\n'
		errorcalc=False

	totalcounts=0
	for specfile in specfiles:
		fraction=int(toolbar_width*filecount/nfiles)
		toolbar_update(fraction,toolbar_width)

		if first:
			#If first, get channels as well as spectrum
			counts,fluxes,exptime,channels=read_spec(specfile,first,bkgcorr)
			first=False
		else:
			counts,fluxes,exptime=read_spec(specfile,first,bkgcorr)

		totalcounts+=sum(counts)*exptime
		counts=rebin(channels,counts,channel_bins)

		if errorcalc:
			# Find perturbed spectra for errors
			rand_spec=[]
			for i in range(0,n_errors):
				spec_i=perturb(fluxes)
				spec_i=[i/exptime for i in spec_i]
				spec_i=rebin(channels,spec_i,channel_bins)
				rand_spec.append(spec_i)

		total=sum(counts)

		if total>0 and exptime>tstep/10.:
			totalspectra+=1
			datalist.append(counts)
			exposure+=(exptime)
			if errorcalc:
				rand_spectra.append(rand_spec)
		else:
			sys.stdout.write('%s' % (' '*(toolbar_width+2)))
			sys.stdout.write('%s' % ('\b'*(toolbar_width+2)))
			sys.stdout.flush()
			if total <=0:
				sys.stderr.write('Warning: No counts in time bin')
				sys.stderr.write('\n')
			elif exptime <=tstep/10.:
				sys.stderr.write('Warning: Exposure time too short')
				sys.stderr.write('\n')
			sys.stderr.write(specfile+ '\n')
			sys.stderr.write('\n')
		filecount+=1

	sys.stdout.write('[%s' % ('=' * toolbar_width))
	sys.stdout.write('\n')
	print 'Exposure time:', exposure
	print 'Total counts:',int(totalcounts),'\n'

	spectrum_array=n.array(datalist)

	return spectrum_array, rand_spectra


def normalize_spectra(spectrum_array,subtract=True):
	means=n.mean(spectrum_array,axis=0)
	print '\nNormalizing spectra...'
	newlist=[]
	for row in spectrum_array:
		if subtract:
			newrow=[(c-m)/m for c,m in zip(row,means)]
		else:
			newrow=[c/m for c,m in zip(row,means)]
		newlist.append(newrow)
	return n.transpose(n.array(newlist)), means


def calc_errors(n_errors,rand_spectra,means,eigenvals):
	print '\nCalculating errors...'
	pca_list=[]
	new_rand_spectra=[]
	eigenvals_ptbd=[]
	for rand_spec in rand_spectra:
		newspec=[]
		for spectrum in rand_spec:
			newspec.append([(c-m)/m for c,m in zip(spectrum,means)])
		new_rand_spectra.append(newspec)
	for i in range(0,n_errors):
		datalist=[]
		for spec_set in new_rand_spectra:
			datalist.append(spec_set[i])
		error_array=n.transpose(datalist)
		U,A,V=svd(error_array)
		U=n.transpose(U)
		if i==0:
			U0=n.copy(U)
		tempvals=[]
		for k in A:
			tempvals.append(k**2)
		total=sum(tempvals)
		eigenvals_ptbd.append([k/total for k in tempvals])
		# eigenvals_ptbd.append(tempvals)
		for j in range(0,len(U)):
			row=U[j]
			if sum([a*b for a,b in zip(U0[j],row)])<0:
				row=[k*-1 for k in row]
				U[j]=row
		pca_list.append(U)
	errors=[]
	eigenvals_ptbd=n.array(eigenvals_ptbd)
	eigenerrs=n.std(eigenvals_ptbd,axis=0)
	for i in range(0,len(pca_list[0])):
		errors.append([])
		mean_flux=means[i]

		#iterate over energies
		for e in range(0,len(pca_list[0][0])):
			vals=[]
			#iterate over spectra
			for j in range(0,n_errors):
				val=pca_list[j][i][e]
				vals.append(val)
			meanval=sum(vals)/len(vals)
			variance=0.
			for val in vals:
				variance+=(val-meanval)**2/float(n_errors-1)
			error=variance**0.5
			errors[i].append(error)
	return errors, eigenerrs


def decompose(spectrum_array):
	print '\nDecomposing spectra using SVD...'
	U,A,V=svd(spectrum_array)
	U=n.transpose(U)
	eigenvals=[i**2 for i in A]
	total=sum(eigenvals)
	eigenvals=[i/total for i in eigenvals]
	return U,eigenvals

def decompose_nmf(spectrum_array,n_spectra):
	print '\nDecomposing spectra using NMF...'
	fctr = nimfa.mf(spectrum_array, method="nmf", max_iter=10000, rank=n_spectra, update='divergence', objective='div')
	fctr_res = nimfa.mf_run(fctr)
	a=n.transpose(n.array(fctr_res.basis()))
	coeffs=n.array(fctr_res.coef())
	return a,coeffs

def decompose_ica(spectrum_array):
	print '\nDecomposing spectra using ICA...'
	ics=n.transpose(mdp.fastica(spectrum_array,whitened=True))
	return ics

def print_eigenvalues(eigenvals,n_spectra,eigenval_errors=[]):

	'''Prints out the fractional variability of components up to n_spectra'''
	print '\nPercentage variability in 1st %s components:\n' % n_spectra
	if len(eigenval_errors)>0:
		for i in range(0,n_spectra):
			print '\tEigenvector %s:' % str(i+1),str(eigenvals[i]*100)[0:6],\
			'+/-', str(eigenval_errors[i]*100)[0:6], '%'
		print '\n\tRemaining variability:',str(sum(eigenvals[n_spectra:])*100)[0:6],'%'

	else:
		for i in range(0,n_spectra):
			print '\tEigenvector %s:' % str(i+1),str(eigenvals[i]*100)[0:6], '%'
		print '\nRemaining variability:',str(sum(eigenvals[n_spectra:])*100)[0:6],'%'

def plot_spectra(spectrum_array,energies):
	ax=p.subplot(111)
	for spectrum in spectrum_array:
		p.plot(energies[:-1],spectrum)
	ax.set_xscale('log')
	ax.set_yscale('log')
	p.show()

def plot_results(n_spectra,pc_array, energies, eigenvalues, spectrum_errors=None, eigenerrors=None):
	'''Plot the component spectra and eigenvalues'''
	print '\nPlotting spectra...'
	colours=['k','r','b','g','m','c','y']
	axisfontsize=14

	spectrumfigure=p.figure(figsize=(14,max([6,2*n_spectra])))
	for specnum in range(0,n_spectra):
		# subplotnum=100*n_spectra+21+specnum
		sub_plot=p.subplot2grid(shape=(n_spectra,2),loc=(specnum,0))
		if spectrum_errors != None:
			if n_spectra<=len(colours):
				p.errorbar(energies[:-1],pc_array[specnum],spectrum_errors[specnum],ls='none',marker='x',color=colours[specnum])
			else:
				p.errorbar(energies[:-1],pc_array[specnum],spectrum_errors[specnum],ls='none',marker='x')
		else:
			if n_spectra<=len(colours):
				p.plot(energies[:-1],pc_array[specnum],marker='x',color=colours[specnum])
			else:
				p.plot(energies[:-1],pc_array[specnum],marker='x')
		sub_plot.set_xscale('log')
		sub_plot.xaxis.set_major_locator(FixedLocator([0.1,0.2,0.5,1,2,5,10,20,50,100]))
		sub_plot.xaxis.set_minor_locator(FixedLocator([0.3,0.4,0.6,0.7,0.8,0.9,3,4,6,7,8,9,15,30,40,60,70,80,90]))
		sub_plot.xaxis.set_major_formatter(ScalarFormatter())
		p.xlim(min(energies),max(energies))
	p.xlabel('Energy (KeV)',fontsize=axisfontsize)
	p.ylabel('Normalised Count Rate',fontsize=axisfontsize)

	sub_plot=p.subplot2grid(shape=(n_spectra,2),loc=(0,1),rowspan=3)
	p.xlabel('Eigenvector',fontsize=axisfontsize)
	p.ylabel('Fractional Variability',fontsize=axisfontsize)
	eigenvalmarker='o'
	for i in range(0,n_spectra):
		if n_spectra<=len(colours):
			p.plot(i+1,[eigenvalues[i]],marker=eigenvalmarker,color=colours[i],ls='-',ms=8)
		else:
			p.plot(i+1,[eigenvalues[i]],marker=eigenvalmarker,ls='-',ms=8)
	p.plot(range(n_spectra+1,len(eigenvalues)),eigenvalues[n_spectra:-1],marker=eigenvalmarker,ls='None',color='y',ms=5)
	sub_plot.set_yscale('log')

	p.show()


def plot_results_nmf(n_spectra,pc_array, energies):
	'''Plot the component spectra and eigenvalues'''
	print '\nPlotting spectra...'
	colours=['k','r','b','g','m','c','y']
	axisfontsize=14

	spectrumfigure=p.figure(figsize=(6,max([6,2*n_spectra])))
	for specnum in range(0,n_spectra):
		# subplotnum=100*n_spectra+21+specnum
		sub_plot=p.subplot2grid(shape=(n_spectra,1),loc=(specnum,0))
		if n_spectra<=len(colours):
			p.plot(energies[:-1],pc_array[specnum],marker='x',color=colours[specnum])
		else:
			p.plot(energies[:-1],pc_array[specnum],marker='x')
		sub_plot.set_xscale('log')
		sub_plot.xaxis.set_major_locator(FixedLocator([0.1,0.2,0.5,1,2,5,10,20,50,100]))
		sub_plot.xaxis.set_minor_locator(FixedLocator([0.3,0.4,0.6,0.7,0.8,0.9,3,4,6,7,8,9,15,30,40,60,70,80,90]))
		sub_plot.xaxis.set_major_formatter(ScalarFormatter())
		p.xlim(min(energies),max(energies))
	p.xlabel('Energy (KeV)',fontsize=axisfontsize)
	p.ylabel('Normalised Count Rate',fontsize=axisfontsize)
	p.show()

def write_spectra(output_stem,n_spectra,pc_array,energies,pc_errors=None):
	print '\nWriting PC spectra...'
	for specnum in range(0,n_spectra):
		specfilename=output_stem+'_spectrum_%s.csv' % str(specnum+1)
		specfile=open(specfilename,'w')
		if pc_errors != None:
			row0='E%s, Counts%s, +- \n' % (str(specnum+1),str(specnum+1))
		else:
			row0='E%s, Counts%s\n' % (str(specnum+1),str(specnum+1))
		specfile.write(row0)
		for i in range(0,len(energies)-1):
			if pc_errors != None:
				row=str(10.**((n.log10(energies[i])+n.log10(energies[i+1]))/2.))+', '+str(pc_array[specnum][i])+', '+str(pc_errors[specnum][i])+'\n'
			else:
				row=str(10.**((n.log10(energies[i])+n.log10(energies[i+1]))/2.))+', '+str(pc_array[specnum][i])+'\n'
			specfile.write(row)


def write_eigenvalues(output_stem,eigenvalues,eigenval_errors):
	print '\nWriting eigenvalues...'
	eigenfilename=output_stem+'_eigenvalues.csv'
	eigenfile=open(eigenfilename,'w')
	if eigenval_errors != None:
		row0='n, eigenvalue, +-\n'
	else:
		row0='n,eigenvale\n'
	eigenfile.write(row0)
	for i in range(0,len(eigenvalues)):
		if eigenval_errors != None:
			row=str(i)+', '+str(eigenvalues[i])+', '+str(eigenval_errors[i])+'\n'
		else:
			row=str(i)+', '+str(eigenvalues[i])+'\n'
		eigenfile.write(row)


def calc_vector_orientations(n_spectra,pc_array,means,energies,ebounds):
	'''Calculates the orientation of each principal component on a count-count plot,
	where the energy bands are defined by abounds'''

	print '\n Calculating vector orientations...'

	if len(ebounds)!=2:
		print 'Error! ebounds should be a list containing tuples, to describe the energy bounds'
		return 0

	for specnum in range(0,n_spectra):
		pc_spectrum=pc_array[specnum]
		pc_spectrum=[a*b for a,b in zip(pc_spectrum,means)]
		vals=[]

		for emin,emax in ebounds:

			indexmin=min(range(len(energies)), key=lambda i: abs(energies[i]-emin))
			indexmax=min(range(len(energies)), key=lambda i: abs(energies[i]-emax))

			vals.append(sum(pc_spectrum[indexmin:indexmax]))

		#Angle from vertical
		angle=n.arctan(vals[0]/vals[1])*180/n.pi
		print [v/vals[1] for v in vals], angle
