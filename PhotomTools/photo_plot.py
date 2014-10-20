#! /usr/bin/python

#=======================================================================
#+
# NAME:
#   photo_plot.py
#
# PURPOSE:
#   Plot the SED and photometry using the .spec outputs of le phare.
#
# 
# INPUTS: 
#		.spec file containing (currently just) one SED spectrum best fit to our photometry
#
# OUTPUTS:
#   	Currently just outputs SED to screen, attempting to match the output of Brian's supermongo script, specBCLedit.sm
#
# BUGS:
#		I'd ideally like to use Ureka python, but it is about twice as slow as my default /usr/bin/python installation
#
# Needed Fixes:
# 		Should be an option whether Y axis should is in flux or AB magnitude. 
#		Make default opening window size/scaling appropriate
#		Reformat model values text to line up in columns
#		Fix what looks like an error bar in the magnitude on the non detections
#		Add P(z=7.42) onto probability plots somewhere	 
#		Treat case where we have spec-z. Will this even make a prob plot?
#
# 		Inevitably, be more general in treating an SED catalog.
#
#		Use functions, maybe in some other modules to make more readable
# REVISION HISTORY:
#   2014-05-06  started by Hoag (UC Davis)
#
#-
#=======================================================================


import numpy as np
import pyfits
import pylab 
import sys
import os

masterfile = sys.argv[1]
def main(masterfile):
	master = [x.partition("\n")[0] for x in open(masterfile).readlines() ]   # Contains a list of the lines in .spec file
	for index,line in enumerate(master):
		if index == 1:
			secondline = [float(x) for x in line.split(' ') if x and x !='#']
			if secondline[1] == secondline[2]:
				havespecz = 1
				specz = secondline[1]
				# print "have specz input=", specz
			else:
				havespecz = 0
				# print "don't have specz as input"
		while 1:
			if 'FILT' in line:
				nfilt = int(line.split(' ')[-1])	
				findex = index
				# print findex
				# print nfilt
			break
		while 1:	
			if 'STAR' in line:
				begindex = index
				firstprobindex = begindex + nfilt + 1
				lastprobindex = firstprobindex + 100  # will need to change if I change my redshift stepsize or range in config/*.para file in lephare_dev
				# print firstprobindex
			break
		
	# ---------------------------------------------- #
	# Get the probability data to plot (only if we don't have specz)
	# ---------------------------------------------- #

	if not havespecz:
		problist = master[firstprobindex:lastprobindex+1]
		zs = [float(0.1*x) for x in range(0,101)]   # horizontal axis of prob plot. Note that the range here will need to change if I change my redshift stepsize or range in config/*.para file in lephare_dev
		probs = [float(s.split(' ')[-1]) for s in problist] # vertical axis of prob plot

		# Figure out P(z>7) and P(z=7.42)
		sum=0
		norm=0
		for ind,z in enumerate(zs):
			norm+=probs[ind]
			if z>=7:
				sum+=probs[ind]
				if z==7.4:
					pz7pt4 = probs[ind]
		pzgt7=sum/norm
		pz7pt4 = pz7pt4/norm

	# ---------------------------------------------- #
	# Get the model spectrum data to plot
	# ---------------------------------------------- #

	mindex = findex + 3
	modelstring = master[mindex]
	quantstring = master[mindex+1]
	modellist = [x for x in modelstring.split(' ') if x and x != '#']
	quantlist = [x for x in quantstring.split(' ') if x]
	mdict = { m: q for m,q in zip(modellist,quantlist) }
	z = 'Zphot' # won't matter for the case of zspec 

	# We only want the spectrum of the most likely galaxy template
	nlinesgal1 = int(mdict['Nline'])
	# print master[firstspec1index]
	# print master[lastspec1index]
	if not havespecz:
		firstspec1index = lastprobindex+1
		# print firstspec1index
	else:
		firstspec1index = firstprobindex	
	lastspec1index = firstspec1index+nlinesgal1
	modelspeclist = master[firstspec1index:lastspec1index]
	# print modelspeclist
	# print [s.split(' ')[1] for s in modelspeclist]
	wave = [ s.split(' ')[1] for s in modelspeclist ]
	modelspec = [ float(s.split(' ')[3]) for s in modelspeclist ]   # Right now these are in AB magnitudes	

	# ---------------------------------------------- #
	# Get the photometric data + errors to plot
	# ---------------------------------------------- #

	photlist = master[begindex+1:begindex+1+nfilt]
	splitlist = [ s.split(' ') for s in photlist]
	photlist=[]

	for filt in splitlist:
		newlist=[]
		for field in filt:
			if field:
				newlist.append(field)
		photlist.append(newlist)

	# First make arrays for detections		
	maglist = [float(x[0]) for x in photlist if float(x[1]) > 0 ]
	errlist = [float(x[1]) for x in photlist if float(x[1]) > 0 ]
	centlist = [float(x[2]) for x in photlist if float(x[1]) > 0 ]
	fwhmlist = [float(x[3])/2 for x in photlist if float(x[1]) > 0 ]  

	# Next do non-detections

	limmaglist = [float(x[0]) for x in photlist if float(x[1]) < 0 ]
	limerrlist = [float(x[1]) for x in photlist if float(x[1]) < 0 ]
	limcentlist = [float(x[2]) for x in photlist if float(x[1]) < 0 ]
	limfwhmlist = [float(x[3])/2 for x in photlist if float(x[1]) < 0 ]

	# print limerrlist

	# ---------------------------------------------- #
	# Plot all the data on a figure + subfigure (probs)
	# ---------------------------------------------- #
	
	if not havespecz:
		fg, axs = pylab.subplots(nrows=2, ncols=1)

		# SED plot
		ax = axs[0]
		ax1 = axs[1]
		ax1.patch.set_facecolor('black')
		ax1.plot(zs,probs,'w')
		# ax.set_ylim(0,1)
		# ax1.set_xlim(1,9)
		ax1.set_xlabel('Redshift, z')
		ax1.set_ylabel('Probability (Unnormalized)')
		# ax1.text(5,0.75, 'P(z>7) = %.2f \n P(z=7.4) = %.2f' % (float(pzgt7),float(pz7pt4)), color='w',horizontalalignment='right',
		     # verticalalignment='bottom')
		fg.suptitle("Best Fit SED to SExtractor ISO photometry scaled to F160W AUTO MAG", fontsize=14)
	else: 
		fg = pylab.figure()
		ax = pylab.subplot(111)
		fg.suptitle('Best Fit SED forced to z=%.2f with SExtractor ISO photometry scaled to MAG_AUTO_F160W' % specz)

	ax.set_xscale('log')
	ax.patch.set_facecolor('black')
	# ax.set_yscale('log')
	ax.plot(wave,modelspec,'b')
	ax.set_ylim(30,18)
	ax.set_xlim(2000,100000)
	ax.set_xlabel('Wavelength (A)')
	ax.set_ylabel('Magnitude (AB)')
	# ax.xaxis.label.set_color('red')
	ax.tick_params(axis='x', colors='red')
	# ax.xaxis.set_ticks(np.arange(1e3,1e5,1e4))
	ax.set_xticks([1e3,1e4,1e5])
	# ax.get_xaxis().set_major_formatter(pylab.ticker.ScalarFormatter())
	ax.xaxis.set_ticks_position('bottom')
	# ax.errorbar(centlist,maglist,xerr=fwhmlist,yerr=errlist,marker='o',markersize=9,linestyle='',color='w')
	ax.errorbar(centlist,maglist,xerr=fwhmlist,yerr=errlist,marker='o',markersize=9,linestyle='',color='w')
	ax.errorbar(limcentlist,limmaglist,xerr=limfwhmlist,yerr=limerrlist,lolims=True,marker='o',markersize=9,linestyle='',color='w')
	# ax.arrow( x, y, dx, dy, **kwargs )
	# ax.invert_yaxis()
	# ax.axvline(x=10245,linestyle='--',color='w',label='Emission line seen in MOSFIRE Y-band')

	# ax.text(0.3,0.95, r'Type # Model $N_{band}$ $\chi_{\nu}^2$ Z E(B-V) $log(L_{IR})$ $log(Age)$ $log(SFR)$ $log(M_{\odot})$ ''\n'' %s %s %s %.2f %.2f %.2f %.2f %.2f %.2f %.2f  ' \
	# 	% (mdict['Type'],mdict['Model'],mdict['Nband'],float(mdict['Chi2']),float(mdict[z]),float(mdict['EB-V']),float(mdict['Lir']),float(mdict['Age']),float(mdict['SFR']),float(mdict['Mass'])),
	#      horizontalalignment='center',
	#      verticalalignment='center',
	#      transform = ax.transAxes,color='w',fontsize = 12)
	ax.text(0.3,0.95, r'$\chi_{\nu}^2$ $\,\,$ Z $\,\,$ E(B-V) $\,\,$ $log(L_{IR})$ $\,\,$ $log(Age)$ $\,\,$ $log(SFR)$ $\,\,$ $log(M_{\odot})$ ''\n'' %.2f $\,\,$ %.2f $\,\,$ %.2f $\,\,$ %.2f $\,\,$ %.2f $\,\,$ %.2f $\,\,$ %.2f  ' \
		% (float(mdict['Chi2']),float(mdict[z]),float(mdict['EB-V']),float(mdict['Lir']),np.log10(float(mdict['Age'])),float(mdict['SFR']),float(mdict['Mass'])),
	     horizontalalignment='left',
	     verticalalignment='center',
	     transform = ax.transAxes,color='w',fontsize = 12)

	# Probability plot 


	
	# ax.legend()
	pylab.show()
main(masterfile)