#!/usr/bin/env python

from numpy import *
import numpy as np

class robust(object):
	def __init__(self, x, k=1.5, mu0=None, sig0=None, eps=1.e-8):
		self.x = x  # data array
		self.k = k  # the scale of the Huber function
		if mu0 == None:
			mu0 = median(x)
		self.mu0 = mu0  # initial guess of location
		if sig0 == None:
			self.sig0 = median(abs(x-median(x)))  # initial guess of scale is MADN
		self.eps = eps
		
	def location_huber(self):
		return location_huber(self.x, k=self.k, eps=self.eps, sig=self.sig0)
		
	def loc_scale(self, maxiter=50):
		loc, scale, weights = loc_scale(self.x, self.mu0, sig0=self.sig0, 
		                                k=self.k, eps=self.eps, maxiter=maxiter)
		self.loc = loc
		self.scale = scale
		self.weights = weights 

def MADN(x):
	"""
	Calculate MADN (median absolute deviation).
	"""
	m = np.median(np.abs(x - np.median(x)))
	return m

def std_MADN(x):
	"""
	Estimate standard deviation from MADN.
	"""
	m = MADN(x)
	return 1.4826 * x
	
def W_huber(y, k=1.5):
	# calculates the weight function
	W = ma.where(ma.abs(y)>k, k/ma.abs(y), 1.0)
	return W

def location_huber(x, k=1.5, eps=0.01, sig=-1.0):
	# estimate the location mu given an array of data x
	# Use Huber's function family
	# k is the parameter in Huber function
	# eps is the convergence factor: mu_j+1 - mu_j <= eps * MADN
	if sig <= 0:
		sig = median(abs(x - median(x)))/0.6745  # use MADN as the initial estimate of scale
	mu_prev = median(x)   # initial guess of location
	y = (x - mu_prev) / sig   
	w = W_huber(y, k=k)  # calculate the weights
	mu_next = sum(w*x) / sum(w)   # the next estimate of location
	niter = 1
	while (abs(mu_next-mu_prev)>(eps*sig)):
		# if not converged yet, update mu_prev and keep going
		mu_prev = mu_next
		y = (x - mu_prev) / sig
		w = W_huber(y, k=k)
		mu_next = sum(w*x) / sum(w)
		niter += 1
		
	#print "Estimated location is %.2f within %d iterations" % (mu_next, niter)
	return mu_next
	
def W_bisquare_scale(x):
	seterr(divide='ignore')   # suppress the divide-by-zero warnings
	y = where(abs(x)>0, minimum(3.0-3.*x**2+x**4, 1./(x**2)), 3.)  # the weight at x=0 is 3 (or 6?)
	#y = minimum(3.0-3.*x**2+x**4, 1./(x**2))
	seterr(divide='print')    # reset the divide-by-zero warnings
	return y

def scale_bisquare(x, mu, sig0=-1., eps=0.01):
	# estimate the scale sigma given an array of data x and an estimate of location mu
	# use the bisquare weight function: W(x) = min{3-3x**2+x**4, x**-2}
	# a commonly-used value of delta is 0.5 (justification?)
	# if not given an initial estimate of scale, use MADN
	delta = 0.5
	if sig0 <= 0:
		sig0 = median(abs(x-median(x)))/0.6745
	# calculate the weights first
	w = W_bisquare_scale((x)/sig0)
	# calculate the next estimate of sigma
	n = len(x)
	sig1 = sqrt((1./(n*delta))*sum(w*(x**2)))
	niter = 1
	#print "sig1=", sig1
	while abs(sig1/sig0 - 1.) >= eps:  # test for convergence
		sig0 = sig1
		w = W_bisquare_scale((x)/sig0)
		sig1 = sqrt((1./(n*delta))*sum(w*(x**2)))
		#print "sig1=", sig1
		niter += 1
	#print "Estimated scale is %.2f within %d iterations." % (sig1, niter)
	# to asymptotically approach SD when x is normally distributed and delta=0.5
	# only works if x is centered at zero!
	return sig1/1.56   
	
def loc_scale(x, mu0, mask=None, sig0=-1., k=1.5, eps=0.01, maxiter=50):
	# iteratively estimate location AND scale
	# input array x should NOT be a masked array
	seterr(divide='ignore')
	if sig0 <= 0:
		sig0 = median(abs(x-median(x)))/0.6745  # use the non-masked version of x
	if mask == None:
		mask = zeros(len(x),'int')  # don't mask anything
	x = ma.array(x, mask=mask)   # create masked array...
	#print x.mask
	# in the below iterations, only use the non-masked values
	r = (x - mu0) / sig0
	w1 = W_huber(r, k=k) #; print w1.mask
	w2 = W_bisquare_scale(r)
	# w1 and w2 will also be masked array
	#print sum(w1)
	mu1 = sum(w1*x) / sum(w1)
	#n = len(x)
	n = x.count()   # number of elements not masked
	sig1 = sqrt((sig0**2/(n*0.5))*sum(w2*r**2))
	niter=1
	while (abs(mu1-mu0)>(eps*sig1)) & (niter<=maxiter):
		mu0 = mu1
		sig0 = sig1
		r = (x - mu0) / sig0
		w1 = W_huber(r, k=k) #; print w1
		w2 = W_bisquare_scale(r)
		mu1 = sum(w1*x) / sum(w1)
		sig1 = sqrt((sig0**2/(n*0.5))*sum(w2*r**2))
		niter+=1
		if niter == maxiter:
			print "Maximum iterations reached."
			break
	#print "niter=", niter
	seterr(divide='print')
	return mu1, sig1/1.56, w1   # also return the weights for location
	
	
	