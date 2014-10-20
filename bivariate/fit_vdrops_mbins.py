#!/usr/bin/env python

from numpy import *
import bivariate_lf as bl
import bivariate_fit as bf
import fit_lbg as fl
import mlutil
import zdist
import os, sys, time
from multiprocessing import Queue, Process, Pool

## initialize necessary things here
parv = array([-1.68527018,-20.50967054,0.8757289,0.85187255,0.26594964])
dlimits1 = array([[24.0,25.0],[-2.0,3.0]])
dlimits2 = dlimits1.copy()
mlimits1 = array([[21.0,26.5],[-2.0,3.0]])
mlimits2 = array([[23.0,28.5],[-2.0,3.0]])
kgrid1 = mlutil.readkgrid('kernel_Z.p')
kgrid2 = mlutil.readkgrid('kernel_Z_udf.p')
zdgrid1 = zdist.read_zdgrid('zdgrid_vdrops.p')
zdgrid2 = zdist.read_zdgrid('zdgrid_vdrops_udf.p')
mc = bl.mconvert('M1500_to_z.txt')
logr0_arr = arange(0.7, 1.0, 0.005)
sigma_arr = arange(1.0, 1.3, 0.005)
logr0_grid, sigma_grid = meshgrid(logr0_arr, sigma_arr)
logr0_grid, sigma_grid = map(ravel, [logr0_grid, sigma_grid])
#print logr0_grid
nprocs = 3
#chunksize = float(niter) / nproc
q_logl = Queue()
q_pars = Queue()
par = parv.copy()
N = len(logr0_arr)*len(sigma_arr)
logl = zeros(N)

chunksize = float(N) / float(nprocs)
procs = []
q_pars = Queue()
q_logl = Queue()

def worker(logr0_wk, sigma_wk, q_logl_wk, q_pars_wk):
	# given an array of logr0 and sigma, calculate the logl
	#print len(logr0_wk)
	logl_wk = zeros(len(logr0_wk))
	for i in range(len(logr0_wk)):
		par[2] = logr0_wk[i]
		par[3] = sigma_wk[i]
		logl_wk[i] = fl.quickmlfunc(par, dlimits1, dlimits2, kgrid1, kgrid2, zdgrid1, zdgrid2, mc,
			drop='v')
		#print logl_wk[i]
	#print "(logr0, sigma, logl):", logr0, sigma, l
		
	q_pars_wk.put([logr0_wk, sigma_wk])
	q_logl_wk.put(logl_wk)

t1 = time.time()
# Using pool or workers
#po = Pool(processes=1)  # create a pool with however many processors there are
#for i in range(len(logr0_arr)):
#	for j in range(len(sigma_arr)):
#		logr0 = logr0_arr[i]
#		sigma = sigma_arr[j]
#		q_pars.put([logr0, sigma])
#		res = po.apply_async(worker, (logr0, sigma))
#		q_logl.put(res.get())

# Create my own "pool" and use queues
for i in range(nprocs):
	i0 = int(round(chunksize)*i)
	i1 = int(round(chunksize)*(i+1))
	if i != nprocs-1:
		logr0_wk = logr0_grid[i0:i1]
		sigma_wk = sigma_grid[i0:i1]
	else:
		logr0_wk = logr0_grid[i0:]
		sigma_wk = sigma_grid[i0:]
	#print logr0_wk
	p = Process(target=worker, args=(logr0_wk,sigma_wk,q_logl,q_pars))
	procs += [p]
	p.start()
		
# collect results
# Use pools
#k = 0
#for i in range(len(logr0_arr)):
#	for j in range(len(sigma_arr)):
#		par_list += [q_pars.get()]
#		logl[k] = q_logl.get()
#		k += 1

	
# Use processes
logr0_list = zeros(0)
sigma_list = zeros(0)
logl_list = zeros(0)
for i in range(nprocs):
	#print "q_pars.get():", q_pars.get()
	#print "q_logl.get():", q_logl.get()
	l = q_pars.get()
	#print l[0]
	logr0_list = concatenate((logr0_list, l[0]))
	sigma_list = concatenate((sigma_list, l[1]))
	logl_list = concatenate((logl_list,q_logl.get()))

print len(logr0_list), len(logl_list), N
for i in range(nprocs):
	procs[i].join()

#po.close()
#po.join()   # wait for all processes to finish
t2 = time.time()

print "Total time: %.2f seconds." % (t2-t1)

f = open('vdrops_mbins_24_25.txt','a')
#f.write('# 1 logr0\n')
#f.write('# 2 sigma\n')
#f.write('# 3 logl\n')
for i in range(len(logl)):
	f.write('%f %f %f \n' % (logr0_list[i], sigma_list[i], logl_list[i]))
f.close()
	
