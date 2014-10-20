#!/usr/bin/env python

import os, glob

def countgal(dirs, root, bands, nperrun=40):
   nsim = {}
   for b in bands:
      nsim[b] = 0
   for d in dirs:
      if d[-1] == '/': d = d[:-1]
      for b in bands:
         ncat = len(glob.glob(d+'/'+root+'_'+b+'/run*.cat'))
         nsim[b] += ncat * nperrun
   for b in bands:
      print "Number in %s-band is %d" % (b, nsim[b])

#bands = ['B','V','I','Z']
#for b in bands:
#    cats = glob.glob('/data/goods202/khuang/goodssim/run*m_%s/run*.cat' % b)
#    Ngal = 0
#    for c in cats:
#        f = open(c)
#        lines = f.readlines()
#        f.close()
#        Ngal += (len(lines)-59)
#    print "%s-band has %d recovered objects" % (b,Ngal)
