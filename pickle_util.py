#!/usr/bin/env python

import os
import cPickle
import gzip

"""
Some utility functions to pickle and compress objects.
"""

def load_pickle(fname):
   """
   Load a regular pickled file.
   """
   f = open(fname, 'rb')
   x = cPickle.load(f)
   f.close()
   return x

def dump_pickle(x, fname):
   """
   Dump an object to a regular pickle file.
   """
   f = open(fname, 'wb')
   cPickle.dump(x, f, 2)
   f.close()

def pickle_gzip(x, fname):
   """
   Pickle the object x, then save to a gzipped pickle file.
   """
   if os.path.splitext(fname)[1] != '.gz':
      fname = fname + '.gz'
   f = gzip.open(fname, 'wb')
   cPickle.dump(x, f)
   f.close()

def pickle2gzip(fname):
   """
   Convert a regular pickled file into a gzipped pickled file.
   """
   print "loading pickled file..."
   x = load_pickle(fname)
   fname_gzip = fname + '.gz'
   fp = gzip.open(fname_gzip, 'wb')
   print "dumping to gzipped file..."
   cPickle.dump(x, fp)
   fp.close()
   print "removing regular pickled file..."
   os.remove(fname)
   print "Done."

