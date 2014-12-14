#!/usr/bin/env python

import numpy as np
import os

"""
Quick functions to compare training set results between various people.
"""

trainDIR = "/Users/khuang/Dropbox/Github/glass/docs/gig_inspectionresults/training"
mine = trainDIR + '/GiGtraining141212_huang.txt'


def load_results(datafile):
   GiGoutput = os.path.join(trainDIR, datafile)
   data = np.genfromtxt(GiGoutput,comments='#',skip_header=2,names=True)
   return data

def print_emline_objs(data):
   selection = data[np.logical_or((data['G102_Emission_Line'] == 1),(data['G141_Emission_Line'] == 1))]
   print '\nIDs and PA of objects with at least 1 emission line:'
   for ii in xrange(len(selection)):
      print 'ID =',str("%.5d" % selection['ID'][ii]),' PA =',str("%.3d" % selection['PA'][ii])
   return selection

def print_multi_emline_objs_141(data):
   selection = data[(data['G141_Emission_Lines_Multiple'] == 1)]
   print '\nIDs and PA of objects with multiple emission lines in G141:'
   for ii in xrange(len(selection)):
      print 'ID =',str("%.5d" % selection['ID'][ii]),' PA =',str("%.3d" % selection['PA'][ii])
   return selection

def print_awesome_objs(data):
   selection = data[(data['AWESOME'] == 1)]
   print '\nIDs and PA of AWESOME objects:'
   for ii in xrange(len(selection)):
      print 'ID =',str("%.5d" % selection['ID'][ii]),' PA =',str("%.3d" % selection['PA'][ii])
   return selection

def get_emline_id(data):
   selection = data[np.logical_or((data['G102_Emission_Line'] == 1),(data['G141_Emission_Line'] == 1))]
   emline_id = np.unique([int(x) for x in selection['ID']])
   return emline_id

