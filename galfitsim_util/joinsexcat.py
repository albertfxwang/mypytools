#!/usr/bin/env python

from numpy import *
from pygoods import *


def joinsexcat(cat1,cat2,joinedcat):
    """Join cat1 and cat2, where cat1 has only fake galaxies,
       and cat2 has all galaxies. Identify fake galaxies in 
       cat2 by matching x and y. Write the fake galaxies on
       top in the joined catalog."""
    c1 = sextractor(cat1)
    c2 = sextractor(cat2)
    c2x = array(c2._colentries)[:,1]
    c2y = array(c2._colentries)[:,2]
    id = c2.number # ID from catalog 2 is what we need
    f = open(joinedcat,'w')
    f.write(c2._header)  # write headers
    matched_cat2 = []
    n1 = len(c1)
    for i in range(len(c1)):
        x1 = c1._colentries[i][1] # string format
        y1 = c1._colentries[i][2]
        crit = (c2x==x1)*(c2y==y1)
        if sum(crit): # if there is a match
            j = compress(crit,arange(len(c2)))[0]
            # the jth object in cat2 is a fake galaxy
            matched_cat2 += [j]
            f.write('%d ' % id[j])  # use the ID number from cat2
            for k in range(1,len(c2._colentries[0])):
                f.write('%s ' % c2._colentries[j][k])
            f.write('\n')
        else: # no match... it's an error
            raise ValueError, 'no match found'
    # Now write the real galaxies in the bottom
    for j in range(len(c2)):
        if j not in matched_cat2:
            n1 = n1+1
            f.write('%d ' % id[j])
            for k in range(1,len(c2._colentries[0])):
                f.write('%s ' % c2._colentries[j][k])
            f.write('\n')
    f.close()

