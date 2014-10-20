#!/usr/bin/env python

import os, sys

def clean_psf(psfname):
   psfroot = os.path.splitext(psfname)[0]
   os.system('rm %s_??_??.fits' % psfroot)


if __name__ == "__main__":
   psf = sys.argv[1]
   print psf
   clean_psf(psf)