#!/usr/bin/env python
# Version 1.1r KH Huang - Re-wrote for galfit 3.0 and added x, y position in the output for 
#                         sersic type objects (on royal cluster)
import os, glob, sys, string, shutil, time
from math import *
from string import *
from pyraf import iraf
from iraf import noao, stsdas, fourier, fitting, artdata, images
from pygoods.sextutils import *
import readgfheader

## NOTE: GALFIT output images should have names obj*_out.fits!

def run(funcname,outname):
    # output as SExtractor format
    #iraf.ls('obj*_out.fits',Stdout='tmplist')
    outlist = glob.glob('obj*_out.fits')
    f = open('tmplist','wb')
    for t in outlist:
      f.write('%s\n' % t)
    f.close()
    if funcname=="sersic":
       f2=open(outname,"w")
       colnames = ['OBJNUM','XCENT','XCENTERR','YCENT','YCENTERR','MAG',
                   'MAGERR','RE','REERR','N','NERR','AR','ARERR','PA','PAERR',
                   'CHISQ','CHISQNU']
       writeheader(f2,colnames)
    #if funcname=="devauc":
    #   f2=open(outname,"w")
    #   colnames=['#OBJNUM','XCENT1','YCENT1','MAG1','MAG1ERR','RE1','RE1ERR',
    #              'AR1','AR1ERR','PA1','PA1ERR','XCENT2','YCENT2','MAG2',
    #              'MAG2ERR','RS2','RS2ERR','AR2','AR2ERR','PA2','PA2ERR'
    #              'CHISQ','CHISQNU']
    #   writeheader(f2,colnames) 
    f=open('tmplist')
    objnum = []
    xcent = []
    xcenterr = []
    ycent = []
    ycenterr = []
    mag = []
    magerr = []
    re = []
    reerr = []
    sersicn = []
    sersicnerr = []
    axratio = []
    axratioerr = []
    pa = []
    paerr = []
    chisq = []
    chisqnu = []
    for line in f.readlines():
        id = line.split()[0]
        imgroot=id+"[2]"
        comp1=funcname
        if comp1=="sersic":
          iraf.imgets.saveParList(filename="imgets.par")
	  # figure out how many components there are
          n = 1
          ifcomp = 1
          while ifcomp:
              compstr = 'COMP_%d' % n
              iraf.imgets(imgroot,compstr)
              if iraf.imgets.value == 'sersic':    # nth fitted component found 
                  xcstrn = '%d_XC' % n       # x center
                  iraf.imgets(imgroot,xcstrn)
                  hxcentn=str(iraf.imgets.value)   
                  hxcentn=hxcentn.split('+/-')
                  if len(hxcentn) == 1:
                      xcentn = hxcentn[0]
                      xcenterrn = '1.e-6'    # there's no x center error information
                  else:
                      xcentn,xcenterrn=hxcentn
                  xcent += [xcentn]
                  xcenterr += [xcenterrn]
                  ycstrn = '%d_YC' % n       # y center
                  iraf.imgets(imgroot,ycstrn)
                  hycentn=str(iraf.imgets.value)
                  hycentn=hycentn.split('+/-')
                  if len(hycentn) == 1:
                      ycentn = hycentn[0]
                      ycenterrn = '1.e-6'    # there's no y center error information
                  else:
                      ycentn,ycenterrn=hycentn
                  ycent += [ycentn]
                  ycenterr += [ycenterrn]
                  magstrn = '%d_MAG' % n     # mag
                  iraf.imgets(imgroot,magstrn)
                  hmagn=str(iraf.imgets.value)
                  magn,magerrn=hmagn.split('+/-')
                  mag += [magn]
                  magerr += [magerrn]
                  restrn = '%d_RE' % n     # effective radius r_e
                  iraf.imgets(imgroot,restrn)
                  hren=str(iraf.imgets.value)
                  ren,reerrn=hren.split('+/-')
                  re += [ren]
                  reerr += [reerrn]
                  sersicnstrn = '%d_N' % n     # sersic n
                  iraf.imgets(imgroot,sersicnstrn)
                  hsersicnn=str(iraf.imgets.value)
                  sersicnn,sersicnerrn=hsersicnn.split('+/-')
                  sersicn += [sersicnn]
                  sersicnerr += [sersicnerrn]
                  axratiostrn = '%d_AR' % n     # axial ratio
                  iraf.imgets(imgroot,axratiostrn)
                  haxration=str(iraf.imgets.value)
                  axration,axratioerrn=haxration.split('+/-')
                  axratio += [axration]
                  axratioerr += [axratioerrn]
                  pastrn = '%d_PA' % n     # position angle
                  iraf.imgets(imgroot,pastrn)
                  hpan=str(iraf.imgets.value)
                  pan,paerrn=hpan.split('+/-')
                  pa += [pan]
                  paerr += [paerrn]
                  iraf.imgets(imgroot,'CHISQ')    # chi-square (same for the same image)
                  hchisq=str(iraf.imgets.value)
                  chisq += [hchisq]
                  iraf.imgets(imgroot,'CHI2NU')   # reduced chi-square (same for the same image)
                  hchisqnu=str(iraf.imgets.value)
                  chisqnu += [hchisqnu]
                  iraf.imgets(imgroot,'OBJNUM')   # object number (same for the same image!)
                  hobjnum=str(iraf.imgets.value)
                  objnum += [hobjnum]

                  # process other attributes...
                  n += 1     # jump to next possible component
              else:
                  ifcomp = 0
          # the last n will be 1 greater than the number of components
          # now start to print out the entries
          for i in range(n-1):
              f2.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n"
                        %(objnum[i],xcent[i],xcenterr[i],ycent[i],ycenterr[i],mag[i],magerr[i],
                          re[i],reerr[i],sersicn[i],sersicnerr[i],axratio[i],axratioerr[i],
                          pa[i],paerr[i],chisq[i],chisqnu[i]))
    f.close()
    f2.close()
      
    iraf.delete(files="tmplist", verify="no")
