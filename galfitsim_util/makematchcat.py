#!/usr/bin/env python

from numpy import *
import os,glob
from pygoods import *
from pyraf import iraf
import pyfits

"""
To join the input object list *.list, SExtractor catalog run*.cat, and the GALFIT simulation output gfit*_s.cat.
Should include all input objects, whether they are detected and successfully measured by GALFIT or not.
Make FITS tables instead of ASCII catalogs.
"""
# columns related to GALFIT measurements
inputcolumns = ['x_in','y_in','mag_in','galaxy_type','re_in','axratio_in','pa_in']
incol_fmt_pf = ['D','D','D','I','D','D','D']  # format code for pyfits columns
gfcolumns = ['run_id','x_out','y_out','x_out_err','y_out_err','x_flag','y_flag','recovered_GALFIT',
           're_out','re_out_err','re_flag','mag_out','mag_out_err','mag_flag',
           'axratio_out','axratio_out_err','axratio_flag','n_out','n_out_err','n_flag','pa_out',
           'pa_out_err','pa_flag','chi2nu','ncomps']
gfcol_value0 = [0,-1.,-1.,-1.,-1.,-1,-1,0,
                -1.,-1.,-1,-1.,-1.,-1,
                -1.,-1.,-1,-1.,-1.,-1,-999.,
                -999.,-1,-1.,-1]  # initial value for numpy arrays
gfcol_fmt_np = ['i','d','d','d','d','i','i','i',
                'd','d','i','d','d','i',
                'd','d','i','d','d','i','d',
                'd','i','d','i']   # format code for numpy arrays
gfcol_fmt_pf = ['J','D','D','D','D','I','I','I',
                'D','D','I','D','D','I',
                'D','D','I','D','D','I','D',
                'D','I','D','I']   # format code for pyfits columns
# additional columns to be included from SE catalogs
secolumns = ['number','mag_auto','magerr_auto','ellipticity','theta_image','background','alpha_j2000','delta_j2000','flux_radius',
             'flux_radius_1','flux_radius_2','imaflags_iso','recovered_SE']
secol_value0 = [-1,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,0,0]
secol_fmt_np = ['i','d','d','d','d','d','d','d','d',
                'd','d','i','i']
secol_fmt_pf = ['I','D','D','D','D','D','D','D','D',
                'D','D','J','I']

def getbandresults(root,band,run):
  
  f = open(root+'_%s/gfit%d_s.cat' % (band,run))
  if len(f.readlines()) == 17: # header only
    f.close()
    return 0
  f.close()
  try:
    secat = sextractor(root+'_%s/run%d.cat' % (band,run))   # SE catalog with fake sources only
    #newcat = sextractor('run%d.newcat' % run)   # SE catalog including all sources in the image
    incat = sextractor(root+'_%s/gal%d.list' % (band,run))
    gfcat = sextractor(root+'_%s/gfit%d_s.cat' % (band,run))
  except: 
    print "Unable to merge matched catalog for iteration %d" % (run)
    return 1  
  #rootnum = root.split('run')[-1]
  #rootnum = int(rootnum.split('m')[0])  # e.g. run1m
  ntot = len(incat)
  # copy input columns 
  incols = {}
  incols['x_in'] = incat._1.copy()
  incols['y_in'] = incat._2.copy()
  incols['mag_in'] = incat._3.copy()
  incols['galaxy_type'] = where(incat._4=='expdisk',0,1)  # 1=devauc, 0=expdisk
  incols['re_in'] = incat._5.copy()
  incols['axratio_in'] = incat._6.copy()
  incols['pa_in'] = incat._7.copy()
  
  # initialize gfcolumns
  gfcols = {}
  gfcols['run_id'] = ones(ntot,'int')*run
  for l in range(1,len(gfcolumns)):
    col = gfcolumns[l]
    gfcols[col] = ones(ntot,gfcol_fmt_np[l]) * gfcol_value0[l]

  # initialize secolumns --- assume they are all float numbers
  secols = {}
  for m in range(len(secolumns)):
    col = secolumns[m]
    secols[col] = ones(ntot,secol_fmt_np[m]) * secol_value0[m]

  # Find matches between input and SE, GALFIT catalogs  
  #xin_se, yin_se = secat.vector_assoc, secat.vector_assoc_1
  #sedetect = zeros(ntot,'int')
  #seindex = ones(ntot,'int')*(-1) 
  # whether SE detects or not, 1=detected, 0=not detected
  
  # Find out which objects are detected in SE, and also
  # the indexes of these objects in SE catalog
  for i in range(ntot):  # i is the index in *.list
    j = -1; k = -1
    if (incat._1[i] in secat.vector_assoc) & (incat._2[i] in secat.vector_assoc_1):
      # determine the index in the SE catalog
      try:
        posmatch = (secat.vector_assoc==incat._1[i]) & (secat.vector_assoc_1==incat._2[i])
        j = arange(ntot)[posmatch==True][0]  # j is the index in run*.cat
      except:
        pass
      id_sex = secat.number[j]
      # determine the index in the GALFIT catalog gfit*_s.cat
      #filename = 'obj%d-wmask-out.fits' % id_sex  # reconstruct the GALFIT output image name
      filename = 'obj%d_out.fits' % id_sex
      if filename in gfcat.filename:
        k = arange(len(gfcat))[gfcat.filename==filename][0]  # the index in gfcat
      # Now start copying GALFIT results
      if j >= 0:
        secols['recovered_SE'][i] = 1
        for col in secolumns:
          if col != 'recovered_SE':
            secols[col][i] = getattr(secat,col)[j]
      if k >= 0:
        gfcols['recovered_GALFIT'][i] = 1
        gfcols['x_out'][i] = gfcat.xout[k]
        gfcols['x_out_err'][i] = gfcat.xout_err[k]
        gfcols['y_out'][i] = gfcat.yout[k]
        gfcols['y_out_err'][i] = gfcat.yout_err[k]
        gfcols['x_flag'][i] = gfcat.xflag[k]
        gfcols['y_flag'][i] = gfcat.yflag[k]
        gfcols['re_out'][i] = gfcat.reout[k]
        gfcols['re_out_err'][i] = gfcat.reout_err[k]
        gfcols['re_flag'][i] = gfcat.reflag[k]
        gfcols['mag_out'][i] = gfcat.magout[k]
        gfcols['mag_out_err'][i] = gfcat.magout_err[k]
        gfcols['mag_flag'][i] = gfcat.magflag[k]
        gfcols['axratio_out'][i] = gfcat.arout[k]
        gfcols['axratio_out_err'][i] = gfcat.arout_err[k]
        gfcols['axratio_flag'][i] = gfcat.arflag[k]
        gfcols['n_out'][i] = gfcat.nout[k]
        gfcols['n_out_err'][i] = gfcat.nout_err[k]
        gfcols['n_flag'][i] = gfcat.nflag[k]
        gfcols['pa_out'][i] = gfcat.paout[k]
        gfcols['pa_out_err'][i] = gfcat.paout_err[k]
        gfcols['pa_flag'][i] = gfcat.paflag[k]
        gfcols['chi2nu'][i] = gfcat.chi2nu[k]
        gfcols['ncomps'][i] = gfcat.ncomps[k]
  # Now start constructing pyfits columns
  fc = []
  for n in range(len(inputcolumns)):
    col = inputcolumns[n]
    fc += [pyfits.Column(name=col,array=incols[col],format=incol_fmt_pf[n])]
  for l in range(len(gfcolumns)):
    col = gfcolumns[l]
    fc += [pyfits.Column(name=col,array=gfcols[col],format=gfcol_fmt_pf[l])]
  for m in range(len(secolumns)):
    col = secolumns[m]
    fc += [pyfits.Column(name=col,array=secols[col],format=secol_fmt_pf[m])]
  coldefs = pyfits.ColDefs(fc)
  tbhdu = pyfits.new_table(coldefs)
  if os.path.exists(root+'_%s/gfit%d_matched.fits'%(band,run)):
    os.remove(root+'_%s/gfit%d_matched.fits'%(band,run))
  tbhdu.writeto(root+'_%s/gfit%d_matched.fits'%(band,run))
  print "Made gfit%d_matched.fits" % run


if __name__ == '__main__':
    # Usage: python makematchcat.py "run*m" I (band)
    curdir = os.getcwd()
    #f = open('/data/raid10/khuang/galfit_sims/galfitsim.dir')
    bands = sys.argv[2:]
    rt = sys.argv[1]
    #dirs = f.readlines()
    #f.close()
    print bands
    h = pyfits.open('si_bright23_seg.fits')
    seg = h[0].data
    h.close()
    for b in bands:
        rcat = glob.glob(rt+"_%s"%b)
        for rc in rcat:
            os.chdir(rc)
            root = rc.split('/')[-1]
            root = root.split('_')[0]
            #if not len(glob.glob('*matched*')):
            scat = glob.glob('*s.cat')
            for sc in scat:
                #os.system('pwd')
                l = sc.split('_')[0]
                run = l.split('gfit')[-1]
                run = int(run)
                #run = int(l[4:])
                getbandresults(root,b,run,seg)
            os.chdir(curdir)

