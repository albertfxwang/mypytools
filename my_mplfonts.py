#!/usr/bin/env python

"""
Store my fonts for matplotlib.
"""
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
fontpath = '/Users/khuang/Library/Fonts'
sysfontpath = '/System/Library/Fonts'
rootfontpath = '/Library/Fonts'
customfp = '/Users/khuang/Dropbox/downloaded_fonts'

def Helvetica(s=12):
	return FontProperties(fname=customfp+'/helvetica/Helvetica.ttf', size=s)

def HelveticaLight(s=12):
   return FontProperties(fname=fontpath+'/HelveticaLight.ttf', size=s)

def OpenSans(s=12):
	return FontProperties(fname=fontpath+'/OpenSans-Regular.ttf',size=s)

def Default(s=12):
	return FontProperties(family='sans-serif',size=s)

def Menlo(s=12):
	return FontProperties(fname=sysfontpath+'/Menlo.ttc',size=s)

def AnonymousPro(s=12):
	return FontProperties(fname=fontpath+'/Anonymous Pro.ttf',size=s)

def Erika(s=12):
   return FontProperties(fname=fontpath+'/Erika Type.ttf',size=s)

def ErikaBold(s=12):
   return FontProperties(fname=fontpath+'/Erika Type_B.ttf',size=s)

def NouveauIBM(s=12):
   return FontProperties(fname=fontpath+'/Nouveau_IBM.ttf',size=s)

def SVBasicManual(s=12):
   return FontProperties(fname=fontpath+'/SVBasicManual.ttf',size=s)

def AmerType(s=12):
   return FontProperties(fname=rootfontpath+'/AmericanTypewriter.ttc',size=s)


#def Arial(s=12):
#   return FontProperties(fname=)