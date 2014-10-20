# Kennicutt relations  Kennicutt 1998 ARAA, 36, 189

from numpy import *

def sfr_from_uv(lnu):
    return 1.4e-28*lnu # lnu in erg s-1 Hz-1

def sfr_from_abmag(ab):
    lnu = 10.**((ab+48.6)/-2.5)
    return sfr_from_uv(lnu)

def sfr_from_halpha(lhalpha):
    return 7.9e-42*lhalpha # loii in erg/s

def sfr_from_oii(loii):
    return 1.4e-41*loii # loii in erg/s

def sfr_from_fir(lfir):
    return 4.5e-44*lfir # fir in erg/s

def uv_from_sfr(sfr):
    return sfr/1.48e-28

def abmg_from_sfr(sfr):
    lnu = uv_from_sfr(sfr)
    return -2.5*log10(lnu)-48.6

def halpha_from_sfr(sfr):
    return sfr/7.9e-42

def oii_from_sfr(oii):
    return sfr/1.4e-41

def fir_from_sfr(oii):
    return sfr/4.5e-44
