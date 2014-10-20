from nirspec import niredux as nir

flats = ['jan09s0010.fits', 'jan09s0011.fits', 'jan09s0012.fits',
         'jan09s0013.fits', 'jan09s0014.fits']
stars = [
    # 'aug14s0018.fits',
    'aug14s0019.fits',
    # 'aug14s0020.fits',
    # 'aug14s0021.fits'
    ]
imgpairs = [
#    ['Raw/aug14s0054.fits','Raw/aug14s0055.fits'], # No signal in these
    ['jan09s0001.fits','jan09s0002.fits']
    ]
nir.niredux(flats,stars,imgpairs,'1938n5',userlow=15068.,userhigh=17890.)

