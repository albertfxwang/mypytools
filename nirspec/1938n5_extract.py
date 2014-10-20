from nirspec import nirspec_extract as ne
import matplotlib.pyplot as plt

ne.nir_extract('1938n5_00.fits','1938n5_00.spec',270,1065,8,238,apmin=-7,
               apmax=7)
plt.show()
ne.nir_extract('1938n5_01.fits','1938n5_01.spec',270,1065,8,238,apmin=-7,
               apmax=7)
plt.show()
