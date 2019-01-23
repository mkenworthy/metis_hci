from hcipy import *
import numpy as np
from matplotlib import pyplot as plt
import pyfits as pf
import scipy.ndimage as snd

ELT = read_fits('ELT_pupil_36m_11.1m_60cm_spiders_243px.fits')

aperture = make_obstructed_circular_aperture(37,11.1/37,6,0.6)

pupil_grid_xl = make_pupil_grid(4096,39.146).shifted((-0.07,-0.07))

aperture_xl = aperture(pupil_grid_xl)

aperture_final = Field(snd.morphology.binary_erosion(aperture_xl.shaped, iterations= 30).ravel(),pupil_grid_xl)

aperture_final = subsample_field(aperture_final,4096//256)

aperture_final = aperture_final.shaped.T.ravel()


pupil_grid2 = make_pupil_grid(243)

aperture_cutout=Field(aperture_final.shaped[7:-6,7:-6].ravel(),pupil_grid2)

write_fits(aperture_final,'METIS_pupil_256_undersized.fits.gz')

write_fits(aperture_cutout,'METIS_pupil_243_undersized.fits.gz')


plt.figure()
imshow_field(aperture_cutout-Field(ELT.ravel(),pupil_grid2))
plt.show()



