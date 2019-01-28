from hcipy import *
import numpy as np
from matplotlib import pyplot as plt



if __name__ == "__main__":


        amp = read_fits('METIS_APP_20jul2018_amp.fits')
        phase = read_fits('METIS_APP_20jul2018_phase.fits')


#        amp = (np.abs(phase)>0)
#        write_fits(amp.astype(int), "L_ELT_binary_pupil.fits")
        aper1 = amp * np.exp(-1j * phase)
        aper2 = amp * np.exp(1j * phase)
        aper3 = amp + 0j
 #       N = 512
        D = 39.2 # meters

        npix = 512 # LM detector size (should be 2048)
        pscale = 0.0052 # LM IMG pixel scale arcsec/pix

        wavelength_0 = 3E-6     # reference wavelength

        pixels_per_lambda = (206265. * wavelength_0 / D) / pscale

        number_of_airy_rings = npix / pixels_per_lambda / 2# should give npix as the FOV

        # This variable is used to determine the scale of a pixel in the simulation (effectively it is focal_length/D the fnumber of course...)
        pixel_size = 18e-6 # microns
        #telescope_focal_length = 100 # meters
        telescope_focal_length = 206265. / (pscale / pixel_size)

        # Make the Telescope
        pupil_grid = make_pupil_grid(aper1.shape, D)

        aper11 = Field(aper1.ravel(), pupil_grid)
        aper22 = Field(aper2.ravel(), pupil_grid)
        aper33 = Field(aper3.ravel(), pupil_grid)

        # This detector grid is in the units of telescope_focal_length
        detector_grid = make_focal_grid(pupil_grid, wavelength=wavelength_0, q=pixels_per_lambda, num_airy=number_of_airy_rings, focal_length=telescope_focal_length)
        prop = FraunhoferPropagator(pupil_grid, detector_grid, wavelength_0=wavelength_0, focal_length=telescope_focal_length)

        # Make the on-sky grid
        pixel_scale = (pscale / pixel_size) # arcseconds per meter, because the telescope focal length is in meters.
        on_sky_grid = detector_grid.scaled(pixel_scale)

        # Propagate different wavelengths
        wavelengths = np.linspace(3E-6, 5E-6, 3)
        for wavelength in wavelengths:

                # make psf A
                wf1 = Wavefront(aper11, wavelength)
                wf1_foc = prop.forward(wf1)
                wf2 = Wavefront(aper22, wavelength)
                wf2_foc = prop.forward(wf2)
                wf3 = Wavefront(aper33, wavelength)
                wf3_foc = prop.forward(wf3)

                fname = 'LMS_APP_PSF_{:.1f}.pdf'.format(wavelength*1e6)
                ti = 'LMS APP at {:.1f} microns'.format(wavelength*1e6)

                plt.clf()
                # Plot it with the on_sky_grid coordinates
                psf_tot = wf1_foc.power + wf2_foc.power + 0.03 * wf3_foc.power
                imshow_field(np.log10(psf_tot/psf_tot.max()), on_sky_grid, vmin=-5)
                plt.colorbar()
                plt.xlabel('RA (arcsec)')
                plt.ylabel('DEC (arcsec)')
                plt.title(ti)
                plt.draw()

                plt.savefig(fname)

                plt.pause(2)
