from hcipy import *
import numpy as np
from matplotlib import pyplot as plt



if __name__ == "__main__":


        amp = read_fits('ELT_pupil_downsampled_Lp.fits')
        phase = read_fits('METIS_vAPP_PDR_L.fits')

        aper = amp * np.exp(1j * phase)

        N = 512
        D = 39.2 # meters

        npix = 101 # LM detector size (should be 2048)
        pscale = 0.0052 # LM IMG pixel scale arcsec/pix

        wavelength_0 = 3E-6     # reference wavelength


        pixels_per_lambda = 4   # number of pixels per lambda over D at reference wavelength
        number_of_airy_rings = 20       # Field of view in number of airy rings of the reference wavelength

        # This variable is used to determine the scale of a pixel in the simulation (effectively it is focal_length/D the fnumber of course...)
        telescope_focal_length = 100 # meters

        # Make the Telescope
        #pupil_grid = make_pupil_grid(N, D) = circular_aperture(D)(pupil_grid)
        pupil_grid = make_pupil_grid(aper.shape)

        aper = Field(aper.ravel(), pupil_grid)

        # This detector grid is in the units of telescope_focal_length
        detector_grid = make_focal_grid(pupil_grid, wavelength=wavelength_0, q=pixels_per_lambda, num_airy=number_of_airy_rings, focal_length=telescope_focal_length)
        prop = FraunhoferPropagator(pupil_grid, detector_grid, wavelength_0=wavelength_0, focal_length=telescope_focal_length)

        # Make the on-sky grid
        pixel_scale = 1000 # arcseconds per meter, because the telescope focal length is in meters.
        on_sky_grid = detector_grid.scaled(pixel_scale)

        # Propagate different wavelengths
        wavelengths = np.linspace(0.5E-6, 2E-6, 11)
        for wavelength in wavelengths:
                wf = Wavefront(aper, wavelength)
                wf_foc = prop.forward(wf)

                plt.clf()
                # Plot it with the on_sky_grid coordinates
                imshow_field(np.log10(wf_foc.power/wf_foc.power.max()), on_sky_grid, vmin=-5)
                plt.xlabel('RA (arcsec)')
                plt.ylabel('DEC (arcsec)')
                plt.draw()
                plt.pause(0.001)
