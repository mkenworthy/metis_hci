from hcipy import *
import numpy as np
import matplotlib.pyplot as plt

fname = 'METIS_APP_20jul2018_'

amp = read_fits(fname + 'amp.fits')
phase = read_fits(fname + 'phase.fits')

for wavelength in [3e-6, 5e-6]:
	ld = 1.5 / (360 * 60 * 60) * (2*np.pi) / (wavelength / 39)
	print(ld)

	pupil_grid = make_pupil_grid(amp.shape)
	focal_grid = make_focal_grid(pupil_grid, 8, 243/2)
	prop = FraunhoferPropagator(pupil_grid, focal_grid)
	focal_grid_view = make_focal_grid(pupil_grid, 16, ld)
	prop_view = FraunhoferPropagator(pupil_grid, focal_grid_view)

	img = prop_view(Wavefront(Field(amp.ravel() * np.exp(1j * phase.ravel()), pupil_grid)))
	psfs = []

	mask = (circular_aperture(16)(focal_grid_view) - circular_aperture(4)(focal_grid_view)) * (focal_grid_view.x > 2) + circular_aperture(4.5)(focal_grid_view) + (1 - circular_aperture(16)(focal_grid_view))
	mask = mask > 1e-3

	l_D = (wavelength / 39) / (2*np.pi) * 360 * 60 * 60
	print(l_D)

	col = ['k', colors.red]

	grating_0 = -30

	plt.figure(1)
	for grating in [30]:
		for i, r in enumerate([1,1.5]):
			wf = Wavefront(Field(amp.ravel(), pupil_grid))
			wf.electric_field *= np.exp(2j * np.pi * pupil_grid.y * (-12.5 - 2.2))
			wf.electric_field *= np.exp(2j * np.pi * pupil_grid.x * -grating)
			wf = prop.forward(wf)
			if r != -1:
				if grating == 0:
					wf.electric_field *= circular_aperture(2*r / l_D, center=[grating_0,0])(focal_grid)
				else:
					wf.electric_field *= circular_aperture(2*r / l_D)(focal_grid)#, center=[grating, -2.2])(focal_grid)
			
			if grating != 0:
				I = wf.intensity / wf.intensity.max()
				imsave_field('PSF_at_beamsplitter_%d_%dum.png' % (int(2*r), wavelength*1e6), np.log10(I), vmin=-5, mask=(I > 1e-20), mask_color='w')
			
			wf = prop.backward(wf)
			
			wf2 = wf.copy()
			wf2.electric_field *= np.exp(1j * phase.ravel())
			wf2 = prop_view.forward(wf2)
			
			I = wf2.intensity / wf2.intensity.max()
			
			if grating != 0:
				wf3 = wf.copy()
				wf3.electric_field *= np.exp(-1j * phase.ravel())
				wf3 = prop_view.forward(wf3)
				
				I += wf3.intensity / wf2.intensity.max()
				
				wf4 = wf.copy()
				wf4 = prop_view.forward(wf4)
				
				I += 1e-2 * wf4.intensity / wf4.intensity.max()
				
			imsave_field('gvAPP_img_post_beamsplitter_%d_%dum.pdf' % (int(2*r), wavelength*1e6), np.log10(I), vmin=-6)
			
			ind = np.unravel_index(np.argmax(I * circular_aperture(50)(focal_grid_view)), I.grid.shape)[1]
			rr = I.grid.y.reshape(I.grid.shape)[:,ind]
			p = I.reshape(I.grid.shape)[:,ind]
			rr -= rr[np.argmax(p)]

			if r == -1:
				plt.plot(rr * l_D, p, label='%dum - Design' % (wavelength*1e6), c=col[i])
			else:
				plt.plot(rr * l_D, p, label='%dum - With beamsplitter' % (wavelength*1e6), c=col[i])
	plt.yscale('log')
	plt.ylim(2e-7, 2)
	plt.xlim(-0.1, 0.7)
	plt.legend()
	plt.xlabel('Angular separation (arcsec)')
	plt.ylabel('Raw contrast')
	plt.savefig('APP_contrast_beamsplitter%dum.pdf' % (wavelength*1e6))
	plt.clf()
	#plt.show()
