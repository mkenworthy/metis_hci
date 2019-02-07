from hcipy import *
import numpy as np
from matplotlib import pyplot as plt
import pyfits as pf


def phase_to_PSF(Keckaper,Keckphase,pupil_grid,focal_grid):
        wf = Wavefront(Keckaper*np.exp(1j*Keckphase))
        wf2 = Wavefront(Keckaper*np.exp(-1j*Keckphase))
        wf3 = Wavefront(Keckaper)

        prop = FraunhoferPropagator(pupil_grid, focal_grid)

        focalplane = 0
        for wl in np.linspace(1,1.075,1):
                wf.wavelength = wl
                wf2.wavelength = wl
                wf3.wavelength = wl
                focalplane += prop(wf).intensity #+prop(wf2).intensity +0.02*prop(wf3).intensity 
                print(wl)

        test = np.gradient(focalplane.shaped)
        ghost = (np.roll(test[0]+test[1],170,axis =0)).ravel()
        ghost = ghost/np.max(ghost)*np.max(focalplane)*0.01
        #(ghost)
        Strehl = prop(wf).intensity.max()/prop(wf3).intensity.max()
        print('Strehl ratio of plate (by w.f.s calc) is: {}'.format(Strehl))
        plt.figure()
        imshow_field(np.log10(focalplane/focalplane.max()),vmin = -5.0,vmax = 0)
        plt.colorbar()

        plt.figure()

        Keckphase +=np.pi
        Keckphase = Keckphase%(2*np.pi)

        imshow_field(Keckphase, cmap = 'RdBu')
        contourf_field(Keckaper, levels = [-1,0], colors = ['k','w'], alpha = 1.)
        plt.colorbar()

        focshape = focalplane.shaped.shape
        
        slicefoc = focalplane.shaped[:,focshape[0]//2]
        plt.figure(figsize=(10,5))
        #rint(focal_grid.x.reshape(1600,1600)[:,0])
        plt.plot(focal_grid.x.reshape(1600,1600)[0,:], np.log10(slicefoc/slicefoc.max()))
        plt.ylim(-7,0)
        plt.xlim(0,25)
        plt.ylabel('Normalised intensity')
        plt.xlabel('y [$\lambda/D$]')

        plt.draw()
        plt.savefig('IMG_LM_APP_raw_contrast.pdf')
        plt.show()





if __name__ == '__main__':

        ### Get holograms ###
        #names = np.loadtxt('coro.list',dtype = str, delimiter  = '\n')
        #print(names)
        name = 'L/METIS_vAPP_PDR_L.fits'
#       name = 'METIS_HR_19_100_19_16_90_2_rad_phase_upscaled_1024.fits'
        vAPP = (read_fits(name))
        Pupil = np.abs(vAPP)>0
        Npix = vAPP.shape[0]

        pupil_grid = make_pupil_grid(Npix,1)    
        focal_grid = make_focal_grid(pupil_grid, 8 ,100)
        aperture = circular_aperture(1)
        vAPP2 = Field(vAPP.ravel(),pupil_grid)
        angle = np.radians(45)
        grating = Field(pupil_grid.x*2*np.pi*13*np.cos(angle)+pupil_grid.y*2*np.pi*13*np.sin(angle),pupil_grid)
        vAPP3 = (vAPP2+grating)%(2*np.pi)*(Pupil.ravel())

        phase_to_PSF(Pupil.ravel()*aperture(pupil_grid),vAPP2%(2*np.pi),pupil_grid,focal_grid)





        
