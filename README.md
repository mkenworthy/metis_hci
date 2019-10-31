# METIS High Contrast Imaging - gvAPP coronagraphs

This repository contains the graving vector Apodizing Phase Plate (gvAPP) designs for the METIS PDR.

The subdirectories contain the phase and amplitude masks for the IMG-APP and LMS-APP, along with
the METIS masked ELT pupil used as a base for their design. Python scripts in the directories can generate
the PSF in the METIS IMG and LMS focal planes for any wavelength and any bandwidth.

The `transmission_curve` directory contains an estimate for the transmission of the gvAPP optics that includes all absorption and reflection losses.

The python scripts require `astropy` and `hcipy` which is available on github.

The designs are by David Doelmand and Emiel Por, and this page is maintained by Matthew Kenworthy
