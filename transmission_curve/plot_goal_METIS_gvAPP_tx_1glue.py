import numpy as np
import matplotlib.pyplot as plt

from astropy.io import ascii

# plot_goal_METIS_gvAPP_tx_1glue - make plot of estimated transmission of METIS gvAPP

# M Kenworthy

f = plt.figure(figsize=(10,5))

data = ascii.read('180327 reference transmittance.csv')
dataA = ascii.read('180327 FTIR transmittance A.csv')

wl = data['col1']/1000.                    # wavelength in microns
caf2_reflect = data['col2']/100.           # CaF2, single 2mm layer, no antireflection coatings (ARC)
lc_and_caf2_reflect = data['col3']/100.    # CaF2, single 2mm layer, no ARC, Liquid Crystal (LC)

wlA = dataA['col1']/1000.
eris_gvapp_A = dataA['col2']/100.          # 8mm caf2, 3 glue layers, LC layer - ERIS gvAPP

# I assume that the bulk absorption of CaF2 is minimal over all wavelength ranges

# to calculate the TX of one glue layer only:

lc_and_caf2_only = lc_and_caf2_reflect / caf2_reflect # 2mm CaF2 and LWP only

glue_3_only = eris_gvapp_A / lc_and_caf2_only

glue_1_only = np.power(glue_3_only, 1./3.)

# gvAPP for METIS should be 1 glue layer and one LC layer

metis_gvapp = lc_and_caf2_only * glue_1_only

# plot out all relevant curves to confirm method

plt.plot(wl, lc_and_caf2_only, label='LC layer') # TX of 2mm substrate with no reflection losses

plt.plot(wl, glue_3_only, label='3 glue') # TX of 2mm substrate with no reflection losses
plt.plot(wl, glue_1_only, label='1 glue') # TX of 2mm substrate with no reflection losses
plt.plot(wl, metis_gvapp, label='METIS gvAPP') # TX of 2mm substrate with no reflection losses

plt.legend()

plt.xlim([3.,5.])
plt.ylim([0.,1.])
plt.xlabel('Wavelength [microns]')
plt.ylabel('Transmission')
plt.title('METIS gvAPP - one glue layer and one ($d=18\mu m$) LQ layer')

plt.draw()
plt.show()

# make output plot for HCI METIS documents

f = plt.figure(figsize=(10,5))

plt.plot(wl, metis_gvapp, label='METIS gvAPP') # TX of 2mm substrate with no reflection losses

plt.plot(wl, lc_and_caf2_only, label='LC layer') # TX of 2mm substrate with no reflection losses

plt.legend()

plt.xlim([3.,5.])
plt.ylim([0.,1.])
plt.xlabel('Wavelength [microns]')
plt.ylabel('Transmission')
plt.title('METIS gvAPP - $d=18\mu m$ LC layer and one glue layer')

plt.draw()
plt.savefig('goal_METIS_gvAPP_tx_1glue.pdf', bbox_inches='tight')
