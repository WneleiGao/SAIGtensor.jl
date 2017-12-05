import numpy as np
from vispy import io, plot as vp
import matplotlib.pyplot as plt
#dat_dir = '../../../'
fname = '/Users/wenlei/Desktop/seismic_cube.npy'
dat_file = fname
data = np.load(dat_file)
nILnXL, nt = data.shape

# Need to know:
nIL = 194  # number of inlines
nXL = 299  # number of crosslines
nt = 463   # number of samples per trace
dt = 0.004 # sample rate in seconds

data = np.reshape(data, newshape = (nIL, nXL, nt))
norm = max((np.abs(np.amax(data)), np.abs(np.amin(data))))

new_data = data.astype('float16')

# A sloppy way to manipulate amplitude values so they
# go from 0-255 as integers to display on graphics card.
new_data /= -2 * norm
new_data = 0.5 + new_data
new_data = 255 * new_data
new_data = new_data[:, :, :450]
vis_data = new_data.astype('int16')

vol_data = np.flipud(np.rollaxis(vis_data, 2))  # data to send to renderer

###
### Make figure
###

fig = vp.Fig()
clim = [130,200]  # This is a hack to scale the colors to something pleasing
vol_pw = fig[0, 0]
vol_pw.volume(vol_data, clim=clim, cmap='grays')
fig.show(run=True)
