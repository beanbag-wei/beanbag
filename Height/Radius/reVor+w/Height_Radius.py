from wrf import getvar, to_np, get_cartopy
from netCDF4 import Dataset
import numpy as np
import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs 
from palettable.lightbartlein.diverging import BlueDarkRed18_5, BlueOrange10_3
from palettable.colorbrewer.diverging import RdGy_11_r
from scipy.fftpack import fft
import copy
import metpy.calc as mpcalc
import matplotlib.ticker
import cmaps

def fun():
    """
    the global parameter
    """
    global numr, numth, subx, suby
    numr = 256     # number of points in polar grid in r direction
    numth = 512    # number of points in polar grid in theta direction

# GET THE AZIMUTHAL WIND AND RATIAL WIND IN CARTESIAN GRID
def polar_wind(ucar, vcar):
    """
    This function is used to calculate the azimuth and radii component of wind
    ucar, vcar are the u and v component of wind in cartesion ordinate
    Mr, Nr are the number of points in rectangle grids.
    """
    num = ucar.shape
    Mr = num[0]
    Nr = num[1]
    uth = np.zeros((Mr, Nr))
    ur = np.zeros((Mr, Nr))
    i = np.arange(-(Mr//2), Mr//2+0.001, 1)
    j = np.arange(-(Nr//2), Nr//2+0.001, 1)
    x, y = np.meshgrid(i, j, indexing='ij')
    th = np.arctan2(y, x)
    uth = np.multiply(vcar, np.cos(th)) - np.multiply(ucar, np.sin(th))
    ur = np.multiply(ucar, np.cos(th)) + np.multiply(vcar, np.sin(th))
    return uth, ur
    
# GET THE U,V COMPONENT OF WIND IN CARTESIAN GRID
def cartesian_wind(uth, ur, Mr, Nr):
	"""
	This function is used to calculate the u,v component of wind
	utheta, ur are the azimuthal and radial wind in Polar ordinate
	Mr, Nr are the number of points in rectangle grids.
	"""
	ucar = np.zeros((Mr, Nr))
	vcar = np.zeros((Mr, Nr))
	i = np.arange(-(Mr//2), Mr//2+0.001, 1)
	j = np.arange(-(Nr//2), Nr//2+0.001, 1)
	x, y = np.meshgrid(i, j, indexing='ij')
	th = np.arctan2(y, x)
	ucar = np.multiply(ur, np.cos(th)) - np.multiply(uth, np.sin(th))
	vcar = np.multiply(ur, np.sin(th)) + np.multiply(uth, np.cos(th))
	return ucar, vcar

## do bilinear interpolation
def InterToPolar(imR, xR, yR):
	"""
	interpolate to Polar ordinate
	"""
	xf = math.floor(xR)
	xc = math.ceil(xR)
	yf = math.floor(yR)
	yc = math.ceil(yR)
	if xf == xc and yc == yf:
		v = imR[yc, xc]
	elif xf == xc:
		v = imR[yf, xf] + (yR - yf) * (imR[yf, xc] - imR[yf, xf])
	elif yf == yc:
		v = imR[yf, xf] + (xR - xf) * (imR[yc, xf] - imR[yf, xf])
	else:
		A = np.mat([[xf, yf, xf*yf, 1], 
					[xf, yc, xf*yc, 1],
					[xc, yf, xc*yf, 1],
					[xc, yc, xc*yc, 1]])
		r = np.mat([imR[yf, xf], imR[yf, xc], imR[yc, xf], imR[yc, xc]]).transpose()
		a = np.linalg.inv(A) * np.float64(r)
		w = np.mat([xR, yR, xR*yR, 1])
		v = np.dot(w, a)
	return v

def InterToCart(imP, r, t, rMin, rMax, M, N, delR, delT):
	"""
	interpolate to Cartesian ordinate
	"""
	ri = (r - rMin) / delR
	ti = t/delT
	rf = math.floor(ri)
	rc = math.ceil(ri)
	tf = math.floor(ti)
	tc = math.ceil(ti)
	
	# these grids are between 0 and 512, it is in forth quadrant
	# so ti of these grids are greater than 511.0, correspondingly
	# the tc is 512. tc would out of the boundary, whose maximum is 511
	# 
	
	if tc > N-1:
		tc = 0

	if rf == rc and tc == tf:
		v = imP[rc, tc]
	elif rf == rc:
		v = imP[rf, tf] + (ti - tf)*(imP[rf, tc] - imP[rf, tf])
	elif tf == tc:
		v = imP[rf, tf] + (ri - rf)*(imP[rc, tf] - imP[rf, tf])
	else:
		A = np.mat([[rf, tf, rf*tf, 1],
					[rf, tc, rf*tc, 1],
					[rc, tf, rc*tf, 1],
					[rc, tc, rc*tc, 1]])
		z = np.mat([imP[rf, tf], imP[rf, tc], imP[rc, tf], imP[rc, tc]]).T
		a = np.linalg.inv(A) * np.float64(z)
		w = [ri, ti, ri*ti, 1]
		v = np.dot(w, a)
	return v

def CarToPolar(imR, rMin, rMax, M, N):
	"""
	M, N is the dimensize of polar grid
	M along the with M points along the r axis and N points along the theta axis
	the imR is the array need to convert
	"""
	num = imR.shape # get the dimension size of input array
	Mr = num[0] # latitude
	Nr = num[1] # longitude
	Om = (Nr-1) / 2
	On = (Mr-1) / 2
	
	sx = (Nr-1) / 2 # scale factors 
	sy = (Mr-1) / 2
	
	imP = np.zeros((M, N), dtype=float)
	
	delR = (rMax - rMin) / (M-1)
	delT = 2*math.pi / N
	for ri in range(0, M):
		for ti in range(0, N):
			r = rMin + ri * delR
			t = ti * delT
			x = r*math.cos(t)
			y = r*math.sin(t)
			xR = x*sx + Om
			yR = y*sy + On
			imP[ri, ti] = InterToPolar(imR, xR, yR)
	return imP

def PolarToCar(imP, rMin, rMax, Mr, Nr):
	"""
	converts polar grid to rectangular grid. 
	imP is the polar grid with M rows and N columns of data (double data 
	between 0 and 1). M is the number of samples along the radius from 
	rMin to rMax (which are between 0 and 1 and rMax > rMin). 
	Mr and Nr are the number of points in the rectangular domain, 
	respectively along the latitude and logitude axis.
	The center of the grid is assumed to be the origin for the polar
	co-ordinates, and half the width of the grid corresponds to r = 1.
	Bilinear interpolation is performed for points not in the imP grid and
	points not between rMin and rMax are rendered as zero. The output is a Mr
	x Nr grid (with double values between 0.0 and 1.0).
	"""
	imR = np.zeros((Mr, Nr))
	num = imP.shape # get the dimension size of input array
	M = num[0]
	N = num[1]
	
	Om = (Nr-1) / 2 # co-ordinates of the center of the grid
	On = (Mr-1) / 2 # Mr along latitude and Nr along longitude
	sx = (Nr-1) / 2 # scale factors 
	sy = (Mr-1) / 2
	
	delR = (rMax - rMin)/(M-1)
	delT = 2*math.pi/N
	
	for yi in range(0, Mr):
		for xi in range(0, Nr):
			x = (xi - Om) / sx
			y = (yi - On) / sy
			r = math.sqrt(x*x + y*y)
			if r >= rMin and r <= rMax:
				t = math.atan2(y, x)
				if t < 0:
					t = t + 2*math.pi
				imR[yi, xi] = InterToCart(imP, r, t, rMin, rMax, M, N, delR, delT)
	return imR

def Asymmetric(ori):
    """
    find the asymmetric part of meteorological. 
    The array of ori has two dimension and is in Polar ordinate,
    and the radii dimension is in the first dimension.
    return the asymmetric part, asy
    """
    asy = np.zeros(ori.shape, dtype=float)
    for ir in range(0, numr):
        ori_mean = np.mean(ori[ir, :])
        asy[ir, :] = ori[ir, :] - ori_mean
    return asy

def Fourier(x, Fs, ReturnFre):
    """
    do the positive fourier transformation
    """
    N = len(x) # get the number of points
    k = np.arange(0, N, 1) # create a vector from 0 to N-1
    T = N/Fs # get the frequency interval
    freq = k/T # create the frequency range
    
    X = fft(x)/N # Fourier transform and normalize the data
    
    # only want the first half of the FFT, since the remainder is the same as the first half
    cutoff = math.ceil(N/2) 
    X = copy.copy(X[:cutoff])
    freq = copy.copy(freq[:cutoff])
    return X[0:ReturnFre], freq[0:ReturnFre]

def WaveDomain(Rm, wavenum, numr, numth):
    """
    After transforming from catesian grid to polar grid, 
    doing Fourier transform
    Rm is the array need to do fourier analysis
    wavenum is the number need to store, often wavenum = 4
    """
    Rmwav = np.zeros((wavenum, numr, numth))
    for ir in range(0, numr):
        # doing Fourier transform focused on Rm
        RmFreqDom, freqrange = Fourier(Rm[ir, :], numth, wavenum) 
    
        mag = 2 * ((RmFreqDom.real**2 + RmFreqDom.imag**2) **0.5) # the magnitude of waaves
        phase = np.arctan2(RmFreqDom.imag, RmFreqDom.real) # the phase of waves
        for iwave in range(0, wavenum):    
            Rmwav[iwave, ir, :] = mag[iwave] * np.cos(np.multiply(freqrange[iwave],
                                                                  np.radians(np.multiply(360/numth, np.arange(0, numth, 1))))\
                                                                 + phase[iwave])
    return Rmwav

#################################################################################################
###################################### The main Function ########################################
fun()
xlim = 2401    # number of grid points in x direction
ylim = 1921    # number of grid points in y direction
xcen = 2172    # 118.73868°E  the center of the tornado (already minus 1)
ycen = 826     # 34.10570°N   the center of the tornado (already minus 1)
loncen = 118.73868
latcen = 34.10570
subx =  30   # 250
suby =  40
# dx = 0.049383 # x resolution in km
# dy = 0.049383 # y resolution in km 
rMin = 0 # define whihc part of circle converting to polar grid 
rMax = 1 
Time = "07:24:00"

# read data
ncfile = Dataset(r"/data/data_s11_3/whuang/data/LES_d05/wrfout_d05_2016-06-23_" + Time + ".nc")
XLAT = getvar(ncfile, "lat")
XLON = getvar(ncfile, "lon")

nlev = 18
xmin = int(max(0, xcen - subx/2))
ymin = int(max(0, ycen - suby/2))
xmax = int(min(xlim-1, xcen + subx/2))
ymax = int(min(ylim-1, ycen + suby/2))

# the domain need to draw
minlon = to_np(XLON[ymin:ymax+1, xmin:xmax+1]).flatten().min()
maxlon = to_np(XLON[ymin:ymax+1, xmin:xmax+1]).flatten().max()
minlat = to_np(XLAT[ymin:ymax+1, xmin:xmax+1]).flatten().min()
maxlat = to_np(XLAT[ymin:ymax+1, xmin:xmax+1]).flatten().max()

# set the domain of subset data
# because the lambert projection, the picture couldn't fill the full picture, subset data need to be small larger than the picture
# it need larger domain than the boundary of drawed picture.
expansion = 5 
data_subx = subx+2*expansion # the total number of subdata
data_suby = suby+2*expansion
# the latitude and longitude 
lat = XLAT[ymin-expansion:ymax+1+expansion, xmin-expansion:xmax+1+expansion]
lon = XLON[ymin-expansion:ymax+1+expansion, xmin-expansion:xmax+1+expansion]
# the minimum and maximum of subset data
xminData = int(max(0, xcen - subx/2 - expansion))
yminData = int(max(0, ycen - suby/2 - expansion))
xmaxData = int(min(xlim-1, xcen + subx/2 + expansion))
ymaxData = int(min(ylim-1, ycen + suby/2 + expansion))
minlonData = to_np(XLON[yminData:ymaxData+1, xminData:xmaxData+1]).flatten().min()
maxlonData = to_np(XLON[yminData:ymaxData+1, xminData:xmaxData+1]).flatten().max()

# The Subset of data
# read u and v
uVer = getvar(ncfile, "ua", timeidx=0)[0:nlev, ymin-expansion:ymax+1+expansion, \
									xmin-expansion:xmax+1+expansion]
vVer = getvar(ncfile, "va", timeidx=0)[0:nlev, ymin-expansion:ymax+1+expansion, \
									xmin-expansion:xmax+1+expansion]
wVer = getvar(ncfile, "wa", timeidx=0)[0:nlev, ymin-expansion:ymax+1+expansion, \
									xmin-expansion:xmax+1+expansion]
height_agl = getvar(ncfile, 'height_agl')[0:nlev, ymin-expansion:ymax+1+expansion, \
									xmin-expansion:xmax+1+expansion]
# height_agl[17, :, :].min() # 18th level . The maximum height of it is 3568.92871094m, 
                                        # min is 3476.81762695m

########################################################################
########################################################################
# Transform to Polar Grid and do Fourier
# height(AGL)
# vertical velocity
# relative vortivity
########################################################################
########################################################################
# convert the AGL height to Polar Grid
heightP = np.zeros((nlev, numr, numth))
for ilev in range(0, nlev):
	heightP[ilev, :, :] = CarToPolar(to_np(height_agl[ilev, :, :]), rMin, rMax, numr, numth)

wavenum = 6
wVerWavP = np.zeros((wavenum, nlev, numr, numth))
for ilev in range(0, nlev):
	wVerP = CarToPolar(to_np(wVer[ilev, :, :]), rMin, rMax, numr, numth)
	wVerWavP[:, ilev, :, :] = WaveDomain(wVerP, wavenum, numr, numth) # do the fourier transform

# calculate the relative vorticity
dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

revorVerWavP = np.zeros((wavenum, nlev, numr, numth))
for ilev in range(0, nlev):
	# calculate the vorticity
	# ?mpcalc.vorticity
	revorVer = mpcalc.vorticity(uVer[ilev, :, :], vVer[ilev, :, :], dx = dx, dy = dy)
	revorVerP = CarToPolar(to_np(revorVer), rMin, rMax, numr, numth)
	revorVerWavP[:, ilev, :, :] = WaveDomain(revorVerP, wavenum, numr, numth) # do the fourier transform

########################################################################
# draw the height-radius plane
# Be careful, when caculating, the sequence is E-N-W-S, as to the rotation direction of quadrant
# For convenience, when drawing, will use ::-1 to make the sequence is E-S-W-N, whose direction is inverse with quadrant
# So the itheta=270 in data is refered to South
########################################################################
# draw the height-azimuth of azimuthal wind 

ir = 50
itheta = 270 # 
delta = 49.383 # m. the space between grids
distance_r = np.multiply((xcen-xminData)/numr*delta, np.arange(0, numr, 1))
dis_r = np.repeat(distance_r, nlev, axis=0).reshape(numr, nlev).T
for iwave in range(0, wavenum):
	fig3 = plt.figure(figsize=(6, 6), dpi=600)
	ax1 = plt.axes([0.1, 0.1, 0.4, 0.6])
	ac1 = ax1.contourf(dis_r, heightP[:, :, itheta], revorVerWavP[iwave, :, :, itheta], \
						cmap=cmaps.WhiteBlueGreenYellowRed, levels=15)

    # Default: 1 color, negatives are dashed# Default: 1 color, negatives are dashed
    # use the loc, lvls, and the linestyles to set different linestyles with different colors 
    # under 0 and above 0
    # if don't use ::-1, the sequence is E-N-W-S, use ::-1, the sequence is E-S-W-N
    # the 0° is in East
	loc = matplotlib.ticker.MaxNLocator(10)
	lvls = loc.tick_values(wVerWavP[iwave, :, :, itheta].min(), wVerWavP[iwave, :, :,  itheta].max())
	ac2 = ax1.contour(dis_r, heightP[:, :, itheta], wVerWavP[iwave, :, :,  itheta], levels=lvls, \
						cmap=BlueOrange10_3.mpl_colormap, linestyles=np.where(lvls >= 0, "-", "--"), \
						linewidths=0.75)
	
	# set the reference line
	ax1.axvline(x=(xcen-xminData)/numr*delta*ir, color='r', linewidth=1)

	# set the extent, as a result that the model height varies with latitude and longitude, 3500m
	ax1.set_ylim(top=3500)

	# set the xy label
	ax1.set_xlabel('Radius (m)', fontsize=8, labelpad=1)
	ax1.set_ylabel('Height(AGL) m', fontsize=8)
	
	# set the xy tick
	ax1.tick_params(axis='x', labelsize=8, pad=2) 
	ax1.tick_params(axis='y', labelsize=8, pad=2) 
	
	# set the attribute of colorbar
	cb = ax1.figure.colorbar(ac1, orientation="horizontal", pad=0.1)
	# cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
	cb.set_label(label='Relative Vorticity')
	cb.ax.tick_params(labelsize=6.5) 
	# fig.colorbar(ac1, orientation='horizontal', pad = 0.1)
	
	# save the figure, pad_inches = 0 and bbox_inches = 'tight' are used to crop the white space
	fig3.savefig('/gpfs3/whuang/Height/Radius/reVor+w/Radius' + Time + '_wave'+str(iwave)+'_reVor+w_in_'+str(itheta)+'.png', \
				pad_inches = 0, bbox_inches = 'tight')
