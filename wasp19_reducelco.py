import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import glob
from photutils import centroid_2dg
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy import wcs
from astropy.table import Table
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy import units as u
import pdb



datadir = '/Users/jspake/projects/inprogress/wasp19/data/lco/all'

def main(datadir=datadir, apsize=7, annsize=[12, 17]):
	# save list of fits files
	os.chdir(datadir)
	flist = glob.glob( datadir + "/*.fits*" )
	flist_file = open("fitsfile_list.txt", "w")
	for item in flist:
		flist_file.write("%s\n" % item)
	flist_file.close()

	# make directory to hold photometry
	phot_dir_name = "/photometry_ap{0}_ann{1}-{2}/".format( apsize, annsize[0], annsize[1] )
	if not os.path.exists( datadir + phot_dir_name ):
		os.makedirs( datadir + phot_dir_name ) 

	# do photometry on first file and get sky coords
	print 'Frame {0} of {1}'.format( 0, len(flist) )
	hdul = fits.open( flist[0] )
	data = hdul[0].data
	header = hdul[0].header
	hdul.close()
	gain = header['GAIN']
	rdnoise = header['RDNOISE']
	error = np.sqrt(data*gain) + rdnoise
	mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5) 
	daofind = DAOStarFinder(fwhm=3.0, threshold=10.*std) 
	sources = daofind(data - median) 
	positions = (sources['xcentroid'], sources['ycentroid'])
	apertures_pix = CircularAperture(positions, r=apsize)
	annulus_pix = CircularAnnulus(positions, r_in=annsize[0], r_out=annsize[1])
	apertures_sky = apertures_pix.to_sky( wcs.WCS(header) )
	annulus_sky = annulus_pix.to_sky( wcs.WCS(header) )
	apers = [ apertures_sky, annulus_sky ]
	norm = ImageNormalize(stretch=SqrtStretch())
	plt.imshow(data)#, cmap='Greys', origin='lower')#, norm=norm)
	apertures_pix.plot(color='white', lw=1, alpha=0.1)
	annulus_pix.plot(color='white', lw=1, alpha=0.1)
	plt.show()
	pdb.set_trace()
	phot_table = aperture_photometry(data, apers, wcs=wcs.WCS(header), error=error)
	bkg_mean = phot_table['aperture_sum_1'] / annulus_pix.area()
	bkg_sum = bkg_mean * apertures_pix.area()
	final_sum = phot_table['aperture_sum_0'] - bkg_sum
	final_error = np.sqrt( phot_table['aperture_sum_err_0']**2 +  phot_table['aperture_sum_err_1']**2 )
	phot_table['residual_aperture_sum'] = final_sum
	phot_table['residual_err'] = final_error
	fname = datadir + phot_dir_name + "phot_table_{0}_{1}s.txt".format(header['MJD-OBS'], header['EXPTIME'])
	phot_table.write(fname, format='ascii')
	replace_comma(fname)

	# pass on sky coords and do photometry on all frames
	for i in range(1,len(flist)):
		print 'Frame {0} of {1}'.format( i, len(flist) )
		do_phot( flist[i], datadir, phot_dir_name, apertures_sky, annulus_sky )



# function to replace comma with space
def replace_comma(fname):
	with open(fname, 'r') as file :
	  filedata = file.read()
	filedata = filedata.replace(',', ' ')
	with open(fname, 'w') as file:
	  file.write(filedata)


# function to do photometry on a given frame
def do_phot(fname, datadir, phot_dir_name, apertures_sky, annulus_sky ):
	apers = [ apertures_sky, annulus_sky ]
	hdul = fits.open( fname )
	if '.fz' in fname:
		data = hdul[1].data
		header = hdul[1].header
		hdul.close()
	else:
		data = hdul[0].data
		header = hdul[0].header
		hdul.close()
	gain = header.get('GAIN')
	if gain == None:
		a = 1
	else:
		rdnoise = header['RDNOISE']
		error = np.sqrt(data*gain) + rdnoise
		mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5) 
		apertures_pix = apertures_sky.to_pixel( wcs.WCS(header) )
		annulus_pix = annulus_sky.to_pixel( wcs.WCS(header) )
		norm = ImageNormalize(stretch=SqrtStretch())
		plt.imshow(data)#, cmap='Greys', origin='lower', norm=norm)
		apertures_pix.plot(color='white', lw=1.5, alpha=0.5)
		annulus_pix.plot(color='white', lw=1.5, alpha=0.5)
		plt.show()
		pdb.set_trace()
		phot_table = aperture_photometry(data, apers, wcs=wcs.WCS(header), error=error)
		bkg_mean = phot_table['aperture_sum_1'] / annulus_pix.area()
		bkg_sum = bkg_mean * apertures_pix.area()
		final_sum = phot_table['aperture_sum_0'] - bkg_sum
		final_error = np.sqrt( phot_table['aperture_sum_err_0']**2 +  phot_table['aperture_sum_err_1']**2 )
		phot_table['residual_aperture_sum'] = final_sum
		phot_table['residual_err'] = final_error
		fname = datadir + phot_dir_name + "phot_table_{0:.7f}_{1}s.txt".format(header['MJD-OBS'], header['EXPTIME'])
		phot_table.write(fname, format='ascii')
		replace_comma(fname)


def convert_phot_files(dirname):
	os.chdir( dirname )
	flist = glob.glob( "*.txt" )
	nframes = len(flist)
	conv_dirname = dirname + '_converted/'
	if not os.path.exists( conv_dirname ):
		os.makedirs( conv_dirname ) 

	for i in range( nframes ):
		mjd = float(flist[i][11:24])
		numb, xcenter, ycenter, ra, dec, aperture_sum_0, \
		aperture_sum_err_0, aperture_sum_1, aperture_sum_err_1, \
		residual_aperture_sum, residual_err = np.genfromtxt( flist[i],\
		 unpack=True, skip_header=1 )
		nstars = len( numb )
		for j in range(nstars):
			line =  mjd, xcenter[j], ycenter[j], ra[j], dec[j], aperture_sum_0[j], \
			aperture_sum_err_0[j], aperture_sum_1[j], aperture_sum_err_1[j], \
			residual_aperture_sum[j], residual_err[j] 
			starfile_name = 'star{0}.txt'.format( j )
			with open( conv_dirname + starfile_name , 'a') as the_file:
				the_file.write( str(line).strip('()') )
				the_file.write( '\n' )

def list_brightness(dirname):
	# make a file that lists the mean brightness over time of each star
	os.chdir( dirname )
	flist = glob.glob( "star*.txt" )
	for i in range(len(flist)):
		ap_sum = np.genfromtxt( flist[i], unpack=True, usecols=5, delimiter=',' )
		ap_median = np.nanmedian( ap_sum )
		with open( "median_brightness.txt" , 'a') as the_file:
			the_file.write( flist[i] + " " + str( ap_median ) )
			the_file.write( '\n' )

def plot_phot(starno):
	# quick plot of photometry, need to be in right directory
	# just give number of star
	fname = "star{0}.txt".format( starno )
	mjd, xcenter, ycenter, ra, dec, aperture_sum_0, \
	aperture_sum_err_0, aperture_sum_1, aperture_sum_err_1, \
	residual_aperture_sum, residual_err = np.genfromtxt( fname, unpack=True,\
	delimiter=',' )
	plt.errorbar( mjd, residual_aperture_sum, yerr=residual_err, fmt='o' )
	plt.show()


# hdul = fits.open(fname)
# data = hdul[0].data
# header = hdul[0].header
# hdul.close()
# gain = header['GAIN']
# rdnoise = header['RDNOISE']
# error = np.sqrt(data*gain) + rdnoise
# mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5) 
 

# daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std) 
# sources = daofind(data - median) 


# positions = (sources['xcentroid'], sources['ycentroid'])
# apertures_pix = CircularAperture(positions, r=7.)
# annulus_pix = CircularAnnulus(positions, r_in=12., r_out=17.)
# apertures_sky = apertures_pix.to_sky( wcs.WCS(header) )
# annulus_sky = annulus_pix.to_sky( wcs.WCS(header) )
# norm = ImageNormalize(stretch=SqrtStretch())
# plt.imshow(data)#, cmap='Greys', origin='lower', norm=norm)
# apertures_pix.plot(color='white', lw=1.5, alpha=0.5)
# annulus_pix.plot(color='white', lw=1.5, alpha=0.5)
# hdul = fits.open(fname2)
# data2 = hdul[0].data
# header2 = hdul[0].header
# hdul.close()
# plt.figure()
# plt.imshow(data2)#, cmap='Greys', origin='lower', norm=norm)
# apertures_pix2 = apertures_sky.to_pixel( wcs.WCS(header2) )
# apertures_pix2.plot(color='white', lw=1.5, alpha=0.5)
# plt.show()
# apers = [ apertures_sky, annulus_sky ]
# phot_table = aperture_photometry(data, apers, wcs=wcs.WCS(header), error=error)
# bkg_mean = phot_table['aperture_sum_1'] / annulus_pix.area()
# bkg_sum = bkg_mean * apertures_pix.area()
# final_sum = phot_table['aperture_sum_0'] - bkg_sum
# final_error = np.sqrt( phot_table['aperture_sum_err_0']**2 +  phot_table['aperture_sum_err_1']**2 )
# phot_table['residual_aperture_sum'] = final_sum
# phot_table['residual_err'] = final_error
# print phot_table

# fname = "table_test.txt"
# phot_table.write(fname, format='ascii')

# # replace comma with space
# with open(fname, 'r') as file :
#   filedata = file.read()
# filedata = filedata.replace(',', ' ')
# with open(fname, 'w') as file:
#   file.write(filedata)



# norm = ImageNormalize(stretch=SqrtStretch())
# plt.imshow(data)#, cmap='Greys', origin='lower', norm=norm)
# apertures_pix.plot(color='white', lw=1.5, alpha=0.5)
# annulus_pix.plot(color='white', lw=1.5, alpha=0.5)


# plt.show()

# Want to first get list of file names
# load fits file
# get time of observation
# find background
# find stars
# convert pixel positions of stars to sky coordinates
# do aperture photometry on sky coord pixels (remove background and get errors too)
# for each data frame save photometry?


