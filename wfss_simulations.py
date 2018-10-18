#! /usr/bin/env python

"""
Simulate direct and dispersed images for the NIRISS/WFSS detector.
All input parameters and globals are stored in wfss_parameters.py
and can be identified as all-caps variables. 
Make sure that rotate_CD_matrix is its own module in the same
directory as this script. 

Usage:
    > python wfss_simulations.py
    
Important parameters to change in wfss_parameters.py are:
    - RA_CENTER & DEC_CENTER (should be selected to include targets of 
      choice in the source catalog)
    - OUTROOT 
    - SOURCE_CATALOG
    - DITHER_X & DITHER_Y
    - FLAT_SPECTRUM
    - SEDS_FILE 
    
To create a new source catalog file with new sources, use the function
pick_sources from ipython or a new script. Sample usage:
    > from wfss_simulations import pick_sources
    > c = pick_sources(nsources=200)
This will create a new catalog file with 200 randomly selected sources.
To change selection criteria, modify wfss_parameters.py parameters:
    - IND_MAX
Additional selection criteria can be added as desired to pick_sources() 
"""

import os
os.environ["WEBBPSF_PATH"] = "/grp/jwst/ote/webbpsf-data"
os.environ["PYSYN_CDBS"] = "/grp/hst/cdbs/"
from collections import OrderedDict
import copy
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt 
#plt.style.use("mystandard")
from astropy.table import Table
from astropy.io import fits
import astropy.wcs as pywcs
import astropy.units as u
from astropy.modeling import models
from astropy.io import ascii
import stsci.convolve
import webbpsf
import reproject
import glob
import grizli
import grizli.multifit
import grizli.jwst
from jwst.resample import gwcs_blot
import asdf
from asdf import fits_embed

from rotate_CD_matrix import rotate_CD_matrix
from wfss_parameters import * # all these imports are global variables (capitals)

#-----------------------------------------------------------------------------#

def pick_sources(catalog=GOODS_CATALOG, nsources=1, only_gal=True, write=True):
    """
    Randomly pick sources from the GOODS-S catalog from the wfss-simulation
    grit repo: https://grit.stsci.edu/NIRISS/wfss-simulation

    Selection criteria are stored in wfss_parameters.py.
    
    Args:
        catalog (str): The (ascii) catalog that contains sources to choose 
            from.
        nsources (int): Number of sources to choose from catalog.
        only_gal (Bool): If True, only galaxies will be selected.
        write (Bool): True if selected sources should be written to new
            catalog file. If so, name will be "{}sources.cat".format(nsources)
    
    Returns:
        out_cat (Table): astropy Table object with randomly selected sources.  
    """
    
    # Include a buffer of the dither size in case large dithers are used.
    max_x_dither = max(DITHER_X)
    max_y_dither = max(DITHER_Y)
    max_ra_dither = (max_x_dither*0.0655)/3600. # Approx. plate scale
    max_dec_dither = (max_x_dither*0.0655)/3600. # Approx. plate scale

    ralo = RA_CENTER - (RA_SPAN/2.) - max_ra_dither
    rahi = RA_CENTER + (RA_SPAN/2.) + max_ra_dither
    declo = DEC_CENTER - (DEC_SPAN/2.) - max_dec_dither
    dechi = DEC_CENTER + (DEC_SPAN/2.) + max_dec_dither

    rows = []
    t = Table.read(catalog, format="ascii.commented_header")
    for i in range(len(t)):
        if ralo < t[i]["ra"] < rahi:
            if declo < t[i]["dec"] < dechi:
                if t[i]["n"] <= IND_MAX:
                    rows.append(i)

    if nsources > len(rows):
        print("Number of sources desired, {}, exceeds number "
              "of sources that match selection criteria. "
              "Using {} instead".format(nsources, len(rows)))
        nsources = len(rows)
    source_rows = np.random.choice(rows, nsources)

    out_cat = Table(rows=[t[x] for x in source_rows], names=t.columns)
    if write is True:
        catname = "{}sources.cat".format(nsources)
        out_cat.write(catname, 
                      format="ascii.commented_header")
        print("Wrote catalog {}".format(catname))
    
    return out_cat

#-----------------------------------------------------------------------------#

def update_catalog(catalog=SOURCE_CATALOG):
    """
    Add MAG_AUTO column to input source catalog.

    Args:
        catalog (str): Name of source catalog ascii table.
    """
    
    cat = Table.read(catalog, format="ascii.commented_header")
    if "MAG_AUTO" not in cat.columns:
        cat["MAG_AUTO"] = cat[MAG_COL]
        Table.write(cat, SOURCE_CATALOG, format="ascii.commented_header", overwrite=True)
        print("Updated catalog {}".format(catalog))

#-----------------------------------------------------------------------------#

def read_rate_im(in_rate=RATE_IMAGE):
    """
    Read the rate image from NIRISS day-in-the-life testing.

    Returns:
        rate_im_orig (HDUList): The rate image file in astropy.io.fits 
            HDUList format.
        aa (AsdfInFits): The rate image file in asdf.fits_embed AsdfInFits
            format.  
    """
    
    rate_im_orig = fits.open(in_rate)
    rate_im_orig[0].header['PUPIL'] = NIS_FILTER.upper()
    
    aa = asdf.open(in_rate)
    
    #Add photometry information to header
    #Possibly ABZP is calculated in model.py - need to double check
    rate_im_orig[0].header['ABZP'] = (ZPs[NIS_FILTER.lower()], 'Filter AB zeropoint')
    rate_im_orig[0].header['PHOTPLAM'] = (PIVOT[NIS_FILTER.lower()], 'Filter pivot wavelength')
    log_photflam = rate_im_orig[0].header['ABZP'] - 18.6921 + 21.1 + 5*np.log10(rate_im_orig[0].header['PHOTPLAM'])
    rate_im_orig[0].header['PHOTFLAM'] = (10**(-0.4*log_photflam), 'Transformation to flambda cgs flux density')
    Jy = 10**(-0.4*(ZPs[NIS_FILTER.lower()]-23.9))*1.e-6
    rate_im_orig[0].header['PHOTFNU'] = (Jy, 'Transformation to Jansky')
    fnu_scale = Jy*u.Jy/(0.0656*u.arcsec)**2
    rate_im_orig[0].header['DN2MJYSR'] = (fnu_scale.to(u.MJy/u.sr).value, 'Transformation to MJy/steradian')

    #Adding PA_V3 just in case this is the reason for crushing
    rate_im_orig[0].header['PA_V3'] = 0
    rate_im_orig['SCI'].header['XPA_V3'] = 0

    return rate_im_orig, aa

#-----------------------------------------------------------------------------#

def make_niriss_image(plot=False):
    """
    Create a direct image with sources that have been convolved with the 
    NIRISS PSF from webbpsf as well as the segmentation (ID) image.
    
    Args:
        plot (Bool): True if plots should be saved.
    """

    # Add only ra-dec, rotation, and wcs to add the sources at the right places
    # `size` is image dimension in arcsec
    hdu = grizli.utils.make_wcsheader(size=NAXIS[0]*0.0656, pixscale=0.0656, get_hdu=True) 
    nis_header = hdu.header
    nis_header["CRVAL1"] = RA_CENTER
    nis_header["CRVAL2"] = DEC_CENTER
    
    # Rotate image to desired PA 
    nis_wcs = pywcs.WCS(nis_header)
    cd_rot = rotate_CD_matrix(nis_wcs.wcs.cd, PA_APER)
    for i in range(2):
        for j in range(2):
            nis_header["CD{0}_{1}".format(i+1,j+1)] = cd_rot[i][j]
    nis_wcs = pywcs.WCS(nis_header)    
    nis_wcs.pscale = 0.0656 #grizli.utils.get_wcs_pscale(nis_wcs)
    nis_header["PA_APER"] = PA_APER
    
    # Load the catalog and compute detector flux
    gfit = Table.read(SOURCE_CATALOG, format='ascii.commented_header')
    object_mag = gfit[MAG_COL]
    ZP = ZPs[NIS_FILTER.lower()]
    object_flux = 10**(-0.4*(object_mag-ZP))
    
    # Determine which objects are within the NIRISS FoV
    xc, yc = nis_wcs.all_world2pix(gfit['ra'], gfit['dec'], 0)
    obj_in_img = (xc > 1) & (yc > 1) & (xc < NAXIS[0]-1) & (yc < NAXIS[1]-1) 
    obj_in_img &= object_flux > 0
    
    # Point sources and faint sources
    stars = obj_in_img & ((gfit['star_flag'] == 1) | (object_mag > MAX_MAG))
    star_idx = np.arange(len(gfit))[stars]
    
    # Galaxies
    if 're' not in gfit.colnames:
        gals = obj_in_img < -100 # False
    else:
        gals = obj_in_img & (~stars) & (object_mag <= MAX_MAG) & (gfit['re'] > 0)
    gal_idx = np.arange(len(gfit))[gals]
    
    # Put sources in the image
    # Initialize model and segmentation images with zeros
    nis_model = np.zeros(NAXIS[::-1], dtype=np.float32)
    nis_seg = np.zeros(NAXIS[::-1], dtype=int)
        
    # pixel indices
    yp, xp = np.indices(NAXIS[::-1])    
        
    # Add star point sources
    xpix = np.cast[int](np.round(xc))
    ypix = np.cast[int](np.round(yc))
    for ix in star_idx:
        nis_model[ypix[ix], xpix[ix]] = object_flux[ix]
        
    if SEG_THRESHOLD > 0:
        Rseg = 5
        for i, ix in enumerate(star_idx):
            print('seg: {0} ({1}/{2})'.format(ix, i, stars.sum()))
            R = np.sqrt((xp-xpix[ix])**2+(yp-ypix[ix])**2)
            clip_seg = (R <= Rseg) & (nis_seg == 0)
            nis_seg[clip_seg] = gfit['id'][ix]
            
    # Add Sersic sources
    for i, ix in enumerate(gal_idx):
        print('ID: {0} ({1}/{2})'.format(ix, i, gals.sum()))
            
        # Effective radius, in pixels
        re = gfit['re'][ix]/nis_wcs.pscale
        se = models.Sersic2D(amplitude=1., r_eff=re, n=gfit['n'][ix],
                             x_0=xc[ix], y_0=yc[ix], ellip=1-gfit['q'][ix],
                             theta=(gfit['pa'][ix]+90-PA_APER)/180*np.pi)
            
        # Normalize to catalog flux
        m = se(xp, yp)
        renorm = object_flux[ix]/m.sum()
        if not np.isfinite(renorm):
            continue
                
        # Add to model image
        nis_model += m*renorm
        if SEG_THRESHOLD > 0:
            # Add to segmentation image
            clip_seg = (m*renorm > SEG_THRESHOLD) & (nis_seg == 0)
            nis_seg[clip_seg] = gfit['id'][ix]
    
    # Check image 
    if plot is True:
        plt.figure(figsize=(10,10))
        plt.subplot(projection=pywcs.WCS(nis_header))
        plt.imshow(np.log10(nis_model))
        #plt.imshow(nis_seg)
        plt.grid(color='black', ls='solid')
        plt.xlabel('Galactic Longitude')
        plt.ylabel('Galactic Latitude')
        plt.title("Sources (logscale)")
        figname = "{}_sources.png".format(OUTROOT)
        plt.savefig(figname)
        print("Wrote {}".format(figname))
        #plt.show()
    
    # Convolve with NIRISS PSF and save image and segmentation map.
    nis = webbpsf.NIRISS()
    # Add -0.5, -0.5 pixel offset to get convolved image in correct place ??
    nis.options['source_offset_r'] = 0.5*nis_wcs.pscale # offset in arcseconds
    nis.options['source_offset_theta'] = 135.   # degrees CCW from +Y
    
    # Get the PSF
    nis.filter = NIS_FILTER.upper()
    psf_hdu = nis.calcPSF(fov_pixels=64) # power of 2 for fast FFT convolution  
    psf = psf_hdu[1].data
    psf /= psf.sum()
    
    # Convolve the model image
    nis_fullsim = stsci.convolve.convolve2d(nis_model, psf, output=None, mode='nearest', cval=0.0, fft=1)
    
    # Mask low fluxes to make images gzip smaller
    #mask = nis_fullsim > CLIP_FLUX
    #nis_fullsim[~mask] = 0
        
    # Save the output image
    filename = '{0}-{1}.fits'.format(OUTROOT, nis.filter.lower())
    fits.writeto(filename, 
                   data=nis_fullsim, header=nis_header,
                   overwrite=True, output_verify='fix')
    print("Wrote {}".format(filename))

    # Save the segmentation image
    if SEG_THRESHOLD > 0:
        filename = '{0}-{1}_seg.fits'.format(OUTROOT, nis.filter.lower())
        fits.writeto(filename,
                       data=nis_seg, header=nis_header,
                       overwrite=True, output_verify='fix')
        print("Wrote {}".format(filename))

    # Check image. 
    if plot is True: 
        plt.subplot(projection=pywcs.WCS(nis_header))
#        plt.imshow(nis_fullsim)
        plt.imshow(np.log10(nis_fullsim))
        #plt.imshow(nis_seg)
        plt.grid(color='black', ls='solid')
        plt.xlabel('Galactic Longitude')
        plt.ylabel('Galactic Latitude')
        plt.title("Distorted sources convolved with PSF (logscale)")
        figname = "{}_sourcepsf.png".format(OUTROOT)
        plt.savefig(figname)
        print("Wrote {}".format(figname))
        #plt.show()

#-----------------------------------------------------------------------------#

def perform_dither_distortion(plot=False):
    """
    Dither the initial direct image and apply detector distortion.

    Args:
        plot (Bool): True if plots should be saved.
    """

    # Make four dithered exposures following the standard template pattern 
    # or whatever specified pattern. It reopens image just created and it 
    # uses a reference rate image to create the right and complete header.

    # Dithers for WFSS (use the sim image just created)
    ref = fits.open('{0}-{1}.fits'.format(OUTROOT, NIS_FILTER.lower()))
    ref_model = grizli.jwst.hdu_to_imagemodel(ref[0])
    blot_parent = gwcs_blot.GWCSBlot(ref_model)
    # Center position on the image. 
    r0 = RA_CENTER * u.deg
    d0 = DEC_CENTER * u.deg

    # Dither pattern
    dx = DITHER_X
    dy = DITHER_Y
    
    if plot is True:
        plt.figure(figsize=(20,20))
    
    #Make the 4 dithers with distortion
    for pos in range(4):
        rate_im, aa = read_rate_im()
        aa['data'] = rate_im['SCI'].data
        aa['dq'] = rate_im['DQ'].data
        wcsinfo = aa['meta']['wcsinfo']
        pointing = aa['meta']['pointing']
        inst = aa['meta']['instrument']
        print('Position {0}, {1}'.format(pos+1, NIS_FILTER))
        inst['filter'] = rate_im[0].header['FILTER'] = 'CLEAR'
        inst['pupil'] = rate_im[0].header['PUPIL'] = NIS_FILTER.upper()
        
        cosd = np.cos(d0/180*np.pi)
        pos_ra = (r0+dx[pos]/cosd*u.arcsec).value
        pos_dec = (d0+dy[pos]*u.arcsec).value
            
        #Update header in rate image and wcf info
        rate_im[1].header['CRVAL1'] = rate_im[1].header['RA_REF'] = pos_ra
        rate_im[1].header['CRVAL2'] = rate_im[1].header['DEC_REF'] = pos_dec
        wcsinfo['crval1'] = wcsinfo['ra_ref'] = pointing['ra_v1'] = pos_ra
        wcsinfo['crval2'] = wcsinfo['dec_ref'] = pointing['dec_v1'] = pos_dec
    
        # Blot reference data to distorted frame
        rate_model = grizli.jwst.img_with_wcs(rate_im)
        blot_child = blot_parent.extract_image(rate_model, interp='poly5')       
        rate_im['SCI'].data = blot_child
    
        # zero-out DQ array
        rate_im['DQ'].data *= 0
        
#        rate_im[0].header['TARGNAME'] = aa['meta']['target'] = OUTROOT
    
        print('{0}-{1}-clear-pos{2}_rate.fits'.format(OUTROOT, NIS_FILTER, pos+1))
        ff = fits_embed.AsdfInFits(rate_im, {'meta':aa['meta']})
        ff.write_to('{0}-{1}-clear-pos{2}_rate.fits'.format(OUTROOT, NIS_FILTER, pos+1), overwrite=True)
                
        # Make grism files (sci extension is direct image for now)
        for grism in ['gr150r', 'gr150c']:
            inst['filter'] = rate_im[0].header['FILTER'] = grism.upper()
            out = '{0}-{1}-{grism}-pos{2}_rate.fits'.format(OUTROOT, NIS_FILTER, pos+1, grism=grism)
            print(out)
            ff = fits_embed.AsdfInFits(rate_im, aa['meta'])
            ff.write_to(out, overwrite=True)
    
        rate_im.close() 
        
        # Show dithered images
        if plot is True:
            plt.subplot(2, 2, pos+1, projection=pywcs.WCS(rate_im[1].header))
            plt.title("Dither=({}, {})".format(dx[pos], dy[pos]))
            plt.imshow(np.log10(rate_im['SCI'].data))
            plt.grid(color='black', ls='solid')
    
    figname = "{}_dithers.png".format(OUTROOT)
    plt.savefig(figname)
    print("Wrote {}".format(figname))
        #plt.colorbar()
        #plt.xlabel('Galactic Longitude')
        #plt.ylabel('Galactic Latitude')
        #plt.show()

#-----------------------------------------------------------------------------#

def disperse_image(plot=False):
    """
    Disperse the dithered, distorted images.

    Args:
        plot (Bool): True if plots should be saved.
    """
    # Compute the grism model and fill the grism exposures with these models.
    # Uses the images just created (seg map, direct image, and the placeholder for the grism).
    # Uses the catalog file.
    
    # It can use:
    # - flat spectra
    # - spectra generated with EAZY (have not tested this yet)
    # - arbitrary spectra (still need to code this one)


    # Use the image just created
    grism_files = glob.glob('{0}-{1}-gr15*pos*_rate.fits'.format(OUTROOT,NIS_FILTER.lower()))
    direct_files = [x.replace("gr150c", "clear") if "gr150c" in x else x.replace("gr150r", "clear") for x in grism_files]*2
    g_seg_file = '{0}-{1}_seg.fits'.format(OUTROOT, NIS_FILTER.lower())

    grp = grizli.multifit.GroupFLT(grism_files=grism_files,
                                       direct_files=direct_files,
                                       ref_file=None,
                                       seg_file=g_seg_file,
                                       catalog=SOURCE_CATALOG, cpu_count=4)
    
    #not using EAZY to get the spectra, just default flat spectra (if it works)
    if FLAT_SPECTRUM is True:
        # Flat f-lambda spectra
        fit_info = None
    else:
        seds = pd.read_table(SEDS_FILE, sep="\s+")
        seds_clip = seds[(seds['wave']>5000.) & (seds['wave']<35000.)]
        sedwave = seds_clip['wave'].values
        sedflux = seds_clip['f_lam'].values
        fit_info = OrderedDict()
        #only one galaxy, so ID=1 (has to match with the ID in the catalog file)
        sources = Table.read(SOURCE_CATALOG, format='ascii.commented_header')
        ids = list(sources["id"])
        fit_info = {ids[x]: {"mag": -1, "spec": [sedwave, sedflux]} for x in range(len(ids))}
    
    # With real spectra
    #fit_info = generate_fit_info(filter=NIS_FILTER)
    
    # Use flat spectra for everything
    print('Compute all objects to AB<chosen mag limit, flat spectrum')
    grp.compute_full_model(fit_info=fit_info, verbose=False, store=False, mag_limit=MAX_MAG, coeffs=[1], cpu_count=0)
    # Add galaxy SEDs for brighter objects
    #print('Compute bright galaxies (N={0})'.format(len(fit_info)))
    #grp.compute_full_model(fit_info=fit_info, verbose=True, store=False, mag_limit=28, coeffs=[1], cpu_count=0)

    # Put the models into the grism FITS exposure files
    for i in range(len(grism_files)):
        if grp.FLTs[i].is_rotated:
            grp.FLTs[i].transform_NIRISS(verbose=False)
                
        im = fits.open(grism_files[i], mode='update')
        pad = grp.FLTs[i].pad
        im['SCI'].data = grp.FLTs[i].model[pad:-pad, pad:-pad]
        im.flush()
        
        # Rotate back
        grp.FLTs[i].transform_NIRISS(verbose=False)
        print("Updated {}".format(grism_files[i]))

    if plot is True:
        maxm = None
        for final_im in grism_files:
            data = fits.getdata(final_im)
            if maxm is None:
                maxm = max(data.flatten())
                vmax = maxm/4.
            plt.figure(figsize=(10,10))
            plt.imshow(data, vmin=-0.1, vmax=vmax)
            plt.colorbar()
            plt.title(final_im)
            plt.xlabel('X [detector]')
            plt.ylabel('Y [detector]')
#            plt.xlim(800,1100)
#            plt.ylim(800,1100)
            figname = "{}_dispersed.png".format(final_im.strip(".fits"))
            plt.savefig(figname)
            print("Saved {}".format(figname))

#-----------------------------------------------------------------------------#

if __name__ == "__main__":
    plot = True
    
    update_catalog()
    make_niriss_image(plot)
    perform_dither_distortion(plot)
    disperse_image(plot) 
