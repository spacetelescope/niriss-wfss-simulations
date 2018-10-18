#! /usr/bin/env python

from matplotlib import pyplot as plt
import matplotlib.patches as patches
import os
import numpy as np
import astropy
from astropy.io import fits
from astropy.table import QTable
import jwst
from jwst.source_catalog import source_catalog_step
from jwst.datamodels import image, DrizProductModel, WavelengthrangeModel, ModelContainer
from jwst import assign_wcs
from jwst.assign_wcs import assign_wcs_step
from jwst.resample.resample_step import ResampleStep

from pipeline_parameters import * #global variables in all caps

print("Using jwst pipeline version: {}\nastropy version: {}".format(jwst.__version__, astropy.__version__ ))

#-----------------------------------------------------------------------------#

def create_models(infiles=INFILES, plot=True):
    allmodels = []
    allmodels_wcs = []
    if plot is True:
            fig = plt.figure(figsize=(12,10))
    infiles.sort()
    for i,item in enumerate(infiles):
        if plot is True:
            # Plot dithered images
            data = fits.getdata(item)
            ax = fig.add_subplot(2,2,i+1)
            ax.set_adjustable('box-forced')
            #pos = re.search("pos\d", item)
            ax.set_title(os.path.basename(item))
            ax.imshow(data, origin="lower", vmin=-0.5, vmax=0.5)

        # Add WCS to the files
        model = image.ImageModel(item)
        allmodels.append(model)
        awcs = assign_wcs.AssignWcsStep()
        model_wcs = awcs(model)
        allmodels_wcs.append(model_wcs)
    
    if plot is True:
        figname = OUTROOT+"_dithered.png"
        plt.savefig(figname)
    model_container = ModelContainer(allmodels_wcs) 
    model_container.meta.filename = RESAMPLED_MOSAIC
    return model_container

#-----------------------------------------------------------------------------#

def resample(model_container, plot=True):
    resample_step = ResampleStep()
    driz_model = resample_step.call(model_container) #default config file
    print(driz_model.data.shape, driz_model.meta.wcsinfo.roll_ref, driz_model.meta.wcsinfo.ra_ref, driz_model.meta.wcsinfo.dec_ref, driz_model.meta.wcsinfo.crpix1, driz_model.meta.wcsinfo.crpix2)
    if plot is True:
        plt.figure(figsize=(5,5))
        plt.imshow(driz_model.data, origin="lower", vmin=-0.5, vmax=0.5)
        plt.title("Resampled Image")
        figname = OUTROOT+"_resample.png"
        plt.savefig(figname)
    return driz_model

#-----------------------------------------------------------------------------#

def create_source_catalog(driz_model, catname=CATNAME, plot=True):
    sc=source_catalog_step.SourceCatalogStep()
    save_results=sc.call(driz_model, save_results=True, output_file=catname, kernel_fwhm=3, npixels=50, deblend=False, snr_threshold=50)
    catalog = QTable.read(catname, format="ascii.ecsv")
    
    if plot is True:
        plt.figure(figsize=(6,6))
        plt.imshow(driz_model.data, origin="lower", vmin=-0.5, vmax=0.5)
        plt.title("Drizzled Image Catalog")
        for obj in catalog:
            plt.plot(obj["xcentroid"].value, obj["ycentroid"].value, "ro", 
                     markersize=3)
        figname = OUTROOT+"_cat.png"
        plt.savefig(figname)

#-----------------------------------------------------------------------------#

if __name__ == "__main__":
    model_container = create_models()

    driz_model = resample(model_container)
    create_source_catalog(driz_model)
