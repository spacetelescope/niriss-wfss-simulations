import glob

#OUTROOT = "niriss-sim-sept18"
#INFILES = glob.glob("/user/jotaylor/NIRISS/wfss_simulation/box/"+OUTROOT+"-f150w-clear*rate.fits")
OUTROOT = "niriss-sim-sept18-stars-2048"
INFILES = glob.glob("/user/jotaylor/NIRISS/wfss_simulation/box/images_2048/"+OUTROOT+"-f150w-clear*rate.fits")
RESAMPLED_MOSAIC = OUTROOT+"-resample.fits"

CATNAME = "sept18_test_cat.ecsv"
