from collections import OrderedDict

# Ab image zeropoints computed to match the ETC as of Feb. 21.
ZPs = OrderedDict([("f090w", 28.27100318778076),
                   ("f115w", 28.329561852474136),
                   ("f150w", 28.183657683252214),
                   ("f200w", 28.260233956567966),
                   ("f140m", 27.39442027680035),
                   ("f158m", 27.458684529328927)])
PIVOT = OrderedDict([("f090w", 9009.324775150702),
                     ("f115w", 11494.953696000099),
                     ("f150w", 14929.433412860893),
                     ("f200w", 19925.71415412438),
                     ("f140m", 14035.440274803774),
                     ("f158m", 15857.162918853444)])

# Define the center of the image.
#RA_CENTER = 53.090861 # ORIG
#DEC_CENTER = -27.956420 # ORIG
RA_CENTER = 53.15 # Center of GOODS-S master catalog (GOODS_CATALOG)
DEC_CENTER = -27.8 # Center of GOODS-S master catalog (GOODS_CATALOG)

# The span of the NIRISS detector in degrees (2.2').
RA_SPAN = 0.0366667 # 2.2' = 0.036667deg
DEC_SPAN = 0.0366667

# Max Sersic index to use when picking targets from GOODS-S.
IND_MAX = .8 
# Catalog to use when choosing GOODS-S targets
GOODS_CATALOG = "/user/jotaylor/NIRISS/wfss_simulation/inputs/goodss_3dhst.v4.1.niriss_galfit.cat"

# Position Angle of orientation on the sky.
PA_APER = 0
NAXIS = (2448, 2448)
NIS_FILTER = "f150w"
MAG_COL = "mag_f150w"
MAX_MAG = 26
SEG_THRESHOLD = 0.02
CLIP_FLUX = 1.e-5

# Rootname of output files.
# Catalog to use that has source information
OUTROOT = "500sources_realspec"
SOURCE_CATALOG = "/user/jotaylor/NIRISS/wfss_simulation/inputs/500sources.cat"
#OUTROOT = "1source" # ORIG
#SOURCE_CATALOG = "/user/jotaylor/NIRISS/wfss_simulation/inputs/single_source.cat" #ORIG

# I don't think this is used vv
TEMPLATE_IMAGE = "/user/jotaylor/NIRISS/wfss_simulation/inputs/jw87600018001_02101_00007_nis_rate.fits"
# Rate image from NIRISS day in the life testing
RATE_IMAGE = "/user/jotaylor/NIRISS/wfss_simulation/inputs/jw00306001001_02101_00001_nis_rate.fits"

# Dither pattern 
# Very large 4-point dither pattern
DITHER_X = [0, 10, 10, -10] #in arcsec [0, 0.3270, 0.1635, -0.1635]
DITHER_Y = [0, 10, -10, -10] #in arcsec [0, 0.1645, 0.3948, 0.2961]

# Spectral model for dispersed objects.
FLAT_SPECTRUM = False
SEDS_FILE = "/user/jotaylor/NIRISS/wfss_simulation/inputs/single_source_big.spec"
