CAT_COLUMNS = [
        "NUMBER", "X_IMAGE", "Y_IMAGE", "ALPHA_J2000", "DELTA_J2000",
        "MAG_AUTO", "MAGERR_AUTO", "MAG_BEST", "MAGERR_BEST",
        "MAG_APER", "MAGERR_APER", "A_WORLD", "ERRA_WORLD",
        "B_WORLD", "ERRB_WORLD", "THETA_J2000", "ERRTHETA_J2000",
        "ISOAREA_IMAGE", "MU_MAX", "FLUX_RADIUS", "FLAGS"
    ]
#   1 NUMBER                 Running object number                                     
#   2 X_IMAGE                Object position along x                                    [pixel]
#   3 Y_IMAGE                Object position along y                                    [pixel]
#   4 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
#   5 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
#   6 MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]
#   7 MAGERR_AUTO            RMS error for AUTO magnitude                               [mag]
#   8 MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
#   9 MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
#  10 MAG_APER               Fixed aperture magnitude vector                            [mag]
#  11 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
#  12 A_WORLD                Profile RMS along major axis (world units)                 [deg]
#  13 ERRA_WORLD             World RMS position error along major axis                  [deg]
#  14 B_WORLD                Profile RMS along minor axis (world units)                 [deg]
#  15 ERRB_WORLD             World RMS position error along minor axis                  [deg]
#  16 THETA_J2000            Position angle (east of north) (J2000)                     [deg]
#  17 ERRTHETA_J2000         J2000 error ellipse pos. angle (east of north)             [deg]
#  18 ISOAREA_IMAGE          Isophotal area above Analysis threshold                    [pixel**2]
#  19 MU_MAX                 Peak surface brightness above background                   [mag * arcsec**(-2)]
#  20 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
#  21 FLAGS                  Extraction flags                                          

FITS_DIR = "fits"
CAT_DIR = "cats"
LSB_CAT_DIR = "cats_lsb"
GRID_DIR = "gridded_cats"
GRID_END = "_gridded.csv"
GRID_CENTER_DIR = "grid_centers.csv"
GRID_CENTER_END = "_grid_centers.csv"
GRID_COUNT_DIR = "grid_counts"
GRID_COUNT_END = "_grid_counts.csv"
GRID_FLAG_DIR = "grid_flags"
GRID_FLAG_END = "_grid_flags.csv"

FILTER_MATCHED_DIR = "filter_matched"
FILTER_MATCHED_OBJ_DIR = "filter_matched_obj"
FILTER_MATCHED_END = "_filter_matched.csv"
FILTER_MATCHED_OBJ_END = "_obj"
DESI_PRIMARY_DIR = "desi_primaries"
DESI_ID = "targetid"
PRIMARY_INFO_END = "_info.csv"
HEADER_DIR = "headers"
HEADER_END = "_header.csv"

TILE_SIZE = 10000

MASK_DIR = "masks"
GOOD_FRACTION_THRESHOLD = 0.6

# # sdss consts
# GAL_RA = "ra"
# GAL_DEC = "dec"
# SDSS_PRIMARY_DIR = "sdss_primaries"
# SDSS_ID = "objID"
# TILE_COLUMN = "tile_name"

# MAG_ULIM = 25
# MAG_LLIM = 22
# SIZE_LIM = 2.7
# MAG_ABOVE = f"MAG_ABOVE_{MAG_LLIM}"
# MAG_BELOW = f"MAG_BELOW_{MAG_ULIM}"
# SIZE_ABOVE = f"SIZE_ABOVE_{SIZE_LIM}"

# kt17 consts
GAL_RA = "RAdeg"
GAL_DEC = "DEdeg"
KT17_PRIMARY_DIR = "kt17_primaries"
KT17_ID = "# PGC"
TILE_COLUMN = "tile"

MAG_ULIM = 24
MAG_LLIM = 19
SIZE_ULIM = 5000
SIZE_LLIM = 3
MAG_ABOVE = f"MAG_ABOVE_{MAG_LLIM}"
MAG_BELOW = f"MAG_BELOW_{MAG_ULIM}"
SIZE_BETWEEN = f"SIZE_BETWEEN"
SIZE_ABOVE = f"SIZE_ABOVE_{SIZE_LLIM}"
SIZE_BELOW = f"SIZE_BELOW_{SIZE_ULIM}"

# # gaia
# GAL_RA = "star_ra"
# GAL_DEC = "star_dec"
# STAR_ID = "number"


