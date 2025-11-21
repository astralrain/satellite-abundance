CAT_COLUMNS = [
        "NUMBER", "X_IMAGE", "Y_IMAGE", "ALPHA_J2000", "DELTA_J2000",
        "MAG_AUTO", "MAGERR_AUTO", "MAG_BEST", "MAGERR_BEST",
        "MAG_APER", "MAGERR_APER", "A_WORLD", "ERRA_WORLD",
        "B_WORLD", "ERRB_WORLD", "THETA_J2000", "ERRTHETA_J2000",
        "ISOAREA_IMAGE", "MU_MAX", "FLUX_RADIUS", "FLAGS"
    ]
CAT_DIR = "cats"
GRID_DIR = "gridded_cats"
GRID_END = "_gridded"
GRID_CENTER_DIR = "grid_centers"
GRID_CENTER_END = "_grid_centers"
GRID_COUNT_DIR = "grid_counts"
GRID_COUNT_END = "_grid_counts"
GRID_FLAG_DIR = "grid_flags"
GRID_FLAG_END = "_grid_flags"
FILTER_MATCHED_DIR = "filter_matched"
FILTER_MATCHED_OBJ_DIR = "filter_matched_obj"
FILTER_MATCHED_END = "_filter_matched"
FILTER_MATCHED_OBJ_END = "_obj"
DESI_PRIMARY_DIR = "desi_primaries"
DESI_ID = "targetid"
SDSS_PRIMARY_DIR = "sdss_primaries"
SDSS_ID = "objID"
PRIMARY_INFO_END = "_info"
HEADER_DIR = "headers"
HEADER_END = "_header"

TILE_SIZE = 10000

MASK_DIR = "masks"
MAG_ULIM = 24.5
MAG_LLIM = 18
SIZE_LIM = 3
MAG_ABOVE = f"MAG_ABOVE_{MAG_LLIM}"
MAG_BELOW = f"MAG_BELOW_{MAG_ULIM}"
SIZE_ABOVE = f"SIZE_ABOVE_{SIZE_LIM}"