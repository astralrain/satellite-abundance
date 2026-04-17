from libs import *
from consts import *

def read_cat(tile: str) -> pd.DataFrame:
    """
    Read a CFIS tile.cat to exclude the comment header and into a
        pandas dataframe.

    Parameters:
        tile (str): Name of the tile of interest. (ie. CFIS_LSB.xxx.yyy.r)

    Returns:
        pd.DataFrame: Catalogue contents with columns CAT_COLUMNS.
    """
    if not os.path.exists(tile):
        raise FileNotFoundError(tile)
    return pd.read_csv(
        tile,
        sep=r"\s+",
        comment="#",
        names=CAT_COLUMNS,
        header=None,
        engine="python"
    )

def get_lsb_cat(tile: str) -> str:
    """
    Acquires the LSB tile.cat file

    Parameters:
        tile (str): Name of the tile of interest

    Returns:
        str: The path of the LSB .cat file.
    """
    vosclient = Client()

    source = f"vos:cfis/tiles_LSB_DR5/{tile}.cat"

    local_cache = os.path.join(os.getcwd(), LSB_CAT_DIR)
    os.makedirs(local_cache, exist_ok=True)

    local_cat = os.path.join(local_cache, f"{tile}.cat")
    
    if not os.path.exists(local_cat):
        vosclient.copy(source, local_cat)
        print(f"Downloaded {tile}.cat")
    else:
        print(f"{tile}.cat Already in folder")

    return local_cat

def get_fits(tile: str) -> str:
    """
    Acquires the tile.fits file

    Parameters:
        tile (str): Name of the tile of interest

    Returns:
        str: The path of the .fits file.
    """
    vosclient = Client()

    source = f"vos:cfis/tiles_LSB_DR5/{tile}.fits"

    local_cache = os.path.join(os.getcwd(), "fits")
    os.makedirs(local_cache, exist_ok=True)

    local_fit = os.path.join(local_cache, f"{tile}.fits")
    
    if not os.path.exists(local_fit):
        vosclient.copy(source, local_fit)
        print(f"Downloaded {tile}.fits")
    else:
        print(f"{tile}.fits Already in folder")

    return local_fit

def save_header(tile: str) -> str:
    """
    Save the FITS primary header of tile to a CSV file (full header 
        with key, value, comment).

    Parameters:
        tile (str): Name of the tile of interest.

    Returns: 
        str: The path of the header file
    """
    fits_path = get_fits(tile)

    with fits.open(fits_path) as hdul:
        header = hdul[0].header

    header_file = os.path.join(HEADER_DIR, f"{tile}_header.csv")

    with open(header_file, mode="w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["KEY", "VALUE", "COMMENT"])
        for card in header.cards:
            writer.writerow([card.keyword, card.value, card.comment])

    print(f"Header saved to {header_file}")
    
    os.remove(fits_path)
    
    return header_file

def load_header(tile: str):
    """
    Load FITS header from CSV and rebuild a FITS Header object.
    
    Parameters:
        tile (str): Name of the tile of interest.

    Returns:
        headers of tile.fits.
    """
    header_file = Path(HEADER_DIR) / f"{tile}_header.csv" 

    if not header_file.exists():
        save_header(tile)

    header = fits.Header()

    with header_file.open(mode="r", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            key = row["KEY"].strip()
            value = row["VALUE"]
            comment = row["COMMENT"]

            try:
                if value.isdigit():
                    value = int(value)
                else:
                    value = float(value)
            except (ValueError, TypeError):
                pass

            if key:
                header[key] = (value, comment)

    wcs = WCS(header)
    nx = header.get("NAXIS1", None)
    ny = header.get("NAXIS2", None)

    return wcs, nx, ny

def get_hdf5(tile: str):
    """
    Acquires the hdf5 mask file for the tile of interest.

    Parameters:
        tile (str): Name of the tile of interest.

    Returns.
        The tile.mask.h5 file.
    """
    source = Path(f"/arc/projects/NearbyDwarfSearching/masks/{tile}/{tile}.mask.h5")
    file = h5py.File(source, 'r')
    
    return file

def get_mask(tile: str):
    """
    Acquires mask values of tile.mask.h5 in array form.
    
    Parameters:
        tile (str): Name of the tile of interest.

    Returns.
        2D array of the contents of MaskData in the .h5 file.
    """
    source = Path(f"/arc/projects/NearbyDwarfSearching/masks/{tile}/{tile}.mask.h5")
    with h5py.File(source, 'r') as file:
        return file["MaskData"][:]

def load_bad_tiles(bad_tiles_file: str) -> set:
    """
    Load the set of bad tile names from a plain-text file.

    Lines beginning with '#' and blank lines are ignored.  Each remaining
    line is stripped and treated as one tile name (e.g. "CFIS_LSB.123.456.r").

    Parameters:
        bad_tiles_file (str): Path to the bad-tile list.

    Returns:
        set[str]: Set of bad tile name strings.
    """
    bad: set[str] = set()
    with open(bad_tiles_file, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                bad.add(line.split(",")[0].strip())
    print(f"Loaded {len(bad)} bad tiles from {bad_tiles_file}.")
    return bad
 