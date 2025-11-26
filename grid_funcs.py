from libs import *
from consts import *
from basic_funcs import *

def grid_tile(tile: str, cat_df: pd.DataFrame, 
              ngrid: int = 10) -> pd.DataFrame:
    """
    Grids a CFIS tile using pixel coordinates.

    Parameters:
        tile (str): Name of the tile/catalogue of interest.
        cat_df (pd.DataFrame): DataFrame containing the catalogue 
            file of the tile.
        ngrid (int): Number of grid divisions in X and Y, default = 10.

    Returns:
        Output CSV path of gridded catalogue
    """

    x_step = TILE_SIZE / ngrid
    y_step = TILE_SIZE / ngrid

    xs = np.asarray(cat_df["X_IMAGE"], dtype=float)
    ys = np.asarray(cat_df["Y_IMAGE"], dtype=float)

    x_idx = np.floor(xs / x_step).astype(int)
    y_idx = np.floor(ys / y_step).astype(int)

    x_idx = np.clip(x_idx, 0, ngrid - 1)
    y_idx = np.clip(y_idx, 0, ngrid - 1)

    cat_df["i"] = x_idx
    cat_df["j"] = y_idx

    output_file = os.path.join(GRID_DIR, f"{tile}_gridded.csv")
    cat_df.to_csv(output_file, index=False)
    
    print(f"Pixel-based gridded catalogue written to {output_file}")

    return output_file

def compute_grid_centers(tile: str, ngrid: int = 10):
    """
    Compute RA/Dec centers of each grid using pixel coordinates and WCS.
    Grids are defined in pixel space, then converted to sky coordinates.

    Parameters:
        tile (str): Name of the tile/catalogue of interest.
        ngrid (int): Number of grid divisions in X and Y, default = 10.
    
    Returns:
        Output CSV path of the grid centers file.
    """

    wcs, nx, ny = load_header(tile)
    
    print(f"Loaded WCS. Image size = {nx} x {ny} pixels")

    px_step_x = nx / ngrid
    px_step_y = ny / ngrid

    grid_data = []

    for i in range(ngrid):
        for j in range(ngrid):
            x0 = i * px_step_x
            x1 = (i + 1) * px_step_x
            y0 = j * px_step_y
            y1 = (j + 1) * px_step_y

            x_center = x0 + px_step_x / 2.0
            y_center = y0 + px_step_y / 2.0

            ra_center, dec_center = wcs.all_pix2world(x_center, y_center, 0)
            
            coord = SkyCoord(ra=ra_center*u.deg, dec=dec_center*u.deg, frame='icrs')
            gal = coord.galactic

            l_center = gal.l.degree
            b_center = gal.b.degree

            grid_data.append({
                "i": i,
                "j": j,
                "x_min": x0, "x_max": x1,
                "y_min": y0, "y_max": y1,
                "x_center": x_center,
                "y_center": y_center,
                "ra_center": ra_center,
                "dec_center": dec_center,
                "l_center": l_center,
                "b_center": b_center,
            })

    grid_info = pd.DataFrame(grid_data)
    output_file = os.path.join(GRID_CENTER_DIR, tile + "_grid_centers.csv")
    grid_info.to_csv(output_file, index=False)

    print(f"Grid centers written to {output_file}")

    return output_file

def flag_regions(tile: str, flag_file: str) -> pd.DataFrame:
    """
    Given a catalogue and its tile, mark which objects fall in good or 
    bad regions using masks.

    Parameters:
        tile (str): Name of tile of interest.
        gridded_file (str): File of the gridded catalogue.

    Returns:
        The pandas DataFrame of the updated catalogue.
    """
    flag_df = pd.read_csv(flag_file)
    
    if f"GOOD_REGION" in flag_df.columns:
        print(f"The catalogue already contains a 'GOOD_REGION' column.")
        return

    mask = get_mask(tile)
    ny, nx = mask.shape

    x = flag_df["X_IMAGE"].astype(int)
    y = flag_df["Y_IMAGE"].astype(int)

    flag_df["GOOD_REGION"] = False

    valid = (x >= 0) & (x < nx) & (y >= 0) & (y < ny)

    flag_df.loc[valid, "GOOD_REGION"] = mask[y[valid], x[valid]]

    flag_df.to_csv(flag_file, index=False)
    
    print(f"Saved updated flag file with mask info to: {flag_file}")
    
    return flag_df

def flag_magllim(gridded_file: str, flag_file: str) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on magnitude threshold.

    Parameters:
        gridded_file (str): Path to the gridded catalogue CSV file.

    Returns:
        pandas.DataFrame: Updated catalogue with a new boolean column 'MAG_ABOVE_18'.
    """
    grid_df = pd.read_csv(gridded_file)

    if "MAG_AUTO" not in grid_df.columns:
        raise KeyError("The catalogue does not contain a 'MAG_AUTO' column.")

    flag_df = pd.read_csv(flag_file)
    
    if MAG_ABOVE in flag_df.columns:
        print(f"The flag file already contains a {MAG_ABOVE} column.")
        return

    flag_df[MAG_ABOVE] = grid_df["MAG_AUTO"] > MAG_LLIM

    flag_df.to_csv(flag_file, index=False)

    print(f"Saved updated flag file with magnitude flags above {MAG_LLIM} to: {flag_file}")

    return flag_df

def flag_magulim(gridded_file: str, flag_file: str) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on magnitude threshold.

    Parameters:
        gridded_file (str): Path to the gridded catalogue CSV file.

    Returns:
        pandas.DataFrame: Updated catalogue with a new boolean column 'MAG_ABOVE_18'.
    """
    grid_df = pd.read_csv(gridded_file)

    if "MAG_AUTO" not in grid_df.columns:
        raise KeyError("The catalogue does not contain a 'MAG_AUTO' column.")

    flag_df = pd.read_csv(flag_file)
    
    if MAG_BELOW in flag_df.columns:
        print(f"The flag file already contains a {MAG_BELOW} column.")
        return

    flag_df[MAG_BELOW] = grid_df["MAG_AUTO"] < MAG_ULIM

    flag_df.to_csv(flag_file, index=False)

    print(f"Saved updated flag file with magnitude flags below {MAG_ULIM} to: {flag_file}")

    return flag_df

def flag_size(gridded_file: str, flag_file: str) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on magnitude threshold.

    Parameters:
        gridded_file (str): Path to the gridded catalogue CSV file.

    Returns:
        pandas.DataFrame: Updated catalogue with a new boolean column 'MAG_ABOVE_18'.
    """
    grid_df = pd.read_csv(gridded_file)

    if "FLUX_RADIUS" not in grid_df.columns:
        raise KeyError("The catalogue does not contain a 'FLUX_RADIUS' column.")

    flag_df = pd.read_csv(flag_file)
    if SIZE_ABOVE in flag_df.columns:
        print(f"The flag file already contains a {SIZE_ABOVE} column.")
        return

    flag_df[SIZE_ABOVE] = grid_df["FLUX_RADIUS"] > SIZE_LIM

    flag_df.to_csv(flag_file, index=False)

    print(f"Saved updated flag file with flux radius flags above {SIZE_LIM} to: {flag_file}")

    return flag_df

def flagging(tile: str, gridded_file: str, flag_file: str):
    df = pd.read_csv(gridded_file)
    base_cols = ["i", "j", "X_IMAGE", "Y_IMAGE"]
    flags_df = df[base_cols].copy()
    
    flags_df.to_csv(flag_file, index=False)
    
    flag_regions(tile, flag_file)
    flag_magllim(gridded_file, flag_file)
    flag_magulim(gridded_file, flag_file)
    flag_size(gridded_file, flag_file)

def add_good_fraction(tile: str, centers_file: str):
    """
    Adds a 'good_fraction' column to centers_file based on the mask.

    Parameters:
        tile (str): Name of tile of interest.
        centers_file (str): Path to the grid centers file. 
            Must contain ['xmin', 'xmax', 'ymin', 'ymax'] columns.

    Returns:
        pd.DataFrame of the updated grid centres file.
    """
    mask = get_mask(tile)
    ny, nx = mask.shape

    df = pd.read_csv(centers_file)
    good_fractions = []

    for _, row in df.iterrows():
        xmin = max(int(row["x_min"]), 0)
        xmax = min(int(row["x_max"]), nx)
        ymin = max(int(row["y_min"]), 0)
        ymax = min(int(row["y_max"]), ny)

        submask = mask[ymin:ymax, xmin:xmax]
        total_pixels = submask.size
        good_pixels = np.sum(submask)
        if total_pixels > 0:
            fraction = good_pixels / total_pixels
        else:
            fraction = np.nan

        good_fractions.append(fraction)

    df["good_fraction"] = good_fractions

    df.to_csv(centers_file, index=False)

    print(f"Added 'good_fraction' column and saved to: {centers_file}")
    return df

def grid_counts(tile: str, flag_file: str, centers_file: str) -> None:
    """
    Add normalized object counts per grid cell into centers_file,
    using the GOOD_REGION flag from gridded_file and correcting 
    for coverage fraction.

    Parameters:
        gridded_file (str): Path to the input catalogue CSV.
            Must have columns 'i', 'j', 'X_IMAGE', 'Y_IMAGE', 'GOOD_REGION'.
        centers_file (str): Path to the existing grid info CSV.
            Must have columns 'i', 'j'.

    Returns:
        None (writes to CSV).
    """
    flag_df = pd.read_csv(flag_file)
    grid_info = pd.read_csv(centers_file)

    required_columns = ["i", "j", "X_IMAGE", "Y_IMAGE", 
                        "GOOD_REGION", MAG_ABOVE, 
                        MAG_BELOW, SIZE_ABOVE,]
    
    for col in required_columns:
        if col not in flag_df.columns:
            raise ValueError(f"Input catalogue must have column '{col}'")

    selected = flag_df[
        (flag_df["GOOD_REGION"]) &
        (flag_df[MAG_ABOVE]) &
        (flag_df[MAG_BELOW]) &
        (flag_df[SIZE_ABOVE])
    ]

    good_counts = (
        selected
        .groupby(["i", "j"])
        .size()
        .reset_index(name="good_count")
        )

    merged = pd.merge(
        grid_info[["i", "j", "x_center", "y_center", 
                   "ra_center", "dec_center", "l_center", "b_center", "good_fraction"]], 
        good_counts, on=["i", "j"], 
        how="left"
    )
    merged["good_count"] = merged["good_count"].fillna(0)

    merged["normalized_count"] = merged["good_count"] / merged["good_fraction"]
    merged.loc[~np.isfinite(merged["normalized_count"]), 
               "normalized_count"] = np.nan
    
    output_file = os.path.join(GRID_COUNT_DIR, f"{tile}_grid_counts.csv")
    merged.to_csv(output_file, index=False)
    print(f"Updated grid counts written to {centers_file}")
