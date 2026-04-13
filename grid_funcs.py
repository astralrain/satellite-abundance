from libs import *
from consts import *
from basic_funcs import *

def grid_tile(tile: str, cat_df: pd.DataFrame, ngrid: int = 10) -> str:
    """
    Grids a CFIS tile using pixel coordinates with its catalogue data.

    Parameters:
        tile (str): Name of the tile/catalogue of interest.
        cat_df (pd.DataFrame): DataFrame containing the catalogue 
            file of the tile of interest.
        ngrid (int): Number of grid divisions in X and Y, 
            default is 10.

    Returns:
        str: Path of gridded catalogue CSV file
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

    output_file = Path(GRID_DIR) / (tile + GRID_END)
    cat_df.to_csv(output_file, index=False)
    
    print(f"Pixel-based gridded catalogue written to {output_file}")

    return output_file

def compute_grid_centers(tile: str, ngrid: int = 10) -> str:
    """
    Compute the centers of each grid using pixel coordinates and WCS. 
        Grids are defined in pixel space, then converted to sky 
        coordinates, RA/Dec and galactic.

    Parameters:
        tile (str): Name of the tile/catalogue of interest.
        ngrid (int): Number of grid divisions in X and Y, default is 10.
    
    Returns:
        str: Path of grid centers CSV file
    """

    wcs, nx, ny = load_header(tile)
    
    print(f"Loaded WCS. Image size = {nx} x {ny} pixels")

    x_step = nx / ngrid
    y_step = ny / ngrid

    grid_data = []

    for i in range(ngrid):
        for j in range(ngrid):
            x0 = i * x_step
            x1 = (i + 1) * x_step
            y0 = j * y_step
            y1 = (j + 1) * y_step

            # pixel centers
            x_center = x0 + x_step / 2.0
            y_center = y0 + y_step / 2.0

            # ra/dec centers
            ra_center, dec_center = wcs.all_pix2world(x_center, y_center, 0)

            # galactic centers
            coord = SkyCoord(
                ra=ra_center*u.deg, 
                dec=dec_center*u.deg, 
                frame='icrs'
                )
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
    output_file = Path(GRID_CENTER_DIR) / (tile + GRID_CENTER_END)
    grid_info.to_csv(output_file, index=False)

    print(f"Grid centers written to {output_file}")

    return output_file

def flag_regions(tile: str, flag_file: Path) -> pd.DataFrame:
    """
    Given a catalogue and its tile, mark which objects fall in good or 
    bad regions using masks where True indicates a good region.

    Parameters:
        tile (str): Name of tile of interest.
        flag_file (Path): File path of the flag file.

    Returns:
        pandas.DataFrame: Updated flag file with a new boolean column.
    """
    flag_df = pd.read_csv(flag_file)
    
    if "GOOD_REGION" in flag_df.columns:
        print(f"The flag file already contains a 'GOOD_REGION' column.")
        return flag_df

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

def flag_magllim(gridded_df: pd.DataFrame, flag_file: Path) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on lower magnitude 
        threshold.

    Parameters:
        gridded_df (pd.DataFrame): Datafram of the gridded catalogue file.
        flag_file (Path): Path to the flag file.

    Returns:
        pandas.DataFrame: Updated flag file with a new boolean column.
    """
    if "MAG_AUTO" not in gridded_df.columns:
        raise KeyError("The catalogue does not contain a 'MAG_AUTO' column.")

    flag_df = pd.read_csv(flag_file)

    if len(flag_df) != len(gridded_df):
        raise ValueError("Flag file and gridded file have different lengths.")
    
    if MAG_ABOVE in flag_df.columns:
        print(f"The flag file already contains a {MAG_ABOVE} column.")
        return flag_df

    flag_df[MAG_ABOVE] = gridded_df["MAG_AUTO"] > MAG_LLIM

    flag_df.to_csv(flag_file, index=False)

    print(f"Saved updated flag file with magnitude flags above {MAG_LLIM} to: {flag_file}")

    return flag_df

def flag_magulim(gridded_df: pd.DataFrame, flag_file: Path) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on upper magnitude 
        threshold.

    Parameters:
        gridded_df (pd.DataFrame): Datafram of the gridded catalogue file.
        flag_file (Path): Path to the flag file.

    Returns:
        pandas.DataFrame: Updated flag file with a new boolean column.
    """
    if "MAG_AUTO" not in gridded_df.columns:
        raise KeyError("The catalogue does not contain a 'MAG_AUTO' column.")

    flag_df = pd.read_csv(flag_file)

    if len(flag_df) != len(gridded_df):
        raise ValueError("Flag file and gridded file have different lengths.")
    
    if MAG_BELOW in flag_df.columns:
        print(f"The flag file already contains a {MAG_BELOW} column.")
        return flag_df

    flag_df[MAG_BELOW] = gridded_df["MAG_AUTO"] < MAG_ULIM

    flag_df.to_csv(flag_file, index=False)

    print(f"Saved updated flag file with magnitude flags below {MAG_ULIM} to: {flag_file}")

    return flag_df

def flag_sizellim(gridded_df: pd.DataFrame, flag_file: Path) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on lower size threshold.

    Parameters:
        gridded_df (pd.DataFrame): Datafram of the gridded catalogue file.
        flag_file (Path): Path to the flag file.

    Returns:
        pandas.DataFrame: Updated flag file with a new boolean column.
    """
    if "FLUX_RADIUS" not in gridded_df.columns:
        raise KeyError("The catalogue does not contain a 'FLUX_RADIUS' column.")

    flag_df = pd.read_csv(flag_file)

    if len(flag_df) != len(gridded_df):
        raise ValueError("Flag file and gridded file have different lengths.")

    if SIZE_ABOVE in flag_df.columns:
        print(f"The flag file already contains a {SIZE_ABOVE} column.")
        return flag_df

    flag_df[SIZE_ABOVE] = gridded_df["FLUX_RADIUS"] > SIZE_LLIM

    flag_df.to_csv(flag_file, index=False)

    print(f"Saved updated flag file with flux radius flags to: {flag_file}")

    return flag_df

def flag_sizeulim(gridded_df: pd.DataFrame, flag_file: Path) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on upper size threshold.

    Parameters:
        gridded_df (pd.DataFrame): Datafram of the gridded catalogue file.
        flag_file (Path): Path to the flag file.

    Returns:
        pandas.DataFrame: Updated flag file with a new boolean column.
    """
    if "FLUX_RADIUS" not in gridded_df.columns:
        raise KeyError("The catalogue does not contain a 'FLUX_RADIUS' column.")

    flag_df = pd.read_csv(flag_file)

    if len(flag_df) != len(gridded_df):
        raise ValueError("Flag file and gridded file have different lengths.")

    if SIZE_BELOW in flag_df.columns:
        print(f"The flag file already contains a {SIZE_BELOW} column.")
        return flag_df
            
    flag_df[SIZE_BELOW] = gridded_df["FLUX_RADIUS"] < SIZE_ULIM

    flag_df.to_csv(flag_file, index=False)

    print(f"Saved updated flag file with flux radius flags to: {flag_file}")

    return flag_df

# def flag_size_kt17(gridded_file: str, flag_file: str) -> pd.DataFrame:
#     """
#     Flags sources in a gridded catalogue based on magnitude threshold.

#     Parameters:
#         gridded_file (str): Path to the gridded catalogue CSV file.

#     Returns:
#         pandas.DataFrame: Updated catalogue with a new boolean column 'MAG_ABOVE_18'.
#     """
#     grid_df = pd.read_csv(gridded_file)

#     if "A_WORLD" not in grid_df.columns:
#         raise KeyError("The catalogue does not contain a 'A_WORLD' column.")
#     if "B_WORLD" not in grid_df.columns:
#         raise KeyError("The catalogue does not contain a 'B_WORLD' column.")
#     if "MAG_AUTO" not in grid_df.columns:
#         raise KeyError("The catalogue does not contain a 'MAG_AUTO' column.")
    
#     flag_df = pd.read_csv(flag_file)
#     if len(flag_df) != len(grid_df):
#         raise ValueError("Flag file and gridded file have different lengths.")
        
#     if SIZE_BETWEEN in flag_df.columns:
#         print(f"The flag file already contains a {SIZE_BETWEEN} column.")
#         return

#     A_arcsec = grid_df["A_WORLD"] * 3600
#     B_arcsec = grid_df["B_WORLD"] * 3600
#     r_rms = np.sqrt(A_arcsec * B_arcsec)

#     mag = grid_df["MAG_AUTO"]
    
#     flag_df[SIZE_BETWEEN] = (mag > (24 - r_rms)) & (mag < (30 - r_rms))

#     flag_df.to_csv(flag_file, index=False)

#     print(f"Saved updated flag file with size flags to: {flag_file}")

#     return flag_df

def flagging(tile: str, gridded_file: Path, flag_file: Path, 
             mask: bool=True, cut: bool=True) -> pd.DataFrame:
    """
    Based on mask and cut, runs programs that adds corresponding 
        flags for each object.

    Parameters:
        tile (str): Name of tile of interest
        gridded_file (Path): Path of gridded catalogue file
        flag_file (Path): Path of flag file
        mask (bool): If the mask is applied, default is True
        cut (bool): If the cut is applied, default is True

    Returns:
        pandas.DataFrame: Returns flag file.
    """
    gridded_df = pd.read_csv(gridded_file)
    
    if not flag_file.exists():
        print(f"Grid flags for {tile} not found, generating files.")
        base_cols = ["i", "j", "X_IMAGE", "Y_IMAGE"]
                     #"A_WORLD", "B_WORLD"]
        flags_df = gridded_df[base_cols].copy()
    
        flags_df.to_csv(flag_file, index=False)

    if mask:
        flag_regions(tile, flag_file)
        
    if cut:
        flag_magllim(gridded_df, flag_file)
        flag_magulim(gridded_df, flag_file)
        flag_sizellim(gridded_df, flag_file)
        flag_sizeulim(gridded_df, flag_file)
        # flag_size_kt17(gridded_file, flag_file)

    return flag_file

def add_good_fraction(tile: str, centers_file: Path) -> pd.DataFrame:
    """
    Adds a 'good_fraction' column to centers_file based on the mask.

    Parameters:
        tile (str): Name of tile of interest.
        centers_file (Path): Path to the grid centers file. 

    Returns:
        pd.DataFrame: Grid centers DataFrame with an added 
            'good_fraction' column.
    """
    mask = get_mask(tile)
    ny, nx = mask.shape

    centers_df = pd.read_csv(centers_file)
    
    required_columns = ["x_min", "x_max", "y_min", "y_max"]
    for col in required_columns:
        if col not in centers_df.columns:
            raise ValueError(f"Input catalogue must have column '{col}'")
    
    good_fractions = []

    for row in centers_df.itertuples(index=False):
        xmin = max(int(row.x_min), 0)
        xmax = min(int(row.x_max), nx)
        ymin = max(int(row.y_min), 0)
        ymax = min(int(row.y_max), ny)
    
        submask = mask[ymin:ymax, xmin:xmax]
        total_pixels = submask.size
        good_pixels = np.sum(submask)
    
        fraction = np.nan
        if total_pixels > 0:
            fraction = good_pixels / total_pixels
    
        good_fractions.append(fraction)

    centers_df["good_fraction"] = good_fractions

    centers_df.to_csv(centers_file, index=False)

    print(f"Added 'good_fraction' column and saved to: {centers_file}")
    
    return centers_df

def grid_counts(tile: str, flag_file: Path, centers_file: Path, 
                mask: bool=True, cut: bool=True,
                mask_exclude:bool=False, cut_exclude:bool=False) -> Path:
    """
    Computes the good counts and normalized counts

    Parameters:
        tile (str): Name of tile of interest
        flag_file (Path): Path of flag file
        centers_file (Path): Path of centers file
        mask (bool): If masks are applied, default is True
        cut (bool): If cuts are applied, default is True
        mask_exclude (bool): If mask selection is the excluded area, default is False
        cut_exclude (bool): If cut selection is the excluded area, default is False

    Returns:
        Path: Path of counts file
    """
    flag_df = pd.read_csv(flag_file)
    grid_info = pd.read_csv(centers_file)

    if mask and "good_fraction" not in grid_info.columns:
        raise ValueError("centers_file must contain 'good_fraction' when mask=True, run add_good_fraction first.")

    if not mask and not cut:
        required_columns = ["i", "j", "X_IMAGE", "Y_IMAGE",]
    elif mask and not cut:
        required_columns = ["i", "j", "X_IMAGE", "Y_IMAGE", "GOOD_REGION"]
    elif cut and not mask:
        required_columns = ["i", "j", "X_IMAGE", "Y_IMAGE", 
                            MAG_ABOVE, MAG_BELOW, SIZE_ABOVE, SIZE_BELOW,]
    else:
        required_columns = ["i", "j", "X_IMAGE", "Y_IMAGE", 
                            "GOOD_REGION", MAG_ABOVE, 
                            MAG_BELOW, SIZE_ABOVE, SIZE_BELOW,]
        
    for col in required_columns:
        if col not in flag_df.columns:
            raise ValueError(f"Input catalogue must have column '{col}'")

    cond = np.ones(len(flag_df), dtype=bool)

    if mask:
        if mask_exclude:
            cond &= ~flag_df["GOOD_REGION"]
        else:
            cond &= flag_df["GOOD_REGION"]

    if cut:
        cut_cond = (
            flag_df[MAG_ABOVE] &
            flag_df[MAG_BELOW] &
            flag_df[SIZE_ABOVE] &
            flag_df[SIZE_BELOW]
        )
        if cut_exclude:
            cond &= ~cut_cond
        else:
            cond &= cut_cond

    selected = flag_df[cond]

    good_counts = (
        selected
        .groupby(["i", "j"])
        .size()
        .reset_index(name="good_count")
        )

    if not mask:
        merged = pd.merge(
            grid_info[["i", "j", "x_center", "y_center", 
                       "ra_center", "dec_center", "l_center", "b_center", ]],
            good_counts, on=["i", "j"], 
            how="left"
        )
        merged["good_count"] = merged["good_count"].fillna(0)
    else:
        merged = pd.merge(
            grid_info[["i", "j", "x_center", "y_center", 
                       "ra_center", "dec_center", "l_center", "b_center", "good_fraction"]], 
            good_counts, on=["i", "j"], 
            how="left"
        )
        merged["good_count"] = merged["good_count"].fillna(0)
    
        merged = merged[merged["good_fraction"] > GOOD_FRACTION_THRESHOLD].copy()
        merged["normalized_count"] = merged["good_count"] / merged["good_fraction"]
    
    output_file = os.path.join(GRID_COUNT_DIR, f"{tile}_grid_counts.csv")
    merged.to_csv(output_file, index=False)
    print(f"Updated grid counts written to {output_file}")

    return output_file