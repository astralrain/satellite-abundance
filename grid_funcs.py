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

            grid_data.append({
                "i": i,
                "j": j,
                "x_min": x0, "x_max": x1,
                "y_min": y0, "y_max": y1,
                "x_center": x_center,
                "y_center": y_center,
                "ra_center": ra_center,
                "dec_center": dec_center
            })

    grid_info = pd.DataFrame(grid_data)
    output_file = os.path.join(GRID_CENTER_DIR, tile + "_grid_centers.csv")
    grid_info.to_csv(output_file, index=False)

    print(f"Grid centers written to {output_file}")

    return output_file

def flag_regions(tile: str, gridded_file: str) -> pd.DataFrame:
    """
    Given a catalogue and its tile, mark which objects fall in good or 
    bad regions using masks.

    Parameters:
        tile (str): Name of tile of interest.
        gridded_file (str): File of the gridded catalogue.

    Returns:
        The pandas DataFrame of the updated catalogue.
    """
    df = pd.read_csv(gridded_file)

    mask = get_mask(tile)
    ny, nx = mask.shape

    x = df["X_IMAGE"].astype(int)
    y = df["Y_IMAGE"].astype(int)

    df["GOOD_REGION"] = False

    valid = (x >= 0) & (x < nx) & (y >= 0) & (y < ny)

    df.loc[valid, "GOOD_REGION"] = mask[y[valid], x[valid]]

    df.to_csv(gridded_file, index=False)
    
    print(f"Saved updated catalogue with mask info to: {gridded_file}")
    
    return df

def flag_magllim(gridded_file: str) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on magnitude threshold.

    Parameters:
        gridded_file (str): Path to the gridded catalogue CSV file.

    Returns:
        pandas.DataFrame: Updated catalogue with a new boolean column 'MAG_ABOVE_18'.
    """
    df = pd.read_csv(gridded_file)

    if "MAG_AUTO" not in df.columns:
        raise KeyError("The catalogue does not contain a 'MAG_AUTO' column.")

    df[f"MAG_ABOVE_{MAG_LLIM}"] = df["MAG_AUTO"] > MAG_LLIM

    df.to_csv(gridded_file, index=False)

    print(f"Saved updated catalogue with magnitude flags to: {gridded_file}")

    return df

def flag_magulim(gridded_file: str) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on magnitude threshold.

    Parameters:
        gridded_file (str): Path to the gridded catalogue CSV file.

    Returns:
        pandas.DataFrame: Updated catalogue with a new boolean column 'MAG_ABOVE_18'.
    """
    df = pd.read_csv(gridded_file)

    if "MAG_AUTO" not in df.columns:
        raise KeyError("The catalogue does not contain a 'MAG_AUTO' column.")

    df[f"MAG_BELOW_{MAG_ULIM}"] = df["MAG_AUTO"] < MAG_ULIM

    df.to_csv(gridded_file, index=False)

    print(f"Saved updated catalogue with magnitude flags to: {gridded_file}")

    return df

def flag_size(gridded_file: str) -> pd.DataFrame:
    """
    Flags sources in a gridded catalogue based on magnitude threshold.

    Parameters:
        gridded_file (str): Path to the gridded catalogue CSV file.

    Returns:
        pandas.DataFrame: Updated catalogue with a new boolean column 'MAG_ABOVE_18'.
    """
    df = pd.read_csv(gridded_file)

    if "FLUX_RADIUS" not in df.columns:
        raise KeyError("The catalogue does not contain a 'FLUX_RADIUS' column.")

    df[f"SIZE_ABOVE_{SIZE_LIM}"] = df["FLUX_RADIUS"] > SIZE_LIM

    df.to_csv(gridded_file, index=False)

    print(f"Saved updated catalogue with flux radius flags to: {gridded_file}")

    return df

def flagging(tile: str, gridded_file: str):
    flag_regions(tile, gridded_file)
    flag_magllim(gridded_file)
    flag_magulim(gridded_file)
    flag_size(gridded_file)

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
        fraction = good_pixels / total_pixels if total_pixels > 0 else np.nan

        good_fractions.append(fraction)

    df["good_fraction"] = good_fractions

    df.to_csv(centers_file, index=False)

    print(f"Added 'good_fraction' column and saved to: {centers_file}")
    return df

def add_grid_counts(gridded_file: str, centers_file: str) -> None:
    """
    Add normalized object counts per grid cell into centers_file,
    using the GOOD_REGION flag from gridded_file and correcting for coverage fraction.

    Parameters:
        gridded_file (str): Path to the input catalogue CSV.
            Must have columns 'i', 'j', 'X_IMAGE', 'Y_IMAGE', 'GOOD_REGION'.
        centers_file (str): Path to the existing grid info CSV.
            Must have columns 'i', 'j'.

    Returns:
        None (writes to CSV).
    """
    gridded_df = pd.read_csv(gridded_file)
    grid_info = pd.read_csv(centers_file)

    required_columns = ["i", "j", "X_IMAGE", "Y_IMAGE", 
                        "GOOD_REGION"]
    for col in required_columns:
        if col not in gridded_df.columns:
            raise ValueError(f"Input catalogue must have column '{col}'")

    good_counts = (
        gridded_df[gridded_df["GOOD_REGION"]]
        .groupby(["i", "j"])
        .size()
        .reset_index(name="good_count")
        )

    merged = pd.merge(grid_info, good_counts, on=["i", "j"], how="left")

    merged["good_count"] = merged["good_count"].fillna(0)
    merged["normalized_count"] = (
        (merged["good_count"] / merged["good_fraction"])
        .replace([np.inf, np.nan], 0)
        .astype(int)
        )

    merged.to_csv(centers_file, index=False)
    
    print(f"Updated grid counts written to {centers_file}")
