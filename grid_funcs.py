from libs import *
from consts import *
from basic_funcs import *

def grid_tile(cat_name: str, cat_df: pd.DataFrame, ngrid: int = 10,
              padding_fraction: float = 0.0) -> pd.DataFrame:
    """
    Grids a CFIS tile using pixel coordinates.

    Parameters:
        cat_name (str): Name of the tile/catalogue of interest.
        cat_df (pd.DataFrame): DataFrame containing the catalogue file of the tile.
        ngrid (int): Number of grid divisions in X and Y, default = 10.
        padding_fraction (float): Padding beyond min/max extents, as a fraction, default = 0.0.

    Returns:
        Output CSV path of gridded catalogue
    """

    x_min, x_max = 0.0, TILE_SIZE
    y_min, y_max = 0.0, TILE_SIZE

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

    output_file = os.path.join(GRID_DIR, f"{cat_name}_gridded.csv")
    cat_df.to_csv(output_file, index=False)
    
    print(f"Pixel-based gridded catalogue written to {output_file}")

    return output_file

def compute_grid_centers(cat_name: str, cat_df: pd.DataFrame, ngrid: int = 10):
    """
    Compute RA/Dec centers of each grid using pixel coordinates and WCS.
    Grids are defined in pixel space, then converted to sky coordinates.

    Parameters:
        cat_name (str): Name of the tile/catalogue of interest.
        cat_df (pd.DataFrame): DataFrame containing the catalogue file of the tile.
        ngrid (int): Number of grid divisions in X and Y, default = 10.
    
    Returns:
        Output CSV path of the grid centers file.
    """

    wcs, nx, ny = load_header(cat_name)
    
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
    output_file = os.path.join(GRID_CENTER_DIR, cat_name + "_grid_centers.csv")
    grid_info.to_csv(output_file, index=False)

    print(f"Grid centers written to {output_file}")

    return output_file

def flag_regions(cat_name: str, gridded_file: str) -> pd.DataFrame:
    """
    Given a catalogue and its tile, mark which objects fall in good or bad regions using masks.

    Parameters:
        cat_name (str): Name of tile of interest.
        gridded_file (str): File of the gridded catalogue.

    Returns:
        The pandas DataFrame of the updated catalogue.
    """
    df = pd.read_csv(gridded_file)

    mask = get_mask(cat_name)
    ny, nx = mask.shape

    x = df["X_IMAGE"].astype(int)
    y = df["Y_IMAGE"].astype(int)

    df["GOOD_REGION"] = False

    valid = (x >= 0) & (x < nx) & (y >= 0) & (y < ny)

    df.loc[valid, "GOOD_REGION"] = mask[y[valid], x[valid]]

    df.to_csv(gridded_file, index=False)
    
    print(f"Saved updated catalogue with mask info to: {gridded_file}")
    
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

    required_columns = ["i", "j", "X_IMAGE", "Y_IMAGE", "GOOD_REGION"]
    for col in required_columns:
        if col not in gridded_df.columns:
            raise ValueError(f"Input catalogue must have column '{col}'")

    x = gridded_df["X_IMAGE"].astype(int)
    y = gridded_df["Y_IMAGE"].astype(int)

    total_counts = gridded_df.groupby(["i", "j"]).size().reset_index(name="total_count")
    good_counts = gridded_df[gridded_df["GOOD_REGION"]].groupby(["i", "j"]).size().reset_index(name="good_count")
    merged_counts = pd.merge(total_counts, good_counts, on=["i", "j"], how="left")
    
    merged_counts["good_count"] = merged_counts["good_count"].fillna(0)

    merged_counts["good_fraction"] = merged_counts["good_count"] / merged_counts["total_count"]
    merged_counts["good_fraction"] = merged_counts["good_fraction"].replace(0, np.nan)

    merged_counts["normalized_count"] = merged_counts["good_count"] / merged_counts["good_fraction"]
    merged_counts["normalized_count"] = merged_counts["normalized_count"].fillna(0).astype(int)

    merged = pd.merge(grid_info, merged_counts[["i", "j", "normalized_count", "good_fraction"]],
                      on=["i", "j"], how="left")

    merged["normalized_count"] = merged["normalized_count"].fillna(0).astype(int)
    merged["good_fraction"] = merged["good_fraction"].fillna(0)

    merged.to_csv(centers_file, index=False)
    
    print(f"Updated grid counts written to {centers_file}")