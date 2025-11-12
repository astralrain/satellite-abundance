from libs import *
from consts import *
from grid_funcs import *
from ang_funcs import *
from basic_funcs import *

def process_galaxy(galaxy_file: str, output_dir: str):
    """
    Takes in a csv galaxy_file containing a list of primaries and goes 
        through each, creating a gridded catalogue for all tiles the 
        galaxies lies on and determines the count and distance.

    Parameters:
        galaxy_file (str): CSV file of galaxies and their information
            - must include ra, dec, z, and tile name
    
    Returns:
        None, writes one or more outputs to CSV files.
    """

    galaxy_df = pd.read_csv(galaxy_file)
    number = 0
    count = 0
    for _, gal in galaxy_df.iterrows():
        tile = gal["tile_name"]
        cat_file = os.path.join(CAT_DIR, tile + ".cat")
        grid_file = Path(GRID_DIR) / f"{tile}_gridded.csv"
        grid_center_file = Path(GRID_CENTER_DIR) / f"{tile}_grid_centers.csv"
        
        if not grid_file.exists():
            print(f"Gridded catalog for {tile} not found, generating files.")
            cat_df = read_cat(cat_file)
            grid_tile(tile, cat_df)
            flagging(tile, grid_file)
        
        if not grid_center_file.exists():
            print(f"Grid centers for {tile} not found, generating files.")
            cat_df = read_cat(cat_file)
            compute_grid_centers(tile)
            add_good_fraction(tile, grid_center_file)
            add_grid_counts(grid_file, grid_center_file)
        
        compute_primary_info(grid_center_file, gal, output_dir)
        number += 1
        print(f"Processed {number}/1860 galaxies")
        count += 1

        if count >= 186:
            break
    return

def compute_primary_info(grid_center_file: str, gal: pd.DataFrame,
                         output_dir: str) -> None:
    """
    Compute angular separation from a primary's ra and dec taken from 
        gal to each of the grid centers in grid_file and write results 
        to a CSV in output_dir.

    Parameters:
        grid_file (str): File of the grid information.
        primary_file (str): File of the primary.
        output_dir (str): Intended directory of output.

    Returns:
        None, writes output to a csv file.
    """
    grid_df = pd.read_csv(grid_center_file)
    results = []

    sep = angular_separation(gal["ra"], gal["dec"],
                             grid_df["ra_center"].values,
                             grid_df["dec_center"].values)

    z = gal["z"]
    ang_d = angular_diameter_distance(z) 

    proj_dist = projected_distance(sep, ang_d)

    for i, (s, projd) in enumerate(zip(sep, proj_dist)):
        row = {**gal.to_dict(),
               **grid_df.iloc[i].to_dict(),
               "angular_separation_deg": s,
               "projected_distance_mpc": projd}
        results.append(row)

    out_df = pd.DataFrame(results)

    # exclude primary
    x = gal["X_IMAGE"]
    y = gal["Y_IMAGE"]

    inside = (
        (out_df["x_min"] <= x) & (x < out_df["x_max"]) &
        (out_df["y_min"] <= y) & (y < out_df["y_max"])
        )

    if inside.any():
        idx = out_df.index[inside][0]
        out_df.loc[idx, "normalized_count"] = (
            out_df.loc[idx, "normalized_count"] - 1
            )
        print(f"Primary galaxy located in grid cell {idx}: count adjusted by -1.")
    else:
        print("Warning: Primary galaxy not found inside any grid cell!")

    primary_name = gal[SDSS_ID]
    output_file = os.path.join(output_dir, f"{primary_name}_info.csv")
    out_df.to_csv(output_file, index=False)
    print(f"Saved separations to {output_file}")

process_galaxy(Path("filter_matched_obj/sdss_with_matches_obj.csv"), SDSS_PRIMARY_DIR)