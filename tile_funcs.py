from libs import *
from basic_funcs import *
from consts import *
from ang_funcs import *

def get_good_tile_fraction(tile, grid: bool = False):
    mask = get_mask(tile)
    total = mask.size
    if total == 0:
        return float("nan")
    return float(np.sum(mask) / total)

def process_tiles(tile_file, tile_info_df):
    if not os.path.exists("tile_information.csv"):
        pd.DataFrame(columns=["tile", "ra", "dec", "good_fraction"]) \
          .to_csv("tile_information.csv", index=False)
    
    tile_list_df = pd.read_csv(tile_file)
    tiles_df = pd.read_csv("cfis_r_tiles.csv")
    bad_tiles = load_bad_tiles("bad_tiles.csv")
    
    tile_frac = dict(zip(tile_info_df["tile"], tile_info_df["good_fraction"]))
    tile_coords = dict(zip(tiles_df["tile"], zip(tiles_df["ra"], tiles_df["dec"])))

    filtered_tiles: set[str] = set()
    new_rows = []

    for tile in tile_list_df["tile"]:
        if tile in bad_tiles:
            print(f"SKIP (bad tile list): {tile}")
            continue
        
        if tile in tile_frac:
            tile_good_frac = tile_frac[tile]
        else:
            tile_good_frac = get_good_tile_fraction(tile)

            if tile not in tile_coords:
                raise ValueError(f"{tile} not found in cfis_r_tiles.csv")
            
            ra, dec = tile_coords[tile]
            new_rows.append({"tile": tile,
                             "ra": ra,
                             "dec": dec,
                             "good_fraction": tile_good_frac
                             })
            
            tile_frac[tile] = tile_good_frac
    
        if tile_good_frac < TILE_GOOD_FRACTION_THRESHOLD:
            print(f"SKIP (good_fraction={tile_good_frac:.3f}): {tile}")
            continue

        filtered_tiles.add(tile)

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        tile_info_df = pd.concat([tile_info_df, new_df], ignore_index=True)

        tile_info_df = tile_info_df.drop_duplicates(subset="tile", keep="last")
        tile_info_df = tile_info_df.sort_values(by=["ra", "dec"])

        tile_info_df.to_csv("tile_information.csv", index=False)
    print("Filtered tiles")       
    return filtered_tiles, tile_info_df

def filter_tile(tile):
    """
    Given a catalogue, keep only objects in GOOD regions (mask == True)
    and save to a new filtered catalogue.

    Parameters:
        tile (str): Tile name
        cat_file (Path): Path to input catalogue
        filtered_cat_dir (Path): Directory to save filtered catalogues

    Returns:
        pd.DataFrame: Filtered catalogue
    """
    cat_path = download_lsb_cat(tile)
    cat_df = read_cat(cat_path)
    os.remove(cat_path)
    mask = get_mask(tile)

    x = cat_df["X_IMAGE"].astype(int)
    y = cat_df["Y_IMAGE"].astype(int)

    valid = (x >= 0) & (x < TILE_SIZE) & (y >= 0) & (y < TILE_SIZE)

    good = pd.Series(False, index=cat_df.index)

    good.loc[valid] = mask[y[valid], x[valid]]

    filtered_df = cat_df[good].copy()

    dir_path = Path(FILTERED_CAT_DIR)
    dir_path.mkdir(parents=True, exist_ok=True)

    out_file = Path(FILTERED_CAT_DIR) / (tile + FILTERED_CAT_END)

    filtered_df.to_csv(out_file, index=False)

    print(f"Saved filtered catalogue: {out_file} ({len(filtered_df)} objects)")

    return filtered_df

def compute_tile_counts(tile, cat_df, tile_info_df, mag_llim, mag_ulim, size_llim, cut=True):
    row = tile_info_df[tile_info_df["tile"] == tile].iloc[0]

    good_frac = row["good_fraction"]
    tile_ra = row["ra"]
    tile_dec = row["dec"]

    cond = np.ones(len(cat_df), dtype=bool)

    if cut:
        cond &= (cat_df["MAG_AUTO"].values > mag_llim)
        cond &= (cat_df["MAG_AUTO"].values < mag_ulim)
        cond &= (cat_df["FLUX_RADIUS"].values > size_llim)

    good_count = int(cond.sum())
    normalized_count = good_count / good_frac if good_frac > 0 else 0

    return {
        "tile": tile,
        "ra": tile_ra,
        "dec": tile_dec,
        "good_fraction": good_frac,
        "good_count": good_count,
        "normalized_count": normalized_count,
    }

def process_tile_counts(tile_file, counts_file, mag_llim, mag_ulim, size_llim):
    tile_info_df = pd.read_csv("tile_information.csv")

    tiles, tile_info_df = process_tiles(tile_file, tile_info_df)
    clear_output(wait=True)

    tile_rows = []

    for tile in tiles:
        cat_path = Path(FILTERED_CAT_DIR) / f"{tile}{FILTERED_CAT_END}"

        if cat_path.exists():
            cat_df = pd.read_csv(cat_path)
        else:
            cat_df = filter_tile(tile)

        row = compute_tile_counts(
            tile,
            cat_df,
            tile_info_df,
            mag_llim,
            mag_ulim,
            size_llim
        )
        print(row)
        tile_rows.append(row)

    counts_df = pd.DataFrame(tile_rows)
    counts_df = counts_df.sort_values(by=["ra", "dec"])
    counts_df.to_csv(counts_file, index=False)
    clear_output(wait=True)

    return counts_df

def match_tiles_to_primaries(valid_tiles_df: pd.DataFrame, primaries_file: str, output_file:str, max_separation_deg: float,):
    """
    Match tiles to nearest primary galaxies within a given angular radius.

    Parameters:
        valid_tiles_df (pd.DataFrame): Tile table (should include ra, dec, tile, etc.)
        primaries_file (str): CSV of primary galaxies
        max_separation_deg (float): Matching radius in degrees
        tile_primary_dir (str): Output directory

    Returns:
        pd.DataFrame: Combined matched tile-primary table
    """
    os.makedirs(PRIMARY_DIR, exist_ok=True)

    primaries_df = pd.read_csv(primaries_file)

    for col in ("tile", "ra", "dec"):
        if col not in valid_tiles_df.columns:
            raise ValueError(f"valid_tiles_df must contain '{col}'")

    for col in (PRIMARY_ID, GAL_RA, GAL_DEC):
        if col not in primaries_df.columns:
            raise ValueError(f"primaries_file must contain '{col}'")
    
    if "Deff" in primaries_df.columns:
        dist_mode = "deff"
    elif "z" in primaries_df.columns:
        dist_mode = "redshift"
    else:
        raise ValueError("primaries_file must contain either 'Deff' or 'z'")

    tile_ras  = valid_tiles_df["ra"].values
    tile_decs = valid_tiles_df["dec"].values

    prim_ras  = primaries_df[GAL_RA].values
    prim_decs = primaries_df[GAL_DEC].values

    print(
        f"Matching {len(valid_tiles_df)} tiles against "
        f"{len(primaries_df)} primaries (radius={max_separation_deg}°)..."
    )

    seps = angular_separation(
        tile_ras[:, np.newaxis],
        tile_decs[:, np.newaxis],
        prim_ras[np.newaxis, :],
        prim_decs[np.newaxis, :],
    )

    nearest_prim_idx = np.argmin(seps, axis=1)
    nearest_sep = seps[np.arange(len(valid_tiles_df)), nearest_prim_idx]

    within_radius = nearest_sep <= max_separation_deg

    print(
        f"{within_radius.sum()} tiles matched | "
        f"{(~within_radius).sum()} dropped"
    )

    tiles_in = valid_tiles_df[within_radius].copy().reset_index(drop=True)
    prim_idx_in = nearest_prim_idx[within_radius]
    sep_in = nearest_sep[within_radius]

    matched_prims = primaries_df.iloc[prim_idx_in].reset_index(drop=True)

    proj_dists = (
            projected_distance(sep_in, matched_prims["Deff"].values)
            if dist_mode == "deff"
            else projected_distance(sep_in, angular_diameter_distance(matched_prims["z"].values))
        )
    tiles_in["primary_id"] = matched_prims[PRIMARY_ID].values
    tiles_in["primary_ra"] = matched_prims[GAL_RA].values
    tiles_in["primary_dec"] = matched_prims[GAL_DEC].values
    tiles_in["ang_sep_deg"] = sep_in
    tiles_in["proj_dist_mpc"] = proj_dists

    combined_chunks = []

    for pid, group in tiles_in.groupby("primary_id"):
        out_path = os.path.join(PRIMARY_DIR, f"{pid}_tiles.csv")
        group.to_csv(out_path, index=False)
        print(f"{pid}: {len(group)} tiles → {out_path}")
        combined_chunks.append(group)

    if not combined_chunks:
        print("WARNING: No matches found.")
        return pd.DataFrame()

    combined_df = pd.concat(combined_chunks, ignore_index=True)

    combined_path = os.path.join(os.getcwd(), output_file)
    combined_df.to_csv(combined_path, index=False)

    print(f"\nCombined file ({len(combined_df)} rows) → {combined_path}")

    return combined_df