from libs import *
from consts import *
from grid_funcs import *
from ang_funcs import *
from basic_funcs import *

def combine_primary_info(input_dir: str, output_file: str) -> None:
    """
    Combine all individual primary information CSV files in a directory
    into one combined CSV.

    Parameters:
        input_dir (str): Directory containing all individual *_info.csv files.
        output_file (str): Path to the output combined CSV.

    Returns:
        None
    """
    csv_files = sorted(Path(input_dir).glob("*_info.csv"))
    if not csv_files:
        print(f"No CSV files found in {input_dir}")
        return

    combined = []
    for file in csv_files:
        try:
            df = pd.read_csv(file)
            combined.append(df)
        except Exception as e:
            print(f"Skipping {file.name}: {e}")

    if not combined:
        print("No valid dataframes to combine.")
        return

    combined_df = pd.concat(combined, ignore_index=True)
    combined_df.to_csv(output_file, index=False)
    print(f"Combined {len(csv_files)} files into {output_file}")

combine_primary_info("sdss_primaries_10", "subset_sdss_primaries_10.csv")
