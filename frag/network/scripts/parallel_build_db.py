import argparse

if __name__ == "__main__":

    # Read in a SD or SMILES file - then write out into a specified directory
    parser = argparse.ArgumentParser(
        description="Convert a SMILES or SDFile to input for Astex Fragment network."
    )
    parser.add_argument("--input")
    parser.add_argument("--base_dir")
    args = parser.parse_args()
    # Now spawn up
    sdf_file_list = ["sdf_one.sdf", "sdf_two.sdf"]
    output_file_list = ["output_one", "output_two"]
