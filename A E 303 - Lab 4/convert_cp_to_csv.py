import os
import pandas as pd

def read_cp_file(filepath):
    """Reads a .cp file from XFOIL and returns a DataFrame with columns x, y, Cp."""
    with open(filepath, 'r') as file:
        lines = file.readlines()[3:]  # Skip first 3 header lines
        data = [list(map(float, line.split())) for line in lines if line.strip()]
    return pd.DataFrame(data, columns=['x', 'y', 'Cp'])

def convert_cp_files_to_csv(folder_path):
    """Converts all .cp files in the folder to .csv format."""
    for filename in os.listdir(folder_path):
        if filename.endswith(".cp"):
            cp_path = os.path.join(folder_path, filename)
            df = read_cp_file(cp_path)

            # Output filename
            csv_name = filename.replace(".cp", ".csv")
            csv_path = os.path.join(folder_path, csv_name)

            df.to_csv(csv_path, index=False)
            print(f"[✓] Converted {filename} → {csv_name}")

if __name__ == "__main__":
    # Set this to the path where your .cp files are
    folder_path = "."  # Current directory; replace with full path if needed
    convert_cp_files_to_csv(folder_path)
