# main.py
import os
import sys
import shutil
from provided_code.solver import linel_fem_solver

INPUT_DIR = "inputs"
OUTPUT_DIR = "outputs"

def ensure_dir(path):
    """Create directory if it doesn't exist."""
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")

def main():
    # 1. Ensure input_files exists
    ensure_dir(INPUT_DIR)

    # 2. List Python files in input_files
    input_files = [f for f in os.listdir(INPUT_DIR) if f.endswith(".py")]
    if not input_files:
        print(f"Error: No .py input files found in '{INPUT_DIR}/'.")
        sys.exit(1)

    # 3. Ensure output directory exists
    ensure_dir(OUTPUT_DIR)

    # 4. Create output subfolders matching input file names
    for file_name in input_files:
        base_name = os.path.splitext(file_name)[0]  # remove .py extension
        output_subdir = os.path.join(OUTPUT_DIR, base_name)
        ensure_dir(output_subdir)

        src_file = os.path.join(INPUT_DIR, file_name)
        dest_file = os.path.join(output_subdir, file_name)

        shutil.copy2(src_file, dest_file)  # copy with metadata
        print(f"Copied {file_name} → {output_subdir}")
        try:
            linel_fem_solver(dest_file)
        except Exception as e:
            print(f"Error running solver for '{file_name}': {e}")

if __name__ == "__main__":
    main()
