import numpy as np
import sys
import os

def convert_file_to_double(file):
    # Read raw binary data as float32
    float_data = np.fromfile(file, dtype=np.float32)
    print(f"Read {len(float_data)} float32 elements.")

    # Convert to float64 (double)
    double_data = float_data.astype(np.float64)

    # Write float64 data to output file
    double_data.tofile(file)
    print(f"Written float64 data to: {file}")

if __name__ == "__main__":
    if len(sys.argv) < 1:
        print("Usage: python float2double.py <file>")
        sys.exit(1)

    input_filename = sys.argv[1]

    convert_file_to_double(input_filename)
