import numpy as np
import sys
import os

def convert_file_to_double(input_file, output_file):
    # Read raw binary data as float32
    float_data = np.fromfile(input_file, dtype=np.float32)
    print(f"Read {len(float_data)} float32 elements.")

    # Convert to float64 (double)
    double_data = float_data.astype(np.float64)

    # Write float64 data to output file
    double_data.tofile(output_file)
    print(f"Written float64 data to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python float2double.py <file.f32> [output_file]")
        sys.exit(1)

    input_filename = sys.argv[1]
    if len(sys.argv) > 2:
        output_filename = sys.argv[2]
    elif input_filename.endswith(".f32"):
        # strip the .f32 suffix: VelocityX.dat.f32 -> VelocityX.dat
        output_filename = input_filename[:-4]
    else:
        output_filename = input_filename

    convert_file_to_double(input_filename, output_filename)
