#!/bin/bash
#SBATCH --time=03:00:00                 # Time limit for the job (adjust as needed)
#SBATCH --job-name=S3D_encoding_time    # Job name
#SBATCH --nodes=1                       # Number of nodes to allocate
#SBATCH --ntasks=16                     # Number of cores to allocate
#SBATCH --partition=                    # Partition/queue to run the job in
#SBATCH -e slurm-%j.err                 # Error file for this job
#SBATCH -o slurm-%j.out                 # Output file for this job
#SBATCH -A                              # Project allocation account name
#SBATCH --mail-type=end                 # Send email when job starts/ends
#SBATCH --mail-user=                    # Email address for notifications

variables=(CH4 CO CO2 H2O O2 Temperature VelocityX VelocityY VelocityZ)

data_dict_path=$1

output_file="./ae_results/S3D_Results/S3D_encoding_time.txt"
tmp_file="./ae_results/S3D_Results/S3D_encoding_time_tmp.txt"

> $output_file
> $tmp_file

for var in "${variables[@]}"; do
    data_file="${data_dict_path}/data/${var}.dat"
    refactor_file="${data_dict_path}/refactor/${var}_refactored"
    ./build/test/test_mdr_fuse_refactor $data_file "-d" 4 60 3 500 500 500 $refactor_file 3 > $tmp_file
    encoding_time=$(grep "Encoding" $tmp_file | head -n 1)
    echo "Method#1, old perbit, Variable=$var, $encoding_time" >> $output_file
    ./build/test/test_mdr_fuse_refactor $data_file "-d" 4 60 3 500 500 500 $refactor_file 2 > $tmp_file
    encoding_time=$(grep "Encoding" $tmp_file | head -n 1)
    echo "Method#2, perbit, Variable=$var, $encoding_time" >> $output_file
done

rm $tmp_file