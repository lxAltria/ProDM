#!/bin/bash
#SBATCH --time=03:00:00                 # Time limit for the job (adjust as needed)
#SBATCH --job-name=S3D_decompose_time   # Job name
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

output_file="./ae_results/S3D_Results/S3D_decompose_time.txt"
tmp_file="./ae_results/S3D_Results/S3D_decompose_time_tmp.txt"

> $output_file
> $tmp_file

for var in "${variables[@]}"; do
    data_file="${data_dict_path}/data/${var}.dat"
    ./external/MGARDx/build/test/test_decompose $data_file 1 4 3 500 500 500 > $tmp_file
    decompose_time=$(grep "Decomposition" $tmp_file | head -n 1)
    echo "Method#1, PMGARD, Variable=$var, $decompose_time" >> $output_file
    ./external/MGARDx/build/test/test_decompose_interleave $data_file 1 4 3 500 500 500 > $tmp_file
    decompose_time=$(grep "Decomposition" $tmp_file | head -n 1)
    echo "Method#2, Fast, Variable=$var, $decompose_time" >> $output_file
done

rm $tmp_file