# Improving Progressive Compression with Adaptive Interpolation and Coefficient Decomposition

This is the code repo for the paper "Improving Progressive Compression with Adaptive Interpolation and Coefficient Decomposition".
<br />

# Installation

One-command compilation using "sh build_script.sh". It will automatically builds ProDM libaray and the dependencies.
```
git clone https://github.com/Linus-Li-1037/ProDM.git
cd ProDM
git switch Predict
sh build_script.sh
```
# Examples
**Template**<br />

```
cd build
Refactor: ./test/two_modes_refactor data_file -[dataType: f/d] target_level num_bitplanes num_dims dim1 dim2 ... dimn output_path -[encoder_option: Nega/XOR/PerBit] -[prior_mode: eb(default)/PSNR] -[CP_or_not: CP/no_CP] (coeff_interp_direction, default tune) 
Retrieval: ./test/two_modes_reconstructor data_file -[dataType: f/d] num_of_tolerance tol1 tol2 ... toln refactored_path -[encoder_option: Nega/XOR/PerBit] -[interpreter_option: Greedy/DP/BFS] -[CP_or_not: CP/no_CP] [Optional: Reconstructed data path]
```

# Artifacts
Please use evaluation.py for reproducing the results in paper.<br />

