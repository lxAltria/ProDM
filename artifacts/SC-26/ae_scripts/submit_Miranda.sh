data_dict_path=$1

variables=(density diffusivity pressure velocityx velocityy velocityz viscocity)

for var in "${variables[@]}"; do
   sbatch ./ae_scripts/Miranda_ablation_study.sh $data_dict_path $var
done