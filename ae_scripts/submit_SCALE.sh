data_dict_path=$1

variables=(PRES QC QG QI QR QS QV RH T U V W)

for var in "${variables[@]}"; do
   sbatch ./ae_scripts/SCALE_ablation_study.sh $data_dict_path $var
done
