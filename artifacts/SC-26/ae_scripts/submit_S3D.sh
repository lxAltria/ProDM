data_dict_path=$1

sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path CH4
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path CO
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path CO2
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path H2O
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path N2
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path O2
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path Pressure
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path Temperature
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path VelocityX
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path VelocityY
sbatch ./ae_scripts/S3D_ablation_study.sh $data_dict_path VelocityZ
