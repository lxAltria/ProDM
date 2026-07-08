data_dict_path=$1

variables=(CLDICE CLDLIQ CLOUD CMFDQR CMFDQ CMFDT CONCLD DCQ DTCOND DTV FICE GCLDLWP ICIMR ICLDIWP ICLDTWP ICWMR OMEGAT OMEGA QC QRL QRS Q RELHUM T UU U VD01 VQ VT VU VV V Z3)

for var in "${variables[@]}"; do
    sbatch ./ae_scripts/CESM_ablation_study.sh $data_dict_path $var
done