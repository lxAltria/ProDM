data_dict_path=$1

sh ./ae_scripts/submit_CESM.sh "$data_dict_path/CESM"
sh ./ae_scripts/submit_Miranda.sh "$data_dict_path/Miranda"
sh ./ae_scripts/submit_SCALE.sh "$data_dict_path/SCALE"
sh ./ae_scripts/submit_S3D.sh "$data_dict_path/S3D"
