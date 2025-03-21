# QPro: An Efficient Framework for Quantity-of-Interest Based Progressive Retrieval with Guaranteed Error Control

This is the code repo for the paper "QPro: An Efficient Framework for Quantity-of-Interest Based Progressive Retrieval with Guaranteed Error Control".

Major authors: Dr. Xin Liang (UK), Dr. Qing Liu (NJIT), Dr. Xubin He (Temple)<br />
Other contributors: Xuan Wu (UK), Qirui Tian (NJIT), Wenbo Li (UK)<br />
Collaborators: Dr. Scott Klasky (ORNL), Dr. Qian Gong (ORNL), Dr. Jill Zhang (LLNL), Dr. Seung-Hoe Ku (PPPL), Dr. Xiaohua Zhang (LLNL), Dr. Jieyang Chen (UAB) etc.<br />

# Installation

One-command compilation using "sh build_script.sh". It will automatically builds ProDM libaray and the dependencies.

git clone https://github.com/Linus-Li-1037/ProDM.git<br />
cd ProDM<br />
sh build_script.sh<br />

# Examples
**Template**<br />
Precision data refactoring using approximators:<br />
```
cd build
Refactor: ./test/refactor_d64 1 1 $dataset $path_to_dataset $max_weight_v $block_size
Retrieval: ./test/qoi_{$target QoI}_d64 1 1 $eb $path_to_dataset
```

Taking Hurricane dataset as an example, Hurricane ISABEL can be downloaded from https://sdrbench.github.io/.<br /> 
Double versions of VelocityX.dat, VelocityY.dat, and VelocityZ.dat are extended from Uf48.bin.f32, Vf48.bin.f32, and Wf48.bin.f32 of Hurricane ISABEL.<br />
We arrange the dataset like followings:<br />
```
Hurricane_d64
├── data
│   ├── VelocityX.dat
│   ├── VelocityY.dat
│   └── VelocityZ.dat
└── refactor
    ├── VelocityX_refactored
    ├── VelocityY_refactored
    └── VelocityZ_refactored
```

Thus, the Template can be modified into the following commands to test with Hurricane ISABEL using V_total as targeted QoI:<br />
```
cd build
Refactor: ./test/refactor_d64 1 1 Hurricane ../Hurricane_d64 7 4
Retrieval: ./test/qoi_Vtot_d64 1 1 $eb ../Hurricane_d64
```

# Q&A

Please address your questions to xliang@uky.edu with subject title ProDM<br />
