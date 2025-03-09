# ProDM: A Unified Progressive Data Management Library

This is the code repo for NSF project "Collaborative Research: Elements: ProDM: Developing A Unified Progressive Data Management Library for Exascale Computational Science". It is a joint collaborative effort from the University of Kentucky (UK), New Jersey Institute of Technology (NJIT), and Temple University. 

Major authors: Dr. Xin Liang (UK), Dr. Qing Liu (NJIT), Dr. Xubin He (Temple)<br />
Other contributors: Xuan Wu (UK), Qirui Tian (NJIT), Wenbo Li (UK)<br />
Collaborators: Dr. Scott Klasky (ORNL), Dr. Qian Gong (ORNL), Dr. Jill Zhang (LLNL), Dr. Seung-Hoe Ku (PPPL), Dr. Xiaohua Zhang (LLNL), Dr. Jieyang Chen (UAB) etc.<br />

# Installation

One-command compilation using "sh build_script.sh". It will automatically builds ProDM libaray and the dependencies.

git clone https://https://github.com/lxAltria/ProDM.git<br />
cd ProDM<br />
sh build_script.sh<br />

# Examples

Multilevel data refactoring using PMGARD:<br />
cd build<br />
mkdir -p refactored_data<br />
Refactor: ./test/test_mdr_refactor $data_file $num_level $num_bitplanes $num_dims $dim0 $dim1 $dim2<br />
Retrieval: ./test/test_pdr_reconstructor $data_file $num_tolerance $tolerance_1 ...<br />

Precision data refactoring using approximators:<br />
cd build<br />
mkdir -p refactored_data<br />
Refactor: ./test/test_pdr_refactor $data_file $num_level $num_bitplanes $num_dims $dim0 $dim1 $dim2<br />
Retrieval: ./test/test_pdr_reconstructor $data_file $num_tolerance $tolerance_1 ...<br />

# Q&A

Please address your questions to xliang@uky.edu with subject title ProDM<br />
