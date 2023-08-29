This group of codes is to generate multiple extended patterns (foci) with aberration compensation based on the PAC-GS algorithm.

First, run WGS_para_noab_mod.m to generate holograms which corresponding to multiple extended patterns and foci at target positions before aberration correction, save the file: focus_ori.bmp, BR_ori.bmp and ori_3pattern.mat.

Then, run Zernike_base_generation.m to generate Zernike base.

Then run main.m, perform AO test and aberration correction based on PAC-GS, and get holograms corresponding to multiple extended patterns and foci at target positions after aberration correction.

JC, 2023/8/29