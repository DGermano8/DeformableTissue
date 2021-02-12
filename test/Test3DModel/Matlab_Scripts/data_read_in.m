clear all;
close all;

addpath /Users/domenicgermano/workspace/Chaste/anim/matlab

addpath /Users/domenicgermano/workspace/Chaste/projects/results/Test_Periodic_with_bend_large_Anoikis_no_noise_1/results_from_time_0

filename = 'results';

nodesfile = [filename, '.viznodes'];
celltypesfile = [filename, '.vizcelltypes'];

nodedata = LoadNonConstantLengthData(nodesfile);
celltypesdata = LoadNonConstantLengthData(celltypesfile);

%%
% data structure is:
% time | node_0 | x_0 | y_0 | z_0 | u_0 | v_0 | w_0 | node_1 | x_1 | y_1 | z_1 | u_1 | v_1 | w_1 | ...
%   1  |    2   |  3  |  4  |  5  |  6  |  7  |  8  |    9   |  10 |  11 |  12 |  13 |  14 |  15 | ... 
nodeVelocityRaw = readtable('nodevelocities.dat');

dimensionData = size(nodeVelocityRaw);

timeData = table2array(nodeVelocityRaw(1:dimensionData(1),1));

for i = 1:dimensionData(1)
    
    node_ind_i = table2array(nodeVelocityRaw(i,2:7:end));
    
    node_x_i = table2array(nodeVelocityRaw(i,3:7:end));
    node_y_i = table2array(nodeVelocityRaw(i,4:7:end));
    node_z_i = table2array(nodeVelocityRaw(i,5:7:end));

    node_u_i = table2array(nodeVelocityRaw(i,6:7:end));
    node_v_i = table2array(nodeVelocityRaw(i,7:7:end));
    node_w_i = table2array(nodeVelocityRaw(i,8:7:end));
    
end
