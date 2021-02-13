clear all;
close all;

% addpath /Users/domenicgermano/workspace/Chaste/anim/matlab
addpath /Users/germanod/workspace/Chaste/anim/matlab

% addpath /Users/domenicgermano/workspace/Chaste/projects/results/Test_Periodic_with_bend_large_Anoikis_no_noise_1/results_from_time_0
addpath /Users/germanod/workspace/Chaste/projects/results/Test_Periodic_with_bend_large_Anoikis_no_noise_vel_o/results_from_time_0

filename = 'results';

nodesfile = [filename, '.viznodes'];
celltypesfile = [filename, '.vizcelltypes'];

nodedata = LoadNonConstantLengthData(nodesfile);
celltypesdata = LoadNonConstantLengthData(celltypesfile);

%%
% input these manually...

% tissue size
width = 10;
hight = 14;

ghost = 1 + 1;
stromal = 1;
epithelial = 1;

initial_cells = (stromal + epithelial)*width*hight;

x_centre = 0.5 * width;
y_centre = 0.5 * sqrt(0.75) * hight;


%%
% data structure is:
% time | node_0 | x_0 | y_0 | z_0 | u_0 | v_0 | w_0 | node_1 | x_1 | y_1 | z_1 | u_1 | v_1 | w_1 | ...
%   1  |    2   |  3  |  4  |  5  |  6  |  7  |  8  |    9   |  10 |  11 |  12 |  13 |  14 |  15 | ... 
nodeVelocityRaw = readtable('nodevelocities.dat');

dimensionData = size(nodeVelocityRaw);

timeData = table2array(nodeVelocityRaw(1:dimensionData(1),1));

t1_centre_dist = zeros(dimensionData(1),(dimensionData(2)-1)/7);
t1_velocity_i = zeros(dimensionData(1),(dimensionData(2)-1)/7);


for i = 1:dimensionData(1)
    
    time = timeData(i);
    
    node_ind_i = table2array(nodeVelocityRaw(i,2:7:end));
    
    node_x_i = table2array(nodeVelocityRaw(i,3:7:end));
    node_y_i = table2array(nodeVelocityRaw(i,4:7:end));
    node_z_i = table2array(nodeVelocityRaw(i,5:7:end));

    node_u_i = table2array(nodeVelocityRaw(i,6:7:end));
    node_v_i = table2array(nodeVelocityRaw(i,7:7:end));
    node_w_i = table2array(nodeVelocityRaw(i,8:7:end));
    
    types = celltypesdata{i}(2:end);
    
%     if(i==1)
%         for j=node_ind_i(~isnan(node_ind_i))
%             matlab_node_numb = j+1;
%             t0_centre_dist(i,matlab_node_numb) = sqrt( (node_x_i(matlab_node_numb) - x_centre).^2 + (node_y_i(matlab_node_numb) - y_centre).^2 );
%         
%             t0_velocity_i(i,matlab_node_numb) = 0.0;
%         end
%         previous_time = 0.0;
%     end
    
    
    for j=node_ind_i(~isnan(node_ind_i))
        
        matlab_node_numb = j+1;
        if (matlab_node_numb > (dimensionData(2)-1)/7)
            matlab_node_numb = matlab_node_numb - initial_cells;
        end
        if (ismember(types(matlab_node_numb),[1 2 3 6]))
            t1_centre_dist(i,matlab_node_numb) = sqrt( (node_x_i(matlab_node_numb) - x_centre).^2 + (node_y_i(matlab_node_numb) - y_centre).^2 );

            t1_velocity_i(i,matlab_node_numb) = sqrt( (node_u_i(matlab_node_numb)).^2 + (node_v_i(matlab_node_numb)).^2 );

            % Determine Cell type
            CellColour = [0 0 0];
            if     (types(matlab_node_numb) == 0)
                CellColour = [0 0 0];
            elseif (types(matlab_node_numb) == 1)
                CellColour = [0 0 1];
            elseif (types(matlab_node_numb) == 2)
                CellColour = [1 0 1];
            elseif (types(matlab_node_numb) == 3)
                CellColour = [1 0 0];
            elseif (types(matlab_node_numb) == 6)
                CellColour = [1 0 0];
            end
        end
%         subplot(1,2,1)
%         hold on;
%         title('distance')
%         plot( [previous_time time], [t0_centre_dist(matlab_node_numb) t1_centre_dist(matlab_node_numb)],'-','color',CellColour)
% 
%         subplot(1,2,2)
%         hold on;
%         title('speed')
%         plot( [previous_time time], [t0_velocity_i(matlab_node_numb) t1_velocity_i(matlab_node_numb)],'-','color',CellColour)
%             
            
    end
    
%     t0_centre_dist = t1_centre_dist;
%     t0_velocity_i = t1_velocity_i;
%     previous_time = time;
end

%%

for i = 1:((dimensionData(2)-1)/7)
    
    subplot(1,2,1)
    hold on;
    title('distance')
    plot(t1_centre_dist(:,i))
    
    subplot(1,2,2)
    hold on;
    title('speed')
    plot(t1_velocity_i(:,i))

end