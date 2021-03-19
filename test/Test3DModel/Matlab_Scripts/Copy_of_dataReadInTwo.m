clear all;
close all;

% addpath /Users/domenicgermano/Workspace/ChasteDom/anim/matlab
addpath /Users/germanod/workspace/ChasteDom/anim/matlab

% addpath /Users/domenicgermano/Workspace/results/FlatWntRadius_0303_06/results_from_time_0
addpath /Users/germanod/Workspace/results/FlatWntRadius_0803_01/results_from_time_0

filename = 'results';

nodesfile = [filename, '.viznodes'];
celltypesfile = [filename, '.vizcelltypes'];

nodedata = LoadNonConstantLengthData(nodesfile);
celltypesdata = LoadNonConstantLengthData(celltypesfile);

%%
% input these manually...

% tissue size
width = 14;
hight = 18;

ghost = 1 + 1;
stromal = 1;
epithelial = 1;

initial_cells = (stromal + epithelial)*width*hight;

x_centre = 0.5 * width;
y_centre = 0.5 * sqrt(0.75) * hight;
% x_centre = 0.5 * width * 0.96;
% y_centre = 0.5 * sqrt(0.75) * hight * 0.96;

max_radius = sqrt((x_centre)^2 + (y_centre)^2);


%%

% data structure is:
% time | node_0 | x_0 | y_0 | z_0 | u_0 | v_0 | w_0 | node_1 | x_1 | y_1 | z_1 | u_1 | v_1 | w_1 | ...
%   1  |    2   |  3  |  4  |  5  |  6  |  7  |  8  |    9   |  10 |  11 |  12 |  13 |  14 |  15 | ... 
nodeVelocityRaw = readtable('nodevelocities.dat');

CellIdRaw = readtable('loggedcell.dat');

dimensionData = size(nodeVelocityRaw);

timeData = table2array(nodeVelocityRaw(1:dimensionData(1),1));

holder_init = str2double(split(table2array(nodeVelocityRaw(end,2))));

CellIdholder = str2double(split(table2array(CellIdRaw(end,2))));

t1_centre_dist = zeros(1,(length(holder_init))/7);
t1_velocity_i = zeros(1,(length(holder_init))/7);

% t1_centre_dist = zeros(1,(dimensionData(2)-1)/7);
% t1_velocity_i = zeros(1,(dimensionData(2)-1)/7);

node_i_radial_velocity_with_time = zeros(CellIdholder(end-3),length(timeData));


number_to_do = 20;
interval = 20;

radius_intervals = zeros(number_to_do,101);
counter_intervals = zeros(number_to_do,101);
count = 1;
for i = (length(timeData)- number_to_do*interval+interval):1:length(timeData);
    
    time = timeData(i);
    
    holder = str2double(split(table2array(nodeVelocityRaw(i,2))));
%     holder = table2array(nodeVelocityRaw(i,2:end));
    holder_CellId = str2double(split(table2array(CellIdRaw(i,2))));

    holder_cell_id = holder_CellId(1:5:end);
    holder_node_id = holder_CellId(2:5:end);
    
    node_ind_i = holder(1:7:end)';
    
    node_x_i = holder(2:7:end);
    node_y_i = holder(3:7:end);
    node_z_i = holder(4:7:end);

    node_u_i = holder(5:7:end);
    node_v_i = holder(6:7:end);
    node_w_i = holder(7:7:end);
    
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
    
    
%     for j=node_ind_i(~isnan(node_ind_i))
    for j=node_ind_i(~isnan(node_ind_i))
        
        matlab_node_numb = j+1;
        if (matlab_node_numb > (length(holder_init))/7)
        %if (matlab_node_numb > ((dimensionData(2)-1)/7) ) 
            matlab_node_numb = matlab_node_numb - initial_cells;
        end
        
        if (ismember(types(matlab_node_numb),[1 2 3 6]))
            xt = node_x_i(matlab_node_numb) - x_centre;
            yt = node_y_i(matlab_node_numb) - y_centre;
            
            t1_centre_dist(matlab_node_numb) = sqrt( (xt).^2 + (yt).^2 );

            t1_velocity_i(matlab_node_numb) = (xt*node_u_i(matlab_node_numb) + yt*node_v_i(matlab_node_numb))/(sqrt( (xt).^2 + (yt).^2 ));
        end
        
        prev_int = 0.0;
        for k=1:100
            this_int = (k/100)*max_radius;
            
            if (t1_centre_dist(matlab_node_numb) >= prev_int && t1_centre_dist(matlab_node_numb) < this_int)
                radius_intervals(count,k) = radius_intervals(count,k) + t1_velocity_i(matlab_node_numb);
                counter_intervals(count,k) = radius_intervals(count,k) + 1;
                
            end
            prev_int = this_int;
        end
        if (t1_centre_dist(matlab_node_numb) >= prev_int )
            radius_intervals(count,101) = radius_intervals(count,101) + t1_velocity_i(matlab_node_numb);
            counter_intervals(count,101) = radius_intervals(count,101) + 1;
        end
       
    end
    
    count = count + 1;
    
end
%%
% scatter(t1_centre_dist,(t1_velocity_i))
scatter(t1_centre_dist,(t1_velocity_i))
xlabel('Distance from centre')
ylabel('Radial velocity')


%%
figure
av_velocity = zeros(1,101);
for k = 1:101
    av_velocity(k) = sum(radius_intervals(:,k))/sum(counter_intervals(:,k));
end
plot(0:((1/100)*max_radius):((1)*max_radius),av_velocity,'ko')

figure
for k = 1:number_to_do
hold on
plot(0:((1/100)*max_radius):((1)*max_radius),radius_intervals(k,:),'o','color',[0.8 0.8 0.8])
end
plot(0:((1/100)*max_radius):((1)*max_radius),av_velocity,'ko','linewidth',2)

%%


for i = 1:(length(holder_init)/7)
    for j=1:(length(t1_centre_dist(:,1))-1)
        if (t1_centre_dist(j,i) == 0)
            t1_centre_dist_post(j,i) = NaN;
        elseif abs(t1_centre_dist(j+1,i) - t1_centre_dist(j,i)) > 0.5
            t1_centre_dist_post(j,i) = NaN;
        else
            t1_centre_dist_post(j,i) = t1_centre_dist(j,i);
        end
        
        if (t1_velocity_i(j,i) == 0)
            t1_velocity_i_post(j,i) = NaN;
        else
            t1_velocity_i_post(j,i) = t1_velocity_i(j,i);
        end
        
    end
    
    j = length(t1_centre_dist(:,1));
        if (t1_centre_dist(j,i) == 0)
            t1_centre_dist_post(j,i) = NaN;
        else
            t1_centre_dist_post(j,i) = t1_centre_dist(j,i);
        end
        
        if (t1_velocity_i(j,i) == 0)
            t1_velocity_i_post(j,i) = NaN;
        else
            t1_velocity_i_post(j,i) = t1_velocity_i(j,i);
        end
    
end

%%
figure;
for i = 1:(length(holder_init)/7)
    
    subplot(1,2,1)
    hold on;
    title('distance')
    plot(t1_centre_dist_post(:,i))
    
    subplot(1,2,2)
    hold on;
    title('speed')
    plot(t1_velocity_i_post(:,i))

end




%%

% data structure is:
% time | node_0 | x_0 | y_0 | z_0 | u_0 | v_0 | w_0 | node_1 | x_1 | y_1 | z_1 | u_1 | v_1 | w_1 | ...
%   1  |    2   |  3  |  4  |  5  |  6  |  7  |  8  |    9   |  10 |  11 |  12 |  13 |  14 |  15 | ... 
nodeVelocityRaw = readtable('nodevelocities.dat');

CellIdRaw = importdata('loggedcell.dat');

dimensionData = size(nodeVelocityRaw);

timeData = table2array(nodeVelocityRaw(1:dimensionData(1),1));

% holder_init = str2double(split(table2array(nodeVelocityRaw(end,2))));
% holder_init = str2double(split(table2array(nodeVelocityRaw(end,2:end))));

CellIdholder = str2double(split(CellIdRaw(end)));

number_of_cells_end = CellIdholder(end-4);

node_i_radial_velocity_with_time = zeros(number_of_cells_end,length(timeData));
node_i_radial_position_with_time = zeros(number_of_cells_end,length(timeData));

for i = 1:length(timeData)
    
    time = timeData(i);
    
%     size_12 = size(nodeVelocityRaw(i,2:end));
%     holder = zeros(1,size_12(2));
%     for ii = 1:size_12(2)
%         hold_val = table2array(nodeVelocityRaw(i,ii));
%         if iscell(hold_val)
%             holder(ii) = NaN;
%         else
%             holder(ii) = hold_val;
%         end
%     end
    holder = str2double(split(table2array(nodeVelocityRaw(i,2:end))));
%     holder = table2array(nodeVelocityRaw(i,2:end));
    holder_CellId = str2double(split(CellIdRaw(i)));

    holder_cell_id = holder_CellId(2:5:end);
%     holder_node_id = holder_CellId(2:5:end);
    
    node_ind_i = holder(1:7:end)';
    
    node_x_i = holder(2:7:end);
    node_y_i = holder(3:7:end);
    node_z_i = holder(4:7:end);

    node_u_i = holder(5:7:end);
    node_v_i = holder(6:7:end);
    node_w_i = holder(7:7:end);
    
    types = celltypesdata{i}(2:end);
    
    for j=1:length(node_ind_i)
        
        node_i = holder_cell_id(j);
        
        if (ismember(types(j),[1 2 3 6]))
            
            xt = node_x_i(j) - x_centre;
            yt = node_y_i(j) - y_centre;
            
            node_i_radial_position_with_time(node_i,i) = sqrt( (xt).^2 + (yt).^2 );

            node_i_radial_velocity_with_time(node_i,i) = (xt*node_u_i(j) + yt*node_v_i(j))/(sqrt( (xt).^2 + (yt).^2 ));
        end
        
    end
end
%%
delta_t = 0.01;
sample_interval = 5;
when_to_sample = 7200;
when_to_stop = 9000;
average_velocity = zeros(number_of_cells_end,length(when_to_sample:sample_interval:when_to_stop));
average_position =  zeros(number_of_cells_end,length(when_to_sample:sample_interval:when_to_stop));
figure;
for j=(hight*width):number_of_cells_end
    count = 0;
    for i=1:length(timeData)
        if node_i_radial_position_with_time(j,i) == 0
            node_i_radial_position_with_time(j,i) = NaN;
        end
    end
    for i=when_to_sample:sample_interval:when_to_stop
        count = count + 1;
        
        if count == 1
            average_velocity(j,count) = (1/(delta_t*sample_interval))*(node_i_radial_position_with_time(j,i) - node_i_radial_position_with_time(j,i-sample_interval));
%             average_velocity(j,count) = node_i_radial_velocity_with_time(j,i);
%             average_velocity(j,count) = mean(node_i_radial_velocity_with_time(j,(i-sample_interval):i));
%             average_position(j,count) = node_i_radial_position_with_time(j,i);
            average_position(j,count) = mean(node_i_radial_position_with_time(j,(i-sample_interval):i));
        else
            average_velocity(j,count) = (1/(delta_t*sample_interval))*(node_i_radial_position_with_time(j,i) - average_position(j,count-1));
%             average_velocity(j,count) = node_i_radial_velocity_with_time(j,i);
%             average_velocity(j,count) = mean(node_i_radial_velocity_with_time(j,(i-sample_interval):i));
%             average_position(j,count) = node_i_radial_position_with_time(j,i);
            average_position(j,count) = mean(node_i_radial_position_with_time(j,(i-sample_interval):i));
        end
%         if mod(i,sample_interval)==0 && i > sample_interval+1
            
            
%             average_velocity(j,count) = (1/(delta_t*sample_interval))*(node_i_radial_position_with_time(j,i) - node_i_radial_position_with_time(j,i-sample_interval));
%             average_position(j,count) = node_i_radial_position_with_time(j,i);
%         elseif i ==1
%             average_velocity(j,1) = NaN;
%             average_position(j,1) = NaN;
%         elseif i ==2
%             average_velocity(j,2) = NaN;
%             average_position(j,2) = NaN;
%         end
    end
    subplot(2,1,1)
    hold on;
    plot(node_i_radial_position_with_time(j,:))
    xlabel('time')
    ylabel('radial distance')
    subplot(2,1,2)
    hold on;
    plot(average_position(j,:),average_velocity(j,:),'o-')
    plot([0 7], [0 0] ,'k-')
    xlabel('radial distance')
    ylabel('radial velocity')
end



%%

region_count = zeros(1,2*ceil(max_radius));
region_velocity =  zeros(1,2*ceil(max_radius));
region_mean = zeros(1,2*ceil(max_radius));
figure;
for j=(hight*width):number_of_cells_end
    
    for i=1:2*ceil(max_radius)
        low_bound = (i-1)/2;
        upp_bound = (i-1)/2+1;
        
        if( (average_position(j,end) > low_bound) && (average_position(j,end) <= upp_bound))
            region_count(i) = region_count(i) + 1;
            region_velocity(i) = region_velocity(i) + (average_velocity(j,end));
        end
    end
    
end

for i=1:2*ceil(max_radius)
    low_bound = (i-1)/2;
    upp_bound = (i-1)/2+1;
    
    region_mean(i) = 0.5*(low_bound+upp_bound);
    if(region_count(i) ~=0)
        region_velocity(i) = region_velocity(i)/region_count(i);
    end
end

hold on
plot(region_mean,region_velocity,'o-')
plot([0 upp_bound], [0 0] ,'k-')
xlabel('radial distance')
ylabel('radial velocity')
hold off

%%
region_count = zeros(length(when_to_sample:sample_interval:when_to_stop),2*ceil(max_radius));
region_velocity =  zeros(length(when_to_sample:sample_interval:when_to_stop),2*ceil(max_radius));
region_mean = zeros(length(when_to_sample:sample_interval:when_to_stop),2*ceil(max_radius));
figure;


for j=(hight*width):number_of_cells_end
    count = 0;
    for time=when_to_sample:sample_interval:when_to_stop
        count = count + 1;
        for i=1:2*ceil(max_radius)
            low_bound = (i-1)/2;
            upp_bound = (i-1)/2+1;

            if( (average_position(j,count) > low_bound) && (average_position(j,count) <= upp_bound))
                region_count(count,i) = region_count(count,i) + 1;
                region_velocity(count,i) = region_velocity(count,i) + (average_velocity(j,count));
            end
        end
    end
    
end

count = 0;
for time=when_to_sample:sample_interval:when_to_stop
    count = count + 1;
    for i=1:2*ceil(max_radius)
        low_bound = (i-1)/2;
        upp_bound = (i-1)/2+1;

        region_mean(count,i) = 0.5*(low_bound+upp_bound);
        if(region_count(count,i) ~=0)
            region_velocity(count,i) = region_velocity(count,i)/region_count(count,i);
        end
    end
end

[X,Y] = meshgrid(1:2*ceil(max_radius),when_to_sample:sample_interval:when_to_stop);
hold on
s=surf(X,Y,region_velocity);
% s.EdgeColor = 'none';
% plot([0 upp_bound], [0 0] ,'k-')
xlabel('radial distance')
ylabel('time')
zlabel('radial velocity')
hold off