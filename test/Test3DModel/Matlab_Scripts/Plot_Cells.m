clear all;
filename = 'results';

addpath /Users/germanod/workspace/Chaste/anim/matlab
% addpath /Users/domenicgermano/workspace/Chaste/anim/matlab

folder = '/Users/germanod/workspace/Chaste/projects/results/Test_writter_7/results_from_time_0';
% folder = '/Users/domenicgermano/workspace/Chaste/projects/results/Test_Periodic_with_bend_large_Anoikis_no_noise_2/results_from_time_0';

addpath /Users/germanod/workspace/Chaste/projects/results/Test_writter_7/results_from_time_0
% addpath /Users/domenicgermano/workspace/Chaste/projects/results/Test_Periodic_with_bend_large_Anoikis_no_noise_2/results_from_time_0

nodesfile = [filename, '.viznodes'];
celltypesfile = [filename, '.vizcelltypes'];
elementsfile = [filename, '.vizelements'];
setupfile = [filename, '.vizsetup'];
cellmutationstatefile = [filename, '.vizmutationstates'];

nodedata = LoadNonConstantLengthData(nodesfile);
celltypesdata = LoadNonConstantLengthData(celltypesfile);
elementsdata = LoadNonConstantLengthData(elementsfile);
% cellmutationstatedata = LoadNonConstantLengthData(cellmutationstatefile);

numtimes = length(nodedata);

x_centre = 0.5 * 10;
y_centre = 0.5 * sqrt(0.75) * 14;

% tissue size
width = 10;
hight = 10;

ghost = 1 + 1;
stromal = 1;
epithelial = 1;

initial_cells = (stromal + epithelial)*width*hight;

ghost_number = ghost * width * hight;

%%
figure;
t0_dist = zeros(1,length(nodedata{1}(2:3:end-2)));
previous_time = 0;
xvals = nodedata{1}(2:3:end-2);
yvals = nodedata{1}(3:3:end-1);
% for j=1:length(nodedata{1}(2:3:end-2)) %[1 50 100 150 200 250 280]
%  
%     t0_dist(j) = sqrt( (xvals(j) - x_centre).^2 + (yvals(j) - y_centre).^2 );
%         
% end
t1_dist = t0_dist;

for i = 1:numtimes
    
    % clf;
    
    time = nodedata{i}(1)  ;   % Gives first element of ith row

    types = celltypesdata{i}(2:end);		% Proliferative state
%     mutation = cellmutationstatedata{i}(2:end);		% Proliferative state
    xvals = nodedata{i}(2:3:end-2);
    yvals = nodedata{i}(3:3:end-1);
    zvals = nodedata{i}(4:3:end);
    NumCells = length(xvals);
    
    if i ==1
        t0_xvals = xvals;
        t0_yvals = yvals;
        
        t1_xvals = xvals;
        t1_yvals = yvals;
    end
    
    
    additional_it = [];
    if (initial_cells + ghost_number < NumCells)
        
        additional_it = (1:NumCells-(initial_cells + ghost_number)) + initial_cells;
        
    end
    
    cell_ids = [1:initial_cells additional_it];
    t2_dist = zeros(1,length(cell_ids));

    for j=cell_ids%[1 50 100 150 200 250 280]  -> 280 is a magic number! this is to do with nodes and cell numbers not alinging
        
        if(j <= length(types))
            CellColour = [0 0 0];
        	if (types(j) == 0)
                CellColour = [0 0 0];        % Black (stem) epithelial cells

            elseif (types(j) == 6)
                CellColour = [1 0 0];        % Yellow (transit) epithelial cells
            elseif (types(j) == 1)
            	CellColour = [0 0 1];
            elseif (types(j) == 2)
            	CellColour = [1 0 1];        % Pink (differentiated) tissue cells
            elseif(types(j) == 3)
            	CellColour = [1 0 0]; 
            end
            Colour = CellColour;
    %           plot3(xvals(j),yvals(j), zvals(j), 'o','Color',Colour)
%             scatter3(xvals(j),yvals(j), zvals(j), 50,'filled','MarkerFaceColor',Colour)
            %hold on
            
        elseif (0) % PlotGhosts
            CellColour = [0.4 0.4 0.4];  % Grey cells - could be ghosts or apoptotic
            Colour = CellColour;
    %           plot3(xvals(j),yvals(j), zvals(j), 'o','Color',Colour)
%             scatter3(xvals(j),yvals(j), zvals(j), 50,'filled','MarkerFaceColor',Colour)
%         	hold on
        end
        
%         if ismember(types(j),[1 2 3 6])
%             hold on;
%             scatter3(xvals(j),yvals(j), zvals(j), 50,'filled','MarkerFaceColor',Colour)
%         end
        
%         t2_dist(j) = sqrt( (xvals(j) - x_centre).^2 + (yvals(j) - y_centre).^2 );
        
        t2_dist(j) = sqrt( (xvals(j) - t1_xvals(j)).^2 + (yvals(j) - t1_yvals(j)).^2 );
        
        % Periodic displacement
        if (abs(t2_dist(j)) > 1.0)
            t2_dist(j) = 0.0;
        end

        
        if (ismember(types(j),[1 2 3 6]) &&  length(t2_dist)==length(t1_dist) && length(t0_dist)==length(t2_dist))
            hold on;
            if ismember(types(j),[1 2 3 6])

            % distance itself
            subplot(1,2,1)
            hold on;
            title('distance')
%             plot( [previous_time time], [t1_dist(j) t2_dist(j)],'-','color',Colour)
            plot( [previous_time time], [t1_dist(j) t2_dist(j)],'-','color',Colour)

            % maybe we need the derivative of the distance?
            subplot(1,2,2)
            hold on;
            title('speed')
            plot( [previous_time time], [t1_dist(j) - t0_dist(j) t2_dist(j) - t1_dist(j)],'-','color',Colour)
        
%             if ismember(types(j),[1])
%                 subplot(2,2,3)
%                 hold on;
%                 plot( [t0_xvals(j) xvals(j)], [t0_yvals(j) yvals(j)],'-','color',Colour)
           end
        end
        
    end
    % drawnow;
    
    t0_dist = t1_dist;
    t1_dist = t2_dist;
    previous_time = time;
    
    t0_xvals = t1_xvals;
    t0_yvals = t1_yvals;
    
    t1_xvals = xvals;
    t1_yvals = yvals;
    
    
end