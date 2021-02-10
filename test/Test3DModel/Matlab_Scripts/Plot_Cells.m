clear all;
filename = 'results';

addpath /Users/germanod/workspace/Chaste/anim/matlab

folder = '/Users/germanod/workspace/Chaste/projects/results/Test_Periodic_with_bend_large_Anoikis_no_noise/results_from_time_0';

addpath /Users/germanod/workspace/Chaste/projects/results/Test_Periodic_with_bend_large_Anoikis_no_noise/results_from_time_0



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

fig=figure;
pause(0.1)

x_centre = 0.5 * 10;
y_centre = 0.5 * sqrt(0.75) * 14;

%%
previous_dist_from_centre = zeros(1,length(nodedata{1}(2:3:end-2)));
previous_time = 0;
xvals = nodedata{1}(2:3:end-2);
yvals = nodedata{1}(3:3:end-1);
for j=1:length(nodedata{1}(2:3:end-2)) %[1 50 100 150 200 250 280]
 
    previous_dist_from_centre(j) = sqrt( (xvals(j) - x_centre).^2 + (yvals(j) - y_centre).^2 );
        
end
    
for i = 1:100
    
    % clf;
    
    time = nodedata{i}(1)  ;   % Gives first element of ith row

    types = celltypesdata{i}(2:end);		% Proliferative state
%     mutation = cellmutationstatedata{i}(2:end);		% Proliferative state
    xvals = nodedata{i}(2:3:end-2);
    yvals = nodedata{i}(3:3:end-1);
    zvals = nodedata{i}(4:3:end);
    NumCells = length(xvals)
    
    dist_from_centre = zeros(1,NumCells);
    
    for j=[1 50 100 150 200 250 280]
            
        if(j <= length(types))
        	if (types(j) == 0)
                CellColour = [0 0 0];        % Black (stem) epithelial cells

            elseif (types(j) == 6)
                CellColour = [0 0 1];        % Yellow (transit) epithelial cells

            elseif (types(j) == 2)
            	CellColour = [1 0 1];        % Pink (differentiated) tissue cells
            elseif(types(j) == 3)
            	CellColour = [1 0 0]; 
            end
            Colour = CellColour;
    %           plot3(xvals(j),yvals(j), zvals(j), 'o','Color',Colour)
            %scatter3(xvals(j),yvals(j), zvals(j), 50,'filled','MarkerFaceColor',Colour)
            %hold on
            
        elseif (0) % PlotGhosts
            CellColour = [0.4 0.4 0.4];  % Grey cells - could be ghosts or apoptotic
            Colour = CellColour;
    %           plot3(xvals(j),yvals(j), zvals(j), 'o','Color',Colour)
            scatter3(xvals(j),yvals(j), zvals(j), 50,'filled','MarkerFaceColor',Colour)
        	hold on
        end
        
        dist_from_centre(j) = sqrt( (xvals(j) - x_centre).^2 + (yvals(j) - y_centre).^2 );
        hold on;
        plot( [previous_time time], [previous_dist_from_centre(j) dist_from_centre(j)])
    
    end
    % drawnow;
    previous_dist_from_centre = dist_from_centre;
    previous_time = time;
    
    
end