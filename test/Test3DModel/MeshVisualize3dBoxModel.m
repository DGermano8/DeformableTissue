%MeshVisualize3dBoxModel(filename)
clear all;
filename = 'results';
%
% function Visualize3dBoxModel(filename)

% Reads in <filename>.viznodes, <filename>.vizc;celltypes, <filename>.vizelements and <filename>.vizsetup
%
% <filename>.viznodes has form:
% (time value 1) (x_1 y_1) (x_2 y_2) ... (x_m y_m)
% (time value 2) (x_1 y_1) (x_2 y_2) ... (x_n y_n)
% .
% .
% N.B. m not necessarily equal to n
%
% <filename>.vizcelltypes has form:
%
% This requires LoadNonConstantLengthData('filename'), which is in the Chaste anim/matlab folder

addpath /Users/germanod/workspace/Chaste/anim/matlab
addpath /Users/germanod/workspace/Chaste/projects/DeformableTissue/results/results_from_time_1

% Visualiser options
RandomColour = false;		% For cells
PlotSprings = true;
PlotTissueCells = true;
PlotCells = false;			% Chinese-lantern type cells
PlotCellCenters = true;
PlotGhosts = true;
NumPtsInCircles = 10;
Defaultcolor = [0.0,0.0,0.7];

nodesfile = [filename, '.viznodes'];
celltypesfile = [filename, '.vizcelltypes'];
elementsfile = [filename, '.vizelements'];
setupfile = [filename, '.vizsetup'];

nodedata = LoadNonConstantLengthData(nodesfile);
celltypesdata = LoadNonConstantLengthData(celltypesfile);
elementsdata = LoadNonConstantLengthData(elementsfile);

if (size(nodedata) ~= size(celltypesdata))
    error('Numbers of timesteps recorded in node and cell type files do not match.')
end

if (RandomColour)
    %Initialise Random Cell colours
    CellColour = rand(1,3);
end

% Setup the unit sphere for drawing cells
[sphereX sphereY sphereZ] = sphere(NumPtsInCircles) ;
sphereX = sphereX;
sphereY = sphereY;
sphereZ = sphereZ;

numtimes = length(nodedata);

%close all
fig=figure;
pause(0.1)
%%

for i = 1;numtimes;        % Timestep is 30 seconds
    clf 

    time = nodedata{i}(1)     % Gives first element of ith row

    types = celltypesdata{i}(2:end);		% Proliferative state
    xvals = nodedata{i}(2:3:end-2);
    yvals = nodedata{i}(3:3:end-1);
    zvals = nodedata{i}(4:3:end);
    NumCells = length(xvals);
    NumNodes = length(types);
    
    %Ele = reshape(elementsdata{i}(2:end)+1,4,length(elementsdata{i}(2:end))/4)';
    Ele = reshape(elementsdata{i}(2:end)+1,4,length(elementsdata{i}(2:end))/4)';
    NumEle = length(Ele(:,1));

    % Plot CellCenters
    
    if (PlotCellCenters)
        for j=1:NumCells
            
            if(j <= length(types))
                if (types(j) == 0)
                    CellColour = [0 0 0];        % Black (stem) epithelial cells

                elseif (types(j) == 1)
                    CellColour = [0 0 1];        % Yellow (transit) epithelial cells

                elseif (types(j) == 2)
                    CellColour = [1 0 1];        % Pink (differentiated) tissue cells
                end
                Colour = CellColour;
    %           plot3(xvals(j),yvals(j), zvals(j), 'o','Color',Colour)
                scatter3(xvals(j),yvals(j), zvals(j), 50,'filled','MarkerFaceColor',Colour)
                hold on
            
            elseif (PlotGhosts)
                CellColour = [0.4 0.4 0.4];  % Grey cells - could be ghosts or apoptotic
                Colour = CellColour;
    %           plot3(xvals(j),yvals(j), zvals(j), 'o','Color',Colour)
                scatter3(xvals(j),yvals(j), zvals(j), 50,'filled','MarkerFaceColor',Colour)
                hold on
            end         
            
%             Colour = CellColour;
% %           plot3(xvals(j),yvals(j), zvals(j), 'o','Color',Colour)
%             scatter3(xvals(j),yvals(j), zvals(j), 50,'filled','MarkerFaceColor',Colour)
%             hold on
        end
    end
    
% Trying to colour the cell centres according to proliferative state    
%     for j = 1:NumNodes
%     
%         if (PlotCellCenters)
%             if (types(j) == 1)
%             plot3(xvals,yvals, zvals, 'o', 'MarkerFaceColor','y', 'MarkerEdgeColor','y')      % Plot epithelial cells in yellow
%             hold on
%             end
%             if (types(j) == 2)
%             plot3(xvals,yvals, zvals, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor','r')      % Plot tissue cells in red
%             hold on
%             end
%         end
%     end

    % Plot Cells
    if (PlotCells)
        CellSeparation=zeros(NumCells,1);
        CellMultiplier=zeros(NumCells,1);
        
        % Calculate cells radius as half average spring length
        for k=1:NumEle
            EdgeIndices = [1,1,1,2,2,3;
                2,3,4,3,4,4];
            for j=1:6
                CellSeparation(Ele(k,EdgeIndices(:,j))) = CellSeparation(Ele(k,EdgeIndices(:,j))) ...
                    + norm([xvals(Ele(k,EdgeIndices(1,j)))-xvals(Ele(k,EdgeIndices(2,j)));
                    yvals(Ele(k,EdgeIndices(1,j)))-yvals(Ele(k,EdgeIndices(2,j)));
                    zvals(Ele(k,EdgeIndices(1,j)))-zvals(Ele(k,EdgeIndices(2,j)))]);
            end
            CellMultiplier(Ele(k,:)) = CellMultiplier(Ele(k,:)) + 3;
        end
        
        for j=1:NumCells
            CellRadius=0.5;%0.5*CellSeparation(j)/CellMultiplier(j);
            CellX = xvals(j) + CellRadius*sphereX;
            CellY = yvals(j) + CellRadius*sphereY;
            CellZ = zvals(j) + CellRadius*sphereZ;
            
            if(RandomColour)
                if (length(CellColour(:,1))<j)
                    CellColour = [CellColour; rand(j - length(CellColour(:,1)),3)];
                end
                Colour = CellColour(j ,:);
            end
            
            % Trying to plot cell colours according to proliferative state
            
            if(j <= length(types))
                if (types(j) == 2)
                    CellColour = [0 0 0];        % Black (stem) epithelial cells

                elseif (types(j) == 1)
                    CellColour = [0 0 1];        % Yellow (transit) epithelial cells

                elseif (types(j) == 0)
                    CellColour = [1 0 1];        % Pink (differentiated) tissue cells
                end
            else
                CellColour = [0.4 0.4 0.4];  % Grey cells - could be ghosts or apoptotic
            end        
            
            Colour = CellColour;
            %if (types(j) <3)
            %if (types(j) == 2 )
                surf (CellX, CellY, CellZ, 'FaceColor', Colour, 'FaceAlpha', 0.25) ;
                hold on
            %end
        end
    end

    % Plot Springs
    if (PlotSprings)
        for k=1:NumEle
            EdgeIndices = [1 2 3 1 4 3 2 4]; 
            plot3(xvals(Ele(k,EdgeIndices)), yvals(Ele(k,EdgeIndices)),...
                  zvals(Ele(k,EdgeIndices)),'b')
        end
    end

    %view(3)
    az = 45;
    el = 5;
    view(az, el);
%    axis([-1,16,-1,16,-1,6]);				
%     axis([-1,7,-1,7,-1,6]);				
    title({sprintf('Time = %2.3f', time);' '})

    %pause(0.1)
    drawnow
    grid off
    M(i) = getframe(fig);
    
    
    pause(0.01);
end
       

% folder = '/Users/germanod/Workspace/Chaste/projects/results_from_time_0/';
% 
% writerObj = VideoWriter( strcat(folder,'BoxModel3DMovie'),'Motion JPEG AVI');
% writerObj.Quality = 100;
% writerObj.FrameRate = 10;
% open(writerObj);
% for i=1:length(M)
%     frame_obj = M(i) ;    
%     writeVideo(writerObj, frame_obj);
% end
% close(writerObj);

%movie2avi(M,'BoxModel3DMovie','fps', 5)