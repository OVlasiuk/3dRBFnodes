% % % % % % % % % MAIN SCRIPT FOR NODE SETTING: PARALLEL CPUS % % % % % % %

% TODO: dimension-agnostic code

%% % % % % % % % % % % % PARAMETERS  % % % % % % % % % % % % % % % % % % %

N = 50;                         % number of boxes per side of the cube
maxNodesPerBox = 40;  
repelSteps = 20;                % the number of iterations of the repel.m routine
densityF = @density;            %     put the handle to your density function here
kValue = 20;                    % number of nearest neighbors used in the repel.m
A = 6;                          % The outer cube sidelength; all will be 
                                %  contained in [-A/2, A/2]^3. 
jitter = 0;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
dim = 3;                        % ATTN: the subsequent code is NOT dimension-independent
repelPower = 5;                 
oct = 2^dim;
delta = 1/(256* N * maxNodesPerBox^(1/dim));
cubeShrink = 1 - maxNodesPerBox^(-1/dim)/64;
r1 = sqrt(2);
r2 = (sqrt(5)-1)/(sqrt(2));
threshold = .7;                 % domain choice threshold (used for strictly positive density)
adjacency = 3^dim;              % the number of nearest boxes to consider
close all;

s = char(mfilename('fullpath'));
cd(s(1:end-12))                         % cd to the mfile folder; 
                                        % The constant 12 depends on the
                                        % length of the filename.
if ~exist('Output','dir')
    mkdir Output;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
try 
    load('./Output/unit_lattice_radius.mat','DELTA', 'CUBE_SHRINK','R1', 'R2')
    if (delta ~= DELTA)...
            || (cubeShrink ~= CUBE_SHRINK)...
            || (r1 ~= R1)...
            || (r2 ~= R2)
        throw(MException('ReadTable:NoFile','I could not find the table of radii.'));
    end
catch
    fprintf('\nLooks like the interpolation table for this number of lattice\n');
    fprintf('nodes is missing or not up to date... Hang on there, I''ll make\n'); 
    fprintf('a new one for you. This may take a few minutes, but we''ll only\n');
    fprintf('do it once.\n');
    lattice_by_count(2000,delta,cubeShrink,r1,r2,'y');
    fprintf('...\nDone.\n\n')
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Populate vertices of the unit cube 
cubeVectors = zeros(dim, oct);                                                                
for i=1:dim
    len = 2 ^ (dim-i);
    for j=0:2^i-1                       
        cubeVectors(i, j*len+1:(j+1)*len) =  mod(j,2)/N;
    end
end
           


%% Main                                               
tic
I=1:N^dim;
corners = A*[rem((I-1), N);  floor(rem(I-1, N^2)/N);  floor((I-1)/N^2)]/N-A/2.0;
[IDX, ~] = knnsearch(corners', corners', 'k', adjacency);
%     Uncomment the following lines to use a domain indicator function.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     [corners_bool, ~] = in_domain(corners(1,:), corners(2,:),  corners(3,:) );    
cornerIndices = true(1,size(corners,2));  %logical(sum(corners_bool(IDX),2));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

cornersUsed = corners(:,cornerIndices);
cornersDensity = densityF(cornersUsed);
cornersAveragedDensity = mean(cornersDensity(IDX'),1);        
cornersRadii = cornersAveragedDensity;
currentNumNodes = num_radius(cornersRadii*N/A);
% % % Note that the maximum number of nodes is capped in the following
% line: not doing it can cause up to a thousand nodes in a single box.
currentNumNodes = min([currentNumNodes; maxNodesPerBox*ones(1,numel(currentNumNodes))],[],1);
nodes = zeros(dim,sum(currentNumNodes));
previousNodes = [0 cumsum(currentNumNodes)];
    
for i=1:N^dim
    J = 1:currentNumNodes(i);
    box = A * cubeShrink * [J/currentNumNodes(i);  mod(r1*J,1);  mod(r2*J,1)]/N;
    box = bsxfun(@plus, corners(:,i)+delta, box);
     
    nodes(:,previousNodes(i)+1:previousNodes(i+1)) = box;   
end
toc

%% Remove nodes outside the density support
% bdry_nodes = nodes(:, :, boundaryBoxes);  
% coordinates of the nodes that belong to boundary boxes
% flat = num2cell(reshape(bdry_nodes, dim, []), 1);                                      
% 
% f_vals = densityF(nodes);                                                        
% values of the density function at those nodes
% bdry_removed_indices = reshape( f_vals > threshold, maxNodesPerBox, []);            
% indices of nodes that have to be removed; hard-coded threshold condition
% nodeIndices(:, boundaryBoxes) = nodeIndices(:, boundaryBoxes) .* ~bdry_removed_indices; 
% indices in the global nodeset
% node_indices_flat = reshape(nodeIndices, 1, []);  
% nodes = reshape (nodes, dim, []);                                                     
% all nodes in a single array;      dim x num possible nodes
% cnf = nodes(:, node_indices_flat);                                                    
% after removing nodes with zero density
cnf = nodes; 


%% Node stats
outtemp = length(cnf);
% fprintf( fileID, '\nNumber of nodes:      %d\n',  outtemp);
fprintf( '\nNumber of nodes:      %d\n',  outtemp)
fprintf( 'Mean number of nodes per box:      %d\n', mean(currentNumNodes ))
fprintf( 'Max number of nodes per box:      %d\n', max(currentNumNodes ))
fprintf( 'Min number of nodes per box:      %d\n', min(currentNumNodes ))
toc
fprintf('\n');
 
%% Repel and save nodes
% fprintf( fileID, 'Performing %d repel steps.\n',  repelSteps);
fprintf( 'Performing %d repel steps.\n',  repelSteps)
cnf = repel(cnf,kValue,repelSteps,A,0,densityF,jitter,repelPower,0);
 

%% Plot the results
pbaspect([1 1 1])
% view([1 1 0])
F = figure(1);
msize = ceil(max(1, 22-3.8*log10(size(cnf,2)) ));
plot3(cnf(1,:), cnf(2,:), cnf(3,:),'.k','MarkerSize',msize);


% savefig(F,'./Output/nodes','compact')
% save('slanttrui.mat', 'cnf')
% dlmwrite('./Output/cnf.txt',cnf','delimiter','\t'); % ,'precision',3