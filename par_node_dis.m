% % % % % % % % % MAIN SCRIPT FOR NODE SETTING: PARALLEL CPUS % % % % % % %

% TODO: dimension-agnostic code
% TODO: even more masks
%% % % % % % % % % % % % PARAMETERS  % % % % % % % % % % % % % % % % % % %

N = 10;                         % number of boxes per side of the cube
maxNodesPerBox = 40;  
repelSteps = 20;               % the number of iterations of the repel.m routine
densityF = @density;                % put the handle to your density function here
kValue = 15;                   % number of nearest neighbors used in the repel.m
A = 2.4;                        % the bigger cube side length
sepRadius = 0.2*A/N;
%%
dim = 3;                        % ATTN: the subsequent code is NOT dimension-independent
repelPower = 5;
bins = 100;
oct = 2^dim;
delta = 1/(32* N * maxNodesPerBox^(1/dim));
cubeShrink = 1 - maxNodesPerBox^(-1/dim)/3;
r1 = sqrt(2);
r2 = sqrt(5);
threshold = .7;   % domain choice threshold (used for strictly positive density)
adjacency = 3^dim;              % the number of nearest boxes to consider
close all;
if ~exist('Output','dir')
    mkdir Output;
end
fileID = fopen('./Output/console.txt','w');
%% populate vertices of the unit cube 
cubeVectors = zeros(dim, oct);                                                                
for i=1:dim
    len = 2 ^ (dim-i);
    for j=0:2^i-1                       
        cubeVectors(i, j*len+1:(j+1)*len) =  mod(j,2)/N;
    end
end

%% node data array and masks
% nodes = zeros(dim, maxNodesPerBox, N^dim); 
% nodeIndices = false(maxNodesPerBox, N^dim);  % logical-styled integer array to be used as mask
% boundaryBoxes = false (N^dim,1);                   % boxes on the boundary
% corner = -ones(dim,1);
           
% corner = [rem((i-1), N);  floor(rem(i-1, N^2)/N);  floor((i-1)/N^2)]/N  ;             
% TODO: this is not dimension-independent



%% main                                               
tic
I=1:N^dim;
corners = A*[rem((I-1), N);  floor(rem(I-1, N^2)/N);  floor((I-1)/N^2)]/N-A/2.0;
[IDX, ~] = knnsearch(corners', corners', 'k', adjacency);
%     Uncomment the following lines to use a domain indicator function.
%     [corners_bool, ~] = in_domain(corners(1,:), corners(2,:),  corners(3,:) );    
cornerIndices = true(1,size(corners,2));  %logical(sum(corners_bool(IDX),2));
cornersUsed = corners(:,cornerIndices);
cornersDensity = densityF(cornersUsed);
cornersAveragedDensity = mean(cornersDensity(IDX'),1) ;        
% averaged density is normalized here to be at most 1

cornersRadii = sepRadius * cornersAveragedDensity;               
currentNumNodes = boxRadiusNum(cornersRadii*N/A);
currentNumNodes = min([currentNumNodes; maxNodesPerBox*ones(1,numel(currentNumNodes))],[],1);
nodes = zeros(dim,sum(currentNumNodes));
previousNodes = [0 cumsum(currentNumNodes)];
    
for i=1:N^dim
%     for j=1:currentNumNodes
%         box(:,j) = A * cubeShrink * [j/currentNumNodes;  mod(r1*j,1);  mod(r2*j,1)]/N; 
%         % the box is shrunk to account for inwards corner % TODO: can this be vectorized?
%         box(:,j) = box(:,j) + corner + delta;
%     end
    J = 1:currentNumNodes(i);
    box = A * cubeShrink * [J/currentNumNodes(i);  mod(r1*J,1);  mod(r2*J,1)]/N;
    box = bsxfun(@plus, corners(:,i)+delta, box);
     
    nodes(:,previousNodes(i)+1:previousNodes(i+1)) = box;   
end
toc

%% remove nodes outside the density support
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
%% node stats
% fprintf( '\nNumber of nodes:      %d\n',  length(cnf))
% fprintf( 'Mean number of nodes per box:      %d\n', mean(sum(nodeIndices,1) ))
% fprintf( 'Max number of nodes per box:      %d\n', max(sum(nodeIndices,1) ))
% fprintf( 'Min number of nodes per box:      %d\n', min(sum(nodeIndices,1) ))
% toc
% fprintf('\n');
% clf;
% clf;
% pbaspect([1 1 1])
% view([1 1 0])
% figure(2);
% plot3(cnf(1,:), cnf(2,:), cnf(3,:),  '.k');
 
% repel and save nodes
% fprintf( fileID, 'Performing %d repel steps.\n',  repelSteps);
fprintf( 'Performing %d repel steps.\n',  repelSteps)
cnf = repel(cnf, kValue, repelSteps, @in_domain, repelPower, densityF, fileID);
toc 


pbaspect([1 1 1])
% view([1 1 0])
F = figure(1);
plot3(cnf(1,:), cnf(2,:), cnf(3,:),  '.k');

% savefig(F,'./Output/nodes','compact')
% save('slanttrui.mat', 'cnf')
% dlmwrite('./Output/cnf.txt',cnf','delimiter','\t'); % ,'precision',3
