% % % % % % % % % MAIN SCRIPT FOR NODE SETTING: PARALLEL CPUS % % % % % % %

% TODO: dimension-agnostic code
% TODO: even more masks
%% % % % % % % % % % % % PARAMETERS  % % % % % % % % % % % % % % % % % % %

N = 20;                         % number of boxes per side of the cube
max_nodes_per_box = 15;  
repel_steps = 10;               % the number of iterations of the repel.m routine
density = @trui;                % put the handle to your density function here
k_value = 40;                   % number of nearest neighbors used in the repel.m
%%
dim = 3;                        % ATTN: the subsequent code is NOT dimension-independent
repel_power = 7;
oct = 2^dim;
delta = 1/(8* N * max_nodes_per_box^(1/dim));
cube_shrink = 1 - max_nodes_per_box^(-1/dim)/2;
r1 = sqrt(2);
r2 = sqrt(5);
threshold = .7;   % domain choice threshold (used for strictly positive density)
if ~exist('Output')
    mkdir Output;
end
fileID = fopen('./Output/console.txt','w');
%% populate vertices of the unit cube 
cube_vectors = zeros(dim, oct);                                                                
for i=1:dim
    len = 2 ^ (dim-i);
    for j=0:2^i-1                       
        cube_vectors(i, j*len+1:(j+1)*len) =  mod(j,2)/N;
    end
end

%% node data array and masks
nodes = zeros(dim, max_nodes_per_box, N^dim); 
node_indices = true(max_nodes_per_box, N^dim);                                          % logical-styled integer array to be used as mask
box_indices = false (N^dim,1);                                                           % boxes on the boundary
corner = -ones(dim,1);
           



%% main parfor                                                
tic

parfor i=1:N^dim
    corner = [rem((i-1), N);  floor(rem(i-1, N^2)/N);  floor((i-1)/N^2)]/N  ;               % TODO: this is not dimension-independent
    l = bsxfun(@plus, corner, cube_vectors);
    box_indices(i) =  (oct-sum(in_domain(l(1,:), l(2,:), l(3,:))) )  ;                 % 
    box = zeros(dim, max_nodes_per_box);    
    for j=1:max_nodes_per_box
        box(:,j) = cube_shrink * [j/max_nodes_per_box;  mod(r1*j,1);  mod(r2*j,1)]/N;       % the box is shrunk to account for inwards corner % TODO: can this be vectorized?
        box(:,j) = box(:,j) + corner + delta;
    end    
    nodes(:,:,i) = box;   
    
                                                            
end
toc

%% remove nodes outside the density support
bdry_nodes = nodes(:, :, box_indices);                                                  % coordinates of the nodes that belong to boundary boxes
flat = reshape(bdry_nodes, dim, []);                                       % 
f_vals = ~in_domain(flat(1,:), flat(2,:), flat(3,:));                                                        % values of the density function at those nodes
bdry_removed_indices = reshape( f_vals, max_nodes_per_box, []);             % indices of nodes that have to be removed; hard-coded threshold condition
node_indices(:, box_indices) = node_indices(:, box_indices) .* ~bdry_removed_indices;   % indices in the global nodeset
node_indices_flat = reshape(node_indices, 1, []);  
nodes = reshape (nodes, dim, []);                                                       % all nodes in a single array;      dim x num possible nodes
cnf = nodes(:, node_indices_flat);                                                      % after removing nodes with zero density
                             
%% node stats
outtemp = length(cnf);
fprintf( fileID, '\nNumber of nodes:      %d\n',  outtemp);
fprintf( '\nNumber of nodes:      %d\n',  outtemp)
outtemp = mean(sum(node_indices,1) );
fprintf( fileID, 'Mean number of nodes per box:      %d\n',  outtemp);
fprintf( 'Mean number of nodes per box:      %d\n',  outtemp )
outtemp = max(sum(node_indices,1) );
fprintf( fileID, 'Max number of nodes per box:      %d\n', outtemp );
fprintf( 'Max number of nodes per box:      %d\n', outtemp )
outtemp = min(sum(node_indices,1) );
fprintf( fileID, 'Min number of nodes per box:      %d\n', outtemp );
fprintf( 'Min number of nodes per box:      %d\n', outtemp )
toc
fprintf('\n')
clf;
clf;
% pbaspect([1 1 1])
% view([1 1 0])
% figure(2);
% plot3(cnf(1,:), cnf(2,:), cnf(3,:),  '.k');
 
%% repel and save nodes
fprintf( fileID, 'Performing %d repel steps.\n',  repel_steps);
fprintf( 'Performing %d repel steps.\n',  repel_steps)
cnf = repel(cnf, k_value, repel_steps, repel_power, fileID);
toc 

pbaspect([1 1 1])
% view([1 1 0])
F = figure(1);
plot3(cnf(1,:), cnf(2,:), cnf(3,:),  '.k');

savefig(F,'./Output/nodes','compact')
% save('slanttrui.mat', 'cnf')
dlmwrite('./Output/cnf.txt',cnf','delimiter','\t'); % ,'precision',3
