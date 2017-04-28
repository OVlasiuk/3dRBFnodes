%NODE_EARTH
% Places a uniform distribution in a thin layer about the Earth surface,
% using the ETOPO1 data and the same irrational-lattices-based approach as
% node_dis

% TODO: dimension-agnostic code
% %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % MAIN SCRIPT FOR NODES IN GEO-SETTING % % % % % % %
N = 80;                         % number of boxes per side of the cube
max_nodes_per_box = 20;          % 
repel_steps = 30;               % the number of iterations of the repel.m routine
density_f = @density_earth;     % put the handle to your density function here
k_value = 30;                   % number of nearest neighbors used in the repel.m
A = 2.4;                        % the bigger cube side length
%%
dim = 3;                        % ATTN: the subsequent code is NOT dimension-independent
repel_power = 5;
bins = 100;
oct = 2^dim;
delta = 1/(256* N * max_nodes_per_box^(1/dim));
cubeShrink = 1 - max_nodes_per_box^(-1/dim)/64;
r1 = sqrt(2);
r2 = (sqrt(5)+1)/(sqrt(2));
adjacency = 3^dim;              % the number of nearest boxes to consider


close all;
s = char(mfilename('fullpath'));
cd(s(1:end-10))                         % cd to the mfile folder; 
                                        % The constant 12 depends on the
                                        % length of the filename.
addpath helpers/                                                                                
if ~exist('Output','dir')
    mkdir Output;
end
fileID = fopen('./Output/console.txt','w');

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
    lattice_by_count(500,delta,cubeShrink,r1,r2,'y');
    fprintf('...\nDone.\n\n')
end

%% populate vertices of the unit cube 
cube_vectors = zeros(dim, oct);                                                                
for i=1:dim
    len = 2 ^ (dim-i);
    for j=0:2^i-1                       
        cube_vectors(i, j*len+1:(j+1)*len) =  mod(j,2)/N;
    end
end

%% node data array and masks
% nodes = zeros(dim, max_nodes_per_box, N^dim); 
% node_indices = false(max_nodes_per_box, N^dim);                                          
% logical-styled integer array to be used as mask
% box_indices = false (N^dim,1);                                                           
% boxes on the boundary
corner = -ones(dim,1);
           
box = zeros(dim, max_nodes_per_box);    
for j=1:max_nodes_per_box
    box(:,j) = A*cubeShrink * [j/max_nodes_per_box;  mod(r1*j,1);  mod(r2*j,1)]/N;       
    % the box is shrunk to account for inwards corner % TODO: can this be vectorized?
    box(:,j) = box(:,j) +  delta;
end    


%% main                                               
tic
I=1:N^dim;
corners = A*[rem((I-1), N);  floor(rem(I-1, N^2)/N);  floor((I-1)/N^2)]/N-A/2.0;
[corners_bool, ~] = in_domain(corners(1,:), corners(2,:),  corners(3,:) );
[IDX, ~] = knnsearch(corners', corners', 'k', adjacency); 
corner_indices = logical(sum(corners_bool(IDX),2));
corners_used = corners(:,corner_indices);
count = size(corners_used,2);
nodes = reshape(repmat(corners_used,max_nodes_per_box,1), dim,[]) + ...
    reshape(repmat(box,1,count), dim, []);                                  
% coordinates of the actual nodes
%% remove nodes outside the density support                                                 
[f_vals, ~] = in_domain(nodes(1,:), nodes(2,:), nodes(3,:));            
% values of the density function at those nodes
cnf = nodes(:, f_vals);                                                      
% after removing nodes with zero density
                         
%% node stats
outtemp = length(cnf);
fprintf( fileID, '\nNumber of nodes:      %d\n',  outtemp);
fprintf( '\nNumber of nodes:      %d\n',  outtemp)

% F = figure(1);
% plot3(cnf(1,:), cnf(2,:), cnf(3,:),  '.k','MarkerSize',1);

toc
fprintf('\n')
clf;
clf;
% pbaspect([1 1 1])
% view([1 1 0])
% figure(2);
% plot3(cnf(1,:), cnf(2,:), cnf(3,:),  '.k');
 
% repel and save nodes
fprintf( fileID, 'Performing %d repel steps.\n',  repel_steps);
fprintf( 'Performing %d repel steps.\n',  repel_steps)
cnf = repel(cnf, k_value, repel_steps, A, @in_domain, density_f, 0, repel_power, 0);
 

pbaspect([1 1 1])
% view([1 1 0])
F = figure(1);
plot3(cnf(1,:), cnf(2,:), cnf(3,:),  '.k','MarkerSize',1);
axis vis3d

%% plotting and diagnostic
% figure(3);
% pbaspect([1 1 1])
% Y = surface_nodes(size(cnf,2));
% [surface_indices, distances] = knnsearch(cnf', Y, 'k', 8);
% surface_cnf = unique(surface_indices);
% plot3(cnf(1,surface_cnf), cnf(2,surface_cnf), cnf(3,surface_cnf),  '.r','MarkerSize',7);
% hold on;
% interior_cnf=setdiff(1:size(cnf,2),surface_cnf);
% plot3(cnf(1,interior_cnf), cnf(2,interior_cnf), cnf(3,interior_cnf),  '.k','MarkerSize',1);
% [~, knn_surface] =  knnsearch(cnf', cnf(:,surface_cnf)', 'k', 4);
% [~, knn_interior] =  knnsearch(cnf', cnf(:,interior_cnf)', 'k', 4);
% figure(4)
% h3 = histogram(knn_surface(:,2),bins);
% h3.FaceColor = [0.9 0.2 0];
% hold on;
% h4 = histogram(knn_interior(:,2),bins);
% h4.FaceColor = [0.1 0.1 0.8];
% figure(5)
% h5 = histogram(knn_surface(:,2)./knn_surface(:,4),bins);
% h5.FaceColor = [0.6 0.6 0];
% hold on;
% h6 = histogram(knn_interior(:,2)./knn_interior(:,4),bins);
% h6.FaceColor = [0.05 0.8 0];


% % % % % % % % % % % % % % % % % % % % 
% savefig(F,'./Output/nodes','compact')
% dlmwrite('./Output/cnf.txt',cnf','delimiter','\t'); % ,'precision',3)
% save('cnf.mat', 'cnf')
