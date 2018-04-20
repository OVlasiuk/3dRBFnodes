function cnf = node_earth(densityF,in_domainF)
%NODE_EARTH
% cnf = node_earth(densityF,in_domainF)
% Places a uniform distribution in a thin layer about the Earth surface,
% using the ETOPO1 data and the same irrational-lattices-based approach as
% node_dis 
% Output:
% cnf -- 3x(#of points)-sized matrix containing point coordinates;
%   Additionally produces 
% Input:
% densityF -- handle to the density function, accepts an array of size 
%   (dim)x(#of points); (currently only dim=3);
% in_domainF -- handle to the point inclusion function, accepts three arrays
%   of corrdinates, [is, radii] = in_domainF(x, y, z); returns a logical array
%   "is" of the same size as x, and array "radii" containing the
%   interpolated ETOPO radii for the spherical angles defined by input
%   values.
% 
%   See also RUNME, NODE_DIS.

% % % % % % % % % MAIN SCRIPT FOR NODES IN GEO-SETTING % % % % % % %
%% % % % % % % % % % % % PARAMETERS  % % % % % % % % % % % % % % % % % % %
N = 80;                         % number of boxes per side of the cube
maxNodesPerBox = 20;          % 
A = 2.4;                        % the bigger cube side length
dim = 3;                        % ATTN: the subsequent code is NOT dimension-independent
repelPower = 4;
oct = 2^dim;
cubeShrink = 1 - maxNodesPerBox^(-1/dim)/8;
delta = (1-cubeShrink)/2;
r1 = sqrt(2);
r2 = (sqrt(5)-1)/(sqrt(2));
adjacency = (dim+1)*2^dim;              % the number of nearest boxes to consider


close all;
s = char(mfilename('fullpath'));
cd(s(1:end-10))                         % cd to the mfile folder; 
                                        % The constant 12 depends on the
                                        % length of the filename.
addpath helpers/                                                                                
if ~exist('output','dir')
    mkdir output;
end
if ~exist('densityF','var')
    densityF=@density_earth;
end
if ~exist('in_domainF','var')
    in_domainF = @in_domain;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
try 
    load('./output/unit_lattice_radius.mat')
    if (cubeShrink ~= CUBE_SHRINK)...
            || (r1 ~= R1)...
            || (r2 ~= R2)
        throw(MException('ReadTable:NoFile','I could not find the table of radii.'));
    end
catch
    fprintf('\nLooks like the interpolation table for this number of lattice\n');
    fprintf('nodes is missing or not up to date... Hang on there, I''ll make\n'); 
    fprintf('a new one for you. This may take a few minutes, but we''ll only\n');
    fprintf('do it once.\n');
    lattice_by_count(maxNodesPerBox,cubeShrink,r1,r2,'y');
    fprintf('...\nDone.\n\n')
end

%% Populate vertices of the unit cube 
cube_vectors = zeros(dim, oct);                                                                
for i=1:dim
    len = 2 ^ (dim-i);
    for j=0:2^i-1                       
        cube_vectors(i, j*len+1:(j+1)*len) =  mod(j,2)/N;
    end
end

%% Populate a sample voxel
% This example uses uniform distribution, so we may just as well
% preallocate one voxel and reuse it for all corners.
j = 1:maxNodesPerBox;
voxel = A*cubeShrink * [mod(j/maxNodesPerBox + .5,1);  mod(r1*j,1);  mod(r2*j,1)]/N;       
voxel = voxel +  A*delta/N;


%% Main                                               
tic
I=1:N^dim;
corners = A*[rem((I-1), N);  floor(rem(I-1, N^2)/N);  floor((I-1)/N^2)]/N-A/2.0;
if ~exist('in_domainF','var')
    cornersUsed = corners;
else
    [corners_bool, ~] = in_domainF(corners(1,:), corners(2,:),  corners(3,:) );
    [IDX, ~] = knnsearch(corners', corners'+A/2/N, 'k', adjacency); 
    cornerIndices = logical(sum(corners_bool(IDX),2));
    cornersUsed = corners(:,cornerIndices);
end
count = size(cornersUsed,2);
nodes = reshape(repmat(cornersUsed,maxNodesPerBox,1), dim,[]) + ...
    reshape(repmat(voxel,1,count), dim, []);                                  

%% Remove nodes outside the density support                                                 
if ~exist('in_domainF','var')
    cnf = nodes;
else
    [f_vals, ~] = in_domainF(nodes(1,:), nodes(2,:), nodes(3,:));            
    % values of the density function
    cnf = nodes(:, f_vals);                                                      
    % after removing nodes with zero density
end
                         
%% Node stats
outtemp = length(cnf);
fprintf( '\nNumber of nodes:      %d\n',  outtemp)
toc
fprintf('\n')
clf;
 
%% Repel and save nodes
kValue =  ceil(max(log10( size(cnf,2) ).^2-5,12));
repelSteps = ceil(max(.9*log10( size(cnf,2) ).^2-5,10));
fprintf( 'Performing %d repel steps using %d nearest neighbors.\n',  repelSteps, kValue)
if ~exist('in_domainF','var')
    in_domainF = 0;
end
cnf = repel(cnf, size(cnf,2), kValue, repelSteps, A, in_domainF, 's', repelPower,'histogram',true);
 
%% Plot the results
pbaspect([1 1 1])
daspect([1 1 1])
figure(1);
msize = ceil(max(1, 22-5*log10(size(cnf,2)) ));
plot3(cnf(1,:), cnf(2,:), cnf(3,:),'.k','MarkerSize',msize);
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',12)
grid on;
axis vis3d

% collectfigs

%% Plotting and diagnostic
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
% dlmwrite('./output/cnf_earth.txt',cnf','delimiter','\t','precision',10); %
