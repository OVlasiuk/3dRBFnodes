function cnf = node_dis(densityF,in_domainF)
%NODE_DIS
% cnf = node_dis(densityF,in_domainF)
% Distributes nodes with the variable density (locally defining the
% distance to the nearest neighbor) given by the handle densityF.
% densityF -- handle to the density function, accepts an array of size 
%   (dim)x(#of points); (currently only dim=3);
% in_domainF -- handle to the point inclusion function, accepts three arrays
% of corrdinates, in_domainF(x,y,z); returns a logical array of the same
% size as x;
% 
% See also RUNME, NODE_EARTH.

% % % % % % % % % MAIN SCRIPT FOR NODE SETTING: VARIABLE DENSITY % % % % % % %
%% % % % % % % % % % % % PARAMETERS  % % % % % % % % % % % % % % % % % % %


N = 100;                         % number of boxes per side of the cube
maxNodesPerBox = 40;
A = 12;                          
jitter = 0;                     % The amount of jitter to add to the repel procedure.
dim = 3;                        % ATTN: the subsequent code is NOT dimension-independent
repelPower = 5;                 
oct = 2^dim;
cubeShrink = 1 - maxNodesPerBox^(-1/dim)/64;
delta = (1-cubeShrink)/2;
r1 = sqrt(2);
r2 = (sqrt(5)-1)/(sqrt(2));
adjacency = (dim+1)*2^dim;              % the number of nearest boxes to consider

close all;
s = char(mfilename('fullpath'));
cd(s(1:end-8))                         % cd to the mfile folder; 
                                        % The constant 12 depends on the
                                        % length of the filename.
addpath helpers/                                        
if ~exist('Output','dir')
    mkdir Output;
end
if ~exist('densityF','var')
    densityF=@density;
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
    lattice_by_count(2*maxNodesPerBox,delta,cubeShrink,r1,r2,'y');
    fprintf('...\nDone.\n\n')
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Populate vertices of the unit cube 
cubeVectors = zeros(dim, oct);                                                                
for i=1:dim
    len = 2 ^ (dim-i);
    for j=0:2^i-1                       
        cubeVectors(i, j*len+1:(j+1)*len) =  A*mod(j,2)/N;
    end
end         

%% Main                                               
tic
I=1:N^dim;
corners = A*[rem((I-1), N);  floor(rem(I-1, N^2)/N);  floor((I-1)/N^2)]/N-A/2.0;
if ~exist('in_domainF','var')
    cornersUsed = corners;
else
    [IDX, ~] = knnsearch(corners', corners'+A/2/N, 'k', adjacency);
    corners_bool = in_domainF(corners(1,:), corners(2,:),  corners(3,:) );    
    cornerIndices = logical(sum(corners_bool(IDX),2));
    cornersUsed = corners(:,cornerIndices);
end

Density = densityF(reshape(bsxfun(@plus,cubeVectors(:),repmat(cornersUsed,oct,1) ),dim,[]));
cornersAveragedDensity = mean(reshape(Density,oct,[]),1);        
cornersRadii = cornersAveragedDensity;
currentNumNodes = num_radius(cornersRadii*N/A);
% % % Note that the maximum number of nodes is capped in the following
% line: not doing it can cause MANY nodes in MANY boxes.
currentNumNodes = min([currentNumNodes; maxNodesPerBox*ones(1,numel(currentNumNodes))],[],1);
nodes = zeros(dim,sum(currentNumNodes));
previousNodes = [0 cumsum(currentNumNodes)];
    
for i=1:N^dim
    J = 1:currentNumNodes(i);
    box = A * cubeShrink * [J/currentNumNodes(i);  mod(r1*J,1);  mod(r2*J,1)]/N;
    box = bsxfun(@plus, corners(:,i)+A*delta/N, box);
     
    nodes(:,previousNodes(i)+1:previousNodes(i+1)) = box;   
end
toc

%% Remove nodes outside the density support
if ~exist('in_domainF','var')
    cnf = nodes;
else
    f_vals = in_domainF(nodes(1,:), nodes(2,:), nodes(3,:));            
    % values of the density function
    cnf = nodes(:, f_vals);                                                      
    % after removing nodes with zero density
end

%% Node stats
outtemp = length(cnf);
fprintf( '\nNumber of nodes:      %d\n',  outtemp)
fprintf( 'Mean number of nodes per box:      %d\n', mean(currentNumNodes ))
fprintf( 'Max number of nodes per box:      %d\n', max(currentNumNodes ))
fprintf( 'Min number of nodes per box:      %d\n', min(currentNumNodes ))
toc
fprintf('\n');
 
%% Repel and save nodes
kValue =  ceil(max(log10( size(cnf,2) ).^2-5,12));
repelSteps = ceil(max(.9*log10( size(cnf,2) ).^2-5,10));
fprintf( 'Performing %d repel steps using %d nearest neighbors.\n',  repelSteps, kValue)
if ~exist('in_domainF','var')
    in_domainF = 0;
end
cnf = repel(cnf,kValue,repelSteps,A,in_domainF,densityF,jitter,repelPower,0);

%% Plot the results
pbaspect([1 1 1])
figure(1);
msize = ceil(max(1, 22-5*log10(size(cnf,2)) ));
plot3(cnf(1,:), cnf(2,:), cnf(3,:),'.k','MarkerSize',msize);
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',12)
az = -101.5;
el = 30;
view(az,el);
grid on;
axis vis3d

figure(3);
[~, D] = knnsearch(cnf', cnf', 'k', 2);
rdens_cnf = D(:,2);
rdens_fun = densityF(cnf);
ratio = rdens_fun./rdens_cnf';
plot(ratio);
set(gca,'FontSize',12)
xlabel('Node {\bf\it{N}}','FontSize',24);
ylabel('\rho({\bf\it{N}})/\Delta({\bf\it{N}})','FontSize',24);
minratio = min(ratio)
maxratio = max(ratio)
meanratio = mean(ratio)
varratio = var(ratio)

% dlmwrite('./Output/cnf.txt',cnf','delimiter','\t','precision',10); % 