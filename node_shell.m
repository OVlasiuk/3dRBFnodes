function cnf = node_shell(cnf_bdry, densityF,in_domainF)
%NODE_SHELL
% cnf = node_dis(densityF,in_domainF)
% Distributes nodes with the variable density (locally defining the
% distance to the nearest neighbor) given by the handle densityF.
% densityF -- handle to the density function, accepts an array of size 
%   (dim)x(#of points); (currently only dim=3);
% in_domainF -- handle to the point inclusion function, accepts three arrays
% of coordinates, in_domainF(x,y,z); returns a logical array of the same
% size as x;
% 
% See also RUNME, NODE_EARTH.

% % % % % % % % % MAIN SCRIPT FOR NODE SETTING: VARIABLE DENSITY % % % % % % %
%% % % % % % % % % % % % PARAMETERS  % % % % % % % % % % % % % % % % % % %
N = 50;                         % number of boxes per side of the cube
maxNodesPerBox = 80;
A = 14;                          
jitter = 0;                     % The amount of jitter to add to the repel procedure.
dim = 3;                        
oct = 2^dim;
cubeShrink = 1 - maxNodesPerBox^(-1/dim);
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
if ~exist('Output','dir')
    mkdir Output;
end
if ~exist('densityF','var')
    densityF=@density_shell_uni;
end
if ~exist('in_domainF','var')
    in_domainF = @(x,y,z) in_shell(x,y,z,1+density_shell_uni(1)/4,6.4046-density_shell_uni(6.4046)/4);
%     in_domainF = @(x,y,z) in_shell(x,y,z,1,6.1062);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
try 
    load('./Output/unit_lattice_radius.mat')
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
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
if ~exist('in_domainF','var') || ~isa(in_domainF,'function_handle')
    cornersUsed = corners;
else
    IDX = knnsearch(corners', corners'+A/2/N, 'k', adjacency);
    corners_bool = in_domainF(corners(1,:), corners(2,:),  corners(3,:) );
    cornerIndices = logical(sum(corners_bool(IDX),2));
    cornersUsed = corners(:,cornerIndices);
end

Density = densityF(reshape(bsxfun(@plus,cubeVectors(:),repmat(cornersUsed,oct,1) ),dim,[]));
cornersAveragedDensity = mean(reshape(Density,oct,[]),1);
currentNumNodes = num_radius(cornersAveragedDensity*N/A);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Put centers into empty boxes
centersEmpty = cornersUsed(:, ~currentNumNodes) + A/2/N;
centersDensity = densityF(centersEmpty);
[sortedEmptyDensity, sortEmpty] = sort(centersDensity);
sortedCentersEmpty = centersEmpty(:,sortEmpty);
centersEmptyFill = false(1,size(centersEmpty,2));

% [eIDX, eD] = knnsearch(sortedCentersEmpty',sortedCentersEmpty','k',100);

for emptyIndex=1:size(sortedCentersEmpty, 2)
    distMatrix = bsxfun(@minus, sortedCentersEmpty(:,emptyIndex),...
        sortedCentersEmpty(:,centersEmptyFill) );
    if isempty(distMatrix) ||...
        (min( sqrt(sum(distMatrix.*distMatrix,1)) ) > sortedEmptyDensity(emptyIndex))
        centersEmptyFill(emptyIndex) = true;
    end
end
currentNumNodes(~currentNumNodes) = double(centersEmptyFill);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Place irrational lattices everywhere
nodes = zeros(dim,sum(currentNumNodes));
previousNodes = [0 cumsum(currentNumNodes)];
    
for i=1:size(cornersUsed,2)
    J = 1:currentNumNodes(i);
%     deltaI =  cornersAveragedDensity(i)/2;
%     cubeShrinkI = 1-2*N*deltaI/A;
    box = A * cubeShrink * [mod(.5+J/currentNumNodes(i),1);  mod(r1*J,1);  mod(r2*J,1)]/N;
    box = bsxfun(@plus, cornersUsed(:,i)+ delta*A/N , box(randperm(dim),:));    
    nodes(:,previousNodes(i)+1:previousNodes(i+1)) = box;   
end
toc

%% Remove nodes outside the density support
if ~exist('in_domainF','var') || ~isa(in_domainF,'function_handle')
    cnf = nodes;
else
    in_check = in_domainF(nodes(1,:), nodes(2,:), nodes(3,:));            
    % values of the density function
    cnf = nodes(:, in_check);                                                      
    % after removing nodes with zero density
end

%% Node stats
outtemp = length(cnf);
fprintf( '\nNumber of nodes:      %d\n',  outtemp)
fprintf( 'Boxes per side of the enclosing cube:      %d\n',  N)
fprintf( 'Mean number of nodes per box:      %d\n', mean(currentNumNodes ))
fprintf( 'Max number of nodes per box:      %d\n', max(currentNumNodes ))
fprintf( 'Min number of nodes per box:      %d\n', min(currentNumNodes ))
toc
fprintf('\n');

%% Repel and save nodes
kValue =  ceil(max(log10( size(cnf,2) ).^2-5,12));
repelSteps = 20;
fprintf( 'Performing %d repel steps using %d nearest neighbors.\n',  repelSteps, kValue)
if ~exist('in_domainF','var')
    in_domainF = 0;
end
if ~exist('cnf_bdry','var')
    cnf_bdry = [];
end
cnf = [cnf cnf_bdry];
cnf = repel(cnf,size(cnf,2)-size(cnf_bdry,2),kValue,repelSteps,A,in_domainF,jitter);

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
diff = abs(rdens_fun - rdens_cnf');
radii = sqrt( sum( cnf.*cnf,1 ) );

plot(radii,ratio,'.k', 'MarkerSize',4)
hold on;
plot(radii,diff,'.g', 'MarkerSize',4)
set(gca,'FontSize',12)
xlabel('Radius {\bf\it{N}}','FontSize',24);
ylabel('\rho({\bf\it{N}})/\Delta({\bf\it{N}})','FontSize',24);

%% Density recovery
maxdiff = max(diff);
meandiff = mean(diff);
minratio = min(ratio);
maxratio = max(ratio);
quantile10 = quantile(ratio,0.1);
quantile90 = quantile(ratio,0.9);
meanratio = mean(ratio);
varratio = var(ratio);
fprintf('\nmaxdiff\t\tmeandiff\n');
fprintf('%3.6f\t%3.6f\n',maxdiff,meandiff)
fprintf('minratio\tmaxratio\tquantile10\tquantile90\n');
fprintf('%3.6f\t%3.6f\t%3.6f\t%3.6f\t\n', minratio,maxratio,quantile10,quantile90)
fprintf('meanratio\tvarratio\n')
fprintf('%3.6f\t%3.6f\n',meanratio,varratio)
figure(4)
plot3(cnf(1,ratio>quantile90),cnf(2,ratio>quantile90),cnf(3,ratio>quantile90),'.k','MarkerSize',msize)
hold on;
plot3(cnf(1,ratio<quantile10),cnf(2,ratio<quantile10),cnf(3,ratio<quantile10),'.r','MarkerSize',msize)
axis vis3d;
% min_bad_radius = min(sqrt(sum(cnf(:,ratio>quantile75).^2,1)))

% dlmwrite('./Output/cnf.txt',cnf','delimiter','\t','precision',10); % 