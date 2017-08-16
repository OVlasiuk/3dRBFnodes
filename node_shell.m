function cnf = node_shell(cnf_bdry, densityF,in_domainF)
%NODE_SHELL
% cnf = node_dis(densityF,in_domainF)
% Distributes nodes with the variable density (locally defining the
% distance to the nearest neighbor) given by the handle densityF: a 
% constant multiple of the norm. The nodes are confined to a spherical
% shell with radii 1 and rcapRad; the latter is determined by a coordinate
% transform from Earth scale to the working scale. In the Earth scale, the
% desired node set is such that between radii 'a' and 'ztop' there are 
% about 'Nr' layers with about 'Ns' nodes in each.
% 
% densityF -- handle to the density function, accepts an array of size 
%   (dim)x(#of points): densityF(cnf); (currently only dim=3);
% in_domainF -- handle to the point inclusion function, accepts three arrays
% of coordinates: in_domainF(x,y,z); returns a logical array of the same
% size as x;
% 
% See also RUNME, NODE_EARTH, NODE_DIS.

% % % % % % % % % MAIN SCRIPT FOR NODE SETTING: VARIABLE DENSITY % % % % % % %
%% % % % % % % % % % % % PARAMETERS  % % % % % % % % % % % % % % % % % % %
N = 65;                         % number of boxes per side of the cube
maxNodesPerBox = 60;
A = 6;          
cutoffLength = 5e3;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Specific to the atmodeling problems:
a = 6371220;
ztop = 12000;
Ns = 12100;
Nr = 30;
ksep = @(Ns) sqrt(8*pi /Ns/sqrt(3)) ;
rcap = @(r) a * exp( sqrt(8*pi/Ns/sqrt(3)) * (r-a) *(Nr-1)/ztop);
rcapRad = rcap(a+ztop) / rcap(a);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

dim = 3;                        
oct = 2^dim;

r1 = 6*pi;
r2 = exp(1)/2;
adjacency = (dim+1)*2^dim;              % the number of nearest boxes to consider
close all;
s = char(mfilename('fullpath'));
cd(s(1:end-10))                         % cd to the mfile folder; 
% The constant depends on the
% length of the filename.
addpath helpers/                                        
if ~exist('Output','dir')
    mkdir Output;
end
if ~exist('densityF','var')
    densityF=@(v)  ksep(Ns) * sqrt(sum(v.*v, 1))   ;
end
if ~exist('in_domainF','var')
    in_domainF = @(x,y,z) in_shell(x,y,z,1,rcapRad);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% MAIN
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Populate vertices of the unit cube 
load('Output/lattice_riesz.mat');
cubeVectors = zeros(dim, oct);                                                                
for i=1:dim
    len = 2 ^ (dim-i);
    for j=0:2^i-1                       
        cubeVectors(i, j*len+1:(j+1)*len) =  A*mod(j,2)/N;
    end
end         

%% Select corners in the support                                              
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
currentNumNodes = num_radius(cornersAveragedDensity*N/A, 'riesz');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Place centers into empty boxes
fillingIndices = ~currentNumNodes & (cornersAveragedDensity*N/A > 1.5);
centersEmpty = cornersUsed(:, fillingIndices) + A/2/N;
centersEmptyDensity = densityF(centersEmpty);
[sortedEmptyDensity, sortEmpty] = sort(centersEmptyDensity,'ascend');
sortedCentersEmpty = centersEmpty(:,sortEmpty);
centersEmptyFill = false(1,size(sortedCentersEmpty,2));

I = 1:size(sortedCentersEmpty, 2);
emptynum = size(sortedCentersEmpty, 2)
cutoff = 0;
dlarge = inf;
new = true;

tic
for emptyIndex=1:emptynum
    if ~any(centersEmptyFill)
        centersEmptyFill(emptyIndex) = true;
        continue
    end
    if (mod(sum(centersEmptyFill), cutoffLength) == 1) && (sum(centersEmptyFill)>1) && new
%         emptyIndex
%         sum(centersEmptyFill)
        cutoff = emptyIndex-1;
        indlarge = centersEmptyFill & [true(1,cutoff),false(1,emptynum-cutoff)];
        tic
        ns = createns(sortedCentersEmpty(:,indlarge)', 'nsmethod','kdtree');
        toc
        new = false;
    end
    indsmall = centersEmptyFill & [false(1,cutoff),true(1,emptynum-cutoff)];
    distMatrix = bsxfun(@minus, sortedCentersEmpty(:,emptyIndex),...
                                sortedCentersEmpty(:, indsmall));
    dsmall = min(sqrt(sum(distMatrix.*distMatrix,1)));           
    if cutoff
        [~, dlarge] = knnsearch(ns,sortedCentersEmpty(:,emptyIndex)');
    end
    if isempty(dsmall)
        dsmall = dlarge;
    end
    if  (min(dsmall, dlarge) > sortedEmptyDensity(emptyIndex))
        centersEmptyFill(emptyIndex) = true;
        new = true;
    end
end
I(sortEmpty) = I;
currentNumNodes(fillingIndices) = double(centersEmptyFill(I));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Place irrational lattices everywhere
cubeShrink = 1  - currentNumNodes.^(-1/dim)/2;
cubeShrink(currentNumNodes > 8) = 1  - currentNumNodes(currentNumNodes > 8).^(-1/dim)/sqrt(2);
delta = (1-cubeShrink)/2;

nodes = zeros(dim,sum(currentNumNodes));
previousNodes = [0 cumsum(currentNumNodes)];    

for i=1:size(cornersUsed,2)
    if ~currentNumNodes(i)
        continue
    end
    box = A  * cubeShrink(i)  * ltable{currentNumNodes(i)}/N;
    box = bsxfun(@plus, cornersUsed(:,i)+ delta(i) *A/N , box);    
    nodes(:,previousNodes(i)+1:previousNodes(i+1)) = box;   
end
toc

%% Remove nodes outside the density support
if ~exist('in_domainF','var') || ~isa(in_domainF,'function_handle')
    cnf = nodes;
else
    in_check = in_domainF(nodes(1,:), nodes(2,:), nodes(3,:));
    sum(~in_check)
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
kValue = 20;% ceil(max(log10( size(cnf,2) ).^2-5,12));
repelSteps = 60;
fprintf( 'Performing %d repel steps using %d nearest neighbors.\n',  repelSteps, kValue)
if ~exist('in_domainF','var')
    in_domainF = 0;
end
if ~exist('cnf_bdry','var')
    cnf_bdry = [];
end
% Right density:
rdensity = @(v)  ksep(Ns) * sqrt(sum(v.*v, 1));
dcompare(cnf, rdensity);

figure(5);
hold on;
rcnf = sqrt(sum(cnf.*cnf,1));
histogram(rcnf((rcnf>1+.0001) & (rcnf< rcapRad - .0001)), 500);
cnf = [cnf cnf_bdry];

clear in_domainF;
in_domainF = @(x,y,z) in_shell(x,y,z,1,rcapRad);
cnf = repel(cnf,size(cnf,2)-size(cnf_bdry,2),kValue,repelSteps,rdensity,in_domainF);

%% Plot the results
figure(1);
msize = ceil(max(1, 22-5*log10(size(cnf,2)) ));
plot3(cnf(1,:), cnf(2,:), cnf(3,:),'.k','MarkerSize',msize);
if ~usejava('desktop')
     print('cnf','-dpdf','-r300','-bestfit')
end
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',12)
az = -101.5;
el = 30;
view(az,el);
pbaspect([1 1 1])
daspect([1 1 1])
grid on;
axis vis3d
% 
%% Density recovery
dcompare(cnf, rdensity, 1);
figure(5);
rcnf = sqrt(sum(cnf.*cnf,1));
histogram(rcnf((rcnf>1+.0001) & (rcnf< rcapRad - .0001)), 500);
if ~usejava('desktop')
    print('radial','-dpdf','-r300','-bestfit')
end

% dlmwrite('./Output/cnf.txt',cnf','delimiter','\t','precision',10); % 
