function cnf = node_cloud(densityF,in_domainF)
%NODE_DIS
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


N = 70;                         % number of boxes per side of the cube
maxNodesPerBox = 200;
A = 2;                          
dim = 3;                        % ATTN: the subsequent code is NOT dimension-independent
oct = 2^dim;
cubeShrink = 1 - maxNodesPerBox^(-1/dim)/3;
delta = (1-cubeShrink)/2;
 r1 = sqrt(2);
 r2 = (sqrt(5)-1)/(sqrt(2));
%r1 = 0.179373654819913;
%r2 = 0.531793804909494;
adjacency = (dim+1)*2^dim;              % the number of nearest boxes to consider
cutoffLength = 5e3;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
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
    densityF = @(x) 0.8*density_cloud(x); % empiric scale adjustment
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
%% MAIN
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Populate vertices of the unit cube 
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
currentNumNodes = num_radius(cornersAveragedDensity*N/A);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Place centers into empty boxes
fillingIndices = (cornersAveragedDensity*N/A > 1.2);

for J = 1:10
    fillingIndices = ~currentNumNodes & fillingIndices;
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
        if  (min(dsmall, dlarge) > (1 + (J-1)/20.0) * sortedEmptyDensity(emptyIndex))
            centersEmptyFill(emptyIndex) = true;
            new = true;
        end
    end
    I(sortEmpty) = I;
    currentNumNodes(fillingIndices) = double(centersEmptyFill(I));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Place irrational lattices everywhere
% % % Note that the maximum number of nodes is capped in the following
% line: not doing it can cause MANY nodes in MANY boxes.
currentNumNodes = min([currentNumNodes; maxNodesPerBox*ones(1,numel(currentNumNodes))],[],1);
nodes = zeros(dim,sum(currentNumNodes));
previousNodes = [0 cumsum(currentNumNodes)];
    
for i=1:size(cornersUsed,2)
    J = 1:currentNumNodes(i);
    box = A * cubeShrink * [mod(J/currentNumNodes(i) + .5,1);  mod(r1*J,1);  mod(r2*J,1)]/N;
    box = bsxfun(@plus, cornersUsed(:,i)+A*delta/N, box(randperm(3),:));   
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
kValue =  30;
repelSteps = 20;
fprintf( 'Performing %d repel steps using %d nearest neighbors.\n',  repelSteps, kValue)
if ~exist('in_domainF','var')
    in_domainF = 0;
end

in_domainF = 0;
plback = @(v) pback(v, A);
cnf = repel(cnf,size(cnf,2),kValue,repelSteps,densityF,in_domainF,'A',A,...
                        'pullback', plback);
close all


rep = 1;   % (quantile(r, .97) > 1) &&
while  rep < 10
    cnf = repel(cnf,size(cnf,2),kValue,repelSteps,densityF,in_domainF,'A',A,...
                        'pullback', plback);
    rep = rep + 1;
    r = dcompare(cnf,densityF);
end

%% Plot the results
figure;
msize = ceil(max(1, 22-5*log10(size(cnf,2)) ));
plot3(cnf(1,:), cnf(2,:), cnf(3,:),'.k','MarkerSize',msize);
xlabel('x')
ylabel('y')
zlabel('z')
pbaspect([1 1 1])
daspect([1 1 1])
set(gca, 'Clipping', 'off')
set(gca,'FontSize',12)
grid on;
axis vis3d


% figure(3);
% [~, D] = knnsearch(cnf', cnf', 'k', 2);
% rdens_cnf = D(:,2);
% rdens_fun = densityF(cnf);
% ratio = rdens_fun./rdens_cnf';
% diff = abs(rdens_fun - rdens_cnf');
% plot(ratio);
% hold on;
% plot(diff)

set(gca,'FontSize',12)
xlabel('Node {\bf\it{N}}','FontSize',24);
ylabel('\rho({\bf\it{N}})/\Delta({\bf\it{N}})','FontSize',24);
legend('Ratio','Difference');


dcompare(cnf,@density_cloud,1);


% dlmwrite('./Output/cnf.txt',cnf','delimiter','\t','precision',10); % 
