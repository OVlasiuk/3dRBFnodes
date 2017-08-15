function pt_analyzer(cnf, in_domainF, figfactor, varargin)
%PT_ANALYZER
% pt_analyzer(cnf, in_domainF, figfactor, varargin)
% Given a (dim)x(number of points)-array, will determine its separation
% distance, distribution of distances to the nearest neighbors, and radii
% of the largest holes. Provided an indicator function of the restricting
% domain, in_domainF, will detect a subset of nodes, that are at most
% (separation_distance/2) away from the boundary in l1-metric. 
% 
% in_domainF is expected as 'in = in_domainF(x, y, z)',
%   where 'in' is a logical array of the same size as 'x', and 'x', 'y', 'z'
%   are the respective coordinates.
% figfactor -- the default figure numbers will be shifted by
%   100*(figfactor-1); must be an integer.
% Pass the string 'holes' as the last argument to perform hole analysis.
% 
%    See also KNNSEARCH, VORONOIN.

path_old = pwd;
path_new = char(mfilename('fullpath'));
cd(path_new(1:end-11))
addpath .
whether_holes = 0;
for i=1:numel(varargin)-1
    if (varargin{i} == 'holes') 
        whether_holes = 1;
        break
    end
end
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if length(size(cnf)) ~= 2
    fprintf( 'Sorry, can only parse a (dim)x(number of points) array.\n')
    return
end
[dim, N] = size(cnf);
if ~exist('figfactor','var') || (figfactor <= 0)
    figfactor = 1;
end

format long;
% bins = 200;
% binwidth = .0002;
adjacency = 13;
    colors = [[138,186,195]
    [86,54,41]
    [169,191,160]
    [78,55,75]
    [205,181,157]
    [37,63,81]
    [207,175,195]
    [38,66,52]
    [162,169,199]
    [65,63,35]
    [192,146,149]
    [57,110,125]
    [159,136,113]
    [88,96,123]
    [123,144,114]
    [135,118,142]
    [70,101,86]
    [125,88,90]
    [104,146,143]
    [110,99,73]]/256;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
fprintf( 'Analyzing a set of cardinality %d in %d-dimensional space.\n', N, dim);
% % % % % % % % % % SEPARATION OF THE WHOLE NODE SET % % % % % % % % % % %
[~, Dcnf] = knnsearch(cnf', cnf', 'k', adjacency+1);
Dcnf = Dcnf(:,2:end);     % the first column contains only zeros
separation_all = min(Dcnf(:,1));
fprintf( 'The separation distance of this node set:\t%3.6f\n',...
        separation_all);
fprintf( '.25 separation distance quantile:\t%3.6f\n',...
    quantile(Dcnf(:,1),.25));
% % % % % % % % % % SEPARATION OF THE SURFACE NODE SET % % % % % % % % % % 
if exist('in_domainF', 'var') && isa(in_domainF,'function_handle')
    CNF = repmat(cnf,2*dim,1);
    e=[eye(dim) -eye(dim)];
    shifted=bsxfun(@plus, separation_all*e(:)/2,CNF);
    shifted = reshape(shifted,3,[]);
    indices = ~in_domainF( shifted(1,:), shifted(2,:), shifted(3,:));
    indices = reshape(indices,6,[]);
    II=sum(indices,1);
    I=logical(II);
    fprintf( 'The number of nodes l1-close to the surface of the domain:\t%d\n',...
        sum(I));
    cnfsurf=cnf(:,I);          
    [~, Dsurf] = knnsearch(cnfsurf', cnfsurf', 'k', adjacency+1);
    Dsurf = Dsurf(:,2:end);     % the first column contains only zeros
else
     in_domainF = 1;
end
% % % % % % % % % % % % % % % % % % HOLE RADII % % % % % % % % % % % % % %
if whether_holes
    [V,~] = voronoin(cnf');
    if isa(in_domainF,'function_handle')
        V = V(in_domainF(V(:,1),V(:,2),V(:,3)),:);       % Voronoi centers inside the shell
    else
        fprintf('No domain indicator function provided, proceeding anyway...\n');
    end
    [~, holedists] = knnsearch(cnf',V);
    fprintf('Number of holes:\t %d \n',size(holedists,1));    
    fprintf('The deepest hole:\t %3.6f \n',max(holedists(:,1)));
    fprintf('Average hole depth:\t %3.6f \n',mean(holedists(:,1)));
    fprintf('.25 hole radius quantile:\t %3.6f \n',quantile(holedists(:,1), .25 ));
%     fprintf('The largest increase from 1st to %d-th hole radius is:\t %3.6f\n'...
%         ,dim+1, max(abs(holedists(:,1)-holedists(:,dim+1))) );
%     fprintf('(should be zero up to roundoff error, if a domain function is used)\n');
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % PLOTTING % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                                    figure(20 + 100* (figfactor - 1 ));    
msize = ceil(max(1, 22-3.5*log10(size(cnf,2)) ));
plot3(cnf(1,:), cnf(2,:), cnf(3,:),  '.b','MarkerSize',msize);
axis vis3d;
pbaspect([1 1 1])
daspect([1 1 1])
title('The nodes (all)')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if isa(in_domainF,'function_handle')
                                figure(30+ 100* (figfactor - 1 ));
    msize = ceil(max(1, 22-3.5*log10(size(cnfsurf,2)) ));
    plot3(cnfsurf(1,:), cnfsurf(2,:), cnfsurf(3,:),  '.b','MarkerSize',msize);
    axis vis3d;
    pbaspect([1 1 1])
    daspect([1 1 1])
    title('Surface nodes')
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Solid color denotes nearest neighbors; contour only denotes surface
% nearest neighbors.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if whether_holes
                                    figure(40+ 100* (figfactor - 1 ));
    pbaspect([1 1 1])
    daspect([1 1 1])
    hold on;
    for i=1:adjacency
        cnf_nearest_neighbors = histogram(Dcnf(:,i));%,bins,'BinWidth',binwidth
        cnf_nearest_neighbors.EdgeAlpha=0;
        cnf_nearest_neighbors.FaceAlpha = .4;
        cnf_nearest_neighbors.Normalization = 'probability';
        cnf_nearest_neighbors.EdgeColor = colors(i,:);
        if isa(in_domainF,'function_handle')
            cnfsurf_nearest_neighbors = histogram(Dsurf(:,i)); %,bins,'BinWidth',binwidth
            cnfsurf_nearest_neighbors.EdgeAlpha=1;   
            cnfsurf_nearest_neighbors.FaceAlpha = .1;
            cnfsurf_nearest_neighbors.Normalization = 'probability';
            cnfsurf_nearest_neighbors.DisplayStyle = 'stairs';  
            cnfsurf_nearest_neighbors.EdgeColor = cnf_nearest_neighbors.EdgeColor;
        end
    end
    deepest_holes = histogram(holedists(:,1)); %,bins,'BinWidth',binwidth
    deepest_holes.EdgeAlpha=1;
    deepest_holes.DisplayStyle='stairs';
    deepest_holes.Normalization = 'probability';
    deepest_holes.EdgeColor = 'black';
    deepest_holes.LineStyle=':';
    deepest_holes.LineWidth=2.0;
    
%             figure(60+ 100* (figfactor - 1 ));
% plot3(cnfsurf(1,:), cnfsurf(2,:), cnfsurf(3,:),  'or','MarkerSize',ceil(msize/2));            
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % RATIOS OF NNS % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
else
                                    figure(40+ 100* (figfactor - 1 ));
    pbaspect([1 1 1])
    daspect([1 1 1])
    hold on;
    for i=2:adjacency
        fprintf('Mean ratio of %d-th neighbor to 1st:\t%3.6f\n',i,mean(Dcnf(:,i) ./ Dcnf(:,1)));
        cnf_nearest_neighbors_ratio = histogram( Dcnf(:,i) ./ Dcnf(:,1) );%,bins,'BinWidth',binwidth
        cnf_nearest_neighbors_ratio.EdgeAlpha=0;
        cnf_nearest_neighbors_ratio.FaceAlpha = .4;
        cnf_nearest_neighbors_ratio.Normalization = 'probability';
        cnf_nearest_neighbors_ratio.EdgeColor = colors(i,:);
        if isa(in_domainF,'function_handle')
            cnfsurf_nearest_neighbors_ratio = histogram( Dsurf(:,i) ./ Dsurf(:,1) ); %,bins,'BinWidth',binwidth
            cnfsurf_nearest_neighbors_ratio.EdgeAlpha=1;   
            cnfsurf_nearest_neighbors_ratio.FaceAlpha = .1;
            cnfsurf_nearest_neighbors_ratio.Normalization = 'probability';
            cnfsurf_nearest_neighbors_ratio.DisplayStyle = 'stairs';  
            cnfsurf_nearest_neighbors_ratio.EdgeColor = cnf_nearest_neighbors_ratio.EdgeColor;
        end
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
set(gca,'FontSize',12)
ylabel('Probability of the characteristic','FontSize',24);
if whether_holes
    xlabel('Distances to the nearest neighbors vs hole radii','FontSize',24);
else
    xlabel('Distances to the nearest neighbors: ratios','FontSize',24);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% BOXPLOT:
% Input data is specified as a numeric vector or numeric matrix. If x is a
% vector, boxplot plots one box. If x is a matrix, boxplot plots one box 
% for each column of x.
% On each box, the central mark indicates the median, and the bottom and top
% edges of the box indicate the 25th and 75th percentiles, respectively. 
% The whiskers extend to the most extreme data points not considered outliers,
% and the outliers are plotted individually using the '+' symbol.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                                figure(50+ 100* (figfactor - 1 ));
if isa(in_domainF,'function_handle')
    Ry = max(max(max(Dsurf)), max(max(Dsurf)));
    subplot(1,2,1)
    hold on;
    boxplot(Dsurf) % ,'PlotStyle','compact'
    set(gca,'FontSize',12)
    figure(50 + 100* (figfactor - 1 ))
    plot(max(Dsurf,[],1));
    plot(min(Dsurf,[],1));
    leg = legend('Maximal distances for the surface node set','Minimal distances for the surface node set');
    leg.FontSize = 16;
    leg.Location = 'southeast';
    xlim([1 adjacency]);
    ylim([0 Ry]);
    %
    figure(50 + 100* (figfactor - 1 ))
    subplot(1,2,2)
    set(gca,'FontSize',12)
    boxplot(Dcnf)
    hold on
    plot(max(Dcnf,[],1));
    plot(min(Dcnf,[],1));
    leg = legend('Maximal distances for the whole node set','Minimal distances for the whole node set');
    leg.FontSize = 16;
    leg.Location = 'southeast';
    xlim([1 adjacency]);
    ylim([0 Ry]);
    figure(50 + 100* (figfactor - 1 ))
else
    hold on;
    set(gca,'FontSize',12)
    boxplot(Dcnf)
    plot(max(Dcnf,[],1));
    plot(min(Dcnf,[],1));
    leg = legend('Maximal distances for the whole node set','Minimal distances for the whole node set');
    leg.FontSize = 16;
    leg.Location = 'southeast';
    xlim([1 adjacency]);
    ylim([0 max(max(Dcnf))]);
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
cd(path_old)
