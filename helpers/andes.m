%ANDES
% Extracts and plots a piece of South-American coast next to the Andes,
% using the output of node_earth.m
% Also, computes and displays some statistics about the distribution of
% distances to the first nearest neighbors (the exact number is adjusted by
% *adjacency* option below), for the whole configuration, and for the 
% surface nodes.
% 
%      See also  NODE_EARTH, IN_DOMAIN

s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-5))
addpath ..
%
cnf = node_earth;
%
format long;
bins = 200;
binwidth = .0002;
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

% % % % % % % % % % SEPARATION OF THE WHOLE NODE SET % % % % % % % % % % %
[~, Dcnf] = knnsearch(cnf', cnf', 'k', adjacency);
num_nodes_total = size(cnf,2)
separation_all = min(Dcnf(:,2))

% % % % % % % % % % SEPARATION OF THE SURFACE NODE SET % % % % % % % % % % 
CNF = repmat(cnf,6,1);
e=[eye(3) -eye(3)];
shifted=bsxfun(@plus, separation_all*e(:)/2,CNF);
shifted = reshape(shifted,3,[]);
indices = ~in_domain( shifted(1,:), shifted(2,:), shifted(3,:));
indices = reshape(indices,6,[]);
II=sum(indices,1);
I=logical(II);
num_surface_nodes = sum(I)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
rads = sum(cnf.^2,1);
inds=rads>1.19;                 % radii squared are close to 1.21 on the 
                                % outer boundary: the upper bound in in_domain 
cnfsurf=cnf(:,I&~inds);          
[~, Dsurf] = knnsearch(cnfsurf', cnfsurf', 'k', adjacency);

% % % % % % % % % % % % % % ANDES % % % % % % % % % % % % % % % % % % % % % 
hemi=cnf(1,:)>0 & cnf(2,:)<0 & cnf(3,:)<0;
cnfhemi=cnf(:,I&~inds&hemi);
rot=[1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1/sqrt(2)];
hemirot=zeros(size(cnfhemi));
hemirot(3,:)=cnfhemi(3,:);
hemirot(1:2,:) = rot*cnfhemi(1:2,:);
r= [2 3 1];
chemi = hemirot(r,:);
Ihemi = chemi(1,:)<0 & chemi(3,:)>0.4;
chemirads = sqrt(sum(chemi.^2,1));
az = -67.5;
el = 32;

% % % % % % % % % % % % % % % % % % HOLE RADII % % % % % % % % % % % % % % 
[V,~] = voronoin(cnf');
S = sqrt(sum(V.*V,2));
[~, radii] = in_domain(V(:,1),V(:,2),V(:,3));
VV = V((S<1.1)&(S>radii'),:);       % Voronoi centers inside the shell
[~, holedists] = knnsearch(cnf',VV,'k',4);
increase_from_first_to_fourth_largest_hole = max(abs(holedists(:,1)-holedists(:,4)))

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
figure(3)
plot3(cnfsurf(1,:), cnfsurf(2,:), cnfsurf(3,:),  '.b','MarkerSize',1);
title('Surface nodes')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
figure(5);
pbaspect([1 1 1])
hold on;
for i=2:adjacency
    hcnf = histogram(Dcnf(:,i),bins,'BinWidth',binwidth);
    hsurf = histogram(Dsurf(:,i),bins, 'BinWidth',binwidth);
    hcnf.EdgeAlpha=0;
    hsurf.EdgeAlpha=1;
    hcnf.FaceAlpha = .4;
    hsurf.FaceAlpha = .1;
    hcnf.Normalization = 'probability';
    hsurf.Normalization = 'probability';
%     hcnf.DisplayStyle = 'stairs';
    hsurf.DisplayStyle = 'stairs';
    hcnf.EdgeColor = colors(i-1,:);
    hsurf.EdgeColor = hcnf.EdgeColor;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
hh = histogram(holedists(:,1),bins, 'BinWidth',binwidth);
hh.EdgeAlpha=1;
hh.DisplayStyle='stairs';
hh.Normalization = 'probability';
hh.EdgeColor = 'black';
hh.LineStyle=':';
hh.LineWidth=1.5;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
set(gca,'FontSize',12)
ylabel('Number of nodes','FontSize',24);
xlabel('Distances to the nearest neighbors vs hole radii','FontSize',24);
separation_surface = min(Dsurf(:,2))

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
figure(6);
subplot(1,2,1)
hold on;
boxplot(Dsurf) % ,'PlotStyle','compact'
figure(6)
plot(max(Dsurf,[],1));
plot(min(Dsurf,[],1));
leg = legend('Maximal distances for the surface node set','Minimal distances for the surface node set');
leg.FontSize = 16;
leg.Location = 'southeast';
xlim([1 adjacency]);
ylim([0 max(max(Dsurf))]);
%
figure(6)
subplot(1,2,2)
boxplot(Dcnf)
hold on
plot(max(Dcnf,[],1));
plot(min(Dcnf,[],1));
leg = legend('Maximal distances for the whole node set','Minimal distances for the whole node set');
leg.FontSize = 16;
leg.Location = 'southeast';
xlim([1 adjacency]);
ylim([0 max(max(Dsurf))]);
figure(6)
set(gca,'FontSize',12)


% % % % % % % % % % % % % % % PLOTTING THE ANDES % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
figure(10);
scatter3(chemi(1,Ihemi), chemi(2,Ihemi), chemi(3,Ihemi),...
    ones(1,size(chemi(Ihemi),2)) * 12,chemirads(Ihemi),'filled');
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',12)
% title('The Western coast of South America with color-encoded heights')
view(az,el); axis vis3d;
% print('andes_scatter','-dpdf','-r300','-bestfit')
figure(11);
plot3(chemi(1,Ihemi), chemi(2,Ihemi), chemi(3,Ihemi),  '.b',...
    'MarkerSize',5);
% title('The Western coast of South America')
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',12)
view(az,el);grid on; axis vis3d;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% print('andes','-dpdf','-r300','-bestfit')
cd(s_old)