%PAPER_APPENDIX
% A script that will recreate the figures, contained in the Appendix
% section, in particular, show the comparison of irrational lattices,
% generated using the following pairs of parameters:
%     sqrt(2)             (sqrt(5)-1)/(sqrt(2));
%     sqrt(3)             sqrt(5);
%     0.785477746971993   0.348797903953764;
%     0.179373654819913   0.531793804909494;
%     0.429397337085805   0.257282784769860;
%     rand(1)             rand(1),
% as well as Riesz periodic minimizers.
% 
% See also LATTICE_BY_COUNT, NUM_RADIUS.
s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-14))
load('../Output/lattice_riesz.mat')
if exist('../Output/mrtable_riesz.mat','file')
    load('../Output/mrtable_riesz.mat');
else
    mtable = ones(1,numel(ltable));
    rtable = ones(1,numel(ltable));
    for N=2:numel(ltable)
        [~, D] = knnsearch(ltable{N}', ltable{N}','k',2);
        mtable(N) = mean(D(:,2));
        rtable(N) = min(D(:,2));
    end
    save('../Output/mrtable_riesz.mat', 'mtable','rtable')
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
dim = 3;
maxNodesPerBox = 80; 
cubeShrink = 1; % - maxNodesPerBox^(-1/dim)/2;
delta = (1-cubeShrink)/2;
num = 200;
extent = 200;
ratios = [
    sqrt(2)             (sqrt(5)-1)/(sqrt(2));
    sqrt(3)             sqrt(5);
    0.785477746971993   0.348797903953764;
    0.179373654819913   0.531793804909494;
    0.429397337085805   0.257282784769860;
    rand(1)             rand(1);
    ];
sratios = {
    'sqrt(2)',   '(sqrt(5)-1)/(sqrt(2));',...
    'sqrt(3)',   'sqrt(5);              ',...
    '0.785477746971993',   '0.348797903953764;'    ,...
    '0.179373654819913',   '0.531793804909494;'    ,...
    '0.429397337085805',   '0.257282784769860;'    ,...
    'rand(1)',   'rand(1);'
    };
markers = [
'ok';
'xk';
'+k';
'*k';
'sk';
'dk';
'vk';
'^k';
'<k';
'>k';
'pk';
            ];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if exist('../Output/lattices_rational.mat','file')
    load('../Output/lattices_rational.mat');
else
    lattices = zeros(size(ratios,1),num);
    minlattices = zeros(size(ratios,1),num);
    parfor I=1:size(lattices,1)
        [lattices(I,:), minlattices(I,:)] = ...
            lattice_by_count(num,cubeShrink,ratios(I,1), ratios(I,2),'n');
    end
    save('../Output/lattices_rational.mat', 'lattices','minlattices')
end
[lmax, ind]=max(lattices,[],1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
close all;
f1 = figure;
f1.PaperType = 'A2';
hold on;
for i=1:size(lattices,1)
    plot(lattices(i,1:extent),markers(i),'MarkerSize',6)
end
i = i+1;
plot(mtable,markers(i),'MarkerSize',8)
ii = 1:extent;
p=plot(ii.^(-1/dim),'-.','LineWidth',3);


set(gca,'FontSize',12)
xlabel('Number of nodes in the unit cube','FontSize',20);
ylabel('Mean separation distance','FontSize',20);
% RR = ratios';
legend_ratios = string(sratios);
legend_string = cell(1,size(ratios,1)+2);
for i=1:size(ratios,1)
    legend_string{i} = legend_ratios(2*i-1)+ string(' and ') + legend_ratios(2*i);
end
legend_string{end-1} = string('Periodic Riesz optimizers');
legend_string{end} = string('n^{-1/d}');

[leg, ico] = legend(legend_string{:});
leg.FontSize = 19;
i = 1;
while isa(ico(i),'matlab.graphics.primitive.Text')
    ico(i).FontSize = 17;
    i=i+1;
end
for j=i:numel(ico)
    if string(ico(j).Marker) ~= string('none')
        ico(j).MarkerSize = 17;
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
f2 = figure;
f2.PaperType = 'A2';
hold on;
for i=1:size(lattices,1)
    plot(lattices(i,1:extent)./minlattices(i,1:extent),markers(i),'MarkerSize',6)
end
i = i+1;
plot(mtable(1:extent)./rtable(1:extent),markers(i),'MarkerSize',12)

set(gca,'FontSize',12)
xlabel('Number of nodes in the unit cube','FontSize',20);
ylabel('Mean to minimal separation ratio','FontSize',20);
% RR = ratios';

% leg=legend(legend_string{:});
% leg.FontSize = 16;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
f3 = figure;
hold on
h1 = histogram(ind,10);
h2 = histogram(ind(1:80),10);

leg=legend('Index of the lattice with the largest separation',...
            'Same, up to 80 nodes per v');
h1.FaceColor = [0 0 0.9];  
h1.FaceAlpha = .2;
h1.DisplayStyle='stairs';

h2.FaceColor = [0.9 0 0];  
h2.FaceAlpha = .3;

% disp('Type ''dbcont'' to continue, ''dbquit'' to terminate:')
% keyboard
% figure(f1)
% ii = 1:extent;
% p=plot(ii.^(-1/dim),'-','LineWidth',5);
% 
% leg.Visible='off';
% l = legend(p,'N^{-1/d}');
% [l, ico] = legend(p,'N^{-1/d}');
% ico(1).FontSize = 16;
% l.FontSize = 20;


currentNumNodes = min(100, extent);
figure;
plot3(ltable{currentNumNodes}(1,:),ltable{currentNumNodes}(2,:),...
    ltable{currentNumNodes}(3,:),'.k','MarkerSize',20)
pbaspect([1 1 1])
daspect([1 1 1])
set(gca, 'Clipping', 'off')
set(gca,'xtick',linspace(0,1,5))
set(gca,'ytick',linspace(0,1,5))
set(gca,'ztick',linspace(0,1,5))
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
zlabel('z','FontSize',20);
az = 135.7; el = 32.4; view(az,el)
grid on
axis vis3d

figure;
J = 1:currentNumNodes;
currentRatio = 4;
r1 = ratios(4,1);
r2 = ratios(4,2);
irrlat = [mod(J/currentNumNodes + .5,1);  mod(r1*J,1);  mod(r2*J,1)]; %;
plot3(irrlat(1,:), irrlat(2,:), irrlat(3,:), '.b','MarkerSize',20)
pbaspect([1 1 1])
daspect([1 1 1])
set(gca, 'Clipping', 'off')
set(gca,'xtick',linspace(0,1,5))
set(gca,'ytick',linspace(0,1,5))
set(gca,'ztick',linspace(0,1,5))
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
zlabel('z','FontSize',20);
az = 135.7; el = 32.4; view(az,el)
grid on
axis vis3d

cd(s_old)