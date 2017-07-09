%LATTICE_PAPER
% A script that will produce and plot an example of lattice separation 
% distances for the four pairs of parameters:
%   'sqrt(2) and (sqrt(5)-1)/(sqrt(2))',...
%   'sqrt(2) and sqrt(5)',...
%   '0.785477746971993 and 0.348797903953764',...
%   '0.179373654819913 and 0.531793804909494'
% 
% See also LATTICE_BY_COUNT, NUM_RADIUS.
s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-13))
if exist('../Output/mrtable_riesz.mat','file')
    load('../Output/mrtable_riesz.mat');
else
    load('../Output/lattice_riesz.mat')
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
cubeShrink = 1 - maxNodesPerBox^(-1/dim)/2;
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
        [lattices(I,:), minlattices(I,:)] = lattice_by_count(num,cubeShrink,ratios(I,1), ratios(I,2),'n');
    end
    save('../Output/lattices_rational.mat', 'lattices','minlattices')
end
[lmax, ind]=max(lattices,[],1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
close all;
f1 = figure;
f1.PaperType = 'A2';
% f1.Pa
hold on;
for i=1:size(lattices,1)
    plot(lattices(i,1:extent),markers(i),'MarkerSize',6)
end
i = i+1;
plot(mtable,markers(i),'MarkerSize',8)

set(gca,'FontSize',12)
xlabel('Number of nodes in a box','FontSize',20);
ylabel('Mean separation distance','FontSize',20);
% RR = ratios';
legend_ratios = string(sratios);
legend_string = cell(1,size(ratios,1)+1);
for i=1:numel(legend_string)-1
    legend_string{i} = legend_ratios(2*i-1)+ string(' and ')+ legend_ratios(2*i);
end
legend_string{end} = string('Periodic Riesz optimizers');
[leg, ico] = legend(legend_string{:});
leg.FontSize = 17;
i = 1;
while isa(ico(i),'matlab.graphics.primitive.Text')
    ico(i).FontSize = 16;
    i=i+1;
end
for j=i:numel(ico)
    if ico(i).Marker ~= 'none'
        ico(i).MarkerSize = 18;
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
f2 = figure;
hold on;
for i=1:size(lattices,1)
    plot(lattices(i,1:extent)./minlattices(i,1:extent),markers(i),'MarkerSize',6)
end
i = i+1;
plot(mtable(1:extent)./rtable(1:extent),markers(i),'MarkerSize',12)

set(gca,'FontSize',12)
xlabel('Number of nodes in a box','FontSize',20);
ylabel('Mean to minimal separation ratio','FontSize',20);
% RR = ratios';

leg=legend(legend_string{:});
leg.FontSize = 16;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
f3 = figure;
hold on
h1 = histogram(ind);
h2 = histogram(ind(1:80));

leg=legend('Index of the lattice with the largest separation','Same, up to 80 nodes per box');
h1.FaceColor = [0 0 0.9];  
h1.FaceAlpha = .2;
h1.DisplayStyle='stairs';

h2.FaceColor = [0.9 0 0];  
h2.FaceAlpha = .3;


cd(s_old)