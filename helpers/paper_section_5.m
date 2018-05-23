%PAPER_SECTION_5
% Reproduces the figures contained in Section 5 of the associated 
% paper; in particular:
% 
% See also LATTICE_BY_COUNT, NUM_RADIUS, LATTICE_RIESZ.
s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-15))

%% Initialize
m  = 5;                     % Power in RBFs
d  = 2;                     % Use up through degree d tems
maxK = 100;                 % Maximal number of nodes in the stencil
% Riesz
nodes_riesz = dlmread('../output/riesz1k.txt');
N = size(nodes_riesz,1);
% Halton
halton_obj = haltonset(3);           % Create Halton nodes throughout unit cube
nodes_halton = halton_obj(1:N,:);
% Cartesian
N = 1e3;
x = linspace(0, 1, ceil(N^(1/3)) );
[X,Y,Z] = meshgrid(x);
nodes_cart = [X(:),Y(:),Z(:)];

%% Build knn-tree
ktree_riesz = createns(nodes_riesz,'nsmethod','kdtree'); 
ktree_halton = createns(nodes_halton,'nsmethod','kdtree'); 
ktree_cart = createns(nodes_cart,'nsmethod','kdtree');

%% Compute weights
rng(5);                     % Specify seed for reproducible results
C = [.5 .5 .5] + randn(1,3)*5e-2; 
condition_numbers = ones(maxK, 3);

for k = 1:maxK
    idx_riesz = knnsearch(ktree_riesz, C,'k',k);
    condition_numbers(k, 1) = RBF_FD_PHS_pol_condnum(nodes_riesz(idx_riesz,:), m, d);
    %
    idx_halton = knnsearch(ktree_halton, C,'k',k);
    condition_numbers(k, 2) = RBF_FD_PHS_pol_condnum(nodes_halton(idx_halton,:), m, d);
    %
    idx_cart = knnsearch(ktree_cart, C,'k',k);
    condition_numbers(k, 3) = RBF_FD_PHS_pol_condnum(nodes_cart(idx_cart,:), m, d);
end




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
markers = [
'v ';
'* ';
's ';
'sk';
'vk';
            ];
        
legend_string = cell(1,3);
legend_string{1} = "Periodic Riesz minimizers";
legend_string{2} = "Halton nodes";
legend_string{3} = "Cartesian nodes";
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
f1 = figure;
f1.PaperType = 'A2';
hold on;

plot(1:maxK, condition_numbers(:,1), markers(1,:),'MarkerSize',6,...
    'MarkerEdgeColor', [0.6350    0.0780    0.1840])
plot(1:maxK, condition_numbers(:,2), markers(2,:),'MarkerSize',6,...
    'MarkerEdgeColor', [0    0.4470    0.7410])
plot(1:maxK, condition_numbers(:,3), markers(3,:),'MarkerSize',6,...
    'MarkerEdgeColor', [0.4660    0.6740    0.1880])

set(gca, 'YScale', 'log')
set(gca,'FontSize',12)
xlabel('Number of nearest nodes in the stencil','FontSize',20);
ylabel('Condition number of RBF-FD PHS-poly','FontSize',20);

[leg, ico] = legend(legend_string{:});
leg.FontSize = 19;
i = 1;
while isa(ico(i),'matlab.graphics.primitive.Text')
    ico(i).FontSize = 17;
    i=i+1;
end
for j=i:numel(ico)
    if string(ico(j).Marker) ~= "none"
        ico(j).MarkerSize = 17;
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

cd(s_old)
