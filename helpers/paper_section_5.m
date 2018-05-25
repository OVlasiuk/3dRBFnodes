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
maxK = 200;                 % Maximal number of nodes in the stencil
% Riesz
nodes_riesz = dlmread('../output/riesz1k.txt');
N = size(nodes_riesz,1);
% Halton
halton_obj = haltonset(3);           % Create Halton nodes throughout unit cube
nodes_halton = halton_obj(1:N,:);
% Cartesian
N = 1e3;
x = 0:N^(-1/3):1;
[X,Y,Z] = meshgrid(x);
nodes_cart = [X(:),Y(:),Z(:)];

%% Build knn-tree
ktree_riesz = createns(nodes_riesz,'nsmethod','kdtree'); 
ktree_halton = createns(nodes_halton,'nsmethod','kdtree'); 
ktree_cart = createns(nodes_cart,'nsmethod','kdtree');

%% Compute weights
rng(5);                     % Specify seed for reproducible results
condition_numbers = zeros(maxK, 3);
IT = 500;
for it=1:IT
C = [.5 .5 .5] + randn(1,3)*5e-2; 
    for k = 1:maxK
        idx_riesz = knnsearch(ktree_riesz, C,'k',k);
        stencil_riesz = [C; nodes_riesz(idx_riesz,:)];
        condition_numbers(k, 1) = condition_numbers(k, 1) + RBF_FD_PHS_pol_condnum(stencil_riesz, m, d);
        %
        idx_halton = knnsearch(ktree_halton, C,'k',k);
        stencil_halton = [C; nodes_halton(idx_halton,:)];
        condition_numbers(k, 2) = condition_numbers(k, 2) + RBF_FD_PHS_pol_condnum(stencil_halton, m, d);
        %
        idx_cart = knnsearch(ktree_cart, C,'k',k);
        stencil_cart = [C; nodes_cart(idx_cart,:)];
        condition_numbers(k, 3) = condition_numbers(k, 3) + RBF_FD_PHS_pol_condnum(stencil_cart, m, d);
    end
end
condition_numbers = condition_numbers / IT;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Legend
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
%% Plotting
close all;
f1 = figure;
f1.PaperType = 'A2';
hold on;
lower_cardinality = 100;

plot(lower_cardinality:maxK, condition_numbers(lower_cardinality:end,1), markers(1,:),'MarkerSize',6,...
    'MarkerEdgeColor', [0.6350    0.0780    0.1840])
plot(lower_cardinality:maxK, condition_numbers(lower_cardinality:end,2), markers(2,:),'MarkerSize',6,...
    'MarkerEdgeColor', [0    0.4470    0.7410])
plot(lower_cardinality:maxK, condition_numbers(lower_cardinality:end,3), markers(3,:),'MarkerSize',6,...
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
