%LATTICE_PAPER
% A script that will produce and plot an example of lattice separation 
% distances for the four pairs of parameters:
%   'sqrt(2) and (sqrt(5)-1)/(sqrt(2))',...
%   'sqrt(2) and sqrt(5)',...
%   '0.785477746971993 and 0.348797903953764',...
%   '0.179373654819913 and 0.531793804909494'
% 
% See also LATTICE_BY_COUNT, NUM_RADIUS.
dim = 3;
maxNodesPerBox = 40; 
delta = 1/(128* maxNodesPerBox^(1/dim));
cubeShrink = 1 - maxNodesPerBox^(-1/dim)/64;
num = 200;
extent = 200;
l = @(r1,r2) lattice_by_count(num,cubeShrink,r1,r2,'n');
ratios = [
    sqrt(2)             (sqrt(5)-1)/(sqrt(2));
    sqrt(3)             sqrt(5);
    0.785477746971993   0.348797903953764;
    0.179373654819913   0.531793804909494;
    0.429397337085805   0.257282784769860;
    rand(1)             rand(1);
    ];
lattices = zeros(size(ratios,1),num);
for i=1:size(lattices,1)
    lattices(i,:) = l(ratios(i,1), ratios(i,2));
end
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
close all;
figure(1);
hold on;
for i=1:size(lattices,1)
    plot(lattices(i,1:extent),markers(i),'MarkerSize',6)
end
[lmax, ind]=max(lattices,[],1);
set(gca,'FontSize',12)
xlabel('Number of nodes in a box','FontSize',20);
ylabel('Minimal separation distance','FontSize',20);
RR = ratios';
legend_ratios = string(num2str(RR(:),'%2.16f\n'))';
legend_string = cell(0);
for i=1:size(ratios,1)
    legend_string = [legend_string, legend_ratios(2*i-1)+ string(' and ')+ legend_ratios(2*i)];
end
leg=legend(legend_string{:});
leg.FontSize = 16;
figure(2)
hold on
h1 = histogram(ind);
h1.FaceColor = [0 0 0.9];  
h1.FaceAlpha = .2;
h1.DisplayStyle='stairs';
h2 = histogram(ind(1:80));
h2.FaceColor = [0.9 0 0];  
h2.FaceAlpha = .3;

