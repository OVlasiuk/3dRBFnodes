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
num = 500;
extent = 200;
l = @(r1,r2) lattice_by_count(num,delta,cubeShrink,r1,r2,'n');
ratios = [
    sqrt(2)             (sqrt(5)-1)/(sqrt(2));
    sqrt(3)             sqrt(5);
    0.785477746971993   0.348797903953764;
    0.179373654819913   0.531793804909494;
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
set(gca,'FontSize',12)
xlabel('Number of nodes in a box','FontSize',20);
ylabel('Minimal separation distance','FontSize',20);
leg=legend('sqrt(2) and (sqrt(5)-1)/(sqrt(2))',...
    'sqrt(2) and sqrt(5)',...
    '0.785477746971993 and 0.348797903953764',...
    '0.179373654819913 and 0.531793804909494');
leg.FontSize = 16;
