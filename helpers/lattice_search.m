addpath ./helpers
its = 1e5;
format long;
max_nodes_per_box = 40;          
num = 500;
dim = 3;                        
cubeShrink = 1 - maxNodesPerBox^(-1/dim)/64;
delta = (1-cubeShrink)/2;

l = @(r1,r2) lattice_by_count(num,cubeShrink,r1,r2,'n');
r1 = 0.785477746971993;
r2 = 0.348797903953764;
lold = l(r1,r2);
improvement = 0;
I = 0;
vG = zeros(1,2);

tic
for i=1:its
%     v = .01*rand(1,2) + [.78 .34];
    v =  rand(1,2);
    lnew = l(v(1),v(2));
    dff = sum((lnew > lold).*(lnew - lold));
    if dff > improvement
        improvement = dff
        vG = v
    end
end
toc
