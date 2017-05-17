%LATTICE_SEARCH
% Will perform a Monte Carlo search for a pair of ratios maximizing
% l1-distance above a separation distance sequence determined by r1, r2.
% Prints found values to console.
% At the time of writing has to be run with the working directory
% containing both this m-file, and lattice_by_count.m
% 
%   See also LATTICE_BY_COUNT, NUM_RADIUS.
r1 = 0.785477746971993;
r2 = 0.348797903953764;

addpath .
its = 1e5;
format long;
max_nodes_per_box = 40;          
num = 500;
dim = 3;                        
cubeShrink = 1 - maxNodesPerBox^(-1/dim)/64;
delta = (1-cubeShrink)/2;
l = @(r1,r2) lattice_by_count(num,cubeShrink,r1,r2,'n');

lold = l(r1,r2);
improvement = 0;
vG = zeros(1,2);

tic
for i=1:its
%   v = .01*rand(1,2) + [.78 .34];  % promising neighborhood 
    v =  rand(1,2);
    lnew = l(v(1),v(2));
    dff = sum((lnew > lold).*(lnew - lold));
    if dff > improvement
        improvement = dff
        vG = v
    end
end
toc
