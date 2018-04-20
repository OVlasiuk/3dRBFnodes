
function [mtable, rtable] = lattice_by_count(COUNT,CUBE_SHRINK, R1, R2, SAVE)
%LATTICE_BY_COUNT 
% rtable = lattice_by_count(COUNT, CUBE_SHRINK, R1, R2, SAVE)
% Given the maximal number of nodes COUNT, the centering term and the box
% shrinking factor, construct irrational lattices with N nodes and
% parameters R1, R2 for each N=1:COUNT, and determine the
% separation distances. The lattices are placed in 27=3^3 adjacent boxes.
% 
% COUNT -- the maximal number of nodes in a single box to be considered;
% CUBE_SHRINK -- nodes inside each box will be scaled by this factor; must
% be < 1.
% R1, R2 -- parameters of the lattice; must be linearly independent over
% the rationals;
% SAVE -- one of 'y', 'n'.
DELTA = (1 - CUBE_SHRINK)/2;
rtable = zeros(1,COUNT);
mtable = zeros(1,COUNT);

N = 3;
dim = 3;        % works only in 3d

I=1:N^dim;
corners = [rem((I-1), N);  floor(rem(I-1, N^2)/N);  floor((I-1)/N^2)];

for nodes_per_box=1:COUNT    
    j=1:nodes_per_box;
    box = CUBE_SHRINK * [mod(j/nodes_per_box+.5,1);  mod(R1*j,1);  mod(R2*j,1)];
    box = box +  DELTA;
    nodes = reshape(repmat(corners,nodes_per_box,1),dim, []) +...
                        reshape( repmat(box,1,size(corners,2)), dim, []);
    [~, D] = knnsearch(nodes', nodes','k',2);
    rtable(nodes_per_box) = min(D(:,2));
    mtable(nodes_per_box) = mean(D(:,2));
end

if SAVE == 'y'
    s_old = pwd;
    s = char(mfilename('fullpath'));
    cd(s(1:end-16))
    save('../output/unit_lattice_radius.mat', 'rtable','mtable','DELTA', 'CUBE_SHRINK',...
    'R1', 'R2')
    cd(s_old)
end
