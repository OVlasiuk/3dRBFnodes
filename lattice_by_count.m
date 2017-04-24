function rtable = lattice_by_count(COUNT, DELTA, CUBE_SHRINK, R1, R2, SAVE)
%LATTICE_BY_COUNT 
% rtable = lattice_by_count(COUNT, DELTA, CUBE_SHRINK, R1, R2, SAVE)
% Given the maximal number of nodes COUNT, the centering term and the box
% shrinking factor, construct irrational lattices with N nodes and
% parameters R1, R2 for each N=1:COUNT, and determine the
% separation distances. The lattices are placed in 27=3^3 adjacent boxes.
% 
% COUNT -- the maximal number of nodes in a single box to be considered;
% DELTA -- nodes inside each box will be translated by this value;
% CUBE_SHRINK -- nodes inside each box will be scaled by this factor; must
% be < 1.
% R1, R2 -- parameters of the lattice; must be linearly independent over
% the rationals;
% SAVE -- one of 'y', 'n'.

rtable = zeros(1,COUNT);

N = 3;
dim = 3;        % works only in 3d

I=1:N^dim;
corners = [rem((I-1), N);  floor(rem(I-1, N^2)/N);  floor((I-1)/N^2)];

parfor nodes_per_box=1:COUNT
    box = zeros(dim, nodes_per_box);    
    for j=1:nodes_per_box
        box(:,j) = CUBE_SHRINK * [j/nodes_per_box;  mod(R1*j,1);  mod(R2*j,1)];
        box(:,j) = box(:,j) +  DELTA;
    end
    nodes = reshape(repmat(corners,nodes_per_box,1),dim, []) +...
                        reshape( repmat(box,1,size(corners,2)), dim, []);
    [~, D] = knnsearch(nodes', nodes','k',2);
    rtable(nodes_per_box) = min(D(:,2));
end

s = char(mfilename('fullpath'));
cd(s(1:end-16))

if SAVE == 'y'
    save('./Output/unit_lattice_radius.mat', 'rtable','DELTA', 'CUBE_SHRINK',...
    'R1', 'R2')
end