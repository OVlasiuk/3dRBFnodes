function rtable = lattice_by_count(COUNT, DELTA, CUBE_SHRINK, SAVE)
%LATTICE_BY_COUNT 
% rtable = lattice_by_count(COUNT, DELTA, CUBE_SHRINK, SAVE:y/n)
% Given the maximal number of nodes COUNT, the centering term and the box
% shrinking factor, construct irrational lattices with N nodes and
% parameters sqrt(2) and sqrt(5) for each N=1:COUNT, and determine the
% separation distances. The number of boxes used is 27=3^3.

rtable = zeros(1,COUNT);

N = 3;
dim = 3;        % works only in 3d
r1 = sqrt(2);
r2 = sqrt(5);

I=1:N^dim;
corners = [rem((I-1), N);  floor(rem(I-1, N^2)/N);  floor((I-1)/N^2)];

parfor nodes_per_box=1:COUNT
    box = zeros(dim, nodes_per_box);    
    for j=1:nodes_per_box
        box(:,j) = CUBE_SHRINK * [j/nodes_per_box;  mod(r1*j,1);  mod(r2*j,1)];
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
    save('./Output/unit_lattice_radius.mat', 'rtable','DELTA', 'CUBE_SHRINK')
end