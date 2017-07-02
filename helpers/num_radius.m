function  num = num_radius(r)
%NUM_RADIUS 
% num = num_radius(r)
% Returns the number of nodes in 3-dimensional unit cube, drawn from an
% irrational lattice, such that the minimal separation of nodes between 27
% such stacked cubes is approximately r (multiple copies are to account for
% distances between neighboring cubes).
% For specifics of how the nodes are picked from the lattice, and for the
% scaling/shift we apply in each cube:
%   see also LATTICE_BY_COUNT, NODE_DIS.
s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-10))

persistent mtable;
if isempty(mtable)
    load('../Output/unit_lattice_radius.mat');    
end
num = interp1(mtable,1:numel(mtable),r,'pchip');
%  We assume that the neighbor boxes have at least one node in them
num = round(num.*(r<=1.0));
num(r<mtable(end)) = size(mtable,2);
% num((1+sqrt(3)/2<r) & (r<=2+1+sqrt(3)/2)) = 1;

% ind = INDEX(uint16(max(num,1)));
cd(s_old)