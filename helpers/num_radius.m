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

persistent rtable;
if isempty(rtable)
    load('../Output/unit_lattice_radius.mat');    
end
num = interp1(rtable,1:numel(rtable),r,'pchip');
%  We assume that the neighbor boxes have at least one node in them
num = round(num.*(r<=1.0));
num(r<rtable(end)) = size(rtable,2);

% ind = INDEX(uint16(max(num,1)));
cd(s_old)