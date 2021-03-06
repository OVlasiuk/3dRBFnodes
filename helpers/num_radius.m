function  num = num_radius(r, varargin)
%NUM_RADIUS 
% num = num_radius(r,TYPE)
% Returns the number of nodes in 3-dimensional unit cube, drawn from an
% irrational lattice or the periodic Riesz minimizers, such that the mean 
% separation of nodes between 27 such stacked cubes is approximately r 
% (multiple copies are to account for distances between neighboring cubes).
% r -- radius to be approximated;
% TYPE -- (optional)is one of the string parameters: 'riesz', 'irrational';
%   if no argument passed, the default 'irrational' is used.
% 
% For specifics of how the nodes are picked, and for the scaling/shift we 
% apply in each cube:
%   see also LATTICE_BY_COUNT, NODE_DIS.

s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-10))

% persistent mtable;
% if isempty(mtable)
if isempty(varargin)
    load('../output/unit_lattice_radius.mat');
elseif any(varargin{:} == 'riesz')
    load('../output/mrtable_riesz.mat');    
elseif any(varargin{:} == 'irrational')
    load('../output/unit_lattice_radius.mat');
end
% end
num = interp1(mtable,1:numel(mtable),r,'pchip');
%  We assume that the neighbor boxes have at least one node in them
num = round(num.*(r<=1.0));
num(r<mtable(end)) = size(mtable,2);
% num((1+sqrt(3)/2<r) & (r<=2+1+sqrt(3)/2)) = 1;

% ind = INDEX(uint16(max(num,1)));
cd(s_old)
