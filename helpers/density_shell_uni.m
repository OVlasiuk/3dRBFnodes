function d = density_shell_uni(v)
%DENSITY
% d = density(coords) OR d = density(x,y,z)
% Return a constant density.
% Takes matrices of size dim x N, where N is the number of vectors.
% v = cell2mat(varargin);
d = .048 * sqrt(sum(v.*v, 1));