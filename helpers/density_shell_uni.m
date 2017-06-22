function d = density_shell_uni(v)
%DENSITY
% d = density(coords) OR d = density(x,y,z)
% Return a constant density.
% Takes matrices of size dim x N, where N is the number of vectors.
% v = cell2mat(varargin);
rs = sqrt(sum(v.*v, 1));
% 0.034629556613235
% d = 0.0174597 * rs;          % 0.0346/0.0215 = 1.6093
d = 0.0215 * rs;
% 2.729878606661392
d(rs > 2.72986 | rs < 1.00001 ) = .15;