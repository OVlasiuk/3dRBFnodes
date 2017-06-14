function [is, radii] = in_shell(x, y, z, r, R)
%IN_SHELL
% [is, radii] = in_shell(x, y, z, r, R)
% r, R -- inner/outer shell radii;
% is -- logical array containing 'true' at the indices for which the
% corresponding (x,y,z)-point is inside the shell;

radii = sqrt(x.*x + y.*y + z.*z);
is = min(max(radii, r), R) == radii;