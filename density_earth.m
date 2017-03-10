function [d] = density_earth(x, y, z)
% returns density above the Earth surface that decreases with height
outer = 1.1;

[~, surface_radii] = in_domain(x,y,z);
[~, ~, radii] = cart2sph(x,y,z);

d = 100000000*(outer-radii).^2./(outer-surface_radii);