function slicehist(x,y,z, a, b)
%% slicehist(x,y,z, a, b)
% Given a collection of points described by their coordinates x, y, z, plots 
% the histogram of those lying in the range [a,b] in z-coordinate.
%%
clf;
indices = logical((z>a).*(z<b));
sum(indices)
histogram2(x(indices),y(indices),50);