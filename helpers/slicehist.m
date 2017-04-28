function slicehist(x,y,z, a, b)
%% slicehist(x,y,z, a, b)
% Given a collection of points described by their coordinates x, y, z, plots 
% the frequency histogram of (x,y)-coordinates of those lying in the range [a,b] in z-coordinate.
%%
bins = 100;
clf;
indices = (z>a)&&(z<b);
sum(indices)
histogram2(x(indices),y(indices),bins)
