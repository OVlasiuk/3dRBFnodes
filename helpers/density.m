function d = density(varargin)
%DENSITY
% d = density(coords) OR d = density(x,y,z)
% Return a density that consists of two Gaussians.
% Takes matrices of size dim x N, where N is the number of vectors.
v = cell2mat(varargin);
% d = ones(1,size(v,2)) * .12;
d = min(exp( sum((v+ones(size(v))).^2,1) - 5  ),...
    exp( sum(( v ).^2,1) - 4  )) ;