function [d] = density(varargin)
% Return a density that consists of two Gaussians.
% Takes matrices of size dim x N, where N is the number of vectors.
v = cell2mat(varargin);
d = .15+.8*exp(sum(v.^2,1));%+.2*exp(sum(-(v-ones(size(v))).^2,1));