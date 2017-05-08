%RUNME
% A demo script showcasing node distribution using irrational lattices and 
% Riesz energy minimization. Serves as a wrapper to the example routines in
% this folder; the filenames of routines must begin with 'node_'.
% 
%   See also NODE_DIS, NODE_EARTH.

s = char(mfilename('fullpath'));
cd(s(1:end-5))
addpath helpers/

fprintf('The following node distribution examples are available:\n');
D = dir;
names = cell(size(D));
for i=1:numel(D)
    names{i} = D(i).name;  
end
chnames = strvcat(names{:});
idx  = strmatch('node_', chnames);
for i=1:numel(idx)
    fprintf('%d     %s\n', i, names{idx(i)});
end
fprintf('\nMost of them will accept some input values, e.g., radial density,\n')
fprintf('point inclusion function, etc.\n')
fprintf('Type *help name_of_the_example* while in Matlab prompt to see details\n'); 
fprintf('and possible input values. Alternatively, type the number of a \n')
fprintf('script you''d like to run now with the default parameters, or press\n');
I=input('return to go back to Matlab prompt.\n');
if (I > 0) & (I <= numel(idx))
    run(names{idx(I)})
else
    return
end