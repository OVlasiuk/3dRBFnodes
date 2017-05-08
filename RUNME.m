function RUNME
%RUNME
% A demo script showcasing node distribution using irrational lattices and 
% Riesz energy minimization. Serves as a wrapper to the example routines in
% this folder; the filenames of routines must begin with 'node_'.
% 
%   See also NODE_DIS, NODE_EARTH.

s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-5))
addpath helpers/

fprintf('The following node distribution examples are available:\n');
D = dir('*.m');
names = cell(size(D));
for i=1:numel(D)
    names{i} = D(i).name;  
end
chnames = strvcat(names{:});
idx  = strmatch('node_', chnames);
I = cell(size(idx));
for i=1:numel(idx)
    fprintf('%d     %s\n', i, names{idx(i)});
    I{i} = num2str( i );
end
fprintf('\nMost of them will accept some input values, e.g., radial density,\n')
fprintf('point inclusion function, etc.\n')
fprintf('Type "help name_of_the_example" while in Matlab prompt to see details\n'); 
fprintf('and possible input values. Alternatively, type the number of a \n')
fprintf('script you''d like to run now with the default parameters, or press\n');
fprintf('Return to go back to Matlab prompt. Type "A" to list all m-files\n');
fprintf('in the current project. For most of these, typing "help name_of_the_file\n');
inp=input('will produce some additional info.\n','s');
if isempty(I)
	d(s_old)
    return
else
    switch inp
        case 'A'
            dir *.m;
            dir helpers*/*.m;
        case I
            run( names{ idx(eval(inp)) } )
    end
end