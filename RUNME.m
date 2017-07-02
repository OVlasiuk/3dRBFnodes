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

fprintf('  The following node distribution examples are available:\n');
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
fprintf('point inclusion function, etc. Type "help name_of_the_example" while\n')
fprintf('in Matlab prompt to see details and possible input values.\n'); 
fprintf('Now you have the following options:\n \t-- enter the number of a script you''d like to run with the \n')
fprintf('default parameters;\n');
fprintf('\t-- press Return to go back to Matlab prompt;\n');
fprintf('\t-- type "A<Return>" to list all m-files in the current project.\n')
fprintf('For most of these, typing "help name_of_the_file" will produce\n');
inp=input('some additional info.\n>> ','s');
if isempty(I)
	d(s_old)
    return
else
    switch inp
        case {'a','A'}
            [s, r] = system('git ls-files | egrep "*.m\>"');
            if s ~=0
                dir *.m;
                dir helpers*/*.m;
            else
                fprintf('\n%s',r)
            end
        case I
            run( names{ idx(eval(inp)) } )
    end
end