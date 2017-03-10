function num = boxRadiusNum(r)
% Returns the number of nodes in 3-dimensional unit cube, drawn from the
% irrational lattice with parameters (sqrt(2), sqrt(5)), such that the
% minimal separation of nodes between 27 such stacked cubes is
% approximately r (multiple copies are to account for distances between 
% neighboring cubes).
% For specifics of how the nodes are picked from the lattice, and for the
% scaling/shift we apply in each cube, see for example par_node


persistent lookup_table;
if isempty(lookup_table)
    load('unit_lattice_radius.mat');    
end
num = round(interp1(lookup_table,2:1000,r,'pchip'));