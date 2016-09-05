% % % % % % % % % MAIN SCRIPT FOR NODE SETTING: PARALLEL CPUS % % % % % % %
% % % % % % % % % % % % % PARAMETERS  % % % % % % % % % % % % % % % % % % %
                                                                                       % TODO replace index counts and such with masks
                                                                                       % TODO column-major adjustment
dim = 3;
oct = 2^dim;
N = 50;                                                                                % number of boxes per side of the cube
max_nodes_per_box = 15;  
k_value = 15;                                                                          % number of nearest neighbors used in the knnsearch
repel_steps = 1;
repel_power = 5;
r1 = sqrt(2);
r2 = sqrt(5);
density = @trui;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


cube_vectors = zeros(dim, oct);                                                          % populate vertices of the unit cube       
for i=1:dim
    len = 2 ^ (dim-i);
    for j=0:2^i-1                       
        cube_vectors(i, j*len+1:(j+1)*len) = mod(j,2);
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
nodes = zeros(dim, max_nodes_per_box, N^dim); 
node_indices = false(max_nodes_per_box, N^dim);                                          % logical-styled integer array to be used as mask
box_indices = false (N^dim,1);                                                            % boxes on the boundary
corners = -ones(dim,N^dim);
                                                                                        

                                                
tic
parfor i=1:N^dim
    corners(:,i) = [rem((i-1), N)  floor(rem(i-1, N^2)/N)  floor((i-1)/N^2)]/N;
%    floor(rem(i-1, N^2)/N)
    % TODO: this is not dimension-independent
                                                                                        %       also, it is quite nasty
                                                                                        %       
    eval_pts = num2cell(bsxfun(@plus, corners(:,i), cube_vectors/N),1);
    fun_values = cellfun(density, eval_pts);               
    current_num_nodes = min(max_nodes_per_box-ceil(max_nodes_per_box * mean(fun_values)), max_nodes_per_box);    
    box = zeros(dim, max_nodes_per_box);    
    for j=1:current_num_nodes
        box(:,j) = [j/current_num_nodes;  frac_part(r1*j);  frac_part(r2*j)]/N;     % TODO: can this be vectorized?
        box(:,j) = box(:,j) + corners(:,i);
    end    
    nodes(:,:,i) = box;   
    if max(fun_values)>=.7                                                               % this definitely can be vectorized easily
        box_indices(i) = true;        
    end
    temp = false(max_nodes_per_box,1);
    temp(1:current_num_nodes) = true;                         
    node_indices(:, i) = temp;                                                           % crazy stuff to make parfor work

end
toc

bdry_nodes = nodes(:, :, box_indices);                                                  % coordinates of the nodes that belong to boundary boxes
flat = num2cell(reshape(bdry_nodes, dim, []), 1);                                       % 
f_vals = cellfun(density, flat);                                                        % values of the density function at those nodes
bdry_removed_indices = reshape( f_vals >.7, max_nodes_per_box, []);                     % indices of nodes that have to be removed
node_indices(:, box_indices) = node_indices(:, box_indices) .* ~bdry_removed_indices;   % indices in the global nodeset
node_indices_flat = reshape(node_indices, 1, []);  
nodes = reshape (nodes, dim, []);
cnf = nodes(:, node_indices_flat);                                                           % after removing all nodes with zero density
                             
% % % % % % % % % % % % % % 
fprintf( '\nNumber of nodes:      %d\n',  length(cnf))
fprintf( 'Mean number of nodes per box:      %d\n', mean(sum(node_indices,1) ))
fprintf( 'Max number of nodes per box:      %d\n', max(sum(node_indices,1) ))
fprintf( 'Min number of nodes per box:      %d\n', min(sum(node_indices,1) ))
toc
fprintf('\n');

% pbaspect([1 1 1])
% view([1 1 0])
% figure(2);
% plot3(cnf(1,:), cnf(2,:), cnf(3,:),  '.k');

fprintf( 'Performing %d repel steps.\n',  repel_steps)
cnf = repel(cnf', k_value, repel_steps, repel_power);
toc 

pbaspect([1 1 1])
% view([1 1 0])
figure(1);
plot3(cnf(:,1), cnf(:,2), cnf(:,3),  '.k');
save('slanttrui.mat', 'cnf')
