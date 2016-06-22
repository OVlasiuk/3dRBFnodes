dim = 3;
oct = 2^dim;
% 
% vertices of the unit cube
cube_vectors = zeros(oct,3);
count=1;
for i=0:1
    for j=0:1
        for k=0:1
            cube_vectors(count,:) = [i j k];
            count = count + 1;
        end
    end
end

N = 15; %    number of boxes per side of the cube

max_nodes_per_box = 40; 
tolerance = 1e-4;


corners = -ones(N^dim,dim);
nodes = zeros(N^dim, max_nodes_per_box+1,dim); % the first element or vector will contain the number 
                    %  of nodes in the current box
                    
tic

for i=1:N^dim
    corners(i,:) = [rem((i-1), N)/N  floor(rem(i-1, N^2)/N)/N floor((i-1)/N^2)/N];
    eval_pts = num2cell(bsxfun(@plus, corners(i,:), cube_vectors/N),2);
    fun_values = cellfun(@density, eval_pts);
    [max_dens, ind1] = max(fun_values);
    [min_dens, ind2] = min(fun_values);
    radius = 1/(1+3/2*min_dens);

    current_box = make_fcc_scaled(corners(i,:), corners(i,:)+1/N, radius, max_dens/min_dens,cell2mat(eval_pts(ind1)),cell2mat(eval_pts(ind2)));
        if(size(current_box,1)>max_nodes_per_box)
            fprintf ('Warning: Nodes in box exceeds maximum. Consider raising maximum or adjusting radius function')
            return
        end
    nodes = writeNodes(nodes, current_box, i);

end



% nodes_sparse = reshape(nodes(:,2:end,:),3,[]);       % this contains
% all the sparsity; we will turn nodes into a full matrix instead

num_nodes = sum(nodes(:,1,1));

cnf = zeros(num_nodes,dim);
count = 1;
for i=1:N^dim
    cur_num = nodes(i,1,1);
    cnf(count:(count+cur_num-1),:) = nodes(i, 2:cur_num+1, :);
    count = count + cur_num;
end
toc 
% % % % % % % % % % % % % 

pbaspect([1 1 1])
figure(1);
fprintf( 'Number of nodes:      %d\n\n',  num_nodes )
plot3(cnf(:,1), cnf(:,2), cnf(:,3),  '.k');


k_value = 15;           % number of nearest neighbors used in the knnsearch
repel_steps = 5;

cnf = repel(cnf, k_value, repel_steps);

pbaspect([1 1 1])
figure(2);
plot3(cnf(:,1), cnf(:,2), cnf(:,3),  '.k');
