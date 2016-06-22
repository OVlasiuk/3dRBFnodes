dim = 3;
oct = 2^dim;

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

max_nodes_per_box = 40;  %00;
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

nodes_consecutive = zeros(num_nodes,dim);
count = 1;
for i=1:N^dim
    cur_num = nodes(i,1,1);
    nodes_consecutive(count:(count+cur_num-1),:) = nodes(i, 2:cur_num+1, :);
    count = count + cur_num;
end
toc

figure(1);
fprintf( 'Number of nodes:      %d\n\n',  num_nodes )
plot3(nodes_consecutive(:,1), nodes_consecutive(:,2), nodes_consecutive(:,3),  '.k');



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % %  REPULSION  % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

k_value = 15;           % number of nearest neighbors used in the knnsearch
repel_steps = 1;

cnf = nodes_consecutive;

IDX = knnsearch(cnf, cnf, 'k', k_value+1);
forces = zeros(size(cnf,1), size(cnf,2));        

for iter=1:repel_steps
    for i=1:size(cnf,1)
        force_i = normr(ones(k_value,1) * cnf(i,:) - cnf(IDX(i,2:end),:));
%         
        factor = 1/(iter*(100 + density(cnf(i,:))));
%         
        forces(i, :) = factor * sum(force_i, 1);
    end
    cnf = cnf + forces;
end



figure(2);
% repel(nodes_consecutive);
plot3(cnf(:,1), cnf(:,2), cnf(:,3),  '.k');
tmp = cnf - nodes_consecutive;