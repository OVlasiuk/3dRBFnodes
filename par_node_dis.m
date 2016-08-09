dim = 3;
oct = 2^dim;
N = 15; %    number of boxes per side of the cube
max_nodes_per_box = 20; 
k_value = 15;           % number of nearest neighbors used in the knnsearch
repel_steps = 20;
r1 = sqrt(2);
r2 = sqrt(5);
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
% % % % % % % % % % % % % % % % % % 
corners = -ones(N^dim,dim);
nodes = zeros(N^dim, max_nodes_per_box+1,dim); % the vector (n,1,:) will contain the number 
                    %  of nodes in the n-th box
                    
tic
parfor i=1:N^dim
    corners(i,:) = [rem((i-1), N)/N  floor(rem(i-1, N^2)/N)/N floor((i-1)/N^2)/N];   % TODO
    eval_pts = num2cell(bsxfun(@plus, corners(i,:), cube_vectors/N),2);
    fun_values = cellfun(@trui, eval_pts);   
%     current_num_nodes = floor(max_nodes_per_box*mean(fun_values)/2);
    current_num_nodes = max_nodes_per_box-ceil(max_nodes_per_box * mean(fun_values));
    
    node = zeros(max_nodes_per_box+1,dim);
    corner1 = corners(i,:);
%   corner2 = corners(i,:)+1/N;     
    for j=1:current_num_nodes
        node(j+1,:) = [j/current_num_nodes/N,  frac_part(r1*j)/N,  frac_part(r2*j)/N];
        node(j+1,:) = node(j+1,:) + corner1;
    end           
    node(1,:) = node(1, :) + current_num_nodes;
    nodes(i,:,:) = node;
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


% % % % % % % % % % % % % 
fprintf( '\nNumber of nodes:      %d\n',  num_nodes )
fprintf( 'Mean number of nodes per box:      %d\n', mean(nodes(:,1,1) ))
fprintf( 'Max number of nodes per box:      %d\n', max(nodes(:,1,1) ))
fprintf( 'Min number of nodes per box:      %d\n', min(nodes(:,1,1) ))
toc
fprintf('\n');

% pbaspect([1 1 1])
% view([1 1 0])
% figure(2);
% plot3(cnf(:,1), cnf(:,2), cnf(:,3),  '.k');

fprintf( 'Performing %d repelling steps.\n',  repel_steps)
cnf = repel(cnf, k_value, repel_steps);
toc 

pbaspect([1 1 1])
% view([1 1 0])
figure(1);
plot3(cnf(:,1), cnf(:,2), cnf(:,3),  '.k');
save('slanttrui.mat', 'cnf')