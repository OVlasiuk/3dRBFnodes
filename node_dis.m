dim = 3;
oct = 2^dim;

cube_vectors = zeros(3,oct);
count=1;
for i=0:1
    for j=0:1
        for k=0:1
            cube_vectors(:,count) = [i j k];
            count = count + 1;
        end
    end
end



N = 15; %    number of boxes per side of the cube

max_nodes_per_box = 5;  %00;
tolerance = 1e-4;


corners = -ones(3, N^3);
nodes = zeros(3, max_nodes_per_box+1, N^dim); % the first element or vector will contain the number 
                    %  of nodes in the current box


                    
tic

for i=1:N^dim
    corners(:,i) = [rem((i-1), N)/N  floor(rem(i-1, N^2)/N)/N floor((i-1)/N^2)/N];
    evaluation_pts = num2cell(bsxfun(@plus, corners(1:3,i), cube_vectors/N),1);
    fun_values = cellfun(@density, evaluation_pts);
%     if max(fun_values)-min(fun_values) > tolerance
        current_box = irrational_nodes(corners(:,i), corners(:,i)+1/N, floor(max_nodes_per_box*mean(fun_values)));
        nodes = writeNodes(nodes, current_box, i);
%     end       %  run this block if density does not vary much, else - use
%                  directed density
end



% nodes_sparse = reshape(nodes(:,2:end,:),3,[]);       % this contains
% all the sparsity; we will turn nodes into a full matrix instead
num_nodes = sum(nodes(1,1,:));


nodes_consecutive = zeros(3, num_nodes);
count = 1;
for i=1:N^dim
    cur_num = nodes(1,1,i);
    nodes_consecutive(:, count:(count+cur_num-1)) = nodes(:, 2:cur_num+1, i);
    count = count + cur_num;
end
toc
fprintf( 'number of nodes:      %d\n\n',  num_nodes )
plot3(nodes_consecutive(1,:), nodes_consecutive(2,:), nodes_consecutive(3,:),  '.k');

% repel((nodes_consecutive));