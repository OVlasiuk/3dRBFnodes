function rad = latradius(COUNT)

% COUNT = 2000;
rad = zeros(1,COUNT);

N = 3;
dim = 3; 
r1 = sqrt(2);
r2 = sqrt(5);

I=1:N^dim;
corners = [rem((I-1), N);  floor(rem(I-1, N^2)/N);  floor((I-1)/N^2)]/N;

parfor max_nodes_per_box=2:COUNT
    delta = 1/(32* N * max_nodes_per_box^(1/dim));
    cube_shrink = 1 - max_nodes_per_box^(-1/dim)/3;
    box = zeros(dim, max_nodes_per_box);    
    for j=1:max_nodes_per_box
        box(:,j) = cube_shrink * [j/max_nodes_per_box;  mod(r1*j,1);  mod(r2*j,1)]/N;       % the box is shrunk to account for inwards corner % TODO: can this be vectorized?
        box(:,j) = box(:,j) +  delta;
    end
    nodes = reshape(repmat(corners,max_nodes_per_box,1),dim, []) +  reshape( repmat(box,1,size(corners,2)), dim, []);
%     nodes = box;
    [~, D] = knnsearch(nodes', nodes','k',2);
    rad(max_nodes_per_box) = min(D(:,2));
end

