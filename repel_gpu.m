function cnf = repel_gpu(cnf, k_value, repel_steps,s)

bins = 50;
riesz_s = @(x)riesz(x,s+1);
dim = size(cnf,1);
pt_num = size(cnf,2);
directions = zeros(size(cnf));        

[IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
IDX = IDX(:,2:end)';                     % drop the trivial first column in IDX

step = min(D(:,2));
fprintf( 'Minimal separation before repel steps:      %f\n', step)
fprintf( 'Mean separation before repel steps:      %f\n\n',  mean(D(:,2)) )
% % % % % % % % % % % % % % % % % % % % 
% % %  histogram
fprintf('\n');
clf;
clf;
figure(2);
h1=histogram(D(:,2),bins);
h1.FaceColor = [0 0 0.9];        % blue
hold on;
% % % % % % % % % % % % % % % % % % % % 

G_neighbors = gpuArray( cnf(:,IDX) );
G_cnf = gpuArray( cnf );

for iter=1:repel_steps 
       cnf_repeated = reshape(repmat(G_cnf,k_value,1), dim, k_value*pt_num); 
       riesz_gradient = cnf_repeated - G_neighbors;                      
%      vectors pointing from each node to its k_value nearest neighbors
       inverse_norms_riesz = sum(riesz_gradient.^2,1).^(-0.5*(s+1));
%      norms of riesz_gradient raised to the power -s
       riesz_gradient = bsxfun(@times,inverse_norms_riesz,riesz_gradient);
       riesz_gradient = sum(reshape(riesz_gradient, dim, k_value, pt_num),2);
       riesz_gradient = reshape(riesz_gradient, dim, pt_num);
%      Riesz gradient for the node configuration       
       inverse_norms = sum(riesz_gradient.^2,1).^(-0.5);
       directions =  bsxfun(@times,inverse_norms,riesz_gradient);
%      normalized Riesz gradient
       
    G_cnf = G_cnf + directions*step/5/iter;
    G_cnf(G_cnf<0) =  -G_cnf(G_cnf<0);
    G_cnf(G_cnf>1) =  2-G_cnf(G_cnf>1);
    G_neighbors = G_cnf(:,IDX);
end

cnf = gather(G_cnf);


[~, D] = knnsearch(cnf', cnf', 'k', k_value+1);   
fprintf( 'Minimal separation after:      %f\n',  min(D(:,2)) )
fprintf( 'Mean separation after:      %f\n',  mean(D(:,2)) )

% % % % % % % % % % % % % % % % % % % % 
% % %  histogram
figure(2);
h2 = histogram(D(:,2),bins);
h2.FaceColor = [0.9 0 0];       % red
saveas(h2,'./Output/histogram.png');
hold off;
