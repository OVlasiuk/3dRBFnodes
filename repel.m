function cnf = repel(cnf, k_value, repel_steps,s)

bins = 50;
riesz_s = @(x)riesz(x,s+1);
dim = size(cnf,1);

[IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
forces = zeros(size(cnf));        
IDXflat = reshape(IDX', 1, []);

cnf_neighbors = cnf(:,IDXflat);
cnf_neighbors = reshape(cnf_neighbors, dim, k_value+1, size(cnf,2));
cnf_neighbors = cnf_neighbors(:,2:end,:);

step = min(D(:,2));
fprintf( 'Minimal separation before repel steps:      %f\n', step  )
fprintf( 'Mean separation before repel steps:      %f\n\n',  mean(D(:,2)) )
% % % % % % % % % % % % % % % % % % % % 
% % %  histogram
figure(2);
h1=histogram(D(:,2),bins);
h1.FaceColor = [0 0 0.9];        % blue
hold on;
% % % % % % % % % % % % % % % % % % % % 

% cnf_full = zeros(dim, k, size(cnf,1) );
% for i=1:size(cnf,1)
%     cnf_full(:,:,i) = 
% end



for iter=1:repel_steps
   parfor i=1:size(cnf,2)      
%        directions =   cnf(:,i)* ones(1,k_value) - cnf(:,IDX(i,2:end));
       directions =   cnf(:,i)* ones(1,k_value) - cnf_neighbors(:,:,i);          %this has dim rows, k_value columns
       force_coeffs = cellfun(riesz_s,  num2cell(directions,1));                 %this has k_value rows  
       forces(:,i) =  normc(directions * force_coeffs' );                        % 1 x k_value * k_value x 3
   end
    cnf = cnf + forces*step/3/iter;
    cnf(cnf<0) =  -cnf(cnf<0);
    cnf(cnf>1) =  2-cnf(cnf>1);
    cnf_neighbors = cnf(:,IDXflat);
    cnf_neighbors = reshape(cnf_neighbors, dim, k_value+1, size(cnf,2));
    cnf_neighbors = cnf_neighbors(:,2:end,:);
end

[~, D] = knnsearch(cnf', cnf', 'k', k_value+1);   
fprintf( 'Minimal separation after:      %f\n',  min(D(:,2)) )
fprintf( 'Mean separation after:      %f\n',  mean(D(:,2)) )

% % % % % % % % % % % % % % % % % % % % 
% % %  histogram
figure(2);
h2 = histogram(D(:,2),bins);
h2.FaceColor = [0.9 0 0];       % red
hold off;