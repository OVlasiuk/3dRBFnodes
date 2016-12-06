function cnf = repel(cnf, k_value, repel_steps,s, outfile)

bins = 100;
offset = 7;         % divides the minimal separation in the main loop
% riesz_s = @(x)riesz(x,s+1);
dim = size(cnf,1);
pt_num = size(cnf,2);
forces = zeros(size(cnf));        

[IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
IDX = IDX(:,2:end)';                     % drop the trivial first column in IDX

step = min(D(:,2));
fprintf( outfile, 'Minimal separation before repel steps:      %f\n', step)
fprintf( outfile, 'Mean separation before repel steps:      %f\n\n',  mean(D(:,2)) )
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


for iter=1:repel_steps
       cnf_neighbors = cnf(:,IDX);
       cnf_repeated = reshape(repmat(cnf,k_value,1), dim, k_value*pt_num); 
       directions = cnf_repeated - cnf_neighbors;
       inverse_norms_riesz = sum(directions.^2,1).^(-0.5*(s+1));      
       directions = bsxfun(@times,inverse_norms_riesz,directions);
       directions = sum(reshape(directions, dim, k_value, pt_num),2);
       directions = reshape(directions, dim, pt_num);
       inverse_norms = sum(directions.^2,1).^(-0.5);
       forces =  bsxfun(@times,inverse_norms,directions); 
    
    cnf_tentative = cnf + forces*step/offset/iter;
    domain_check = in_domain( cnf_tentative(1,:), cnf_tentative(2,:), cnf_tentative(3,:) );
       
    cnf = cnf + bsxfun(@times,domain_check,forces)*step/5/iter;
    
%     cnf(cnf<0) =  -cnf(cnf<0);
%     cnf(cnf>1) =  2-cnf(cnf>1);
end


[~, D] = knnsearch(cnf', cnf', 'k', k_value+1);   
fprintf( outfile, 'Minimal separation after:      %f\n',  min(D(:,2)) )
fprintf( outfile, 'Mean separation after:      %f\n',  mean(D(:,2)) )

% % % % % % % % % % % % % % % % % % % % 
% % %  histogram
figure(2);
h2 = histogram(D(:,2),bins);
h2.FaceColor = [0.9 0 0];       % red
saveas(h2,'./Output/histogram.png');
hold off;
