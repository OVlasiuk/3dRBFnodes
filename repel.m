function cnf = repel(cnf, k_value, repel_steps)

bins = 50;
% LARGE_NUMBER = 8000;            % we will divide by it; really arbitrary

[IDX, D] = knnsearch(cnf, cnf, 'k', k_value+1);
forces = zeros(size(cnf));        

step = min(D(:,2));
% max_step = max(D(:,2));
fprintf( 'Minimal separation before repel steps:      %f\n', step  )
fprintf( 'Mean separation before repel steps:      %f\n\n',  mean(D(:,2)) )
% % % % % % % % % % % % % % % % % % % % 
% % %  histogram
figure(2);
h1=histogram(D(:,2),bins);
h1.FaceColor = [0 0 0.9];
hold on;
% % % % % % % % % % % % % % % % % % % % 
for iter=1:repel_steps
   parfor i=1:size(cnf,1)
       directions = ones(k_value,1) * cnf(i,:) - cnf(IDX(i,2:end),:);      %this has k_value rows
       force_coeffs = cellfun(@riesz,  num2cell(directions,2));            %this has k_value rows  
       forces(i, :) =   normr(force_coeffs' * normr(directions));         % 1 x k_value * k_value x 3
   end
    cnf = cnf + forces*step/3/iter;
    cnf(cnf<0) =  -cnf(cnf<0);
    cnf(cnf>1) =  2-cnf(cnf>1);
end

[~, D] = knnsearch(cnf, cnf, 'k', k_value+1);   
fprintf( 'Minimal separation after:      %f\n',  min(D(:,2)) )
fprintf( 'Mean separation after:      %f\n',  mean(D(:,2)) )

% % % % % % % % % % % % % % % % % % % % 
% % %  histogram
figure(2);
h2 = histogram(D(:,2),bins);
h2.FaceColor = [0.9 0 0];
hold off;