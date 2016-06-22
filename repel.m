function cnf = repel(cnf, k_value, repel_steps)

% k_value = 15;           % number of nearest neighbors used in the knnsearch
% repel_steps = 1;
bins = 50;

[IDX, D] = knnsearch(cnf, cnf, 'k', k_value+1);
forces = zeros(size(cnf));        

fprintf( 'Minimal separation before repulsion:      %f\n',  min(D(:,2)) )
fprintf( 'Mean separation before repulsion:      %f\n',  mean(D(:,2)) )
figure(3);
h1=histogram(D(:,2),50);
h1.FaceColor = [0 0 0.8];
hold on;

% force_i = zeros(k_value, 3);              

for iter=1:repel_steps
    for i=1:size(cnf,1)
        force_i = normr(ones(k_value,1) * cnf(i,:) - cnf(IDX(i,2:end),:));
%         
        factor = 1/(iter*(1600 + density(cnf(i,:))));
%         
        forces(i, :) = factor * sum(force_i, 1);
    end
    cnf = cnf + forces;
    cnf(cnf<0) =  -cnf(cnf<0);
    cnf(cnf>1) =  2-cnf(cnf>1);
end

[~, D] = knnsearch(cnf, cnf, 'k', k_value+1);   

fprintf( 'Minimal separation after:      %f\n',  min(D(:,2)) )
fprintf( 'Mean separation after:      %f\n',  mean(D(:,2)) )
h2 = histogram(D(:,2),50);
h2.FaceColor = [0.8 0 0];