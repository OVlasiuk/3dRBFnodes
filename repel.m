function cnf = repel(cnf, k_value, repel_steps)

bins = 50;
LARGE_NUMBER = 4000;            % we will divide by it; really arbitrary

[IDX, D] = knnsearch(cnf, cnf, 'k', k_value+1);
forces = zeros(size(cnf));        

fprintf( 'Minimal separation before repulsion:      %f\n',  min(D(:,2)) )
fprintf( 'Mean separation before repulsion:      %f\n',  mean(D(:,2)) )
figure(3);
h1=histogram(D(:,2),bins);
h1.FaceColor = [0 0 0.9];
hold on;

for iter=1:repel_steps
    for i=1:size(cnf,1)
        force_i = normr(ones(k_value,1) * cnf(i,:) - cnf(IDX(i,2:end),:));
%         
        factor = 1/(iter*(LARGE_NUMBER + density(cnf(i,:))));
%         
        forces(i, :) = factor * sum(force_i, 1);
    end
    cnf = cnf + forces;
    cnf(cnf<0) =  -cnf(cnf<0);
    cnf(cnf>1) =  2-cnf(cnf>1);
end
% forces(1:10,:)

[~, D] = knnsearch(cnf, cnf, 'k', k_value+1);   

fprintf( 'Minimal separation after:      %f\n',  min(D(:,2)) )
fprintf( 'Mean separation after:      %f\n',  mean(D(:,2)) )
h2 = histogram(D(:,2),bins);
h2.FaceColor = [0.9 0 0];
hold off;