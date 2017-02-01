function cnf = repel(cnf, k_value, repel_steps,s, outfile)

bins = 100;
offset = 9;         % divides the minimal separation in the main loop
dim = size(cnf,1);
pt_num = size(cnf,2);   


[IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
IDX = IDX(:,2:end)';                     % drop the trivial first column in IDX
step = min(D(:,2));
% cutoff = (k_value*step)^2;
var(D(:,2))
fprintf( outfile, 'Minimal separation before repel steps:      %f\n', step);
fprintf( 'Minimal separation before repel steps:      %f\n', step)
outtemp = mean(D(:,2));
fprintf( outfile, 'Mean separation before repel steps:      %f\n\n',   outtemp);
fprintf(   'Mean separation before repel steps:      %f\n\n',   outtemp)
% % % % % % % % % % % % % % % % % % % % 
% % %  histogram
fprintf('\n')
clf;
clf;
figure(2);
h1=histogram(D(:,2),bins);
h1.FaceColor = [0 0 0.9];        % blue
hold on;
% D_old = D;       % uncomment this line to compare distribution below /
% invoke plot3

% % % % % % % % % % % % % % % % % % % % 
for iter=1:repel_steps
       cnf_neighbors = cnf(:,IDX);
       cnf_repeated = reshape(repmat(cnf,k_value,1), dim, []); 
       riesz_gradient = cnf_repeated - cnf_neighbors;
%      vectors pointing from each node to its k_value nearest neighbors
       norms_riesz = sum(riesz_gradient.^2,1);
       inverse_norms_riesz=norms_riesz.^(-0.5*(s+1)); % .*(norms_riesz<cutoff)
%      norms of riesz_gradient raised to the power -s-1
       riesz_gradient = bsxfun(@times,inverse_norms_riesz,riesz_gradient);
       riesz_gradient = sum(reshape(riesz_gradient, dim, k_value, pt_num),2);
       riesz_gradient = reshape(riesz_gradient, dim, pt_num);
%      Riesz gradient for the node configuration       
       inverse_norms = sum(riesz_gradient.^2,1).^(-0.5);
       directions =  bsxfun(@times,inverse_norms,riesz_gradient); 
%      normalized Riesz gradient
    cnf_tentative = cnf + directions*step/offset/iter;
    domain_check = in_domain( cnf_tentative(1,:), cnf_tentative(2,:), cnf_tentative(3,:));   
    cnf(domain_check) = cnf_tentative(domain_check);
end
% % % % % % % % % % % % % % % % % % % % 




[~, D] = knnsearch(cnf', cnf', 'k', k_value+1);   
outtemp = min(D(:,2));
fprintf( outfile, 'Minimal separation after:      %f\n',  outtemp );
fprintf(   'Minimal separation after:      %f\n',  outtemp );
outtemp =  mean(D(:,2));
var(D(:,2))
fprintf( outfile, 'Mean separation after:      %f\n',  outtemp);
fprintf( 'Mean separation after:      %f\n',  outtemp)
% % % % % % % % % % % % % % % % % % % % 
% uncomment these lines to highlight nodes for which the distance to the
% nearest neighbor has decreased/increased/remained the same
% clf;
% figure(3)
% pbaspect([1 1 1])
% hold on;
% CG=cnf(:,D_old(:,2)<D(:,2));
% CB=cnf(:,D_old(:,2)>D(:,2));
% CE=cnf(:,D_old(:,2)==D(:,2));
% plot3(CG(1,:), CG(2,:), CG(3,:),  '.g','MarkerSize',1);
% plot3(CB(1,:), CB(2,:), CB(3,:),  '.r','MarkerSize',1);
% plot3(CE(1,:), CE(2,:), CE(3,:),  '.k','MarkerSize',1);
% % % % % % % % % % % % % % % % % % % % 
% % %  histogram
figure(2);
h2 = histogram(D(:,2),bins);
h2.FaceColor = [0.9 0 0];       % red
saveas(h2,'./Output/histogram.png');
hold off;
