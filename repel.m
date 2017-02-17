function cnf = repel(cnf, k_value, repel_steps, in_domain, s, outfile)
% repel(cnf, k_value, repel_steps,s, outfile)
% tries to distribute the configuration cnf of size dim x N, repelling 
% it in the direction of the Riesz gradient by a constant multiple of the 
% distance to the nearest neighbor. Uses the domain checker function in_domain
% cnf - dim x N matrix with the points
% k_value - the number of nearest neighbors used in the repel algorithm
% repel_steps - iterations of the step process to be made
% in_domain - domain checker; must take (x,y,z) as coordinates and return a
% boolean array of answers "in domain/not in the domain" for each point

dim = size(cnf,1);
bins = 100;
offset = 9;         % divides the minimal separation in the main loop
dim = size(cnf,1);
pt_num = size(cnf,2);   


[IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
IDX = IDX(:,2:end)';                     % drop the trivial first column in IDX
step = min(D(:,2));
% cutoff = (k_value*step)^2;

% % % % % % % % % % % % % % % % % % % % 
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
% D_old = D;       % uncomment this line to compare distribution/invoke plot3 below 


for iter=1:repel_steps
    cnf_repeated = reshape(repmat(cnf,k_value,1),dim,[]);
    knn_differences = cnf_repeated - cnf(:,IDX);
%       vectors pointing from each node to its (quasi-) k_value nearest 
%       neighbors
    knn_norms = sum(knn_differences.^2,1).^(0.5);
    riesz_weights = knn_norms.^(-s-1);
    directions = bsxfun(@times,riesz_weights,knn_differences);
    directions = reshape(directions, dim, k_value, []);
    directions = sum(directions,2);
%       sum along the dimension that contains (quasi-) k_value nearest
%       neighbors
    directions = reshape(directions, dim, []);
%       Riesz gradient for this node configuration 
    norms = sum(directions.^2,1);
    normals = directions.*(norms.^(-0.5));
    
    cnf_tentative = cnf + normals.*knn_norms(k_value*(0:pt_num-1)+1)/iter/offset;
    [domain_check, ~] = in_domain( cnf_tentative(1,:), cnf_tentative(2,:), cnf_tentative(3,:));
    cnf(:,domain_check) = cnf_tentative(:,domain_check); 
end
 

[~, D] = knnsearch(cnf', cnf', 'k', k_value+1);   
outtemp = min(D(:,2));
fprintf( outfile, 'Minimal separation after:      %f\n',  outtemp );
fprintf(   'Minimal separation after:      %f\n',  outtemp );
outtemp =  mean(D(:,2));
fprintf( outfile, 'Mean separation after:      %f\n',  outtemp);
fprintf( 'Mean separation after:      %f\n',  outtemp)

% % % % % % % % % % % % % % % % % % % % 
% % uncomment these lines to highlight nodes for which the distance to the
% % nearest neighbor has decreased/increased/remained the same
% clf;
% figure(5)
% pbaspect([1 1 1])
% hold on;
% CG=cnf(:,D_old(:,2)<D(:,2));
% CB=cnf(:,D_old(:,2)>D(:,2));
% CE=cnf(:,D_old(:,2)==D(:,2));
% fprintf(   'The number of nodes that have improved separation is:      %f\n',  sum(D_old(:,2)<D(:,2)) );
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
