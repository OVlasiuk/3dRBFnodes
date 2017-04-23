function cnf = repel(cnf,...
                    k_value,...
                    repel_steps,...
                    A,...
                    in_domainF,...
                    densityF,...
                    jitter,...
                    s,...
                    outfile)
%REPEL 
% cnf = repel(cnf,...
%             k_value,...
%             repel_steps,...
%             A,...
%             in_domainF,...
%             densityF,...
%             jitter,...
%             s,...
%             outfile)
% Tries to distribute the configuration cnf of size dim x N, repelling 
% it in the direction of the Riesz gradient by a constant multiple of the 
% distance to the nearest neighbor. Uses the domain checker function
% in_domainF.
% 
% cnf -- dim x N matrix with the points;
% k_value -- the number of nearest neighbors used in the repel algorithm;
% repel_steps -- iterations of the step process to be made;
% A -- sidelength of the outer bounding cube: the repel process will be
%   restricted to the cube [-A/2, A/2]^3.
% in_domainF -- domain checker; must take (x,y,z) as (arrays of) coordinates
%   and return a boolean array of answers "in domain/not in the domain" for 
%   each point. Pass 0 or nothing to not perform any domain checks.
% densityF -- determines radial distance to the nearest node;
% jitter -- a number between 0 and 1 serving as a factor of a random summand
%   for the direction of repulsion.
% s -- the exponent used in the Riesz kernel;
%   It is HIGHLY recommended to use either s=5.0 or s=0.5, as these are 
%   pre-coded, or to modify the source code. Otherwise you'll be using the 
%   Matlab's power function, which turns out to be not that great.
% outfile -- the log is printed to this file. Pass 0 to only print to
%   console;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
dim = size(cnf,1);
pt_num = size(cnf,2);   
bins = 200;
offset = 3;         % divides the minimal separation in the main loop
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if ~exist('s', 'var')
    s = 5.0;
end
if s==5.0
    compute_weights = @(x) 1./x./x./x;
else
    if s==0.5
        compute_weights = @(x) 1./x./sqrt(x);
    else
        compute_weights = @(x) sqrt(x).^(-s-1);
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
fprintf( '\nEntering the repel.m subroutine.\n\n')
if ~exist('jitter', 'var') || jitter==0
    noise = @() 0;
else 
    noise = @() jitter*normc(randn(dim,pt_num));
end
if ~exist('outfile', 'var')
    outfile = 0;
end
if ~exist('in_domainF', 'var') || ~isa(in_domainF,'function_handle')
     in_domainF = @(x,y,z) ones(size(x));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
tic
[IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
IDX = IDX(:,2:end)';          % drop the trivial first column in IDX
step = min(D(:,2));
% cutoff = (k_value*step)^2;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if outfile ~= 0
    fprintf( outfile, 'Minimal separation before repel steps:      %f\n', step);
end    
fprintf( 'Minimal separation before repel steps:      %f\n', step)

outtemp = mean(D(:,2));
if outfile ~=0
    fprintf( outfile, 'Mean separation before repel steps:      %f\n\n',   outtemp);
end
fprintf(   'Mean separation before repel steps:      %f\n\n',   outtemp)
toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% %  histogram

fprintf('\n')
clf;
clf;
figure(2);
h1=histogram(D(:,2),bins);
h1.FaceColor = [0 0 0.9];        % blue
hold on;
% D_old = D;       % uncomment this line to compare distribution/invoke plot3 below 

%% Main loop
for iter=1:repel_steps
    
    cnf_repeated = reshape(repmat(cnf,k_value,1),dim,[]);
%       vectors pointing from each node to its (quasi-) k_value nearest 
%       neighbors
    knn_differences = cnf_repeated - cnf(:,IDX);
    knn_norms = sum(knn_differences.*knn_differences,1);               
    riesz_weights = compute_weights(knn_norms);     
    directions = bsxfun(@times,riesz_weights,knn_differences);
    directions = reshape(directions, dim, k_value, []);
    directions = reshape(sum(directions,2), dim, []);
%     Add jitter:
    directions = directions + noise() * mean(sqrt(sum(directions.*directions,1)));

% % % % % % %    Uncomment to use external field term.    % % % % % % % % %
%     cnf_for_gradient = reshape(repmat(cnf,dim,1),dim,[]);
%     finite_difference1 = repmat(eye(3), 1, size(cnf,2));
%     finite_difference2 = repmat(-eye(3), 1, size(cnf,2));
%     cnf_for_gradient1 = cnf_for_gradient + step*finite_difference1;
%     cnf_for_gradient2 = cnf_for_gradient + step*finite_difference2;
%     density_gradient = (densityF(cnf_for_gradient1)-densityF(cnf_for_gradient2))/step/2.0;
%     density_gradient = reshape(density_gradient, dim,[]);
%     directions = directions + pt_num^(s/dim) * density_gradient;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % %   Riesz gradient for this node configuration 
    normals = directions./sqrt(sum(directions.*directions,1));
    cnf_tentative = cnf + normals.*step/offset/iter;
    domain_check = ~sum(cnf_tentative<-A/2.0 + cnf_tentative>A/2.0,1) &...
         in_domainF( cnf_tentative(1,:), cnf_tentative(2,:), cnf_tentative(3,:));
    cnf(:,domain_check) = cnf_tentative(:,domain_check); 
end
 
%% New separation
[~, D] = knnsearch(cnf', cnf', 'k', k_value+1);   
outtemp = min(D(:,2));
if outfile ~=0
    fprintf( outfile, 'Minimal separation after:      %f\n',  outtemp );
end
fprintf(   'Minimal separation after:      %f\n',  outtemp );
outtemp =  mean(D(:,2));
if outfile ~=0
    fprintf( outfile, 'Mean separation after:      %f\n',  outtemp);
end
fprintf( 'Mean separation after:      %f\n',  outtemp)
toc

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % Histogram % % % % % % % % % % % % % % % % % %  
figure(2);
h2 = histogram(D(:,2),bins);
h2.FaceColor = [0.9 0 0];       % red
% saveas(h2,'./Output/histogram.png');
hold off;
