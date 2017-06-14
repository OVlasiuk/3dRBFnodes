function cnf = repel(cnf,...
                    N_moving,...
                    k_value,...
                    repel_steps,...
                    ~,...
                    in_domainF,...
                    jitter,...
                    s,...
                    outfile,...
                    pullbackF...
                )
%REPEL 
% cnf = repel(cnf,...
%             N_moving,
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
% N_moving -- move only the first N_moving points in cnf;
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
switch s
    case 5.0
        compute_weights = @(x) 1./x./x./x;
    case 2.5
        compute_weights = @(x) 1./x./sqrt(x)./sqrt(sqrt(x));
    case 0.5
        compute_weights = @(x) 1./sqrt(x)./sqrt(sqrt(x));
    otherwise
        compute_weights = @(x) sqrt(x).^(-s-1);     
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if ~exist('jitter', 'var') || jitter==0
    noise = 0;
else 
    noise = @() jitter*normc(randn(dim,pt_num));
end
if ~exist('outfile', 'var')
    outfile = 0;
end
if ~exist('in_domainF', 'var') || ~isa(in_domainF,'function_handle')
     in_domainF = @(x, y, z) ones(size(x));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
[IDX, D] = knnsearch(cnf', cnf(:,1:N_moving)', 'k', k_value+1);
IDX = IDX(:,2:end)';          % drop the trivial first column in IDX
step = min(D(:,2));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
fprintf( '\nEntering the repel.m subroutine; a timer starts.\n\n')
tic
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if outfile ~= 0
    fprintf( outfile, 'Minimal separation before repel steps:      %3.8f\n', step);
end    
fprintf( 'Minimal separation before repel steps:      %3.8f\n', step)
outtemp = mean(D(:,2));
if outfile ~=0
    fprintf( outfile, 'Mean separation before repel steps:      %3.8f\n\n',   outtemp);
end
fprintf(   'Mean separation before repel steps:      %3.8f\n\n',   outtemp)
toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% %  histogram
fprintf('\n')
figure(2);
h1=histogram(D(:,2),bins);
h1.FaceColor = [0 0 0.9];        % blue
h1.EdgeAlpha=.1;
hold on;
% D_old = D;       % uncomment this line to compare distribution/invoke plot3 below 

%% Main loop
for iter=1:repel_steps
    cnf_repeated = reshape(repmat(cnf(:,1:N_moving),k_value,1),dim,[]);
%       vectors pointing from each node to its (quasi-) k_value nearest 
%       neighbors
    knn_differences = cnf_repeated - cnf(:,IDX);
    knn_norms_squared = sum(knn_differences.*knn_differences,1);               
    riesz_weights = compute_weights(knn_norms_squared);     
    gradient = bsxfun(@times,riesz_weights,knn_differences);
    gradient = reshape(gradient, dim, k_value, []);
% % % %   Riesz gradient for this node configuration 
    gradient = reshape(sum(gradient,2), dim, []);
%     Add jitter:
    if isa(noise,'function_handle')
        gradient = gradient + noise() * mean(sqrt(sum(gradient.*gradient,1)));
    end
    directions = gradient./sqrt(sum(gradient.*gradient,1));
%     step = sqrt(min(reshape(knn_norms_squared,k_value,[]),[],1));
    cnf_tentative = cnf(:,1:N_moving) +...
                        directions(:,1:N_moving).*step/(offset+iter-1);  
    domain_check = in_domainF( cnf_tentative(1,:), cnf_tentative(2,:), cnf_tentative(3,:));
% %                     ~sum((cnf_tentative<-A/2.0) + (cnf_tentative>A/2.0),1)
    if exist('pullbackF', 'var') && isa(pullbackF,'function_handle')
        cnf(:,~domain_check) = pullbackF(cnf_tentative(:,~domain_check)); 
    end    
    cnf(:,domain_check) = cnf_tentative(:,domain_check); 
end
toc

 
%% New separation
[~, D] = knnsearch(cnf', cnf', 'k', k_value+1);   
outtemp = min(D(:,2));
if outfile ~=0
    fprintf( outfile, 'Minimal separation after:      %3.8f\n',  outtemp );
end
fprintf(   'Minimal separation after:      %3.8f\n',  outtemp );
outtemp =  mean(D(:,2));
if outfile ~=0
    fprintf( outfile, 'Mean separation after:      %3.8f\n',  outtemp);
end
fprintf( 'Mean separation after:      %3.8f\n',  outtemp)
toc
% % % % % % % % % % % % % % % Histogram % % % % % % % % % % % % % % % % % %  
figure(2);
h2 = histogram(D(:,2),bins);
h2.FaceColor = [0.9 0 0];       % red
h2.EdgeAlpha=.1;
set(gca,'FontSize',12)
ylabel('Number of nodes','FontSize',16);
xlabel('Distance to the nearest neighbor','FontSize',16);
% saveas(h2,'./Output/histogram.png');
hold off;
