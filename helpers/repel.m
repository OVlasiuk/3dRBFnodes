function cnf = repel(cnf,...
                    N_moving,...
                    k_value,...
                    repel_steps,...
                    densityF,...
                    in_domainF,...                   
                    varargin...
                )
%REPEL 
% cnf = repel(cnf, N_moving, k_value, repel_steps, densityF, in_domainF,...                   
%                 'NAME1',VALUE1,...,'NAMEN',VALUEN)
% Tries to distribute the configuration cnf of size dim x N, repelling 
% it in the direction of the Riesz gradient by a constant multiple of the 
% distance to the nearest neighbor. Uses the domain checker function
% in_domainF.
% 
% cnf -- dim x N matrix with the points;
% N_moving -- move only the first N_moving points in cnf;
% k_value -- the number of nearest neighbors used in the repel algorithm;
% repel_steps -- iterations of the step process to be made;
% densityF -- the target radial density: desired distance to the nearest 
%   node; may be unstable if the function isn't Lip-1;
% in_domainF -- domain checker; must take (x,y,z) as (arrays of) coordinates
%   and return a boolean array of answers "in domain/not in the domain" for 
%   each point. Pass 0 or nothing to not perform any domain checks.
% 
% Optional argument name/value pairs:
% Name          Value
% 
% 'jitter'      a number between 0 and 1 serving as a factor of a 
%               random summand for the direction of repulsion.
% 
% 'pullback'    a mapping to use when pulling back nodes, pushed outside
%               the domain by repulsion.
% 
% 'A'           sidelength of the outer bounding cube: the repel process 
%               will be restricted to the cube [-A/2, A/2]^dim; default:
%               100;
% 
% 's'           the exponent used in the Riesz kernel; default: 4.0;
% 
% 'instats'     a logical value for whether to compute the nearest-neighbor
%               statistics for input; since the output is always processed,
%               useful to disable in a loop, etc; default: true;
% 
% 'histogram'   a logical value indicating whether to produce histograms of
%               the initial and post-repel distributions of the
%               nearest-neighbor distances; requires 'instats' to be true;
%               default: false;
% 
% 'bins'        the number of bins to use in the histogram; default: 200;
% 
% 'offset'      offset in the arithmetic progression, dividing the nearest
%               separation, that determines the step length; default: 18; 
% 
%   It is HIGHLY recommended to use either s=4.0 or s=0.5, as these are 
%   pre-coded, or to modify the source code. Otherwise you'll be using the 
%   Matlab's power function, which turns out to be not that great.

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Initialize variables
dim = size(cnf,1);
C = 0.034629556613235;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
pnames = { 'jitter' 'pullback' 'A'     's'     'histogram' 'bins' 'offset' 'instats'};
dflts =  { 0            []      100.0   4.0     false       200    18       true};
[jitter, pullbackF, A, s, htrue, bins, offset, instats, ~] =...
     internal.stats.parseArgs(pnames, dflts, varargin{:});
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if jitter==0
    noise = 0;
else 
    noise = @() jitter*normc(randn(dim,N_moving));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

switch s
    case 4.0
        compute_riesz = @(x) 1./x./x;
        compute_weights = @(x) x.*x.*x.*x;
    case 2.0
        compute_riesz = @(x) 1./x;
        compute_weights = @(x) x.*x;
    case 0.5
        compute_riesz = @(x) 1./sqrt(sqrt(x));
        compute_weights = @(x) sqrt(x);
    otherwise
        compute_riesz = @(x) sqrt(x).^(-s);
        compute_weights = @(x) x.^s;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

fprintf( '\nEntering the repel.m subroutine; a timer starts.\n\n')
tic
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if instats
    [IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
    IDX = IDX(:,2:end)';          % drop the trivial first column in IDX
    step = D(:,2);   
    fprintf( 'Minimal separation before repel steps:\t%3.8f\n', min(step))
    outtemp = mean(D(:,2));
    fprintf(   'Mean separation before repel steps:\t%3.8f\n\n',   outtemp)
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Histogram
if htrue && instats
    fprintf('\n')
    F = figure;
    h1=histogram(D(:,2),bins,'Normalization','probability');
    h1.FaceColor = [0 0 0.9];        % blue
    h1.EdgeAlpha=.1;
    hold on;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Main loop
for iter=1:repel_steps
    if mod(iter,10) == 1
        [IDX, ~] = knnsearch(cnf', cnf(:,1:N_moving)', 'k', k_value+1);
        IDX = IDX(:,2:end)';
    end
%% Vectors from nearest neighbors    
    cnf_repeated = reshape(repmat(cnf(:,1:N_moving),k_value,1),dim,[]);
    cnf_repeated_concentric = cnf_repeated./sqrt(sum(cnf_repeated.*cnf_repeated,1));
    knn_cnf = cnf(:,IDX);
    knn_differences = cnf_repeated - knn_cnf;    
    knn_norms_squared = sum(knn_differences.*knn_differences,1); 
%% Weights using radial density
    riesz_weights = compute_riesz(knn_norms_squared);
    if isa(densityF,'function_handle')
        knn_density =  densityF(knn_cnf);   
        density_weights = compute_weights(knn_density);
        weights = s*density_weights .* riesz_weights ./ ...
                                                        knn_norms_squared;
    else
        weights = s*riesz_weights./knn_norms_squared;
    end
%% Sum up over the nearest neighbors    
    gradient = bsxfun(@times,weights,knn_differences);% -...
%                       s*riesz_weights .* density_weights./knn_density .* cnf_repeated_concentric; 
    gradient = reshape(gradient, dim, k_value, []);
    gradient = reshape(sum(gradient,2), dim, []);
%% Add noise and renormalize    
    if isa(noise,'function_handle')
        gradient = gradient + noise() * mean(sqrt(sum(gradient.*gradient,1)));
    end
    directions = gradient./sqrt(sum(gradient.*gradient,1)); 
    step = sqrt(min(reshape(knn_norms_squared,k_value,[]),[],1));
    cnf_tentative = cnf(:,1:N_moving) +...
                            directions(:,1:N_moving).*step/(offset+iter-1); 
%% Detect the domain                        
    if exist('in_domainF', 'var') && isa(in_domainF,'function_handle')
        domain_check = in_domainF( cnf_tentative(1,:), cnf_tentative(2,:), cnf_tentative(3,:));
    else
        domain_check = ~any((cnf_tentative<-A/2.0) | (cnf_tentative>A/2.0),1);
    end
    if isa(pullbackF,'function_handle')
        cnf(:,~domain_check) = pullbackF(cnf_tentative(:,~domain_check)); 
    end    
    cnf(:,domain_check) = cnf_tentative(:,domain_check); 
end
toc

 
%% New separation
[~, D] = knnsearch(cnf', cnf', 'k', k_value+1);   
fprintf(   'Minimal separation after:\t\t%3.8f\n',  min(D(:,2)));
fprintf( 'Mean separation after:\t\t\t%3.8f\n',  mean(D(:,2)))

%% % % % % % % % % % % % % % % Histogram % % % % % % % % % % % % % % % % % 
if htrue    
    figure(F);
    h2 = histogram(D(:,2),bins,'Normalization','probability');
    h2.FaceColor = [0.9 0 0];       % red
    h2.EdgeAlpha=.1;
    set(gca,'FontSize',12)
    ylabel('Probability','FontSize',24);
    xlabel('Distance to the nearest neighbor','FontSize',24);
    hold off;
end
