function ratios = dcompare(pts,densityF, plotit)
%DCOMPARE
% ratios = dcompare(pts,densityF, plotit)
% Display statistics about how the radial density (aka distance to the
% nearest neighbor) compares with the values of function 'densityF' at the
% array of column vectors 'pts'; must be of size (dim)x(N).
% ratios -- array of values (densityF(pts) ./ D), where D are the distances
%   to nearest neighbors for columns of 'pts'.
% Input:
% pts -- 3x(numpts) array
% densityF -- function handle; must accepts 'pts' as argument.
% plotit -- pass 'y' or 1, etc., to plot the configuration 'pts'.
%   See also PT_ANALYZER
if size(pts,1) ~= 3
    pts = pts';
end

[~, D] = knnsearch(pts', pts', 'k', 2);
rdens_cnf = D(:,2);
rdens_fun = densityF(pts);
ratios = rdens_fun./rdens_cnf';
diff = abs(rdens_fun - rdens_cnf');
radii = sqrt( sum( pts .*pts,1 ) );
maxdiff = max(diff);
meandiff = mean(diff);
minratio = min(ratios);
maxratio = max(ratios);
quantile5 = quantile(ratios,0.05);
quantile95 = quantile(ratios,0.95);
meanratio = mean(ratios);
varratio = var(ratios);
fprintf('\t\t= DCOMPARE =\n');
fprintf('"diff" = abs(rho-D) \n');
fprintf('"ratio" = rho/D \n');
fprintf('\nmaxdiff\t\tmeandiff\n');
fprintf('%3.6f\t%3.6f\n\n',maxdiff,meandiff)
fprintf('minratio\tmaxratio\tquantile5\tquantile95\n');
fprintf('%3.6f\t%3.6f\t%3.6f\t%3.6f\t\n\n', minratio,maxratio,quantile5,quantile95)
fprintf('meanratio\tvarratio\n')
fprintf('%3.6f\t%3.6f\n',meanratio,varratio)

if exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
    msize = ceil(max(1, 22-5*log10(size(pts,2)) ));
    figure;
    plot(radii,ratios,'.k', 'MarkerSize',4)
    hold on;
    plot(radii,diff,'.g', 'MarkerSize',4)
    set(gca,'FontSize',12)
    xlabel('Radius {\bf\it{N}}','FontSize',24);
    ylabel('\rho({\bf\it{N}})/\Delta({\bf\it{N}})','FontSize',24);
    if ~usejava('desktop')
        print('ratio','-dpdf','-r300','-bestfit')
    end
    if (size(pts,1)==3) 
        figure   
        plot3(pts(1,ratios>quantile95),pts(2,ratios>quantile95),pts(3,ratios>quantile95),'.k','MarkerSize',msize)
        hold on;
        plot3(pts(1,ratios<quantile5),pts(2,ratios<quantile5),pts(3,ratios<quantile5),'.r','MarkerSize',msize)
        axis vis3d;
        daspect([1 1 1]);
        pbaspect([1 1 1]);
        warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
        [leg, ico] = legend('Large ratio rho/Delta','Small ratio rho/Delta');
        leg.Location = 'best';
        ico(4).MarkerSize = min(15, 10*msize);
        ico(6).MarkerSize = min(15, 10*msize);
        if ~usejava('desktop')
            print('error_location','-dpdf','-r300','-bestfit')
        end
    end
end
