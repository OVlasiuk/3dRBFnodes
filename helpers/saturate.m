function cnf = saturate(cnf, rdensity, in_domainF)
%SATURATE
% cnf = saturate(cnf, rdensity, in_domainF)
% Find all the Voronoi centers of the configuration 'cnf', contained inside
% the domain with the indicator function 'in_domainF', for which the hole
% depth multiplied by 0.9 is larger than the value of 'rdensity' at that
% respective center.
% cnf -- must be of the size 3x(num_pts)
% rdensity -- density, must accept arrays of size 3x(num_pts) and return
%   vectors of length num_pts
% in_domainF -- must accept 3 vectors (x, y, z) and return a logical vector
%   of the same dimension as, say, x.

[V,~] = voronoin(cnf');
V = V(in_domainF( V(:,1),V(:,2),V(:,3)) , : ); 
[~, holedepths] = knnsearch(cnf',V);
good_inds = (rdensity(V') < .9 * holedepths');
Vg = V(good_inds,:);

[sortedDists, sortHoles] = sort(holedepths(good_inds),'descend');
sortedVg = Vg(sortHoles,:);
sortedDensity = rdensity(sortedVg');
numHoles = size(sortedVg,1)
holeFill = false(1, numHoles);
new = true;
cutoff = 0;
sumHoleFill = 0;

for holeInd =1:numHoles
    if ~any(holeFill)
        holeFill(holeInd) = true;
        continue
    end
     if (mod(sumHoleFill, 5e3) == 1) && (sumHoleFill>1) && new
%         holeInd
%         sumHoleFill
        cutoff = holeInd-1;
        indlarge = holeFill & [true(1,cutoff),false(1,numHoles-cutoff)];
        tic
        ns = createns(sortedVg(indlarge,:)', 'nsmethod','kdtree');
        toc
        new = false;
     end
    indsmall = holeFill & [false(1,cutoff),true(1,numHoles-cutoff)];
    distMatrix = bsxfun(@minus, sortedVg(holeInd,:), sortedVg(indsmall,:))';
    dsmall = min(sqrt(sum(distMatrix.*distMatrix,1))); 
    if cutoff
        [~, dlarge] = knnsearch(ns,sortedVg(holeInd,:));
    else
        dlarge = sortedDists(holeInd);
    end
    if isempty(dsmall)
        dsmall = sortedDists(holeInd);
    end
    d = min([dsmall, dlarge, sortedDists(holeInd)]);
    if  (sortedDensity(holeInd) <   d)
       holeFill(holeInd) = true;
       sumHoleFill = sumHoleFill + 1;
       new = true;
    end
end
sumHoleFill 
Vnew = sortedVg(holeFill,:);
cnf = [cnf Vnew'];
% r = dcompare(cnf, rdensity);