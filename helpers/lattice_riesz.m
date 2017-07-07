function lattice_riesz
%LATTICE_RIESZ
% Finds and saves the periodic Riesz minimizers in the unit cube; by
% default, the Riesz exponent is s=4.0. Also, generates the tables of mean
% and minimal separation distances for the minimizers; these are saved as 
% 'mrtable_riesz.mat' in the 'Output/' folder.
% The table of (coordinates of) minimizers is saved as a Matlab cell array
% 'lattice_riesz.mat'.


options = optimoptions(@fmincon,'Algorithm','interior-point','SpecifyObjectiveGradient',true,...
'Display','final','MaxFunctionEvaluations',10000,'MaxIterations',2000,'UseParallel',1);
ltable = cell(1,30);

parfor N=1:20
    A = [-eye(3*N); eye(3*N)];
    b = [zeros(3*N,1); ones(3*N,1)];
    Xinitial = rand(3,N);
    [xfinal fval exitflag output] = fmincon(@penergy,Xinitial,...
    A,b,[],[],[],[],[],options);
    ltable{N} = xfinal;
end

mtable = ones(1,numel(ltable));
rtable = ones(1,numel(ltable));
for N=2:numel(ltable)
    [~, D] = knnsearch(ltable{N}', ltable{N}','k',2);
    mtable(N) = mean(D(:,2));
    rtable(N) = min(D(:,2));
end

s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-13))
save('../Output/lattice_riesz.mat', 'ltable')
save('../Output/mrtable_riesz.mat', 'mtable','rtable')
cd(s_old)



% plot3(xfinal(1,:),xfinal(2,:),xfinal(3,:),'.b','MarkerSize',20)
% axis vis3d