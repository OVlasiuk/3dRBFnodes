%ANDES
% Extracts and plots a piece of South-American coast next to the Andes,
% using the output of node_earth.m
% 
%      See also  NODE_EARTH
s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-5))
run('../node_earth')
format long;
bins = 200;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
[~, Dcnf] = knnsearch(cnf', cnf', 'k', 2);
figure(5);
pbaspect([1 1 1])
hcnf = histogram(Dcnf(:,2),bins);
hcnf.FaceColor = [0 0.9 0];
hcnf.EdgeAlpha=.1;
set(gca,'FontSize',12)
ylabel('Number of nodes','FontSize',24);
xlabel('Distance to the nearest neighbor','FontSize',24);
hold on;
separation_all = min(Dcnf(:,2))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
CNF = repmat(cnf,6,1);
e=[eye(3) -eye(3)];
shifted=bsxfun(@plus, separation_all*e(:)/2,CNF);
shifted = reshape(shifted,3,[]);
indices = ~in_domain( shifted(1,:), shifted(2,:), shifted(3,:));
indices = reshape(indices,6,[]);
II=sum(indices,1);
I=logical(II);
num_surface_nodes = sum(I)

[~, D] = knnsearch(cnfsurf', cnfsurf', 'k', 2);
figure(5);
hsurf = histogram(D(:,2),bins);
hsurf.FaceColor = [0.9 0 0];
hsurf.EdgeAlpha=.1;

figure(3)
rads = sum(cnf.^2,1);
inds=rads>1.19;                 % radii squared are close to 1.21 on the 
                                % outer boundary: the upper bound in in_domain 
cnfsurf=cnf(:,I&~inds);          
plot3(cnfsurf(1,:), cnfsurf(2,:), cnfsurf(3,:),  '.b','MarkerSize',1);
title('Surface nodes')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
separation_surface = min(D(:,2))


hemi=cnf(1,:)>0 & cnf(2,:)<0 & cnf(3,:)<0;
cnfhemi=cnf(:,I&~inds&hemi);

rot=[1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1/sqrt(2)];
hemirot=zeros(size(cnfhemi));

hemirot(3,:)=cnfhemi(3,:);
hemirot(1:2,:) = rot*cnfhemi(1:2,:);
r= [2 3 1];
chemi = hemirot(r,:);


Ihemi = chemi(1,:)<0 & chemi(3,:)>0.4;
chemirads = sqrt(sum(chemi.^2,1));
az = -67.5;
el = 32;
figure(10);
scatter3(chemi(1,Ihemi), chemi(2,Ihemi), chemi(3,Ihemi),...
    ones(1,size(chemi(Ihemi),2)) * 12,chemirads(Ihemi),'filled');
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',12)
% title('The Western coast of South America with color-encoded heights')
view(az,el); axis vis3d;
% print('andes_scatter','-dpdf','-r300','-bestfit')
figure(11);
plot3(chemi(1,Ihemi), chemi(2,Ihemi), chemi(3,Ihemi),  '.b',...
    'MarkerSize',5);
% title('The Western coast of South America')
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',12)
view(az,el);grid on; axis vis3d;
% print('andes','-dpdf','-r300','-bestfit')
cd(s_old)