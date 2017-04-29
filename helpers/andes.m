%ANDES
% Extracts and plots a piece of South-American coast next to the Andes,
% using the output of node_earth.m
% 
%      See also  NODE_EARTH
s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-5))
run('../node_earth')

CNF = repmat(cnf,6,1);
e=[eye(3) -eye(3)];
shifted=bsxfun(@plus, 0.007156*e(:),CNF);
shifted = reshape(shifted,3,[]);
indices = ~in_domain( shifted(1,:), shifted(2,:), shifted(3,:));
indices = reshape(indices,6,[]);
II=sum(indices,1);
I=logical(II);
num_surface_nodes = sum(I)

rads = sum(cnf.^2,1);
inds=rads>1.19;
cnfsurf=cnf(:,I&~inds);
plot3(cnfsurf(1,:), cnfsurf(2,:), cnfsurf(3,:),  '.b','MarkerSize',1);
bins = 100;
[~, D] = knnsearch(cnfsurf', cnfsurf', 'k', 2);
figure(5);
h2 = histogram(D(:,2),bins);
h2.FaceColor = [0.9 0 0];  
hold on;
[~, Dcnf] = knnsearch(cnf', cnf', 'k', 2);
figure(5);
hcnf = histogram(Dcnf(:,2),bins);
hcnf.FaceColor = [0 0.9 0]; 

separation_surface = min(D(:,2))
separation_all = min(Dcnf(:,2))

hemi=cnf(1,:)>0 & cnf(2,:)<0 & cnf(3,:)<0;
cnfhemi=cnf(:,I&~inds&hemi);

rot=[1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1/sqrt(2)];
hemirot=zeros(size(cnfhemi));

hemirot(3,:)=cnfhemi(3,:);
hemirot(1:2,:) = rot*cnfhemi(1:2,:);
figure(10);plot3(hemirot(1,:), hemirot(2,:), hemirot(3,:),  '.b','MarkerSize',5);
r= [2 3 1];
chemi = hemirot(r,:);

figure(10);plot3(chemi(1,:), chemi(2,:), chemi(3,:),  '.b','MarkerSize',5);
xlabel('x')
ylabel('y')
zlabel('z')
Ihemi = chemi(1,:)<0 & chemi(3,:)>0.4;

xlabel('x')
ylabel('y')
zlabel('z')
figure(10);plot3(chemi(1,Ihemi), chemi(2,Ihemi), chemi(3,Ihemi),'.b','MarkerSize',5);

chemirads = sqrt(sum(chemi.^2,1));
az = -67.5;
el = 32;
figure(10);scatter3(chemi(1,Ihemi), chemi(2,Ihemi), chemi(3,Ihemi),...
    ones(1,size(chemi(Ihemi),2)) * 13,chemirads(Ihemi),'filled');view(az,el); axis vis3d;
 print('andes_scatter','-dpdf','-r300','-bestfit')
figure(11);plot3(chemi(1,Ihemi), chemi(2,Ihemi), chemi(3,Ihemi),  '.b',...
    'MarkerSize',5);view(az,el);grid on; axis vis3d;
print('andes','-dpdf','-r300','-bestfit')
cd(s_old)