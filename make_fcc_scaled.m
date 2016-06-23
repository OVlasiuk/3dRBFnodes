function [M] = make_fcc_scaled(A,B,r,s,C,D)
%Returns fcc lattice of variable density on cube with opposing vertices
%A,B row vectors. Density scales linearly from separation radius r at vertex
%C to radius s*r at vertex D.
%Output is Nx3 array of coordinates.
%Rescale r
r = (1+s)/2*r;
len = abs(B(1)-A(1));
ilim = ceil(len/r);
jlim = ceil(len/(sqrt(3)/2*r));
klim = ceil(len/(sqrt(2/3)*r));
[X,Y,Z] = meshgrid(-ilim:ilim,-jlim:jlim,-klim:klim);
siz = (2*ilim+1)*(2*jlim+1)*(2*klim+1);

e1 = r*[1;0;0]; e2 = r*[1/2; sqrt(3)/2;0]; e3= r*[1/2; 1/sqrt(12);sqrt(2/3)];
M=ones(3,siz);
for m = 1:siz
M(:,m) = [X(m);Y(m);Z(m)];
end

%Remove vectors outside of cube.
M = [e1,e2,e3]*M;
M(:,unique(arrayfun(@ceil,find(M>len)/3))) = [];
M(:,unique(arrayfun(@ceil,find(M<0)/3))) = [];

%Rescale lattice in direction of largest density change.
V = D-C;
N = find(V>0);
M(N,:) = (M(N,:)+M(N,:).^2*(s-1)/(2*len))*2/(1+s);

N = find(V<0);
M(N,:) = (M(N,:)+M(N,:).^2*(s-1)/(2*len))*2/(1+s);
M(N,:) = len-M(N,:);


%Translate to vertex.
[~,sizeM] = size(M); 
M = M+bsxfun(@min,A,B)'*ones(1,sizeM);
M = M';

end