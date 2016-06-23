function r=slanttrui(xyz)

n = [1,1,1]/sqrt(3);
proj = (xyz-dot(xyz,n)*n)';
S = null(n);

T = [S(1,1),S(1,1);S(2,1), S(3,1); S(3,1),S(2,1)]/[1/2,-1/2,;-1/2,-1/2,;0,1];
proj = T*(proj - [0;1/2;-1/2]);

r = trui(frac_part(S\proj)');