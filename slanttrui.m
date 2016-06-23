function r=slanttrui(xyz)

n = [1,1,1]/sqrt(3);
proj = (xyz-dot(xyz,n)*n)';
S = null(n);

r = trui(frac_part(S\proj)');