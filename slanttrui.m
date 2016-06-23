function r=slanttrui(xyz)

n = [1,1,0]/sqrt(2);
proj = (xyz-dot(xyz,n)*n)';
S = null(n);
proj = proj.*[sqrt(2); sqrt(2);1];
r = trui(frac_part(S\proj)');