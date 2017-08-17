function x = pback(x, varargin)
%PBACK
% Takes matrices of size (dim x N)


pnames = { 'shape'  'A'     'r'     'rcap'  };
dflts =  {  []      1.0     1.0     1.2  };
[ shape, A, r, R,  ~] =...
     internal.stats.parseArgs(pnames, dflts, varargin{:});
    
if isempty(shape)
    disp('No shape specified; trying to proceed with the cube...')
    shape = 'cube';
end    

switch shape
    case 'cube'
        x(x> A/2.0) =  A - x(x> A/2.0);
        x(x<-A/2.0) = -A - x(x<-A/2.0);
    case 'sphere'
        rs = sum(x.*x,1);
        x(rs > r^2) =  r*x(rs > r^2)./sqrt(rs(rs > r^2));
    case 'shell'
        rs = sqrt(sum(x.*x,1));
        x = x ./ rs;
        rs = max(min(rs,R),r);
        x = rs .* x;
end

