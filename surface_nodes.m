function [X] = surface_nodes(k)
%Distribute n uniformily distributed Fibonacci nodes around the Earth using ETOPO1 data and linear
%interpolation between data points.

%Separation between the nodes is approximately sphere radius*3.0921*k^(-1/2)

Y = getFibonacciNodes(k);

persistent Z;
if isempty(Z)
    load('z_transp.mat');    
end

inner = .9;

[m, n] = size(Z);
nn = numel(Z);
delta = pi/n;
[p1, p2, p3] = cart2sph(Y(:,1),Y(:,2),Y(:,3));

azdelta = mod(p1,delta);
eldelta = mod(p2,delta);

gridaz = uint16( (p1-azdelta)/delta + m/2 +1);
gridel = uint16( (p2-eldelta)/delta + n/2 +1);
gridind = sub2ind(size(Z), gridaz, gridel);

% Instead of dividing by the Earth radius in meters, i.e. 6,378,000, we exaggerate the
% radial coordinate.
R = 0.9+[Z(gridind) Z(mod(gridind+m,nn)) Z(mod(gridind+1,nn)) Z(mod(gridind+m+1,nn))]/63780;


indices =  (azdelta + eldelta) <= delta;

    r_interpolated_lower = R(indices,1) + ( R(indices,2)-R(indices,1) ) .*...
        eldelta(indices,1)/delta +...
        ( R(indices,2)-R(indices,1) ) .* azdelta(indices,1)/delta;
    
    r_interpolated_upper = R(~indices,4) + (R(~indices,2)-R(~indices,4)).*...
        (1-azdelta(~indices,1)/delta) +...
        (R(~indices,3)-R(~indices,4)).*(1-eldelta(~indices,1)/delta);

heights(indices) =  (r_interpolated_lower > inner) .* r_interpolated_lower;
heights(~indices) =  (r_interpolated_upper > inner) .* r_interpolated_upper;

X = [heights',heights',heights'].*Y;
end