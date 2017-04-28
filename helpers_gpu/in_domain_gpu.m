function is = in_domain_gpu(x, y, z)

% G_cnf = gpuArray( cnf );
% cnf = gather(G_cnf);
%persistent Z;

persistent Z_gpu;
if isempty(Z_gpu)
    load('z_transp.mat');    
end

 Z_gpu = gpuArray( Z );
 x_gpu = gpuArray( reshape(x,[],1) );
 y_gpu = gpuArray( reshape(y,[],1) );
 z_gpu = gpuArray( reshape(z,[],1) );

  %x_gpu = reshape(x,[],1);
  %y_gpu = reshape(y,[],1);
  %z_gpu = reshape(z,[],1);


[m, n] = size(Z);
delta = pi/n;
[p1, p2, p3] = cart2sph(x_gpu,y_gpu,z_gpu);

azdelta = mod(p1,delta);
eldelta = mod(p2,delta);

gridaz = uint16( (p1-azdelta)/delta + m/2 );
gridel = uint16( (p2-eldelta)/delta + n/2 );
gridind = sub2ind(size(Z), gridaz, gridel);

R = 1+[Z(gridind) Z(gridind+m) Z(gridind+1) Z(gridind+m+1)]/6378000;


indices =  (azdelta + eldelta) <= delta;

    r_interpolated = R(indices,1) + ( R(indices,2)-R(indices,1) ) .*...
        eldelta(indices,1)/delta +...
        ( R(indices,2)-R(indices,1) ) .* azdelta(indices,1)/delta;
% else 
    r_not_interpolated = R(~indices,4) + (R(~indices,2)-R(~indices,4)).*...
        (1-azdelta(~indices,1)/delta) +...
        (R(~indices,3)-R(~indices,4)).*(1-eldelta(~indices,1)/delta);
% end
is_gpu = zeros(1,numel(x_gpu),'gpuArray');
%is_gpu =zeros(1,numel(x));

is_gpu(indices) =  (p3(indices,1) < 1.1) .* (r_interpolated < p3(indices,1));
is_gpu(~indices) =  (p3(~indices,1) < 1.1) .* (r_not_interpolated < p3(~indices,1));
is = gather( is_gpu );
