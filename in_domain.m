function is = in_domain_gpu(x, y, z)

% G_neighbors = gpuArray( cnf(:,IDX) );
% G_cnf = gpuArray( cnf );
% 
% for iter=1:repel_steps k
%        cnf_repeated = reshape(repmat(G_cnf,k_value,1), dim, k_value*pt_num); 
%        directions = cnf_repeated - G_neighbors;
%        inverse_norms_riesz = sum(directions.^2,1).^(-0.5*(s+1));
%        directions = bsxfun(@times,inverse_norms_riesz,directions);
%        directions = sum(reshape(directions, dim, k_value, pt_num),2);
%        directions = reshape(directions, dim, pt_num);
%        inverse_norms = sum(directions.^2,1).^(-0.5);
%        forces =  bsxfun(@times,inverse_norms,directions); 
%        
%     G_cnf = G_cnf + forces*step/5/iter;
%     G_cnf(G_cnf<0) =  -G_cnf(G_cnf<0);
%     G_cnf(G_cnf>1) =  2-G_cnf(G_cnf>1);
%     G_neighbors = G_cnf(:,IDX);
% end
% 
% cnf = gather(G_cnf);
persistent Z;
if isempty(Z)
    load('z_transp.mat');    
end

% persistent Z_gpu = gpuArray( Z );
% x_gpu = gpuArray( x );
% y_gpu = gpuArray( y );
% z_gpu = gpuArray( z );

  x_gpu = reshape(x,[],1);
  y_gpu = reshape(y,[],1);
  z_gpu = reshape(z,[],1);


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
% is = zeros(numel(x_gpu),'gpuArray');
is =zeros(1,numel(x));

is(indices) =  (p3(indices,1) < 1.1) .* (r_interpolated < p3(indices,1));
is(~indices) =  (p3(~indices,1) < 1.1) .* (r_not_interpolated < p3(~indices,1));