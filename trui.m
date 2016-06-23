
function r = trui(xy)
persistent F
% Return grain radius at location (x,y)
if isempty(F)
    % Read in the trui image once
    A = double(imread('trui.png','PNG')); 
    A = flipud(A(:,:,1));
    rf = @(s) 0.002+0.006*s+0.012*s.^8; % Conversion of brightness to grain radius
    F =  A/255;
    
end
ixy = round(255*xy);
r = F(1+ ixy(:,2)+ 256*ixy(:,1) ); 

% Given a location (x,y), evaluate
% the corresponding grain radius
