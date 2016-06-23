
function r = trui(xy)
persistent F
% Return grain radius at location (x,y)
if isempty(F)
    % Read in the trui image once
    A = double(imread('trui.png','PNG')); 
    A = flipud(A(:,:,1));
    F =  A/255; %Brightness
    
end
ixy = round(255*xy);
r = (F(1+ ixy(:,2)+ 256*ixy(:,1) )); 

% Given a location (x,y), evaluate
% the corresponding grain radius
