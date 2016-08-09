function r=slanttrui(xyz)

<<<<<<< HEAD
n = normr([1,1,1]);
proj = (xyz-dot(xyz,n)*n)';
=======
n = normr([1 1 0]);
>>>>>>> origin/HEAD
S = null(n);
proj = (xyz-dot(xyz,n)*n)';
proj = proj.*[sqrt(2); sqrt(2);1];

r = trui(frac_part(S\proj)');

%ixy = round(255*r);
%persistent F
%if isempty(F)
    % Read in the trui image once
    %A = double(imread('trui.png','PNG')); 
    %A = flipud(A(:,:,1));
    %F =  A/255; %Brightness    
%end


%f = F(1+ ixy(:,2)+ 256*ixy(:,1) );
