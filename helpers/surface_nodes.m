function X = surface_nodes(k,plotit)
%SURFACE_NODES
% X = surface_nodes(k,plotit)
% Distribute k uniformily distributed Fibonacci nodes on the Earth bedrock
% using ETOPO1 data and linear interpolation between data points.
% k -- number of nodes;
% plotit -- pass 'y' or 1, etc., to plot the produced configuration.

%Separation between the Fibonacci nodes is approximately sphere radius*3.0921*k^(-1/2)
if ~exist('plotit','var')
    plotit = 0;
end
fib = getFibonacciNodes(k);

[is, heights] = in_domain(fib(:,1),fib(:,2),fib(:,3));
% heights(~is)=1;
X = bsxfun(@times, heights',fib);

msize = ceil(max(1, 22-4*log10(k) ));
if exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
    pbaspect([1 1 1]);
    colormap(winter)
    grid on;
    plot3(X(:,1),X(:,2),X(:,3), '.b', 'MarkerSize',msize)
    axis vis3d
end

%  plot3(X(:,1),X(:,2),X(:,3), '.k', 'MarkerSize',1)
% plot3(X(is,1),X(is,2),X(is,3), '.k', 'MarkerSize',1)
% plot3(X(~is,1),X(~is,2),X(~is,3), '.b', 'MarkerSize',1)

end
