function [X] = surface_nodes(k)
%Distribute n uniformily distributed Fibonacci nodes around the Earth using ETOPO1 data and linear
%interpolation between data points.

%Separation between the nodes is approximately sphere radius*3.0921*k^(-1/2)

fib = getFibonacciNodes(k);

[is, heights] = in_domain(fib(:,1),fib(:,2),fib(:,3));
heights(~is)=1;
X = repmat((heights)',1,3).*fib;
%  X = bsxfun(@times, heights', fib);
% hold on;
% pbaspect([1 1 1]);
% plot3(X(is,1),X(is,2),X(is,3), '.r', 'MarkerSize',1)
% plot3(X(~is,1),X(~is,2),X(~is,3), '.k', 'MarkerSize',1)
end