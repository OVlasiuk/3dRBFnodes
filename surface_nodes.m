function [X] = surface_nodes(k)
%Distribute n uniformily distributed Fibonacci nodes on the Earth bedrock
% using ETOPO1 data and linear interpolation between data points.

%Separation between the Fibonacci nodes is approximately sphere radius*3.0921*k^(-1/2)
fib = getFibonacciNodes(k);

[is, heights] = in_domain(fib(:,1),fib(:,2),fib(:,3));
% heights(~is)=1;
X = bsxfun(@times, heights',fib);

%% Uncomment the following lines to plot X.
%  clf;
%  hold on;
%  pbaspect([1 1 1]);
%  plot3(X(:,1),X(:,2),X(:,3), '.k', 'MarkerSize',1)
% plot3(X(is,1),X(is,2),X(is,3), '.k', 'MarkerSize',1)
% plot3(X(~is,1),X(~is,2),X(~is,3), '.b', 'MarkerSize',1)
end
