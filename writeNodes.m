function X = writeNodes(X, box, index)

l = size(box,1)+1;
X(index,1,:) = X(index, 1, :) + size(box,1);
X(index, 2:l, :) = box;
