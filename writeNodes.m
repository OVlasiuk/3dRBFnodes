function X = writeNodes(X, box, index)

l = size(box,2)+1;
X(:,1,index) = X(:,1,index) + size(box,2);
X(:,2:l,index) = box;
