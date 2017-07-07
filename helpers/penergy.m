function [en, grad]  = penergy(cnf,s)
%PENERGY
% en  = penergy(cnf,s)
% Periodic Riesz energy in the unit 3-dimensional cube; accepts short fat
% matrices.
if ~exist('s', 'var')
    s = 4.0;
end
if size(cnf,1) ~= 3
    cnf = reshape(cnf,3,[]);
    tr_back = 1;
else
    tr_back = 0;
end
switch s
    case 4.0
        compute_riesz = @(x) 1./x./x;
        compute_grad = @(x) 1./x./x./x;
    case 2.0
        compute_riesz = @(x) 1./x;
        compute_grad = @(x) 1./x./x;
    case 0.5
        compute_riesz = @(x) 1./sqrt(sqrt(x));
        compute_grad = @(x) 1./sqrt(sqrt(x))./x;
    otherwise
        compute_riesz = @(x) sqrt(x).^(-s);
        compute_grad = @(x) sqrt(x).^(-s-1);
end
N = numel(cnf)/3;       % number of 3-vectors
D = diag(true(1,N));

cnfd1 = cnf(1,:)-cnf(1,:)';
cnf_d1a = abs(cnfd1);
[cnf_d1a, idx1] = min([cnf_d1a(:),  1-cnf_d1a(:)],[],2);
cnf_d1a = reshape(cnf_d1a,[],N);
gr1 = [cnfd1(:), cnfd1(:)-sign(cnfd1(:))];
gr1 = reshape( gr1(sub2ind(size(gr1),1:N^2,idx1')), [], N );

cnfd2 = cnf(2,:)-cnf(2,:)';
cnf_d2a = abs(cnfd2);
[cnf_d2a, idx2] = min([cnf_d2a(:),  1-cnf_d2a(:)],[],2);
cnf_d2a = reshape(cnf_d2a,[],N);
gr2 = [cnfd2(:), cnfd2(:)-sign(cnfd2(:))];
gr2 = reshape( gr2(sub2ind(size(gr2),1:N^2,idx2')), [], N );

cnfd3 = cnf(3,:)-cnf(3,:)';
cnf_d3a = abs(cnfd3);
[cnf_d3a, idx3] = min([cnf_d3a(:),  1-cnf_d3a(:)],[],2);
cnf_d3a = reshape(cnf_d3a,[],N);
gr3 = [cnfd3(:), cnfd3(:)-sign(cnfd3(:))];
gr3 = reshape( gr3(sub2ind(size(gr3),1:N^2,idx3')), [], N );

en1      =  cnf_d1a.*cnf_d1a +...
            cnf_d2a.*cnf_d2a +...
            cnf_d3a.*cnf_d3a...
            ;
en = compute_riesz(en1);

grad_weight = compute_grad(en1);
grad = -s * [sum(gr1.*grad_weight,1,'omitnan');...
            sum(gr2.*grad_weight,1,'omitnan');...
            sum(gr3.*grad_weight,1,'omitnan')];
en = sum(sum(en( ~D )));
if tr_back
    grad = grad(:);
end