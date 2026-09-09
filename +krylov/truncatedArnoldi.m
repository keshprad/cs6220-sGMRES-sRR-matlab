function [V,AV,info] = truncatedArnoldi(A,v1,d,truncation)
%TRUNCATEDARNOLDI Build a Krylov basis using recent-vector orthogonalization.

arguments
    A
    v1 (:,1) {mustBeNumeric,mustBeFinite}
    d (1,1) double {mustBeInteger,mustBePositive}
    truncation (1,1) double {mustBeInteger,mustBePositive}
end

if ~isfloat(v1)
    error("krylov:truncatedArnoldi:InvalidStart", ...
        "v1 must be single or double precision.");
end
n = numel(v1);
if d > n
    error("krylov:truncatedArnoldi:InvalidDimension", ...
        "d must not exceed the ambient dimension.");
end
beta = norm(v1);
if beta == 0
    error("krylov:truncatedArnoldi:ZeroStart", ...
        "v1 must be nonzero.");
end

apply = krylov.asOperator(A,n);
V = zeros(n,d,"like",v1);
AV = zeros(n,d,"like",v1);
V(:,1) = v1/beta;
breakdown = false;
m = 0;

for j = 1:d
    image = apply(V(:,j));
    AV(:,j) = image;
    m = j;
    if j == d
        break
    end

    w = image;
    first = max(1,j-truncation+1);
    active = V(:,first:j);
    for pass = 1:2
        w = w - active*(active'*w);
    end

    threshold = 10*eps(class(real(w)))*max(1,norm(image));
    remainderNorm = norm(w);
    if remainderNorm <= threshold
        breakdown = true;
        break
    end
    V(:,j+1) = w/remainderNorm;
end

V = V(:,1:m);
AV = AV(:,1:m);
info = struct( ...
    Dimension=m, ...
    Breakdown=breakdown, ...
    Truncation=min(truncation,m));
end
