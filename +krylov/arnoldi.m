function [V,H,info] = arnoldi(A,v1,k)
%ARNOLDI Build an orthonormal Krylov basis with two Gram--Schmidt passes.

arguments
    A
    v1 (:,1) {mustBeNumeric,mustBeFinite}
    k (1,1) double {mustBeInteger,mustBePositive}
end

if ~isfloat(v1)
    error("krylov:arnoldi:InvalidStart", ...
        "v1 must be single or double precision.");
end
n = numel(v1);
if k > n
    error("krylov:arnoldi:InvalidDimension", ...
        "k must not exceed the ambient dimension.");
end
beta = norm(v1);
if beta == 0
    error("krylov:arnoldi:ZeroStart", ...
        "v1 must be nonzero.");
end

apply = krylov.asOperator(A,n);
V = zeros(n,k+1,"like",v1);
H = zeros(k+1,k,"like",v1);
V(:,1) = v1/beta;
breakdown = false;
m = 0;

for j = 1:k
    image = apply(V(:,j));
    w = image;
    active = V(:,1:j);
    for pass = 1:2
        coefficients = active'*w;
        H(1:j,j) = H(1:j,j) + coefficients;
        w = w - active*coefficients;
    end

    remainderNorm = norm(w);
    H(j+1,j) = remainderNorm;
    m = j;
    threshold = 10*eps(class(real(w)))*max(1,norm(image));
    if remainderNorm <= threshold
        breakdown = true;
        break
    end
    V(:,j+1) = w/remainderNorm;
end

if breakdown
    V = V(:,1:m);
    H = H(1:m,1:m);
else
    V = V(:,1:m+1);
    H = H(1:m+1,1:m);
end
activeBasis = V(:,1:m);
info = struct( ...
    Dimension=m, ...
    Breakdown=breakdown, ...
    Orthogonality=norm(activeBasis'*activeBasis-eye(m,"like",H),"fro"));
end
