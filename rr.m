function [theta,U,info] = rr(A,v1,d,options)
%RR Rayleigh--Ritz extraction from a full Arnoldi basis.
%
%   [THETA,U,INFO] = RR(A,V1,D) returns up to five Ritz pairs whose Ritz
%   values have largest real part. A may be a matrix or a function handle.
%
%   [THETA,U,INFO] = RR(...,NumEigenpairs=Q) returns Q Ritz pairs, or all
%   available pairs if happy breakdown produces a smaller subspace.

arguments
    A
    v1 (:,1) {mustBeNumeric,mustBeFinite}
    d (1,1) double {mustBeInteger,mustBePositive}
    options.NumEigenpairs = []
end

if ~isfloat(v1)
    error("rr:InvalidStart", ...
        "v1 must be single or double precision.");
end
n = numel(v1);
if d > n
    error("rr:InvalidDimension", ...
        "d must not exceed the ambient dimension.");
end
requested = validatedCount(options.NumEigenpairs,d);

[V,H,arnoldiInfo] = krylov.arnoldi(A,v1,d);
m = arnoldiInfo.Dimension;
count = min(requested,m);
projected = H(1:m,1:m);
[theta,U,Y] = krylov.ritzPairs(projected,V(:,1:m),count);

images = V*H;
residuals = images*Y-U.*theta.';
info = struct( ...
    Dimension=m, ...
    Breakdown=arnoldiInfo.Breakdown, ...
    RitzResiduals=vecnorm(residuals)./vecnorm(U), ...
    BasisOrthogonality=arnoldiInfo.Orthogonality);
end

function count = validatedCount(candidate,d)
if isempty(candidate)
    count = min(5,d);
elseif isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && ...
        candidate == fix(candidate) && candidate >= 1 && candidate <= d
    count = double(candidate);
else
    error("rr:InvalidEigenpairCount", ...
        "NumEigenpairs must be an integer between 1 and d.");
end
end
