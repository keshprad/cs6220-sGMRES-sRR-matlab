function [theta,U,info] = rr(A,v1,k,options)
%RR Rayleigh--Ritz extraction from a full Arnoldi basis.
%
%   [THETA,U,INFO] = RR(A,V1,K) returns up to five Ritz pairs whose Ritz
%   values have largest real part. A may be a matrix or a function handle.
%
%   [THETA,U,INFO] = RR(...,NumEigenpairs=Q) returns Q Ritz pairs, or all
%   available pairs if happy breakdown produces a smaller subspace.

arguments
    A
    v1 (:,1) {mustBeNumeric,mustBeFinite}
    k (1,1) double {mustBeInteger,mustBePositive}
    options.NumEigenpairs = []
end

totalTimer = tic;
if ~isfloat(v1)
    error("rr:InvalidStart", ...
        "v1 must be single or double precision.");
end
n = numel(v1);
if k > n
    error("rr:InvalidDimension", ...
        "k must not exceed the ambient dimension.");
end
requested = validatedCount(options.NumEigenpairs,k);

[V,H,arnoldiInfo] = krylov.arnoldi(A,v1,k);
m = arnoldiInfo.Dimension;
count = min(requested,m);
projected = H(1:m,1:m);
[theta,U,Y] = krylov.ritzPairs(projected,V(:,1:m),count);
coreTime = toc(totalTimer);

images = V*H;
residuals = images*Y-U.*theta.';
info = struct( ...
    Dimension=m, ...
    Breakdown=arnoldiInfo.Breakdown, ...
    RitzResiduals=vecnorm(residuals)./vecnorm(U), ...
    BasisOrthogonality=arnoldiInfo.Orthogonality, ...
    CoreTime=coreTime, ...
    WallTime=toc(totalTimer));
end

function count = validatedCount(candidate,k)
if isempty(candidate)
    count = min(5,k);
elseif isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && ...
        candidate == fix(candidate) && candidate >= 1 && candidate <= k
    count = double(candidate);
else
    error("rr:InvalidEigenpairCount", ...
        "NumEigenpairs must be an integer between 1 and k.");
end
end
