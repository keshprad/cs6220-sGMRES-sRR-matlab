function [theta,U,info] = srr(A,v1,k,options)
%SRR Sketched Rayleigh--Ritz extraction from a truncated Arnoldi basis.
%
%   [THETA,U,INFO] = SRR(A,V1,K) uses truncation length ell=2 and a DCT-IV
%   sketch with MIN(4*K,NUMEL(V1)) rows.
%
%   Name-value options are SketchSize, Truncation, NumEigenpairs, and Seed.

arguments
    A
    v1 (:,1) {mustBeNumeric,mustBeFinite}
    k (1,1) double {mustBeInteger,mustBePositive}
    options.SketchSize = []
    options.Truncation (1,1) double {mustBeInteger,mustBePositive} = 2
    options.NumEigenpairs = []
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end

totalTimer = tic;
if ~isfloat(v1)
    error("srr:InvalidStart", ...
        "v1 must be single or double precision.");
end
n = numel(v1);
if k > n
    error("srr:InvalidDimension", ...
        "k must not exceed the ambient dimension.");
end
ell = options.Truncation;
sketchSize = validatedSketchSize(options.SketchSize,k,n);
requested = validatedCount(options.NumEigenpairs,k);

[V,AV,arnoldiInfo] = krylov.truncatedArnoldi( ...
    A,v1,k,ell);
m = arnoldiInfo.Dimension;
sketch = krylov.srdct(n,sketchSize,4,options.Seed);
sketchedBasis = sketch.Apply(V);
sketchedImages = sketch.Apply(AV);
[Q,R,p] = krylov.pivotedQR(sketchedBasis);
permutedProjection = R\(Q'*sketchedImages);
projected = zeros(m,m,"like",permutedProjection);
projected(p,:) = permutedProjection;

count = min(requested,m);
[theta,U,Y] = krylov.ritzPairs(projected,V,count);
coreTime = toc(totalTimer);
residuals = AV*Y-U.*theta.';
sketchedVectors = sketchedBasis*Y;
sketchedResiduals = sketchedImages*Y-sketchedVectors.*theta.';
info = struct( ...
    Dimension=m, ...
    Breakdown=arnoldiInfo.Breakdown, ...
    RitzResiduals=vecnorm(residuals)./vecnorm(U), ...
    SketchedRitzResiduals= ...
        vecnorm(sketchedResiduals)./vecnorm(sketchedVectors), ...
    SketchSize=sketchSize, ...
    Truncation=ell, ...
    Seed=options.Seed, ...
    ReducedCondition=cond(R), ...
    CoreTime=coreTime, ...
    WallTime=toc(totalTimer));
end

function sketchSize = validatedSketchSize(candidate,k,n)
if isempty(candidate)
    sketchSize = min(4*k,n);
elseif isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && ...
        candidate == fix(candidate)
    sketchSize = double(candidate);
else
    error("srr:InvalidSketchSize", ...
        "SketchSize must be an integer between k and numel(v1).");
end
if sketchSize < k || sketchSize > n
    error("srr:InvalidSketchSize", ...
        "SketchSize must be an integer between k and numel(v1).");
end
end

function count = validatedCount(candidate,k)
if isempty(candidate)
    count = min(5,k);
elseif isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && ...
        candidate == fix(candidate) && candidate >= 1 && candidate <= k
    count = double(candidate);
else
    error("srr:InvalidEigenpairCount", ...
        "NumEigenpairs must be an integer between 1 and k.");
end
end
