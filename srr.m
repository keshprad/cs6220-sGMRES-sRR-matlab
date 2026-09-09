function [theta,U,info] = srr(A,v1,d,options)
%SRR Sketched Rayleigh--Ritz extraction from a truncated Arnoldi basis.
%
%   [THETA,U,INFO] = SRR(A,V1,D) uses truncation 2 and a DCT-IV sketch
%   with MIN(4*D,NUMEL(V1)) rows.
%
%   Name-value options are SketchSize, Truncation, NumEigenpairs, and Seed.

arguments
    A
    v1 (:,1) {mustBeNumeric,mustBeFinite}
    d (1,1) double {mustBeInteger,mustBePositive}
    options.SketchSize = []
    options.Truncation (1,1) double {mustBeInteger,mustBePositive} = 2
    options.NumEigenpairs = []
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end

if ~isfloat(v1)
    error("srr:InvalidStart", ...
        "v1 must be single or double precision.");
end
n = numel(v1);
if d > n
    error("srr:InvalidDimension", ...
        "d must not exceed the ambient dimension.");
end
sketchSize = validatedSketchSize(options.SketchSize,d,n);
requested = validatedCount(options.NumEigenpairs,d);

[V,AV,arnoldiInfo] = krylov.truncatedArnoldi( ...
    A,v1,d,options.Truncation);
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
    Truncation=options.Truncation, ...
    Seed=options.Seed, ...
    ReducedCondition=cond(R));
end

function sketchSize = validatedSketchSize(candidate,d,n)
if isempty(candidate)
    sketchSize = min(4*d,n);
elseif isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && ...
        candidate == fix(candidate)
    sketchSize = double(candidate);
else
    error("srr:InvalidSketchSize", ...
        "SketchSize must be an integer between d and numel(v1).");
end
if sketchSize < d || sketchSize > n
    error("srr:InvalidSketchSize", ...
        "SketchSize must be an integer between d and numel(v1).");
end
end

function count = validatedCount(candidate,d)
if isempty(candidate)
    count = min(5,d);
elseif isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && ...
        candidate == fix(candidate) && candidate >= 1 && candidate <= d
    count = double(candidate);
else
    error("srr:InvalidEigenpairCount", ...
        "NumEigenpairs must be an integer between 1 and d.");
end
end
