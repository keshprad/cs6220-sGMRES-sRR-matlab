function [x,info] = sgmres(A,b,d,options)
%SGMRES Fixed-dimension sketched GMRES using truncated Arnoldi.
%
%   [X,INFO] = SGMRES(A,B,D) uses truncation 4 and a DCT-II sketch with
%   MIN(2*(D+1),NUMEL(B)) rows.
%
%   Name-value options are InitialGuess, SketchSize, Truncation, and Seed.

arguments
    A
    b (:,1) {mustBeNumeric,mustBeFinite}
    d (1,1) double {mustBeInteger,mustBePositive}
    options.InitialGuess = []
    options.SketchSize = []
    options.Truncation (1,1) double {mustBeInteger,mustBePositive} = 4
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end

if ~isfloat(b)
    error("sgmres:InvalidRightHandSide", ...
        "b must be single or double precision.");
end
n = numel(b);
if d > n
    error("sgmres:InvalidDimension", ...
        "d must not exceed the ambient dimension.");
end
sketchSize = validatedSketchSize(options.SketchSize,d,n);
apply = krylov.asOperator(A,n);
x0 = initialGuess(options.InitialGuess,b);
r0 = b-apply(x0);
denominator = max(norm(b),eps(class(real(b))));

if norm(r0) == 0
    x = x0;
    info = struct( ...
        Dimension=0, ...
        Breakdown=true, ...
        RelativeResidual=0, ...
        SketchedRelativeResidual=0, ...
        SketchSize=sketchSize, ...
        Truncation=options.Truncation, ...
        Seed=options.Seed, ...
        ReducedCondition=1);
    return
end

[V,AV,arnoldiInfo] = krylov.truncatedArnoldi( ...
    apply,r0,d,options.Truncation);
m = arnoldiInfo.Dimension;
sketch = krylov.srdct(n,sketchSize,2,options.Seed);
sketchedResidual = sketch.Apply(r0);
sketchedImages = sketch.Apply(AV);
[Q,R,p] = krylov.pivotedQR(sketchedImages);
permutedCoefficients = R\(Q'*sketchedResidual);
y = zeros(m,1,"like",permutedCoefficients);
y(p) = permutedCoefficients;
x = x0+V*y;

trueResidual = b-apply(x);
sketchedRightHandSide = sketch.Apply(b);
sketchedDenominator = max(norm(sketchedRightHandSide), ...
    eps(class(real(sketchedRightHandSide))));
info = struct( ...
    Dimension=m, ...
    Breakdown=arnoldiInfo.Breakdown, ...
    RelativeResidual=norm(trueResidual)/denominator, ...
    SketchedRelativeResidual= ...
        norm(sketchedResidual-sketchedImages*y)/sketchedDenominator, ...
    SketchSize=sketchSize, ...
    Truncation=options.Truncation, ...
    Seed=options.Seed, ...
    ReducedCondition=cond(R));
end

function sketchSize = validatedSketchSize(candidate,d,n)
if isempty(candidate)
    sketchSize = min(2*(d+1),n);
elseif isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && ...
        candidate == fix(candidate)
    sketchSize = double(candidate);
else
    error("sgmres:InvalidSketchSize", ...
        "SketchSize must be an integer between d and numel(b).");
end
if sketchSize < d || sketchSize > n
    error("sgmres:InvalidSketchSize", ...
        "SketchSize must be an integer between d and numel(b).");
end
end

function x0 = initialGuess(candidate,b)
if isempty(candidate)
    x0 = zeros(size(b),"like",b);
elseif isnumeric(candidate) && isfloat(candidate) && ...
        isequal(size(candidate),size(b)) && all(isfinite(candidate),"all")
    x0 = candidate;
else
    error("sgmres:InvalidInitialGuess", ...
        "InitialGuess must be a finite floating-point vector the size of b.");
end
end
