function [x,info] = sgmres(A,b,k,options)
%SGMRES Fixed-dimension sketched GMRES using truncated Arnoldi.
%
%   [X,INFO] = SGMRES(A,B,K) uses truncation length ell=4 and a DCT-II
%   sketch with MIN(2*(K+1),NUMEL(B)) rows.
%
%   Name-value options are InitialGuess, SketchSize, Truncation, and Seed.

arguments
    A
    b (:,1) {mustBeNumeric,mustBeFinite}
    k (1,1) double {mustBeInteger,mustBePositive}
    options.InitialGuess = []
    options.SketchSize = []
    options.Truncation (1,1) double {mustBeInteger,mustBePositive} = 4
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end

totalTimer = tic;
if ~isfloat(b)
    error("sgmres:InvalidRightHandSide", ...
        "b must be single or double precision.");
end
n = numel(b);
if k > n
    error("sgmres:InvalidDimension", ...
        "k must not exceed the ambient dimension.");
end
ell = options.Truncation;
sketchSize = validatedSketchSize(options.SketchSize,k,n);
apply = krylov.asOperator(A,n);
x0 = initialGuess(options.InitialGuess,b);
r0 = b-apply(x0);
denominator = max(norm(b),eps(class(real(b))));

if norm(r0) == 0
    x = x0;
    coreTime = toc(totalTimer);
    info = struct( ...
        Dimension=0, ...
        Breakdown=true, ...
        RelativeResidual=0, ...
        SketchedRelativeResidual=0, ...
        SketchSize=sketchSize, ...
        Truncation=ell, ...
        Seed=options.Seed, ...
        ReducedCondition=1, ...
        CoreTime=coreTime, ...
        WallTime=toc(totalTimer));
    return
end

[V,AV,arnoldiInfo] = krylov.truncatedArnoldi( ...
    apply,r0,k,ell);
m = arnoldiInfo.Dimension;
sketch = krylov.srdct(n,sketchSize,2,options.Seed);
sketchedResidual = sketch.Apply(r0);
sketchedImages = sketch.Apply(AV);
[Q,R,p] = krylov.pivotedQR(sketchedImages);
permutedCoefficients = R\(Q'*sketchedResidual);
y = zeros(m,1,"like",permutedCoefficients);
y(p) = permutedCoefficients;
x = x0+V*y;
coreTime = toc(totalTimer);

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
    Truncation=ell, ...
    Seed=options.Seed, ...
    ReducedCondition=cond(R), ...
    CoreTime=coreTime, ...
    WallTime=toc(totalTimer));
end

function sketchSize = validatedSketchSize(candidate,k,n)
if isempty(candidate)
    sketchSize = min(2*(k+1),n);
elseif isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && ...
        candidate == fix(candidate)
    sketchSize = double(candidate);
else
    error("sgmres:InvalidSketchSize", ...
        "SketchSize must be an integer between k and numel(b).");
end
if sketchSize < k || sketchSize > n
    error("sgmres:InvalidSketchSize", ...
        "SketchSize must be an integer between k and numel(b).");
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
