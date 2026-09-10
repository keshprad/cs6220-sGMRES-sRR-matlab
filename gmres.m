function [x,info] = gmres(A,b,k,options)
%GMRES Fixed-dimension GMRES using a full Arnoldi basis.
%
%   [X,INFO] = GMRES(A,B,K) performs exactly K Arnoldi steps unless a
%   happy breakdown occurs. A may be a square matrix or a function handle.
%
%   [X,INFO] = GMRES(...,InitialGuess=X0) starts from X0.

arguments
    A
    b (:,1) {mustBeNumeric,mustBeFinite}
    k (1,1) double {mustBeInteger,mustBePositive}
    options.InitialGuess = []
end

totalTimer = tic;
if ~isfloat(b)
    error("gmres:InvalidRightHandSide", ...
        "b must be single or double precision.");
end
n = numel(b);
if k > n
    error("gmres:InvalidDimension", ...
        "k must not exceed the ambient dimension.");
end
apply = krylov.asOperator(A,n);
x0 = initialGuess(options.InitialGuess,b);
r0 = b-apply(x0);
beta = norm(r0);
denominator = max(norm(b),eps(class(real(b))));

if beta == 0
    x = x0;
    coreTime = toc(totalTimer);
    info = struct( ...
        Dimension=0, ...
        Breakdown=true, ...
        RelativeResidual=0, ...
        EstimatedRelativeResidual=0, ...
        BasisOrthogonality=0, ...
        CoreTime=coreTime, ...
        WallTime=toc(totalTimer));
    return
end

[V,H,arnoldiInfo] = krylov.arnoldi(apply,r0,k);
m = arnoldiInfo.Dimension;
smallRightHandSide = zeros(size(H,1),1,"like",H);
smallRightHandSide(1) = beta;
y = H\smallRightHandSide;
x = x0+V(:,1:m)*y;
coreTime = toc(totalTimer);

trueResidual = b-apply(x);
info = struct( ...
    Dimension=m, ...
    Breakdown=arnoldiInfo.Breakdown, ...
    RelativeResidual=norm(trueResidual)/denominator, ...
    EstimatedRelativeResidual=norm(smallRightHandSide-H*y)/denominator, ...
    BasisOrthogonality=arnoldiInfo.Orthogonality, ...
    CoreTime=coreTime, ...
    WallTime=toc(totalTimer));
end

function x0 = initialGuess(candidate,b)
if isempty(candidate)
    x0 = zeros(size(b),"like",b);
elseif isnumeric(candidate) && isfloat(candidate) && ...
        isequal(size(candidate),size(b)) && all(isfinite(candidate),"all")
    x0 = candidate;
else
    error("gmres:InvalidInitialGuess", ...
        "InitialGuess must be a finite floating-point vector the size of b.");
end
end
