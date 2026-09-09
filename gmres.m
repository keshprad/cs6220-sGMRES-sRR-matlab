function [x,info] = gmres(A,b,d,options)
%GMRES Fixed-dimension GMRES using a full Arnoldi basis.
%
%   [X,INFO] = GMRES(A,B,D) performs exactly D Arnoldi steps unless a
%   happy breakdown occurs. A may be a square matrix or a function handle.
%
%   [X,INFO] = GMRES(...,InitialGuess=X0) starts from X0.

arguments
    A
    b (:,1) {mustBeNumeric,mustBeFinite}
    d (1,1) double {mustBeInteger,mustBePositive}
    options.InitialGuess = []
end

if ~isfloat(b)
    error("gmres:InvalidRightHandSide", ...
        "b must be single or double precision.");
end
n = numel(b);
if d > n
    error("gmres:InvalidDimension", ...
        "d must not exceed the ambient dimension.");
end
apply = krylov.asOperator(A,n);
x0 = initialGuess(options.InitialGuess,b);
r0 = b-apply(x0);
beta = norm(r0);
denominator = max(norm(b),eps(class(real(b))));

if beta == 0
    x = x0;
    info = struct( ...
        Dimension=0, ...
        Breakdown=true, ...
        RelativeResidual=0, ...
        EstimatedRelativeResidual=0, ...
        BasisOrthogonality=0);
    return
end

[V,H,arnoldiInfo] = krylov.arnoldi(apply,r0,d);
m = arnoldiInfo.Dimension;
smallRightHandSide = zeros(size(H,1),1,"like",H);
smallRightHandSide(1) = beta;
y = H\smallRightHandSide;
x = x0+V(:,1:m)*y;

trueResidual = b-apply(x);
info = struct( ...
    Dimension=m, ...
    Breakdown=arnoldiInfo.Breakdown, ...
    RelativeResidual=norm(trueResidual)/denominator, ...
    EstimatedRelativeResidual=norm(smallRightHandSide-H*y)/denominator, ...
    BasisOrthogonality=arnoldiInfo.Orthogonality);
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
