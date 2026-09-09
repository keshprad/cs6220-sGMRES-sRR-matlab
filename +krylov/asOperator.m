function apply = asOperator(A,n)
%ASOPERATOR Return a matrix-vector function for a matrix or function handle.

arguments
    A
    n (1,1) double {mustBeInteger,mustBePositive}
end

if isnumeric(A)
    if ~isfloat(A) || ~ismatrix(A) || ~isequal(size(A),[n n])
        error("krylov:asOperator:InvalidMatrix", ...
            "A must be an n-by-n floating-point matrix.");
    end
    if any(~isfinite(nonzeros(A)))
        error("krylov:asOperator:NonfiniteMatrix", ...
            "A must contain only finite values.");
    end
    apply = @(X) A*X;
elseif isa(A,"function_handle")
    apply = @(X) applyFunction(A,X,n);
else
    error("krylov:asOperator:InvalidOperator", ...
        "A must be a numeric matrix or a function handle.");
end
end
function Y = applyFunction(A,X,n)
Y = A(X);
if ~isnumeric(Y) || ~isfloat(Y) || ...
        size(Y,1) ~= n || size(Y,2) ~= size(X,2)
    error("krylov:asOperator:InvalidOutput", ...
        "The operator output must have the same size as its input.");
end
if any(~isfinite(Y),"all")
    error("krylov:asOperator:NonfiniteOutput", ...
        "The operator returned a nonfinite value.");
end
end
