function [Q,R,p] = pivotedQR(B)
%PIVOTEDQR Compute an economical pivoted QR and check numerical rank.

arguments
    B {mustBeNumeric,mustBeFinite}
end

if ~isfloat(B) || ~ismatrix(B) || size(B,1) < size(B,2) || isempty(B)
    error("krylov:pivotedQR:InvalidMatrix", ...
        "B must be a nonempty, tall floating-point matrix.");
end

[Q,R,p] = qr(B,"econ","vector");
diagonal = abs(diag(R));
threshold = eps(class(real(R)))*max(size(B))*max(diagonal);
if any(diagonal <= threshold)
    error("krylov:pivotedQR:RankDeficient", ...
        "The reduced matrix is numerically rank deficient.");
end
end
