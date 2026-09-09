function [theta,U,Y] = ritzPairs(projected,V,numEigenpairs)
%RITZPAIRS Extract normalized Ritz pairs with largest real parts.

arguments
    projected {mustBeNumeric,mustBeFinite}
    V {mustBeNumeric,mustBeFinite}
    numEigenpairs (1,1) double {mustBeInteger,mustBePositive}
end

m = size(projected,1);
if ~isfloat(projected) || ~isequal(size(projected),[m m]) || ...
        ~isfloat(V) || size(V,2) ~= m
    error("krylov:ritzPairs:InvalidDimensions", ...
        "projected must be square and match the column count of V.");
end
if numEigenpairs > m
    error("krylov:ritzPairs:TooManyPairs", ...
        "numEigenpairs must not exceed the projected dimension.");
end

[allVectors,allValues] = eig(projected,"vector");
[~,order] = sort(real(allValues),"descend");
selected = order(1:numEigenpairs);
theta = allValues(selected);
Y = allVectors(:,selected);
U = V*Y;
scales = vecnorm(U);
Y = Y./scales;
U = V*Y;
end
