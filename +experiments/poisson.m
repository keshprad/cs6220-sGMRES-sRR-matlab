function [A,b,details] = poisson(gridSize,seed)
%POISSON Build a sparse five-point Poisson system with a known solution.

arguments
    gridSize (1,1) double {mustBeInteger,mustBePositive}
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end

if gridSize < 2
    error("experiments:poisson:InvalidGridSize", ...
        "gridSize must be at least 2.");
end

A = gallery("poisson",gridSize);
matrixSize = gridSize^2;
stream = RandStream("mt19937ar",Seed=seed);
exact = randn(stream,matrixSize,1);
exact = exact/norm(exact);
b = A*exact;

details = struct( ...
    MatrixSize=matrixSize, ...
    GridSize=gridSize, ...
    Exact=exact);
end
