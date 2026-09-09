function sketch = srdct(n,sketchSize,dctType,seed)
%SRDCT Construct a subsampled randomized discrete cosine transform.

arguments
    n (1,1) double {mustBeInteger,mustBePositive}
    sketchSize (1,1) double {mustBeInteger,mustBePositive}
    dctType (1,1) double {mustBeInteger}
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end

if sketchSize > n
    error("krylov:srdct:InvalidSketchSize", ...
        "sketchSize must lie between 1 and n.");
end
if ~ismember(dctType,[2 4])
    error("krylov:srdct:InvalidDctType", ...
        "dctType must be 2 or 4.");
end

stream = RandStream("mt19937ar",Seed=seed);
signs = 2*(rand(stream,n,1) >= 0.5)-1;
rows = randperm(stream,n,sketchSize).';
scale = sqrt(n/sketchSize);

sketch = struct( ...
    InputSize=n, ...
    SketchSize=sketchSize, ...
    Type=dctType, ...
    Seed=seed, ...
    Rows=rows, ...
    Signs=signs, ...
    Scale=scale, ...
    Apply=@(X) applySketch(X,signs,rows,scale,dctType,n));
end

function Y = applySketch(X,signs,rows,scale,dctType,n)
if ~isnumeric(X) || ~isfloat(X) || size(X,1) ~= n
    error("krylov:srdct:InvalidInput", ...
        "The sketched input must be a floating-point array with n rows.");
end
if any(~isfinite(X),"all")
    error("krylov:srdct:NonfiniteInput", ...
        "The sketched input must contain only finite values.");
end
transformed = dct(signs.*X,[],1,Type=dctType);
Y = scale*transformed(rows,:);
end
