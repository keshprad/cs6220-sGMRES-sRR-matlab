function raw = collectRaw(cfg)
%COLLECTRAW Run balanced paired timings and return one row per method call.

arguments
    cfg (1,1) struct
end

if ~license("test","Signal_Toolbox")
    error("experiments:collectRaw:MissingSignalToolbox", ...
        "The Signal Processing Toolbox is required for DCT sketches.");
end

rows = repmat(emptyRow(),0,1);
for suiteIndex = 1:numel(cfg.Suites)
    suite = cfg.Suites(suiteIndex);
    for gridIndex = 1:numel(cfg.GridSizes)
        gridSize = cfg.GridSizes(gridIndex);
        problemSeed = derivedSeed(cfg.ProblemSeed,gridSize);
        problem = buildProblem(suite,gridSize,problemSeed);

        for dimensionIndex = 1:numel(cfg.Dimensions)
            dimension = cfg.Dimensions(dimensionIndex);
            for sketchTrial = 1:cfg.SketchTrials
                sketchSeed = derivedSeed( ...
                    cfg.SketchSeed,gridSize,dimension,sketchTrial);
                methods = methodNames(suite);

                if dimension > problem.MatrixSize
                    for repetition = 1:cfg.Repetitions
                        order = orderedMethods(methods,suiteIndex,gridIndex, ...
                            dimensionIndex,sketchTrial,repetition);
                        for position = 1:2
                            row = baseRow(problem,dimension,order(position), ...
                                sketchTrial,repetition,position,problemSeed, ...
                                sketchSeed,cfg);
                            row.Status = "skipped-dimension";
                            rows(end+1,1) = row; %#ok<AGROW>
                        end
                    end
                    continue
                end

                warmMethods(problem,dimension,sketchSeed,cfg,methods);
                for repetition = 1:cfg.Repetitions
                    order = orderedMethods(methods,suiteIndex,gridIndex, ...
                        dimensionIndex,sketchTrial,repetition);
                    for position = 1:2
                        row = runOne(problem,dimension,order(position), ...
                            sketchTrial,repetition,position,problemSeed, ...
                            sketchSeed,cfg);
                        rows(end+1,1) = row; %#ok<AGROW>
                    end
                end
            end
        end
    end
end

raw = struct2table(rows);
end

function problem = buildProblem(suite,gridSize,problemSeed)
if suite == "gmres"
    [A,b,details] = experiments.poisson(gridSize,problemSeed);
    problem = struct( ...
        Suite=suite, ...
        Name="poisson-2d", ...
        GridSize=gridSize, ...
        MatrixSize=details.MatrixSize, ...
        MeshPeclet=NaN, ...
        A=A, ...
        B=b, ...
        Start=[], ...
        Oracle=emptyOracle());
else
    [A,~,details] = experiments.poisson(gridSize,problemSeed);
    v1 = details.Exact;
    problem = struct( ...
        Suite=suite, ...
        Name="poisson-2d", ...
        GridSize=gridSize, ...
        MatrixSize=details.MatrixSize, ...
        MeshPeclet=NaN, ...
        A=A, ...
        B=[], ...
        Start=v1, ...
        Oracle=computeOracle(A,v1));
end
end

function oracle = computeOracle(A,v1)
oracle = emptyOracle();
try
    count = min(6,numel(v1)-2);
    options = struct(tol=1e-12,v0=v1, ...
        maxit=max(1000,10*numel(v1)));
    [vectors,diagonal,flag] = eigs(A,count,"largestreal",options);
    if flag ~= 0
        error("experiments:collectRaw:OracleDidNotConverge", ...
            "eigs returned convergence flag %d.",flag);
    end
    values = diag(diagonal);
    residuals = A*vectors-vectors.*values.';
    relativeResiduals = vecnorm(residuals)./vecnorm(vectors);
    largestReal = max(real(values));
    tolerance = 1e-12+1e-9*abs(largestReal);
    cluster = abs(real(values)-largestReal) <= tolerance;
    oracle.Available = true;
    oracle.Values = values(cluster);
    oracle.Residual = max(relativeResiduals(cluster));
catch exception
    oracle.ErrorIdentifier = string(exception.identifier);
    oracle.ErrorMessage = string(exception.message);
end
end

function oracle = emptyOracle()
oracle = struct( ...
    Available=false, ...
    Values=complex(zeros(0,1)), ...
    Residual=NaN, ...
    ErrorIdentifier="", ...
    ErrorMessage="");
end

function methods = methodNames(suite)
if suite == "gmres"
    methods = ["GMRES" "sGMRES"];
else
    methods = ["RR" "sRR"];
end
end

function order = orderedMethods(methods,suiteIndex,gridIndex, ...
        dimensionIndex,sketchTrial,repetition)
cellParity = mod(suiteIndex+gridIndex+dimensionIndex+sketchTrial,2);
if mod(repetition+cellParity,2) == 0
    order = methods;
else
    order = fliplr(methods);
end
end

function warmMethods(problem,dimension,sketchSeed,cfg,methods)
for warmup = 1:cfg.Warmups
    for method = methods
        try
            invoke(problem,dimension,method,sketchSeed,cfg);
        catch
            % The recorded calls retain any repeatable method error.
        end
    end
end
end

function row = runOne(problem,dimension,method,sketchTrial,repetition, ...
        runOrder,problemSeed,sketchSeed,cfg)
row = baseRow(problem,dimension,method,sketchTrial,repetition, ...
    runOrder,problemSeed,sketchSeed,cfg);
externalTimer = tic;
try
    outcome = invoke(problem,dimension,method,sketchSeed,cfg);
    row.WallTime = toc(externalTimer);
    info = outcome.Info;
    row.AchievedDimension = info.Dimension;
    row.CoreTime = info.CoreTime;
    row.InstrumentedWallTime = info.WallTime;
    row.RelativeResidual = outcome.RelativeResidual;
    row.SketchedRelativeResidual = outcome.SketchedRelativeResidual;
    row.RitzResidual = outcome.RitzResidual;
    row.SketchedRitzResidual = outcome.SketchedRitzResidual;
    row.ReducedCondition = outcome.ReducedCondition;
    row.RitzReal = real(outcome.RitzValue);
    row.RitzImag = imag(outcome.RitzValue);
    if problem.Oracle.Available && isfinite(row.RitzReal)
        row.OracleDistance = min(abs( ...
            outcome.RitzValue-problem.Oracle.Values));
    end
    row.Status = "ok";
catch exception
    row.WallTime = toc(externalTimer);
    row.Status = "error";
    row.ErrorIdentifier = string(exception.identifier);
    row.ErrorMessage = string(exception.message);
end
end

function outcome = invoke(problem,dimension,method,sketchSeed,cfg)
outcome = emptyOutcome();
switch method
    case "GMRES"
        [~,info] = gmres(problem.A,problem.B,dimension);
        outcome.Info = info;
        outcome.RelativeResidual = info.RelativeResidual;
    case "sGMRES"
        [~,info] = sgmres(problem.A,problem.B,dimension, ...
            Truncation=cfg.GmresTruncation,Seed=sketchSeed);
        outcome.Info = info;
        outcome.RelativeResidual = info.RelativeResidual;
        outcome.SketchedRelativeResidual = info.SketchedRelativeResidual;
        outcome.ReducedCondition = info.ReducedCondition;
    case "RR"
        [theta,~,info] = rr( ...
            problem.A,problem.Start,dimension,NumEigenpairs=1);
        outcome.Info = info;
        outcome.RitzValue = theta(1);
        outcome.RitzResidual = info.RitzResiduals(1);
    case "sRR"
        [theta,~,info] = srr(problem.A,problem.Start,dimension, ...
            Truncation=cfg.RrTruncation,NumEigenpairs=1,Seed=sketchSeed);
        outcome.Info = info;
        outcome.RitzValue = theta(1);
        outcome.RitzResidual = info.RitzResiduals(1);
        outcome.SketchedRitzResidual = info.SketchedRitzResiduals(1);
        outcome.ReducedCondition = info.ReducedCondition;
end
end

function outcome = emptyOutcome()
outcome = struct( ...
    Info=struct(), ...
    RelativeResidual=NaN, ...
    SketchedRelativeResidual=NaN, ...
    RitzResidual=NaN, ...
    SketchedRitzResidual=NaN, ...
    ReducedCondition=NaN, ...
    RitzValue=complex(NaN,NaN));
end

function row = baseRow(problem,dimension,method,sketchTrial,repetition, ...
        runOrder,problemSeed,sketchSeed,cfg)
row = emptyRow();
row.Suite = problem.Suite;
row.Problem = problem.Name;
row.GridSize = problem.GridSize;
row.MatrixSize = problem.MatrixSize;
row.MeshPeclet = problem.MeshPeclet;
row.Dimension = dimension;
row.Method = method;
row.SketchTrial = sketchTrial;
row.Repetition = repetition;
row.RunOrder = runOrder;
row.ProblemSeed = problemSeed;
row.SketchSeed = sketchSeed;
row.OracleAvailable = problem.Oracle.Available;
row.OracleResidual = problem.Oracle.Residual;
if startsWith(method,"s")
    if problem.Suite == "gmres"
        row.SketchSize = min(2*(dimension+1),problem.MatrixSize);
        row.Truncation = cfg.GmresTruncation;
    else
        row.SketchSize = min(4*dimension,problem.MatrixSize);
        row.Truncation = cfg.RrTruncation;
    end
end
end

function row = emptyRow()
row = struct( ...
    Suite="", ...
    Problem="", ...
    GridSize=NaN, ...
    MatrixSize=NaN, ...
    MeshPeclet=NaN, ...
    Dimension=NaN, ...
    Method="", ...
    SketchTrial=NaN, ...
    Repetition=NaN, ...
    RunOrder=NaN, ...
    ProblemSeed=NaN, ...
    SketchSeed=NaN, ...
    SketchSize=NaN, ...
    Truncation=NaN, ...
    AchievedDimension=NaN, ...
    CoreTime=NaN, ...
    InstrumentedWallTime=NaN, ...
    WallTime=NaN, ...
    RelativeResidual=NaN, ...
    SketchedRelativeResidual=NaN, ...
    RitzResidual=NaN, ...
    SketchedRitzResidual=NaN, ...
    ReducedCondition=NaN, ...
    RitzReal=NaN, ...
    RitzImag=NaN, ...
    OracleDistance=NaN, ...
    OracleResidual=NaN, ...
    OracleAvailable=false, ...
    Status="", ...
    ErrorIdentifier="", ...
    ErrorMessage="");
end

function seed = derivedSeed(baseSeed,varargin)
modulus = 2^32;
seed = mod(double(baseSeed),modulus);
for index = 1:numel(varargin)
    coordinate = double(varargin{index});
    seed = mod(seed*1664525+coordinate*1013904223+1013904223,modulus);
end
seed = floor(seed);
end
