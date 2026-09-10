function cfg = config(options)
%CONFIG Resolve and validate one crossover experiment configuration.

arguments
    options.Preset (1,1) string = "quick"
    options.Suites string = strings(1,0)
    options.GridSizes double = []
    options.Dimensions double = []
    options.Repetitions double = []
    options.SketchTrials double = []
    options.Warmups double = []
    options.ProblemSeed double = []
    options.SketchSeed double = []
    options.GmresTruncation double = []
    options.RrTruncation double = []
    options.ResidualRatioLimit double = []
    options.ResidualTolerance double = []
    options.OracleFactor double = []
    options.OutputDirectory (1,1) string = ""
    options.Force (1,1) logical = false
    options.MakePlots (1,1) logical = true
end

preset = lower(options.Preset);
switch preset
    case "smoke"
        defaults = struct(GridSizes=[8 12],Dimensions=[4 8], ...
            Repetitions=2,SketchTrials=2,Warmups=1);
    case "quick"
        defaults = struct(GridSizes=[24 48],Dimensions=[8 16 32 64], ...
            Repetitions=4,SketchTrials=3,Warmups=1);
    case "lecture"
        defaults = struct(GridSizes=[64 128], ...
            Dimensions=[16 32 64 128 256],Repetitions=6, ...
            SketchTrials=5,Warmups=1);
    otherwise
        error("experiments:config:InvalidPreset", ...
            "Preset must be smoke, quick, or lecture.");
end

suites = options.Suites;
if isempty(suites)
    suites = ["gmres" "rr"];
else
    suites = lower(reshape(suites,1,[]));
end
if any(~ismember(suites,["gmres" "rr"])) || ...
        numel(unique(suites)) ~= numel(suites)
    error("experiments:config:InvalidSuites", ...
        "Suites must be a nonempty, duplicate-free subset of gmres and rr.");
end

gridSizes = vectorOption(options.GridSizes,defaults.GridSizes,"GridSizes");
dimensions = vectorOption(options.Dimensions,defaults.Dimensions,"Dimensions");
repetitions = scalarOption(options.Repetitions,defaults.Repetitions, ...
    "Repetitions",1);
if mod(repetitions,2) ~= 0
    error("experiments:config:InvalidRepetitions", ...
        "Repetitions must be an even positive integer.");
end
sketchTrials = scalarOption(options.SketchTrials,defaults.SketchTrials, ...
    "SketchTrials",1);
warmups = scalarOption(options.Warmups,defaults.Warmups,"Warmups",0);
problemSeed = scalarOption(options.ProblemSeed,1729,"ProblemSeed",0);
sketchSeed = scalarOption(options.SketchSeed,2718,"SketchSeed",0);
gmresTruncation = scalarOption(options.GmresTruncation,4, ...
    "GmresTruncation",1);
rrTruncation = scalarOption(options.RrTruncation,2,"RrTruncation",1);
residualRatioLimit = positiveScalarOption(options.ResidualRatioLimit,10, ...
    "ResidualRatioLimit");
residualTolerance = positiveScalarOption(options.ResidualTolerance,0.1, ...
    "ResidualTolerance");
oracleFactor = positiveScalarOption(options.OracleFactor,10,"OracleFactor");

outputDirectory = options.OutputDirectory;
if strlength(outputDirectory) == 0
    outputDirectory = fullfile("results",preset);
end

cfg = struct( ...
    Preset=preset, ...
    Suites=suites, ...
    GridSizes=gridSizes, ...
    Dimensions=dimensions, ...
    Repetitions=repetitions, ...
    SketchTrials=sketchTrials, ...
    Warmups=warmups, ...
    ProblemSeed=problemSeed, ...
    SketchSeed=sketchSeed, ...
    GmresTruncation=gmresTruncation, ...
    RrTruncation=rrTruncation, ...
    ResidualRatioLimit=residualRatioLimit, ...
    ResidualTolerance=residualTolerance, ...
    OracleFactor=oracleFactor, ...
    OutputDirectory=outputDirectory, ...
    Force=options.Force, ...
    MakePlots=options.MakePlots);
end

function value = vectorOption(candidate,defaultValue,name)
if isempty(candidate)
    value = defaultValue;
else
    value = reshape(candidate,1,[]);
end
if isempty(value) || any(~isfinite(value)) || any(value < 1) || ...
        any(value ~= fix(value)) || numel(unique(value)) ~= numel(value)
    error("experiments:config:Invalid"+name, ...
        "%s must contain unique positive integers.",name);
end
end

function value = scalarOption(candidate,defaultValue,name,minimum)
if isempty(candidate)
    value = defaultValue;
else
    value = candidate;
end
if ~isscalar(value) || ~isfinite(value) || value ~= fix(value) || ...
        value < minimum
    error("experiments:config:Invalid"+name, ...
        "%s must be an integer greater than or equal to %d.",name,minimum);
end
end

function value = positiveScalarOption(candidate,defaultValue,name)
if isempty(candidate)
    value = defaultValue;
else
    value = candidate;
end
if ~isscalar(value) || ~isfinite(value) || value <= 0
    error("experiments:config:Invalid"+name, ...
        "%s must be a positive finite scalar.",name);
end
end
