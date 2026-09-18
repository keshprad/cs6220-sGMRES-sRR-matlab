function [raw,summary,paths] = runCrossover(options)
%RUNCROSSOVER Run the MATLAB-native GMRES and RR crossover experiment.

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

cfg = experiments.config( ...
    Preset=options.Preset, ...
    Suites=options.Suites, ...
    GridSizes=options.GridSizes, ...
    Dimensions=options.Dimensions, ...
    Repetitions=options.Repetitions, ...
    SketchTrials=options.SketchTrials, ...
    Warmups=options.Warmups, ...
    ProblemSeed=options.ProblemSeed, ...
    SketchSeed=options.SketchSeed, ...
    GmresTruncation=options.GmresTruncation, ...
    RrTruncation=options.RrTruncation, ...
    ResidualRatioLimit=options.ResidualRatioLimit, ...
    ResidualTolerance=options.ResidualTolerance, ...
    OracleFactor=options.OracleFactor, ...
    OutputDirectory=options.OutputDirectory, ...
    Force=options.Force, ...
    MakePlots=options.MakePlots);

repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
cfg.OutputDirectory = prepareOutputDirectory( ...
    cfg.OutputDirectory,cfg.Force,repositoryRoot);

raw = experiments.collectRaw(cfg);
summary = experiments.summarize(raw,cfg);
paths = experiments.writeOutputs(raw,summary,cfg,cfg.OutputDirectory);
if cfg.MakePlots
    paths.Figures = experiments.plotCrossover(raw,summary,cfg.OutputDirectory);
else
    paths.Figures = strings(0,1);
end
end

function target = prepareOutputDirectory(requested,force,repositoryRoot)
requested = char(requested);
if isempty(strtrim(requested))
    error("experiments:runCrossover:UnsafeOutputDirectory", ...
        "OutputDirectory must not be empty after configuration resolution.");
end

requestedFile = java.io.File(requested);
if ~requestedFile.isAbsolute()
    requestedFile = java.io.File(repositoryRoot,requested);
end
target = char(requestedFile.getCanonicalPath());
repositoryRoot = char(java.io.File(repositoryRoot).getCanonicalPath());
temporaryRoot = char(java.io.File(tempdir).getCanonicalPath());

insideRepository = isWithin(target,repositoryRoot);
insideTemporary = isWithin(target,temporaryRoot);
if (~insideRepository && ~insideTemporary) || ...
        strcmp(target,repositoryRoot) || strcmp(target,temporaryRoot)
    error("experiments:runCrossover:UnsafeOutputDirectory", ...
        "OutputDirectory must be a child of the repository or tempdir.");
end

if isfile(target)
    error("experiments:runCrossover:UnsafeOutputDirectory", ...
        "OutputDirectory names an existing file.");
end
if isfolder(target)
    contents = dir(target);
    names = string({contents.name});
    nonempty = any(names ~= "." & names ~= "..");
    if nonempty && ~force
        error("experiments:runCrossover:OutputExists", ...
            "OutputDirectory already exists and is nonempty. Use Force=true to replace it.");
    end
    if force
        [removed,message] = rmdir(target,"s");
        if ~removed
            error("experiments:runCrossover:OutputRemovalFailed", ...
                "Could not replace OutputDirectory: %s",message);
        end
    end
end
if ~isfolder(target)
    [created,message] = mkdir(target);
    if ~created
        error("experiments:runCrossover:OutputCreationFailed", ...
            "Could not create OutputDirectory: %s",message);
    end
end
target = string(target);
end

function result = isWithin(candidate,parent)
result = startsWith(candidate,[parent filesep]);
end
