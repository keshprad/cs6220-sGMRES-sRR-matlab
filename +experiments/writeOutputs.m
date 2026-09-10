function paths = writeOutputs(raw,summary,cfg,outputDirectory)
%WRITEOUTPUTS Write experiment tables and reproducibility metadata.

arguments
    raw table
    summary table
    cfg (1,1) struct
    outputDirectory (1,1) string
end

if ~isfolder(outputDirectory)
    error("experiments:writeOutputs:MissingOutputDirectory", ...
        "The validated output directory does not exist.");
end

rawPath = fullfile(outputDirectory,"raw.csv");
summaryPath = fullfile(outputDirectory,"summary.csv");
metadataPath = fullfile(outputDirectory,"metadata.json");
writetable(raw,rawPath);
writetable(summary,summaryPath);

repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
metadata = environmentMetadata(cfg,repositoryRoot);
writeText(metadataPath,jsonencode(metadata,PrettyPrint=true));

paths = struct( ...
    Raw=string(rawPath), ...
    Summary=string(summaryPath), ...
    Metadata=string(metadataPath));
end

function metadata = environmentMetadata(cfg,repositoryRoot)
metadata = struct( ...
    Config=cfg, ...
    MatlabRelease=string(version("-release")), ...
    MatlabVersion=string(version), ...
    Architecture=string(computer("arch")), ...
    Toolboxes=ver, ...
    Blas="", ...
    TimestampUtc=string(datetime("now",TimeZone="UTC", ...
        Format="yyyy-MM-dd'T'HH:mm:ss.SSSXXX")), ...
    GitCommit="", ...
    WorkingTreeClean=false);

try
    metadata.Blas = string(version("-blas"));
catch
    metadata.Blas = "unavailable";
end

commitCommand = sprintf('git -C "%s" rev-parse HEAD',repositoryRoot);
[commitStatus,commit] = system(commitCommand);
if commitStatus == 0
    metadata.GitCommit = strtrim(string(commit));
end
[treeStatus,tree] = system( ...
    sprintf('git -C "%s" status --porcelain',repositoryRoot));
metadata.WorkingTreeClean = treeStatus == 0 && strlength(strtrim(tree)) == 0;
end

function writeText(path,text)
[file,openMessage] = fopen(path,"w","n","UTF-8");
if file < 0
    error("experiments:writeOutputs:OpenFailed", ...
        "Could not open metadata file: %s",openMessage);
end
cleanup = onCleanup(@() fclose(file));
count = fprintf(file,"%s",text);
if count ~= strlength(text)
    error("experiments:writeOutputs:WriteFailed", ...
        "Metadata output was incomplete.");
end
clear cleanup
end
