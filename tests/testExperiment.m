function tests = testExperiment
%TESTEXPERIMENT Tests for the MATLAB-native crossover experiment.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repositoryRoot);
temporaryRoot = tempname;
mkdir(temporaryRoot);
testCase.TestData.RepositoryRoot = repositoryRoot;
testCase.TestData.TempFolder = temporaryRoot;
end

function teardownOnce(testCase)
rmpath(testCase.TestData.RepositoryRoot);
if isfolder(testCase.TestData.TempFolder)
    rmdir(testCase.TestData.TempFolder,"s");
end
end

function testQuickPresetAndOverrides(testCase)
cfg = experiments.config(Preset="quick",GridSizes=[5 7],MakePlots=false);

verifyEqual(testCase,cfg.GridSizes,[5 7]);
verifyEqual(testCase,cfg.Dimensions,[8 16 32 64]);
verifyEqual(testCase,cfg.Repetitions,4);
verifyEqual(testCase,cfg.OutputDirectory,fullfile("results","quick"));
verifyFalse(testCase,cfg.MakePlots);
verifyError(testCase,@() experiments.config(Repetitions=3), ...
    "experiments:config:InvalidRepetitions");
end

function testPoissonConstruction(testCase)
[A,b,details] = experiments.poisson(2,19);
expected = [4 -1 -1 0; -1 4 0 -1; -1 0 4 -1; 0 -1 -1 4];

verifyTrue(testCase,issparse(A));
verifyEqual(testCase,full(A),expected);
verifySize(testCase,A,[4 4]);
verifyEqual(testCase,b,A*details.Exact,AbsTol=1e-14);
verifyEqual(testCase,details.MatrixSize,4);
verifyEqual(testCase,details.GridSize,2);
end

function testPoissonSeedIsReproducibleAndLocal(testCase)
before = rng;
[~,~,first] = experiments.poisson(5,19);
[~,~,second] = experiments.poisson(5,19);
after = rng;

verifyEqual(testCase,after,before);
verifyEqual(testCase,first.Exact,second.Exact);
verifyEqual(testCase,norm(first.Exact),1,AbsTol=1e-14);
end

function testGmresSweepUsesPoissonProblem(testCase)
cfg = experiments.config( ...
    Suites="gmres",GridSizes=4,Dimensions=3, ...
    Repetitions=2,SketchTrials=1,Warmups=0,MakePlots=false);

raw = experiments.collectRaw(cfg);

verifyEqual(testCase,unique(raw.Problem),"poisson-2d");
verifyTrue(testCase,all(isnan(raw.MeshPeclet)));
end

function testRrSweepUsesExplicitPoissonProblem(testCase)
cfg = experiments.config( ...
    Suites="rr",GridSizes=4,Dimensions=3, ...
    Repetitions=2,SketchTrials=1,Warmups=0,MakePlots=false);

raw = experiments.collectRaw(cfg);

verifyEqual(testCase,unique(raw.Problem),"poisson-2d");
verifyEqual(testCase,unique(raw.MatrixSize),16);
verifyTrue(testCase,all(raw.OracleAvailable));
end

function testTinySweepIsPairedAndDeterministic(testCase)
cfg = experiments.config( ...
    Suites=["gmres" "rr"],GridSizes=4,Dimensions=3, ...
    Repetitions=2,SketchTrials=1,Warmups=0, ...
    OutputDirectory=fullfile(testCase.TestData.TempFolder,"first"), ...
    MakePlots=false);
before = rng;
raw = experiments.collectRaw(cfg);
after = rng;

verifyEqual(testCase,after,before);
verifyEqual(testCase,height(raw),8);
verifyEqual(testCase,sort(unique(raw.Method)), ...
    sort(["GMRES";"sGMRES";"RR";"sRR"]));
verifyTrue(testCase,all(raw.CoreTime >= 0));
verifyTrue(testCase,all(raw.WallTime >= raw.CoreTime));
verifyTrue(testCase,all(raw.Status == "ok"));

groups = findgroups(raw.Suite,raw.SketchTrial,raw.Repetition);
verifyEqual(testCase,splitapply(@numel,raw.Method,groups),2*ones(4,1));
for group = 1:max(groups)
    members = groups == group;
    verifyEqual(testCase,sort(raw.RunOrder(members)),[1;2]);
    verifyEqual(testCase,numel(unique(raw.SketchSeed(members))),1);
end
end

function testSummarySeparatesQualityGates(testCase)
raw = syntheticRows( ...
    Suite=repmat("gmres",4,1), ...
    Method=["GMRES";"sGMRES";"GMRES";"sGMRES"], ...
    Repetition=[1;1;2;2], ...
    CoreTime=[2;1;4;1], ...
    WallTime=[3;1.5;5;1.5], ...
    RelativeResidual=[0.01;0.02;0.001;0.1]);
cfg = experiments.config(GridSizes=4,Dimensions=3, ...
    Repetitions=2,SketchTrials=1,MakePlots=false);

summary = experiments.summarize(raw,cfg);

verifyEqual(testCase,summary.FixedWorkPairs,2);
verifyEqual(testCase,summary.ComparablePairs,1);
verifyEqual(testCase,summary.ConvergedPairs,1);
verifyEqual(testCase,summary.ComparableCoreSpeedup,2,AbsTol=1e-14);
verifyEqual(testCase,summary.ComparableCoreSpeedupQ25,2,AbsTol=1e-14);
verifyEqual(testCase,summary.ComparableCoreSpeedupQ75,2,AbsTol=1e-14);
verifyEqual(testCase,summary.CellStatus,"ok");
end

function testRrSummaryUsesOracleTarget(testCase)
raw = syntheticRows( ...
    Suite=repmat("rr",2,1), ...
    Method=["RR";"sRR"], ...
    Repetition=[1;1], ...
    CoreTime=[3;1], ...
    WallTime=[4;2], ...
    RitzResidual=[0.01;0.02], ...
    RitzReal=[5;5.01], ...
    RitzImag=[1;1.02], ...
    OracleDistance=[0.02;0.03], ...
    OracleResidual=[0.001;0.001]);
cfg = experiments.config(Suites="rr",GridSizes=4,Dimensions=3, ...
    Repetitions=2,SketchTrials=1,MakePlots=false);

summary = experiments.summarize(raw,cfg);

verifyEqual(testCase,summary.ComparablePairs,1);
verifyEqual(testCase,summary.ConvergedPairs,1);
verifyEqual(testCase,summary.ComparableCoreSpeedup,3,AbsTol=1e-14);
verifyEqual(testCase,summary.EigenvalueDifferenceMedian, ...
    abs((5+1i)-(5.01+1.02i)),AbsTol=1e-14);
end

function testSummaryMarksIncomparableCell(testCase)
raw = syntheticRows( ...
    Suite=repmat("gmres",2,1), ...
    Method=["GMRES";"sGMRES"], ...
    Repetition=[1;1], ...
    CoreTime=[2;1], ...
    WallTime=[3;2], ...
    RelativeResidual=[1e-4;0.1]);
cfg = experiments.config(Suites="gmres",GridSizes=4,Dimensions=3, ...
    Repetitions=2,SketchTrials=1,MakePlots=false);

summary = experiments.summarize(raw,cfg);

verifyEqual(testCase,summary.FixedWorkPairs,1);
verifyEqual(testCase,summary.ComparablePairs,0);
verifyEqual(testCase,summary.CellStatus,"not-comparable");
verifyTrue(testCase,isnan(summary.ComparableCoreSpeedup));
verifyTrue(testCase,isnan(summary.ComparableCoreSpeedupQ25));
verifyTrue(testCase,isnan(summary.ComparableCoreSpeedupQ75));
end

function testRunCrossoverWritesProtectedBundle(testCase)
output = fullfile(testCase.TestData.TempFolder,"bundle");
[raw,summary,paths] = experiments.runCrossover( ...
    Suites=["gmres" "rr"],GridSizes=4,Dimensions=3, ...
    Repetitions=2,SketchTrials=1,Warmups=0, ...
    OutputDirectory=output,MakePlots=true);

verifyEqual(testCase,numel(paths.Figures),4);
verifyTrue(testCase,all(isfile(paths.Figures)));
verifyEqual(testCase,height(raw),8);
verifyEqual(testCase,height(summary),2);
verifyTrue(testCase,isfile(paths.Raw));
verifyTrue(testCase,isfile(paths.Summary));
verifyTrue(testCase,isfile(paths.Metadata));
verifyGreaterThan(testCase,dir(paths.Raw).bytes,0);
verifyGreaterThan(testCase,dir(paths.Summary).bytes,0);
verifyGreaterThan(testCase,dir(paths.Metadata).bytes,0);
metadata = jsondecode(fileread(paths.Metadata));
verifyEqual(testCase,string(metadata.Config.Preset),"quick");

verifyError(testCase,@() experiments.runCrossover( ...
    Suites="gmres",GridSizes=4,Dimensions=3,Repetitions=2, ...
    SketchTrials=1,Warmups=0,OutputDirectory=output,MakePlots=false), ...
    "experiments:runCrossover:OutputExists");

sentinel = fullfile(testCase.TestData.TempFolder,"keep.txt");
writelines("keep",sentinel);
experiments.runCrossover( ...
    Suites="gmres",GridSizes=4,Dimensions=3,Repetitions=2, ...
    SketchTrials=1,Warmups=0,OutputDirectory=output, ...
    Force=true,MakePlots=false);
verifyTrue(testCase,isfile(sentinel));
verifyEqual(testCase,string(strtrim(fileread(sentinel))),"keep");
end

function testRunCrossoverRejectsUnsafeOutputTarget(testCase)
verifyError(testCase,@() experiments.runCrossover( ...
    Suites="gmres",GridSizes=4,Dimensions=3,Repetitions=2, ...
    SketchTrials=1,Warmups=0, ...
    OutputDirectory=testCase.TestData.RepositoryRoot,MakePlots=false), ...
    "experiments:runCrossover:UnsafeOutputDirectory");
end

function testLinePlotsExportIndependentSketchSpread(testCase)
raw = plotFixture();
rrRaw = raw;
rrRaw.Suite(:) = "rr";
rrRaw.Method(raw.Method == "GMRES") = "RR";
rrRaw.Method(raw.Method == "sGMRES") = "sRR";
rrRaw.RitzResidual = raw.RelativeResidual;
rrRaw.RitzReal(:) = 4;
rrRaw.RitzImag(:) = 0;
rrRaw.OracleDistance(:) = 0;
rrRaw.OracleResidual(:) = 1e-14;
raw = [raw;rrRaw];
cfg = experiments.config(GridSizes=4,Dimensions=3,MakePlots=false);
summary = experiments.summarize(raw,cfg);
summary.CellStatus(summary.Suite == "rr") = "unconverged";
missing = summary(1,:);
missing.Dimension = 6;
missing.CellStatus = "not-comparable";
missing.ComparableCoreSpeedup = NaN;
summary = [summary;missing];
output = fullfile(testCase.TestData.TempFolder,"line-figures");
mkdir(output);
figuresBefore = numel(findall(groot,Type="figure"));

[paths,plotData] = experiments.plotCrossover(raw,summary,output);

verifyEqual(testCase,numel(paths),4);
verifyTrue(testCase,all(isfile(paths)));
for path = reshape(paths,1,[])
    verifyGreaterThan(testCase,dir(path).bytes,0);
end
verifyEqual(testCase,numel(findall(groot,Type="figure")),figuresBefore);
data = plotData([plotData.Suite] == "gmres");
verifyEqual(testCase,data.Speedup,[2 NaN]);
verifyEqual(testCase,data.SpeedQuartileLow,[2 NaN]);
verifyEqual(testCase,data.SpeedQuartileHigh,[2 NaN]);
verifyEqual(testCase,data.ResidualRatio,[3 NaN]);
verifyEqual(testCase,data.QuartileLow,[2 NaN]);
verifyEqual(testCase,data.QuartileHigh,[4 NaN]);
verifyEqual(testCase,data.Minimum,[1 NaN]);
verifyEqual(testCase,data.Maximum,[8 NaN]);
verifyEqual(testCase,data.Samples{1},[1 2 3 4 8]);
verifyEmpty(testCase,data.Samples{2});
verifyEqual(testCase,data.XScale,"log");
verifyEqual(testCase,data.SpeedScale,"log");
verifyEqual(testCase,data.ResidualScale,"linear");
verifyEqual(testCase,data.SpeedTicks,[1 2]);
rrData = plotData([plotData.Suite] == "rr");
verifyEqual(testCase,rrData.Statuses,"unconverged");
verifyEqual(testCase,rrData.Samples{1},[1 2 3 4 8]);
end

function testResidualSamplesDoNotWeightTimingRepetitions(testCase)
raw = plotFixture();
for repetition = 3:9
    extra = raw(1:2,:);
    extra.Repetition(:) = repetition;
    raw = [raw;extra]; %#ok<AGROW>
end
values = experiments.residualSamples(raw,"gmres",3);
verifyEqual(testCase,values,[1 2 3 4 8]);
verifyEqual(testCase,median(values),3);
end

function testResidualSamplesExcludeFailedAndIncompletePairs(testCase)
raw = plotFixture();
raw.Status(raw.SketchTrial == 2) = "error";
raw.AchievedDimension(raw.SketchTrial == 3) = 2;
raw(raw.SketchTrial == 4 & raw.Method == "sGMRES",:) = [];
raw.RelativeResidual(raw.SketchTrial == 1) = 0;
verifyEqual(testCase,experiments.residualSamples(raw,"gmres",3),[1 8]);
verifyError(testCase,@() experiments.residualSamples([raw;raw(1,:)],"gmres",3), ...
    "experiments:plotCrossover:DuplicatePair");
raw.RelativeResidual(raw.SketchTrial == 5 & raw.Method == "GMRES") = 0;
verifyError(testCase,@() experiments.residualSamples(raw,"gmres",3), ...
    "experiments:plotCrossover:NonfiniteResidualRatio");
end

function testSpeedQuartilesUsePairedQualityGatedRatios(testCase)
raw = syntheticRows(Suite=repmat("gmres",10,1), ...
    Method=repmat(["GMRES";"sGMRES"],5,1), ...
    Repetition=repelem((1:5)',2), ...
    CoreTime=[1;1;4;2;9;3;16;4;500;5], ...
    WallTime=ones(10,1)*1000, ...
    RelativeResidual=[0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;1]);
cfg = experiments.config(GridSizes=4,Dimensions=3,MakePlots=false);
summary = experiments.summarize(raw,cfg);
verifyEqual(testCase,summary.ComparablePairs,4);
verifyEqual(testCase,summary.ComparableCoreSpeedup,2.5);
verifyEqual(testCase,summary.ComparableCoreSpeedupQ25,1.75);
verifyEqual(testCase,summary.ComparableCoreSpeedupQ75,3.25);
end

function raw = plotFixture()
ratios = [1 2 3 4 8];
raw = table();
for trial = 1:5
    pair = syntheticRows(Suite=repmat("gmres",4,1), ...
        Method=["GMRES";"sGMRES";"GMRES";"sGMRES"], ...
        Repetition=[1;1;2;2],CoreTime=[2;1;2;1],WallTime=[3;2;3;2], ...
        RelativeResidual=0.01*[1;ratios(trial);1;ratios(trial)]);
    pair.SketchTrial(:) = trial;
    raw = [raw;pair]; %#ok<AGROW>
end
end

function raw = syntheticRows(options)
arguments
    options.Suite (:,1) string
    options.Method (:,1) string
    options.Repetition (:,1) double
    options.CoreTime (:,1) double
    options.WallTime (:,1) double
    options.RelativeResidual (:,1) double = []
    options.RitzResidual (:,1) double = []
    options.RitzReal (:,1) double = []
    options.RitzImag (:,1) double = []
    options.OracleDistance (:,1) double = []
    options.OracleResidual (:,1) double = []
end

count = numel(options.Method);
raw = table( ...
    options.Suite,4*ones(count,1),32*ones(count,1),3*ones(count,1), ...
    options.Method,ones(count,1),options.Repetition, ...
    repmat("ok",count,1),3*ones(count,1),options.CoreTime, ...
    options.WallTime,filled(options.RelativeResidual,count), ...
    filled(options.RitzResidual,count), ...
    filled(options.RitzReal,count),filled(options.RitzImag,count), ...
    filled(options.OracleDistance,count), ...
    filled(options.OracleResidual,count),NaN(count,1),NaN(count,1), ...
    VariableNames=["Suite" "GridSize" "MatrixSize" "Dimension" ...
    "Method" "SketchTrial" "Repetition" "Status" ...
    "AchievedDimension" "CoreTime" "WallTime" "RelativeResidual" ...
    "RitzResidual" "RitzReal" "RitzImag" "OracleDistance" ...
    "OracleResidual" "ReducedCondition" "SketchSize"]);
end

function values = filled(candidate,count)
if isempty(candidate)
    values = NaN(count,1);
else
    values = candidate;
end
end
