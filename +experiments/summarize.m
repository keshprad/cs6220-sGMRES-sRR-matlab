function summary = summarize(raw,cfg)
%SUMMARIZE Aggregate paired calls into crossover and quality statistics.

arguments
    raw table
    cfg (1,1) struct
end

required = ["Suite" "GridSize" "MatrixSize" "Dimension" "Method" ...
    "SketchTrial" "Repetition" "Status" "AchievedDimension" ...
    "CoreTime" "WallTime" "RelativeResidual" "RitzResidual" ...
    "RitzReal" "RitzImag" "OracleDistance" "OracleResidual" ...
    "ReducedCondition" "SketchSize"];
if any(~ismember(required,string(raw.Properties.VariableNames)))
    error("experiments:summarize:InvalidRawTable", ...
        "The raw table does not contain the required experiment columns.");
end

cells = unique(raw(:,["Suite" "GridSize" "MatrixSize" "Dimension"]), ...
    "rows","stable");
rows = repmat(emptySummaryRow(),height(cells),1);
for cellIndex = 1:height(cells)
    key = cells(cellIndex,:);
    members = raw.Suite == key.Suite & ...
        raw.GridSize == key.GridSize & ...
        raw.MatrixSize == key.MatrixSize & ...
        raw.Dimension == key.Dimension;
    rows(cellIndex) = summarizeCell(raw(members,:),cfg);
end

summary = struct2table(rows);
summary = sortrows(summary,["Suite" "GridSize" "Dimension"]);
end

function result = summarizeCell(cellRows,cfg)
suite = cellRows.Suite(1);
dimension = cellRows.Dimension(1);
if suite == "gmres"
    classicalName = "GMRES";
    sketchedName = "sGMRES";
else
    classicalName = "RR";
    sketchedName = "sRR";
end

pairKeys = unique(cellRows(:,["SketchTrial" "Repetition"]), ...
    "rows","stable");
pairCount = height(pairKeys);
fixed = false(pairCount,1);
comparable = false(pairCount,1);
converged = false(pairCount,1);
classicalCore = NaN(pairCount,1);
sketchedCore = NaN(pairCount,1);
classicalWall = NaN(pairCount,1);
sketchedWall = NaN(pairCount,1);
classicalResidual = NaN(pairCount,1);
sketchedResidual = NaN(pairCount,1);
residualRatio = NaN(pairCount,1);
eigenvalueDifference = NaN(pairCount,1);
classicalOracleDistance = NaN(pairCount,1);
sketchedOracleDistance = NaN(pairCount,1);
oracleResidual = NaN(pairCount,1);
sketchedCondition = NaN(pairCount,1);
sketchSize = NaN(pairCount,1);

for pairIndex = 1:pairCount
    pairMembers = cellRows.SketchTrial == pairKeys.SketchTrial(pairIndex) & ...
        cellRows.Repetition == pairKeys.Repetition(pairIndex);
    pair = cellRows(pairMembers,:);
    classical = pair(pair.Method == classicalName,:);
    sketched = pair(pair.Method == sketchedName,:);
    if height(classical) ~= 1 || height(sketched) ~= 1
        continue
    end

    classicalCore(pairIndex) = classical.CoreTime;
    sketchedCore(pairIndex) = sketched.CoreTime;
    classicalWall(pairIndex) = classical.WallTime;
    sketchedWall(pairIndex) = sketched.WallTime;
    sketchedCondition(pairIndex) = sketched.ReducedCondition;
    sketchSize(pairIndex) = sketched.SketchSize;

    if suite == "gmres"
        classicalResidual(pairIndex) = classical.RelativeResidual;
        sketchedResidual(pairIndex) = sketched.RelativeResidual;
        finiteResult = all(isfinite([classical.RelativeResidual, ...
            sketched.RelativeResidual]));
    else
        classicalResidual(pairIndex) = classical.RitzResidual;
        sketchedResidual(pairIndex) = sketched.RitzResidual;
        classicalValue = complex(classical.RitzReal,classical.RitzImag);
        sketchedValue = complex(sketched.RitzReal,sketched.RitzImag);
        eigenvalueDifference(pairIndex) = abs(classicalValue-sketchedValue);
        classicalOracleDistance(pairIndex) = classical.OracleDistance;
        sketchedOracleDistance(pairIndex) = sketched.OracleDistance;
        oracleResidual(pairIndex) = max( ...
            classical.OracleResidual,sketched.OracleResidual);
        finiteResult = all(isfinite([classical.RitzResidual, ...
            sketched.RitzResidual,classical.RitzReal,classical.RitzImag, ...
            sketched.RitzReal,sketched.RitzImag]));
    end

    fixed(pairIndex) = classical.Status == "ok" && ...
        sketched.Status == "ok" && ...
        classical.AchievedDimension == dimension && ...
        sketched.AchievedDimension == dimension && finiteResult && ...
        all(isfinite([classical.CoreTime,sketched.CoreTime, ...
        classical.WallTime,sketched.WallTime])) && ...
        all([classical.CoreTime,sketched.CoreTime, ...
        classical.WallTime,sketched.WallTime] >= 0);
    if ~fixed(pairIndex)
        continue
    end

    residualRatio(pairIndex) = safeRatio( ...
        sketchedResidual(pairIndex),classicalResidual(pairIndex));
    comparable(pairIndex) = residualRatio(pairIndex) <= ...
        cfg.ResidualRatioLimit;
    if suite == "rr"
        classicalTarget = oracleTargetPass(classical,cfg.OracleFactor);
        sketchedTarget = oracleTargetPass(sketched,cfg.OracleFactor);
        comparable(pairIndex) = comparable(pairIndex) && ...
            classicalTarget && sketchedTarget;
    end
    converged(pairIndex) = comparable(pairIndex) && ...
        classicalResidual(pairIndex) <= cfg.ResidualTolerance && ...
        sketchedResidual(pairIndex) <= cfg.ResidualTolerance;
end

result = emptySummaryRow();
result.Suite = suite;
result.GridSize = cellRows.GridSize(1);
result.MatrixSize = cellRows.MatrixSize(1);
result.Dimension = dimension;
result.TotalPairs = pairCount;
result.FixedWorkPairs = sum(fixed);
result.ComparablePairs = sum(comparable);
result.ConvergedPairs = sum(converged);
result.FixedWorkCoreSpeedup = ratioMedian( ...
    classicalCore,sketchedCore,fixed);
result.ComparableCoreSpeedup = ratioMedian( ...
    classicalCore,sketchedCore,comparable);
result.ConvergedCoreSpeedup = ratioMedian( ...
    classicalCore,sketchedCore,converged);
result.FixedWorkWallSpeedup = ratioMedian( ...
    classicalWall,sketchedWall,fixed);
result.ComparableWallSpeedup = ratioMedian( ...
    classicalWall,sketchedWall,comparable);
result.ConvergedWallSpeedup = ratioMedian( ...
    classicalWall,sketchedWall,converged);
result.ClassicalCoreMedian = selectedMedian(classicalCore,comparable);
result.SketchedCoreMedian = selectedMedian(sketchedCore,comparable);
result.ClassicalCoreMad = selectedMad(classicalCore,comparable);
result.SketchedCoreMad = selectedMad(sketchedCore,comparable);
result.CoreSpeedupMad = selectedMad( ...
    classicalCore./sketchedCore,comparable);
result.ClassicalWallMedian = selectedMedian(classicalWall,comparable);
result.SketchedWallMedian = selectedMedian(sketchedWall,comparable);
result.WallSpeedupMad = selectedMad( ...
    classicalWall./sketchedWall,comparable);
result.ClassicalResidualMedian = selectedMedian( ...
    classicalResidual,fixed);
result.SketchedResidualMedian = selectedMedian(sketchedResidual,fixed);
result.ResidualRatioMedian = selectedMedian(residualRatio,fixed);
result.EigenvalueDifferenceMedian = selectedMedian( ...
    eigenvalueDifference,fixed);
result.ClassicalOracleDistanceMedian = selectedMedian( ...
    classicalOracleDistance,fixed);
result.SketchedOracleDistanceMedian = selectedMedian( ...
    sketchedOracleDistance,fixed);
result.OracleResidualMedian = selectedMedian(oracleResidual,fixed);
result.SketchedConditionMedian = selectedMedian(sketchedCondition,fixed);
result.SketchSize = selectedMedian(sketchSize,fixed);

if all(cellRows.Status == "skipped-dimension")
    result.CellStatus = "skipped-dimension";
elseif ~any(fixed)
    result.CellStatus = "no-fixed-work";
elseif ~any(comparable)
    result.CellStatus = "not-comparable";
elseif ~any(converged)
    result.CellStatus = "unconverged";
else
    result.CellStatus = "ok";
end
end

function passed = oracleTargetPass(row,oracleFactor)
passed = isfinite(row.OracleDistance) && isfinite(row.OracleResidual) && ...
    row.OracleDistance <= oracleFactor*max( ...
    row.RitzResidual,row.OracleResidual);
end

function value = safeRatio(numerator,denominator)
if denominator == 0
    if numerator == 0
        value = 1;
    else
        value = Inf;
    end
else
    value = numerator/denominator;
end
end

function value = ratioMedian(classical,sketched,selection)
ratios = classical(selection)./sketched(selection);
ratios = ratios(isfinite(ratios));
value = medianOrNan(ratios);
end

function value = selectedMedian(values,selection)
selected = values(selection & isfinite(values));
value = medianOrNan(selected);
end

function value = selectedMad(values,selection)
selected = values(selection & isfinite(values));
if isempty(selected)
    value = NaN;
else
    center = median(selected);
    value = median(abs(selected-center));
end
end

function value = medianOrNan(values)
if isempty(values)
    value = NaN;
else
    value = median(values);
end
end

function row = emptySummaryRow()
row = struct( ...
    Suite="", ...
    GridSize=NaN, ...
    MatrixSize=NaN, ...
    Dimension=NaN, ...
    TotalPairs=0, ...
    FixedWorkPairs=0, ...
    ComparablePairs=0, ...
    ConvergedPairs=0, ...
    FixedWorkCoreSpeedup=NaN, ...
    ComparableCoreSpeedup=NaN, ...
    ConvergedCoreSpeedup=NaN, ...
    FixedWorkWallSpeedup=NaN, ...
    ComparableWallSpeedup=NaN, ...
    ConvergedWallSpeedup=NaN, ...
    ClassicalCoreMedian=NaN, ...
    SketchedCoreMedian=NaN, ...
    ClassicalCoreMad=NaN, ...
    SketchedCoreMad=NaN, ...
    CoreSpeedupMad=NaN, ...
    ClassicalWallMedian=NaN, ...
    SketchedWallMedian=NaN, ...
    WallSpeedupMad=NaN, ...
    ClassicalResidualMedian=NaN, ...
    SketchedResidualMedian=NaN, ...
    ResidualRatioMedian=NaN, ...
    EigenvalueDifferenceMedian=NaN, ...
    ClassicalOracleDistanceMedian=NaN, ...
    SketchedOracleDistanceMedian=NaN, ...
    OracleResidualMedian=NaN, ...
    SketchedConditionMedian=NaN, ...
    SketchSize=NaN, ...
    CellStatus="");
end
