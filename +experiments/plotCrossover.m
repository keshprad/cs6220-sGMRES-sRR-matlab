function [paths,plotData] = plotCrossover(raw,summary,outputDirectory)
%PLOTCROSSOVER Export speedup curves and residual spread across sketch trials.
% Raw measurements are required: summary medians cannot recover quartiles.

arguments
    raw table
    summary table
    outputDirectory (1,1) string
end

required = ["Suite" "GridSize" "MatrixSize" "Dimension" ...
    "ComparableCoreSpeedup" "ComparableCoreSpeedupQ25" ...
    "ComparableCoreSpeedupQ75" "CellStatus"];
if any(~ismember(required,string(summary.Properties.VariableNames)))
    error("experiments:plotCrossover:InvalidSummary", ...
        "Recompute the summary with experiments.summarize to include timing quartiles.");
end
requiredRaw = ["Suite" "GridSize" "Dimension" "Method" "SketchTrial" ...
    "Repetition" "Status" "AchievedDimension" "CoreTime" "WallTime" ...
    "RelativeResidual" "RitzResidual" "RitzReal" "RitzImag"];
if any(~ismember(requiredRaw,string(raw.Properties.VariableNames)))
    error("experiments:plotCrossover:InvalidRaw", ...
        "Raw paired measurements are required for residual spread.");
end
if ~isfolder(outputDirectory)
    error("experiments:plotCrossover:MissingOutputDirectory", ...
        "The output directory does not exist.");
end

paths = strings(0,1);
plotData = struct([]);
visible = ismember(summary.CellStatus,["ok" "unconverged"]);
upper = summary.ComparableCoreSpeedupQ75(visible);
upper = upper(isfinite(upper));
speedLimit = max(17,ceil(max([upper;0]))+1);
for suite = ["gmres" "rr"]
    suiteSummary = summary(summary.Suite == suite,:);
    if isempty(suiteSummary)
        continue
    end
    [newPaths,data] = plotSuite(raw(raw.Suite == suite,:), ...
        suiteSummary,suite,outputDirectory,speedLimit);
    paths = [paths;newPaths]; %#ok<AGROW>
    if isempty(plotData)
        plotData = data;
    else
        plotData(end+1,1) = data; %#ok<AGROW>
    end
end
end

function [paths,data] = plotSuite(raw,summary,suite,outputDirectory,speedLimit)
grids = unique(summary.GridSize,"sorted");
dimensions = unique(summary.Dimension,"sorted")';
speedup = NaN(numel(grids),numel(dimensions));
speedLow = NaN(size(speedup));
speedHigh = NaN(size(speedup));
residualRatio = NaN(size(speedup));
quartileLow = NaN(size(speedup));
quartileHigh = NaN(size(speedup));
minimum = NaN(size(speedup));
maximum = NaN(size(speedup));
samples = cell(size(speedup));
statuses = repmat("missing",size(speedup));
matrixSizes = NaN(numel(grids),1);
for row = 1:height(summary)
    g = find(grids == summary.GridSize(row),1);
    k = find(dimensions == summary.Dimension(row),1);
    if statuses(g,k) ~= "missing"
        error("experiments:plotCrossover:DuplicateSetting", ...
            "Summary contains more than one row for a grid and dimension.");
    end
    statuses(g,k) = summary.CellStatus(row);
    matrixSizes(g) = summary.MatrixSize(row);
    if ~ismember(statuses(g,k),["ok" "unconverged"])
        continue
    end
    speedup(g,k) = summary.ComparableCoreSpeedup(row);
    speedLow(g,k) = summary.ComparableCoreSpeedupQ25(row);
    speedHigh(g,k) = summary.ComparableCoreSpeedupQ75(row);
    members = raw.GridSize == grids(g) & raw.Dimension == dimensions(k);
    values = experiments.residualSamples(raw(members,:),suite,dimensions(k));
    if isempty(values)
        error("experiments:plotCrossover:MissingSamples", ...
            "A plotted setting has no valid paired residual measurements.");
    end
    samples{g,k} = values;
    residualRatio(g,k) = median(values);
    % Linear interpolation at ranks 1+(N-1)*p matches the lecture plots.
    sorted = sort(values);
    if numel(sorted) == 1
        quartiles = [sorted sorted];
    else
        quartiles = interp1(1:numel(sorted),sorted, ...
            1+(numel(sorted)-1)*[0.25 0.75]);
    end
    quartileLow(g,k) = quartiles(1);
    quartileHigh(g,k) = quartiles(2);
    minimum(g,k) = min(values);
    maximum(g,k) = max(values);
end

figureHandle = figure(Visible="off",Color="white", ...
    Position=[100 100 1320 480]);
cleanup = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle,1,2,Padding="compact",TileSpacing="loose");
speedAxis = nexttile(layout);
residualAxis = nexttile(layout);
hold(speedAxis,"on");
hold(residualAxis,"on");
colors = [0 0.4470 0.6980;0.8350 0.3690 0];
if numel(grids) > 2
    colors = lines(numel(grids));
end
markers = ["o" "s" "^" "d" "v"];
handles = gobjects(numel(grids),1);
for g = 1:numel(grids)
    color = colors(g,:);
    faintColor = 0.55*color+0.45;
    marker = markers(mod(g-1,numel(markers))+1);
    label = sprintf("g = %g (n = %g)",grids(g),matrixSizes(g));
    validSpeed = isfinite(speedLow(g,:)) & isfinite(speedHigh(g,:));
    starts = find(diff([false validSpeed false]) == 1);
    stops = find(diff([false validSpeed false]) == -1)-1;
    for segment = 1:numel(starts)
        idx = starts(segment):stops(segment);
        fill(speedAxis,[dimensions(idx) fliplr(dimensions(idx))], ...
            [speedLow(g,idx) fliplr(speedHigh(g,idx))],color, ...
            FaceAlpha=0.13,EdgeColor="none",HandleVisibility="off");
    end
    handles(g) = errorbar(speedAxis,dimensions,speedup(g,:), ...
        speedup(g,:)-speedLow(g,:),speedHigh(g,:)-speedup(g,:), ...
        Color=color,Marker=marker,MarkerFaceColor=color, ...
        LineWidth=1.5,MarkerSize=5,CapSize=7,DisplayName=label);
    % Draw each contiguous band separately so missing settings remain gaps.
    valid = isfinite(quartileLow(g,:));
    starts = find(diff([false valid false]) == 1);
    stops = find(diff([false valid false]) == -1)-1;
    for segment = 1:numel(starts)
        idx = starts(segment):stops(segment);
        fill(residualAxis,[dimensions(idx) fliplr(dimensions(idx))], ...
            [quartileLow(g,idx) fliplr(quartileHigh(g,idx))],color, ...
            FaceAlpha=0.13,EdgeColor="none",HandleVisibility="off");
    end
    offset = (g-(numel(grids)+1)/2)*0.045;
    errorbar(residualAxis,dimensions*2^offset,residualRatio(g,:), ...
        residualRatio(g,:)-minimum(g,:),maximum(g,:)-residualRatio(g,:), ...
        LineStyle="none",Color=faintColor,CapSize=6, ...
        LineWidth=0.9,HandleVisibility="off");
    for k = 1:numel(dimensions)
        values = samples{g,k};
        if isempty(values)
            continue
        end
        if numel(values) == 1
            jitter = 0;
        else
            jitter = linspace(-0.06,0.06,numel(values));
        end
        scatter(residualAxis,dimensions(k)*2.^(offset+jitter),values, ...
            22,color,marker,"filled",MarkerFaceAlpha=0.3, ...
            MarkerEdgeAlpha=0.3,HandleVisibility="off");
    end
    plot(residualAxis,dimensions,residualRatio(g,:), ...
        Color=color,Marker=marker,MarkerFaceColor=color, ...
        LineWidth=1.8,MarkerSize=5,HandleVisibility="off");

end

for ax = [speedAxis residualAxis]
    ax.XScale = "log";
    ax.XTick = dimensions;
    ax.XTickLabel = string(dimensions);
    ax.XLim = [min(dimensions)/1.12 max(dimensions)*1.12];
    ax.YGrid = "on";
    ax.XMinorTick = "off";
    ax.FontSize = 11;
    ax.Color = "white";
    ax.XColor = [0.1 0.1 0.1];
    ax.YColor = [0.1 0.1 0.1];
    box(ax,"off");
    yline(ax,1,"--",Color=[0.35 0.35 0.35],HandleVisibility="off");
    xlabel(ax,"Krylov dimension k",Color="black");
end
speedAxis.YScale = "linear";
if any(isfinite(speedup(:)) & speedup(:) <= 0)
    error("experiments:plotCrossover:NonpositiveSpeedup", ...
        "Speedup ratios must be positive.");
end
speedAxis.YTick = 0:max(2,2*ceil(speedLimit/18)):speedLimit;
speedAxis.YTickLabel = compose("%g",speedAxis.YTick);
speedAxis.YLim = [0 speedLimit];
if isfinite(speedup(end,end))
    text(speedAxis,dimensions(end)/1.07,speedup(end,end)+0.04*speedLimit, ...
        sprintf("%.1fx",speedup(end,end)),HorizontalAlignment="right", ...
        Color=colors(end,:),FontSize=12,FontWeight="bold");
end
residualValues = [minimum(:);maximum(:);1];
residualValues = residualValues(isfinite(residualValues));
padding = max(0.05,0.1*(max(residualValues)-min(residualValues)));
residualAxis.YLim = [max(0,min(residualValues)-padding) max(residualValues)+padding];
title(speedAxis,"Median speedup",Color="black");
ylabel(speedAxis,"Classical / sketched core time",Color="black");
title(residualAxis,"Residual impact across sketches",Color="black");
ylabel(residualAxis,"Sketched / classical residual",Color="black");
key = legend(speedAxis,handles,Location="northoutside", ...
    Orientation="horizontal",Box="off",TextColor="black");
key.Layout.Tile = "north";
if suite == "gmres"
    comparison = "GMRES versus sGMRES";
else
    comparison = "RR versus sRR";
end
title(layout,comparison,Color="black",FontWeight="bold");

pngPath = fullfile(outputDirectory,suite+"_crossover.png");
pdfPath = fullfile(outputDirectory,suite+"_crossover.pdf");
exportgraphics(figureHandle,pngPath,Resolution=180);
exportgraphics(figureHandle,pdfPath,ContentType="vector");
paths = [string(pngPath);string(pdfPath)];
data = struct(Suite=suite,GridSizes=grids,Dimensions=dimensions, ...
    Speedup=speedup,ResidualRatio=residualRatio, ...
    SpeedQuartileLow=speedLow,SpeedQuartileHigh=speedHigh, ...
    QuartileLow=quartileLow,QuartileHigh=quartileHigh, ...
    Minimum=minimum,Maximum=maximum,Statuses=statuses, ...
    XScale=string(speedAxis.XScale),SpeedScale=string(speedAxis.YScale), ...
    ResidualScale=string(residualAxis.YScale),SpeedTicks=speedAxis.YTick, ...
    SpeedLimits=speedAxis.YLim);
data.Samples = samples;
clear cleanup
end
