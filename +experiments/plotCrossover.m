function [paths,plotData] = plotCrossover(summary,outputDirectory)
%PLOTCROSSOVER Export annotated timing and residual heatmaps.

arguments
    summary table
    outputDirectory (1,1) string
end

required = ["Suite" "GridSize" "MatrixSize" "Dimension" ...
    "ComparableCoreSpeedup" "ResidualRatioMedian" "CellStatus"];
if any(~ismember(required,string(summary.Properties.VariableNames)))
    error("experiments:plotCrossover:InvalidSummary", ...
        "The summary table does not contain the required plot columns.");
end
if ~isfolder(outputDirectory)
    error("experiments:plotCrossover:MissingOutputDirectory", ...
        "The output directory does not exist.");
end

paths = strings(0,1);
plotData = repmat(struct( ...
    Suite="",GridSizes=[],Dimensions=[],Speedup=[], ...
    ResidualRatio=[],ResidualLabels=[],Statuses=[],SpeedColorLimits=[], ...
    ResidualColorLimits=[],SpeedColorScale="",ResidualColorScale="", ...
    SpeedColorTicks=[],SpeedColorTickLabels=[],SpeedColormap=[], ...
    ResidualColormap=[]),0,1);
for suite = ["gmres" "rr"]
    suiteSummary = summary(summary.Suite == suite,:);
    if isempty(suiteSummary)
        continue
    end
    [newPaths,data] = plotSuite(suiteSummary,suite,outputDirectory);
    paths = [paths;newPaths]; %#ok<AGROW>
    plotData(end+1,1) = data; %#ok<AGROW>
end
end

function [paths,data] = plotSuite(summary,suite,outputDirectory)
grids = unique(summary.GridSize,"sorted");
dimensions = unique(summary.Dimension,"sorted");
speedup = NaN(numel(grids),numel(dimensions));
residualRatio = NaN(size(speedup));
statuses = repmat("missing",size(speedup));
matrixSizes = NaN(numel(grids),1);

for row = 1:height(summary)
    gridIndex = find(grids == summary.GridSize(row),1);
    dimensionIndex = find(dimensions == summary.Dimension(row),1);
    status = summary.CellStatus(row);
    statuses(gridIndex,dimensionIndex) = status;
    if status == "ok" || status == "unconverged"
        speedup(gridIndex,dimensionIndex) = ...
            summary.ComparableCoreSpeedup(row);
        residualRatio(gridIndex,dimensionIndex) = ...
            summary.ResidualRatioMedian(row);
    end
    matrixSizes(gridIndex) = summary.MatrixSize(row);
end

figureHandle = figure(Visible="off",Color="white", ...
    Position=[100 100 1320 500]);
cleanup = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle,1,2, ...
    Padding="compact",TileSpacing="loose");

speedAxis = nexttile(layout);
drawImage(speedAxis,speedup);
speedLimits = powerOfTwoRatioLimits(speedup);
speedAxis.ColorScale = "log";
clim(speedAxis,speedLimits);
neutralPosition = log2(1/speedLimits(1))/ ...
    log2(speedLimits(2)/speedLimits(1));
speedColors = divergingMap(257,neutralPosition);
colormap(speedAxis,speedColors);
speedBar = colorbar(speedAxis);
speedExponents = round(log2(speedLimits(1))): ...
    round(log2(speedLimits(2)));
speedTicks = 2.^speedExponents;
speedBar.Ticks = speedTicks;
speedBar.TickLabels = compose("%g",speedTicks);
speedBar.Label.String = "classical/sketched core-time ratio";
speedBar.Color = [0.1 0.1 0.1];
title(speedAxis,"Speedup ratio",Color="black");

residualAxis = nexttile(layout);
drawImage(residualAxis,residualRatio);
residualLimits = ratioLimits(residualRatio);
clim(residualAxis,residualLimits);
neutralPosition = (1-residualLimits(1))/diff(residualLimits);
residualColors = divergingMap(257,neutralPosition);
colormap(residualAxis,residualColors);
residualBar = colorbar(residualAxis);
residualBar.Label.String = "sketched/classical residual ratio";
residualBar.Color = [0.1 0.1 0.1];
title(residualAxis,"Residual-quality ratio",Color="black");

yLabels = compose("%g (n=%g)",grids,matrixSizes);
formatAxis(speedAxis,dimensions,yLabels);
formatAxis(residualAxis,dimensions,yLabels);
residualAxis.YTickLabel = [];
hold(speedAxis,"on");
hold(residualAxis,"on");
residualLabels = annotateAxes( ...
    speedAxis,residualAxis,speedup,residualRatio,statuses);
xlabel(speedAxis,"Krylov dimension k",Color="black");
xlabel(residualAxis,"Krylov dimension k",Color="black");
ylabel(speedAxis,"Grid size (matrix dimension)",Color="black");

if suite == "gmres"
    comparison = "GMRES vs. sGMRES";
else
    comparison = "RR vs. sRR";
end
sgtitle(layout,comparison+ ...
    ": speedup and residual impact", ...
    Color="black",FontWeight="bold");

pngPath = fullfile(outputDirectory,suite+"_crossover.png");
pdfPath = fullfile(outputDirectory,suite+"_crossover.pdf");
exportgraphics(figureHandle,pngPath,Resolution=180);
exportgraphics(figureHandle,pdfPath,ContentType="vector");
paths = [string(pngPath);string(pdfPath)];
data = struct( ...
    Suite=suite, ...
    GridSizes=grids, ...
    Dimensions=dimensions, ...
    Speedup=speedup, ...
    ResidualRatio=residualRatio, ...
    ResidualLabels=residualLabels, ...
    Statuses=statuses, ...
    SpeedColorLimits=speedLimits, ...
    ResidualColorLimits=residualLimits, ...
    SpeedColorScale=string(speedAxis.ColorScale), ...
    ResidualColorScale=string(residualAxis.ColorScale), ...
    SpeedColorTicks=reshape(speedBar.Ticks,1,[]), ...
    SpeedColorTickLabels=reshape(string(speedBar.TickLabels),1,[]), ...
    SpeedColormap=speedColors, ...
    ResidualColormap=residualColors);
clear cleanup
end

function drawImage(axisHandle,values)
imageHandle = imagesc(axisHandle,values);
imageHandle.AlphaData = isfinite(values);
axisHandle.Color = [0.82 0.82 0.82];
axisHandle.YDir = "normal";
box(axisHandle,"on");
end

function formatAxis(axisHandle,dimensions,yLabels)
axisHandle.XTick = 1:numel(dimensions);
axisHandle.XTickLabel = string(dimensions);
axisHandle.YTick = 1:numel(yLabels);
axisHandle.YTickLabel = yLabels;
axisHandle.TickLength = [0 0];
axisHandle.XColor = [0.1 0.1 0.1];
axisHandle.YColor = [0.1 0.1 0.1];
end

function residualLabels = annotateAxes( ...
        speedAxis,residualAxis,speedup,residualRatio,statuses)
residualLabels = strings(size(residualRatio));
for row = 1:size(statuses,1)
    for column = 1:size(statuses,2)
        status = statuses(row,column);
        if isfinite(speedup(row,column))
            text(speedAxis,column,row,sprintf("%.2gx",speedup(row,column)), ...
                HorizontalAlignment="center",FontWeight="bold", ...
                Color="black");
        end
        if isfinite(residualRatio(row,column))
            label = text(residualAxis,column,row, ...
                sprintf("%.2fx",residualRatio(row,column)), ...
                HorizontalAlignment="center",FontWeight="bold", ...
                Color="black");
            residualLabels(row,column) = string(label.String);
        end
        if status == "unconverged"
            markerColumn = column+0.32;
            markerRow = row-0.32;
            plot(speedAxis,markerColumn,markerRow,"ko", ...
                MarkerSize=8,LineWidth=1.5);
            plot(residualAxis,markerColumn,markerRow,"ko", ...
                MarkerSize=8,LineWidth=1.5);
        elseif status ~= "ok"
            text(speedAxis,column,row,"x",HorizontalAlignment="center", ...
                FontSize=15,FontWeight="bold",Color=[0.2 0.2 0.2]);
            text(residualAxis,column,row,"x",HorizontalAlignment="center", ...
                FontSize=15,FontWeight="bold",Color=[0.2 0.2 0.2]);
        end
    end
end
end

function limits = ratioLimits(values)
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    limits = [0 2];
else
    limits = [min([finiteValues;1]) max([finiteValues;1])];
    if limits(1) == limits(2)
        width = max(0.5,0.1*abs(limits(1)));
        limits = [max(0,limits(1)-width) limits(2)+width];
    end
end
end

function limits = powerOfTwoRatioLimits(values)
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    limits = [0.5 2];
    return
end
if any(finiteValues <= 0)
    error("experiments:plotCrossover:NonpositiveSpeedup", ...
        "Speedup ratios must be positive for logarithmic color scaling.");
end
exponents = [floor(log2(min([finiteValues(:);1]))) ...
    ceil(log2(max([finiteValues(:);1])))];
if exponents(1) == exponents(2)
    exponents = exponents+[-1 1];
end
limits = 2.^exponents;
end

function colors = divergingMap(count,neutralPosition)
blue = [0.15 0.35 0.70];
white = [1 1 1];
red = [0.75 0.15 0.20];
if neutralPosition <= 0
    colors = interpolateColors(white,red,count);
elseif neutralPosition >= 1
    colors = interpolateColors(blue,white,count);
else
    lowerCount = round(neutralPosition*(count-1))+1;
    upperCount = count-lowerCount+1;
    lower = interpolateColors(blue,white,lowerCount);
    upper = interpolateColors(white,red,upperCount);
    colors = [lower;upper(2:end,:)];
end
end

function colors = interpolateColors(first,last,count)
colors = [linspace(first(1),last(1),count)', ...
    linspace(first(2),last(2),count)', ...
    linspace(first(3),last(3),count)'];
end
