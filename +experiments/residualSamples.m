function values = residualSamples(raw,suite,dimension)
%RESIDUALSAMPLES One residual ratio per independent sketch, not timing repeat.
% Use the same successful fixed-work pairs as the residual summary. Repeated
% measurements within a sketch are collapsed to their median paired ratio.

if suite == "gmres"
    classicalName = "GMRES";
    sketchedName = "sGMRES";
    field = "RelativeResidual";
else
    classicalName = "RR";
    sketchedName = "sRR";
    field = "RitzResidual";
end
keys = unique(raw(:,["SketchTrial" "Repetition"]),"rows");
ratios = NaN(height(keys),1);
for j = 1:height(keys)
    pair = raw(raw.SketchTrial == keys.SketchTrial(j) & ...
        raw.Repetition == keys.Repetition(j),:);
    classical = pair(pair.Method == classicalName,:);
    sketched = pair(pair.Method == sketchedName,:);
    if height(classical) > 1 || height(sketched) > 1
        error("experiments:plotCrossover:DuplicatePair", ...
            "A timing pair contains duplicate method observations.");
    end
    if height(classical) ~= 1 || height(sketched) ~= 1
        continue
    end
    rows = [classical;sketched];
    times = [rows.CoreTime;rows.WallTime];
    residuals = rows.(field);
    valid = all(rows.Status == "ok") && ...
        all(rows.AchievedDimension == dimension) && ...
        all(isfinite(times) & times >= 0) && ...
        all(isfinite(residuals) & residuals >= 0);
    if suite == "rr"
        valid = valid && all(isfinite([rows.RitzReal;rows.RitzImag]));
    end
    if ~valid
        continue
    end
    if residuals(1) == 0 && residuals(2) == 0
        ratios(j) = 1;
    else
        ratios(j) = residuals(2)/residuals(1);
    end
end
trials = unique(keys.SketchTrial);
values = NaN(1,numel(trials));
for j = 1:numel(trials)
    selected = ratios(keys.SketchTrial == trials(j));
    selected = selected(~isnan(selected));
    if ~isempty(selected)
        values(j) = median(selected);
    end
end
values = values(~isnan(values));
if any(~isfinite(values))
    error("experiments:plotCrossover:NonfiniteResidualRatio", ...
        "An infinite residual ratio cannot be represented by finite whiskers.");
end
end
