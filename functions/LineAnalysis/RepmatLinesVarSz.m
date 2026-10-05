function [output1] = RepmatLinesVarSz(filler, target, targetType, regionLabels, user, targetLength)


for i = 1:length(target)

    if length(target{i}) > targetLength
        target{i} = target{i}(1:50);
        disp(   [' !! Truncated FA line data for FA#' num2str(i)])
    else
    end

    template = zeros([1 targetLength - length(target{i})]);    
    target{i} = [target{i} template];
end

if targetType == "array"
    tmp = zeros(size(target{i}));
elseif targetType == "vector"
    tmp = zeros(length(regionLabels), 1);
else
end

k = 1;

for i = 1:length(regionLabels)
    if regionLabels(i) == 0
        if targetType == "array"
            tmp(i, 1:end) = repmat([filler], 1, user.line.interpolationSampleSize);
        elseif targetType == "vector"
            tmp(i) = filler;
        else
        end
    else
        if targetType == "array"
            tmp(i, 1:end) = target{k};
        elseif targetType == "vector"
            tmp(i, 1:end) = target(k);
        else
        end
        k = k + 1;
    end
end

output1 = tmp;
end