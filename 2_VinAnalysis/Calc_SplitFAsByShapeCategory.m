function T = Calc_SplitFAsByShapeCategory(T, categories, options)
% Classifies FAs into a variable number of categories, each defined by
% AND-combined filter rules. FAs are assigned to the first matching
% category. Unmatched FAs go to a 'Rest' category with a warning.
%
% --- Inputs ---
%   T          : table with FA data
%   categories : cell array of categories. Each category is a cell array
%                of rule triplets: {varName, sign, threshold, ...}
%                e.g. {'MajAx', '>', 30, 'MinMajRatio', '<', 0.5}
%
% --- Options ---
%   FRET            : logical, default false
%   order           : integer array to rearrange output categories.
%                     e.g. order=[3 1 2 4] means output Cat1 = processing
%                     category 3, Cat2 = processing category 1, etc.
%                     MaskClass pixel values follow output order.
%                     Default: no reordering.
%   includePatterns : cell of strings, default {}
%   excludePatterns : cell of strings, default {}
%
% --- Example ---
%   categories = {
%       {'MajAx', '>', 30, 'MinMajRatio', '<', 0.5},  % processed 1st
%       {'MajAx', '>', 15, 'MinMajRatio', '<', 0.5},  % processed 2nd
%       {'MinMajRatio', '>', 0.5, 'MajAx', '<', 15},  % processed 3rd
%       {'MinMajRatio', '>', 0.5, 'MajAx', '>', 15},  % processed 4th
%   };
%   T = Calc_SplitFAsByShapeCategory(T, categories, ...
%           FRET=true, order=[3 4 2 1]);
%   % Output Cat1 = processing cat 3 → pixel value 1 in MaskClass
%   % Output Cat2 = processing cat 4 → pixel value 2 in MaskClass
%   % Output Cat3 = processing cat 2 → pixel value 3 in MaskClass
%   % Output Cat4 = processing cat 1 → pixel value 4 in MaskClass

    arguments
        T                       table
        categories              cell
        options.FRET            logical = false
        options.order           double  = []
        options.includePatterns cell    = {}
        options.excludePatterns cell    = {}
    end

    nCats = numel(categories);

    %% --- Validate: each category must have triplets only ---
    for i = 1:nCats
        rules = categories{i};
        if mod(numel(rules), 3) ~= 0
            error(['Category %d has %d elements — must be a multiple of 3 ' ...
                   '(varName, sign, threshold).'], i, numel(rules));
        end
        for ri = 1:3:numel(rules)
            if ~ischar(rules{ri}) && ~isstring(rules{ri})
                error('Category %d, rule %d: variable name must be a string.', i, (ri+2)/3);
            end
            if ~ismember(rules{ri+1}, {'<','>','<=','>='})
                error('Category %d, rule %d: invalid operator "%s".', i, (ri+2)/3, rules{ri+1});
            end
            if ~isnumeric(rules{ri+2})
                error('Category %d, rule %d: threshold must be numeric.', i, (ri+2)/3);
            end
        end
    end

    %% --- Validate: order argument ---
    if isempty(options.order)
        procToOut = 1:nCats;   % processing index → output index
    else
        if numel(options.order) ~= nCats
            error('order has %d entries but there are %d categories.', ...
                  numel(options.order), nCats);
        end
        if ~isequal(sort(options.order(:)'), 1:nCats)
            error('order must be a permutation of 1:%d.', nCats);
        end
        % options.order(i) = which processing category becomes output Cat i
        % We need procToOut(p) = output index for processing category p
        procToOut = zeros(1, nCats);
        for i = 1:nCats
            procToOut(options.order(i)) = i;
        end
    end

    %% --- Validate: warn if a category's rules are subset of an earlier one ---
    % For each pair (i,j) where i < j, check if cat j's rules are a subset
    % of cat i's rules (meaning cat j can never be reached)
    for i = 1:nCats-1
        for j = i+1:nCats
            if isSubsetRules(categories{i}, categories{j})
                warning(['Category %d rules are a strict subset of Category %d rules. ' ...
                         'Category %d may never be reached — check your category order.'], ...
                         j, i, j);
            end
        end
    end

    %% --- Processing category names (internal) ---
    % Always Cat1..CatN in processing order; Rest for unmatched
    procCatNames = "Cat" + string(1:nCats);
    allProcNames = [procCatNames, "Rest"];

    %% --- Output category names (after reordering) ---
    % procToOut maps processing index → output index
    % output Cat1..CatN + Rest
    outCatNames = "Cat" + string(1:nCats);

    %% --- Filter rows ---
    T = Filter_IncludeExclude(T, options.includePatterns, options.excludePatterns);

    %% --- Helper: apply comparison operator ---
    function out = applySign(val, sign, thr)
        switch sign
            case '<',   out = val <  thr;
            case '>',   out = val >  thr;
            case '<=',  out = val <= thr;
            case '>=',  out = val >= thr;
        end
        out(~isfinite(val)) = false;
    end

    %% --- Helper: check one FA against all rules of one category ---
    function match = matchesCategory(rules, metricValues, k)
        match = true;
        for ri = 1:3:numel(rules)
            vec = metricValues(char(rules{ri}));
            if ~applySign(vec(k), rules{ri+1}, rules{ri+2})
                match = false;
                return;
            end
        end
    end

    nRows = height(T);
    totalUnmatched = 0;

    % initialize summary columns
    MIN_FA = 5;
    summaryInt  = nan(nRows, nCats);
    summaryFRET = nan(nRows, nCats);
    summaryFrac = nan(nRows, nCats);
    summaryIntPx  = nan(nRows, nCats);
    summaryFRETPx = nan(nRows, nCats);

    fprintf('Calc_SplitFAsByShapeCategory: processing %d rows...\n', nRows);

    for r = 1:nRows

        if mod(r, 50) == 0
            fprintf('  row %d / %d\n', r, nRows);
        end

        %% --- Get mask & images ---
        if options.FRET
            F    = T.FRETimageFA{r};
            F(isnan(F)) = 0;
            if isempty(F), continue; end
        else
        end

        Mask = T.Mask{r};
        Int  = T.Int{r};
        Int(~Mask) = 0;

        %% --- Connected components ---
        CC   = bwconncomp(Mask, 8);
        nFAs = CC.NumObjects;

        %% --- Collect unique variable names across all categories ---
        allVarNames = string.empty;
        for i = 1:nCats
            rules = categories{i};
            for ri = 1:3:numel(rules)
                allVarNames(end+1) = string(rules{ri}); %#ok<AGROW>
            end
        end
        allVarNames = unique(allVarNames);

        %% --- Load metric vectors for this row ---
        metricValues = containers.Map('KeyType','char','ValueType','any');
        for i = 1:numel(allVarNames)
            vName = char(allVarNames(i));
            vec   = T.(vName){r};
            if numel(vec) ~= nFAs
                error(['Row %d: FA count mismatch for variable "%s" ' ...
                       '(%d FAs in metric vs %d in connected components).'], ...
                       r, vName, numel(vec), nFAs);
            end
            metricValues(vName) = vec;
        end

        %% --- Initialise output fields (in output order) ---
        for i = 1:numel(allProcNames)   % includes Rest
            c = allProcNames(i);
            T.("Int_" + c){r} = [];
            if options.FRET
                T.("FRET_" + c){r} = [];
            end
        end

        %% --- Class mask (in output pixel values) ---
        if options.FRET
            ClassMask = zeros(size(F));
        else
            ClassMask = zeros(size(Mask));
        end

        %% --- Classify each FA ---
        rowUnmatched = 0;

        % per object
        faInt   = nan(nFAs, 1);
        faFRET  = nan(nFAs, 1);
        faClass = zeros(nFAs, 1);

        % for all pixels
        catPixInt  = cell(nCats, 1);   % accumulate raw pixel intensities per proc category
        catPixFRET = cell(nCats, 1);
        for i = 1:nCats
            catPixInt{i}  = [];
            catPixFRET{i} = [];
        end

        for k = 1:nFAs

            pix      = CC.PixelIdxList{k};
            meanInt  = mean(Int(pix), 'omitnan');
            faInt(k) = meanInt;
            if options.FRET
                meanF = mean(F(pix), 'omitnan');
                faFRET(k) = meanF;
            end
            
            % Find first matching processing category
            procIdx = 0;
            for i = 1:nCats
                if matchesCategory(categories{i}, metricValues, k)
                    procIdx = i;
                    break;
                end
            end

            faClass(k) = procIdx;   % 0 = Rest, 1..nCats = processing category index;

            if procIdx == 0
                % No match → Rest
                assignedName = "Rest";
                pixelVal     = nCats + 1;
                rowUnmatched = rowUnmatched + 1;
            else
                % Map processing index → output index
                outIdx       = procToOut(procIdx);
                assignedName = "Cat" + string(procIdx);  % stored under proc name
                pixelVal     = outIdx;                    % pixel value = output index
                % for all pixels
                catPixInt{procIdx}  = [catPixInt{procIdx};  Int(pix)];
                if options.FRET
                    catPixFRET{procIdx} = [catPixFRET{procIdx}; F(pix)];
                end
            end

            ClassMask(pix) = pixelVal;
            T.("Int_" + assignedName){r}(end+1,1) = meanInt;
            if options.FRET
                T.("FRET_" + assignedName){r}(end+1,1) = meanF;
            end

        end

        if rowUnmatched > 0
            totalUnmatched = totalUnmatched + rowUnmatched;
        end

        % make summaryInt, summaryFRET and summaryFrac
        validMask  = faClass >= 1 & faClass <= nCats;
        nValid     = sum(validMask);
        faClassBin = faClass(validMask);
        
        for c = 1:nCats
            sel = faClass == c;
            if sum(sel) < MIN_FA, continue; end
        
            summaryInt(r, c)   = mean(faInt(sel), 'omitnan');
            summaryIntPx(r, c) = mean(catPixInt{c}, 'omitnan');
        
            if options.FRET
                summaryFRET(r, c)   = mean(faFRET(sel), 'omitnan');
                summaryFRETPx(r, c) = mean(catPixFRET{c}, 'omitnan');
            end
        
            if nValid >= MIN_FA
                summaryFrac(r, c) = sum(sel & validMask) / nValid;
            end
        end
        % make image with index per category
        T.MaskClass{r} = ClassMask;
        %T.FA_summaryInt  = summaryInt;
        T.FA_summaryIntPx  = summaryIntPx;
        T.FA_summaryFrac = summaryFrac;
        if options.FRET
            %T.FA_summaryFRET = summaryFRET;
            T.FA_summaryFRETPx = summaryFRETPx;
        end
    end

    %% --- Rename columns from processing order to output order ---
    % e.g. if procToOut(3)=1, then Int_Cat3 → Int_Cat1
    if ~isempty(options.order)
        T = renameCatColumns(T, procToOut, nCats, options.FRET);
    end

    %% --- Final warning for unmatched FAs ---
    if totalUnmatched > 0
        warning(['Calc_SplitFAsByShapeCategory: %d FA(s) across all rows did not ' ...
                 'match any category and were placed in "Rest". Check your category ' ...
                 'definitions.'], totalUnmatched);
    end

    fprintf('Done.\n');
end


%% =========================================================================
function T = renameCatColumns(T, procToOut, nCats, doFRET)
% Rename Int_Cat<proc> → Int_Cat<out> (and FRET_Cat<proc> → FRET_Cat<out>)
% using a temporary name to avoid collisions during renaming.

    prefixes = "Int_";
    if doFRET
        prefixes = [prefixes, "FRET_"];
    end

    for p = 1:nCats
        o = procToOut(p);
        for pr = prefixes
            oldName = pr + "Cat" + string(p);
            tmpName = pr + "CatTMP" + string(p);
            if ismember(oldName, T.Properties.VariableNames)
                T = renamevars(T, oldName, tmpName);
            end
        end
    end

    for p = 1:nCats
        o = procToOut(p);
        for pr = prefixes
            tmpName = pr + "CatTMP" + string(p);
            newName = pr + "Cat"    + string(o);
            if ismember(tmpName, T.Properties.VariableNames)
                T = renamevars(T, tmpName, newName);
            end
        end
    end

    % Sort columns so Cat1..CatN appear in order (keep other columns in place)
    varNames = T.Properties.VariableNames;
    catCols  = strings(1, nCats);
    for pr = prefixes
        for i = 1:nCats
            catCols(i) = pr + "Cat" + string(i);
        end
        % Find position of first cat column and reorder from there
        idx = find(ismember(varNames, catCols), 1);
        if ~isempty(idx)
            otherCols = varNames(~ismember(varNames, catCols));
            before    = otherCols(1:idx-1);
            after     = otherCols(idx:end);
            T = T(:, [before, catCols, after]);
        end
    end
end


%% =========================================================================
function result = isSubsetRules(rulesA, rulesB)
% Returns true if every rule in rulesB is also present in rulesA
% (same variable, same operator, same threshold), meaning any FA matching
% B also matches A — so B can never be reached after A.

    result = false;
    nB = numel(rulesB) / 3;
    nA = numel(rulesA) / 3;
    if nB == 0 || nA == 0, return; end

    matchCount = 0;
    for bi = 1:nB
        bVar  = string(rulesB{(bi-1)*3 + 1});
        bSign = rulesB{(bi-1)*3 + 2};
        bThr  = rulesB{(bi-1)*3 + 3};
        for ai = 1:nA
            aVar  = string(rulesA{(ai-1)*3 + 1});
            aSign = rulesA{(ai-1)*3 + 2};
            aThr  = rulesA{(ai-1)*3 + 3};
            if aVar == bVar && strcmp(aSign, bSign) && aThr == bThr
                matchCount = matchCount + 1;
                break;
            end
        end
    end

    % B is a subset of A if all of B's rules appear in A
    result = (matchCount == nB);
end

% function T = Calc_SplitFAsByShapeCategory(T, ...
%         MidVar, MidSign, MidThresh, ...         % Line vs Round
%         BelowVar, BelowSign, BelowThresh, ...   % ThinLine vs EllipsoidLine
%         AboveVar, AboveSign, AboveThresh, FRET, ...   % SmallRound vs BigRound
%         includePatterns, excludePatterns)
% % --- Inputs ---
% % - T.FRETimageFA
% % - T.Area
% % - T.MinMajRatio
% %
% % --- Outputs ---
% % - T.Area_Cat1/2/3/4: List of areas within the category
% % - T.FRET_Cat1/2/3/4: List of FRETav within the categoryù
% % - T.MaskClass: 512x512 mask where Cat1/2/3/4 has value 1/2/3/4 (bg = 0)
% % 
% % This function:
% %   - Uses connected components of FRETimageFA (non‑NaN = FA pixels)
% %   - Uses precomputed FA-level metrics stored in table columns:
% %         T.(MidVar){r}(k)
% %         T.(BelowVar){r}(k)
% %         T.(AboveVar){r}(k)
% %   - Classifies each FA into:
% %         Cat1 (1), Cat2 (2), Cat3 (3), Cat4 (4)
% %   - Builds one 512×512 mask per row with these values
% %   - Computes area & FRET per FA in each category
% %   - Stores PixelIdxList per FA for each category
% 
% % Filter if necessary 
% if nargin < 12 || isempty(includePatterns)
%     includePatterns = {};
% end
% if nargin < 13 || isempty(excludePatterns)
%     excludePatterns = {};
% end
% 
% T = Filter_IncludeExclude(T, includePatterns, excludePatterns);
% 
%     %% Helper: comparison operator
%     function out = applySign(val, sign, thr)
%         switch sign
%             case '<',   out = val <  thr;
%             case '>',   out = val >  thr;
%             case '<=',  out = val <= thr;
%             case '>=',  out = val >= thr;
%             otherwise, error('Invalid comparison operator: %s', sign);
%         end
%         out(~isfinite(val)) = false;
%     end
% 
%     nRows = height(T);
% 
%     fprintf('Calc_SplitFAsByShapeCategory: processing %d rows...\n', nRows);
% 
%     for r = 1:nRows
% 
%         if mod(r, 50) == 0
%             fprintf('  row %d / %d\n', r, nRows);
%         end
% 
%         % Get Mask
%         if FRET == true
%             % Get FRET & mask
%             Mask = ~isnan(T.FRETimageFA{r});
%             F = T.FRETimageFA{r};
%             F(isnan(F)) = 0; % 512×512 double
%             if isempty(F), continue; end
%         else
%             Mask = T.Mask{r};
%         end
% 
%         % Get Int
%         Int = T.Int{r};
%         Int(~Mask) = 0; % to apply mask to intensity image
% 
%         %% 1) Build FA mask & components
%         %BW = ~isnan(F);
%         BW = Mask;
%         CC = bwconncomp(BW, 4);
%         nFAs = CC.NumObjects;
% 
%         %% 2) Retrieve precomputed metrics
%         midVar   = T.(MidVar){r};     % Nx1 double
%         belowVar = T.(BelowVar){r};   % Nx1 double
%         aboveVar = T.(AboveVar){r};   % Nx1 double
% 
%         % Safety check
%         if length(midVar) ~= nFAs
%             error('FA count mismatch between CC and table variables.');
%         end
% 
%         %% 3) Prepare category outputs
%         cats = ["Cat1","Cat2","Cat3","Cat4"];
%         for c = cats
%             if FRET
%                 T.("FRET_" + c){r} = [];
%             else
%             end
%             T.("Int_" + c){r} = [];
%             %T.("Idx_"  + c){r} = {};
%         end
% 
%         %% 4) Prepare unified class mask
%         if FRET
%             ClassMask = zeros(size(F));
%         else
%             ClassMask = zeros(size(Mask));
%         end
% 
%         %% 5) Loop over FAs
%         for k = 1:nFAs
% 
%             pix = CC.PixelIdxList{k};
% 
%             % % FRET per FA (ignore NaNs)
%             % vals = F(pix);
%             % vals = vals(~isnan(vals));
%             % meanF = mean(vals, 'omitnan');
% 
%             if FRET
%                 % F image
%                 meanF = mean(F(pix), 'omitnan');
%             else
%             end
%             % Int image
%             meanInt = mean(Int(pix), 'omitnan');
% 
%             %% Classification using provided metrics
% 
%             % 1) Cat1-2 vs Cat3-4
%             isMiddle = applySign(midVar(k), MidSign, MidThresh);
%             % 2) Cat1 vs Cat2
%             isFirstQuart = applySign(belowVar(k), BelowSign, BelowThresh);
%             % 3) Cat3 vs Cat4
%             isThirdQuart = applySign(aboveVar(k), AboveSign, AboveThresh);
% 
%             % After
%             if isMiddle && isFirstQuart
%                 classCode = 1; c = "Cat1";
%             elseif isMiddle && ~isFirstQuart
%                 classCode = 2; c = "Cat2";
%             elseif ~isMiddle && isThirdQuart
%                 classCode = 3; c = "Cat3";
%             else
%                 classCode = 4; c = "Cat4";
%             end
% 
%             %% Fill class mask
%             ClassMask(pix) = classCode;
% 
%             %% Append outputs
%             T.("Int_" + c){r}(end+1,1) = meanInt;
%             if FRET
%                 T.("FRET_" + c){r}(end+1,1) = meanF;
%             else
%             end
%             %T.("Idx_"  + c){r}{end+1}   = pix;
% 
%         end
% 
%         %% 6) Store unified classification mask
%         T.MaskClass{r} = ClassMask;
% 
%     end
%     fprintf('Done.\n');
% end