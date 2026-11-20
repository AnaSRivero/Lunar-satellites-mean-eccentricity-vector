function outputModelName = buildSelectiveLunarGravityModel(orekitDataDir, maxDegree, maxOrder, extraTerms)
%BUILDSELECTIVELUNARGRAVITYMODEL Write a filtered copy of the lunar gravity field.
%
%   outputModelName = buildSelectiveLunarGravityModel(orekitDataDir, ...
%                         maxDegree, maxOrder, extraTerms)
%
%   Reads 'lunar_grail.gfc' from orekitDataDir and writes, in the same
%   folder, a gravity field file that keeps
%     - all the zonal harmonics (order 0),
%     - the non-zonal harmonics of odd degree,
%     - the non-zonal harmonics of even degree listed in extraTerms
%       (both C_lm and S_lm),
%   up to (maxDegree, maxOrder), and sets every other coefficient to zero.
%   The harmonics kept by default are those retained in the analytical
%   solution of the paper.
%
%   extraTerms - N-by-2 array of [degree, order] pairs, e.g. [2 2; 4 3];
%                use zeros(0, 2) for none.
%
%   Returns the name of the filtered file (not the full path), which is
%   the name to be given to Orekit's ICGEMFormatReader. If the file
%   already exists it is reused.

    if isempty(extraTerms)
        extraTerms = zeros(0, 2);
    end
    if size(extraTerms, 2) ~= 2
        error('buildSelectiveLunarGravityModel:badTerms', ...
            'extraTerms must have two columns: [degree, order].');
    end
    extraTerms = unique(round(extraTerms), 'rows');

    % File name encoding the selection
    if isempty(extraTerms)
        tag = 'none';
    else
        pieces = arrayfun(@(k) sprintf('%d-%d', extraTerms(k, 1), extraTerms(k, 2)), ...
            1:size(extraTerms, 1), 'UniformOutput', false);
        tag = strjoin(pieces, '_');
    end
    outputModelName = sprintf('lunar_grail_reduced_d%d_o%d_extra_%s.gfc', maxDegree, maxOrder, tag);

    sourcePath = fullfile(orekitDataDir, 'lunar_grail.gfc');
    outputPath = fullfile(orekitDataDir, outputModelName);
    if isfile(outputPath)
        return;
    end

    fidIn = fopen(sourcePath, 'r');
    if fidIn == -1
        error('buildSelectiveLunarGravityModel:noSource', ...
            'Could not open the source gravity field: %s', sourcePath);
    end
    closeIn = onCleanup(@() fclose(fidIn));

    fidOut = fopen(outputPath, 'w');
    if fidOut == -1
        error('buildSelectiveLunarGravityModel:noOutput', ...
            'Could not create the filtered gravity field: %s', outputPath);
    end
    closeOut = onCleanup(@() fclose(fidOut));

    zeroCoeff = '+0.0000000000000000E+00';

    while true
        line = fgetl(fidIn);
        if ~ischar(line)
            break;
        end

        if startsWith(strtrim(line), 'max_degree')
            fprintf(fidOut, 'max_degree                  %d\n', maxDegree);
            continue;
        end

        tokens = split(string(strtrim(line)));
        if numel(tokens) < 5 || ~strcmp(tokens(1), "gfc")
            fprintf(fidOut, '%s\n', line);   % header lines are copied unchanged
            continue;
        end

        degree = str2double(tokens(2));
        order  = str2double(tokens(3));
        if isnan(degree) || isnan(order) || degree > maxDegree || order > maxOrder
            continue;
        end

        keep = order == 0 ...                                   % zonal
            || mod(degree, 2) == 1 ...                          % odd-degree non-zonal
            || any(extraTerms(:, 1) == degree & extraTerms(:, 2) == order);

        if keep
            coeffC = char(tokens(4));
            coeffS = char(tokens(5));
        else
            coeffC = zeroCoeff;
            coeffS = zeroCoeff;
        end

        fprintf(fidOut, 'gfc %4d %4d %s %s\n', degree, order, coeffC, coeffS);
    end
end
