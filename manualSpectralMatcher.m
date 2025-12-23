classdef manualSpectralMatcher < handle
    properties
        Files
        CurrentFileIdx = 1
        DataFolder

        RefLambda
        RefIntensity
        ObsSignal
        ObsMetadata

        CurrentCenterLambda = 5000;
        CurrentDispersion = 0.5;
        CurrentOrder = 1;

        Fig
        Ax
        PlotObs
        PlotRef

        lblInfo
        sldCenter
        sldDisp
        edtCenter
        edtDisp
        btnPrev
        btnNext
        btnSave
    end

    methods
        function obj = manualSpectralMatcher(dataFolder, thPath, arPath)
            obj.DataFolder = dataFolder;

            fprintf('Loading reference data...\n');
            [obj.RefLambda, obj.RefIntensity] = obj.loadAndMerge(thPath, arPath);

            obj.Files = dir(fullfile(dataFolder, '*.fit*'));
            if isempty(obj.Files)
                error('No .fit or .fits files found in %s', dataFolder);
            end

            obj.buildGUI();
            obj.loadFile(1);
        end

        function buildGUI(obj)
            obj.Fig = uifigure('Name', 'ThAr Spectral Matcher', ...
                'Position', [100 100 1200 700], ...
                'WindowKeyPressFcn', @(s,e) obj.onKeyPress(e));

            % --- Grid Layout ---
            g = uigridlayout(obj.Fig, [3, 1]);
            g.RowHeight = {50, '1x', 120};

            % --- Header ---
            pnlHeader = uipanel(g);
            pnlHeader.Layout.Row = 1;
            hb = uigridlayout(pnlHeader, [1, 4]);
            hb.ColumnWidth = {100, 100, '1x', 150};

            obj.btnPrev = uibutton(hb, 'Text', '< Prev', 'ButtonPushedFcn', @(s,e) obj.changeFile(-1));
            obj.btnNext = uibutton(hb, 'Text', 'Next >', 'ButtonPushedFcn', @(s,e) obj.changeFile(1));

            obj.lblInfo = uilabel(hb, 'Text', 'Loading...', 'FontWeight', 'bold');

            obj.btnSave = uibutton(hb, 'Text', 'Print Params to Command', 'BackgroundColor', [0.6 0.8 1], ...
                'ButtonPushedFcn', @(s,e) obj.printParams());

            % --- Plot Area ---
            obj.Ax = uiaxes(g);
            obj.Ax.Layout.Row = 2;
            obj.Ax.XGrid = 'on';
            obj.Ax.YGrid = 'on';

            title(obj.Ax, 'Top: Observed (Red) | Bottom: Reference (Blue)');

            % --- Controls ---
            pnlControls = uipanel(g, 'Title', 'Calibration Controls (Use Arrow Keys)');
            pnlControls.Layout.Row = 3;
            cgrid = uigridlayout(pnlControls, [2, 4]);
            cgrid.ColumnWidth = {'fit', '1x', 80, 20};

            uilabel(cgrid, 'Text', 'Center Lambda (A):');
            obj.sldCenter = uislider(cgrid, 'Limits', [2500 8000], ...
                'ValueChangedFcn', @(s,e) obj.updateFromSlider('center'));
            obj.edtCenter = uieditfield(cgrid, 'numeric', ...
                'ValueChangedFcn', @(s,e) obj.updateFromEdit('center'));
            uilabel(cgrid, 'Text', 'A');

            uilabel(cgrid, 'Text', 'Dispersion (A/px):');
            obj.sldDisp = uislider(cgrid, 'Limits', [0.05 1.5], ...
                'ValueChangedFcn', @(s,e) obj.updateFromSlider('disp'));
            obj.edtDisp = uieditfield(cgrid, 'numeric', ...
                'ValueChangedFcn', @(s,e) obj.updateFromEdit('disp'));
            uilabel(cgrid, 'Text', 'A/px');
        end

        function onKeyPress(obj, event)
            if ismember('shift', event.Modifier)
                dCenter = 1.0;
                dDisp = 0.001;
            else
                dCenter = 10.0;
                dDisp = 0.01;
            end

            switch event.Key
                case 'leftarrow'
                    obj.CurrentCenterLambda = obj.CurrentCenterLambda - dCenter;
                case 'rightarrow'
                    obj.CurrentCenterLambda = obj.CurrentCenterLambda + dCenter;
                case 'downarrow'
                    obj.CurrentDispersion = obj.CurrentDispersion - dDisp;
                case 'uparrow'
                    obj.CurrentDispersion = obj.CurrentDispersion + dDisp;
                otherwise
                    return; % Ignore other keys
            end

            % Enforce Limits (prevent crashing the slider)
            limC = obj.sldCenter.Limits;
            limD = obj.sldDisp.Limits;

            obj.CurrentCenterLambda = max(limC(1), min(limC(2), obj.CurrentCenterLambda));
            obj.CurrentDispersion = max(limD(1), min(limD(2), obj.CurrentDispersion));

            % Sync GUI elements
            obj.sldCenter.Value = obj.CurrentCenterLambda;
            obj.edtCenter.Value = obj.CurrentCenterLambda;
            obj.sldDisp.Value = obj.CurrentDispersion;
            obj.edtDisp.Value = obj.CurrentDispersion;

            % Redraw
            obj.updatePlot();
        end

        function loadFile(obj, idx)
            % Boundary checks
            if idx < 1 || idx > length(obj.Files), return; end
            obj.CurrentFileIdx = idx;

            fname = obj.Files(idx).name;
            [~, metadata, data] = scripts.getMetadata(obj.DataFolder, fname);


            obj.ObsMetadata = metadata;

            obj.ObsSignal = prctile(data, 95, 1);
            obj.ObsSignal = obj.ObsSignal - min(obj.ObsSignal); % Baseline removal
            obj.ObsSignal = obj.ObsSignal / max(obj.ObsSignal); % Normalize

            obj.guessParameters();

            obj.lblInfo.Text = sprintf('File: %s | Order: %d | Angle: %.2f | Position: %4.0f', ...
                fname, metadata.order, metadata.gratingAng, metadata.gratingPos);

            obj.updatePlot();
        end

        function guessParameters(obj)
            % Use the user provided logic to estimate Center Wavelength
            angle = obj.ObsMetadata.gratingAng;
            order = obj.ObsMetadata.order;

            % Linear approximations derived from user prompt:
            % 34deg -> 4000(O2) / 8000(O1)
            % 36.6deg -> 4500(O2) / 9000(O1)

            if order == 2 % Blue (3700 - 5100)
                m = (5000 - 4000) / (39.33 - 34);
                c = 4000 - m * 34;
                estLambda = m * angle + c;
                obj.CurrentOrder = 2;
                obj.sldCenter.Limits = [300 5500];
                obj.sldDisp.Value = 0.2; % Lower dispersion for Blue usually
            else % Red (5100 - 9100) or default
                m = (10000 - 8000) / (39.33 - 34);
                c = 8000 - m * 34;
                estLambda = m * angle + c;
                obj.CurrentOrder = 1;
                obj.sldCenter.Limits = [5000 8000];
                obj.sldDisp.Value = 0.4; % Higher dispersion for Red
            end

            obj.CurrentCenterLambda = estLambda;
            obj.CurrentDispersion = obj.sldDisp.Value;

            obj.sldCenter.Value = estLambda;
            obj.edtCenter.Value = estLambda;
            obj.edtDisp.Value = obj.CurrentDispersion;
        end

        function updatePlot(obj)
            nPixels = length(obj.ObsSignal);
            pixelAxis = 1:nPixels;
            centerPixel = nPixels / 2;

            lambdaAxis = obj.CurrentCenterLambda + ...
                obj.CurrentDispersion * (pixelAxis - centerPixel);

            minL = min(lambdaAxis) - 50;
            maxL = max(lambdaAxis) + 50;

            idx = obj.RefLambda > minL & obj.RefLambda < maxL;
            localRefLam = obj.RefLambda(idx);
            localRefInt = obj.RefIntensity(idx);

            if ~isempty(localRefInt)
                localRefInt = localRefInt / max(localRefInt);
                localRefInt = -1 * localRefInt;
            end

            hold(obj.Ax, 'off');
            plot(obj.Ax, lambdaAxis, obj.ObsSignal, 'Color', '#D95319', 'LineWidth', 1.0);
            hold(obj.Ax, 'on');

            if ~isempty(localRefLam)
                stem(obj.Ax, localRefLam, localRefInt, 'Color', '#0072BD', 'Marker', 'none', 'BaseValue', 0);
            end

            yline(obj.Ax, 0, 'k-'); % Zero line
            xlim(obj.Ax, [min(lambdaAxis), max(lambdaAxis)]);
            ylim(obj.Ax, [-1.1, 1.1]);
        end

        function updateFromSlider(obj, type)
            if strcmp(type, 'center')
                val = obj.sldCenter.Value;
                obj.edtCenter.Value = val;
                obj.CurrentCenterLambda = val;
            else
                val = obj.sldDisp.Value;
                obj.edtDisp.Value = val;
                obj.CurrentDispersion = val;
            end
            obj.updatePlot();
        end

        function updateFromEdit(obj, type)
            if strcmp(type, 'center')
                val = obj.edtCenter.Value;
                obj.sldCenter.Value = val;
                obj.CurrentCenterLambda = val;
            else
                val = obj.edtDisp.Value;
                obj.sldDisp.Value = val;
                obj.CurrentDispersion = val;
            end
            obj.updatePlot();
        end

        function changeFile(obj, dir)
            newIdx = obj.CurrentFileIdx + dir;
            if newIdx >= 1 && newIdx <= length(obj.Files)
                obj.loadFile(newIdx);
            end
        end

        function printParams(obj)
            fprintf('\n--- Manual Match Parameters ---\n');
            fprintf('File: %s\n', obj.Files(obj.CurrentFileIdx).name);
            fprintf('Grating Angle: %.2f\n', obj.ObsMetadata.gratingAng);
            fprintf('Center Wavelength: %.4f A\n', obj.CurrentCenterLambda);
            fprintf('Dispersion: %.5f A/pixel\n', obj.CurrentDispersion);
            fprintf('-------------------------------\n');
        end

        function [mergedLam, mergedInt] = loadAndMerge(~, thPath, arPath)
            [thInt, thLam] = cleanParse(thPath);
            [arInt, arLam] = cleanParse(arPath);
            allInt = [thInt; arInt];
            allLam = [thLam; arLam];
            [mergedLam, idx] = sort(allLam);
            mergedInt = allInt(idx);
            [mergedLam, uIdx] = unique(mergedLam);
            mergedInt = mergedInt(uIdx);
        end
    end
end


function [intensity, lambda] = cleanParse(filePath)
str = fileread(filePath);
str = regexprep(str, '\', '');
str = regexprep(str, '[a-zA-Z]+', ' ');
nums = sscanf(str, '%f');
if mod(length(nums), 2) ~= 0, nums = nums(1:end-1); end
data = reshape(nums, 2, []).';
intensity = data(:, 1);
lambda = data(:, 2);
end

% app = manualSpectralMatcher("data/thar", "data/thoriumtable2_a.txt", "data/argontable2_a.txt");