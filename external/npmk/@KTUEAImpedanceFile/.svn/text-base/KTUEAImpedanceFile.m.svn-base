classdef KTUEAImpedanceFile

% Plots a graphical representation of the impedance data collected by
% Cerebus and saved in a text file.
%
% Usage:
%
%   plotImpedances(mapfile)
%
%   WHERE:
%
%   mapfile:    Is the mapfile structure obtained with mapfileLoader
%               function that maps the electrodes to channel numbers. If
%               mapfile is not provided then a default map will be used.
%
%   To see a list of methods type methods(plotImpedances).
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   Version 1.1.0
%
%   Ideas borrowed from Sergey Stavisky, The BrainGate Project

    properties (Hidden, SetAccess = private, GetAccess = private)
        chanCount    = 96;
        impedanceDataValues = NaN(1, 96);
        Mapfile
        fileHandle
        pathHandle
    end
    methods(Hidden)
        function obj = KTUEAImpedanceFile(Mapfile)
            if ~exist('Mapfile', 'var')
                obj.Mapfile = KTUEAMapFile;
            else
                obj.Mapfile = Mapfile;
            end
            if ~obj.Mapfile.isValid; return; end;
            
            if exist('getFile.m', 'file') == 2
                [obj.fileHandle obj.pathHandle] = getFile('*.txt', 'Open an impedance file...');
            else
                [obj.fileHandle obj.pathHandle] = uigetfile('*.txt', 'Open an impedance file...');
            end
            if ~obj.fileHandle; disp('No file was selected.'); return; end;
            impedanceDataCell = importdata([obj.pathHandle obj.fileHandle], ' ', 200);
            impedanceDataCell(1:9) = [];
            impedanceDataCell(97:end) = [];
            impedanceDataCellParsed = regexp(impedanceDataCell, '\t', 'split');
            for i = 1:size(impedanceDataCellParsed, 1)
                impedanceSingleCell = impedanceDataCellParsed{i,:}(2);
                impedanceSingleText = impedanceSingleCell{:};
                obj.impedanceDataValues(i) = str2double(impedanceSingleText(1:end-4));
            end
        end
    end
    methods
        function FilePath = PathName(obj)
            FilePath = obj.pathHandle;
        end
        function FileName = Filename(obj)
            FileName = obj.fileHandle;
        end
        function validFlag = isValid(obj)
            if any(isnan(obj.getChannelImpedances)) 
                validFlag = 0;
            else
                validFlag = 1;
            end
        end
        function impedanceDataValues = getImpedances(obj)
            impedanceDataValues = obj.impedanceDataValues;
        end
        function impedanceDataValue = getChannelImpedance(obj, chanNum)
            if ~exist('chanNum', 'var')
                disp('Channel number is a required argument.');
                return;
            end
            impedanceDataValue = obj.impedanceDataValues(chanNum);
        end
        function analyzedStat = getImpedanceMean(obj)
            analyzedStat = mean(obj.impedanceDataValues);
        end
        function analyzedStat = getImpedanceMedian(obj)
            analyzedStat = median(obj.impedanceDataValues);
        end
        function analyzedStat = getImpedanceRange(obj)
            analyzedStat = obj.getImpedanceMax-obj.getImpedanceMin;
        end
        function analyzedStat = getImpedanceSTD(obj)
            analyzedStat = std(obj.impedanceDataValues);
        end
        function analyzedStat = getImpedance975Perc(obj)
            sortedImpedances = sort(obj.impedanceDataValues, 'ascend');
            analyzedStat = sortedImpedances(round(obj.chanCount*0.975));
        end
        function analyzedStat = getImpedance750Perc(obj)
            sortedImpedances = sort(obj.impedanceDataValues, 'ascend');
            analyzedStat = sortedImpedances(round(obj.chanCount*0.75));
        end
        function analyzedStat = getImpedance250Perc(obj)
            sortedImpedances = sort(obj.impedanceDataValues, 'ascend');
            analyzedStat = sortedImpedances(round(obj.chanCount*0.25));
        end
        function analyzedStat = getImpedance25Perc(obj)
            sortedImpedances = sort(obj.impedanceDataValues, 'ascend');
            analyzedStat = sortedImpedances(round(obj.chanCount*0.025));
        end
        function analyzedStat = getImpedanceMin(obj)
            analyzedStat = min(obj.impedanceDataValues);
        end
        function analyzedStat = getImpedanceMax(obj)
            analyzedStat = max(obj.impedanceDataValues);
        end
        function plotStatistics(obj)
            figure('Name', 'Array Statistics');
            title('Array Impedance Values Statistics');
            set(gcf, 'Position', [657 540 288 287]);
            axis off;
            text(0,0.9, 'Mean');                                                                                                             text(0.6,0.9, num2str(obj.getImpedanceMean, '%5.1f'));
            text(0,0.8, 'Median');            text(0.6,0.8, num2str(obj.getImpedanceMean, '%5.1f'));
            text(0,0.7, 'Range');             text(0.6,0.7, num2str(obj.getImpedanceRange, '%5.1f'));
            text(0,0.6, 'STD');               text(0.6,0.6, num2str(obj.getImpedanceSTD, '%5.1f'));
            text(0,0.5, '97.5th Percentile'); text(0.6,0.5, num2str(obj.getImpedance975Perc, '%5.1f'));
            text(0,0.4, '75th Percentile');   text(0.6,0.4, num2str(obj.getImpedance750Perc, '%5.1f'));
            text(0,0.3, '25th Percentile');   text(0.6,0.3, num2str(obj.getImpedance250Perc, '%5.1f'));
            text(0,0.2, '2.5th Percentile');  text(0.6,0.2, num2str(obj.getImpedance25Perc, '%5.1f'));
            text(0,0.1, 'Min');               text(0.6,0.1, num2str(obj.getImpedanceMin, '%5.1f'));
            text(0,0.0, 'Max');               text(0.6,0.0, num2str(obj.getImpedanceMax, '%5.1f'));
        end
        function plotHistogram(obj, binCount, maxRange)
            if ~exist('binCount', 'var')
                binCount = 20;
            end
            if ~exist('maxRange', 'var')
                if max(obj.impedanceDataValues) > 10000
                    maxRange = 10000;
                else
                    maxRange = max(obj.impedanceDataValues);
                end
            end
            hist(obj.impedanceDataValues, min(obj.impedanceDataValues):...
                                          maxRange/binCount:...
                                          maxRange,1);
        title('Impedance Values Histogram');
        xlabel('Impedance Values (kOhm)');
        ylabel('Impedance Occurances');
        end
        function plotBoxPlot(obj, boxStyle)
            if ~exist('boxStyle', 'var');
                boxStyle = 'traditional';
            end
            sortedImpedances = sort(obj.impedanceDataValues);

            boxplot(obj.impedanceDataValues);
            title(['Box Plot (' boxStyle ')']);
            ylabel('Impedance Values (kOhm)');
            set(gca, 'xTickLabel', ' ');
            xlim([0.9, 1.1]);
            set(gcf, 'Position', [632 356 235 398]);             
        end
        function plotImpedances(obj, yelThreshold, redThreshold)
            if ~exist('yelThreshold', 'var')
                yelThreshold = 100;
            end
            if ~exist('redThreshold', 'var')
                redThreshold = 800;
            end
            plotFigure = KTFigure;
            plotFigure.EnlargeFigure;
            plotFigure.MakeBackgroundWhite;
            hold on;
            for channelIDX = 1:obj.chanCount
                if obj.impedanceDataValues(channelIDX) > redThreshold
                    spColor = [1,0,0];
                    fColor  = [1,1,1];
                elseif obj.impedanceDataValues(channelIDX) < redThreshold && obj.impedanceDataValues(channelIDX) > yelThreshold
                    spColor = [0,0,0];
                    fColor  = [1,1,1];
                elseif obj.impedanceDataValues(channelIDX) < yelThreshold
                    spColor = [1,1,0];
                    fColor  = [0,0,0];
                end
                obj.Mapfile.GenerateChannelSubplot(channelIDX);
                obj.Mapfile.GenerateChannelSubplotNames(channelIDX, fColor);
                text(0.05, 0.7, [num2str(obj.impedanceDataValues(channelIDX)) ' kOhm'], 'Color', fColor, 'FontSize', 10, 'FontWeight', 'bold');
                axis on;
                obj.Mapfile.setAxisBackgroundColor(spColor);
                obj.Mapfile.setAxesColor(spColor);
            end
            hold off;
            plotFigure.ShowFigure;
        end
        function plotImpedancesColorBackground(obj)
            for channelIDX = 1:obj.chanCount
                obj.Mapfile.GenerateChannelSubplot(channelIDX);
                axis on;
                obj.Mapfile.setAxisBackgroundColor([0,obj.impedanceDataValues(channelIDX)/max(max(obj.impedanceDataValues)),0]);        
            end
        end
        function plotImpedancesComparison(obj)
            changeThreshold = 0.2;
            yelThreshold = 100;
            secImpedanceFile = KTUEAImpedanceFile(obj.Mapfile);
            secImpedanceValues = secImpedanceFile.getChannelImpedances;
            plotFigure = KTFigure;
            plotFigure.EnlargeFigure;
            plotFigure.MakeBackgroundWhite;
            hold on;
            for channelIDX = 1:obj.chanCount
                percentChange = (obj.impedanceDataValues(channelIDX) - secImpedanceValues(channelIDX))/obj.impedanceDataValues(channelIDX);
                if abs(percentChange)  > changeThreshold && percentChange > 0
                    % Lime Green -- If impedance is increased by more 
                    % than threshold.
                    if abs(percentChange) > 2 * changeThreshold
                        spColor = [34,139,34] ./255;
                        fColor  = [1,1,1];
                    else
                        spColor = [154,205,50] ./255;
                        fColor  = [0,0,0];
                    end
                elseif abs(percentChange) > changeThreshold && percentChange < 0
                    % Chocolate orange -- If impedance is dropped by more 
                    % than threshold.
                    if abs(percentChange) > 2 * changeThreshold
                        spColor = [210,105,30] ./255;
                        fColor  = [1,1,1];
                    else
                        spColor = [245,222,179] ./255;
                        fColor  = [0,0,0];
                    end
                else
                    % Honeydew -- If threshold change is less than
                    % threshold.
                    spColor = [240,255,240] ./255;
                    fColor  = [0,0,0];
                end
                obj.Mapfile.GenerateChannelSubplot(channelIDX);
                patchHandle = patch([0.5, 0.5, 0, 0], [1, 0, 0, 1], spColor);
                obj.Mapfile.GenerateChannelSubplotNames(channelIDX, fColor);
                text(0.05, 0.7, ['1: ' num2str(obj.impedanceDataValues(channelIDX)) ' kOhm'], 'Color', fColor, 'FontSize', 10, 'FontWeight', 'bold');
                text(0.05, 0.5, ['2: ' num2str(secImpedanceValues(channelIDX)) ' kOhm'], 'Color', fColor, 'FontSize', 10, 'FontWeight', 'bold');
            end
            hold off;
            plotFigure.ShowFigure;
        end
    end
end