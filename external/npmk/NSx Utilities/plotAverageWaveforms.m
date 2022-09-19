% plotAverageWaveforms
%
% Plots average waveforms for all channels and specified units across the
% array. It will prompt for a CMP file to map each channel to its electrode
% representation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Use plotAverageWaveforms(NEV, units)
% 
% NOTE: All input arguments are optional. Input arguments may be in any order.
%
%   NEV:          NEV structure containing all waveforms.
%                 DEFAULT: Will prompt for NEV.
%
%   'units':      Specify the units to plot (i.e. 2 or 0:5 for all).
%                 DEFAULT: will plot all units.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   USAGE EXAMPLE: 
%   
%   openNEV(NEV, 0:3);
%
%   In the example above, the waveforms in the NEV structure coming from
%   units 0 through 3 will be plotted.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   kian@blackrockmicro.com
%   Blackrock Microsystems
%   Salt Lake City, UT
%   
%   Version 2.1.0.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
% 2.1.0.0:
%   - Added the ability to plot scaled or non-scaled plots. By default, the
%   plots are scaled.
%   - The color of the first unit is always selected at random, so when a
%   channel only has 1 unit, that 1 unit will be a different color
%   (visualization only).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotAverageWaveforms(NEV, units, plottype)

mapfile = 1;
colors = repmat({'magenta', 'cyan', 'red', 'green', 'blue', 'black'}, 1, 3);

if ~exist('NEV', 'var')
    NEV = openNEV('read');
end

if ~exist('units', 'var')
    units = 0:5;
end

if ~exist('plottype', 'var')
    plottype = 'scaled';
end

numberOfElectrodes = length(unique(NEV.Data.Spikes.Electrode));

disp('Select a map file for your array. If none, press cancel or ESC.');
mapFile = KTUEAMapFile;
if ~mapFile.isValid
    disp('Map File was not selected.');
    mapfile = 0;
    plotHeight = floor(sqrt(numberOfElectrodes));
    plotWidth = ceil(numberOfElectrodes/floor(sqrt(numberOfElectrodes)));
end

figure = KTFigure;

figure.EnlargeFigure;
figure.MakeBackgroundWhite;
figure.SetActive;
title('Spikes Average Waveforms for Individual Channels');

for chanIDX = 1:numberOfElectrodes
    
    if mapfile
try
        mapFile.GenerateChannelSubplot(chanIDX);
catch
    disp('e');
end
    else
        subplot(plotHeight, plotWidth, chanIDX);
    end
    
    hold on;
    initialColor = round(rand(1)*5);
    for unitIDX = units
        [~, allWaveforms] = findEventTimes(NEV, chanIDX, unitIDX);
        if ~isempty(allWaveforms)
            averageWaveform = mean(allWaveforms, 2);
            plot(averageWaveform, 'color', colors{unitIDX+1+initialColor}, 'LineWidth',2);
        end
    end
    hold off;
    
    if strcmpi(plottype, 'scaled')
        ylim([-1000, 1000]);
    end
%     if size(allWaveforms,2) ~= 0
%         text(1, 0, num2str(size(allWaveforms,2)), 'color', 'red');
%     end
    mapFile.setAxisBackgroundColor('white');
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end