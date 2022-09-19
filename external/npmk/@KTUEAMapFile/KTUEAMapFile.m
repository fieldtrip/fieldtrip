classdef KTUEAMapFile

% UEAMAPFILE -- Defines a class that contains information on a UEA
% electrode array using CMP files provided by Blackrock Microsystems.
%
% To load a CMP file type 'MapName = UEAMapFile' where MapName is the name
% of the variable that will be created and will contain the map class.
%
% Example: myArrayMap = UEAMapFile;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METHODS:
%
% There are several available methods for mapfiles. Here's a list and a
% description of what they perform.
%
% Electrode2Channel(Electrode): 
%   Will return the Channel number corresponding to a passed Electrode.
%
% Channel2Electrode(Channel):   
%   Will return the Electrode number corresponding to a passed Channel.
%
% GetChannelBankID(Channel):
%   Returns the BankID for a passed Channel.
%
% GetChannelPin(Channel):
%   Returns the Pin for a passed Channel.
%
% GetChannelLabel(Channel):
%   Returns the Label for a passed Channel.
%
% GenerateChannelSubplot(Channel):
%   Returns a subplot handle for the passed Channel to be used in plotting
%   UEA-like maps.
%
% GenerateChannelSubplotNames(Channel):
%   Plots the names of channels and electrodes on the subplot corresponding
%   to the passed Channel.
%
% GetChannelColumnRow(Channel):
%   Returns the Column and Row positions of the passed number that is used
%   to plot UEA-like maps.
%
% PlotCMP:
%   Will plot a map of the CMP Map with all the electrode and channel
%   names.
%
%   Kian Torab
%   Blackrock Microsystems
%   ktorab@blackrockmicro.com
%
%   Version 1.6.1.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
% 1.6.1.0:
%   - Minor bug fix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
    properties (SetAccess = private)
        ElecNum
        ChanNum
        Column
        Row
        Bank
        Pin
        Label
        pathHandle
        fileHandle
    end
    methods (Hidden)
        function obj = KTUEAMapFile(fileName)
            if ~exist('fileName', 'var')
                if exist('getFile.m', 'file') == 2
                    [obj.fileHandle obj.pathHandle] = getFile('*.cmp', 'Open a CMP Mapfile...');
                else
                    [obj.fileHandle obj.pathHandle] = uigetfile('*.cmp', 'Open a CMP Mapfile...');
                end
                if ~obj.fileHandle
                    disp('No file was selected.'); 
                    return; 
                end
                mapfileDataCell = importdata([obj.pathHandle obj.fileHandle], ' ', 200);
            else
                mapfileDataCell = importdata(fileName, ' ', 200);
            end
            mapfileDataCellParsed = regexp(mapfileDataCell, '\t', 'split');
            if size(mapfileDataCellParsed, 2) == 1
                mapfileDataCellParsed = regexp(mapfileDataCell, '\s+', 'split');
            end
            lastJunkIDX = 1;
            while 1
                if size(mapfileDataCellParsed{lastJunkIDX}, 2) == 6 || size(mapfileDataCellParsed{lastJunkIDX}, 2) == 5
                    if ~isempty(str2num(mapfileDataCellParsed{lastJunkIDX}{1})) && ...
                    ~isempty(str2num(mapfileDataCellParsed{lastJunkIDX}{2})) && ...
                    ~isempty(str2num(mapfileDataCellParsed{lastJunkIDX}{4}))
                    lastJunkIDX = lastJunkIDX - 1;
                        break;
                    end
                end
                lastJunkIDX = lastJunkIDX + 1;
            end            
            mapfileDataCellParsed(1:lastJunkIDX) = [];
            for i = 1:size(mapfileDataCellParsed, 1)
                if ~all(isspace(mapfileDataCellParsed{i}))
                    obj.Column(i) = str2num(mapfileDataCellParsed{i,:}{1});
                    obj.Row(i)    = str2num(mapfileDataCellParsed{i,:}{2});
                    obj.Bank(i)   = mapfileDataCellParsed{i,:}{3}-'@';
                    obj.Pin(i)    = str2num(mapfileDataCellParsed{i,:}{4});
                    ElectNum      = str2num(mapfileDataCellParsed{i,:}{5}(5:end));
                    if isempty(ElectNum)
                        obj.ElecNum(i) = str2num(mapfileDataCellParsed{i,:}{5}(2:end-4));
                    else
                        obj.ElecNum(i) = ElectNum;
                    end
                    obj.ChanNum(i)  = (obj.Bank(i) - 1) * 32 + obj.Pin(i);
                    obj.Label{i}    = mapfileDataCellParsed{i,:}{5};
                end
            end
        end
    end
    methods
        function validFlag = isValid(obj)
            if ~obj.getFilename
                validFlag = 0;
            else
                validFlag = 1;
            end
        end
        function FilePath = getPathName(obj)
            FilePath = obj.pathHandle;
        end
        function FileName = getFilename(obj)
            FileName = obj.fileHandle;
        end
        function Channel = Electrode2Channel(obj, Electrode)
            if ~exist('Electrode', 'var'); disp('Electrode is a required input. Refer to help for more information.'); return; end;
            if isempty(obj.ChanNum(obj.ElecNum == Electrode(1)))
                disp('Electrode number is invalid.');
                Channel = [];
                return;
            end
            for ElecIDX = 1:length(Electrode)
                Channel(ElecIDX) = obj.ChanNum(obj.ElecNum == Electrode(ElecIDX));
            end
        end
        function Electrode = Channel2Electrode(obj, Channel)
            if isempty(obj.ElecNum(obj.ChanNum == Channel(1)))
                disp('Channel number is invalid.');
                Electrode = [];
                return;
            end            
            if ~exist('Channel', 'var'); disp('Channel is a required input. Refer to help for more information.'); return; end;
            for ChanIDX = 1:length(Channel)
                Electrode(ChanIDX) = obj.ElecNum(obj.ChanNum == Channel(ChanIDX));
            end
        end
        function BankID = getChannelBankID(obj, Channel)
            BankID = obj.Bank(obj.ChanNum == Channel);
        end
        function Label = getChannelLabel(obj, Channel)
            Label = char(obj.Label{obj.ChanNum == Channel});
        end
        function Pin = getChannelPin(obj, Channel)
            Pin = obj.Pin(obj.ChanNum == Channel);
        end
        function spHandle = GenerateChannelSubplot(obj, Channel)
            if ~exist('Channel', 'var'); disp('Channel is a required input. Refer to help for more information.'); return; end;
            [x, y] = getChannelColumnRow(obj, Channel);
            spHandle=subplot('position', [x/10+0.004, y/10+0.004, 1/10-0.004, 1/10-0.004]);
            axis off;
            box on;
            set(spHandle,'XTickLabel', []);
            set(spHandle,'YTickLabel', []);
            %%
        end
        function GenerateChannelSubplotNames(obj, Channel, fColor)
            if ~exist('Channel', 'var'); disp('Channel is a required input. Refer to help for more information.'); return; end;
            if ~exist('fColor', 'var'); fColor  = [1,1,1]; end
            Electrode = Channel2Electrode(obj, Channel);
            text(0.13, 0.23, ['elec ', num2str(Electrode)], 'Color', fColor, 'FontSize', 8);
            text(0.13, 0.10, ['chan ', num2str(Channel)], 'Color', fColor, 'FontSize', 8);            
        end
        function [x, y] = getChannelColumnRow(obj, Channel)
            if ~exist('Channel', 'var'); disp('Channel is a required input. Refer to help for more information.'); return; end;
            x = obj.Column(obj.ChanNum == Channel);
            y = obj.Row(obj.ChanNum == Channel);
        end
        function PlotCMP(obj)
            figure;
            for Channel = 1:96
                GenerateChannelSubplot(obj, Channel);
                GenerateChannelSubplotNames(obj, Channel, [0,0,0]);
            end
        end
        function setAxisBackgroundColor(~, curColor)
            if isa(curColor, 'double') && (length(curColor) == 3)
                if all(curColor <= 1) && all(curColor >= 0)
                    set(gca, 'Color', curColor);
                else
                    disp('The color values have to be between 0 and 1.');
                end
            elseif isa(curColor, 'char')
                try
                	set(gca, 'Color', curColor);
                catch
                    disp('Color is not valid. Type ''doc color'' for more information.');
                end
            else                    
                disp('The input needs to be a valid color input.');
                disp('A valid color input is a vector of length 3 of values between 0 and 1.');
            end
        end
        function setXAxisColor(~, curColor)
            if isa(curColor, 'double') && (length(curColor) == 3)
                if all(curColor <= 1) && all(curColor >= 0)
                    set(gca, 'XColor', curColor);
                else
                    disp('The color values have to be between 0 and 1.');
                end
            elseif isa(curColor, 'char')
                try
                    set(gca, 'XColor', curColor);
                catch
                    disp('Color is not valid. Type ''doc color'' for more information.');
                end
            else                    
                disp('The input needs to be a valid color input.');
                disp('A valid color input is a vector of length 3 of values between 0 and 1.');
            end
        end
        function setYAxisColor(~, curColor)
            if isa(curColor, 'double') && (length(curColor) == 3)
                if all(curColor <= 1) && all(curColor >= 0)
                    set(gca, 'YColor', curColor);
                else
                    disp('The color values have to be between 0 and 1.');
                end
            elseif isa(curColor, 'char')
                try
                    set(gca, 'YColor', curColor);
                catch
                    disp('Color is not valid. Type ''doc color'' for more information.');
                end
            else                    
                disp('The input needs to be a valid color input.');
                disp('A valid color input is a vector of length 3 of values between 0 and 1.');
            end
        end
        function setAxesColor(~, curColor)
            if isa(curColor, 'double') && (length(curColor) == 3)
                if all(curColor <= 1) && all(curColor >= 0)
                    set(gca, 'YColor', curColor);
                    set(gca, 'XColor', curColor);
                else
                    disp('The color values have to be between 0 and 1.');
                end
            elseif isa(curColor, 'char')
                try
                    set(gca, 'YColor', curColor);
                    set(gca, 'XColor', curColor);
                catch
                    disp('Color is not valid. Type ''doc color'' for more information.');
                end
            else                    
                disp('The input needs to be a valid color input.');
                disp('A valid color input is a vector of length 3 of values between 0 and 1.');
            end
        end        
    end
end