% =============================================================
% Simple demo program calling functions from the
% SGLX_readMeta class.
%
function DemoReadSGLXData()

% Ask user for binary file
[binName,path] = uigetfile('*.bin', 'Select Binary File');

% Parse the corresponding metafile
meta = SGLX_readMeta.ReadMeta(binName, path);

% Get first one second of data
nSamp = floor(1.0 * SGLX_readMeta.SampRate(meta));
dataArray = SGLX_readMeta.ReadBin(0, nSamp, meta, binName, path);
dataType = 'A';         %set to 'A' for analog, 'D' for digital data

% For an analog channel: gain correct saved channel ch (1-based for MATLAB).
ch = 1;

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;

% Read these lines in dw (0-based).
% For 3B2 imec data: the sync pulse is stored in line 6.
% May be 1 or more line indices.
dLineList = [0,1,6];

if dataType == 'A'
    if strcmp(meta.typeThis, 'imec')
        dataArray = SGLX_readMeta.GainCorrectIM(dataArray, [ch], meta);
    else
        dataArray = SGLX_readMeta.GainCorrectNI(dataArray, [ch], meta);
    end
    plot(1e6*dataArray(ch,:));
else
    digArray = SGLX_readMeta.ExtractDigital(dataArray, meta, dw, dLineList);
    for i = 1:numel(dLineList)
        plot(digArray(i,:));
        hold on
    end
    hold off
end
end % DemoReadSGLXData